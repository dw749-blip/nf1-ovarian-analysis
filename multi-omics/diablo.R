library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(mixOmics)
library(caret)
library(IlluminaHumanMethylation27kanno.ilmn12.hg19)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(patchwork)
library(stringr)
setwd('E:/NF1/TCGA/diablo')

# ===============================
# Load NF1 groups (only NF1-null and NF1-wildtype)
# ===============================
df <- read.delim("NF1_expression_functional_groups.tsv", sep = "\t", stringsAsFactors = FALSE)
df <- df %>%
  tibble::rownames_to_column(var = "sample")

# NF1-null
NF1_null_samples <- df$sample[
  df$NF1_functional_status == "NF1-null"]

# Wild-type
NF1_wildtype_samples <- df$sample[
  df$NF1_functional_status == "NF1-wildtype" &
    df$group != "GTEx"]

selected_samples <- c(NF1_null_samples, NF1_wildtype_samples)
group_labels <- c(
  rep("NF1_null", length(NF1_null_samples)),
  rep("NF1_wildtype", length(NF1_wildtype_samples))
)

# ===============================
# Load data
# ===============================
# ------------- Gene expression -------------
rna_raw <- read.delim("TcgaGtex_ovary_gene_expected_counts.tsv", stringsAsFactors = FALSE)
colnames(rna_raw) <- gsub("\\.", "-", colnames(rna_raw))
rna_raw <- rna_raw[,selected_samples]

# Convert Ensembl IDs to HUGO symbols
ensembl_id_clean <- sub("\\..*$", "", rownames(rna_raw))
rownames(rna_raw) <- ensembl_id_clean
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_id_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL")
)
gene_map <- gene_map %>% filter(!is.na(SYMBOL) & SYMBOL != "")
rna_raw <- rna_raw[rownames(rna_raw) %in% gene_map$ENSEMBL, ]
gene_map <- gene_map[match(rownames(rna_raw), gene_map$ENSEMBL), ]
rna_raw$SYMBOL <- gene_map$SYMBOL
rna_raw <- rna_raw %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  dplyr::select(-ENSEMBL) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")
rna_raw <- rna_raw %>% column_to_rownames("SYMBOL")

# Normalize and filter to top 10000 variable genes
coldata <- data.frame(
  row.names = colnames(rna_raw),
  NF1_group = group_labels[match(colnames(rna_raw), selected_samples)]
)
dds <- DESeqDataSetFromMatrix(
  countData = rna_raw,
  colData = coldata,
  design = ~ NF1_group
)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
rna_vst <- assay(vsd)
gene_var <- apply(rna_vst, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:10000]
rna_var <- rna_vst[top_genes, ]

# ------------- DNA methylation -------------
methy_raw <- read.delim("TCGA-OV.methylation27.tsv.gz", stringsAsFactors = FALSE)
colnames(methy_raw) <- gsub("\\.", "-", colnames(methy_raw))
colnames(methy_raw) <- gsub("([0-9]{2})[A-Z]$", "\\1", colnames(methy_raw))
methy_raw <- methy_raw %>% column_to_rownames("Composite-Element-REF")
methy_raw <- methy_raw[, intersect(colnames(methy_raw), selected_samples)]

# Filter to top 5000 variable probes with no missing values
methy_complete <- methy_raw[complete.cases(methy_raw), ]
#methy_var_values <- apply(methy_complete, 1, var)
#top_methy <- names(sort(methy_var_values, decreasing = TRUE))[1:5000]
#methy_var <- methy_complete[top_methy, ]

# ------------- miRNA -------------
mirna_raw <- read.delim("TCGA-OV.mirna.tsv.gz", stringsAsFactors = FALSE)
colnames(mirna_raw) <- gsub("\\.", "-", colnames(mirna_raw))
colnames(mirna_raw) <- gsub("([0-9]{2})[A-Z]$", "\\1", colnames(mirna_raw))
mirna_raw <- mirna_raw %>% column_to_rownames("miRNA_ID")
mirna_raw <- mirna_raw[, intersect(colnames(mirna_raw), selected_samples)]

# Keep only miRNAs with non-zero variance
mirna_var_values <- apply(mirna_raw, 1, var)
mirna_var <- mirna_raw[mirna_var_values > 0, ]

# ------------- Protein expression -------------
protein_raw <- read.delim("TCGA-OV.protein.tsv", stringsAsFactors = FALSE)
colnames(protein_raw) <- gsub("\\.", "-", colnames(protein_raw))
colnames(protein_raw) <- gsub("([0-9]{2})[A-Z]$", "\\1", colnames(protein_raw))
protein_raw <- protein_raw %>% column_to_rownames("peptide_target")
protein_raw <- protein_raw[, intersect(colnames(protein_raw), selected_samples)]

# ===============================
# Filter samples to those with all data types
# ===============================
common_samples <- Reduce(intersect, list(colnames(rna_var), colnames(methy_complete), colnames(mirna_var), colnames(protein_raw)))
length(common_samples)

# Count samples per NF1 group
group_df <- data.frame(sample = selected_samples, NF1_group = group_labels, stringsAsFactors = FALSE)
group_common <- group_df %>% 
  filter(sample %in% common_samples) %>% 
  arrange(match(sample, common_samples))
table(group_common$NF1_group)

# Create outcome variable (NF1 group)
Y <- factor(
  group_common$NF1_group,
  levels = c("NF1_wildtype", "NF1_null")
)
names(Y) <- group_common$sample
table(Y)

# Make all matrices have same sample order
rna_block <- t(rna_var[, group_common$sample])
methy_block <- t(methy_complete[, group_common$sample])
mirna_block <- t(mirna_var[, group_common$sample])
protein_block <- t(protein_raw[, group_common$sample])

# ===============================
# Build DIABLO input
# ===============================
data <- list(
  RNA = rna_block,
  Methylation = methy_block,
  miRNA = mirna_block,
  Protein = protein_block
)

# Remove near-zero variance features
remove_nzv <- function(df) {
  nzv <- nearZeroVar(df, saveMetrics = TRUE)
  df[, !nzv$nzv, drop = FALSE]
}
data$RNA <- remove_nzv(data$RNA)
data$Methylation <- remove_nzv(data$Methylation)
data$miRNA <- remove_nzv(data$miRNA)
data$Protein <- remove_nzv(data$Protein)

# ===============================
# Initial pairwise PLS comparisons
# ===============================
list.keepX = c(25, 25)
list.keepY = c(25, 25)

# Generate pairwise PLS models
pls1 <- spls(data[["RNA"]], data[["Methylation"]], 
             keepX = list.keepX, keepY = list.keepY) 
pls2 <- spls(data[["RNA"]], data[["miRNA"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["RNA"]], data[["Protein"]], 
             keepX = list.keepX, keepY = list.keepY)
pls4 <- spls(data[["Methylation"]], data[["miRNA"]], 
             keepX = list.keepX, keepY = list.keepY)
pls5 <- spls(data[["Methylation"]], data[["Protein"]], 
             keepX = list.keepX, keepY = list.keepY)
pls6 <- spls(data[["miRNA"]], data[["Protein"]], 
             keepX = list.keepX, keepY = list.keepY)
plotVar(pls1, cutoff = 0.5, title = "(a) RNA vs Methylation", 
        legend = c("RNA", "Methylation"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
plotVar(pls2, cutoff = 0.5, title = "(b) RNA vs miRNA", 
        legend = c("RNA", "miRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
plotVar(pls3, cutoff = 0.5, title = "(c) RNA vs Protein", 
        legend = c("RNA", "Protein"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
plotVar(pls4, cutoff = 0.5, title = "(c) Methylation vs miRNA", 
        legend = c("Methylation", "miRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
plotVar(pls5, cutoff = 0.5, title = "(c) Methylation vs Protein", 
        legend = c("Methylation", "Protein"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
plotVar(pls6, cutoff = 0.5, title = "(c) miRNA vs Protein", 
        legend = c("miRNA", "Protein"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))
cor(pls1$variates$X, pls1$variates$Y) 
cor(pls2$variates$X, pls2$variates$Y) 
cor(pls3$variates$X, pls3$variates$Y) 
cor(pls4$variates$X, pls4$variates$Y) 
cor(pls5$variates$X, pls5$variates$Y) 
cor(pls6$variates$X, pls6$variates$Y) 

# ===============================
# Build DIABLO model for high predictive ability (design value = 0.1)
# ===============================
# Design matrix
design <- matrix(
  0.1,
  ncol = length(data),
  nrow = length(data),
  dimnames = list(names(data), names(data))
)
diag(design) <- 0

# Tune the number of components
basic_diablo_model <- block.splsda(X = data, Y = Y, ncomp = 5, design = design)
perf_diablo = perf(basic_diablo_model, validation = 'Mfold', folds = 5, nrepeat = 5)
plot(perf_diablo)
perf_diablo$choice.ncomp$WeightedVote
ncomp <- 2

# Tuning of parameters
test.keepX <- list(
  RNA = c(10, 25, 50, 100, 150, 200),
  Methylation = c(10, 25, 50, 100, 150, 200),
  miRNA = c(10, 20, 30, 50, 100),
  Protein = c(5, 10, 25, 50)
)
tune_diablo <- tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                                 test.keepX = test.keepX, design = design,
                                 validation = "Mfold", folds = 5, nrepeat = 5,
                                 dist = "centroids.dist")
list.keepX <- tune_diablo$choice.keepX
list.keepX

# Build model
diablo_model <- block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
diablo_model$design

selectVar(diablo_model, block = 'RNA', comp = 1)$RNA$name
selectVar(diablo_model, block = 'Methylation', comp = 1)$Methylation$name
selectVar(diablo_model, block = 'miRNA', comp = 1)$miRNA$name
selectVar(diablo_model, block = 'Protein', comp = 1)$Protein$name

# ===============================
# Save DIABLO objects for reproducibility
# ===============================
saveRDS(diablo_model, file = "diablo_model_changed.rds")
saveRDS(tune_diablo, file = "tune_diablo_changed.rds")
saveRDS(perf_diablo, file = "perf_diablo_changed.rds")
saveRDS(list.keepX, file = "selected_keepX_changed.rds")

# Read in objects to save
diablo_model <- readRDS("diablo_model_changed.rds")
tune_diablo <- readRDS("tune_diablo_changed.rds")
perf_diablo <- readRDS("perf_diablo_changed.rds")
list.keepX <- readRDS("selected_keepX_changed.rds")
ncomp <- 2

# ===============================
# Clean and save all selected features
# ===============================
# Compute group means for each block to identify higher expression
group_vec <- Y[common_samples]
get_higher_expression <- function(feature, block_matrix) {
  if(!(feature %in% rownames(block_matrix))) return(NA)
  
  values <- block_matrix[feature, ]
  values <- suppressWarnings(as.numeric(values))
  group_vec_aligned <- group_vec[colnames(block_matrix)]
  wt_mean <- mean(values[group_vec_aligned == "NF1_wildtype"], na.rm = TRUE)
  null_mean <- mean(values[group_vec_aligned == "NF1_null"], na.rm = TRUE)
  
  if(is.na(wt_mean) | is.na(null_mean)) return(NA)
  if(null_mean > wt_mean) return("NF1_null")
  if(wt_mean > null_mean) return("NF1_wildtype")
  return("Equal")
}

components <- 1:2
feature_list <- lapply(components, function(comp){
  # RNA
  rna_genes <- selectVar(diablo_model, block = "RNA", comp = comp)$RNA$name
  rna_loadings <- diablo_model$loadings$RNA[, comp]
  
  # Methylation
  cpgs <- selectVar(diablo_model, block = "Methylation", comp = comp)$Methylation$name
  cpg_loadings <- diablo_model$loadings$Methylation[, comp]
  anno27_subset <- anno27[cpgs, ]
  cpg_genes <- anno27_subset$Symbol
  
  # Protein
  proteins <- selectVar(diablo_model, block = "Protein", comp = comp)$Protein$name
  protein_loadings <- diablo_model$loadings$Protein[, comp]
  proteins_clean <- gsub("_.*", "", proteins)
  proteins_clean <- gsub("BETACATENIN", "CTNNB1", proteins_clean)
  proteins_clean <- gsub("Lasu1", "HUWE1", proteins_clean)
  proteins_clean <- gsub("GCN5L2", "KAT2A", proteins_clean)
  proteins_clean <- gsub("EVI1", "MECOM", proteins_clean)
  proteins_clean <- gsub("CD20", "MS4A1", proteins_clean)
  proteins_clean <- gsub("PKA-a", "PRKACA", proteins_clean)
  proteins_clean <- gsub("b-Actin", "ACTB", proteins_clean)
  proteins_clean <- gsub("PAI1", "SERPINE1", proteins_clean)
  proteins_clean <- gsub("ECADHERIN", "CDH1", proteins_clean)
  proteins_clean <- gsub("HMHA1", "ARHGAP45", proteins_clean)
  proteins_clean <- gsub("ERALPHA", "ESR1", proteins_clean)
  proteins_clean <- gsub("FIBRONECTIN", "FN1", proteins_clean)
  proteins_clean <- gsub("BAP1C4", "BAP1", proteins_clean)
  proteins_clean <- gsub("Tau", "MAPT", proteins_clean)
  proteins_clean <- gsub("IGFRb", "IGFBP1", proteins_clean)
  proteins_clean <- gsub("ANNEXIN1", "ANXA1", proteins_clean)
  
  # miRNA
  mirna_features <- selectVar(diablo_model, block = "miRNA", comp = comp)$miRNA$name
  mirna_loadings <- diablo_model$loadings$miRNA[, comp]
  
  # Dataframes
  df_rna <- data.frame(
    Component = comp,
    OmicsType = "RNA",
    Feature = rna_genes,
    GeneSymbol = rna_genes,
    Loading = rna_loadings[rna_genes],
    stringsAsFactors = FALSE)
  df_cpg <- data.frame(
    Component = comp,
    OmicsType = "CpG",
    Feature = cpgs,
    GeneSymbol = cpg_genes,
    Loading = cpg_loadings[cpgs],
    stringsAsFactors = FALSE)
  df_protein <- data.frame(
    Component = comp,
    OmicsType = "Protein",
    Feature = proteins,
    GeneSymbol = proteins_clean,
    Loading = protein_loadings[proteins],
    stringsAsFactors = FALSE)
  df_mirna <- data.frame(
    Component = comp,
    OmicsType = "miRNA",
    Feature = mirna_features,
    GeneSymbol = mirna_features,
    Loading = mirna_loadings[mirna_features],
    stringsAsFactors = FALSE)
  
  df_all <- bind_rows(df_rna, df_cpg, df_protein, df_mirna)
  
  # Phenotype associated with loading weight
  df_all$AssociatedPhenotype <- ifelse(
    df_all$Loading < 0,
    "NF1_null",
    "NF1_wildtype"
  )
  
  # Higher expression group
  df_all$HigherExpression <- NA
  
  df_all$HigherExpression[df_all$OmicsType == "RNA"] <-
    sapply(df_all$Feature[df_all$OmicsType == "RNA"],
           get_higher_expression,
           block_matrix = rna_var)
  
  df_all$HigherExpression[df_all$OmicsType == "CpG"] <-
    sapply(df_all$Feature[df_all$OmicsType == "CpG"],
           get_higher_expression,
           block_matrix = methy_complete)
  
  df_all$HigherExpression[df_all$OmicsType == "Protein"] <-
    sapply(df_all$Feature[df_all$OmicsType == "Protein"],
           get_higher_expression,
           block_matrix = protein_raw)
  
  df_all$HigherExpression[df_all$OmicsType == "miRNA"] <-
    sapply(df_all$Feature[df_all$OmicsType == "miRNA"],
           get_higher_expression,
           block_matrix = mirna_var)
  
  return(df_all)
})

top_features_all <- bind_rows(feature_list)

# Save
write.csv(top_features_all, file = "diablo_all_top_features.csv", row.names = FALSE)

# ===============================
# Model performance
# ===============================
perf.diablo = perf(diablo_model, validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 
perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate

# ===============================
# Plot results
# ===============================
pdf("sample_clusters.pdf", width = 10, height = 8)
group <- diablo_model$Y
cols <- c("NF1_wildtype" = "#FF7F0E", "NF1_null" = "#1F77B4")
shapes <- c("NF1_wildtype" = 17, "NF1_null" = 16)
plotIndiv(diablo_model, col = cols, pch = shapes,
          ind.names = FALSE, legend = TRUE, legend.title = "Phenotype")
dev.off()

# ------------- Circos plot -------------
pdf("circos_plot.pdf", width = 8, height = 8)
circosPlot(
  diablo_model,
  comp = 1,
  cutoff = 0.7,
  line = TRUE,
  size.variables = 0.5
)
dev.off()

# ------------- Heatmap -------------
pdf("diablo_heatmap.pdf", width = 8, height = 6)
row_order <- c("NF1_wildtype", "NF1_null")
group_colors <- c(
  "NF1_wildtype" = "#1F77B4",
  "NF1_null"  = "#FF7F0E"
)
block_order <- c("mRNA", "CpGs", "miRNA", "Proteins")
block_colors <- c(
  "mRNA"     = "#377eb8",
  "CpGs"     = "#4daf4a",
  "miRNA"    = "#e41a1c",
  "Proteins" = "#FF7F0E"
)
cimDiablo(
  diablo_model,
  comp = 1,
  cluster = "column",
  legend = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  size.legend = 1
)
legend(
  "right",
  legend = row_order,
  fill   = group_colors[row_order],
  border = NA,
  bty    = "n",
  title = "Rows",
  cex    = 1
)
legend(
  x = "topright",
  inset = c(0.02, 0),
  legend = block_order,
  fill = block_colors[block_order],
  border = NA,
  bty = "n",
  title = "Columns",
  cex = 1
)
dev.off()

# ------------- Feature plots -------------
# Manual protein name corrections
protein_map <- c(
  "BETACATENIN" = "CTNNB1",
  "Lasu1" = "HUWE1",
  "GCN5L2" = "KAT2A",
  "EVI1" = "MECOM",
  "CD20" = "MS4A1",
  "PKA-a" = "PRKACA",
  "b-Actin" = "ACTB",
  "PAI1" = "SERPINE1",
  "ECADHERIN" = "CDH1",
  "HMHA1" = "ARHGAP45",
  "ERALPHA" = "ESR1",
  "FIBRONECTIN" = "FN1",
  "BAP1C4" = "BAP1",
  "Tau" = "MAPT",
  "IGFRb" = "IGFBP1",
  "ANNEXIN1" = "ANXA1"
)

msig_all <- msigdbr(species = "Homo sapiens")

ovarian_sets <- msig_all %>% filter(grepl("OVARIAN", gs_name, ignore.case = TRUE))
ovarian_msigdb_genes <- ovarian_sets %>% pull(gene_symbol) %>% unique()
ovarian_msigdb_genes <- toupper(trimws(ovarian_msigdb_genes))

nf1_sets <- msig_all %>% filter(grepl("NF1|NEUROFIBROMIN", gs_name, ignore.case = TRUE))
nf1_genes <- unique(nf1_sets$gene_symbol)
nf1_genes <- toupper(trimws(nf1_genes))

mircancer <- read.csv(file = "mircancer.txt", sep = "\t", header = FALSE)
ovarian_mir <- mircancer %>%
  filter(grepl("ovarian", V2, ignore.case = TRUE))
ovarian_mirnas <- unique(ovarian_mir$V1)

blocks <- setdiff(names(diablo_model$loadings), "Y")
feature_df <- do.call(rbind, lapply(blocks, function(block) {
  vars <- selectVar(diablo_model, block = block, comp = 1)[[block]]
  features <- vars$name
  loadings <- abs(vars$value)
  
  # Annotate methylation features
  if(block == "Methylation"){
    anno27 <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
    anno_subset <- anno27[features, ]
    gene_symbols <- anno_subset$Symbol
    gene_symbols[gene_symbols == "" | is.na(gene_symbols)] <- "NA"
    features <- paste0(features, " (", gene_symbols, ")")
  }
  # Annotate protein features
  if(block == "Protein"){
    features <- gsub("_.*", "", features)
    gene_symbols <- protein_map[features]
    gene_symbols[is.na(gene_symbols)] <- features[is.na(gene_symbols)]
    features <- paste0(features, " (", gene_symbols, ")")
  }
  data.frame(
    Feature = features,
    Loading = loadings,
    Block = block,
    stringsAsFactors = FALSE
  )
}))

feature_df <- feature_df[order(feature_df$Block, -feature_df$value.var), ]
feature_df$Feature <- factor(feature_df$Feature,
                             levels = unique(feature_df$Feature))
feature_df$Feature <- as.character(feature_df$Feature)
feature_df$Block <- factor(feature_df$Block,
                           levels = c("RNA", "Methylation", "miRNA", "Protein"))

# Add overlap with databases
feature_df$GeneSymbol <- NA
# RNA block
feature_df$GeneSymbol[feature_df$Block == "RNA"] <- 
  feature_df$Feature[feature_df$Block == "RNA"]
# Methylation block (extract gene inside parentheses)
feature_df$GeneSymbol[feature_df$Block == "Methylation"] <- sub(".*\\((.*)\\)", "\\1", feature_df$Feature[feature_df$Block == "Methylation"])
# Protein block (extract gene inside parentheses)
feature_df$GeneSymbol[feature_df$Block == "Protein"] <- sub(".*\\((.*)\\)", "\\1", feature_df$Feature[feature_df$Block == "Protein"])
feature_df$GeneSymbol[feature_df$Block == "Protein" & feature_df$GeneSymbol %in% names(protein_map)] <-
  protein_map[feature_df$GeneSymbol[feature_df$Block == "Protein" &
                                      feature_df$GeneSymbol %in% names(protein_map)]]
# miRNA block
feature_df$GeneSymbol[feature_df$Block == "miRNA"] <-
  feature_df$Feature[feature_df$Block == "miRNA"]
feature_df$In_miRCancer <- FALSE
feature_df$In_miRCancer[feature_df$Block == "miRNA"] <-
  feature_df$GeneSymbol[feature_df$Block == "miRNA"] %in% ovarian_mirnas

feature_df$GeneSymbol <- toupper(trimws(feature_df$GeneSymbol))
feature_df$Known_OC <- feature_df$GeneSymbol %in% ovarian_msigdb_genes
feature_df$In_NF1 <- feature_df$GeneSymbol %in% nf1_genes

# Remove OvCa-related, NF1-related classifications
make_block_plot <- function(data, block_name) {
  title_map <- c(
    "RNA" = "mRNA",
    "Methylation" = "CpGs",
    "miRNA" = "miRNA",
    "Protein" = "Proteins"
  )
  df_block <- data %>% 
    filter(Block == block_name)
  ggplot(df_block, aes(x = reorder(Feature, value.var), y = value.var)) +
    geom_segment(aes(xend = reorder(Feature, value.var), y = 0, yend = value.var),
                 linewidth = 0.8, color = "black") +
    geom_point(size = 3, color = "black") +
    coord_flip() +
    labs(
      title = title_map[block_name],
      y = "Absolute Loading Weight",
      x = "Features") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none")}

make_block_plot <- function(data, block_name) {
  title_map <- c(
    "RNA" = "mRNA",
    "Methylation" = "CpGs",
    "miRNA" = "miRNA",
    "Protein" = "Proteins"
  )
  
  df_block <- data %>% 
    filter(Block == block_name)
  # Create database label
  if(block_name == "miRNA"){
    df_block$Database <- "No overlap"
    df_block$Database[df_block$In_miRCancer] <- "miRCancer"
    df_block$Database <- factor(
      df_block$Database,
      levels = c("No overlap", "miRCancer"))
    colors <- c(
      "No overlap" = "black",
      "miRCancer"  = "green")
  } else {
    df_block$Database <- "No overlap"
    df_block$Database[df_block$Known_OC] <- "OvCa-related"
    df_block$Database[df_block$In_NF1] <- "NF1-related"
    df_block$Database[df_block$Known_OC & df_block$In_NF1] <- "Both"
    df_block$Database <- factor(
      df_block$Database,
      levels = c("No overlap", "OvCa-related", "NF1-related", "Both"))
    colors <- c(
      "No overlap" = "black",
      "OvCa-related" = "red",
      "NF1-related" = "blue",
      "Both" = "purple"
    )
  }
  
  ggplot(df_block,
         aes(x = reorder(Feature, value.var), y = value.var)) +
    geom_segment(aes(xend = reorder(Feature, value.var),
                     y = 0,
                     yend = value.var,
                     color = Database),
                 linewidth = 0.8) +
    geom_point(aes(color = Database),
               size = 3) +
    coord_flip() +
    scale_color_manual(
      name = "Pathway overlap",
      values = colors,
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        byrow = TRUE,
        keyheight = unit(0.4, "cm"),
        keywidth  = unit(0.4, "cm")
      )
    ) +
    labs(
      title = title_map[block_name],
      y = "Absolute Loading Weight",
      x = "Features"
    ) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.key = element_blank(),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

plot_RNA <- make_block_plot(feature_df, "RNA")
plot_Methylation <- make_block_plot(feature_df, "Methylation")
plot_miRNA <- make_block_plot(feature_df, "miRNA")
plot_Protein <- make_block_plot(feature_df, "Protein")

combined_plot <- (plot_RNA | plot_Methylation) / (plot_miRNA | plot_Protein)
combined_plot
ggsave("feature_plots.pdf",
       combined_plot,
       width = 10,
       height = 10)

# ------------- Perform enrichment of selected features -------------
# Do only features associated with NF1-null (negative loading weight, component = 1)
get_negative_features <- function(model, block){
  loadings_vec <- model$loadings[[block]][, 1]
  neg_features <- names(loadings_vec)[loadings_vec < 0]
  return(neg_features)
}

#get_negative_features <- function(model, block){
#  loadings_mat <- model$loadings[[block]][, 1:2]
#  neg_features <- rownames(loadings_mat)[apply(loadings_mat, 1, function(x) any(x < 0))]
#  return(neg_features)
#}

rna_neg_features <- get_negative_features(diablo_model, "RNA")
cpg_neg_features <- get_negative_features(diablo_model, "Methylation")
protein_neg_features <- get_negative_features(diablo_model, "Protein")

# Convert to gene symbols
rna_genes_neg <- rna_neg_features

anno27 <- getAnnotation(IlluminaHumanMethylation27kanno.ilmn12.hg19)
anno27_subset <- anno27[cpg_neg_features, ]
cpg_genes_neg <- anno27_subset$Symbol
cpg_genes_neg <- unique(na.omit(cpg_genes_neg))

protein_genes_neg <- protein_map[protein_neg_features]
protein_genes_neg[is.na(protein_genes_neg)] <- protein_neg_features[is.na(protein_genes_neg)]
protein_genes_neg <- unique(protein_genes_neg)

combined_genes_neg <- unique(c(rna_genes_neg, cpg_genes_neg, protein_genes_neg))

# Hallmark pathways
msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)
hallmark_list <- split(msig_h$gene_symbol, msig_h$gs_name)

run_hallmark_enrichment <- function(gene_vector){
  gene_vector <- unique(na.omit(gene_vector))
  enricher(
    gene = gene_vector,
    TERM2GENE = msig_h[, c("gs_name", "gene_symbol")],
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}
hallmark_rna <- run_hallmark_enrichment(rna_genes_neg)
hallmark_cpg <- run_hallmark_enrichment(cpg_genes_neg)
hallmark_protein <- run_hallmark_enrichment(protein_genes_neg)
hallmark_combined <- run_hallmark_enrichment(combined_genes_neg)

write.csv(as.data.frame(hallmark_rna), "hallmark_enrichment_rna.csv", row.names = FALSE)
write.csv(as.data.frame(hallmark_cpg), "hallmark_enrichment_cpg.csv", row.names = FALSE)
write.csv(as.data.frame(hallmark_protein), "hallmark_enrichment_protein.csv", row.names = FALSE)
write.csv(as.data.frame(hallmark_combined), "hallmark_enrichment_combined.csv", row.names = FALSE)

make_hallmark_barplot <- function(enrich_obj, title_text){
  df <- as.data.frame(enrich_obj)
  # Order by p.adjust
  df <- df %>%
    arrange(p.adjust)
  # Add significance category
  df$Significant <- ifelse(df$p.adjust < 0.05, "Significant", "Not significant")
  # Convert to factor for plotting
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df$log10_padj <- -log10(df$p.adjust)
  ggplot(df, aes(x = log10_padj, y = Description, fill = Significant)) +
    geom_col() +
    scale_fill_manual(
      values = c(
        "Significant" = "red",
        "Not significant" = "gray")) +
    labs(
      title = title_text,
      x = "-log10(p-adjust)",
      y = "Pathway"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9))
}

plot_rna <- make_hallmark_barplot(hallmark_rna, "mRNA")
plot_cpg <- make_hallmark_barplot(hallmark_cpg, "CpGs")
plot_protein <- make_hallmark_barplot(hallmark_protein, "Proteins")
plot_combined <- make_hallmark_barplot(hallmark_combined, "Hallmark Pathway Enrichment")
plot_combined
combined_plot <- plot_rna / plot_cpg / plot_protein / plot_combined
combined_plot
ggsave("hallmark_enrichment_combined_plot.pdf", plot_combined, width = 10, height = 8)

# GO BP terms (use Hallmark pathways above instead)
# RNA
go_rna <- enrichGO(
  gene = rna_genes_neg,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
go_df <- as.data.frame(go_rna)
sig_terms <- go_df %>%
  dplyr::filter(p.adjust < 0.05)
write.csv(sig_terms, file = "go_enrichment_rna.csv", row.names = FALSE)

# Methylation
go_cpg <- enrichGO(
  gene = cpg_genes_neg,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
go_df_cpg <- as.data.frame(go_cpg)
sig_terms_cpg <- go_df_cpg %>%
  dplyr::filter(p.adjust < 0.05)
write.csv(sig_terms_cpg, file = "go_enrichment_cpg.csv", row.names = FALSE)

# Protein
go_protein <- enrichGO(
  gene = protein_genes_neg,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
go_df_protein <- as.data.frame(go_protein)
sig_terms_protein <- go_df_protein %>%
  dplyr::filter(p.adjust < 0.05)
write.csv(sig_terms_protein, file = "go_enrichment_protein.csv", row.names = FALSE)

# RNA, methylation, protein combined
go_combined <- enrichGO(
  gene = combined_genes_neg,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
go_df_combined <- as.data.frame(go_combined)
sig_terms_combined <- go_df_combined %>%
  dplyr::filter(p.adjust < 0.05)
write.csv(sig_terms_combined, file = "go_enrichment_combined.csv", row.names = FALSE)

# GO barplots
prepare_go_df <- function(go_obj, modality_label) {
  df <- data.frame(as.data.frame(go_obj), stringsAsFactors = FALSE)
  df$GeneRatio_numeric <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  df$Modality <- modality_label
  df <- df[df$p.adjust < 0.05, ]
  df <- df[order(df$p.adjust), ]
  head(df, 20)
}

rna_df <- prepare_go_df(go_rna, "mRNA")
cpg_df <- prepare_go_df(go_cpg, "CpGs")
protein_df <- prepare_go_df(go_protein, "Proteins")
combined_df <- prepare_go_df(go_combined, "Proteins")
rna_df$log10_padj <- -log10(rna_df$p.adjust)
cpg_df$log10_padj <- -log10(cpg_df$p.adjust)
protein_df$log10_padj <- -log10(protein_df$p.adjust)
combined_df$log10_padj <- -log10(combined_df$p.adjust)

make_dotplot <- function(df, title_text) {
  ggplot(df,
         aes(x = GeneRatio_numeric,
             y = reorder(Description, GeneRatio_numeric),
             size = Count,
             color = log10_padj)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#4575b4",
      high = "#d73027",
      name = "-log10(p-adj)",
    ) +
    labs(title = title_text,
         x = "Gene Ratio",
         y = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
}
plot_rna <- make_dotplot(rna_df, "mRNA")
plot_cpg <- make_dotplot(cpg_df, "CpGs")
plot_protein <- make_dotplot(protein_df, "Proteins")
plot_combined <- make_dotplot(combined_df, "GO BP Pathway Enrichment")
combined_plot <- plot_rna / plot_cpg / plot_protein / plot_combined
combined_plot
plot_combined
ggsave("go_enrichment_combined_plot.pdf", plot_combined, width = 8, height = 7)

# ------------- Loading plots and tables (to get which features drive which NF1 group) -------------
plotLoadings(diablo_model, comp = 1, contrib = 'max', method = 'median')

# Extract component scores
scores <- diablo_model$X$RNA[, 1]
# Get group labels
group <- diablo_model$Y
# Combine
df_scores <- data.frame(score = scores, group = group)
aggregate(score ~ group, data = df_scores, median)
positive_group <- "NF1_wildtype"
negative_group <- "NF1_null"

extract_loadings_table <- function(model, block_name, comp = 1,
                                   positive_group, negative_group) {
  # Get loadings
  load_vec <- model$loadings[[block_name]][, comp]
  
  # Keep only selected (non-zero)
  load_vec <- load_vec[load_vec != 0]
  
  # Create dataframe
  df <- data.frame(
    Block = block_name,
    Feature = names(load_vec),
    Component = comp,
    Loading = load_vec,
    Direction = ifelse(load_vec > 0, "Positive", "Negative"),
    Associated_Group = ifelse(load_vec > 0,
                              positive_group,
                              negative_group),
    stringsAsFactors = FALSE
  )
  
  # Sort by absolute contribution
  df <- df[order(abs(df$Loading), decreasing = TRUE), ]
  
  return(df)
}

rna_tab     <- extract_loadings_table(diablo_model, "RNA",     1,
                                      positive_group, negative_group)
mirna_tab   <- extract_loadings_table(diablo_model, "miRNA",   1,
                                      positive_group, negative_group)
cpg_tab     <- extract_loadings_table(diablo_model, "Methylation",     1,
                                      positive_group, negative_group)
protein_tab <- extract_loadings_table(diablo_model, "Protein", 1,
                                      positive_group, negative_group)