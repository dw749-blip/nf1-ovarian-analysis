library(DESeq2)
library(WGCNA)
library(tidyverse)
library(dplyr)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(ggrepel)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(grid)
setwd('E:/NF1/TCGA/')
allowWGCNAThreads()

# ===============================
# Load RNA-seq counts
# ===============================
# Load filtered raw counts
raw_counts <- read.delim("TcgaGtex_ovary_gene_expected_counts.tsv", sep = "\t", stringsAsFactors = FALSE)
colnames(raw_counts) <- gsub("\\.", "-", colnames(raw_counts))

# ===============================
# Load NF1 groups
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

meta_df <- data.frame(
  sample = selected_samples,
  NF1_group = group_labels,
  stringsAsFactors = FALSE
)
rownames(meta_df) <- meta_df$sample

# ===============================
# Subset count matrix
# ===============================
raw_counts <- raw_counts[, colnames(raw_counts) %in% selected_samples]
raw_counts <- raw_counts[, meta_df$sample]

# ===============================
# DESeq2 normalization (VST)
# ===============================
dds <- DESeqDataSetFromMatrix(
  countData = round(raw_counts),
  colData = meta_df,
  design = ~ NF1_group
)

# Filter low-expression genes (counts < 10 in more than 90% of samples)
count_mat <- counts(dds)
keep_genes <- rowMeans(count_mat >= 10) >= 0.10
dds <- dds[keep_genes, ]

vsd <- vst(dds, blind = FALSE)
expr_vst <- assay(vsd)

# ===============================
# Dataset filtering
# ===============================
gene_var <- apply(expr_vst, 1, mad)
sorted_genes <- sort(gene_var, decreasing = TRUE)
top_genes <- head(names(sorted_genes), 10000) # Keep top 10k most variable genes
expr_vst <- expr_vst[top_genes, ]

# ===============================
# Prepare data for WGCNA
# ===============================
expr_df <- as.data.frame(t(expr_vst))
expr_df[1:5, 1:10]

gsg <- goodSamplesGenes(expr_df, verbose = 3)
if (!gsg$allOK) {
  expr_df <- expr_df[gsg$goodSamples, gsg$goodGenes]
}

# ===============================
# Pick soft-threshold power
# ===============================
powers <- c(1:20)
sft <- pickSoftThreshold(expr_df, powerVector = powers, verbose = 5, networkType = "signed")

# Plot scale independence and mean connectivity graphs to determine ideal SFT
par(mar=c(1,1,1,1))
plot(sft$fitIndices[,1],sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],sft$fitIndices[,2],col="red")
abline(h=0.90,col="red")

par(mar=c(1,1,1,1))
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels= sft$fitIndices[,1],col="red")

softPower <- 9

# ===============================
# Build network & detect modules
# ===============================
network <- blockwiseModules(
    expr_df,
    power = softPower,
    networkType = "signed",
    TOMType = "signed",
    minModuleSize = 30,
    mergeCutHeight = 0.25,
    saveTOMs = TRUE,
    saveTOMFileBase = "NF1",
    pamRespectsDendro = FALSE,
    reassignThreshold = 0,
    numericLabels = TRUE,
    verbose = 3
)

# ===============================
# Trait matrix
# ===============================
traitData <- meta_df[rownames(expr_df), , drop = FALSE]

traitData$NF1_status <- ifelse(traitData$NF1_group == "NF1_null", 1, 0)
traitMatrix <- traitData[, "NF1_status", drop = FALSE]

# ===============================
# Heatmap of module–trait relationships
# ===============================
moduleColors <- labels2colors(network$colors)
MEs <- orderMEs(network$MEs)
moduleTraitCor <- cor(MEs, traitMatrix, use = "p")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(expr_df))
# Adjust p-values
moduleTraitFDR <- matrix(
  p.adjust(as.vector(moduleTraitP), method = "BH"),
  nrow = nrow(moduleTraitP),
  ncol = ncol(moduleTraitP)
)
rownames(moduleTraitFDR) <- rownames(moduleTraitP)
colnames(moduleTraitFDR) <- colnames(moduleTraitP)

pdf("wgcna_heatmap.pdf", width = 5, height = 8)
# Round correlations
cor_rounded <- round(moduleTraitCor, 2)
# Create text matrix with significance stars
get_stars <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else ""
}
text_matrix <- matrix("", nrow = nrow(moduleTraitCor), ncol = ncol(moduleTraitCor))
for (i in 1:nrow(moduleTraitCor)) {
  for (j in 1:ncol(moduleTraitCor)) {
    stars <- get_stars(moduleTraitFDR[i, j])
    text_matrix[i, j] <- paste0(cor_rounded[i, j], stars)
  }
}
# Create color scale/palette
wgcna_cols <- WGCNA::blueWhiteRed(100)
col_fun <- colorRamp2(seq(-1, 1, length.out = 100), wgcna_cols)
# Plot
wgcna_heatmap <- Heatmap(
  moduleTraitCor,
  name = "Correlation",
  col = col_fun,
  border = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_title = "Module",
  heatmap_legend_param = list(
    legend_height = unit(8, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (moduleTraitFDR[i, j] < 0.05) {
      grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 12, fontface = "bold"))
    } else {
      grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 12))
    }})
draw(wgcna_heatmap)
dev.off()

moduleLabels <- network$colors
names(moduleLabels) <- colnames(expr_df)

# Module membership (MM)
geneMM <- cor(expr_df, MEs, use = "p")
colnames(geneMM) <- colnames(MEs)

# Gene significance (GS)
geneGS <- cor(expr_df, traitMatrix, use = "p")
colnames(geneGS) <- paste0("GS_", colnames(traitMatrix))

# ===============================
# Function to identify significant modules and hub genes for any NF1 group
# ===============================
get_modules_and_hubs <- function(trait_col, moduleTraitFDR, moduleLabels,
                                 geneMM, geneGS, orgDb,
                                 kME_thresh = 0.8, GS_thresh = 0.3,
                                 min_genes_for_enrichment = 10) {
  # Significant modules for this group
  sig_modules <- rownames(moduleTraitFDR)[moduleTraitFDR[, trait_col] < 0.05]
  if (length(sig_modules) == 0) return(NULL)
  
  # Identify hub genes per module
  hub_list <- list()
  for (mod in sig_modules) {
    mod_label <- as.numeric(sub("ME", "", mod))
    mod_genes <- names(moduleLabels)[moduleLabels == mod_label]
    
    # Columns to use
    kME_col <- mod
    GS_col  <- "GS_NF1_status"
    
    hubs <- mod_genes[
      abs(geneMM[mod_genes, kME_col]) >= kME_thresh &
        abs(geneGS[mod_genes, GS_col]) >= GS_thresh
    ]
    
    if (length(hubs) > 0) hub_list[[mod]] <- hubs
  }
  
  # Convert hubs to a tidy data frame with HUGO symbols
  if (length(hub_list) > 0) {
    hub_df <- tibble(
      module = rep(names(hub_list), lengths(hub_list)),
      ensembl_id = unlist(hub_list)
    ) %>%
      mutate(ensembl_id_clean = sub("\\..*$", "", ensembl_id))
    
    # Keep only valid keys
    valid_keys <- AnnotationDbi::keys(orgDb, keytype = "ENSEMBL")
    hub_df <- hub_df %>%
      filter(ensembl_id_clean %in% valid_keys)
    
    if (nrow(hub_df) == 0) {
      hub_df_final <- tibble(module = character(), ensembl_id = character(), SYMBOL = character())
    } else {
      gene_map <- AnnotationDbi::select(
        orgDb,
        keys = unique(hub_df$ensembl_id_clean),
        keytype = "ENSEMBL",
        columns = c("SYMBOL", "GENENAME")
      )
      hub_df_final <- hub_df %>%
        left_join(gene_map, by = c("ensembl_id_clean" = "ENSEMBL")) %>%
        arrange(module)
    }
  } else {
    hub_df_final <- tibble(module = character(), ensembl_id = character(), SYMBOL = character())
  }
  
  # Prepare gene lists per module for enrichment
  module_gene_list <- lapply(sig_modules, function(mod) {
    mod_label <- as.numeric(sub("ME", "", mod))
    genes <- names(moduleLabels)[moduleLabels == mod_label]
    sub("\\..*$", "", genes) # clean Ensembl IDs
  })
  names(module_gene_list) <- sig_modules
  
  # Map to Entrez IDs
  safe_select <- function(orgDb, keys, keytype = "ENSEMBL", columns = c("SYMBOL", "GENENAME")) {
    tryCatch(
      AnnotationDbi::select(orgDb, keys = keys, keytype = keytype, columns = columns),
      error = function(e) {
        message("Error in select()!")
        message("Number of keys: ", length(keys))
        message("First 20 keys: ", paste(head(keys, 20), collapse = ", "))
        message("Any NAs? ", any(is.na(keys)))
        stop(e)  # re-throw the error after printing info
      }
    )
  }
  gene2entrez <- safe_select(
    orgDb,
    keys = unique(unlist(module_gene_list)),
    keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  )
  map_to_entrez <- function(genes) {
    gene2entrez %>%
      filter(ENSEMBL %in% genes, !is.na(ENTREZID)) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID)
  }
  
  # GO enrichment
  go_results <- lapply(names(module_gene_list), function(mod) {
    entrez <- map_to_entrez(module_gene_list[[mod]])
    if (length(entrez) < min_genes_for_enrichment) return(NULL)
    enrichGO(
      gene = entrez,
      OrgDb = orgDb,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.10,
      readable = TRUE
    )
  })
  names(go_results) <- names(module_gene_list)
  
  return(list(
    sig_modules = sig_modules,
    hub_genes = hub_df_final,
    module_gene_list = module_gene_list,
    go_results = go_results
  ))
}

# ===============================
# Apply to all NF1 groups
# ===============================
results <- get_modules_and_hubs(
  trait_col = "NF1_status",
  moduleTraitFDR = moduleTraitFDR,
  moduleLabels = moduleLabels,
  geneMM = geneMM,
  geneGS = geneGS,
  orgDb = org.Hs.eg.db
)

# ===============================
# Helper functions to extract top terms and hub genes
# ===============================
get_top_terms <- function(enrich_list, module, n = 3) {
  if (!module %in% names(enrich_list)) return(NA_character_)
  enrich_obj <- enrich_list[[module]]
  if (is.null(enrich_obj)) return(NA_character_)
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(NA_character_)
  paste(
    df %>%
      arrange(p.adjust) %>%
      slice_head(n = n) %>%
      pull(Description),
    collapse = "; "
  )
}

get_top_hub_genes <- function(module, trait, moduleLabels, geneMM, geneGS, gene_map, n = 3) {
  ensembl2symbol <- gene_map %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENSEMBL, SYMBOL) %>%
    deframe()
  
  mod_label <- as.numeric(sub("ME", "", module))
  mod_genes <- names(moduleLabels)[moduleLabels == mod_label]
  
  if (length(mod_genes) == 0) return(NA_character_)
  
  kME_col <- module
  GS_col  <- paste0("GS_", trait)
  
  hub_df <- tibble(
    ensembl = mod_genes,
    kME  = geneMM[mod_genes, kME_col],
    GS   = geneGS[mod_genes, GS_col]
  ) %>%
    drop_na() %>%
    mutate(
      hub_score = abs(kME) * abs(GS),
      SYMBOL = ensembl2symbol[sub("\\..*$", "", ensembl)]
    ) %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(desc(hub_score))
  
  if (nrow(hub_df) == 0) return(NA_character_)
  
  paste(
    head(hub_df$SYMBOL, n),
    collapse = "; "
  )
}

# ===============================
# Function to build a module table for one NF1 group
# ===============================
build_module_table <- function(moduleTraitCor, moduleTraitP, moduleTraitFDR,
                               go_results, moduleLabels, geneMM, geneGS, gene_map,
                               p_cutoff = 0.05, n_terms = 3, n_hub = 3) {
  sig_modules <- rownames(moduleTraitFDR)[moduleTraitFDR[, "NF1_status"] < p_cutoff]
  if (length(sig_modules) == 0) return(NULL)
  tibble(
    module = sig_modules,
    association = ifelse(moduleTraitCor[sig_modules, "NF1_status"] > 0, "NF1_null", "NF1_wildtype"),
    correlation = moduleTraitCor[sig_modules, "NF1_status"],
    p_adj = moduleTraitFDR[sig_modules, "NF1_status"],
    GO_BP = map_chr(sig_modules, ~ get_top_terms(go_results, .x, n_terms)),
    top_hub_genes = map_chr(sig_modules, ~ get_top_hub_genes(.x, "NF1_status", moduleLabels, geneMM, geneGS, gene_map, n_hub))
  ) %>%
    arrange(p_adj)}

# ===============================
# Build combined table for all NF1 groups
# ===============================
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(unlist(results$module_gene_list)),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)
final_table <- build_module_table(
  moduleTraitCor = moduleTraitCor,
  moduleTraitP = moduleTraitP,
  moduleTraitFDR = moduleTraitFDR,
  go_results = results$go_results,
  moduleLabels = moduleLabels,
  geneMM = geneMM,
  geneGS = geneGS,
  gene_map = gene_map,
  p_cutoff = 0.05,
  n_terms = 3,
  n_hub = 3
)

write.csv(final_table, file = "wgcna_go_enrichment_hub_genes.csv", row.names = FALSE)

# ===============================
# Plot dotplot with correlation, significance, and enrichment results
# ===============================
get_top_terms_simple <- function(enrich_list, module, n = 2) {
  enrich_obj <- enrich_list[[module]]
  if (is.null(enrich_obj)) return(NA_character_)
  df <- as.data.frame(enrich_obj)
  if (nrow(df) == 0) return(NA_character_)
  terms <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    pull(Description)
  paste(terms, collapse = " | ")
}

manual_labels <- c(
  "ME2" = "oxidative phosphorylation | cellular respiration",
  "ME9" = "small GTPase-mediated signal transduction",
  "ME4" = "extracellular matrix organization",
  "ME6" = "cilium organization | microtubule-based movement",
  "ME1" = "leukocyte cell-cell adhesion | leukocyte migration",
  "ME12" = "cytoplasmic translation | ribosome biogenesis",
  "ME16" = "immunoglobulin mediated immune response",
  "ME18" = "oligosaccharide metabolic process",
  "ME20" = "No enriched GO terms",
  "ME14" = "No enriched GO terms"
)

plot_df <- tibble(
  module = rownames(moduleTraitCor),
  correlation = round(moduleTraitCor[, "NF1_status"], 2),
  p_value = moduleTraitP[, "NF1_status"],
  p_adj = moduleTraitFDR[, "NF1_status"]
) %>%
  mutate(
    neg_log10_p = -log10(p_adj),
    NF1_association = ifelse(correlation > 0, "NF1_null", "NF1_wildtype"))

plot_df <- plot_df %>%
  mutate(
    GO_terms = map_chr(module, ~ get_top_terms_simple(results$go_results, .x, n = 2))
  )

plot_df <- plot_df %>%
  mutate(
    GO_manual = manual_labels[module],
    label = ifelse(!is.na(GO_manual), paste0(module, ": ", GO_manual), module))

plot_df_sig <- plot_df %>%
  filter(p_adj < 0.05)

pdf("wgcna_dotplot.pdf", width = 11, height = 8)
ggplot(plot_df_sig, aes(x = correlation, y = reorder(label, correlation))) +
  geom_point(
    aes(size = neg_log10_p, fill = correlation),
    shape = 21, color = "black", stroke = 0.3
  ) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    limits = c(-0.5, 0.5),
    breaks = c(-0.5, 0, 0.5)
  ) +
  scale_size(range = c(3, 10)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Correlation with NF1 status",
    y = "Module (Top GO terms)",
    size = "-log10(p-adj)",
    fill = "Correlation with NF1-null"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 14),
    legend.position = "right"
  )
dev.off()