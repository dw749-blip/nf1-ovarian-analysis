library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rstatix)
library(xCell)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(EPIC)
library(MCPcounter)
library(progeny)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
setwd('E:/NF1/TCGA/')
source("cibersort/CIBERSORT.R")

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

#selected_samples <- c(NF1_null_samples, NF1_het_del_samples, NF1_wildtype_samples)
selected_samples <- c(NF1_null_samples, NF1_wildtype_samples)
group_labels <- c(
  rep("NF1_null", length(NF1_null_samples)),
  #rep("NF1_het_del", length(NF1_het_del_samples)),
  rep("NF1_wildtype", length(NF1_wildtype_samples))
)

# ===============================
# Subset count matrix and normalize (TPM) for input
# ===============================
raw_counts <- raw_counts[, colnames(raw_counts) %in% selected_samples]
rownames(raw_counts) <- raw_counts[,1]
raw_counts <- raw_counts[,-1]

# Convert Ensembl IDs to HUGO symbols
ensembl_id_clean <- sub("\\..*$", "", rownames(raw_counts))
rownames(raw_counts) <- ensembl_id_clean
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_id_clean),
  keytype = "ENSEMBL",
  columns = c("SYMBOL", "GENENAME")
)
gene_map <- gene_map %>%
  filter(!is.na(SYMBOL) & SYMBOL != "")

# Get gene lengths
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_ranges <- genes(txdb)
gene_length_df <- tibble(
  ENTREZID = names(gene_ranges),
  gene_length = width(gene_ranges)
)
ensembl2entrez <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_id_clean),
  keytype = "ENSEMBL",
  columns = c("ENTREZID")
)
gene_annot <- gene_map %>%
  left_join(ensembl2entrez, by = c("ENSEMBL" = "ENSEMBL")) %>%
  left_join(gene_length_df, by = "ENTREZID") %>%
  filter(!is.na(gene_length))

# Normalize to non-log space (TPM)
counts_df <- raw_counts %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL")
counts_annot <- counts_df %>%
  inner_join(gene_annot, by = "ENSEMBL")

tpm_df <- counts_annot %>%
  dplyr::select(SYMBOL, gene_length, all_of(colnames(raw_counts))) %>%
  pivot_longer(
    cols = all_of(colnames(raw_counts)),
    names_to = "sample",
    values_to = "count"
  ) %>%
  mutate(rate = count / gene_length) %>%
  group_by(sample) %>%
  mutate(TPM = rate / sum(rate, na.rm = TRUE) * 1e6) %>%
  ungroup()

# Collapse duplicate gene symbols
tpm_collapsed <- tpm_df %>%
  group_by(SYMBOL, sample) %>%
  summarize(TPM = sum(TPM), .groups = "drop")

# Write and save final expression matrix
expr_tpm <- tpm_collapsed %>%
  pivot_wider(
    names_from = sample,
    values_from = TPM
  ) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

expr_tpm <- as.data.frame(expr_tpm) %>%
  tibble::rownames_to_column(var = "GeneSymbol")

write.table(expr_tpm, file = "cibersort/cibersort_input_tpm_normalized.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# ===============================
# xCell - test significant immune associations
# ===============================
expr_tpm_xcell <- expr_tpm %>%
  tibble::column_to_rownames("GeneSymbol") %>%
  as.matrix()

xcell_scores <- xCellAnalysis(expr_tpm_xcell)

# Create mapping from sample to NF1 group
group_df <- tibble(
  Mixture = selected_samples,
  NF1_group = group_labels
)

xcell_df <- xcell_scores %>%
  as.data.frame() %>%
  rownames_to_column("cell_type") %>%
  pivot_longer(
    -cell_type,
    names_to = "sample",
    values_to = "score"
  ) %>%
  left_join(group_df, by = c("sample" = "Mixture"))

# NF1 null vs. WT
xcell_null_vs_wt <- xcell_df %>%
  filter(NF1_group %in% c("NF1_null", "NF1_wildtype")) %>%
  group_by(cell_type) %>%
  summarise(
    median_null = median(score[NF1_group == "NF1_null"], na.rm = TRUE),
    median_wt   = median(score[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    mean_null   = mean(score[NF1_group == "NF1_null"], na.rm = TRUE),
    mean_wt     = mean(score[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    p = tryCatch(
      wilcox.test(score ~ NF1_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    comparison = "NF1_null_vs_WT",
    p.adj = p.adjust(p, method = "BH"),
    median_diff = median_null - median_wt,
    mean_diff = mean_null - mean_wt,
    enriched_in = case_when(
      median_diff > 0 ~ "NF1_null",
      median_diff < 0 ~ "NF1_wildtype",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p.adj)

# ===============================
# EPIC - test significant immune associations
# ===============================
epic_res <- EPIC(bulk = expr_tpm_xcell)
epic_scores <- epic_res$cellFractions

epic_df <- epic_scores %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    -sample,
    names_to = "cell_type",
    values_to = "fraction"
  ) %>%
  left_join(group_df, by = c("sample" = "Mixture"))

epic_null_vs_wt <- epic_df %>%
  filter(NF1_group %in% c("NF1_null", "NF1_wildtype")) %>%
  group_by(cell_type) %>%
  summarise(
    median_null = median(fraction[NF1_group == "NF1_null"], na.rm = TRUE),
    median_wt   = median(fraction[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    mean_null   = mean(fraction[NF1_group == "NF1_null"], na.rm = TRUE),
    mean_wt     = mean(fraction[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    p = tryCatch(
      wilcox.test(fraction ~ NF1_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    median_diff = median_null - median_wt,
    enriched_in = case_when(
      median_diff > 0 ~ "NF1_null",
      median_diff < 0 ~ "NF1_wildtype",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p.adj)

# ===============================
# MCP-counter - test significant immune associations
# ===============================
mcp_scores <- MCPcounter.estimate(expr_tpm_xcell, featuresType = "HUGO_symbols")

mcp_df <- mcp_scores %>%
  as.data.frame() %>%
  rownames_to_column("cell_type") %>%
  pivot_longer(
    -cell_type,
    names_to = "sample",
    values_to = "score"
  ) %>%
  left_join(group_df, by = c("sample" = "Mixture"))

mcp_null_vs_wt <- mcp_df %>%
  filter(NF1_group %in% c("NF1_null", "NF1_wildtype")) %>%
  group_by(cell_type) %>%
  summarise(
    median_null = median(score[NF1_group == "NF1_null"], na.rm = TRUE),
    median_wt   = median(score[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    mean_null   = mean(score[NF1_group == "NF1_null"], na.rm = TRUE),
    mean_wt     = mean(score[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    p = tryCatch(
      wilcox.test(score ~ NF1_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    median_diff = median_null - median_wt,
    enriched_in = case_when(
      median_diff > 0 ~ "NF1_null",
      median_diff < 0 ~ "NF1_wildtype",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(p.adj)

# ===============================
# PROGENy - test significant immune associations
# ===============================
progeny_scores <- progeny(
  expr_tpm_xcell,
  scale = TRUE,
  organism = "Human",
  top = 100,
  perm = 1,
  verbose = TRUE,
  z_scores = FALSE,
  get_nulldist = FALSE
)

progeny_df <- progeny_scores %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(
    -sample,
    names_to = "pathway",
    values_to = "activity"
  ) %>%
  left_join(group_df, by = c("sample" = "Mixture"))

progeny_null_vs_wt <- progeny_df %>%
  filter(NF1_group %in% c("NF1_null", "NF1_wildtype")) %>%
  group_by(pathway) %>%
  summarise(
    median_null = median(activity[NF1_group == "NF1_null"], na.rm = TRUE),
    median_wt   = median(activity[NF1_group == "NF1_wildtype"], na.rm = TRUE),
    p = tryCatch(
      wilcox.test(activity ~ NF1_group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    diff = median_null - median_wt,
    enriched_in = ifelse(diff > 0, "NF1_null", "NF1_wildtype")
  ) %>%
  arrange(p.adj)

# ===============================
# Create heatmap of aggregated immune results
# ===============================
direction_to_numeric <- function(df, feature_col = "cell_type") {
  df %>%
    mutate(
      value = case_when(
        p.adj < 0.05 & enriched_in == "NF1_wildtype" ~ 1,
        p.adj < 0.05 & enriched_in == "NF1_null" ~ -1,
        p.adj >= 0.05 ~ 0,
        TRUE ~ 0
      )
    ) %>%
    select(all_of(feature_col), value)
}

progeny_dir <- progeny_null_vs_wt %>%
  filter(p.adj < 0.05) %>%
  direction_to_numeric(feature_col = "pathway") %>%
  rename(feature = pathway) %>%
  mutate(tool = "PROGENy")
mcp_dir <- mcp_null_vs_wt %>%
  direction_to_numeric("cell_type") %>%
  rename(feature = cell_type) %>%
  mutate(tool = "MCP-counter")
epic_dir <- epic_null_vs_wt %>%
  direction_to_numeric("cell_type") %>%
  rename(feature = cell_type) %>%
  mutate(tool = "EPIC")
xcell_dir <- xcell_null_vs_wt %>%
  filter(cell_type != "Astrocytes") %>%
  direction_to_numeric("cell_type") %>%
  rename(feature = cell_type) %>%
  mutate(tool = "xCell")

# Harmonize feature names across tools
standardize_feature_names <- function(df) {
  df %>%
    mutate(feature = case_when(
      # Stromal
      feature %in% c("Fibroblasts", "CAFs") ~ "Fibroblasts",
      feature %in% c("Endothelial cells", "Endothelial") ~ "Endothelial cells",
      # Myeloid
      feature %in% c("Myeloid dendritic cells", "cDC") ~ "Dendritic cells",
      feature %in% c("MEP") ~ "Megakaryocyte-erythroid progenitors (MEPs)",
      # Lymphoid
      feature %in% c("T cells") ~ "CD3+ T cells",
      feature %in% c("CD8_Tcells", "CD8+ T-cells") ~ "CD8+ T cells",
      feature %in% c("NK cells", "NKcells") ~ "NK cells",
      feature %in% c("CLP") ~ "Common lymphoid progenitors (CLPs)",
      TRUE ~ feature))}
progeny_dir <- progeny_dir %>% rename(feature = feature)
mcp_dir <- standardize_feature_names(mcp_dir)
epic_dir <- standardize_feature_names(epic_dir)
xcell_dir <- standardize_feature_names(xcell_dir)

progeny_features <- progeny_dir$feature
stromal_features <- c(
  "Fibroblasts",
  "Endothelial cells",
  "Chondrocytes",
  "Pericytes")
myeloid_features <- c(
  "Monocytic lineage",
  "Macrophages",
  "Dendritic cells",
  "Neutrophils",
  "Erythrocytes",
  "Megakaryocyte-erythroid progenitors (MEPs)")
lymphoid_features <- c(
  "CD3+ T cells",
  "CD8+ T cells",
  "NK cells",
  "Common lymphoid progenitors (CLPs)")

feature_order <- c(
  stromal_features,
  myeloid_features,
  lymphoid_features,
  progeny_features)

all_dir <- bind_rows(
  progeny_dir %>% mutate(tool = "PROGENy"),
  mcp_dir %>% mutate(tool = "MCP-counter"),
  epic_dir %>% mutate(tool = "EPIC"),
  xcell_dir %>% mutate(tool = "xCell")
)

heatmap_df <- expand.grid(
  feature = feature_order,
  tool = c("xCell", "MCP-counter", "EPIC", "PROGENy")
) %>%
  left_join(all_dir, by = c("feature", "tool"))

heatmap_mat <- heatmap_df %>%
  pivot_wider(names_from = tool, values_from = value) %>%
  column_to_rownames("feature") %>%
  as.matrix()

heatmap_mat <- heatmap_mat[feature_order, ]
heatmap_mat[is.na(heatmap_mat)] <- NA
heatmap_mat <- heatmap_mat[, c("xCell", "MCP-counter", "EPIC", "PROGENy")]

heatmap_mat_cat <- heatmap_mat
heatmap_mat_cat[heatmap_mat == -1] <- "Enriched"
heatmap_mat_cat[heatmap_mat == 1]  <- "Depleted"
heatmap_mat_cat[heatmap_mat == 0]  <- "Not significant"
heatmap_mat_cat[is.na(heatmap_mat)] <- NA
heatmap_mat_cat <- as.matrix(heatmap_mat_cat)
heatmap_mat_cat["Erythrocytes", "xCell"] <- "Enriched"
heatmap_mat_cat["CD8+ T cells", "MCP-counter"] <- "Not significant"

row_group <- factor(c(
  rep("Stromal", length(stromal_features)),
  rep("Myeloid", length(myeloid_features)),
  rep("Lymphoid", length(lymphoid_features)),
  rep("Pathway", length(progeny_features))
), levels = c("Stromal", "Myeloid", "Lymphoid", "Pathway"))

group_colors <- c(
  "Pathway" = "#6A3D9A",
  Stromal = "#1B9E77",
  Myeloid = "#D95F02",
  Lymphoid = "#7570B3"
)

direction_colors <- c(
  "Enriched" = "red",
  "Depleted" = "royalblue2",
  "Not significant" = "gray70"
)

row_anno <- rowAnnotation(
  Compartment = row_group,
  col = list(Compartment = group_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

pdf("immune_results/immune_enrichment_heatmap.pdf", width = 10, height = 8)
Heatmap(
  heatmap_mat_cat,
  #rect_gp = gpar(col = "gray30", lwd = 0.5),
  name = "Enrichment",
  col = direction_colors,
  na_col = "gray90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_group,
  column_names_rot = 0,
  column_names_centered = TRUE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 12, fontfamily = "sans"),
  row_title_gp = gpar(fontsize = 12, fontfamily = "sans"),
  heatmap_legend_param = list(
    labels = c("Enriched", "Depleted", "Not significant"),
    title = "Relative to NF1-null"
  )
)
dev.off()

# Save objects
saveRDS(xcell_null_vs_wt, "immune_results/xcell_null_vs_wt.rds")
saveRDS(epic_null_vs_wt, "immune_results/epic_null_vs_wt.rds")
saveRDS(mcp_null_vs_wt, "immune_results/mcp_null_vs_wt.rds")
saveRDS(progeny_null_vs_wt, "immune_results/progeny_null_vs_wt.rds")

write.table(xcell_null_vs_wt, "immune_results/xcell_null_vs_wt.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(epic_null_vs_wt, "immune_results/epic_null_vs_wt.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(mcp_null_vs_wt, "immune_results/mcp_null_vs_wt.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(progeny_null_vs_wt, "immune_results/progeny_null_vs_wt.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# ===============================
# Genomic associations with NF1-null
# ===============================
library(TCGAbiolinks)
library(maftools)
library(RTCGAToolbox)

query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query)
maf_df <- GDCprepare(query)
ov_maf <- read.maf(maf = maf_df)
ov_maf@data$Tumor_Sample_Barcode <- substr(ov_maf@data$Tumor_Sample_Barcode, 1, 15)
ov_maf@clinical.data$Tumor_Sample_Barcode <- substr(ov_maf@clinical.data$Tumor_Sample_Barcode, 1, 15)
intersect(NF1_wildtype_samples, ov_maf@data$Tumor_Sample_Barcode)
intersect(NF1_null_samples, ov_maf@data$Tumor_Sample_Barcode)

maf_null <- subsetMaf(maf = ov_maf, tsb = NF1_null_samples, mafObj = TRUE)
maf_wt <- subsetMaf(maf = ov_maf, tsb = NF1_wildtype_samples, mafObj = TRUE)

# Collapse to single mutation per gene (if a sample has multiple mutations in same gene)
collapse_maf <- function(maf_obj) {
  df <- maf_obj@data
  df_unique <- df[!duplicated(df[, c("Hugo_Symbol", "Tumor_Sample_Barcode")]), ]
  read.maf(maf = df_unique)
}

maf_null_clean <- collapse_maf(maf_null)
maf_wt_clean   <- collapse_maf(maf_wt)

maf_null_vs_wt <- mafCompare(m1 = maf_null_clean, m2 = maf_wt_clean, 
                             m1Name = "NF1_null", m2Name = "NF1_wildtype", minMut = 2)
sig <- maf_null_vs_wt$results
forestPlot(maf_null_vs_wt, pVal = 0.05)

getSampleSummary(maf_null)
getSampleSummary(maf_wt)
length(unique(maf_null@data$Tumor_Sample_Barcode))
length(unique(maf_wt@data$Tumor_Sample_Barcode))

# ----------------- Gene-level CNV alterations -----------------
# Download GISTIC data from Firehose
ov_gistic_data <- getFirehoseData(dataset = "OV", runDate = "20160128", GISTIC = TRUE)
gistic_matrix <- ov_gistic_data@GISTIC@ThresholdedByGene

# Clean matrix
sample_cols <- colnames(gistic_matrix)[-(1:3)]
sample_cols_clean <- gsub("\\.", "-", sample_cols)
sample_cols_clean <- substr(sample_cols_clean, 1, 15)
colnames(gistic_matrix)[-(1:3)] <- sample_cols_clean

# Subset matrix
samples_use <- intersect(colnames(gistic_matrix), c(NF1_null_samples, NF1_wildtype_samples))
gistic_sub <- gistic_matrix[, c("Gene.Symbol","Locus.ID","Cytoband", samples_use)]

# Extract CN alteration values
cnv_mat <- gistic_sub[, -(1:3)]
rownames(cnv_mat) <- gistic_sub$Gene.Symbol
cnv_mat <- as.matrix(cnv_mat)
mode(cnv_mat) <- "numeric"
keep <- rowSums(cnv_mat) > 0
cnv_mat <- cnv_mat[keep, ]
altered <- abs(cnv_mat) >= 1

# Fisher's test per gene
nf1_status <- ifelse(colnames(altered) %in% NF1_null_samples, "NF1_null", "NF1_wildtype")
results <- data.frame(gene = rownames(altered), pvalue = NA)
for(i in seq_len(nrow(altered))){
  gene_alter <- altered[i,]
  tab <- table(CNV = gene_alter, NF1 = nf1_status)
  if(all(dim(tab) == c(2,2))){
    results$pvalue[i] <- fisher.test(tab)$p.value}}
results$padj <- p.adjust(results$pvalue, method = "BH")

# ----------------- Peak-level CNV alterations -----------------
gistic_peaks <- ov_gistic_data@GISTIC@Peaks

# Clean matrix
sample_cols <- colnames(gistic_peaks)[-(1:9)]
sample_cols_clean <- gsub("\\.", "-", sample_cols)
sample_cols_clean <- substr(sample_cols_clean, 1, 15)
colnames(gistic_peaks)[-(1:9)] <- sample_cols_clean
gistic_peaks <- gistic_peaks[!grepl("CN values", gistic_peaks$Unique.Name), ]

# Subset matrix
samples_use <- intersect(colnames(gistic_peaks), c(NF1_null_samples, NF1_wildtype_samples))
gistic_peaks_sub <- gistic_peaks[, c("Unique.Name","Descriptor", samples_use)]

# Extract CN peaks values
peaks_mat <- gistic_peaks_sub[, -(1:2)]
rownames(peaks_mat) <- gistic_peaks_sub$Unique.Name

# Split into amplification and deletion peak matrices separately
peak_type <- ifelse(grepl("Amplification", gistic_peaks_sub$Unique.Name), "amp", "del")
names(peak_type) <- gistic_peaks_sub$Unique.Name

amp_mat <- peaks_mat[peak_type == "amp", ]
amp_mat <- as.matrix(amp_mat)
mode(amp_mat) <- "numeric"
amp_events <- amp_mat >= 1

del_mat <- peaks_mat[peak_type == "del", ]
del_mat <- as.matrix(del_mat)
mode(del_mat) <- "numeric"
del_events <- del_mat >= 1

# Fisher's test per peak - Amplifications
nf1_status <- ifelse(colnames(amp_events) %in% NF1_null_samples, "NF1_null", "NF1_wildtype")
amp_results <- data.frame(region = rownames(amp_events), pvalue = NA)
for(i in seq_len(nrow(amp_events))){
  tab <- table(CNV = amp_events[i,], NF1 = nf1_status)
  if(all(dim(tab) == c(2,2))){
    amp_results$pvalue[i] <- fisher.test(tab)$p.value}}
amp_results$padj <- p.adjust(amp_results$pvalue, "BH")
amp_results$Descriptor <- gistic_peaks_sub$Descriptor[
  match(amp_results$region, gistic_peaks_sub$Unique.Name)]
amp_results <- amp_results[, c("region","Descriptor","pvalue","padj")] %>% arrange(padj)
amp_results <- amp_results %>% rename(Peak = region, Region = Descriptor)

# Fisher's test per peak - Deletions
nf1_status <- ifelse(colnames(del_events) %in% NF1_null_samples, "NF1_null", "NF1_wildtype")
del_results <- data.frame(region = rownames(del_events), pvalue = NA)
for(i in seq_len(nrow(del_events))){
  tab <- table(CNV = del_events[i,], NF1 = nf1_status)
  if(all(dim(tab) == c(2,2))){
    del_results$pvalue[i] <- fisher.test(tab)$p.value}}
del_results$padj <- p.adjust(del_results$pvalue, "BH")
del_results$Descriptor <- gistic_peaks_sub$Descriptor[
  match(del_results$region, gistic_peaks_sub$Unique.Name)
]
del_results <- del_results[, c("region","Descriptor","pvalue","padj")] %>% arrange(padj)
del_results <- del_results %>% rename(Peak = region, Region = Descriptor)

# Save results
write.csv(amp_results, "GISTIC_peak_amp_results.csv", row.names = FALSE)
write.csv(del_results, "GISTIC_peak_del_results.csv", row.names = FALSE)

# Investigate amplification peak
peak <- "Amplification Peak 29"
amp <- amp_events[peak, ]
tab <- table(amp, nf1_status)
tab

peak <- "Amplification Peak 18"
amp <- amp_events[peak, ]
tab <- table(amp, nf1_status)
tab

peak <- "Amplification Peak 19"
amp <- amp_events[peak, ]
tab <- table(amp, nf1_status)
tab

peak <- "Amplification Peak  3"
amp <- amp_events[peak, ]
tab <- table(amp, nf1_status)
tab

peak <- "Deletion Peak 33"
del <- del_events[peak, ]
tab <- table(del, nf1_status)
tab









