library(DESeq2)
library(dplyr)
library(fgsea)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(readr)
library(tidyverse)
library(stats)
library(RRHO2)
library(RRHO)
library(EnhancedVolcano)
library(variancePartition)
library(BiocParallel)
setwd('E:/NF1/TCGA/')

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

# ===============================
# Subset count matrix
# ===============================
filtered_counts <- raw_counts[, selected_samples]

# Remove genes with no counts in any samples
filtered_counts <- filtered_counts[rowSums(filtered_counts) > 0, ]

colData <- data.frame(
  row.names = selected_samples,
  NF1_status = factor(group_labels, levels = c("NF1_wildtype", "NF1_null"))
)

# ===============================
# Add covariates to DESeq2 model
# ===============================
colData$sample_id <- rownames(colData)
colData$patient <- substr(colData$sample_id, 1, 12)

# Sequencing batch/center, age, tumor stage
clinical_data <- read.delim('TCGA-OV.clinical.tsv', header = TRUE, sep = '\t')
clinical_data <- clinical_data %>%
  mutate(patient = substr(sample, 1, 12))
clinical_subset <- clinical_data[, c("patient", "name.tissue_source_site")]
colnames(clinical_subset) <- c("patient", "batch_center")
colData <- merge(colData, clinical_subset, by = "patient", all.x = TRUE)

# Organize colData
colData <- colData[!duplicated(colData$sample_id), ]
colData <- colData[match(colnames(filtered_counts), colData$sample_id), ]
rownames(colData) <- colData$sample_id
colData$batch_center <- factor(colData$batch_center)

# Clean factor names to avoid spaces or symbols
levels(colData$batch_center) <- make.names(levels(colData$batch_center))

# ===============================
# PCA to visualize clustering of samples
# ===============================
dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = colData,
  design = ~ batch_center + NF1_status
)

# VST
vsd <- vst(dds, blind = FALSE)

# PCA
pcaData <- plotPCA(vsd, intgroup = "NF1_status", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot using ggplot2
ggplot(pcaData, aes(PC1, PC2, color = NF1_status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()

# ===============================
# Run DESeq2
# ===============================
dds <- DESeq(dds)

res_list <- list(
  "NF1_null_vs_wildtype" = results(dds, contrast = c("NF1_status", "NF1_null", "NF1_wildtype"))
)

# Save results
names(res_list)
for (comparison_name in names(res_list)) {
  #comparison_name <- "NF1_het_del_vs_wildtype"
  cat("Processing comparison ", comparison_name, "\n")
  flush.console()
  res <- as.data.frame(res_list[[comparison_name]])
  res <- res[order(res$padj, decreasing = FALSE), ]
  write.table(res, paste0("rna_seq_analysis_results/DESeq2_", comparison_name, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
}

# Add converted gene symbols to DESeq2 results
files <- c(
  "rna_seq_analysis_results/DESeq2_NF1_null_vs_wildtype.tsv"
)

# Get valid ENSEMBL keys
valid_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
for (file in files) {
  df <- read.delim(file, row.names = 1)
  df <- df %>% rownames_to_column(var = "ENSEMBL")
  df <- df %>% mutate(ENSEMBL_clean = sub("\\..*$", "", ENSEMBL))
  valid_ids <- intersect(unique(df$ENSEMBL_clean), valid_keys)
  if (length(valid_ids) > 0) {
    gene_map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = valid_ids,
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    ) %>% distinct(ENSEMBL, .keep_all = TRUE)
    df <- df %>% left_join(gene_map, by = c("ENSEMBL_clean" = "ENSEMBL"))
  } else {
    df$SYMBOL <- NA
  }
  df <- df %>% dplyr::select(ENSEMBL, SYMBOL, everything(), -ENSEMBL_clean)
  # Save
  new_name <- sub(".tsv$", "_annotated.tsv", file)
  write.table(
    df,
    sub(".tsv$", "_annotated.tsv", file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
  
# ===============================
# GSEA
# ===============================
setwd('E:/')
pathways_H  <- gmtPathways("mh.all.v2024.1.Mm.symbols.gmt")
pathways_C2 <- gmtPathways("m2.all.v2024.1.Mm.symbols.gmt")
pathways_C5 <- gmtPathways("m5.all.v2024.1.Mm.symbols.gmt")
all_pathways <- c(pathways_H, pathways_C2, pathways_C5)
writeGmtPathways(all_pathways, "mm.all.v2024.1.Mm.symbols.gmt")

# Load pathway gene sets
pathways_H  <- gmtPathways("h.all.v2025.1.Hs.symbols.gmt")
pathways_C2 <- gmtPathways("c2.all.v2025.1.Hs.symbols.gmt")
pathways_C5 <- gmtPathways("c5.all.v2025.1.Hs.symbols.gmt")
all_pathways <- c(pathways_H, pathways_C2, pathways_C5)

# Loop through each pairwise comparison
for (comparison_name in names(res_list)) {
  #comparison_name <- "NF1_het_del_vs_wildtype"
  cat("Running GSEA for ", comparison_name, "\n")
  flush.console()
  
  # Extract the DESeq2 result
  res_df <- as.data.frame(res_list[[comparison_name]]) %>%
    tibble::rownames_to_column("gene_id")
  
  # Map Ensembl IDs to gene symbols
  res_df$gene_id <- gsub("\\.\\d+$", "", res_df$gene_id)
  symbol_map <- AnnotationDbi::select(org.Hs.eg.db,
                                      keys = res_df$gene_id,
                                      columns = c("SYMBOL"),
                                      keytype = "ENSEMBL")
  symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]
  res_df <- merge(res_df, symbol_map, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)
  
  # Remove duplicates, rank by logFC
  res_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(log2FoldChange)) %>%
    dplyr::select(SYMBOL, log2FoldChange) %>%
    group_by(SYMBOL) %>%
    summarize(log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))]) %>%
    ungroup() %>%
    arrange(desc(log2FoldChange))
  
  # Create .rnk file
  rnk_file <- paste0("rna_seq_analysis_results/", comparison_name, ".rnk")
  write_tsv(res_df, rnk_file, col_names = FALSE)
  ranks <- fread(rnk_file, header = FALSE)
  colnames(ranks) <- c("gene", "rank")
  ranks_vector <- deframe(ranks)
  
  # Run fgsea on H pathways only
  fgsea_H <- fgseaMultilevel(pathways = pathways_H, stats = ranks_vector, nproc = 1, nPermSimple = 10000)
  
  # Run fgsea on all pathways
  fgsea_all <- fgsea(pathways = all_pathways, stats = ranks_vector, nproc = 1, nperm = 1000, minSize = 15, maxSize = 500)
  
  # Save results
  output_gsea_H_file <- paste0("rna_seq_analysis_results/fgsea_hallmark_", comparison_name, ".tsv")
  write.table(subset(fgsea_H, select = -leadingEdge), output_gsea_H_file, sep = "\t", quote = FALSE, row.names = FALSE)
  output_gsea_all_file <- paste0("rna_seq_analysis_results/fgsea_all_", comparison_name, ".tsv")
  write.table(subset(fgsea_all, select = -leadingEdge), output_gsea_all_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # --- Plot Hallmark pathways ---
  ggplot(fgsea_H, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<0.05)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title = paste0("Hallmark GSEA: ", comparison_name)) +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"), name = "padj < 0.05") +
    theme_minimal(base_size = 13)
  # Save
  output_hallmark_file <- paste0("rna_seq_analysis_results/fgsea_hallmark_plot_", comparison_name, ".png")
  ggsave(output_hallmark_file, width = 10, height = 10, bg = "white")
  output_hallmark_file <- paste0("figures/fgsea_hallmark_plot_", comparison_name, ".png")
  ggsave(output_hallmark_file, width = 10, height = 10, bg = "white")
  
  # --- Plot top and bottom 25 pathways from all sets ---
  top_bottom_50 <- fgsea_all %>%
    arrange(desc(NES)) %>%
    slice_head(n = 25) %>%
    bind_rows(fgsea_all %>% arrange(NES) %>% slice_head(n = 25)) %>%
    mutate(pathway_short = str_trunc(pathway, 50))
  
  ggplot(top_bottom_50, aes(reorder(pathway_short, NES), NES)) +
    geom_col(aes(fill = padj < 0.05)) +
    coord_flip() +
    labs(title = paste0("Top 25 Up/Down Pathways: ", comparison_name),
         x = "Pathway", y = "Normalized Enrichment Score") +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey70"), name = "padj < 0.05") +
    theme_minimal(base_size = 13)
  # Save
  output_top50_file <- paste0("rna_seq_analysis_results/fgsea_top50_plot_", comparison_name, ".png")
  ggsave(output_top50_file, width = 10, height = 10, bg = "white")
  output_top50_file <- paste0("figures/fgsea_top50_plot_", comparison_name, ".png")
  ggsave(output_top50_file, width = 10, height = 10, bg = "white")
  
  # --- Prepare input for EnrichmentMap ---
  # Generic result file format
  em_generic <- fgsea_all %>%
    dplyr::transmute(
      GO.ID = pathway,
      Description = NA,
      p.Val = pval,
      FDR = padj,
      NES = NES,
      Phenotype = ifelse(NES >= 0, 1, -1)
    ) %>%
    arrange(FDR)
  
  output_em_file <- paste0("figures/GSEA/EM_input_fgsea_generic_", comparison_name, ".tsv")
  write_tsv(em_generic, output_em_file)

  # .rnk file
  rnk_df <- tibble(
    gene = names(ranks_vector),
    rank = as.numeric(ranks_vector)
  )
  output_ranks_file <- paste0("figures/GSEA/EM_input_ranks_", comparison_name, ".rnk")
  write_tsv(rnk_df, output_ranks_file, col_names = FALSE)
  
  # GMT file
  writeGmtPathways(all_pathways, "figures/GSEA/EM_input_gmt_pathways.gmt")
}


# -------------------- PROTEIN EXPRESSION ANALYSIS ----------------------
library(limma)
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

# Load normalized counts
protein_raw <- read.delim("TCGA-OV.protein.tsv", stringsAsFactors = FALSE)
colnames(protein_raw) <- gsub("\\.", "-", colnames(protein_raw))
colnames(protein_raw) <- gsub("([0-9]{2})[A-Z]$", "\\1", colnames(protein_raw))
protein_raw <- protein_raw %>% column_to_rownames("peptide_target")
protein_mat <- protein_raw[, intersect(colnames(protein_raw), selected_samples)]

colData <- data.frame(
  row.names = selected_samples,
  NF1_status = factor(group_labels, levels = c("NF1_wildtype", "NF1_null"))
)

# ===============================
# Add covariates to DESeq2 model
# ===============================
colData$sample_id <- rownames(colData)
colData$patient <- substr(colData$sample_id, 1, 12)

# Sequencing batch/center, age, tumor stage
clinical_data <- read.delim('TCGA-OV.clinical.tsv', header = TRUE, sep = '\t')
clinical_data <- clinical_data %>%
  mutate(patient = substr(sample, 1, 12))
clinical_subset <- clinical_data[, c("patient", "name.tissue_source_site")]
colnames(clinical_subset) <- c("patient", "batch_center")
colData <- merge(colData, clinical_subset, by = "patient", all.x = TRUE)

# Organize colData
colData <- colData[!duplicated(colData$sample_id), ]
colData <- colData[match(colnames(protein_mat), colData$sample_id), ]
rownames(colData) <- colData$sample_id
colData$batch_center <- factor(colData$batch_center)

# Clean factor names to avoid spaces or symbols
levels(colData$batch_center) <- make.names(levels(colData$batch_center))

colData <- colData[match(colnames(protein_mat), colData$sample_id), ]

# ===============================
# Run limma
# ===============================
design <- model.matrix(~ batch_center + NF1_status, data = colData)
fit <- lmFit(protein_mat, design)
fit <- eBayes(fit)

# Save results
res <- topTable(fit, coef = "NF1_statusNF1_null", number = Inf, sort.by = "P")
res$protein <- rownames(res)
write.table(res, "rna_seq_analysis_results/limma_NF1_null_vs_wildtype.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano plot
res <- read.delim("rna_seq_analysis_results/limma_NF1_null_vs_wildtype.tsv")

res <- res %>%
  mutate(
    negLog10FDR = -log10(adj.P.Val),
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "Upregulated",
      adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Label top proteins
top_labels <- res %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 15)

# Volcano plot
ggplot(res, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = significance),
             alpha = 0.8,
             size = 2) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "grey50") +
  geom_text_repel(
    data = top_labels,
    aes(label = protein),
    size = 4,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "firebrick",
      "Downregulated" = "#0072B2",
      "Not Significant" = "grey70"
    )
  ) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("p-adj")),
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 14)
  )

ggsave("rna_seq_analysis_results/protein_volcano_plot.pdf", width = 7, height = 7)

# ===============================
# Re-make Hallmark GSEA results plot for better figure
# ===============================
# Remove non-significant pathways for less clutter
# Shorten pathway names (no "Hallmark" or underscores) for clarity
library(stringr)
fgsea <- read.delim("rna_seq_analysis_results/fgsea_hallmark_NF1_null_vs_wildtype.tsv")

fgsea_sig <- fgsea %>%
  filter(padj < 0.05) %>%
  mutate(
    pathway_clean = pathway %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " "),
    direction = ifelse(NES > 0,
                       "Enriched in NF1-null",
                       "Depleted in NF1-null")
  ) %>%
  arrange(NES)

# Plot
ggplot(
  fgsea_sig,
  aes(
    x = reorder(pathway_clean, NES),
    y = NES,
    fill = direction
  )
) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Enriched in NF1-null" = "firebrick",
      "Depleted in NF1-null" = "#0072B2"
    )
  ) +
  labs(
    x = NULL,
    y = "Normalized Enrichment Score (NES)",
    fill = NULL
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 16),
    legend.position = "top",
    legend.text = element_text(size = 13)
  )

ggsave("rna_seq_analysis_results/fgsea_hallmark_plot_remake_NF1_null_vs_wildtype.pdf", width = 10, height = 10)