library(DESeq2)
library(ggplot2)
library(ggpubr)
library(readxl)
library(mclust)
library(dplyr)
setwd('E:/NF1/TCGA/')

# ===============================
# Assess BRCA1 expression levels with respect to NF1 functional group
# ===============================
# Load filtered raw counts
raw_counts <- read.delim("TcgaGtex_ovary_gene_expected_counts.tsv", sep = "\t", stringsAsFactors = FALSE)
colnames(raw_counts) <- gsub("\\.", "-", colnames(raw_counts))

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

# Subset the count matrix
filtered_counts <- raw_counts[, selected_samples]

# Create colData
colData <- data.frame(
  row.names = selected_samples,
  NF1_status = factor(group_labels, levels = c("NF1_wildtype", "NF1_null"))
)

# Import count matrix and metadata
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ NF1_status)

# Normalize with DESeq (estimate size factors)
dds <- DESeq(dds)

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Log-transform counts
log_norm_counts <- log2(norm_counts + 1)

# NF1 gene ID
nf1_gene <- "ENSG00000196712"
matching_nf1 <- grep(paste0("^", nf1_gene, "(\\.|$)"), rownames(log_norm_counts), value = TRUE)
NF1_expression <- log_norm_counts[matching_nf1[1], ]
centered_NF1 <- NF1_expression - mean(NF1_expression)

# Normalize the expression of BRCA1 around 0
gene_name <- 'ENSG00000012048'
# Find all rows that start with this gene ID (before the ".")
matching_rows <- grep(paste0("^", gene_name, "(\\.|$)"), rownames(log_norm_counts), value = TRUE)
BRCA1_expression <- log_norm_counts[matching_rows[1], ]
centered_BRCA1 <- BRCA1_expression - mean(BRCA1_expression)

combined_df <- data.frame(
  sample = colnames(log_norm_counts),
  BRCA1 = as.numeric(centered_BRCA1),
  NF1 = as.numeric(centered_NF1),
  NF1_status = colData$NF1_status
)

long_df <- combined_df %>%
  pivot_longer(
    cols = c(BRCA1, NF1),
    names_to = "gene",
    values_to = "expression"
  )

pvals <- long_df %>%
  group_by(gene) %>%
  summarise(p = wilcox.test(expression ~ NF1_status)$p.value)

p <- ggplot(long_df, aes(x = gene, y = expression, fill = NF1_status)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2), size = 1, alpha = 0.5) +
  geom_text(
    data = pvals,
    aes(x = gene, y = max(long_df$expression) + 0.5, label = paste0("p = ", signif(p, 3))),
    inherit.aes = FALSE
  ) +
  labs(
    x = "Gene",
    y = "Normalized Expression"
  ) +
  theme_classic()
ggsave("nf1_brca1_expression.png", p)

# ===============================
# Assess tumor purity with respect to NF1 functional group
# ===============================
purity_df <- read.delim("purity_ploidy_data.txt", sep = "\t", stringsAsFactors = FALSE)
purity_df$patient_id <- purity_df$array
colData$patient_id <- rownames(colData)
merged_df <- colData %>%
  mutate(sample = rownames(colData)) %>%
  left_join(purity_df, by = "patient_id")

p <- ggplot(merged_df, aes(x = NF1_status, y = purity, fill = NF1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    x = "NF1 Status",
    y = "Tumor Purity"
  ) +
  theme_classic() +
  theme(legend.position = "none")
ggsave("tumor_purity_vs_nf1_expression.png", p)

# ===============================
# Make NF1 expression box plot
# ===============================
# Read in RNA-seq counts
gz_connection <- gzfile('TcgaTargetGtex_gene_expected_count.gz', 'rt')
data <- read.delim(gz_connection)
close(gz_connection)
data

# Convert back to raw integer counts
rownames(data) <- data[[1]]
data <- data[ , -1]
raw_counts <- 2^data - 1
raw_counts <- round(raw_counts)

# Convert column names to dashes
colnames(raw_counts) <- gsub("\\.", "-", colnames(raw_counts))

# Load ASCAT3 metadata
library(jsonlite)
metadata <- fromJSON("ascat3_metadata.json")
file_map <- metadata %>%
  dplyr::select(file_name, submitter_id)
file_map$sample_id <- substr(file_map$submitter_id, 1, 15)

cn_genes <- list(
  NF1 = list(
    chr = "chr17",
    start = 31094977,
    end = 31377675
  )
)

# Classify NF1 copy-number state
classify_cn <- function(total_cn, minor_cn, ploidy){
  
  if(is.na(ploidy) || length(ploidy)==0)
    ploidy <- 2
  
  cn_ratio <- total_cn / ploidy
  
  if(total_cn == 0){
    return("hom_del")
  } else if(total_cn == 1 & minor_cn == 0){
    return("het_del")
  } else if(minor_cn == 0 & cn_ratio < 0.75){
    return("clloh")
  } else if(minor_cn == 0 & cn_ratio <= 1.25){
    return("cnloh")
  } else if(minor_cn == 0 & cn_ratio > 1.25){
    return("cgloh")
  } else{
    return("wt")
  }
}

# Read ASCAT files
ascn_files <- list.files(
  "ascat3_allele_specific_data/",
  pattern="^TCGA",
  recursive=TRUE,
  full.names=TRUE
)

ploidy_df <- read.delim(
  "purity_ploidy_data.txt",
  sep="\t",
  stringsAsFactors=FALSE
)

ploidy_df$patient_id <- ploidy_df$array

gene_states_list <- list()

library(data.table)
for(file in ascn_files){
  
  file_name <- basename(file)
  
  sample_id <- file_map$sample_id[
    file_map$file_name==file_name
  ]
  
  df <- fread(file)
  
  coords <- cn_genes$NF1
  
  seg <- df[
    Chromosome==coords$chr &
      Start <= coords$end &
      End >= coords$start
  ]
  
  if(nrow(seg)>0){
    
    seg$overlap <-
      pmin(seg$End,coords$end) -
      pmax(seg$Start,coords$start)
    
    seg_best <- seg[which.max(seg$overlap),]
    
    total_cn <- seg_best$Copy_Number
    minor_cn <- seg_best$Minor_Copy_Number
    
    ploidy <- ploidy_df$ploidy[
      ploidy_df$patient_id==sample_id
    ]
    
    state <- classify_cn(
      total_cn,
      minor_cn,
      ploidy
    )
    
  }else{
    
    state <- "unknown"
    
  }
  
  gene_states_list[[sample_id]] <-
    data.frame(
      sample_id=sample_id,
      NF1_cn=state
    )
  
}

gene_states <- bind_rows(gene_states_list)

# Read in genotype groups of sample IDs
NF1_genotype_files <- list(
  het_del = "NF1_het_del.txt",
  hom_del = "NF1_hom_del.txt",
  frameshift = "NF1_frameshift_variant.txt",
  missense = "NF1_missense_variant.txt",
  nonsense = "NF1_stop_gained.txt",
  inframe_del = "NF1_inframe_deletion.txt",
  splice_acceptor = "NF1_splice_acceptor_variant.txt",
  synonymous = "NF1_synonymous_variant.txt",
  gtex = "gtex_ovary_sample_ids.txt"
)
NF1_genotype_lists <- lapply(NF1_genotype_files, function(file) {
  readLines(file)
})

# Match and filter columns
selected_columns <- character()
NF1_genotype_labels <- character()
already_selected <- character()

#-----------------------------
# Helper function
#-----------------------------
add_group <- function(samples, label) {
  
  matches <- setdiff(
    intersect(colnames(raw_counts), samples),
    already_selected
  )
  
  selected_columns <<- c(selected_columns, matches)
  
  NF1_genotype_labels <<- c(
    NF1_genotype_labels,
    rep(label, length(matches))
  )
  
  already_selected <<- c(already_selected, matches)
}

#-----------------------------
# Add deletion groups first
#-----------------------------
add_group(
  NF1_genotype_lists$hom_del,
  "hom_del"
)

add_group(
  NF1_genotype_lists$het_del,
  "het_del"
)

#-----------------------------
# Add mutation groups
#-----------------------------
# Missense
add_group(
  NF1_genotype_lists$missense,
  "missense_loh"
)

# Truncating
truncating_samples <- unique(c(
  NF1_genotype_lists$frameshift,
  NF1_genotype_lists$nonsense,
  NF1_genotype_lists$splice_acceptor
))

add_group(
  truncating_samples,
  "truncating_loh"
)

#-----------------------------
# Add LOH groups
#-----------------------------
for(state in c("cgloh","cnloh","clloh")){
  
  samples <- gene_states$sample_id[
    gene_states$NF1_cn == state
  ]
  
  add_group(samples, state)
}

#-----------------------------
# Finally add GTEx
#-----------------------------
add_group(
  NF1_genotype_lists$gtex,
  "gtex"
)

table(NF1_genotype_labels)

# Subset the count matrix
filtered_counts <- raw_counts[, selected_columns]

# Create colData
colData <- data.frame(row.names = selected_columns, 
                      group = factor(NF1_genotype_labels))

# Import count matrix and metadata
dds <- DESeqDataSetFromMatrix(countData = filtered_counts,
                              colData = colData,
                              design = ~ group)

# Normalize with DESeq (estimate size factors)
dds <- DESeq(dds)

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Log-transform counts
log_norm_counts <- log2(norm_counts + 1)

# Normalize the expression of NF1 around 0
gene_name <- 'ENSG00000196712'
# Find all rows that start with this gene ID (before the ".")
matching_rows <- grep(paste0("^", gene_name, "(\\.|$)"), rownames(log_norm_counts), value = TRUE)
NF1_expression <- log_norm_counts[matching_rows[1], ]
centered_NF1 <- NF1_expression - mean(NF1_expression)

# Extract expression of NF1
NF1_expression_df <- data.frame(
  expression = centered_NF1,
  group = colData$group
)

# Remove wildtype outlier
NF1_expression_df <- NF1_expression_df[row.names(NF1_expression_df) != "TCGA-25-1870-01", ]

# Order groups for plotting
NF1_expression_df$group <- recode(
  NF1_expression_df$group,
  "gtex" = "GTEx",
  
  "cgloh" = "Copy-Gain LOH",
  "cnloh" = "Copy-Neutral LOH",
  "clloh" = "Copy-Loss LOH",
  
  "missense_loh" = "Missense + LOH",
  "truncating_loh" = "Truncating + LOH",
  
  "het_del" = "Het. Del (TCN = 1)",
  "hom_del" = "Hom. Del"
)
NF1_expression_df$group <- factor(
  NF1_expression_df$group,
  levels = c("GTEx", "Copy-Neutral LOH", "Copy-Gain LOH", "Copy-Loss LOH", 
             "Missense + LOH", "Truncating + LOH", "Het. Del (TCN = 1)", "Hom. Del"))

# Color boxes manually
box_colors <- c(
  "GTEx" = "grey70",
  
  "Copy-Gain LOH" = "green",
  "Copy-Neutral LOH" = "forestgreen",
  "Copy-Loss LOH" = "yellow",

  "Missense + LOH" = "blue",
  "Truncating + LOH" = "purple",
  
  "Het. Del (TCN = 1)" = "orange",
  "Hom. Del" = "red"
)

# Quantile-based approach
low_cut <- quantile(NF1_expression_df$expression, 0.10, na.rm = TRUE)
high_cut <- quantile(NF1_expression_df$expression, 0.90, na.rm = TRUE)
NF1_expression_df <- NF1_expression_df %>%
  mutate(
    NF1_functional_status = case_when(
      group %in% c("hom_del") |
        expression < low_cut ~ "NF1-null",
      group %in% expression >= low_cut & expression < high_cut ~ "NF1-hypomorph",
      TRUE ~ "NF1-wildtype"
    )
  )

# Statistical comparisons of interest
comparisons <- list(
  c("Copy-Neutral LOH", "Het. Del (TCN = 1)"),
  c("Copy-Neutral LOH", "Hom. Del"),
  c("Het. Del (TCN = 1)", "Hom. Del")
)

# Count samples per group
group_counts <- table(NF1_expression_df$group)
group_labels <- paste0(names(group_counts), " (n=", group_counts, ")")
unique(NF1_expression_df$group)

# Save to file
write.table(NF1_expression_df, file = "NF1_expression_functional_groups_redo.tsv",
            sep = "\t", quote = FALSE, row.names = TRUE)

# Plot expression box plot
ggboxplot(NF1_expression_df, x = "group", y = "expression", fill = "group") +
  stat_compare_means(method = "t.test", comparisons = comparisons, label = "p.signif") +
  scale_x_discrete(labels = group_labels) +
  scale_fill_manual(values = box_colors) +
  labs(
    x = "NF1 Genotype",
    y = "Normalized NF1 mRNA expression") + 
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_line(
          color = "grey92",
          linewidth = 0.25
        ),
        panel.grid.minor = element_blank())
ggsave("nf1_expression_boxplots_redo.png", width = 10, height = 10, bg = "white")
ggsave("figures/nf1_expression_boxplots_redo.png", width = 10, height = 10, bg = "white")

# Visualize expression distribution
ggplot(NF1_expression_df, aes(x = expression, fill = group)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "NF1 expression distribution by genotype",
       x = "Normalized NF1 expression", y = "Density")