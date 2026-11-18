library(data.table)
library(dplyr)
library(jsonlite)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(tidyr)
setwd("E:/NF1/TCGA/")

# ===================================
# Gene coordinates
# ===================================
cn_genes <- list(
  NF1 = list(chr = "chr17", start = 31094977, end = 31377675),
  TP53 = list(chr = "chr17", start = 7668421, end = 7687490),
  BRCA1 = list(chr = "chr17", start = 43044295, end = 43170327)
)
other_genes <- c("BRCA1", "BRCA2", "PTEN")

# ===================================
# Load mutation data
# ===================================
mut <- fread("TCGA-OV.somaticmutation_wxs.tsv")
mut$sample_id <- substr(mut$sample, 1, 15)
samples_with_mut_data <- unique(mut$sample_id)
genes_of_interest <- c("NF1", "TP53", other_genes)
mut <- mut[gene %in% genes_of_interest]

# ===================================
# Load ASCAT3 metadata
# ===================================
metadata <- fromJSON("ascat3_metadata.json")
file_map <- metadata %>% dplyr::select(file_name, submitter_id)
file_map$sample_id <- substr(file_map$submitter_id, 1, 15)

# ===================================
# Function to classify CN state
# ===================================
classify_cn <- function(total_cn, minor_cn, ploidy) {
  if (is.na(ploidy) || length(ploidy) == 0) {
    ploidy <- 2
  }
  cn_ratio <- total_cn / ploidy
  if (total_cn == 0) {
    return("hom_del")
  } else if (total_cn == 1 & minor_cn == 0) {
    return("het_del")
  } else if (minor_cn == 0 & cn_ratio < 0.75) {
    return("clloh")
  } else if (minor_cn == 0 & cn_ratio >= 0.75 & cn_ratio <= 1.25) {
    return("cnloh")
  } else if (minor_cn == 0 & cn_ratio > 1.25) {
    return("cgloh")
  } else {
    return("wt")
  }
}

# ===================================
# Read ASCAT3 files and extract NF1 + TP53 states
# ===================================
ascn_files <- list.files(
  "ascat3_allele_specific_data/",
  pattern = "^TCGA",
  recursive = TRUE,
  full.names = TRUE
)

gene_states_list <- list()

ploidy_df <- read.delim("purity_ploidy_data.txt", sep = "\t", stringsAsFactors = FALSE)
ploidy_df$patient_id <- ploidy_df$array

for (file in ascn_files) {
  file_name <- basename(file)
  sample_id <- file_map$sample_id[
    file_map$file_name == file_name
  ]
  df <- fread(file)
  sample_states <- list(sample_id = sample_id)
  for (gene in names(cn_genes)) {
    coords <- cn_genes[[gene]]
    # Get largest overlapping segment
    seg <- df[
      Chromosome == coords$chr &
      Start <= coords$end &
      End >= coords$start
    ]
    if (nrow(seg) > 0) {
      # Compute overlap size
      seg$overlap <- pmin(seg$End, coords$end) - pmax(seg$Start, coords$start)
      # choose segment with most overlap
      seg_best <- seg[which.max(seg$overlap), ]
      total_cn <- seg_best$Copy_Number
      minor_cn <- seg_best$Minor_Copy_Number
      ploidy <- ploidy_df$ploidy[ploidy_df$patient_id == sample_id]
      state <- classify_cn(total_cn, minor_cn, ploidy)
    } else {
      state <- "unknown"
    }
    
    sample_states[[paste0(gene, "_cn")]] <- state
  }
  gene_states_list[[sample_id]] <- sample_states
}

gene_states <- bind_rows(gene_states_list)

# ===================================
# Build Oncoprint matrix for both CN and mutation data
# ===================================
# Remove samples with missing mutation data
common_samples <- intersect(
  gene_states$sample_id,
  mut$sample_id
)

# Build CN matrix
gene_states <- gene_states %>%
  filter(sample_id %in% common_samples)
gene_states$NF1_cn <- factor(
  gene_states$NF1_cn,
  levels = c("wt", "cnloh", "cgloh", "clloh", "het_del", "hom_del"))
gene_states <- gene_states %>%
  arrange(NF1_cn)
column_order <- gene_states$sample_id

cn_df <- gene_states %>%
  select(sample_id, NF1_cn, TP53_cn, BRCA1_cn)
cn_mat <- as.matrix(
  tibble::column_to_rownames(
    cn_df,
    "sample_id"))
cn_mat <- t(cn_mat)

cn_mat[cn_mat == "hom_del"] <- "Hom Del"
cn_mat[cn_mat == "het_del"] <- "Het Del"
cn_mat[cn_mat == "clloh"] <- "Copy-Loss LOH"
cn_mat[cn_mat == "cnloh"] <- "Copy-Neutral LOH"
cn_mat[cn_mat == "cgloh"] <- "Copy-Gain LOH"
cn_mat[cn_mat == "wt"] <- "Wild-type"
rownames(cn_mat) <- gsub("_cn$", "", rownames(cn_mat))

col_list <- c(
  "Wild-type" = "blue",
  "Copy-Neutral LOH" = "forestgreen",
  "Copy-Gain LOH" = "green",
  "Copy-Loss LOH" = "yellow",
  "Het Del" = "orange",
  "Hom Del" = "red"
)

# Build mutation matrix
mut <- mut %>% mutate(effect = sub(";.*", "", effect))   # Keep only first annotation
mut <- mut %>%
  mutate(
    effect = case_when(
      effect == "frameshift_variant" ~ "Frameshift",
      
      effect %in% c(
        "inframe_deletion",
        "inframe_insertion"
      ) ~ "Inframe indel",
      
      effect == "missense_variant" ~ "Missense",
      
      effect %in% c(
        "splice_acceptor_variant",
        "splice_donor_variant"
      ) ~ "Splice region",
      
      effect == "stop_gained" ~ "Nonsense",
      
      TRUE ~ NA_character_))

sort(unique(mut$effect))

mut_types <- mut %>%
  filter(!is.na(effect)) %>%
  select(sample_id, gene, effect) %>%
  distinct() %>%
  group_by(sample_id, gene) %>%
  summarise(
    effect = first(effect),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = gene,
    values_from = effect,
    values_fill = ""
  )
mut_types <- tibble(sample_id = column_order) %>%
  left_join(mut_types, by = "sample_id")
mut_types <- mut_types %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, column_order))
mut_types[mut_types == "NA"] <- ""

mut_df <- mut_types %>%
  select(sample_id, TP53, NF1, BRCA1, BRCA2, PTEN)

mut_mat <- as.matrix(
  tibble::column_to_rownames(
    mut_df,
    "sample_id"))

mut_mat <- t(mut_mat)

mut_col <- c(
  "Nonsense" = "#D62728",
  "Frameshift" = "#FF7F0E",
  "Splice region" = "#9467BD",
  "Inframe indel" = "#2CA02C",
  "Missense" = "#4C72B0"
)

# CN Oncoprint colors
cn_alter_fun <- list(
  background = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill="grey95", col=NA)),
  
  `Wild-type` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Wild-type"], col=NA)),
  
  `Copy-Neutral LOH` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Copy-Neutral LOH"], col=NA)),
  
  `Copy-Gain LOH` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Copy-Gain LOH"], col=NA)),
  
  `Copy-Loss LOH` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Copy-Loss LOH"], col=NA)),
  
  `Het Del` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Het Del"], col=NA)),
  
  `Hom Del` = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill=col_list["Hom Del"], col=NA))
)

# Mutation Oncoprint colors
mut_alter_fun <- list(
  background = function(x,y,w,h)
    grid.rect(x,y,w,h*0.9,
              gp=gpar(fill="grey95", col=NA))
)

for(i in seq_along(mut_col)) {
  
  this_col <- mut_col[i]
  
  mut_alter_fun[[names(mut_col)[i]]] <-
    eval(
      substitute(
        function(x,y,w,h)
          grid.rect(
            x,y,w,h*0.9,
            gp=gpar(fill=COL, col=NA)
          ),
        list(COL = this_col)
      )
    )
}

# Create Oncoprints
cn_onco <- oncoPrint(
  cn_mat,
  alter_fun = cn_alter_fun,
  col = col_list,
  column_order = column_order,
  
  row_title = "Copy Number",
  row_names_gp = gpar(fontsize = 14),
  row_title_gp = gpar(fontsize = 18, fontface = "bold"),
  
  show_column_names = FALSE,
  show_pct = FALSE,
  
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  
  top_annotation = NULL
)

mut_onco <- oncoPrint(
  mut_mat,
  alter_fun = mut_alter_fun,
  col = mut_col,
  column_order = column_order,
  
  row_title = "Mutations",
  row_names_gp = gpar(fontsize = 14),
  row_title_gp = gpar(fontsize = 18, fontface = "bold"),
  
  show_column_names = FALSE,
  show_pct = FALSE,
  
  remove_empty_columns = FALSE,
  remove_empty_rows = FALSE,
  
  top_annotation = NULL
)

# Combine Oncoprints
combined_onco <- cn_onco %v% mut_onco

pdf("zygosity_clonality_LOH/oncoprint_redo.pdf", width = 12, height = 8)

lgd_cn <- Legend(
  title = "CN Status",
  at = names(col_list),
  legend_gp = gpar(fill = col_list),
  title_gp = gpar(fontsize = 14, fontface = "bold"),
  labels_gp = gpar(fontsize = 12)
)

lgd_mut <- Legend(
  title = "Variant Type",
  at = names(mut_col),
  legend_gp = gpar(fill = mut_col),
  title_gp = gpar(fontsize = 14, fontface = "bold"),
  labels_gp = gpar(fontsize = 12)
)

draw(
  combined_onco,
  ht_gap = unit(12, "mm"),
  show_heatmap_legend = FALSE,
  heatmap_legend_side = "right",
  annotation_legend_list = list(lgd_cn, lgd_mut)
)

dev.off()