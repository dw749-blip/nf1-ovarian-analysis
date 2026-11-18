library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(survival)
library(survminer)
setwd('E:/NF1/TCGA/')

clinical_data <- read.delim('TCGA-OV.clinical.tsv', header = TRUE, sep = '\t')
clinical_data <- clinical_data %>%
  mutate(patient = substr(sample, 1, 12))
nf1_df <- read.delim("NF1_expression_functional_groups.tsv", sep = "\t", stringsAsFactors = FALSE)
nf1_df <- nf1_df %>%
  tibble::rownames_to_column(var = "sample") %>%
  mutate(patient = substr(sample, 1, 12))

# Download clinical data
query <- GDCquery(project = "TCGA-OV",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR XML")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")
clinical <- clinical %>%
  mutate(OS_time  = ifelse(!is.na(days_to_death), days_to_death, days_to_last_followup),
         OS_event = ifelse(!is.na(days_to_death), 1, 0))
clinical <- clinical %>%
  rename(patient = bcr_patient_barcode,
         age = age_at_initial_pathologic_diagnosis,
         race = race_list)

# Merge with tumor stage column from other file
clinical_data_subset <- clinical_data %>%
  dplyr::select(patient, figo_stage.diagnoses)
clinical_merged <- clinical %>%
  left_join(clinical_data_subset, by = "patient")
clinical_merged <- clinical_merged %>%
  rename(tumor_stage = figo_stage.diagnoses)

# Merge clinical and NF1 data
merged <- inner_join(clinical_merged, nf1_df, by = "patient")

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

ggplot(long_df, aes(x = gene, y = expression, fill = NF1_status)) +
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

# ===============================
# Assess tumor purity with respect to NF1 functional group
# ===============================
purity_df <- read.delim("purity_ploidy_data.txt", sep = "\t", stringsAsFactors = FALSE)
purity_df$patient_id <- purity_df$array
colData$patient_id <- rownames(colData)
merged_df <- colData %>%
  mutate(sample = rownames(colData)) %>%
  left_join(purity_df, by = "patient_id")

ggplot(merged_df, aes(x = NF1_status, y = purity, fill = NF1_status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    x = "NF1 Status",
    y = "Tumor Purity"
  ) +
  theme_classic() +
  theme(legend.position = "none")

# ===============================
# All NF1-null vs. WT
# ===============================
# All samples with NF1 genotype null (hom. del, short mutation)
short_mutation_ids <- read.delim("NF1_short_mutation_IDs.txt", header = TRUE, stringsAsFactors = FALSE)
subset_data <- clinical_merged %>%
  left_join(nf1_df, by = "patient") %>%
  mutate(
    NF1_expression = case_when(
      group != "GTEx" & NF1_functional_status == "NF1-wildtype" ~ "High",
      (patient %in% short_mutation_ids$submitter_id) | (group == "Hom. Del") ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(NF1_expression)) %>%
  mutate(
    NF1_expression = factor(NF1_expression, levels = c("High", "Low"))
  )

subset_data <- merged %>%
  mutate(
    NF1_expression = case_when(
      group != "GTEx" & NF1_functional_status == "NF1-wildtype" ~ "High",
      NF1_functional_status == "NF1-null" ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(NF1_expression)) %>%
  mutate(
    NF1_expression = factor(NF1_expression, levels = c("High", "Low"))
  )

# Add BRCA1 expression and tumor purity as covariates
brca1_df <- data.frame(
  sample = names(BRCA1_expression),
  BRCA1_expression = as.numeric(BRCA1_expression))
subset_data <- subset_data %>%
  left_join(brca1_df, by = c("sample" = "sample"))
purity_df_clean <- purity_df %>%
  select(array, purity)
subset_data <- subset_data %>%
  left_join(purity_df_clean, by = c("sample" = "array"))

# Simplify tumor stage
subset_data <- subset_data %>%
  mutate(
    tumor_stage_raw = tumor_stage,
    tumor_stage_simple = case_when(
      grepl("iv", tumor_stage_raw, ignore.case = TRUE)  ~ "IV",
      grepl("iii", tumor_stage_raw, ignore.case = TRUE) ~ "III",
      grepl("ii", tumor_stage_raw, ignore.case = TRUE)  ~ "II",
      grepl("i", tumor_stage_raw, ignore.case = TRUE)   ~ "I",
      TRUE ~ NA_character_
    ),
    tumor_stage_simple = factor(tumor_stage_simple, levels = c("I","II","III","IV"))
  )

# Clean final df
model_df <- subset_data %>%
  mutate(
    age = as.numeric(age),
    OS_time = as.numeric(OS_time),
    OS_event = as.numeric(OS_event),
    BRCA1_expression = as.numeric(BRCA1_expression),
    purity = as.numeric(purity)
  ) %>%
  filter(!is.na(OS_time) & !is.na(OS_event) & !is.na(age) & !is.na(tumor_stage_simple) & !is.na(BRCA1_expression) & !is.na(purity))
# Make "Low" expression the reference for more interpretable hazard ratio
model_df$NF1_expression <- relevel(model_df$NF1_expression, ref = "Low")

table(model_df$NF1_functional_status)
table(model_df$NF1_expression)
table(model_df$tumor_stage_simple)
table(model_df$BRCA1_expression)
table(model_df$purity)
model_df[, c("group", "expression", "NF1_functional_status", "NF1_expression", "OS_time", "BRCA1_expression", "purity")]
print(colnames(model_df))
levels(model_df$NF1_functional_status)
levels(model_df$NF1_expression)

# Run Cox proportional hazards model
cox_model <- coxph(Surv(OS_time, OS_event) ~ NF1_expression + age + tumor_stage_simple + BRCA1_expression + purity,
                   data = model_df)
summary(cox_model)

# Plot adjusted survival curves directly from the Cox model
model_df$NF1_expression <- factor(model_df$NF1_expression,
                                  levels = c("Low", "High"))
surv_df <- data.frame(
  NF1_expression = factor(c("Low", "High"),
                          levels = levels(model_df$NF1_expression)),
  age = rep(mean(model_df$age, na.rm = TRUE), 2),
  tumor_stage_simple = factor("III", levels = levels(model_df$tumor_stage_simple)),
  BRCA1_expression = rep(mean(model_df$BRCA1_expression, na.rm = TRUE), 2),
  purity = rep(mean(model_df$purity, na.rm = TRUE), 2)
)
surv_df
fit <- survfit(cox_model, newdata = surv_df)
p_adj <- ggsurvplot(
  fit,
  data = model_df,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  palette = c("#0072B2", "#D55E00"),
  legend.title = "NF1 Expression",
  legend.labs = c("Low", "High"),
  ggtheme = theme_minimal(),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  size = 0.8
)

p_adj

# Annotate hazard ratios and p-values
hr <- summary(cox_model)$coefficients
hr_text <- paste0(
  levels(model_df$group2)[2], "HR = ",
  round(hr[1, "exp(coef)"], 2),
  ", p = ",
  signif(hr[1, "Pr(>|z|)"], 3)
)

# Add annotation to plot
p <- p_adj$plot
p <- p +
  annotate("text",
           x = Inf, y = Inf,
           label = hr_text,
           hjust = 1.1, vjust = 1.5,
           size = 4, fontface = "italic")
p_adj$plot <- p
p_adj
ggsave("survival_curve_cox_adjusted_expression_alone_vs_WT.png", p_adj)
ggsave("figures/survival_curve_cox_adjusted_expression_alone_vs_WT.png", p_adj)

# Plot KM survival curve
fit <- survfit(Surv(OS_time, OS_event) ~ NF1_expression, data = model_df)
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
p <- ggsurvplot(
  fit,
  data = model_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#0072B2", "#D55E00"),
  xlab = "Time (days)",
  ylab = "Survival Probability",
  legend.title = "NF1 expression",
  legend.labs = c("Low", "High") # Low is first because model_df$NF1_expression was releveled so "Low" is the reference
)
p
ggsave("survival_curve_expression_alone_vs_WT.png", p)
ggsave("figures/survival_curve_expression_alone_vs_WT.png", p)