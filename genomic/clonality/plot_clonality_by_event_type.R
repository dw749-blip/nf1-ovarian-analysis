library(ggplot2)
library(dplyr)
library(readr)
setwd("E:/NF1/TCGA/NF1_LOH_clonality_ABSOLUTE/")

# ===================================
# Read counts (already computed by build_full_cohort_clonality.py /
# plot_clonality_by_event_type.py)
# ===================================
counts <- read_tsv("clonality_by_event_type_counts.tsv", show_col_types = FALSE)

order <- c("Copy-Neutral LOH", "Copy-Gain LOH", "Copy-Loss LOH",
           "Mutation + LOH", "Het. Del (TCN = 1)", "Hom. Del")
base_labels <- c("Copy-Neutral LOH", "Copy-Gain LOH", "Copy-Loss LOH",
                  "Mutation + LOH", "Het. Del. (TCN = 1)", "Hom. Del.")

counts <- counts %>%
  mutate(event_type = factor(event_type, levels = order)) %>%
  arrange(event_type) %>%
  mutate(label = paste0(base_labels, " (n=", total, ")"),
         label = factor(label, levels = label))

# ===================================
# Long format for stacked bars
# ===================================
plot_df <- counts %>%
  select(label, pct_clonal, pct_subclonal) %>%
  tidyr::pivot_longer(cols = c(pct_clonal, pct_subclonal),
                       names_to = "clonality", values_to = "pct") %>%
  mutate(clonality = recode(clonality,
                             pct_clonal = "Clonal",
                             pct_subclonal = "Subclonal"),
         clonality = factor(clonality, levels = c("Subclonal", "Clonal")))

# ===================================
# Plot
# ===================================
p <- ggplot(plot_df, aes(x = label, y = pct, fill = clonality)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("Clonal" = "#3182bd", "Subclonal" = "#f0b429"),
                     breaks = c("Clonal", "Subclonal")) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "NF1 Genotype", y = "% of samples", fill = NULL) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p
ggsave("clonality_by_event_type_R.png", plot = p, width = 7.5, height = 6, dpi = 200)
