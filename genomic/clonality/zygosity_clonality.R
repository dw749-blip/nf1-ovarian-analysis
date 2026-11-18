library(ggplot2)
library(dplyr)
library(readr)
setwd("E:/NF1/TCGA/zygosity_clonality_LOH/")

# ===================================
# Compute allele ratio
# ===================================
df <- read_csv("zygosity_clonality_data.csv")
df <- na.omit(df)

df <- df %>%
  mutate(
    allele_ratio = mutant_cn / major_cn,
    allele_ratio = ifelse(allele_ratio > 1, 1, allele_ratio),
    major_cn = as.factor(major_cn),
    mutation_type = as.factor(mutation_type)
  )

# ===================================
# Define colors and shapes
# ===================================
# Colors for mutation type
color_values <- c(
  "nonsense" = "#d62728",
  "frameshift" = "purple",
  "missense" = "blue",
  "splice acceptor" = "orange",
  "inframe del" = "pink",
  "synonymous" = "#2ca02c"
)

# Shapes for major_cn value
shape_values <- c(
  "1" = 16,  # circle
  "2" = 17,  # triangle
  "3" = 15   # square
)

p <- ggplot(df, aes(x = allele_ratio,
                    y = cellular_prevalence,
                    color = mutation_type,
                    shape = major_cn)) +
  geom_jitter(width = 0, height = 0, size = 3, alpha = 0.85) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    name = "Allele Ratio (mutant CN / major CN)"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    name = "Clonality (CCF)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey95", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.2)
  ) +
  labs(
    color = "Mutation Type",
    shape = "Major Copy Number"
  )

p
ggsave("zygosity_clonality_plot.pdf", plot = p, width = 9, height = 6)

