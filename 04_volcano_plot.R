
# Lesson 6: Volcano Plot

library(ggplot2)

# Load differential expression results
results_df <- read.csv("03_differential_expression_results.csv")

# Remove genes without adjusted p-values
results_df <- results_df[!is.na(results_df$padj), ]

# Create categories for plotting
results_df$category <- "Not significant"

results_df$category[
  results_df$padj < 0.05 &
  results_df$log2FoldChange >= 1
] <- "Upregulated"

results_df$category[
  results_df$padj < 0.05 &
  results_df$log2FoldChange <= -1
] <- "Downregulated"

# Avoid infinite values when calculating -log10(padj)
results_df$padj_for_plot <- pmax(
  results_df$padj,
  .Machine$double.xmin
)

results_df$minus_log10_padj <- -log10(
  results_df$padj_for_plot
)

# Draw the volcano plot
volcano_plot <- ggplot(
  results_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_padj,
    color = category
  )
) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Expression: Airway RNA-seq",
    x = "Log2 Fold Change (treated vs untreated)",
    y = "-Log10 Adjusted P-value",
    color = "Gene category"
  ) +
  theme_bw()

print(volcano_plot)

ggsave(
  "04_volcano_plot.pdf",
  plot = volcano_plot,
  width = 8,
  height = 6
)

cat("Volcano plot saved successfully.\n")

