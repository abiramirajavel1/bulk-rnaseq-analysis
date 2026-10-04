
# Lesson 5: Differential Expression Analysis
# Compare dexamethasone-treated vs untreated samples

library(airway)
library(SummarizedExperiment)
library(DESeq2)

# Load dataset
data("airway", package = "airway")

# Extract raw counts and sample metadata
count_matrix <- assay(airway, "counts")
sample_metadata <- as.data.frame(colData(airway))

# Set untreated samples as the reference group
sample_metadata$dex <- factor(
  sample_metadata$dex,
  levels = c("untrt", "trt")
)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ dex
)

# Filter genes with low counts
min_group_size <- min(table(sample_metadata$dex))

keep <- rowSums(counts(dds) >= 10) >= min_group_size
dds <- dds[keep, ]

cat("Genes retained for analysis:", nrow(dds), "\n")

# Run differential expression analysis
dds <- DESeq(dds)

# Compare treated against untreated samples
res <- results(
  dds,
  contrast = c("dex", "trt", "untrt"),
  alpha = 0.05
)

# Convert results to a data frame
results_df <- as.data.frame(res)
results_df$gene_id <- rownames(results_df)

# Sort by adjusted p-value
results_df <- results_df[
  order(results_df$padj, na.last = TRUE),
]

# Identify significant genes
significant <- !is.na(results_df$padj) &
  results_df$padj < 0.05

upregulated <- significant &
  results_df$log2FoldChange > 0

downregulated <- significant &
  results_df$log2FoldChange < 0

cat("Total genes tested:", nrow(results_df), "\n")
cat("Significant genes:", sum(significant), "\n")
cat("Upregulated genes:", sum(upregulated), "\n")
cat("Downregulated genes:", sum(downregulated), "\n")

cat("\nTop 10 results:\n")
print(head(results_df, 10))

# Save all results
write.csv(
  results_df,
  "03_differential_expression_results.csv",
  row.names = FALSE
)

# Save significant genes only
write.csv(
  results_df[significant, ],
  "03_significant_genes.csv",
  row.names = FALSE
)

cat("\nAnalysis complete. Results saved.\n")



