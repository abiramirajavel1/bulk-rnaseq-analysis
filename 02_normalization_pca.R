
# Lesson 4: RNA-seq Normalization and PCA
# Dataset: Airway RNA-seq dataset
# Purpose: Normalize counts and explore sample relationships

# --------------------------------------------------
# 1. Install required packages if missing
# --------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("airway", quietly = TRUE)) {
  BiocManager::install("airway")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

# --------------------------------------------------
# 2. Load packages
# --------------------------------------------------

library(airway)
library(SummarizedExperiment)
library(DESeq2)
library(ggplot2)

# --------------------------------------------------
# 3. Load the dataset
# --------------------------------------------------

data("airway", package = "airway")

# Extract count matrix and sample metadata
count_matrix <- assay(airway, "counts")
sample_metadata <- as.data.frame(colData(airway))

# Ensure sample IDs match
stopifnot(
  all(colnames(count_matrix) == rownames(sample_metadata))
)

# Set treatment as a categorical variable
sample_metadata$dex <- factor(sample_metadata$dex)

cat("Count matrix dimensions:\n")
print(dim(count_matrix))

cat("\nSamples per treatment group:\n")
print(table(sample_metadata$dex))

# --------------------------------------------------
# 4. Create a DESeq2 dataset
# --------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ dex
)

cat("\nDESeq2 dataset:\n")
print(dds)

# --------------------------------------------------
# 5. Filter genes with very low counts
# --------------------------------------------------

# Retain genes with at least 10 counts in at least
# as many samples as the smallest treatment group.
min_group_size <- min(table(sample_metadata$dex))

keep <- rowSums(counts(dds) >= 10) >= min_group_size
dds <- dds[keep, ]

cat("\nGenes remaining after filtering:\n")
print(nrow(dds))

# --------------------------------------------------
# 6. Estimate size factors for normalization
# --------------------------------------------------

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

cat("\nSize factors:\n")
print(sizeFactors(dds))

cat("\nFirst six normalized genes and four samples:\n")
print(round(normalized_counts[1:6, 1:4], 2))

# --------------------------------------------------
# 7. Variance-stabilizing transformation
# --------------------------------------------------

vsd <- vst(dds, blind = TRUE)

# Extract transformed expression values
vst_counts <- assay(vsd)

cat("\nDimensions of transformed expression matrix:\n")
print(dim(vst_counts))

# --------------------------------------------------
# 8. Plot normalized count distributions
# --------------------------------------------------

pdf("02_normalized_count_distributions.pdf")

boxplot(
  log2(normalized_counts + 1),
  las = 2,
  main = "Normalized RNA-seq count distributions",
  ylab = "log2(normalized count + 1)",
  xlab = "Samples"
)

dev.off()

# --------------------------------------------------
# 9. Perform PCA
# --------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "dex",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = dex,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_variance[1], "% variance")) +
  ylab(paste0("PC2: ", percent_variance[2], "% variance")) +
  ggtitle("PCA of Airway RNA-seq Samples") +
  theme_bw()

print(pca_plot)

ggsave(
  filename = "02_pca_plot.pdf",
  plot = pca_plot,
  width = 7,
  height = 5
)

# --------------------------------------------------
# 10. Save normalized counts
# --------------------------------------------------

write.csv(
  normalized_counts,
  file = "02_normalized_counts.csv"
)

cat("\nLesson 4 completed successfully!\n")
cat("Output files:\n")
cat("- 02_normalized_count_distributions.pdf\n")
cat("- 02_pca_plot.pdf\n")
cat("- 02_normalized_counts.csv\n")

