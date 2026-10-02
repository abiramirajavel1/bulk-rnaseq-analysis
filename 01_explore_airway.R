
# Lesson 3: Explore RNA-seq counts and metadata

# Load required packages
library(airway)
library(SummarizedExperiment)

# Load the dataset
data("airway", package = "airway")

# Extract count matrix and sample metadata
count_matrix <- assay(airway, "counts")
sample_metadata <- as.data.frame(colData(airway))

# Inspect the count matrix
cat("\n--- Count matrix dimensions ---\n")
print(dim(count_matrix))

cat("\n--- First six genes and four samples ---\n")
print(count_matrix[1:6, 1:4])

# Inspect metadata
cat("\n--- Metadata column names ---\n")
print(colnames(sample_metadata))

cat("\n--- Full sample metadata ---\n")
print(sample_metadata)

# Count samples by treatment, if the dex column exists
if ("dex" %in% colnames(sample_metadata)) {
  cat("\n--- Samples per treatment group ---\n")
  print(table(sample_metadata$dex))
}

# Check sample ID alignment
cat("\n--- Sample IDs match ---\n")
print(all(colnames(count_matrix) == rownames(sample_metadata)))

# Calculate library sizes
library_sizes <- colSums(count_matrix)

cat("\n--- Total counts per sample ---\n")
print(library_sizes)

# Plot raw count distributions
boxplot(
  log10(count_matrix + 1),
  las = 2,
  main = "Raw RNA-seq count distributions",
  ylab = "log10(count + 1)",
  xlab = "Samples"
)

# Plot total assigned counts
barplot(
  library_sizes,
  las = 2,
  main = "Total assigned counts per sample",
  ylab = "Total counts"
)
