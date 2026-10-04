
# Lesson 7: Gene Annotation
# Convert Ensembl gene IDs to gene symbols and descriptions

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(AnnotationDbi)
library(org.Hs.eg.db)

# Load differential expression results
results_df <- read.csv(
  "03_differential_expression_results.csv"
)

# Remove Ensembl version suffixes, if present
results_df$ENSEMBL <- sub(
  "\\..*$",
  "",
  results_df$gene_id
)

# Retrieve gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = unique(results_df$ENSEMBL),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Retrieve gene descriptions
gene_descriptions <- mapIds(
  org.Hs.eg.db,
  keys = unique(results_df$ENSEMBL),
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Add annotations to the results table
results_df$gene_symbol <- unname(
  gene_symbols[results_df$ENSEMBL]
)

results_df$gene_description <- unname(
  gene_descriptions[results_df$ENSEMBL]
)

# Display the top 10 annotated results
cat("\n--- Top 10 annotated genes ---\n")

print(
  head(
    results_df[
      order(results_df$padj, na.last = TRUE),
      c(
        "gene_id",
        "gene_symbol",
        "gene_description",
        "log2FoldChange",
        "pvalue",
        "padj"
      )
    ],
    10
  )
)

# Save the annotated results
write.csv(
  results_df,
  "05_annotated_differential_expression_results.csv",
  row.names = FALSE
)

cat("\nGene annotation completed successfully.\n")

