
# Lesson 8: Gene Ontology Enrichment Analysis

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# Load annotated differential expression results
results_df <- read.csv(
  "05_annotated_differential_expression_results.csv"
)

# Select significantly differentially expressed genes
significant_genes <- results_df$ENSEMBL[
  !is.na(results_df$padj) &
  results_df$padj < 0.05
]

# Define the background as all genes tested in the analysis
background_genes <- results_df$ENSEMBL[
  !is.na(results_df$padj)
]

# Convert Ensembl IDs to Entrez IDs
sig_map <- bitr(
  unique(significant_genes),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_map <- bitr(
  unique(background_genes),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Perform Gene Ontology enrichment for biological processes
go_results <- enrichGO(
  gene = unique(sig_map$ENTREZID),
  universe = unique(background_map$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Display results
go_df <- as.data.frame(go_results)

cat("Number of enriched biological processes:",
    nrow(go_df), "\n")

if (nrow(go_df) > 0) {
  cat("\nTop enriched biological processes:\n")
  print(head(
    go_df[, c("Description", "GeneRatio",
              "Count", "p.adjust")],
    10
  ))

  write.csv(
    go_df,
    "06_GO_enrichment_results.csv",
    row.names = FALSE
  )

  pdf("06_GO_enrichment_dotplot.pdf",
      width = 10, height = 7)

  print(dotplot(
    go_results,
    showCategory = 15
  ) + ggtitle("GO Biological Process Enrichment"))

  dev.off()

  cat("\nEnrichment table and dot plot saved.\n")
} else {
  cat("No GO terms passed the selected significance thresholds.\n")
}

