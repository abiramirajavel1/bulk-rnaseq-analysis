
# Chapter 9: Direction-specific GO enrichment
# Compare biological processes associated with
# upregulated and downregulated genes.

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# 1. Load annotated differential expression results
results_df <- read.csv(
  "05_annotated_differential_expression_results.csv",
  stringsAsFactors = FALSE
)

# 2. Select significant genes
significant <- subset(
  results_df,
  !is.na(padj) & padj < 0.05 & !is.na(log2FoldChange)
)

# 3. Separate upregulated and downregulated genes
upregulated <- subset(significant, log2FoldChange > 0)
downregulated <- subset(significant, log2FoldChange < 0)

cat("Significant genes:", nrow(significant), "\n")
cat("Upregulated genes:", nrow(upregulated), "\n")
cat("Downregulated genes:", nrow(downregulated), "\n")

# 4. Define the background: all genes tested in DESeq2
background <- subset(
  results_df,
  !is.na(padj) & !is.na(log2FoldChange)
)

# 5. Convert Ensembl IDs to Entrez IDs
# Remove Ensembl version suffixes, if present.
clean_ids <- function(ids) {
  sub("\\..*$", "", ids)
}

up_ids <- unique(clean_ids(upregulated$ENSEMBL))
down_ids <- unique(clean_ids(downregulated$ENSEMBL))
background_ids <- unique(clean_ids(background$ENSEMBL))

up_map <- bitr(
  up_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_map <- bitr(
  down_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_map <- bitr(
  background_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 6. Remove duplicated IDs
up_entrez <- unique(up_map$ENTREZID)
down_entrez <- unique(down_map$ENTREZID)
background_entrez <- unique(background_map$ENTREZID)

# 7. Run GO enrichment separately for each direction
run_go <- function(gene_ids, universe_ids) {
  if (length(gene_ids) < 5) {
    return(NULL)
  }

  enrichGO(
    gene = gene_ids,
    universe = universe_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

up_go <- run_go(up_entrez, background_entrez)
down_go <- run_go(down_entrez, background_entrez)

# 8. Save enrichment tables and plots
if (!is.null(up_go) && nrow(as.data.frame(up_go)) > 0) {
  write.csv(
    as.data.frame(up_go),
    "07_upregulated_GO_results.csv",
    row.names = FALSE
  )

  pdf("07_upregulated_GO_dotplot.pdf", width = 10, height = 7)
  print(dotplot(up_go, showCategory = 15) +
          ggtitle("GO enrichment: upregulated genes"))
  dev.off()
}

if (!is.null(down_go) && nrow(as.data.frame(down_go)) > 0) {
  write.csv(
    as.data.frame(down_go),
    "07_downregulated_GO_results.csv",
    row.names = FALSE
  )

  pdf("07_downregulated_GO_dotplot.pdf", width = 10, height = 7)
  print(dotplot(down_go, showCategory = 15) +
          ggtitle("GO enrichment: downregulated genes"))
  dev.off()
}

# 9. Summarize results
cat(
  "Upregulated GO terms:",
  if (is.null(up_go)) 0 else nrow(as.data.frame(up_go)),
  "\n"
)

cat(
  "Downregulated GO terms:",
  if (is.null(down_go)) 0 else nrow(as.data.frame(down_go)),
  "\n"
)

cat("Directional GO enrichment completed.\n")

# Lesson 9: Direction-specific GO enrichment

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

results_df <- read.csv(
  "05_annotated_differential_expression_results.csv"
)

# Significant genes in each direction
up_genes <- unique(results_df$ENSEMBL[
  !is.na(results_df$padj) &
  results_df$padj < 0.05 &
  results_df$log2FoldChange > 0
])

down_genes <- unique(results_df$ENSEMBL[
  !is.na(results_df$padj) &
  results_df$padj < 0.05 &
  results_df$log2FoldChange < 0
])

# Use all tested genes as the background
background <- unique(results_df$ENSEMBL[
  !is.na(results_df$padj)
])

# Convert IDs to Entrez IDs
up_map <- bitr(
  up_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_map <- bitr(
  down_genes,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_map <- bitr(
  background,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Function to run GO enrichment
run_go <- function(gene_ids, background_ids) {
  enrichGO(
    gene = unique(gene_ids),
    universe = unique(background_ids),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

up_go <- run_go(
  up_map$ENTREZID,
  background_map$ENTREZID
)

down_go <- run_go(
  down_map$ENTREZID,
  background_map$ENTREZID
)

up_df <- as.data.frame(up_go)
down_df <- as.data.frame(down_go)

cat("Significant GO terms for upregulated genes:",
    nrow(up_df), "\n")

cat("Significant GO terms for downregulated genes:",
    nrow(down_df), "\n")

write.csv(
  up_df,
  "07_upregulated_GO_results.csv",
  row.names = FALSE
)

write.csv(
  down_df,
  "07_downregulated_GO_results.csv",
  row.names = FALSE
)

if (nrow(up_df) > 0) {
  pdf("07_upregulated_GO_dotplot.pdf", width = 10, height = 7)
  print(dotplot(up_go, showCategory = 15) +
          ggtitle("GO: Upregulated Genes"))
  dev.off()
}

if (nrow(down_df) > 0) {
  pdf("07_downregulated_GO_dotplot.pdf", width = 10, height = 7)
  print(dotplot(down_go, showCategory = 15) +
          ggtitle("GO: Downregulated Genes"))
  dev.off()
}

cat("Direction-specific enrichment analysis completed.\n")

