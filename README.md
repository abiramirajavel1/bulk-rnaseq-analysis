# Bulk RNA-seq Analysis Using R and DESeq2

## Project Overview

This project demonstrates a reproducible bulk RNA-seq analysis workflow using R, Bioconductor, and DESeq2. The analysis uses the publicly available `airway` dataset to investigate gene expression differences between dexamethasone-treated and untreated human airway smooth muscle cell samples.

The project covers exploratory data analysis, normalization, principal component analysis (PCA), differential gene expression, gene annotation, and Gene Ontology (GO) functional enrichment.

## Objectives

- Explore RNA-seq count data and sample metadata.
- Examine sequencing library sizes and gene expression distributions.
- Normalize gene expression counts and assess sample variation using PCA.
- Identify differentially expressed genes between treatment groups.
- Visualize differential expression using a volcano plot.
- Annotate gene identifiers with gene symbols and descriptions.
- Perform GO Biological Process enrichment analysis.
- Compare functional enrichment patterns among upregulated and downregulated genes.

## Dataset

- **Dataset:** `airway` from Bioconductor
- **Organism:** Human (*Homo sapiens*)
- **Sample type:** Airway smooth muscle cells
- **Experimental comparison:** Dexamethasone-treated versus untreated samples
- **Number of samples:** 8
- **Initial number of genes:** 63,677

The dataset is used for learning and demonstrating RNA-seq analysis methods.

## Analysis Workflow

### 1. Data Exploration
- Load the RNA-seq count matrix and sample metadata.
- Examine matrix dimensions, sample identifiers, and library sizes.
- Explore gene expression count distributions.

### 2. Normalization and PCA
- Estimate size factors using DESeq2.
- Filter low-count genes for downstream analysis.
- Apply variance-stabilizing transformation for visualization.
- Use PCA to explore sample-level expression patterns.

### 3. Differential Gene Expression
- Fit a DESeq2 model using treatment status as the explanatory variable.
- Compare dexamethasone-treated samples against untreated samples.
- Identify significant genes using adjusted p-values.
- Examine log2 fold changes to determine expression direction.

### 4. Visualization
- Generate a volcano plot to visualize statistical significance and effect sizes.

### 5. Gene Annotation
- Map gene identifiers to gene symbols and descriptions using `org.Hs.eg.db`.

### 6. Functional Enrichment
- Perform Gene Ontology Biological Process enrichment analysis using `clusterProfiler`.
- Analyse upregulated and downregulated gene sets separately.
- Visualize enriched biological processes using dot plots.

## Key Results

Using the specified filtering and statistical criteria, the analysis identified:

| Metric | Result |
|---|---:|
| Genes retained for differential expression testing | 16,139 |
| Significantly differentially expressed genes | 2,773 |
| Upregulated genes | 1,568 |
| Downregulated genes | 1,205 |
| Enriched GO terms for upregulated genes | 186 |
| Enriched GO terms for downregulated genes | 18 |

Differential expression was assessed using an adjusted p-value threshold of 0.05. The upregulated and downregulated gene counts reflect the direction of the estimated log2 fold change.

The upregulated gene set showed enrichment for biological processes including angiogenesis, blood vessel morphogenesis, peptide hormone responses, and cell-substrate adhesion. The downregulated gene set showed enrichment for developmental and neuronal-related annotations, including ear development and axon guidance.

GO enrichment identifies overrepresented annotations; it does not by itself demonstrate that a biological process is activated or inhibited. Some enriched terms may reflect genes with roles in multiple tissues or developmental contexts.

## Repository Structure

| File | Description |
|---|---|
| `01_explore_airway.R` | Initial count matrix and metadata exploration |
| `02_normalization_pca.R` | Normalization, expression distributions, and PCA |
| `03_differential_expression.R` | Differential expression analysis using DESeq2 |
| `04_volcano_plot.R` | Volcano plot generation |
| `05_gene_annotation.R` | Gene identifier annotation |
| `06_functional_enrichment.R` | GO Biological Process enrichment |
| `07_directional_enrichment.R` | Separate enrichment analyses for upregulated and downregulated genes |
| `03_differential_expression_results.csv` | Differential expression results |
| `03_significant_genes.csv` | Significant differentially expressed genes |
| `05_annotated_differential_expression_results.csv` | Annotated differential expression results |
| `06_GO_enrichment_results.csv` | Overall GO enrichment results |
| `07_upregulated_GO_results.csv` | GO enrichment results for upregulated genes |
| `07_downregulated_GO_results.csv` | GO enrichment results for downregulated genes |
| `04_volcano_plot.pdf` | Differential expression volcano plot |
| `06_GO_enrichment_dotplot.pdf` | GO enrichment dot plot |
| `07_upregulated_GO_dotplot.pdf` | Upregulated gene enrichment dot plot |
| `07_downregulated_GO_dotplot.pdf` | Downregulated gene enrichment dot plot |

## Technologies and Tools

- **Language:** R
- **Statistical analysis:** DESeq2
- **Data infrastructure:** Bioconductor, SummarizedExperiment
- **Gene annotation:** AnnotationDbi, org.Hs.eg.db
- **Functional enrichment:** clusterProfiler
- **Visualization:** ggplot2 and enrichment plotting functions
- **Version control:** Git and GitHub

## Reproducibility

The analysis uses the Bioconductor `airway` dataset and R scripts to document the major steps of the workflow.

To reproduce the analysis:

1. Install R and an appropriate version of Bioconductor.
2. Install the required R packages.
3. Clone or download this repository.
4. Run the scripts in numerical order, beginning with data exploration and continuing through directional enrichment.

The scripts depend on the required packages being installed and on the input data and output paths being available as expected. Package versions and software environments should be documented when preparing the workflow for fully reproducible execution.

## Limitations

- The project uses a small, predefined dataset for educational and portfolio purposes.
- PCA is exploratory and does not establish causation.
- Differential expression results depend on the statistical model, filtering criteria, and significance thresholds.
- Gene annotation and enrichment analyses may be affected by identifier mapping losses, gene-set definitions, and overlapping GO categories.
- Enrichment results should be interpreted as statistical overrepresentation rather than direct evidence of biological activation or suppression.

## Learning Outcomes

Through this project, I practised:

- RNA-seq data handling and quality exploration.
- Count normalization and exploratory multivariate analysis.
- Differential expression modelling and multiple-testing correction.
- Gene identifier annotation.
- Functional interpretation using GO enrichment.
- Data visualization, scripting, and version control with Git and GitHub.

This project forms part of my ongoing development in bioinformatics, transcriptomics, and reproducible computational biology.
