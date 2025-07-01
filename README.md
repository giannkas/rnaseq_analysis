## README: RNAseq_Transcriptomics_Workshop

This repository contains R scripts for various analyses in an RNA-seq transcriptomics workshop. The scripts cover essential steps from data preprocessing and quality control to differential expression analysis, clustering, functional enrichment, and integration with metabolic pathways. Each script is designed to be run sequentially, building upon the results of previous steps.

### Project Overview

The primary goal of this workshop is to provide a practical understanding of RNA-seq data analysis. We analyze *Arabidopsis thaliana* leaf genes differentially expressed upon infection by *Pseudomonas syringae DC3000* after 7 days, comparing infected plants to mock-treated plants.

### Scripts Description

The repository contains the following R scripts:

1.  `pca_rnaseq.R`
2.  `de_analysis.R`
3.  `cluster_analysis.R`
4.  `enrichment_analysis.R`
5.  `transcriptomic_metabolomic_integration.R`

A brief overview of each script is provided below, with a more detailed explanation of the steps performed within each.

---

### 1. `pca_rnaseq.R`: Principal Component Analysis and Normalization

**Goal:** This script focuses on data preprocessing, normalization, and dimensionality reduction using Principal Component Analysis (PCA) to visualize and interpret complex RNA-seq data. PCA helps identify patterns and groupings among samples.

**Steps:**

* **Load Libraries:** Loads `dplyr`, `tidyverse`, and `DESeq2` for data manipulation, visualization, and differential expression analysis.
* **Load Data:**
    * Reads `raw_counts.csv` into a data frame, setting `Geneid` as row names.
    * Reads `samples_to_conditions.csv` which contains metadata about the samples, including treatment conditions.
* **Custom PCA Function (`mypca`)**: Defines a function to perform PCA. This function takes a data matrix, centers and scales it, removes constant variables, and then performs Singular Value Decomposition (SVD). It returns a list containing principal component scores, loadings, and explained variances.
* **Data Preparation for DESeq2**:
    * Reorders `raw_counts` columns to align with the `xp_design` file.
    * Checks that sample names match and are in the same order between the count data and experimental design.
    * Converts `infected` and `dpi` columns in `xp_design` to factors for proper handling by `DESeq2`.
* **DESeqDataSet Creation**: Creates a `DESeqDataSet` object from the raw counts and experimental design, specifying the design formula `~ growth + infected + dpi`.
* **Variance-Mean Relationship Visualization (Pre-VST)**: Generates a plot of the mean versus standard deviation for raw gene counts (log10 scale) to visualize the variance-mean relationship before variance stabilization.
* **Variance Stabilization Transformation (VST)**:
    * Estimates size factors and dispersions using `estimateSizeFactors()` and `estimateDispersions()` from `DESeq2`.
    * Applies `varianceStabilizingTransformation()` to the `DESeqDataSet` to stabilize variance across the range of mean expression values.
    * Extracts the matrix of variance-stabilized counts.
* **Variance-Mean Relationship Visualization (Post-VST)**: Generates a plot of the mean versus standard deviation for variance-stabilized counts to demonstrate the effect of VST.
* **PCA Calculation**: Transposes the variance-stabilized counts matrix and performs PCA using the custom `mypca()` function.
* **Explained Variance Plot**: Plots the explained variance per principal component, showing the contribution of each component to the total variance.
* **Cumulative Explained Variance**: Calculates and displays the cumulative sum of explained variance, indicating how many principal components are needed to explain more than 50% of the variance.
* **PCA Score Plot**:
    * Extracs PCA scores and joins them with experimental conditions from `xp_design`.
    * Creates a PCA score plot (PC1 vs. PC2) with `dpi` (days post-infection) overlaid to visualize sample clustering based on experimental conditions.
* **Size Factor Normalization**:
    * Creates a `DESeqDataSet` specifically for normalization, with the design `~ infected`.
    * Estimates size factors for each sample, which are used to account for differences in library size.
    * Visualizes the size factors using a line plot.
    * Extracts and displays the normalized counts.

---

### 2. `de_analysis.R`: Differential Expression Analysis

**Goal:** This script identifies genes in *Arabidopsis thaliana* that are differentially expressed in response to infection by *Pseudomonas syringae DC3000* after 7 days, comparing infected plants to mock-treated plants.

**Steps:**

* **Load Libraries:** Loads `DESeq2`, `tidyverse`, `apeglm`, and `EnhancedVolcano`.
* **Load Data:**
    * Imports `samples_to_conditions.csv` containing sample metadata.
    * Filters the experimental design to keep only "mock" and "infected P. syringae at 7 dpi" conditions.
    * Imports `raw_counts.csv` and sets `Geneid` as row names.
    * Reorders and filters raw count columns to match the selected samples.
* **DESeqDataSet Creation**: Creates a `DESeqDataSet` object using the filtered raw counts and experimental design, with the design formula `~ infected`.
* **Differential Expression Analysis**:
    * Runs the `DESeq` pipeline to perform differential expression analysis.
    * Extracts results, specifically comparing "Pseudomonas_syringae_DC3000" vs. "mock" conditions.
* **Results Exploration**:
    * Displays a summary of the `DESeqResults` object.
    * Filters genes based on adjusted p-value (`padj < 0.01` and `padj < 0.001`).
    * Plots histograms of adjusted and non-adjusted p-values to visualize their distributions.
    * Extracts differentially expressed genes (padj < 0.01) and sorts them by log2 fold change.
* **Log Fold Change Shrinkage**: Applies `lfcShrink()` to the results to improve the accuracy of log2 fold change estimates, especially for genes with low counts.
* **Volcano Plots**:
    * Generates a basic volcano plot using `EnhancedVolcano` to visualize differential expression.
    * Creates a customized volcano plot with specific `xlim`, `ylim`, `pCutoff`, and `FCcutoff` values, and custom legend labels.
* **Median of Ratios Normalization**:
    * Defines a custom `mor_normalization` function to perform median of ratios normalization. This function takes raw count data, calculates geometric means, ratios to a pseudo-reference, and then uses sample medians to derive scaling factors for normalization.
    * Applies this normalization to the filtered raw count data.
* **Normalization Comparison Visualization**:
    * Converts both scaled and raw counts to a "long" format.
    * Combines them and creates a violin plot to visually compare the distributions of raw and scaled gene counts.
* **Heatmap Visualization**:
    * Normalizes counts, filters for differentially expressed genes, and prepares the data for heatmap visualization.
    * Creates heatmaps of log2 transformed normalized counts, both without clustering and with clustering applied, showing the expression patterns of differentially expressed genes.
    * Scales the normalized counts to Z-scores for better visualization in heatmaps.
    * Plots the distribution of Z-score values to confirm that scaling was applied correctly (mean around zero).
    * Generates heatmaps with sample annotations and custom colors to visualize clustering by experimental conditions.
    * Filters for genes with a high log2 fold change (greater than 4) and creates a heatmap specifically for these genes.
* **MA Plots**: Generates MA plots (Mean-Average plots) before and after log fold change shrinkage to assess the effect of shrinkage on the spread of log fold changes.

---

### 3. `cluster_analysis.R`: Gene Co-Expression Clustering

**Goal:** This script performs hierarchical and k-means clustering on gene expression data to identify clusters of co-expressed genes, inferring potential functional relationships and regulatory mechanisms.

**Steps:**

* **Load Libraries:** Loads `tidyverse`, `gplots`, `dendextend`, `reshape2`, and `cluster`.
* **Load Data:**
    * Imports `samples_to_conditions.csv` (experimental design metadata).
    * Filters the experimental design to keep only samples at 7 days post-infection (dpi).
    * Imports `raw_counts.csv` containing raw gene count data.
* **Data Preparation**:
    * Extracts gene names and count data, setting gene names as row names.
    * Reorders count columns to match the filtered experimental design.
    * Performs checks to ensure sample names align correctly between count data and metadata.
* **Normalization and Differential Expression Selection Function (`normalize`)**:
    * Defines a function `normalize()` that takes count data, experimental design, fold change threshold, and p-value threshold.
    * It prepares `coldata` for `DESeq2`, creates a `DESeqDataSet`, and runs the `DESeq` pipeline.
    * Iterates through all possible condition pairs to identify differentially expressed genes based on specified thresholds for adjusted p-value and log2 fold change.
    * Returns normalized counts for the unique set of differentially expressed genes.
* **Apply Normalization and Selection**: Calls the `normalize()` function with the raw counts and filtered experimental design to obtain differentially expressed genes.
* **Data Scaling**: Converts the differentially expressed gene data to a matrix and scales it using `t(scale(t(DEgenes)))` for clustering.
* **Hierarchical Clustering - Samples**: Performs hierarchical clustering on samples using Spearman correlation to create a sample tree. Plots the sample dendrogram.
* **Hierarchical Clustering - Genes**: Performs hierarchical clustering on genes using Pearson correlation to create a gene tree. Plots the gene dendrogram.
* **Heatmap Visualization**: Creates a heatmap of the differentially expressed genes, with rows and columns ordered by the hierarchical clustering results.
* **Dendrogram Cutting and Visualization (Hierarchical)**:
    * Cuts the gene hierarchical tree at different heights (1.5, 1.0, 0.5) to form gene clusters.
    * Displays the head of the clusters formed at height 0.5 and counts the number of unique clusters.
    * Plots the gene tree with colored bars indicating the clusters at different cut heights, and adds lines to show the cut heights.
    * Alternatively, cuts the tree to form a specified number of clusters (e.g., k=5).
    * Visualizes these k=5 clusters on the gene dendrogram.
* **Heatmap with Cluster Bars**: Generates a heatmap with colored side bars indicating the gene clusters obtained from hierarchical clustering (e.g., k=5) to help visualize if the chosen number of clusters aligns with visual patterns in the heatmap.
* **Core Expression Pattern Calculation and Plotting**:
    * Defines a function `clust.core()` to calculate the mean expression pattern for each cluster.
    * Applies this function to obtain core expression patterns for each cluster.
    * Reshapes the core expression patterns and plots them to show the average expression trend of genes within each cluster across samples.
* **Individual Cluster Visualization**:
    * Subsets data for a specific cluster (e.g., cluster 4).
    * Reshapes the data for plotting.
    * Plots the individual gene expression patterns within cluster 4, overlaid with the core expression pattern of that cluster.
* **Determining Optimal Number of Clusters (Elbow Method)**:
    * Calculates the within-groups sum of squares (WSS) for a range of cluster numbers (2 to 20).
    * Plots the WSS, which helps identify the "elbow" point, suggesting an optimal number of clusters.
* **Determining Optimal Number of Clusters (Silhouette Method)**:
    * Calculates the average silhouette width for a range of cluster numbers (2 to 20).
    * Plots the average silhouette width and adds a vertical line to indicate the number of clusters with the maximum average silhouette width, suggesting an optimal number of clusters.
* **K-means Clustering**: Performs k-means clustering with the optimal number of clusters (e.g., 4, as indicated by silhouette analysis).
* **Comparing Hierarchical and K-means Clusters**: Plots the gene tree with colored bars comparing the clusters obtained from hierarchical tree cutting and k-means clustering.

---

### 4. `enrichment_analysis.R`: Functional Enrichment Analysis

**Goal:** This script performs functional enrichment analysis, specifically Over-Representation Analysis (ORA) for Gene Ontology (GO) terms, to identify biological pathways, processes, or functions overrepresented in the set of differentially expressed genes.

**Steps:**

* **Load Libraries:** Loads `DESeq2`, `tidyverse`, `apeglm`, `biomartr`, `clusterProfiler`, `enrichplot`, `org.At.tair.db`, `biomaRt`, and `splitstackshape`.
* **Load Differentially Expressed Genes**: Reads `differential_genes.csv`.
* **Quantile Calculation and Filtering**:
    * Calculates quantiles of log2 fold changes for upregulated genes.
    * Filters genes with log2 fold changes above the 75th percentile and saves them to `.tsv` files (`diff_genes_for_agrigo.tsv` and `diff_genes_for_agrigo_page.tsv`) for potential use with AgriGO.
* **BioMart Queries for Arabidopsis Annotation**:
    * Checks available datasets for *Arabidopsis thaliana* in BioMart.
    * Retrieves available attributes for *A. thaliana* in BioMart, specifically for the `athaliana_eg_gene` dataset.
    * Queries BioMart to retrieve TAIR symbols and Entrez Gene IDs for the initial set of differentially expressed genes.
* **Background Gene Set Preparation**:
    * Reads all *Arabidopsis* genes from `raw_counts.csv` to create a background/universe set.
    * Queries BioMart to retrieve annotations (TAIR symbol, UniProt ID, Entrez Gene ID) for all *Arabidopsis* genes in the background.
    * Ensures Entrez Gene IDs in the background set are character type for compatibility with `enrichGO`.
* **Differentially Expressed Gene Annotation**: Queries BioMart again to retrieve annotations for the specifically differentially expressed genes.
* **Over-Representation Analysis (ORA) with `enrichGO`**:
    * Performs ORA for Gene Ontology Biological Process (BP) terms using `enrichGO()`.
    * Uses the Entrez Gene IDs of differentially expressed genes as the gene set and all annotated *Arabidopsis* genes as the universe.
    * Specifies `org.At.tair.db` for annotation, "ENTREZID" as `keyType`, "BP" for ontology, "BH" (Benjamini-Hochberg) for p-value adjustment, and a `qvalueCutoff` of 0.05.
* **ORA Results Processing and Visualization**:
    * Writes the ORA results to a `.tsv` file (`go_results.tsv`).
    * Simplifies the ORA results using `clusterProfiler::simplify()` to remove redundant GO terms.
    * Generates a dotplot to visualize the simplified ORA results, showing enriched GO terms.
    * Calculates pairwise term similarity for the simplified ORA results using `pairwise_termsim()`.
    * Creates an `emapplot()` to visualize the relationships and similarity between enriched GO terms.

---

### 5. `transcriptomic_metabolomic_integration.R`: Transcriptomic and Metabolomic Integration

**Goal:** This script integrates transcriptomic data with metabolic pathways to gain deeper insights into potential metabolic changes and regulatory mechanisms. It prepares data for visualization using tools like the Interactive Pathways Explorer (iPath).

**Steps:**

* **Load Libraries:** Loads `biomartr`, `tidyverse`, and `biomaRt`.
* **Load Differentially Expressed Genes**: Reads `differential_genes.csv`.
* **Quantile Calculation**: Calculates quantiles of log2 fold changes for upregulated genes.
* **BioMart Query for All Differentially Expressed Genes**:
    * Defines attributes to retrieve from BioMart: TAIR symbol, UniProt ID, and Entrez Gene ID.
    * Queries BioMart (`plants_mart`, `athaliana_eg_gene`) to retrieve these annotations for all differentially expressed genes.
* **Prepare UniProt IDs for iPath (All DE Genes)**:
    * Filters the annotated genes to keep only those with UniProt IDs.
    * Ensures unique entries.
    * Formats the UniProt IDs by prepending "UNIPROT:" to create an ID that iPath can use.
    * Selects only the formatted ID column and writes it to a `.tsv` file (`diff_genes_swissprot.tsv`). This file can be copy-pasted into iPath for visualization.
* **Filter for Highly Upregulated Genes**:
    * Filters the `diff_genes` data to include only those with `log2FoldChange` greater than the 75th percentile of upregulated genes.
* **BioMart Query for Highly Upregulated Genes**: Queries BioMart again to retrieve annotations (TAIR symbol, UniProt ID, and Entrez Gene ID) specifically for these highly upregulated genes.
* **Prepare UniProt IDs for iPath (Highly Upregulated Genes)**:
    * Filters, makes unique, and formats UniProt IDs for these highly upregulated genes, similar to the previous step.
    * Writes these formatted IDs to a separate `.tsv` file (`diff_genes_swissprot_2.tsv`). This file is intended for a more focused pathway visualization in iPath.
