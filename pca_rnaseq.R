# The goal of RNA-seq experiments is to quantify the transcripts present in a biological sample.
# This involves several steps, including data preprocessing, normalization, and analysis.
# One common challenge in RNA-seq data analysis is dealing with high-dimensional data.
# Principal Component Analysis (PCA) is a powerful technique for dimensionality reduction,
# helping to visualize and interpret the complex data by reducing it to fewer dimensions while retaining most of the variance.

library(dplyr)
library(tidyverse)
library(DESeq2)

# read raw counts data and set gene IDs as row names
raw_counts <- read.csv("tutorial/raw_counts.csv", 
                       header = T, 
                       stringsAsFactors = T) %>% 
  column_to_rownames("Geneid")

# read experimental design data
xp_design <- read.csv("tutorial/samples_to_conditions.csv", 
                        header = T, 
                        stringsAsFactors = T)

# define a custom R function called "mypca()" for Principal Component Analysis.
# PCA is used to reduce the dimensionality of the data while preserving as much variance as possible.
# This helps in visualizing the data and identifying patterns or groupings among samples.
mypca <- function(data_matrix, center = TRUE, scale = TRUE) {
  # ensure samples are in rows and variables in columns
  
  # remove constant variables
  standard_deviations <- apply(data_matrix, 2, sd)
  reduced_data_matrix <- data_matrix[, standard_deviations > 0]
  
  # perform Singular Value Decomposition (SVD)
  svd_result <- svd(scale(reduced_data_matrix, center = center, scale = scale))
  
  # create scores data frame
  principal_component_scores <- as.data.frame(svd_result$u %*% diag(svd_result$d))
  rownames(principal_component_scores) <- rownames(data_matrix)
  colnames(principal_component_scores) <- paste0("PC", 1:ncol(principal_component_scores))
  
  # create loadings data frame
  principal_component_loadings <- as.data.frame(svd_result$v)
  colnames(principal_component_loadings) <- paste0("PC", 1:ncol(principal_component_loadings))
  rownames(principal_component_loadings) <- colnames(reduced_data_matrix)
  
  # create data frame for explained variances
  explained_variances <- as.data.frame(round((svd_result$d^2) / sum(svd_result$d^2) * 100, digits = 1))
  rownames(explained_variances) <- paste0("PC", 1:ncol(principal_component_loadings))
  colnames(explained_variances) <- "explained_variance"
  
  # return result
  return(list("scores" = principal_component_scores,
              "loadings" = principal_component_loadings,
              "explained_var" = explained_variances))
}

# reorder raw_counts columns according to the experimental design file
counts <- raw_counts[, xp_design$sample]

# first five rows and five columns
counts[1:5, 1:5]

# check that sample names match between the two files and are in the same order
all(colnames(counts) %in% xp_design$sample)
all(colnames(counts) == xp_design$sample)

# convert columns to factors for proper handling in the DESeq2 model
xp_design$infected <- factor(xp_design$infected)
xp_design$dpi <- factor(xp_design$dpi)

# create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = xp_design, 
                              design = ~ growth + infected + dpi)

# Plot of mean - sd comparison
# Variance - mean plot for all genes
p_mean_sd_scaled <- 
  raw_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = log10(gene_average), y = log10(gene_stdev))) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (log10 scale)",
       y = "Gene count standard deviation (log10 scale)") +
  ggtitle("Mean - Standard deviation relationship\n(no variance stabilisation ")
p_mean_sd_scaled


# Variance stabilisation
# estimation of size factors and dispersions are required before performing the variance stabilisation

dds = estimateSizeFactors(dds)
dds = estimateDispersions(object = dds, fitType = "parametric", quiet = TRUE)
vsd = varianceStabilizingTransformation(object = dds, 
                                        blind = TRUE, # do not take the design formula into account, best practice for sample-level QC
                                        fitType = "parametric")

# extract the matrix of variance-stabilised counts
variance_stabilised_counts <- assay(vsd)

# create the mean-standard deviation plot for variance-stabilized counts
p_mean_sd_vst <- 
  variance_stabilised_counts %>% 
  as.data.frame() %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(cols = - gene, names_to = "sample", values_to = "counts") %>% 
  group_by(gene) %>% 
  summarise(gene_average = mean(counts), gene_stdev = sd(counts)) %>% 
  ungroup() %>% 
  ggplot(., aes(x = gene_average, y = gene_stdev)) +
  geom_point(alpha = 0.5, fill = "grey", colour = "black") +
  labs(x = "Gene count average (variance stabilised)", 
       y = "Gene count standard deviation (variance stabilised)") +
  ggtitle("Mean - Standard deviation relationship\n(after variance stabilisation ")
p_mean_sd_vst

# transpose the data because in variance_stabilised_counts the rows are the variables and the columns correspond to the samples
t_variance_stabilised_counts <- t(variance_stabilised_counts)

# before computing the PCA, check that samples are in rows and genes in columns
pca_results <- mypca(t_variance_stabilised_counts, 
                     center = TRUE, 
                     scale = TRUE)

# plot explained variance per principal component
ggplot(pca_results$explained_var, 
       aes(x = seq(from = 1, to = nrow(pca_results$explained_var)), 
           y = explained_variance)) +
  ylab('explained variance (%)') + 
  ggtitle('Explained variance per component') + 
  geom_bar(stat = "identity") +
  labs(x = "Principal Component number") +
  scale_x_continuous(breaks = seq(
    from = 1, 
    to = nrow(pca_results$explained_var)))

# calculate the cumulative sum of explained variance
cumsum(pca_results$explained_var)[2,1]

# determine how many principal components are needed to explain more than 50% of the variance
cumsum(pca_results$explained_var) %>% 
  as.data.frame() %>% 
  filter(explained_variance > 50) %>% 
  head(n = 1)

# extract PCA scores
scores <- pca_results$scores

# first 5 rows and columns
scores[1:5,1:5]

# join PCA scores with experimental conditions for visualization
scores_with_conditions <- 
  scores %>% 
  rownames_to_column("sample") %>% # to prepare to join on the "sample" column
  left_join(x = .,                 # this means that we are passing the 'scores' dataframe 
            y = xp_design,         # this dataframe contains the sample to condition correspondence
            by = "sample")

# shows the first 5 rows and the last 4 columns  
scores_with_conditions[1:5, 48:52]

# extract explained variance values for plotting
explained_variance <- 
  pca_results$explained_var %>% 
  pull("explained_variance")

# create a PCA score plot with infection condition overlaid
ggplot(scores_with_conditions, 
       aes(PC1, PC2, color = dpi)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ",explained_variance[1],"% variance")) +
  ylab(paste0("PC2: ",explained_variance[2],"% variance")) + 
  coord_fixed(ratio = 1) +
  ggtitle("PCA score plot with the infection condition overlaid")

# create a DESeqDataSet object for normalization
dds_count_norm <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = xp_design, 
                              design = ~ infected) 

# estimate size factors for normalization
dds_count_norm <- estimateSizeFactors(dds_count_norm)
sizeFactors(dds_count_norm)

# create a tibble for size factors
size_factors_df <- tibble(
  sample = names(sizeFactors(dds_count_norm)), 
  correction_factor = sizeFactors(dds_count_norm)
)

# line plot to connect the different size factor values
p <- ggplot(size_factors_df, aes(x = sample, y = correction_factor, group = 1)) +
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle = 90, size = 5)) +
  scale_y_continuous(limits = c(0.5,2))

p

# extract the normalised counts
counts_normalised = counts(dds_count_norm, normalized = TRUE)

# display the first five rows and columns of the normalized counts
counts_normalised[1:5,1:5]


