# Differential expression analysis is a crucial step in understanding how gene expression changes
# in response to different conditions. We are interested in identifying genes
# in Arabidopsis that are differentially expressed in response to infection by the pathogenic
# bacteria Pseudomonas syringae DC3000 after 7 days. This involves comparing the transcriptome
# of infected plants to that of mock-treated plants.

library(DESeq2) # For differential expression analysis
library(tidyverse) # For data manipulation and visualization
library(apeglm) # For log fold change shrinkage
library(EnhancedVolcano) # For creating volcano plots

# import the samples to conditions correspondence
# this file contains metadata about the samples, including treatment conditions
xp_design <- read.csv("tutorial/samples_to_conditions.csv", 
                      header = T, 
                      stringsAsFactors = F, 
                      colClasses = rep("character",4))

# filter design file to keep only "mock" and the "infected P. syringae at 7 dpi" conditions
# this step ensures we only analyze the relevant conditions for our comparison
xp_design_mock_vs_infected = xp_design %>% 
  filter(growth == "MgCl2" & dpi == "7")

# import the gene raw counts
# this file contains the raw count data for each gene in each sample
raw_counts <- read.csv("tutorial/raw_counts.csv", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("Geneid")


# reorder counts columns according to the complete list of samples
# this ensures that the count data aligns with the metadata.
raw_counts <- raw_counts[ , xp_design$sample]

# filter count file accordingly (to keep only samples present in the filtered xp_design file)
raw_counts_filtered = raw_counts[, colnames(raw_counts) %in% xp_design_mock_vs_infected$sample]

# creation of the DESeqDataSet
# DESeqDataSet is the primary object used in DESeq2 for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = raw_counts_filtered, 
                              colData = xp_design_mock_vs_infected, 
                              design = ~ infected)

# check the levels of the 'infected' factor
dds$infected

# run the DESeq2 pipeline to perform differential expression analysis
dds <- DESeq(dds)

# extract the results of the differential expression analysis
res <- results(dds)

# have a peek at the DESeqResults object
# this object contains the log2 fold changes and p-values for each gene
res

# extract results for all genes with a specific contrast
# this specifies the comparison of interest: infected vs. mock
all_genes_results <- results(dds, contrast = c("infected",                      # name of the factor
                                               "Pseudomonas_syringae_DC3000",    # name of the numerator level for fold change
                                               "mock"))                          # name of the denominator level
# check if the results are the same
all.equal(res, all_genes_results)

head(all_genes_results)

# example calculation of log2 fold change
untreated = 230
treated = 750
log2(treated/untreated)

# display metadata columns of the results
mcols(all_genes_results)

# filter genes with a corrected p-value (padj) less than 0.01
all_genes_results %>% 
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  dim()

# threshold of p = 0.001
all_genes_results %>% 
  as.data.frame() %>% 
  filter(padj < 0.001) %>% 
  dim()

# distribution of adjusted p-values
hist(all_genes_results$padj, col="lightblue", main = "Adjusted p-value distribution")

# distribution of non-adjusted p-values
hist(all_genes_results$pvalue, col="grey", main = "Non-adjusted p-value distribution")

# extract differentially expressed genes with padj < 0.01 and sort by log2 fold change
diff_genes = all_genes_results %>% 
  as.data.frame() %>% 
  rownames_to_column("genes") %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))
head(diff_genes)

# apply log fold change shrinkage to improve accuracy of fold change estimates
resLFC <- lfcShrink(dds = dds, 
                    res = all_genes_results,
                    type = "normal",
                    coef = "infected_Pseudomonas_syringae_DC3000_vs_mock") # name or number of the coefficient (LFC) to shrink

resultsNames(dds)

# create a volcano plot to visualize differential expression results
EnhancedVolcano(toptable = resLFC,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(resLFC)
)

# create a customized volcano plot
EnhancedVolcano(toptable = resLFC,
                x = "log2FoldChange",
                y = "padj",
                lab = rownames(resLFC),
                xlim = c(-10, +10),
                ylim = c(0,100),
                pCutoff = 1e-06,
                pointSize = 2.0,
                FCcutoff = 2, 
                title = "Pseudomonas syringae DC3000 versus mock \n (fold change cutoff = 2, p-value cutoff = 1e-06)",
                legendLabels = c(
                  'Not significant',
                  'Log2 fold-change (but do not pass p-value cutoff)',
                  'Pass p-value cutoff',
                  'Pass both p-value & Log2 fold change')
)

# median of ratios normalization method
mor_normalization = function(data){
  
  # take the log
  log_data = log(data) 
  
  # convert to a data frame and calculate the geometric mean for each gene
  log_data = log_data %>% 
    rownames_to_column('gene') %>% 
    mutate (gene_averages = rowMeans(log_data)) %>% 
    filter(gene_averages != "-Inf")
  
  # the last columns is the pseudo-reference column 
  pseudo_column = ncol(log_data)
  
  # where to stop before the pseudo column 
  before_pseduo = pseudo_column - 1
  
  # find the ratio of the log data to the pseudo-reference
  ratios = sweep(log_data[,2:before_pseduo], 1, log_data[,pseudo_column], "-")
  
  # find the median of the ratios
  sample_medians = apply(ratios, 2, median)
  
  # convert the median to a scaling factor
  scaling_factors = exp(sample_medians)
  
  # use scaling factors to scale the original data
  manually_normalized = sweep(data, 2, scaling_factors, "/")
  
  # return the normalized data
  return(manually_normalized)
}

# apply median of ratios normalization to the filtered count data
scaled_counts <- mor_normalization(raw_counts_filtered)
head(scaled_counts)

# convert scaled counts to long format for visualization
long_scaled_counts = 
  scaled_counts %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>% 
  mutate(scaled = "yes")

# convert raw counts to long format for visualization
long_raw_counts = 
  raw_counts %>% 
  rownames_to_column("gene") %>% 
  pivot_longer(-gene, names_to = "sample", values_to = "counts") %>% 
  mutate(scaled = "no")

# combine raw and scaled counts for comparison
long_raw_and_scaled_counts = bind_rows(long_raw_counts, long_scaled_counts)

# create a violin plot to compare raw and scaled counts
ggplot(long_raw_and_scaled_counts, 
       aes(x = scaled, y = counts + 1, fill = scaled)) +
  geom_violin() +
  scale_y_log10() +
  labs(x = "Gene counts scaled/normalised?", y = "Gene counts (raw or scaled)")

# normalize counts and filter for differentially expressed genes
counts_normalised_only_diff_genes = 
  mor_normalization(raw_counts) %>%             # normalize the counts using our custom function
  rownames_to_column("genes") %>%               
  pivot_longer(- genes,                         
               names_to = "sample", 
               values_to = "counts") %>% 
  filter(genes %in% diff_genes$genes) %>%             
  pivot_wider(names_from = "sample",            
              values_from = "counts")  %>%      
  column_to_rownames("genes")                   # the gene column is converted back to row names to create a matrix usable with pheatmap

# check the dimensions of the filtered and normalized count data
dim(counts_normalised_only_diff_genes)          # check that you have the expected number of rows and columns

# create a heatmap of log2 transformed normalized counts for differentially expressed genes
pheatmap(log2(counts_normalised_only_diff_genes+1), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "row",
         show_rownames = FALSE, 
         show_colnames = TRUE)

# scale the normalized counts for heatmap visualization
counts_scaled = 
  counts_normalised_only_diff_genes %>% 
  t(.) %>%                              # transpose to have the genes in columns 
  scale() %>%                           # scale(x, center = TRUE, scale = TRUE) 
  t(.)                                  # back in original shape

# plot the distribution of Z-score values: the majority of the values should be around zero
apply(counts_scaled, MARGIN = 1, mean) %>%                          # calculate the mean per row
  hist(., main = "", xlab = "Z-score values", col = "dodgerblue2")  

# create a heat map of the scaled counts
pheatmap(counts_scaled, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = FALSE, 
         scale = "none",            # already done "manually"
         show_colnames = TRUE)

# applying clustering
pheatmap(counts_scaled, 
         cluster_rows = TRUE,                      
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         main = "Clustering on")

# filter counts for only the samples of interest
counts_scaled_filtered = 
  counts_scaled %>% 
  as.data.frame() %>%
  dplyr::select(xp_design_mock_vs_infected$sample) # keep the 8 samples

# display the first few rows of the filtered and scaled counts
head(counts_scaled_filtered)

# prepare annotation information for the heatmap
anno_col_info = xp_design_mock_vs_infected %>% column_to_rownames("sample")

# definning colors
anno_info_colors = list(
  seed = c(MgCl2 = "#d8b365"),
  infected = c(mock = "lightgrey", 
               Pseudomonas_syringae_DC3000 = "black"),
  dpi = c("7" = "dodgerblue4")
)

# heatmap with colors and annotations
pheatmap(counts_scaled_filtered, 
         cluster_rows = TRUE,                       
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = TRUE,
         annotation_col = anno_col_info,
         annotation_colors = anno_info_colors,
         main = "Clustering with ward method")

# filter genes with a fold change higher than 4 (positive)
genes_differential_fold = 
  res %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%  
  filter(log2FoldChange > 4) %>% 
  select(gene) 

# filter scaled counts for genes with high fold change
counts_scaled_filtered_high_fold_change = 
  counts_scaled_filtered[row.names(counts_scaled_filtered) %in% genes_differential_fold$gene, ]

# heatmap with clustering for genes with high fold change
pheatmap(counts_scaled_filtered_high_fold_change, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         annotation_col = anno_col_info,
         clustering_method = "average",
         annotation_colors = anno_info_colors,
         main = "Clustering of genes and samples (high fold only)")

# crate MA plot to visualize differential expression results
plotMA(dds, alpha = 0.01)

# apply log fold change shrinkage to improve accuracy 
resLFC <- lfcShrink(dds = dds, 
                    res = res,
                    type = "normal",
                    coef = 2) # corresponds to "infected_Pseudomonas_syringae_DC3000_vs_mock" comparison

plotMA(resLFC, alpha = 0.01)
