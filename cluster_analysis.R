# Clustering is a powerful technique for generating hypotheses and exploring 
# data in biological research. The underlying assumption is that genes with 
# similar expression patterns (co-expressed genes) are often controlled by the same regulatory mechanisms (co-regulated genes).
# By identifying clusters of co-expressed genes, we can infer potential functional 
# relationships and regulatory mechanisms. This script performs hierarchical clustering 
# and k-means clustering on gene expression data to identify clusters of co-expressed genes.

library(tidyverse) # For data manipulation and visualization
library(gplots) # For heatmap visualization
library(dendextend) # For dendrogram manipulation
library(reshape2) # For data reshaping
library(cluster) # For clustering algorithms

# import the samples to conditions correspondence
# this file contains metadata about the samples, including treatment conditions
xp_design <- read.csv("tutorial/samples_to_conditions.csv", 
                      header = T, 
                      stringsAsFactors = F, 
                      colClasses = rep("character",4))

# filter design file to keep only samples at 7 days post-infection (dpi)
# this step ensures we only analyze the relevant conditions for our comparison
xp_design_dpi = xp_design %>% 
  filter(dpi == "7")

# import the gene raw counts
# this file contains the raw count data for each gene in each sample
raw_counts <- read.csv("tutorial/raw_counts.csv", header = T, stringsAsFactors = F)

# extract gene names and count data
genes <- raw_counts[,1]
counts <- raw_counts[,-1]
row.names(counts) <- genes

# reorder counts columns according to the complete list of samples 
# this ensures that the count data aligns with the metadata
counts <- counts[ , xp_design_dpi$sample]

# check samples names match between the two files and are in the same order
# this is crucial for ensuring that the data and metadata align correctly
all(colnames(counts) %in% xp_design_dpi$sample)
all(colnames(counts) == xp_design_dpi$sample)

# function to normalize and select differentially expressed genes
normalize <- function(countdata, xp_design, fold_change_threshold, p_value_threshold) {
  # extract all column names
  col.names <- colnames(xp_design)
  
  # columns to discard
  columns.to.discard <- c("sample", "fq1", "fq2", "fq")
  
  # only keep the columns of interest
  colsForConditions <- col.names[!col.names %in% columns.to.discard]
  
  # determine conditions based on the number of condition columns
  if (length(colsForConditions) == 1) {
    condition <- factor(xp_design[, colsForConditions])
  } else {
    # combine multiple condition columns into a single condition column
    xp_design$conditions <- apply(xp_design[, colsForConditions], 1, paste, collapse = ".")
    condition <- factor(x = xp_design[, "conditions"], levels = unique(xp_design[, "conditions"]))
  }
  
  # create a coldata frame and instantiate the DESeqDataSet
  coldata <- data.frame(row.names = colnames(countdata), condition)
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
  
  # run the DESeq pipeline
  dds <- DESeq(dds)
  
  # create dataframe containing normalized counts
  resdata <- as.data.frame(counts(dds, normalized = TRUE))
  
  # initialize list to store differentially expressed genes
  total_selec <- list()
  
  # iterate through the list of conditions to create differential expression values for all possible pairs
  condition_levels <- levels(condition)
  for (i in 1:(length(condition_levels) - 1)) {
    for (j in (i + 1):length(condition_levels)) {
      res <- results(dds, contrast = c("condition", condition_levels[i], condition_levels[j]))
      res$genenames <- rownames(res)
      resul <- as.data.frame(res)
      significantDE <- resul %>% filter(padj < p_value_threshold &
                                          (log2FoldChange > fold_change_threshold |
                                             log2FoldChange < -fold_change_threshold))
      selec <- as.list(significantDE$genenames)
      total_selec <- append(total_selec, selec)
    }
  }
  
  # extract unique differentially expressed genes
  total_selec <- unique(unlist(total_selec))
  
  # return normalized counts for differentially expressed genes
  selection <- resdata[total_selec, ]
  return(selection)
}

# normalize the count data and select differentially expressed genes
DEgenes <- normalize(counts,xp_design_dpi,2,0.01)

# check the dimensions of the differentially expressed gene matrix
dim(DEgenes)

# convert DEgenes to a matrix
DEgenes <- as.matrix(DEgenes)

# scale the data for clustering
scaledata <- t(scale(t(DEgenes)))

# Hierarchical tree of the samples using Spearman correlation
hs <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
sampleTree = as.dendrogram(hs, method="average")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")

# Hierarchical tree of the genes using Pearson correlation
hg <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
geneTree = as.dendrogram(hg, method="average")
plot(geneTree,
     leaflab = "none",             
     main = "Gene Clustering",
     ylab = "Height",
     xlab = "Genes")

# create a heatmap of the differentially expressed genes
heatmap.2(DEgenes,
          Rowv=as.dendrogram(hg), 
          Colv=as.dendrogram(hs),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none")

# cut the hierarchical tree at different heights to form clusters
hclusth1.5 = cutree(hg, h=1.5) #cut tree at height of 1.5
hclusth1.0 = cutree(hg, h=1.0) #cut tree at height of 1.0
hclusth0.5 = cutree(hg, h=0.5) #cut tree at height of 0.5

# display the first few elements of the clusters formed at height 0.5
head(hclusth0.5)

# count the number of unique clusters formed at height 0.5
length(unique(hclusth0.5))

#plot the tree
plot(geneTree,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

# combine cluster vectors
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)

# add colored bars to the dendrogram to indicate clusters
colored_bars(the_bars, geneTree, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("h=0.5","h=1.0","h=1.5"),cex.rowLabels=0.7)

# add lines showing the cut heights
abline(h=1.5, lty = 2, col="grey")
abline(h=1.0, lty = 2, col="grey")
abline(h=0.5, lty = 2, col="grey")

# alternatively it is also posible to set the number of cluster you want.
hclustk5 = cutree(hg, k=5)

plot(geneTree,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

# add colored bars to the dendrogram to indicate clusters formed by cutting the tree into 5 clusters
colored_bars(hclustk5, geneTree, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("k=5"),cex.rowLabels=0.7)

# It is also posible to visualize the bar indicating the clusters in combination with a heatmap.
# This will make it more easy to see if you choose the right number of clusters.
clustColBar <- rainbow(length(unique(hclustk5)), start=0.1, end=0.9)
clustColBar <- clustColBar[as.vector(hclustk5)]

heatmap.2(DEgenes,
          Rowv=as.dendrogram(hg), 
          Colv=as.dendrogram(hs),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none",
          RowSideColors=clustColBar,
          key = FALSE)

# function to calculate the core expression pattern of a cluster
clust.core = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

# calculate the core expression patterns for each cluster
clusters <- hclustk5
cores <- sapply(unique(clusters), clust.core, scaledata, clusters)

# reshape the core expression patterns for plotting
moltenCores <- melt(cores) ## get the data in a "long" format
colnames(moltenCores) <- c('sample','cluster','value')

# plot the core expression patterns for each cluster
ggplot(moltenCores, aes(x=sample,y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Sample") +
  ylab("Expression") +
  labs(title= "Cluster Expression of the samples",color = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# subset the data for cluster 4
clust4 <- t(scaledata[hclustk5==4,])

# get the data frame into long format for plotting
clust4Molten <- melt(clust4, id.vars = "Sample")
colnames(clust4Molten) <- c('sample','gene','value')

# subset the cores molten dataframe so we can plot core4
core <- moltenCores[moltenCores$cluster==4,]

# plot the expression patterns for clusters 4
ggplot(clust4Molten, aes(x=sample,y=value)) + 
  geom_line(color="grey", aes(color="grey", group=gene)) +
  #this adds the core 
  geom_line(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
  geom_point(data=core, aes(sample,value, group=cluster), color="blue",inherit.aes=FALSE) +
  xlab("Samples") +
  ylab("Expression") +
  labs(title= paste0("Cluster 4 consisting ", nrow(clust4), " genes"),color = "Score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# calculate the within groups sum of squares for different number of clusters
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
for (i in 2:20) wss[i] <- sum(kmeans(scaledata,
                                     centers=i)$withinss)

# plot the within groups sum of squares
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# calculate the average silhouette width for different numbers of clusters
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}

# plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")

# add a vertical line to indicate the optimal number of clusters
abline(v = which.max(sil), lty = 2)

# print the optimal number of clusters based on the average silhouette width
cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")

# perform k-means clustering with the optimal number of clusters
set.seed(20)
kClust <- kmeans(scaledata, centers=4, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster

# plot the hierarchical tree of the genes
plot(geneTree,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

# combine the cluster vectors from hierarchical clustering and k-means clustering
the_bars <- cbind(hclustk5, kClusters)

# add colored bars to the dendrogram to indicate clusters from hierarchical clustering and k-means clustering
colored_bars(the_bars, geneTree, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("Treecut",'K-means'),cex.rowLabels=0.7)
