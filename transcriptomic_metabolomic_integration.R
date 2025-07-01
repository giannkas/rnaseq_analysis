# While we have analyzed transcriptomic data in isolation, integrating 
# this data with metabolic pathways can provide deeper insights.
# The Interactive Pathways Explorer (iPath) is a web-based tool that allows 
# visualization, analysis, and customization of various pathway maps.
# By mapping our transcriptomic data onto these pathways, we can explore 
# potential metabolic changes and regulatory mechanisms.

library(biomartr) # For querying biological databases
library(tidyverse) # For data manipulation and visualization
library(biomaRt) # For querying Ensembl databases

diff_genes <- read_csv(file = "tutorial/differential_genes.csv")

# Define the attributes to retrieve from BioMart: TAIR symbol, UniProt ID, and Entrez Gene ID
attributes_to_retrieve <- c("tair_symbol", "uniprotswissprot", "entrezgene_id")

# query BioMart to retrieve annotations for differentially expressed genes
diff_arabidopsis_genes_annotated <- biomartr::biomart(genes = diff_genes$genes,
                                                      mart       = "plants_mart",                 
                                                      dataset    = "athaliana_eg_gene",           
                                                      attributes = attributes_to_retrieve,        
                                                      filters =  "ensembl_gene_id" ) 

# filter genes with UniProt IDs, ensure unique entries, and format IDs for iPath
diff_arabidopsis_genes_annotated %>% 
  filter(uniprotswissprot != "") %>%                                       # to remove genes with no matching Uniprot entries
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>%   # to create an ID that iPath can use
  dplyr::select(id_for_ipath) %>%                                          # we keep only the relevant ID for further copy-pasting 
  write.table(., 
              file = "tutorial/diff_genes_swissprot.tsv", 
              row.names = FALSE, 
              quote = FALSE)

# calculate quantiles of log2 fold changes for upregulated genes
diff_genes %>% 
  filter(log2FoldChange > 0) %>% 
  with(.,quantile(log2FoldChange, c(0.5,0.75,0.9)))

# filter genes with log2 fold changes for upregulated genes
diff_genes_filtered = 
  diff_genes %>% 
  filter(log2FoldChange > quantile(log2FoldChange, 0.75)) 

# we query Ensembl again to retrieve the attributes
diff_arabidopsis_genes_annotated_2 <- biomartr::biomart(genes = diff_genes_filtered$genes,
                                                        mart       = "plants_mart",                 
                                                        dataset    = "athaliana_eg_gene",           
                                                        attributes = attributes_to_retrieve,        
                                                        filters =     "ensembl_gene_id" )  

# filter genes with UniProt IDs, ensure unique entries, and format IDs for iPath
diff_arabidopsis_genes_annotated_2 %>% 
  filter(uniprotswissprot != "") %>% 
  unique() %>% 
  mutate(id_for_ipath = paste("UNIPROT",uniprotswissprot,sep = ":")) %>% 
  dplyr::select(id_for_ipath) %>% 
  write.table(., file = "tutorial/diff_genes_swissprot_2.tsv", row.names = FALSE, quote = FALSE)
