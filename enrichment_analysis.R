# Functional enrichment analysis is a crucial step in interpreting the biological 
# significance of differentially expressed genes. This process helps identify biological
# pathways, processes, or functions that are overrepresented in a set of genes.
# This example is about analyzing Arabidopsis leaf genes that are differentially 
# expressed upon infection by Pseudomonas syringae DC3000 after 7 days.
# Over-Representation Analysis (ORA) is performed to identify enriched Gene 
# Ontology (GO) terms, which can provide insights into the biological processes 
# involved in the response to infection.

library(DESeq2) # For differential expression analysis
library(tidyverse) # For data manipulation and visualization
library(apeglm) # For log fold change shrinkage
library(biomartr) # For querying biological databases
library(clusterProfiler) # For functional enrichment analysis
library(enrichplot) # For visualizing enrichment results
suppressPackageStartupMessages(library("org.At.tair.db")) # A. thaliana genome database
library(biomaRt)  # only use to remove cache bug
library(splitstackshape) # for data manipulation

# read the differentially expressed genes data
diff_genes <- read_csv(file = "tutorial/differential_genes.csv")

# calculate quantiles of log2 fold changes for upregulated genes
diff_genes %>% 
  filter(log2FoldChange > 0) %>% 
  with(.,quantile(log2FoldChange, c(0.5,0.75,0.9)))

# filter genes with log2 fold changes above the 75th percentile and save to a file for AgriGO analysis
diff_genes %>% 
  filter(log2FoldChange > quantile(log2FoldChange, c(0.75))) %>% # keeping fold changes above the 75th percentile
  dplyr::select(genes) %>% 
  write.table(., file = "diff_genes_for_agrigo.tsv", row.names = FALSE, quote = FALSE)

# filter genes with log2 fold changes above the 75th percentile and save to a file for AgriGO PAGE analysis
diff_genes %>% 
  filter(log2FoldChange > quantile(log2FoldChange, c(0.75))) %>% 
  dplyr::select(genes, log2FoldChange) %>% 
  write.table(., file = "diff_genes_for_agrigo_page.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# check available datasets for Arabidopsis thaliana in BioMart
biomartr::organismBM(organism = "Arabidopsis thaliana")

# retrieve available attributes for A. thaliana in BioMart
arabido_attributes = 
  biomartr::organismAttributes("Arabidopsis thaliana") %>% 
  filter(dataset == "athaliana_eg_gene")
arabido_attributes

# specify attributes to retrieve: TAIR symbol and Entrez Gene ID
attributes_to_retrieve = c("tair_symbol", "entrezgene_id")


# query BioMart to retrieve TAIR symbols and Entrez Gene IDs for differentially 
# expressed genes
result_BM <- biomartr::biomart( genes      = diff_genes$genes, # genes were retrieved using biomartr::getGenome()
                                mart       = "plants_mart", # marts were selected with biomartr::getMarts()
                                dataset    = "athaliana_eg_gene", # datasets were selected with biomartr::getDatasets()
                                attributes = attributes_to_retrieve, # attributes were selected with biomartr::getAttributes()
                                filters =   "ensembl_gene_id" )# query key
head(result_BM)  

# read all Arabidopsis genes from the raw counts file to build the background/universe
all_arabidopsis_genes <- read.delim("tutorial/raw_counts.csv", sep = ",", header = T, stringsAsFactors = F)[,1] # directly selects the gene column

# we want the correspondence of TAIR/Ensembl symbols with NCBI Entrez gene ids
attributes_to_retrieve = c("tair_symbol", "uniprotswissprot","entrezgene_id")

# Query BioMart to retrieve annotations for all Arabidopsis genes
all_arabidopsis_genes_annotated <- biomartr::biomart(genes = all_arabidopsis_genes,
                                                     mart       = "plants_mart",                 
                                                     dataset    = "athaliana_eg_gene",           
                                                     attributes = attributes_to_retrieve,        
                                                     filters =  "ensembl_gene_id" )

# for compatibility with enrichGO universe
# genes in the universe need to be characters and not integers (Entrez gene id)
all_arabidopsis_genes_annotated$entrezgene_id = as.character(
  all_arabidopsis_genes_annotated$entrezgene_id) 

# retrieving NCBI Entrez gene id for our genes called differential
diff_arabidopsis_genes_annotated <- biomartr::biomart(genes = diff_genes$genes,
                                                      mart       = "plants_mart",                 
                                                      dataset    = "athaliana_eg_gene",           
                                                      attributes = attributes_to_retrieve,        
                                                      filters =  "ensembl_gene_id" ) 


# perform Over-Representation Analysis (ORA) for Gene Ontology Biological Process class
ora_analysis_bp <- enrichGO(gene = diff_arabidopsis_genes_annotated$entrezgene_id, # differentially expressed genes
                            universe = all_arabidopsis_genes_annotated$entrezgene_id, # background set of genes
                            OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                            keyType = "ENTREZID", # type of gene identifiers
                            ont = "BP", # either "BP" (biological process), "CC" (celullar component) or "MF" (molecular function),
                            pAdjustMethod = "BH", # Benjamini-Hochberg method for adjusting p-values
                            qvalueCutoff = 0.05, # cutoff for adjusted p-values
                            readable = TRUE, # return readable gene IDs
                            pool = FALSE # do not pool results
                            ) 

# write ORA-tsv file
write_tsv(x = as.data.frame(ora_analysis_bp_simplified@result), 
            path = "tutorial/results/go_results.tsv")

# Simplify the ORA results to remove redundant GO terms
ora_analysis_bp_simplified <- clusterProfiler::simplify(ora_analysis_bp)

# dotplot to visualize the simplified ORA results
dotplot(ora_analysis_bp_simplified)

# calculare pairwise term similarity for the simplified ORA results
ora_analysis_bp <- pairwise_termsim(ora_analysis_bp_simplified, method = "JC")

# create an emap plot to visualize the similarity between GO terms
emapplot(ora_analysis_bp, color = "qvalue")

# InterProScan

# interpro <- read.delim("tutorial/results/interproscan_results_web.tsv", header = F, stringsAsFactors = F, row.names = NULL)
# 
# new_colnames <- c("protein_id",
#                   "md5",
#                   "seq_len",
#                   "analysis",
#                   "signature_accession",
#                   "signature_description",
#                   "start",
#                   "stop",
#                   "score",
#                   "status",
#                   "date",
#                   "interpro_accession",
#                   "interpro_description",
#                   "go")
# colnames(interpro) <- new_colnames
# 
# interpro_go <- 
#   interpro %>% 
#   dplyr::select(protein_id, go) %>% 
#   dplyr::filter(go != "-") %>% 
#   dplyr::filter(go != "")
# tail(interpro_go, n = 10)
# 
# # split by pipes in go terms
# 
# splitted_interpro_go <- cSplit(indt = interpro_go, splitCols = "go", sep = "|", direction = "long")
# dedup_go  <- splitted_interpro_go %>% distinct() %>% 
#   dplyr::rename("genes" = "protein_id") %>% # change column name
#   mutate(genes = substr(genes, 1, 9)) # remove gene versioning
# tail(dedup_go)
# 
# write.csv(x = dedup_go, 
#           file = "tutorial/results/gene_ontologies_all_genes_web.csv", 
#           row.names = F)
# 
# diff_genes <- read_csv(file = "tutorial/differential_genes.csv")
# diff_genes_go <- inner_join(x = dedup_go, y = diff_genes)
# head(diff_genes_go)
# 
# diff_genes_go %>% 
#   dplyr::select(genes, go) %>% 
#   write.csv(file = "tutorial/results/gene_ontologies_diff_genes_web.csv", row.names = F)
