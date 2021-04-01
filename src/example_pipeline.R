# load relevant packages
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)
library(glue)

# load all package functions
load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# call wrapper
results <- network_pipeline(deg = myc_de,
                            edge_conf_score_min = 950,
                            logFC_min = 1.5,
                            pvalue_min = 0.05,
                            method = 'centrality',
                            causal_gene_symbol = 'MYC',
                            export_network = FALSE,
                            n_sim = 9999)

View(results)

# annotate with gene symbols
top_genes <- results$top_genes


