# load relevant packages
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)
library(glue)

# load all package functions
devtools::load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# string_db <- readRDS('data/string_db_v11.RDS')
ppi <- readRDS('data/string_ppi_v11.RDS')
sim <- readRDS('data/string_ppi_v11_jacc_sim_list_dense.RDS')
id_xref <- readRDS('data/biomart_xreference_ppi_genes.RDS')

# call wrapper
results <- centrality_pipeline(deg = myc_de,
                               # string_db = string_db,
                               ppi = ppi,
                               sim = sim,
                               id_xref = id_xref,
                            edge_conf_score_min = 950,
                            logFC_min = 1.5,
                            pvalue_max = 0.05,
                            method = 'betweenness',
                            causal_gene_symbol = 'MYC',
                            weighted = TRUE,
                            connected_filter = TRUE)

View(results)

# annotate with gene symbols
top_genes <- results$top_genes


# plotting
set.seed(4)
plot_graph(results[['network']], method = 'weighted_score', gene_list = c('MYC'))
