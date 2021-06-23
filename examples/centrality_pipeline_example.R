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

deg = myc_de
edge_conf_score_min = 950
logFC_min = 2.0
pvalue_max = 0.05
method = 'betweenness'
causal_gene_symbol = 'MYC'
export_network = FALSE
sim_method = 'jaccard'
n_sim = 9999
weighted = TRUE
ppi = NULL
string_db = NULL


# call wrapper
results <- centrality_pipeline(deg = myc_de,
                            edge_conf_score_min = 950,
                            logFC_min = 2.0,
                            pvalue_max = 0.05,
                            method = 'betweenness',
                            causal_gene_symbol = 'MYC',
                            export_network = FALSE,
                            sim_method = 'jaccard',
                            n_sim = 9999,
                            weighted = TRUE,
                            connected_filter = FALSE)

View(results)

# annotate with gene symbols
top_genes <- results$top_genes


# plotting
set.seed(4)
plot_graph(results[['network']], method = 'weighted_score', gene_list = c('MYC'))
