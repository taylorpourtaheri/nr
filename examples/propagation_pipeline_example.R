# load relevant packages
library(STRINGdb)
library(igraph)
library(ggplot2)
library(dplyr)
library(glue)
library(devtools)

# load all package functions
load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# call wrapper
results <- propagation_pipeline(deg = myc_de,
                            edge_conf_score_min = 950,
                            logFC_min = 1.5,
                            pvalue_max = 0.05,
                            method = 'raw',
                            min_diff_score = 0.15,
                            causal_gene_symbol = 'MYC',
                            export_network = FALSE,
                            sim_method = 'jaccard',
                            n_sim = 9999)

# plot output
set.seed(4)
plot_graph(results[['network']], method = 'diffusion_score', gene_list = c('MYCN'))

# ggsave('test.png', width = 15, height = 15)
