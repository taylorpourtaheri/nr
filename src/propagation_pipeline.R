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
results <- centrality_pipeline(deg = myc_de,
                            edge_conf_score_min = 950,
                            logFC_min = 1.5,
                            pvalue_max = 0.05,
                            method = 'betweenness',
                            causal_gene_symbol = 'MYC',
                            export_network = FALSE,
                            n_sim = 9999)

View(results)

# annotate with gene symbols
top_genes <- results$top_genes


# define plotting network
ggn <- ggnetwork(results$network)

ggplot(ggn, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(color = logFC, size = betweenness), alpha = 0.65) +
    geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    geom_nodelabel_repel(data=subset(ggn, Symbol == 'MYC'), aes(label=Symbol)) +
    scale_color_gradient(low = 'blue', high = 'red') +
    scale_size_continuous(range = c(5, 25)) +
    theme_blank()
# ggsave('test.png', width = 15, height = 15)
