# load relevant packages
library(STRINGdb)
library(devtools)
library(magrittr)
library(igraph)
library(ggplot2)
library(dplyr)
library(glue)

# load all package functions
load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# list of subtypes - fix to match true gene names (need to confirm in publication)
names(de_string)
causal_genes <- c('MYC', 'BCAT2', 'E2F3', 'HRAS', 'SRC')

# call wrapper
subtype_results <- purrr::map2(de_string, causal_genes,
    ~centrality_pipeline(deg = .x,
                         edge_conf_score_min = 950,
                         logFC_min = 1.0,
                         pvalue_max = 0.05,
                         method = 'betweenness',
                         causal_gene_symbol = .y,
                         export_network = FALSE,
                         sim_method = 'jaccard',
                         n_sim = 9999,
                         weighted = TRUE,
                         connected_filter = TRUE))
saveRDS(subtype_results, 'results/all_subtypes/subtype_results_full_1.RDS')


View(subtype_results)

gene_results <- purrr::map(subtype_results, ~.[['top_genes']])
purrr::map(gene_results, ~nrow(.))

performance_results <- purrr::map(subtype_results, ~.[['performance']])
performance_results %<>% bind_rows(.id = 'subtype')
saveRDS(performance_results, 'results/all_subtypes/subtype_results_performance_1.RDS')
openxlsx::write.xlsx(performance_results, 'results/all_subtypes/subtype_results_performance_1.xlsx')

# Observations from DEA results:

# MYC DEA results show good results for MYC (high log2FC/low p)

# BCAT DEA results show negative log2FC/high p for both BCAT1 and BCAT2
    # ranked in the 3000-4500 range by adjusted p-value
    # end up being removed in the network filtering step
    # only 11 genes in the final filtered giant component network
    # also only has 24 non-zero similarity scores in the entire full ppi network??
        # how would this change when adjusting the edge_conf_score_min?
plot(subtype_results$BCAT$network)

# E2F3 DEA results show a negative log2FC/high p for E2F3,
    # but log2FC of 2 with low p for E2F7

# RAS DEA results show several RAS isoforms are significant (HRAS listed in paper)
#  (high log2FC/low p)

# SRC DEA results show good results for SRC (high log2FC/low p)


network <- subtype_results$MYC$network

ggplot(ggnetwork(network), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_nodes() +
    geom_edges()

