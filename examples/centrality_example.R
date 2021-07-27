# load relevant packages
# library(STRINGdb)
# library(ggnetwork)
# library(ggplot2)
# library(igraph)
# library(dplyr)
# library(glue)

# load all package functions
devtools::load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# read gene ID xref
id_xref <- readRDS('data/biomart_xreference_ppi_genes.RDS')

# read similarity score list
sim <- readRDS('data/string_ppi_v11_jacc_sim_list_dense.RDS')

# read ppi
ppi <- readRDS('data/string_ppi_v11.RDS')

# define globals
deg <- myc_de
edge_conf_score_min <- 950
logFC_min <- 1.5
pvalue_max <- 0.05
causal_gene_symbol <- 'MYC'
method <- 'betweenness'
connected_filter <- TRUE
export_network <- FALSE
n_sim <- 9999
weighted <- TRUE


# network generation ------------------------------------------------------

# generate protein association network
#
# string_db <- STRINGdb::STRINGdb$new(version="11",
#                                    species=9606,
#                                    score_threshold=edge_conf_score_min)
# ppi <- readRDS('data/string_ppi_v11.RDS')

# map DEA results onto ppi network
ppi_painted <- df_to_vert_attr(graph=ppi, df=deg, common="STRING_id",
                               attr_name = c("Symbol", "ID", "logFC", "AveExpr",
                                             "t", "P.Value", "adj.P.Val", "B"))

# subset the graph to only include nodes that meet thresholds
ppi_painted_filt <- attribute_filter(ppi_painted,
                                     abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_max)

# select the connected subgraph
ppi_painted_filt_giant <- connected_subgraph(ppi_painted_filt)

# calculate centrality
ppi_painted_filt_giant <- calc_centrality(ppi_painted_filt_giant, method = method, bt=T, len = -1)
# add argument 'method' that calls this function if value = 'centrality'

# write final graph
# igraph::write_graph(ppi_painted_filt_giant,
#                     file=glue("data/MYC_DE_network_example_{edge_conf_score_min}.graphml"),
#                     format = "graphml")
# ^ this might be a good place for a logical argument 'export_graph'

# network scoring ---------------------------------------------------------

# generate network scores
scoring_output <- structural_sim(network = ppi_painted_filt_giant,
                                 ppi = ppi_painted,
                                 id_xref = id_xref,
                                 # string_db = string_db,
                                 method = method,
                                 sim_method = 'jaccard',
                                 causal_gene_symbol = causal_gene_symbol,
                                 weighted = weighted)

# evaluate scoring
performance <- evaluate_performance(target = causal_gene_symbol,
                                            network_df = scoring_output$network_df,
                                            causal_sim = scoring_output$causal_sim,
                                            method = method,
                                            n_sim = n_sim,
                                            weighted = weighted)

performance_results <- performance[['performance_results']]
simulation_scores <- performance[['simulation_scores']]


# save results
final_results[['network']] <- scoring_output$network
final_results[['top_genes']] <- scoring_output$network_df
final_results[['performance']] <- performance_results
final_results[['simulation_scores']] <- simulation_scores

lobstr::mem_used()

# plot results
network <- final_results[['network']]

# plotting
set.seed(4)
plot_graph(network, method = 'weighted_score', gene_list = c('MYC'))

