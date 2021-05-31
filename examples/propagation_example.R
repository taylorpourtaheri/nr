
# load relevant packages
library(STRINGdb)
library(ggnetwork)
library(ggplot2)
library(igraph)
library(dplyr)
library(glue)

# load all package functions
devtools::load_all()

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# define globals
deg <- myc_de
edge_conf_score_min <- 950
logFC_min <- 1.5
pvalue_max <- 0.05
causal_gene_symbol <- 'MYC'
final_results <- c()
export_network <- FALSE
export_dir = NULL
n_sim <- 9999
method <- 'raw'
min_diff_score <- 0.15
sim_method = 'jaccard'
weighted = TRUE
ppi = NULL
string_db = NULL

# in the DEG results:
MYC_stringID <-'9606.ENSP00000479618'

# network generation ------------------------------------------------------

# generate protein association network
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=edge_conf_score_min)
ppi <- string_db$get_graph()

# filter to only include genes that are also in the DEG results - check this
ppi2 <- attribute_filter(ppi, name %in% deg$STRING_id)

# map DEA results onto ppi network
ppi_painted <- df_to_vert_attr(graph=ppi2, df=deg, common="STRING_id",
                               attr_name = c("Symbol", "ID", "logFC", "AveExpr",
                                             "t", "P.Value", "adj.P.Val", "B"))

# seed the graph based on the fold-change and pvalue thresholds
ppi_painted_filt <- attribute_seed(ppi_painted,
                                     abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_max)

# select the connected subgraph
ppi_painted_filt_giant <- connected_subgraph(ppi_painted_filt)

# calculate diffusion scores
ppi_painted_filt_giant <- calc_diffusion(graph = ppi_painted_filt_giant,
                                         method = method)

# filter network to only include scores above threshold
final_network <- attribute_filter(ppi_painted_filt_giant,
                                  diffusion_score > min_diff_score)
# remove duplicated edges
final_network_simple <- igraph::simplify(final_network)

# View(igraph::as_data_frame(ppi_painted_filt_giant, what = 'vertices'))
# MYC is filtered out because of low diffusion score
# View(igraph::as_data_frame(final_network_simple, what = 'vertices'))

# generate network scores
scoring_output <- structural_sim(network = final_network_simple,
                                 string_db = string_db,
                                 ppi = ppi,
                                 method = 'diffusion_score',
                                 sim_method = sim_method,
                                 causal_gene_symbol = causal_gene_symbol,
                                 weighted = weighted)

# evaluate scoring
performance_results <- evaluate_performance(network = scoring_output$network,
                                            network_df = scoring_output$network_df,
                                            causal_sim = scoring_output$causal_sim,
                                            method = 'diffusion_score',
                                            n_sim = n_sim,
                                            weighted = weighted)

# save results
final_results[['network']] <- scoring_output$network
final_results[['top_genes']] <- scoring_output$network_df
final_results[['performance']] <- performance_results


# plotting
set.seed(4)
network <- final_results[['network']]
plot_graph(network, method = 'diffusion_score', gene_list = c('MYCN'))


