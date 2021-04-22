
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
n_sim <- 9999

# new
kernel <- 'raw'
min_diff_score <- 0.15
grid <- FALSE

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

# simulate diffusion
scores <- as_data_frame(vertex_attr(ppi_painted_filt_giant))$seed
names(scores) <- as_data_frame(vertex_attr(ppi_painted_filt_giant))$name
diffusion_scores <- diffuStats::diffuse(graph = ppi_painted_filt_giant,
                                              scores = scores,
                                              method = kernel)
ppi_painted_filt_giant <- igraph::set_vertex_attr(ppi_painted_filt_giant,
                                    name = 'diffusion_score',
                                    value = diffusion_scores)

# filter network to only include scores above threshold
final_network <- attribute_filter(ppi_painted_filt_giant,
                                  diffusion_score > min_diff_score)
# remove duplicated edges
final_network_simple <- igraph::simplify(final_network)

View(igraph::as_data_frame(ppi_painted_filt_giant, what = 'vertices'))
# MYC is filtered out because of low diffusion score
View(igraph::as_data_frame(final_network_simple, what = 'vertices'))

# test plotting
set.seed(4)
ggn <- ggnetwork(final_network_simple)
ggplot(ggn, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(color = logFC, size = diffusion_score), alpha = 0.65) +
    geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    geom_nodelabel_repel(data=subset(ggn, Symbol == 'MYCN'), aes(label=Symbol)) +
    scale_color_gradient(low = 'blue', high = 'red') +
    scale_size_continuous(range = c(5, 25)) +
    theme_blank()

# dataframe of final graph results
network_df <- igraph::as_data_frame(final_network_simple, what = 'vertices')
network_df <- dplyr::arrange(network_df, -diffusion_score)
rownames(network_df) <- NULL

# network scoring ---------------------------------------------------------

# find the STRING ID for the causal gene
xref <- data.frame(symbol = causal_gene_symbol)
xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

# calculate similarity of each node and slice out the causal gene
sim <- igraph::similarity(ppi)
index <- which(igraph::V(ppi)$name == xref$STRING_id)
causal_sim <- sim[index,]

# make scores a named vector
names(causal_sim) <- igraph::V(ppi)$name

# get the scores associated with the subnetwork
pred_scores <- causal_sim[igraph::V(ppi_painted_filt_giant)$name]
mean_pred_score <- mean(pred_scores)

# create a key mapping STRING id to gene symbol for all genes
key <- data.frame(symbol = deg$Symbol)
key <- string_db$map(key, "symbol", removeUnmappedRows=T, quiet=T)

# select top genes and annotate for readability
top_genes <- sort(pred_scores, decreasing = TRUE)
top_genes_df <- data.frame(score = top_genes,
                           STRING_id  = names(top_genes))
top_genes_df <- dplyr::left_join(top_genes_df, key)
rownames(top_genes_df) <- NULL

# estimate uncertainty with a random draw of the full ppi graph
n_sim <- 9999
n_draws <- length(igraph::V(ppi_painted_filt_giant)) #151
samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))
sample_means <- sapply(samples, mean)

# calculate p
score_pval <- sum(sample_means > mean_pred_score) / n_sim

# save results
final_results[['network']] <- ppi_painted_filt_giant
final_results[['top_genes']] <- top_genes_df
final_results[['mean_score']] <- mean_pred_score
final_results[['pvalue']] <- score_pval



# plot results
network <- ppi_painted_filt_giant

set.seed(4)

# define plotting network
ggn <- ggnetwork(network)

ggplot(ggn, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes(color = logFC, size = betweenness), alpha = 0.65) +
    geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    geom_nodelabel_repel(data=subset(ggn, Symbol == 'MYC'), aes(label=Symbol)) +
    scale_color_gradient(low = 'blue', high = 'red') +
    scale_size_continuous(range = c(5, 25)) +
    theme_blank()
# ggsave('test.png', width = 15, height = 15)




