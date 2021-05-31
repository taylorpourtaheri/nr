devtools::load_all()

set.seed(4)

plots <- c()

# import data and define parameters ---------------------------------------

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
method <- 'betweenness'
final_results <- c()
export_network <- FALSE
n_sim <- 9999

gene_list <- c('MYC')

###### for visualization purposes only
deg <- dplyr::filter(deg, abs(logFC) < 10)


# generate ppi ------------------------------------------------------------

# generate protein association network
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=edge_conf_score_min)
ppi <- string_db$get_graph()

####
ggn1 <- ggnetwork::ggnetwork(ppi)
as.name(method)
plots[[1]] <- ggplot2::ggplot(ggn1, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = 'grey') +
    geom_nodes(color = 'black', size = 6, alpha = 0.5) +
    # geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    # geom_nodelabel_repel(data=subset(ggn, Symbol %in% gene_list), aes(label=Symbol)) +
    scale_color_gradient2(low = 'blue', high = 'red') +
    # scale_size_continuous(range = c(5, 25)) +
    theme_blank()
####

# map DEA results onto ppi network
ppi_painted <- df_to_vert_attr(graph=ppi, df=deg, common="STRING_id",
                               attr_name = c("Symbol", "ID", "logFC", "AveExpr",
                                             "t", "P.Value", "adj.P.Val", "B"))

####
ggn2 <- ggnetwork::ggnetwork(ppi_painted)
as.name(method)
plots[[2]] <- ggplot2::ggplot(ggn2, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = 'grey') +
    geom_nodes(aes_string(color = 'logFC'), size = 6, alpha = 0.75) +
    # geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    # geom_nodelabel_repel(data=subset(ggn, Symbol %in% gene_list), aes(label=Symbol)) +
    scale_color_gradient2(low = 'blue', high = 'red') +
    # scale_size_continuous(range = c(5, 25)) +
    theme_blank()
####

# subset the graph to only include nodes that meet thresholds
ppi_painted_filt <- attribute_filter(ppi_painted,
                                     abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_max)

# select the connected subgraph
ppi_painted_filt_giant <- connected_subgraph(ppi_painted_filt)

####
ggn3 <- ggnetwork::ggnetwork(ppi_painted_filt_giant)
as.name(method)
plots[[3]] <- ggplot2::ggplot(ggn3, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = 'grey') +
    geom_nodes(aes_string(color = 'logFC'), size = 6, alpha = 0.75) +
    # geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    # geom_nodelabel_repel(data=subset(ggn3, Symbol %in% gene_list), aes(label=Symbol)) +
    scale_color_gradient2(low = 'blue', high = 'red') +
    # scale_size_continuous(range = c(5, 25)) +
    theme_blank()
####

# calculate centrality
ppi_painted_filt_giant <- calc_centrality(ppi_painted_filt_giant, method = method, bt=T, len = -1)

####
ggn4 <- ggnetwork::ggnetwork(ppi_painted_filt_giant)
as.name(method)
plots[[4]] <- ggplot2::ggplot(ggn4, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes_string(color = 'logFC', size = method), alpha = 0.65) +
    # geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    # geom_nodelabel_repel(data=subset(ggn4, Symbol %in% gene_list), aes(label=Symbol)) +
    scale_color_gradient2(low = 'blue', high = 'red') +
    scale_size_continuous(range = c(5, 25)) +
    theme_blank()
####

# network and method evaluation -------------------------------------------

# generate network scores
scoring_output <- structural_sim(network = ppi_painted_filt_giant,
                                 ppi = ppi,
                                 string_db = string_db,
                                 method = 'betweenness',
                                 sim_method = 'jaccard',
                                 causal_gene_symbol = causal_gene_symbol,
                                 weighted = TRUE)


# evaluate scoring
performance_results <- evaluate_performance(network = scoring_output$network,
                                            network_df = scoring_output$network_df,
                                            causal_sim = scoring_output$causal_sim,
                                            n_sim = n_sim,
                                            method = 'betweenness',
                                            weighted = TRUE)


# examine results ---------------------------------------------------------

# plotting
# plots[[5]] <- plot_graph(scoring_output$network, method = 'weighted_score', gene_list = c('MYC'))


ggn5 <- ggnetwork::ggnetwork(scoring_output$network)
method <- 'weighted_score'
as.name(method)
plots[[5]] <- ggplot2::ggplot(ggn5, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges() +
    geom_nodes(aes_string(color = 'logFC', size = method), alpha = 0.65) +
    geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
    geom_nodelabel_repel(data=subset(ggn5, Symbol %in% gene_list), aes(label=Symbol)) +
    scale_color_gradient2(low = 'blue', high = 'red') +
    scale_size_continuous(range = c(5, 25)) +
    theme_blank()

assayr2::export_pngs(plots, 'presentations/ds6011/plots',
                    width = 8, height = 8)

