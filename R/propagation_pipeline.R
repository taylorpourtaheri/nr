#' @title Pipeline function to run complete noderank network propagation workflow.
#' @description Given a dataframe of differential gene expression results and a set
#' of parameters, '\code{centrality_pipeline()}' will generate network of the most
#' important genes, ranked by the specified centrality measure.
#' @param deg Object of class '\code{dataframe}'. Differential gene expression analysis
#' results including log fold-changes and p-values,
#' i.e. the output of \code{topTable()} from limma or \code{results()} from DESeq2.
#' @param edge_conf_score_min Numeric. A The minimum confidence score between edges in the
#' STRING protein-protein interaction network.
#' @param logFC_min Numeric. The minimum log2 fold-change value for a gene
#' to be included in the final network.
#' @param pvalue_max Numeric. The maximum p-value value for a gene
#' to be included in the final network.
#' @param method A string. One of the following centrality methods:
#' \itemize{
#' \item strength
#' \item degree
#' \item avg_strength
#' \item degree_frac
#' \item strength_scaled
#' \item avg_strength_scaled
#' \item evcent_w
#' \item evcent_uw
#' \item betweenness
#' }
#' @param causal_gene_symbol A string. The gene symbol associated with the
#' causal gene.
#' @param export_network If TRUE, the network object will be returned.
#' @param sim_method A string. The method for calculating the similarity between
#' each gene in the final network and the causal gene. One of the following:
#' #' \itemize{
#' \item jaccard
#' \item dice
#' \item invlogweighted
#' }
#' @param n_sim
#' @return A list of length 4:
#' \describe{
#'   \item{network}{Object of class '\code{igraph}'. The network of important genes.}
#'   \item{top_genes}{An annotated dataframe of genes ranked by centrality method.}
#'   \item{mean_score}{The mean score of the network, which is equal to the average
#'   similarity score between each node in the network and the causal gene.}
#'   \item{pvalue}{The p-value associated with the \code{mean_score}.}
#' }
#' @export
propagation_pipeline <- function(deg,
                             edge_conf_score_min, logFC_min, pvalue_max,
                             method = 'raw', min_diff_score = 0.15,
                             causal_gene_symbol, export_network = FALSE,
                             sim_method = 'jaccard', n_sim = 9999){
    # internal check
    print(causal_gene_symbol)

    # list to store results
    final_results <- c()

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

    # write final graph
    if (export_network){
        igraph::write_graph(final_network_simple,
                            file=glue::glue("data/network_result_{edge_conf_score_min}_{method}.graphml"),
                            format = "graphml")
    }

    # dataframe of final graph results
    network_df <- igraph::as_data_frame(final_network_simple, what = 'vertices')
    network_df <- dplyr::arrange(network_df, -diffusion_score)
    rownames(network_df) <- NULL

    # find the STRING ID for the causal gene
    xref <- data.frame(symbol = causal_gene_symbol)
    xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

    # calculate similarity of each node and slice out the causal gene row
    sim <- igraph::similarity(ppi, method = sim_method)
    index <- which(igraph::V(ppi)$name == xref$STRING_id)
    causal_sim <- sim[index,]

    # make scores a named vector
    names(causal_sim) <- igraph::V(ppi)$name

    # get the scores associated with the subnetwork
    pred_scores <- causal_sim[igraph::V(final_network_simple)$name]
    mean_pred_score <- mean(pred_scores)

    # estimate uncertainty with a random draw of the full ppi graph
    n_draws <- length(igraph::V(final_network_simple))
    samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))
    sample_means <- sapply(samples, mean)

    # calculate p
    score_pval <- sum(sample_means > mean_pred_score) / n_sim

    # save results
    final_results[['network']] <- final_network_simple
    final_results[['top_genes']] <- network_df
    final_results[['mean_score']] <- mean_pred_score
    final_results[['pvalue']] <- score_pval

    return(final_results)

}
