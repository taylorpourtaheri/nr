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
#' @param method String. The method passed to the  diffusion() function
#' from the diffuStats() package. One of the following methods:
#' \itemize{
#' \item raw
#' \item ml
#' \item gm
#' \item ber_s
#' \item mc
#' \item ber_p
#' \item z
#' }
#' @param min_diff_score Numeric. The diffusion score threshold for filtering
#' the final network.
#' @param causal_gene_symbol String. The gene symbol associated with the
#' causal gene.
#' @param export_network Logical. If TRUE, the network object will be returned.
#' @param export_dir String. If export_network = TRUE, the network object will be
#' saved using the file path provided in export_dir.
#' @param sim_method String. The method for calculating the similarity between
#' each gene in the final network and the causal gene. One of the following:
#' #' \itemize{
#' \item jaccard
#' \item dice
#' \item invlogweighted }
#' @param n_sim Numeric. The number of simulations passed to \code{evaluate_performance()}.
#' @return A list of length 4:
#' \describe{
#'   \item{network}{Object of class '\code{igraph}'. The network of important genes.}
#'   \item{top_genes}{An annotated dataframe of genes ranked by diffusion score.}
#'   \item{mean_score}{The mean score of the network, which is equal to the average
#'   similarity score between each node in the network and the causal gene.}
#'   \item{pvalue}{The p-value associated with the \code{mean_score}.}
#' }
#' @export
propagation_pipeline <- function(deg, ppi = NULL, string_db = NULL,
                             edge_conf_score_min, logFC_min, pvalue_max,
                             method = 'raw', min_diff_score = 0.15, causal_gene_symbol,
                             export_network = FALSE, export_dir = NULL,
                             sim_method = 'jaccard', n_sim = 9999, weighted = FALSE,
                             ...){

    # internal check
    print(causal_gene_symbol)

    # list to store results
    final_results <- c()

    if (is.null(ppi) & is.null(string_db)){

        # generate protein association network
        string_db <- STRINGdb::STRINGdb$new(version="11",
                                            species=9606,
                                            score_threshold=edge_conf_score_min)
        ppi <- string_db$get_graph()
    }

    # filter to only include genes that are also in the DEG results
        # removed after conversation with Steve
    # ppi2 <- attribute_filter(ppi, name %in% deg$STRING_id)

    # map DEA results onto ppi network
    ppi_painted <- df_to_vert_attr(graph=ppi, df=deg, common="STRING_id",
                                   attr_name = c("Symbol", "ID", "logFC", "AveExpr",
                                                 "t", "P.Value", "adj.P.Val", "B"))

    # seed the graph based on the fold-change and pvalue thresholds
    ppi_painted_filt <- attribute_seed(ppi_painted,
                                       abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_max)

    # select the connected subgraph
    ppi_painted_filt_giant <- connected_subgraph(ppi_painted_filt)

    # calculate diffusion scores

    # diffusion with diffuStats:
    diffuStats_methods <- c('raw', 'ml', 'gm', 'ber_s', 'mc', 'ber_p', 'z')

    if (method %in% diffuStats_methods){

        ppi_painted_filt_giant <- calc_diffusion(graph = ppi_painted_filt_giant,
                                                 method = method)
    }

    # random walk with dnet:
    if (method == 'random_walk'){

        ppi_painted_filt_giant <- calc_random_walk(graph = ppi_painted_filt_giant,
                                                   ...)
    }

    # # filter network to only include scores above threshold
    # final_network <- attribute_filter(ppi_painted_filt_giant,
    #                                   propagation_score > min_diff_score)
    # remove duplicated edges
    final_network <- ppi_painted_filt_giant

    final_network_simple <- igraph::simplify(final_network)

    # write final graph
    if (export_network){
        igraph::write_graph(final_network_simple,
                            file = export_dir,
                            format = "graphml")
    }

    # generate network scores
    scoring_output <- structural_sim(network = final_network_simple,
                                     string_db = string_db,
                                     ppi = ppi,
                                     method = 'propagation_score',
                                     sim_method = sim_method,
                                     causal_gene_symbol = causal_gene_symbol,
                                     weighted = weighted)

    # evaluate scoring
    performance_results <- evaluate_performance(network_df = scoring_output$network_df,
                                                causal_sim = scoring_output$causal_sim,
                                                method = 'propagation_score',
                                                n_sim = n_sim,
                                                weighted = weighted)

    # save results
    final_results[['network']] <- scoring_output$network
    final_results[['top_genes']] <- scoring_output$network_df
    final_results[['performance']] <- performance_results


    return(final_results)

}
