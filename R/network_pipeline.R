network_pipeline <- function(deg,
                             edge_conf_score_min, logFC_min, pvalue_min,
                             method = 'centrality', causal_gene_symbol,
                             export_network = FALSE, sim_method = 'jaccard', n_sim = 9999){
    # internal check
    print(causal_gene_symbol)

    # list to store results
    final_results <- c()

    # generate protein association network
    string_db <- STRINGdb::STRINGdb$new(version="10",
                                        species=9606,
                                        score_threshold=edge_conf_score_min)
    ppi <- string_db$get_graph()

    # map DEA results onto ppi network
    ppi_painted <- df_to_vert_attr(graph=ppi, df=deg, common="STRING_id",
                                   attr_name = c("Symbol", "ID", "logFC", "AveExpr",
                                                 "t", "P.Value", "adj.P.Val", "B"))

    # subset the graph to only include nodes that meet thresholds
    ppi_painted_filt <- attribute_filter(ppi_painted,
                                         abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_min)

    # select the connected subgraph
    ppi_painted_filt_giant <- connected_subgraph(ppi_painted_filt)

    # calculate centrality
    if (method == 'centrality'){

        ppi_painted_filt_giant <- calc_centrality(ppi_painted_filt_giant,
                                                  bt=TRUE, len = -1)
    }

    # write final graph
    if (export_network){
        igraph::write_graph(ppi_painted_filt_giant,
                            file=glue("data/network_result_{edge_conf_score_min}.graphml"),
                            format = "graphml")
    }

    # find the STRING ID for the causal gene
    xref <- data.frame(symbol = causal_gene_symbol)
    xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

    # calculate similarity of each node and slice out the causal gene row
    sim <- igraph::similarity(ppi, method = sim_method)
    index <- which(V(ppi)$name == xref$STRING_id)
    causal_sim <- sim[index,]

    # make scores a named vector
    names(causal_sim) <- V(ppi)$name

    # get the scores associated with the subnetwork
    pred_scores <- causal_sim[V(ppi_painted_filt_giant)$name]
    mean_pred_score <- mean(pred_scores)

    # create a key mapping STRING id to gene symbol for all genes
    key <- data.frame(symbol = deg$Symbol)
    key <- string_db$map(key, "symbol", removeUnmappedRows=T, quiet=T)

    # select top genes and annotate for readability
    top_genes <- sort(pred_scores, decreasing = TRUE)
    top_genes_df <- data.frame(score = top_genes,
                               STRING_id  = names(top_genes))
    top_genes_df <- left_join(top_genes_df, key)
    rownames(top_genes_df) <- NULL

    # estimate uncertainty with a random draw of the full ppi graph
    n_draws <- length(V(ppi_painted_filt_giant))
    samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))
    sample_means <- sapply(samples, mean)

    # calculate p
    score_pval <- sum(sample_means > mean_pred_score) / n_sim

    # save results
    final_results[['network']] <- ppi_painted_filt_giant
    final_results[['top_genes']] <- top_genes_df
    final_results[['mean_score']] <- mean_pred_score
    final_results[['pvalue']] <- score_pval

    return(final_results)

}
