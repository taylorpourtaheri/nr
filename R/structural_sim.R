structural_sim <- function(network, string_db, ppi, method, sim_method, causal_gene_symbol, weighted = FALSE){

    output_list <- c()

    # dataframe of final graph results
    network_df <- igraph::as_data_frame(network, what = 'vertices')
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
    pred_scores <- causal_sim[igraph::V(network)$name]

    sim_scores <- data.frame('causal_similarity' = pred_scores,
                             'name' = igraph::V(network)$name)
    network_df <- dplyr::left_join(network_df, sim_scores)

    # add vertex attributes to network
    network <- igraph::set_vertex_attr(network, 'causal_similarity',
                                       value = pred_scores)


    if(weighted == TRUE){

        # create final weighted scoring metric and output table
        network_df$weighted_score <- network_df[,method]*network_df[ ,'causal_similarity']

        # add vertex attribute
        network <- igraph::set_vertex_attr(network, 'weighted_score',
                                           value = network_df$weighted_score)

        # sorted output
        network_df <- dplyr::arrange(network_df, -weighted_score)
    }

    output_list[['network']] <- network
    output_list[['causal_similarity']] <- causal_sim
    output_list[['network_df']] <- network_df


    return(output_list)
}
