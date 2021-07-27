# sim - a named list of named similarity score vectors, each element corresponding
    # to a gene. A similarity matrix that has been split row-wise. Missing values
    # will be assigned a similarity score of zero.
structural_sim <- function(network, ppi, id_xref, method, sim_method, causal_gene_symbol,
                           sim = NULL, weighted = FALSE){

    output_list <- c()

    # dataframe of final graph results
    network_df <- igraph::as_data_frame(network, what = 'vertices')
    rownames(network_df) <- NULL

    # find the STRING ID for the causal gene
    xref <- dplyr::filter(id_xref, symbol == causal_gene_symbol)
    # xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

    # calculate similarity of each node and slice out the causal gene row
    if(is.null(sim)){

        sim <- igraph::similarity(ppi, method = sim_method)
        index <- which(igraph::V(ppi)$name == xref$STRING_id)
        causal_sim <- sim[index,]

        # make scores a named vector
        names(causal_sim) <- igraph::V(ppi)$name
    }

    # if list of dense similarity scores is given, select casual gene
    if (is.list(sim)){

        causal_vec <- sim[[xref$STRING_id]]

        add_zeros <- function(sim_vec, gene_id){
            if(!(gene_id %in% names(sim_vec))){
                0
            } else{
                sim_vec[gene_id]
            }
        }

        causal_sim <- unlist(purrr::map(igraph::V(ppi)$name, ~add_zeros(causal_vec, .)))
        names(causal_sim) <- igraph::V(ppi)$name
    }


    # get the scores associated with the subnetwork
    pred_scores <- causal_sim[igraph::V(network)$name]

    sim_scores <- data.frame('causal_similarity' = pred_scores,
                             'name' = igraph::V(network)$name)
    network_df <- dplyr::left_join(network_df, sim_scores)

    # sorted output
    network_df <- dplyr::arrange(network_df, -eval(as.name(method)))

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
