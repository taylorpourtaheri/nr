#' @title Calculate network propagation by random walk
#' @description Calculate various measures of random walk
#' @import dnet, igraph
#' @param graph Graph of class '\code{igraph}'. Must include node
#' attributes '\code{name}' and '\code{seed}'.
#' @param method String. Method for random walk. Passed to \code{dRWR()}
#' from the \code{dnet} package.
#' @return Graph of class '\code{igraph}' with new attribute, '\code{propagation_score}'.
#' @export
calc_random_walk <- function(graph, restart_value = 0.75, norm_method = 'column') {

    dnet_seeds <- as.matrix(igraph::vertex.attributes(graph)$seed)
    rownames(dnet_seeds) <- igraph::vertex.attributes(graph)$name

    affinity_scores <- dnet::dRWR(g = graph,
                                     setSeeds = dnet_seeds,
                                     restart = restart_value,
                                     normalise = norm_method)

    affinity_scores <- as.numeric(as.matrix(affinity_scores))

    graph <- igraph::set_vertex_attr(graph,
                                     name = 'propagation_score',
                                     value = affinity_scores)

    return(graph)
}






