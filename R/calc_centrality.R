#' @title Network centrality
#' @description Calculate various measures of network centrality
#' @param graph Graph of class '\code{igraph}'
#' @param method Method of class '\code{string}'
#' @param bt Logical. Calculate betweenness centrality.
#' @param len Path length to estimate betweenness centrality. Ignored if \code{bt=FALSE}.
#' @export
calc_centrality <- function(graph, method = 'strength', bt = FALSE, len = 3) {

    if(method == 'strength'){
        igraph::V(graph)$strength <- igraph::graph.strength(graph, vids = igraph::V(graph), loops=F)
    }
    if (method == 'degree'){
        igraph::V(graph)$degree <- igraph::degree(graph, v = igraph::V(graph), loops=F)
    }
    if (method == 'avg_strength'){
        igraph::V(graph)$avg_strength <- igraph::V(graph)$strength / igraph::V(graph)$degree
    }
    if ("baseline_degree" %in% igraph::vertex_attr_names(graph)) {
        if(method == 'degree_frac'){
            igraph::V(graph)$degree_frac <- igraph::V(graph)$degree / igraph::V(graph)$baseline_degree
        }
        if (method == 'strength_scaled'){
            igraph::V(graph)$strength_scaled <- igraph::V(graph)$strength * igraph::V(graph)$degree_frac
        }
        if (method == 'avg_strength_scaled'){
            igraph::V(graph)$avg_strength_scaled <- igraph::V(graph)$avg_strength * igraph::V(graph)$degree_frac
        }
    }

    if (method == 'evcent_w'){
        igraph::V(graph)$evcent_w <- igraph::evcent(graph, scale = F)$vector
    }
    if (method == 'evcent_uw'){
        igraph::V(graph)$evcent_uw <- igraph::evcent(graph, scale = F, weights = NA)$vector
    }
    if (bt) {
        igraph::V(graph)$betweenness <- igraph::estimate_betweenness(graph, directed = F, cutoff = len)
    }
    return(graph)
}
