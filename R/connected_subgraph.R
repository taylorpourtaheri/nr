#' @title Extract the connected subcomponents of a graph
#' @description Extract the connected components of a graph based on a size threshold
#' @param G Graph of class '\code{igraph}'
#' @param thresh Numeric. The threshold size for a connected component. If
#' \code{NA} then only the giant (i.e., the largest) component is returned. All ties are returned.
#' @param method Either "less" or "greater". The latter returns all connected
#' components with greater than or equal to \code{thresh} nodes. Ignored if \code{thresh}
#' is \code{NA}.
#' @export
connected_subgraph <- function(graph, thresh=NA, method="greater") {

    methods <- c("less", "greater")
    if (!(method %in% methods)) {
        stop("Error: method must be either 'less' or 'greater'")
    }

    cl <- igraph::clusters(graph)
    if (is.na(thresh)) {
        return(igraph::induced_subgraph(graph, which(cl$membership %in% which(cl$csize == max(cl$csize))),
                                        impl = "copy_and_delete"))
    } else if (method == "greater") {
        return(igraph::induced_subgraph(graph, which(cl$membership %in% which(cl$csize >= thresh)),
                                        impl = "copy_and_delete"))
    } else {
        return(igraph::induced_subgraph(graph, which(cl$membership %in% which(cl$csize <= thresh)),
                                        impl = "copy_and_delete"))
    }
}
