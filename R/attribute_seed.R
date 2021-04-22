#' @title Return seeded graph.
#' @description Use `attribute_seed()` to seed network nodes with diffusion points
#' based on given vertex attributes.
#' @param graph Graph of class '\code{igraph}'
#' @param attr_expression A logical expression defined in terms of the attribute names of \code{graph}.
#' only vertices where the expression evaluates to \code{TRUE} are kept.
#' @export
attribute_seed <- function(graph, attr_expression) {
    vdf <- as.data.frame(igraph::vertex_attr(graph), stringsAsFactors = FALSE)
    ind <- eval(substitute(attr_expression), vdf, parent.frame())
    seed <- as.numeric(ind)
    seed[is.na(seed)] <- 0 # convert na values to zeroes - check this
    g_seeded <- igraph::set_vertex_attr(graph,
                                        name = 'seed',
                                        value = seed)
    return(g_seeded)
}
