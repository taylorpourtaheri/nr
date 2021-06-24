#' @title Return seeded graph.
#' @import igraph
#' @description Use `attribute_seed()` to seed network nodes with diffusion points
#' based on given vertex attributes.
#' @param graph Graph of class '\code{igraph}'
#' @param attr_expression A logical expression defined in terms of the attribute names of \code{graph}.
#' only vertices where the expression evaluates to \code{TRUE} are kept.
#' @examples
#' # generate a graph
#' vertex_meta <- data.frame('name' = c('a','b','c','d'), value = c(1, 5, 10, 3))
#' edge_list <- data.frame('v1' = c('a', 'a', 'd', 'd'), 'v2' = c('b', 'c', 'a', 'c'))
#' graph <- igraph::graph_from_data_frame(d = edge_list, vertices = vertex_meta, directed = FALSE)
#'
#' # examine initial graph
#' igraph::vertex_attr(graph)
#'
#' # seed graph by node value
#' seeded_graph <- attribute_seed(graph = graph, value > 4)
#'
#' # examine seeded graph
#' igraph::vertex_attr(seeded_graph)
#'
#' # there is now a node attribute 'seed', which is equal to 1 for nodes b and c,
#' # and 0 for nodes a and d.
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
