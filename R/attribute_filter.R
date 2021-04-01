#' @title Return subgraph where attributes match conditions
#' @description Use `attribute_filter()` to find cases where conditions are true given vertex attributes.
#'
#' @param graph Graph of class '\code{igraph}'
#' @param attr_expression A logical expression defined in terms of the attribute names of \code{graph}.
#' only vertices where the expression evaluates to \code{TRUE} are kept.
#' @export
attribute_filter <- function(graph, attr_expression) {
    vdf <- as.data.frame(igraph::vertex_attr(graph), stringsAsFactors = FALSE)
    ind <- eval(substitute(attr_expression), vdf, parent.frame())
    g_filt <- igraph::induced_subgraph(graph, which(ind), impl = "copy_and_delete")
    return(g_filt)
}
