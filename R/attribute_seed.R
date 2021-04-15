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
    g_filt <- igraph::induced_subgraph(graph, which(ind), impl = "copy_and_delete")

    temp_list <- c()
    temp_list[[1]]<-vdf
    temp_list[[2]]<-ind

    return(temp_list)
}

a <- attribute_seed(ppi_painted,
                    abs(logFC) > log2(logFC_min) & adj.P.Val < pvalue_max)
