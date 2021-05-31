#' @title Calculate network propagation by diffusion
#' @description Calculate various measures of diffusion
#' @import diffuStats
#' @param graph Graph of class '\code{igraph}'. Must include node
#' attributes '\code{name}' and '\code{seed}'.
#' @param method String. Kernel for diffusion. Passed to \code{diffusion()}
#' from the \code{diffuStats} package.
#' @return Graph of class '\code{igraph}' with new attribute, '\code{propagation_score}'.
#' @export
calc_diffusion <- function(graph, method = 'raw') {

    # simulate diffusion
    scores <- as.data.frame(igraph::vertex_attr(graph))$seed
    names(scores) <- as.data.frame(igraph::vertex_attr(graph))$name

    diffusion_scores <- diffuStats::diffuse(graph = graph,
                                            scores = scores,
                                            method = method)

    graph <- igraph::set_vertex_attr(graph,
                                     name = 'propagation_score',
                                     value = diffusion_scores)

    return(graph)
}
