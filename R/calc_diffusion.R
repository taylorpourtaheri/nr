#' @title Network propagation
#' @description Calculate various measures of network propagation
#' @import diffuStats
#' @param graph Graph of class '\code{igraph}'. Must include node
#' attributes '\code{name}' and '\code{seed}'.
#' @param method String. Kernel for diffusion. Passed to \code{diffusion()}
#' from the \code{diffuStats} package.
#' @return Graph of class '\code{igraph}' with new attribute, '\code{diffusion_score}'.
#' @export
calc_diffusion <- function(graph, method = 'raw') {

    # simulate diffusion
    scores <- as_data_frame(vertex_attr(graph))$seed
    names(scores) <- as_data_frame(vertex_attr(graph))$name
    diffusion_scores <- diffuStats::diffuse(graph = graph,
                                            scores = scores,
                                            method = method)
    graph <- igraph::set_vertex_attr(graph,
                                     name = 'diffusion_score',
                                     value = diffusion_scores)

    return(graph)
}
