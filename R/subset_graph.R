#' @title Subset the loaded igraph based on edge confidence
#' @description Given an iGraph object from StringDb, reduce the
#' full graph to a subset based on edge confidence scores. This 
#' function will improve the performance of multiple calls to 
#' the iGraph object for the parameter optimization.
#' @param graph Graph of class '\code{igraph}'
#' @param score_threshold Integer value of the edge confidence'

subset_graph <- function(graph, score_threshold){
  sub_ppi <- igraph::subgraph.edges(graph, E(graph)[E(graph)$combined_score >= score_threshold], del=F)
  sub_ppi <- noderank::connected_subgraph(sub_ppi)
  
  return(sub_ppi)
}