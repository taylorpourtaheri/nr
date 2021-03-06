% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_centrality.R
\name{calc_centrality}
\alias{calc_centrality}
\title{Network centrality}
\usage{
calc_centrality(graph, method = "strength", bt = FALSE, len = 3)
}
\arguments{
\item{graph}{Graph of class '\code{igraph}'}

\item{method}{String indicating centrality method to use (i.e. 'strength', 'degree')}

\item{bt}{Logical. Calculate betweenness centrality.}

\item{len}{Path length to estimate betweenness centrality. Ignored if \code{bt=FALSE}.}
}
\description{
Calculate various measures of network centrality
}
