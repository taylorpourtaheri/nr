% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_grid.R
\name{parameter_grid}
\alias{parameter_grid}
\title{Implement pipeline workflow over a grid of parameters.}
\usage{
parameter_grid(
  deg,
  target,
  grid,
  pipeline_method = NULL,
  connected_filter = NULL,
  n_cores = 1,
  ...
)
}
\arguments{
\item{deg}{Object of class '\code{dataframe}'. Differential gene expression
analysis results including log fold-changes and p-values,
i.e. the output of \code{topTable()} from limma or \code{results()} from DESeq2.}

\item{target}{A string. The gene symbol associated with the causal gene.}

\item{grid}{Object of class '\code{dataframe}', where each row corresponds to
a combination of parameter values. i.e. the output of \code{expand.grid()}.}

\item{n_cores}{Numeric. Number of cores to run on. Default value is 1.}

\item{pipeline}{String. The workflow pipeline to be used. The value can be either
'centrality' or 'propagation', corresponding the \code{centrality_pipeline()}
and \code{propagation_pipeline()} functions.}

\item{method}{A string specifying the network analysis method for the pipeline.
See ?\code{centrality_pipeline()} or ?\code{propagation_pipeline()} for more.}
}
\value{
Object of class '\code{dataframe}', with columns corresponding to the
parameter combination implemented and the
}
\description{
Given a dataframe of differential gene expression results,
this function will implement the specified pipeline with each combination of
parameters and return a dataframe of performance evaluation metrics.
}
