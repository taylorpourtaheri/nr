#' @title Parameter Grid Heatmap
#' @description Create a heatmap of parameter grid results.
#' @param results  Object of class '\code{dataframe}'. Parameter grid results.
#' @param plot_margin Object of class '\code{list}'. Allows adjustment of plot margins.
#' @return Object of class '\code{image}' showing z-score statistics derived from 
#' parameter_grid results
#' @export

param_heatmap <- function(results, plot_margin = c(8,5)){

  
  results_clip <- results %>% select(pipeline_method, logFC, Adj.P, z_score)
  
  results_mat <- results_clip %>%
    # Convert long-form to wide-form
    pivot_wider(names_from = c(pipeline_method, logFC), values_from = z_score) %>%
    column_to_rownames('Adj.P') %>%
    as.matrix
  
  results_mat<-results_mat[,order(colnames(results_mat))]
  
  
  return(heatmap(results_mat, 
                    Colv = NA, 
                    Rowv = NA, 
                    scale = "column",
                    margins = plot_margin,
                    main = "Z-score Heatmap",
                    xlab = "Method @ logFC",
                    ylab = "Adj-P"))
  }
  