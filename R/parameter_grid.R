#' @title Implement pipeline workflow over a grid of parameters.
#' @description Given a dataframe of differential gene expression results,
#' this function will implement the specified pipeline with each combination of
#' parameters and return a dataframe of performance evaluation metrics.
#' @param deg  Object of class '\code{dataframe}'. Differential gene expression
#' analysis results including log fold-changes and p-values,
#' i.e. the output of \code{topTable()} from limma or \code{results()} from DESeq2.
#' @param target A string. The gene symbol associated with the causal gene.
#' @param grid Object of class '\code{dataframe}', where each row corresponds to
#' a combination of parameter values. i.e. the output of \code{expand.grid()}.
#' @param pipeline String. The workflow pipeline to be used. The value can be either
#' 'centrality' or 'propagation', corresponding the \code{centrality_pipeline()}
#' and \code{propagation_pipeline()} functions.
#' @param method A string specifying the network analysis method for the pipeline.
#' See ?\code{centrality_pipeline()} or ?\code{propagation_pipeline()} for more.
#' @param n_cores Numeric. Number of cores to run on. Default value is 1.
#' @return Object of class '\code{dataframe}', with columns corresponding to the
#' parameter combination implemented and the
#' @export
parameter_grid <- function(deg, target, grid, pipeline_method = NULL,
                           connected_filter = NULL, n_cores = 1, ...){

    if(n_cores > parallel::detectCores()){
        warn('n_cores exceeds number of cores; using maximum number of cores instead.')
        n_cores <- parallel::detectCores()
    }

    # import a list of STRINGdb objects and generate list of ppis
    thresholds <- unique(grid$threshold)
    threshold_string <- as.character(thresholds)

    string_list <- purrr::map(thresholds, ~ STRINGdb::STRINGdb$new(version = '11',
                                                                   species = 9606,
                                                                   score_threshold = .))
    ppi_list <- purrr::map(string_list, ~ .$get_graph())

    names(string_list) <- threshold_string
    names(ppi_list) <- threshold_string

    print('ppi_list imported')

    # convert grid to list of rows
    grid_list <- asplit(grid, 1)

    centrality_methods = c('strength','degree','avg_strength','degree_frac',
                           'strength_scaled','avg_strength_scaled','evcent_w',
                           'evcent_uw','betweenness')
    
    propagation_methods = c('raw', 'ml', 'gm', 'ber_s', 'mc', 'ber_p', 'z','random_walk')
    
    # run pipeline
    results_list <- parallel::mclapply(grid_list, FUN = function(x){

        print(x)

        pipeline_method <- x['pipeline_method']
        
        #not sure what the point of this was, seemed like it would not be updated at each iteration
        #if(is.null(pipeline_method)){
        #    pipeline_method <- x['pipeline_method']
        #}

        if(is.null(connected_filter)){
            connected_filter <- as.logical(x['connected_filter'])
        }

        # how to specify method in a way that won't allow centrality methods
        # to be passed to the propagation pipeline and vice versa?
        # if(is.null(method)){
        #     method <- x['method']
        # }

        if (pipeline_method %in% centrality_methods){
            pipeline_output <- centrality_pipeline(deg = deg,
                                                   ppi = ppi_list[[as.character(x['threshold'])]],
                                                   string_db = string_list[[as.character(x['threshold'])]],
                                                   logFC_min = as.numeric(x['logFC']),
                                                   pvalue_max = as.numeric(x['Adj.P']),
                                                   method = as.character(x['pipeline_method']),
                                                   causal_gene_symbol =  target,
						   connected_filter = connected_filter,
                                                   ...)
            pipeline <- 'centrality'
        }

        else if (pipeline_method %in% propagation_methods){
            pipeline_output <- propagation_pipeline(deg = deg,
                                                    ppi = ppi_list[[as.character(x['threshold'])]],
                                                    string_db = string_list[[as.character(x['threshold'])]],
                                                    logFC_min = as.numeric(x['logFC']),
                                                    pvalue_max = as.numeric(x['Adj.P']),
                                                    method = as.character(x['pipeline_method']),
                                                    min_diff_score = 0.15,
                                                    causal_gene_symbol =  target,
                                                    ...)
            pipeline <- 'propagation'
        }


        x_df <- data.frame('pipeline' = pipeline,
                           pipeline_method = pipeline_method,
                           'target' = target,
                           connected_filter = connected_filter,
                           'threshold' = x['threshold'],
                           'logFC' = x['logFC'],
                           'Adj.P' = x['Adj.P'],
                           row.names = NULL)

        dplyr::bind_cols(x_df, pipeline_output[['performance']])

    }, mc.cores = n_cores)

    final_output <- dplyr::bind_rows(results_list)

    return(final_output)

}
