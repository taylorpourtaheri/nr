#' @title Subset the loaded igraph based on edge confidence
#' @description Given an iGraph object from StringDb, reduce the
#' full graph to a subset based on edge confidence scores. This
#' function will improve the performance of multiple calls to
#' the iGraph object for the parameter optimization.
#' @param graph Graph of class '\code{igraph}'
#' @param score_threshold Integer value of the edge confidence'

devtools::load_all()


parameter_grid <- function(deg, target, grid){

    # set globals

    thresholds <- unique(grid$threshold)
    threshold_string <- as.character(thresholds)

    string_list <- purrr::map(thresholds, ~ STRINGdb::STRINGdb$new(version = '11',
                                                                   species = 9606,
                                                                   score_threshold = .))
    ppi_list <- purrr::map(string_list, ~ .$get_graph())

    names(string_list) <- threshold_string
    names(ppi_list) <- threshold_string

    # define function to iterate over grid
    grid_workflow <- function(var, string_list, ppi_list){

        print(var)

        subset <- dplyr::filter(grid, threshold == var)
        string_db_temp <- string_list[[var]]
        ppi_temp <- ppi_list[[var]]

        var_output <- c()

        for (i in 1:nrow(subset)){

            print(subset[i,])

            out_temp <- centrality_pipeline(deg = deg,
                                            ppi = ppi_temp,
                                            string_db = string_db_temp,
                                            logFC_min = subset[i,]$logFC,
                                            pvalue_max = subset[i,]$Adj.P,
                                            method = 'betweenness',
                                            causal_gene_symbol =  target,
                                            sim_method = 'jaccard',
                                            n_sim = 9999,
                                            weighted = TRUE)

            results_df <- dplyr::bind_cols(subset[i,], out_temp$performance)

            var_output[[i]] <- results_df
        }

        return(dplyr::bind_rows(var_output))

    }

    # run grid workflow
    output_list <- purrr::map(threshold_string, ~ grid_workflow(., string_list, ppi_list))
    output <- dplyr::bind_rows(output_list)

    return(output)

}


# test case -----------------------------------------------------------------

# Load differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

deg <- myc_de
target <- 'MYC'

threshold <- c(850,875,900,925,950)
logFC <- c(1.50, 1.75)
Adj.P <- c(0.025, 0.05)

grid <- expand.grid(threshold=threshold,logFC=logFC,Adj.P=Adj.P)

# issues loading the full package due to namespace issues? diffuStats
results <- parameter_grid(myc_de,'MYC', grid)
saveRDS(results, 'results/parameter_grid_dev_results2.RDS')
openxlsx::write.xlsx(results, 'results/parameter_grid_dev_results2.xlsx')


# explore results ---------------------------------------------------------

library(ggplot2)

results$logFC <- factor(results$logFC)

ggplot(results, aes(x = threshold, y = mean_score)) +
    geom_point(aes(color = Adj.P, shape = logFC),
               size = 5)
ggsave('results/parameter_grid_dev_results2.png', width = 10, height = 6)







