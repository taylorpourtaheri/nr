
library(tictoc) # to assess performance of parallel processing

devtools::load_all()

# Load differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC
deg <- myc_de
target <- 'MYC'

# generate parameter grid
pipeline_vec <- c('centrality', 'propagation')
threshold_vec <- c(950)
logFC_vec <- c(1.75)
Adj.P_vec <- c(0.025)

parameter_grid <- expand.grid(pipeline = pipeline_vec,
                              threshold=threshold_vec,
                              logFC=logFC_vec,
                              Adj.P=Adj.P_vec)

# run function

# sequential - 186.026s / 360.622s
# tic()
# results <- parameter_grid(deg = myc_de,
#                           target = 'MYC',
#                           grid = parameter_grid,
#                           parallel = FALSE,
#                           weighted = TRUE)
# toc()

# parallel - 125.848s / 306.146s
tic()
results <- parameter_grid(deg = myc_de,
                          target = 'MYC',
                          grid = parameter_grid,
                          n_cores = 2,
                          weighted = TRUE)
toc()

saveRDS(results, 'results/parameter_grid_dev_results6.RDS')
openxlsx::write.xlsx(results, 'results/parameter_grid_dev_results6.xlsx')

# explore results ---------------------------------------------------------

library(ggplot2)
library(dplyr)

comparison <- group_by(results, pipeline) %>% summarize(score_min = min(mean_score),
                                                        score_max = max(mean_score))



results$logFC <- factor(results$logFC)

# ggplot(results, aes(x = threshold, y = mean_score)) +
#     geom_point(aes(color = Adj.P, shape = logFC),
#                size = 5)
# ggsave('results/parameter_grid_dev_results4.png', width = 10, height = 6)

ggplot(results, aes(x = threshold, y = mean_score)) +
    geom_point(aes(color = Adj.P, shape = logFC),
               size = 5) +
    facet_grid(.~pipeline + logFC, scales = 'free') +
    scale_y_log10() +
    annotation_logticks(sides = 'l')
ggsave('results/parameter_grid_dev_results6.png', width = 10, height = 6)






