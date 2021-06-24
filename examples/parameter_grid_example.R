library(magrittr)
library(textshape)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tictoc) # to assess performance of parallel processing

devtools::load_all()

# Load differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select HRAS condition as an example
#RAS_de <- de_string$RAS
myc_de <- de_string$MYC
#deg <- RAS_de
deg <- myc_de
target <- 'MYC'
#target <- 'HRAS'

# generate parameter grid
pipeline_vec <- c('strength','degree','avg_strength','evcent_w',
                  'evcent_uw','betweenness')
connected_filter <- c('TRUE')
threshold_vec <- c(950)
Adj.P_vec <- seq(0.05, 0.25, 0.05)
logFC_vec <- seq(0.0, 2.0, 0.5)
#Adj.P_vec <- seq(0.005, 0.05, 0.005)

# single example
# pipeline_vec <- c('centrality')
# connected_filter <- c('TRUE', 'FALSE')
# threshold_vec <- c(950)
# logFC_vec <- 2.0
# Adj.P_vec <- 0.05

parameter_grid <- expand.grid(pipeline_method = pipeline_vec,
                              connected_filter = connected_filter,
                              threshold=threshold_vec,
                              logFC=logFC_vec,
                              Adj.P=Adj.P_vec)

# run function for centrality

tic()
results <- parameter_grid(deg = deg,
                          target = target,
                          grid = parameter_grid,
                          #n_cores = 2,
                          #n_cores = 4,
                          weighted = TRUE)
toc()

saveRDS(results, 'results/parameter_grid_centrality_full_MYC.RDS')
openxlsx::write.xlsx(results, 'results/HRAS/parameter_grid_pvalue_logFC.xlsx')

# explore results ---------------------------------------------------------

# best <- group_by(results, pipeline) %>% top_n(1, z_score)
#
#
# ggplot(results, aes(x = Adj.P, y = z_score, group = pipeline)) +
#     geom_line(aes(color = pipeline)) +
#     geom_point(size = 1.75) +
#     geom_point(data = best, color = 'red', size = 1.75) +
#     facet_grid(. ~logFC, scales = 'free') +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
#     scale_color_manual(values = c('#F8766D', '#00BFC4'),
#                        labels = c("Centrality", "Propagation")) +
#     labs(title = 'Z-score vs p-value, by method pipeline',
#          x = 'Adjusted p-value threshold',
#          y = 'Z-score',
#          color = 'Pipeline')
# ggsave('results/parameter_grid_pvalue.png', width = 14, height = 5)
#
# ggplot(results %>% filter(pipeline == 'centrality'),
#        aes(x = Adj.P, y = z_score, group = pipeline)) +
#     geom_line(color = '#F8766D') +
#     geom_point(size = 1.75) +
#     geom_point(data = best %>% filter(pipeline == 'centrality'),
#                color = 'red', size = 1.75) +
#     facet_grid(. ~logFC, scales = 'free') +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
#     labs(title = 'Z-score vs p-value, network centrality pipeline',
#          subtitle = 'Metric: betweenness',
#          x = 'Adjusted p-value threshold',
#          y = 'Z-score')
# ggsave('results/parameter_grid_pvalue_centrality.png', width = 15, height = 4)
#
# ggplot(results %>% filter(pipeline == 'propagation'),
#        aes(x = Adj.P, y = z_score, group = pipeline)) +
#     geom_line(color = '#00BFC4') +
#     geom_point(size = 1.75) +
#     geom_point(data = best %>% filter(pipeline == 'propagation'),
#                color = 'red', size = 1.75) +
#     facet_grid(. ~logFC, scales = 'free') +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
#     labs(title = 'Z-score vs p-value, network propagation pipeline',
#          subtitle = 'Metric: diffusion (raw)',
#          x = 'Adjusted p-value threshold',
#          y = 'Z-score')
# ggsave('results/parameter_grid_pvalue_propagation.png', width = 15, height = 4)


# -------------------------------------------------------------------------

# with heatmap()

results_clip <- results %>% select(logFC, Adj.P, z_score)

results_mat <- results_clip %>%
    # Convert long-form to wide-form
    spread(key = logFC, value = z_score) %>%
    column_to_rownames('Adj.P') %>%
    as.matrix


heatmap(results_mat, Colv = NA, Rowv = NA, scale = "column")

# with ggplot

ggplot(results %>% filter(connected_filter == FALSE), aes(x = logFC, y = Adj.P)) +
    geom_tile(aes(fill = z_score)) +
    geom_text(aes(label = round(score_pval, 3))) +
    scale_fill_continuous(low = 'white', high = '#0072B2') +
    labs(title = 'Z-score vs. log fold-change and adjusted p-value, centrality pipeline',
         subtitle = 'Metric: betweenness, connected_filter = FALSE',
         x = 'Log Fold-Change',
         y = 'Adjusted p-value',
         fill = 'Z-score')

ggsave('results/HRAS/parameter_grid_pvalue_logFC_heatmap_FALSE.png', width = 15, height = 8)

