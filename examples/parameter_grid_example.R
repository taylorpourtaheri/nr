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
#deg <- RAS_de
#target <- 'HRAS'

# select MYC condition as an example
myc_de <- de_string$MYC
deg <- myc_de
target <- 'MYC'


# generate propagation parameter grid
pipeline_vec <- c('raw','ml','gm','ber_s',
                  'mc','ber_p', 'z', 'random_walk')
# pipeline_vec <- c('gm')

threshold_vec <- c(950)
Adj.P_vec <- seq(0.05, 0.25, 0.05)
logFC_vec <- seq(0.0, 2.0, 0.5)

# single example
# pipeline_vec <- c('centrality')
# connected_filter <- c('TRUE', 'FALSE')
# threshold_vec <- c(950)
# logFC_vec <- 2.0
# Adj.P_vec <- 0.05

parameter_grid <- expand.grid(pipeline_method = pipeline_vec,
                              threshold=threshold_vec,
                              logFC=logFC_vec,
                              Adj.P=Adj.P_vec)

# run function for centrality

tic()
results <- parameter_grid(deg = deg,
                          target = target,
                          grid = parameter_grid,
                          n_cores = 4,
                          weighted = TRUE)
toc()

# saveRDS(results, 'results/MYC/all_methods/parameter_grid_propagation_full_MYC.RDS')
# openxlsx::write.xlsx(results, 'results/MYC/all_methods/parameter_grid_propagation_full_MYC.xlsx')

results <- readRDS('results/MYC/all_methods/parameter_grid_propagation_full_MYC.RDS')

# explore results ---------------------------------------------------------

# scale and center z-score
results$scaled_z_score <- scale(results$z_score)

best <- group_by(results, pipeline_method) %>% top_n(1, z_score)


ggplot(results, aes(x = logFC, y = scaled_z_score, group = pipeline_method)) +
    geom_line(aes(color = pipeline_method)) +
    geom_point(size = 1.75) +
    geom_point(data = best, color = 'red', size = 1.75) +
    facet_grid(pipeline_method ~Adj.P, scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
          ,legend.position = 'none'
          ) +
    # scale_color_manual(values = c('#F8766D', '#00BFC4'),
    #                    labels = c("Centrality", "Propagation")) +
    labs(title = 'Scaled Z-score vs log fold-change, by method pipeline',
         x = 'Log fold-change threshold',
         y = 'Scaled Z-score',
         color = 'Propagation method')
# ggsave('results/MYC/all_methods/parameter_grid_propagation_full_MYC_logFC.png',
#        width = 10, height = 6)
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

ggplot(results %>% filter(pipeline_method == 'raw'), aes(x = logFC, y = Adj.P)) +
    geom_tile(aes(fill = z_score)) +
    geom_text(aes(label = round(score_pval, 3))) +
    scale_fill_continuous(low = 'white', high = '#0072B2') +
    labs(title = 'Z-score vs. log fold-change and adjusted p-value, propagation pipeline',
         subtitle = 'Metric: raw, connected_filter = FALSE',
         x = 'Log Fold-Change',
         y = 'Adjusted p-value',
         fill = 'Z-score')

ggsave('results/MYC/parameter_grid_pvalue_logFC_heatmap_FALSE.png', width = 15, height = 8)

