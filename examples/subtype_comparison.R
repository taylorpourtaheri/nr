library(magrittr)
library(textshape)
library(ggplot2)
library(dplyr)
library(tidyr)

myc_results <- readRDS('results/MYC/parameter_grid_pvalue_logFC.RDS')
ras_results <- readRDS('results/HRAS/parameter_grid_pvalue_logFC.RDS')

ggplot(myc_results, aes(x = logFC, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))

ggplot(ras_results, aes(x = logFC, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))




ggplot(myc_results, aes(x = Adj.P, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))

ggplot(ras_results, aes(x = Adj.P, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))




ggplot(myc_results, aes(x = target_neighbors_in_final, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))

ggplot(ras_results, aes(x = target_neighbors_in_final, y = mean_score, group=connected_filter)) +
    facet_wrap(~connected_filter, scales = "free") +
    geom_line(aes(color=connected_filter)) +
    geom_line(aes(y = sample_mean, color='sample_mean')) +
    geom_line(aes(y = sample_sd, color='sample_sd'))


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
