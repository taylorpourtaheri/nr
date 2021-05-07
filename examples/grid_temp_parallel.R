

devtools::load_all()


# Load differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

deg <- myc_de
target <- 'MYC'

threshold <- c(850,875,900,925,950)
logFC <- c(1.0, 1.25, 1.50, 1.75)
Adj.P <- c(0.01, 0.025, 0.05)

parameter_grid <- expand.grid(threshold=threshold,logFC=logFC,Adj.P=Adj.P)

library(tictoc)

tic()
results <- parameter_grid(myc_de,'MYC', grid, parallel = FALSE)
toc()

# tic()
# results <- parameter_grid(myc_de,'MYC', grid, parallel = TRUE) # 739.779 sec
# toc()

saveRDS(results, 'results/parameter_grid_dev_results4.RDS')
openxlsx::write.xlsx(results, 'results/parameter_grid_dev_results4.xlsx')

# deg <- myc_de
# target <- 'MYC'
# parallel <- FALSE
# grid <- grid[1,]

# explore results ---------------------------------------------------------

library(ggplot2)

results$logFC <- factor(results$logFC)

ggplot(results, aes(x = threshold, y = mean_score)) +
    geom_point(aes(color = Adj.P, shape = logFC),
               size = 5)
ggsave('results/parameter_grid_dev_results4.png', width = 10, height = 6)

ggplot(results, aes(x = threshold, y = score_pval)) +
    geom_point(aes(color = Adj.P, shape = logFC),
               size = 5) +
    scale_y_log10() +
    annotation_logticks(sides = 'l')
ggsave('results/parameter_grid_dev_results4p.png', width = 10, height = 6)






