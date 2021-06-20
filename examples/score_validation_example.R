library(magrittr)
library(igraph)
library(STRINGdb)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tictoc) # to assess performance of parallel processing

devtools::load_all()

causal_gene_symbol <- 'MYC'
sim_method <- 'jaccard'
pvalue_max <- 0.05
logFC_min <- 2


# load differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')

# select MYC condition as an example
myc_de <- de_string$MYC

# select significant genes
sig_genes <- dplyr::filter(myc_de, ((abs(logFC) > log2(logFC_min)) & (adj.P.Val < pvalue_max)))

# generate protein association network
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=950)
ppi <- string_db$get_graph()

# find the STRING ID for the causal gene
xref <- data.frame(symbol = causal_gene_symbol)
xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

# calculate similarity of each node and slice out the causal gene row
sim <- igraph::similarity(ppi, method = sim_method)
index <- which(igraph::V(ppi)$name == xref$STRING_id)
causal_sim <- sim[index,]

# make scores a named vector
names(causal_sim) <- igraph::V(ppi)$name

# get the scores associated with the significant genes
sig_genes$myc_sim_score <- causal_sim[sig_genes$STRING_id]

# arrange by myc_sim_score
sig_genes %<>% dplyr::arrange(desc(myc_sim_score))

# to do:
# write a function that collects causal gene neighbors (1 degree and 2 degree)
# from ppi graph, and checks the final graph for whether they're present

# read in Pat's sig genes df
sig_genes_pat <- readRDS('data/sig_genes_test.rds')

View(sig_genes == sig_genes_pat)
