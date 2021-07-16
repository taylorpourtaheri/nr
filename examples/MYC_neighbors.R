
library(magrittr)
library(dplyr)
library(igraph)

# read differential expression data (annotated with gene symbols)
# de_string <- readRDS('data/de_string_v11.RDS')
# de_MYC <- de_string$MYC

# generate protein association network
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=950)
ppi <- string_db$get_graph()

# find gene ID
gene_symbol <-'MYC'
xref <- data.frame(symbol = gene_symbol)
xref <- string_db$map(xref, "symbol", removeUnmappedRows=T, quiet=T)

# select MYC neighbors
neighbors <- igraph::neighbors(ppi, xref$STRING_id, 'all')

# double check edge list
ppi_edges <- igraph::as_data_frame(ppi, 'edges')
gene_edges <- dplyr::filter(ppi_edges,
                            (from == xref$STRING_id |to == xref$STRING_id))
nrow(gene_edges) == length(neighbors) # TRUE


# subset network
sub <- induced_subgraph(ppi, neighbors)
df <- igraph::as_data_frame(sub, 'vertices')
rownames(df) <- NULL
colnames(df) <- 'STRING_id'

# join with DE results
df %<>% left_join(de_MYC)

# save
saveRDS(df, 'data/MYC_neighbors_dea_results_v11.RDS')


