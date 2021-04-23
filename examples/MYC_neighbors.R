
library(magrittr)
library(dplyr)
library(igraph)

# read differential expression data (annotated with gene symbols)
de_string <- readRDS('data/de_string_v11.RDS')
de_MYC <- de_string$MYC

# in the DEG results:
MYC_stringID <-'9606.ENSP00000479618'

# generate protein association network
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=edge_conf_score_min)
ppi <- string_db$get_graph()

# select MYC neighbors
neighbors <- igraph::neighbors(ppi, MYC_stringID, 'all')

# subset network
sub <- induced_subgraph(ppi, neighbors)
df <- igraph::as_data_frame(sub, 'vertices')
rownames(df) <- NULL
colnames(df) <- 'STRING_id'

# join with DE results
df %<>% left_join(de_MYC)

# save
saveRDS(df, 'data/MYC_neighbors_dea_results_v11.RDS')


