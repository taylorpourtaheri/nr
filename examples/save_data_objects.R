
# save STRING db object
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold = 950)
saveRDS(string_db, 'data/string_db_v11.RDS')

# save PPI object
ppi <- string_db$get_graph()
saveRDS(ppi, 'data/string_ppi_v11.RDS')

# save similarity matrix (jaccard is the default method)
sim <- igraph::similarity(ppi)
saveRDS(sim, 'data/string_ppi_v11_jacc_sim_mat.RDS')

# convert sim matrix to named list and remove all zeros
sim_list <- as.list(as.data.frame(sim))
names(sim_list) <- igraph::V(ppi)$name

condense <- function(sim_vec){
    names(sim_vec) <- igraph::V(ppi)$name
    ind <- which(sim_vec != 0)
    dense_vec <- sim_vec[ind]
}

sim_list_dense <- purrr::map(sim_list, ~condense(.))
saveRDS(sim_list_dense, 'data/string_ppi_v11_jacc_sim_list_dense.RDS')

