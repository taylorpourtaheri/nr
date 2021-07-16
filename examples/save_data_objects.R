
# save STRING db object
string_db <- STRINGdb::STRINGdb$new(version="11",
                                    species=9606,
                                    score_threshold=edge_conf_score_min)
saveRDS(string_db, 'data/string_db_v11.RDS')

# save PPI object
ppi <- string_db$get_graph()
saveRDS(ppi, 'data/string_ppi_v11.RDS')

# save similarity matrix (jaccard is the default method)
sim <- igraph::similarity(ppi)
saveRDS(sim, 'data/string_ppi_v11_jacc_sim_mat.RDS')
