evaluate_performance <- function(network, network_df, causal_sim,
                                 method, n_sim, weighted = FALSE){

    if (weighted == FALSE){

        mean_score <- mean(network_df$causal_similarity)

        # estimate uncertainty with a random draw of the full ppi graph
        n_draws <- length(igraph::V(network))
        samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))

        sample_means <- sapply(samples, mean)

        # calculate p
        score_pval <- sum(sample_means > mean_score) / n_sim

    }

    if (weighted == TRUE){

        mean_score <- mean(network_df$weighted_score)
        weights <- network_df[,method]

        # estimate uncertainty with a random draw of the full ppi graph
        n_draws <- length(igraph::V(network))
        samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))

        weighted_samples <- purrr::map(samples, ~ .*weights)
        sample_means <- sapply(weighted_samples, mean)

        # find z-score
        simulation_mean <- mean(sample_means)
        simulation_sd <- sd(sample_means)
        network_z_score <- (mean_score - simulation_mean) / simulation_sd

        # calculate p
        score_pval <- sum(sample_means > mean_score) / n_sim

    }

    output_df <- data.frame('mean_score' = mean_score,
                            'z_score' = network_z_score,
                            'score_pval' = score_pval)

    return(output_df)

}
