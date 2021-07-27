evaluate_performance <- function(target, network_df, causal_sim,
                                 method, n_sim, weighted = FALSE){

    if (weighted == FALSE){

        target_in_network <- any(network_df$Symbol==target)

        if(target_in_network){
            rank <- as.numeric(rownames(network_df[match(target,network_df$Symbol),]))
        }
        else{
            rank <- NULL
        }

        mean_score <- mean(network_df$causal_similarity)

        # estimate uncertainty with a random draw of the full ppi graph
        n_draws <- nrow(network_df)
        samples <- lapply(1:n_sim, function(x) sample(causal_sim, n_draws))

        sample_means <- sapply(samples, mean)

        # find z-score
        simulation_mean <- mean(sample_means)
        simulation_sd <- sd(sample_means)
        network_z_score <- (mean_score - simulation_mean) / simulation_sd

        # calculate p
        score_pval <- sum(sample_means > mean_score) / n_sim

    }

    if (weighted == TRUE){

        target_in_network <- any(network_df$Symbol==target)

        if(target_in_network){
            rank <- as.numeric(rownames(network_df[match(target,network_df$Symbol),]))
        }
        else{
            rank <- NULL
        }

        mean_score <- mean(network_df$weighted_score)
        weights <- network_df[,method]

        # estimate uncertainty with a random draw of the full ppi graph
        n_draws <- nrow(network_df)
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

    target_in_network <- any(network_df$Symbol==target)

    output_df <- data.frame('mean_score' = mean_score,
                            'sample_mean' = simulation_mean,
                            'sample_sd' = simulation_sd,
                            'z_score' = network_z_score,
                            'score_pval' = score_pval,
                            'target_in_network' = target_in_network,
                            'rank' = rank,
                            'trimmed_network' = nrow(network_df),
                            'full_network' = length(causal_sim))

    output_list <- list()
    output_list[['performance_results']] <- output_df
    output_list[['simulation_scores']] <- data.frame('simulation_mean_score' = sample_means)

    return(output_list)

}
