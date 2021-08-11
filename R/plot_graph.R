#' @title Plot final network.
#' @description Plot network with standard visualization specifications.
#' @import ggplot2
#' @import ggnetwork
#' @param graph Graph of class '\code{igraph}'. Must include node
#' attributes '\code{name}' and '\code{seed}'.
#' @param method String. Name of node attribute that should correspond to node size.
#' @param gene_list List of strings. Symbols of primary genes of interest.
#' @export
plot_graph <- function(graph, method, gene_list = NULL) {

    ggn <- ggnetwork::ggnetwork(graph)

    as.name(method)

    if(!is.null(gene_list)){
        plot <- ggplot2::ggplot(ggn, aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_edges() +
            geom_nodes(aes_string(color = 'logFC', size = method), alpha = 0.65) +
            geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
            geom_nodelabel_repel(data=subset(ggn, Symbol %in% gene_list), aes(label=Symbol)) +
            scale_color_gradient(low = 'blue', high = 'red') +
            scale_size_continuous(range = c(5, 25)) +
            theme_blank()
    }else{

    plot <- ggplot2::ggplot(ggn, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges() +
        geom_nodes(aes_string(color = 'logFC', size = method), alpha = 0.65) +
        geom_nodetext_repel(aes(label = Symbol), size = 2.5) +
        # geom_nodelabel_repel(data=subset(ggn, Symbol %in% gene_list), aes(label=Symbol)) +
        scale_color_gradient(low = 'blue', high = 'red') +
        scale_size_continuous(range = c(5, 25)) +
        theme_blank()
    }

    return(plot)
}
