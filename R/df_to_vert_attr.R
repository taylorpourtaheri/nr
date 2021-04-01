#' @title Set the vertex attributes of a graph
#' @description Given a graph and a data frame where vertex
#' names correspond to a column in the data frame, map data
#' in the data frame to attributes of the vertices in the graph.
#' @param graph Graph of class '\code{igraph}'
#' @param df A data frame that has a column that contains vertex IDs
#' and at least one column that has data that will become a vertex attribute of \code{graph}.
#' @param common A character string providing the column name that contains
#' the vertex IDs (i.e., the 'name' vertex attribute of the graph).
#' @param attr_name A character vector providing the name(s) of the
#' column of \code{df} that will be added as vertex attributes to \code{graph}
#' @param new_attr_name A character vector providing the name(s) of the
#' new vertex attribute(s). If \code{NULL} then the name of the
#' new attribute(s) will match \code{attr_name}.
#' @export
df_to_vert_attr <- function(graph, df, common, attr_name, new_attr_name=NULL) {

    # some checks
    if (!is.null(new_attr_name) & (length(attr_name) != length(new_attr_name))) {
        stop("Length of attr_name must be the same as new_attr_name.")
    }
    if (!is.atomic(attr_name)) {
        stop("attr_name is not atomic.")
    }
    if (!is.atomic(new_attr_name)) {
        stop("new_attr_name is not atomic.")
    }
    if (!is.character(attr_name)) {
        warning("attr_name is not a character vector; attempting coersion.")
        attr_name <- as.character(attr_name)
    }
    if (!is.character(new_attr_name) & !is.null(new_attr_name)) {
        warning("new_attr_name is not a character vector; attempting coersion.")
        new_attr_name <- as.character(new_attr_name)
    }
    ind <- match(igraph::V(graph)$name, df[[common]])
    df <- df[ind, ]
    if (nrow(df) == 0 | all(is.na(ind))) {
        stop(glue::glue("None of values in the column '{common}' match vertex names."))
    }

    # set new names if new_attr_name is NULL
    if (is.null(new_attr_name)) {
        new_attr_name <- attr_name
    }

    # set vertex attributes
    for (i in 1:length(attr_name)) {
        graph <- igraph::set_vertex_attr(graph, new_attr_name[i], value = df[[attr_name[i]]])
    }

    return(graph)

}
