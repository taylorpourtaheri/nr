#' db_create
#' 
#' Download protein-protein interaction (ppi) network and return in the form of
#' an iGraph object for processing.
#' 
#' @param vers A string representing the stringDB version.
#' @param spec A number specifying the species code (default is 9606 for humans.)
#' @param thresh A number specifying score_threshold (a measure of interaction strength.)
#' @return An iGraph object representing the ppi network.
#' @examples
#' \dontrun{
#' ppi <- db_create()
#' ppi_mouse <- db_create(spec = 10090, thresh = 975)
#' }
#' @export


db_create <- function(vers = "11", spec = 9606, thresh = 950){

  string_db <- STRINGdb::STRINGdb$new(version=vers,
                                    species=spec,
                                    score_threshold=thresh)

  ppi <- string_db$get_graph()
  
  return(ppi)
  }