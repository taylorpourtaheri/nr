#' db_create
#' 
#' add two numbers
#' 
#' @param vers stingDB version
#' @param spec species code #9606 for human
#' @param thresh score_threshold
#' @export


db_create <- function(vers = "11", spec = 9606, thresh = 950){

  string_db <- STRINGdb::STRINGdb$new(version=vers,
                                    species=spec,
                                    score_threshold=thresh)

  ppi <- string_db$get_graph()
  
  return(ppi)
  }
