#' @title Get Number of Unique sequences
#'
#' @description Determines the amount of unique sequences in a dada object
#'              (generated with dada2) or a list of dada objects
#’
#' @param x a dada object or a list of dada objects
#' @importFrom dada2 getUniques
#' @return
#'
#' @examples
#'
#' @export


getN <- function(x){
  N <- sum(dada2::getUniques(x))
  return(N)
}
