#' @title
#'
#' @description
#’
#' @param
#'
#' @return
#'
#' @examples
#'
#' @export


shinymouseover <- function(shinyelement,
                           shinytextID,
                           MOtext,
                           defaulttext){
  shinyjs::onevent("mouseenter",
                   shinyelement,
                   shinytextID <- renderText({
                     MOtext
                   }))
  shinyjs::onevent("mouseleave",
                   shinyelement,
                   shinytextID <- renderText({
                     defaulttext
                   }))
}

