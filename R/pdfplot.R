#' @title Save plot as pdf
#'
#' @description A simple internal helperfunction to save the current plot as a
#'              pdf with the option to add a dataframe to a new page of the pdf
#’
#' @param infotab a dataframe
#'                the given datafram will be added to the second page of the pdf
#'                default is set to NULL
#' @param PATH the path where the pdf shall be saved
#' @param name the name that will be given to the pdf - if a pdf with that name
#'             already exists within the dedicated directory - the pdf will be
#'             overwritten
#' @importFrom  grid grid.newpage
#' @importFrom gridExtra tableGrob
#' @importFrom grDevices dev.copy dev.off
#' @return
#'
#' @examples
#'
#' @export
#'


pdfplot <- function(infotab = NULL,
                    # the dataframe that shall be appended
                    PATH,
                    # path of the file
                    name
                    # the name of the file
                    )
  {
  dev.copy(pdf, paste0(PATH, name, ".pdf"))
  grid.newpage()
  grid.table(infotab)
  dev.off()
  }
