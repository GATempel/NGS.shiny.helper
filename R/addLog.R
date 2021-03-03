#' @title Create Log Data
#'
#' @description Adds log data to an existing log dataframe
#’
#' @param statement a character string
#'        The statement that will be written into the file next to the time
#' @param dataframe a dataframe
#'        the dataframe needs at least two columns
#'        the statement will be added to the first column
#'        the Sys.time will be added to the second column
#' @return
#'
#' @examples
#'
#' @export

addLog <- function(statement, dataframe)
{
 if (!is.data.frame(dataframe))
 {
   stop(paste0(deparse(substitute(dataframe)), " is not a Dataframe."))
 } else {
   dataframe[nrow(dataframe) + 1,] <- list(
     statement,
     as.character(Sys.Date()),
     strftime(Sys.time(), format = "%H:%M:%S")
   )
   message(paste0("Finished: ",
                  statement))
   return(dataframe)
 }

}
