#' @title Create Folder
#'
#' @description Checks if the given folder name exists in the provided directory
#'              and creates it if not
#’
#' @param input_path a character string
#'        takes a directory path to check if folder_name exists within
#' @param folder_name a character string
#'        the supplied string will be the name of the created folder
#'        in default with recursive = TRUE it will accept a PATH as well as
#'        long as the PATH does NOT start with a slash
#' @param variable_name a character string
#'        the supplied string will serve as the name of the variable that will
#'        be assigned to the combined string of input_path and folder_name
#' @param variable_name logical
#'        will be parsed to the recursive argument of the dir.create function
#'        if TRUE elements other than just the last in the combined string will
#'        be created
#'
#' @return No direct return - but the given variable_name will be assigned to
#'         the combined input_path/folder_name string
#'
#' @examples
#' create_folder("~", "testFolder", "tsPath")
#'
#' @export

create_folder <- function(input_path,
                          folder_name,
                          variable_name,
                          recursive = TRUE)
  {
  if (is.character(input_path) == F |
      is.character(folder_name) == F |
      is.character(variable_name) == F)
  {
    stop("At least one of the supplied inputs it not a character string.")
  }
  assign(variable_name,
         paste0(input_path,
                "/",
                folder_name),
         envir = .GlobalEnv)
  if (file.exists(get(variable_name)) == F)
    {
    dir.create(get(variable_name),
               recursive = recursive)
    }
  }
