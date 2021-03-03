#' @title plot barplot
#'
#' @description An internal helper function to create a ggplot barplot from the
#'              data provided to the function. Takes data that has been created
#'              with the IlluminaAnalysis function.
#'
#' @param genera a list data object containing a dataframe called "Abundance"
#'        genera lists are created with the IlluminaAnalysis function.
#'        Defines the dataframe from which to take the data to plot.
#' @param samples a character string or a vector/list containing multiple
#'        strings. Each string must represent a column name within in the
#'        genera$Abundance dataframe.
#'        Defines which data to plot.
#' @param cutoff a integer value to determine how many of the most abundant
#'        values (rownames of the genera$Abundance) of the selected samples are
#'        to be displayed.
#'        Defines how much of the data is plotted
#' @param modus a character string. Either "total", "hide_others" or "absolute"
#'        everything else will result in the default "relative" modus.
#'        Defines how the data is plotted.
#' @return plot - a ggplot barplot
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import plotly
#'
#' @examples
#'
#' @export

plot_barplot <- function(
  genera,
  samples,
  cutoff,
  modus
){
  if (length(samples) == 0){
    messagetext <- list(
      x = 1,
      y = 1,
      text = "No sample selected \nPlease select a sample",
      showarrow = FALSE,
      font = list(size = 28)
    )
    plot <- plot_ly(as.data.frame(NULL))
    plot <- plot %>% layout(annotations = messagetext,
                            yaxis = list(visible = FALSE),
                            xaxis = list(visible = FALSE))
  } else {
    plotdata <- dplyr::select(genera$Abundance, all_of(samples))
    if (ncol(plotdata) > 1) {
      plotdata$Total <- rowSums(plotdata)
      circ = F
    } else if (ncol(plotdata) == 1) {
      plotdata$copy <- plotdata[1]
      plotdata$Total <- rowSums(plotdata) / 2
      circ = T
    }
    plotdata <- plotdata[order(plotdata$Total,
                               decreasing = T),]
    plotdata <- rbind(plotdata[1:cutoff,],
                      colSums(plotdata[(cutoff + 1):nrow(plotdata),]))
    row.names(plotdata)[nrow(plotdata)] <- "#others"
    if (circ == T | modus == "total"){
      plotdata <- dplyr::select(plotdata, Total)
      plotdata$Name <- row.names(plotdata)
      plot <- plotly::plot_ly(plotdata,
                              labels = ~Name,
                              values = ~Total,
                              type = "pie")
    } else {
      plotdata <- dplyr::select(plotdata, -Total)
      if (modus == "hide_others"){
        plotdata <- dplyr::slice_head(plotdata, n = (nrow(plotdata) - 1))
      }
      datapoints <- colnames(plotdata)
      plotdata$Name <- row.names(plotdata)
      plotdata <- tidyr::pivot_longer(plotdata,
                                      datapoints,
                                      names_to = "Sample",
                                      values_to = "Abundance")
      if (modus %in% c("absolute", "hide_others")){
        plot <- ggplot2::ggplot(plotdata,
                                ggplot2::aes(fill = Name,
                                             y = Abundance,
                                             x = Sample))
        plot <- plot + ggplot2::geom_bar(position = "stack",
                                         stat = "identity")
        plot <- plot + ggplot2::theme_classic()
        plot <- plotly::ggplotly(plot)
      } else {
        plot <- ggplot2::ggplot(plotdata,
                                ggplot2::aes(fill = Name,
                                             y = Abundance,
                                             x = Sample))
        plot <- plot + ggplot2::geom_bar(position = "fill",
                                         stat = "identity")
        plot <- plot + ggplot2::theme_classic()
        plot <- plotly::ggplotly(plot)
      }
    }
  }


  return(plot)
}
