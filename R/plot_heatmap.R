#' @title plot heatmap
#'
#' @description An internal helper function to create a ggplot heatmap from the
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
#' @return plot - a ggplot heatmap
#' @import ggplot2
#' @import dplyr
#' @import plotly
#'
#' @examples
#'
#' @export

plot_heatmap <- function(
  genera,
  samples,
  lowercutoff,
  uppercutoff,
  modus = "default"
)
{
  if (length(samples) == 0){
    messagetext <- list(
      x = 1,
      y = 1,
      text = "No sample selected \nPlease select a sample",
      showarrow = FALSE,
      font = list(size = 28)
    )
    plot <- plotly::plot_ly(as.data.frame(NULL))
    plot <- plot %>% layout(annotations = messagetext,
                            yaxis = list(visible = FALSE),
                            xaxis = list(visible = FALSE))
  } else {
    plotdata <- dplyr::filter(genera$Tidy, Sample %in% samples)
    if (modus == "total")
    {
      cutNames <- genera$Sorted_Names[lowercutoff:uppercutoff]
    } else
    {
      cutNames <- dplyr::select(genera$Abundance, all_of(samples))
      cutNames$Total <- rowSums(cutNames)
      cutNames <- cutNames[order(cutNames$Total,
                                 decreasing = T),]
      cutNames$Names <- row.names(cutNames)
      cutNames <- cutNames$Names
      cutNames <- cutNames[lowercutoff:uppercutoff]
    }
    plotdata <- dplyr::filter(plotdata, Name %in% cutNames)
    plot <- ggplot2::ggplot(plotdata, ggplot2::aes(x = Sample,
                                                   y = Name,
                                                   fill = Abundance))
    plot <- plot + ggplot2::geom_raster()
    plot <- plot + ggplot2::theme_classic()
    plot <- plot + ggplot2::ylab(as.character(substitute(genera)))
    plot <- plotly::ggplotly(plot)
  }

  return(plot)
}
