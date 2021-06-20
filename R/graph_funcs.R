#' The file contains all the graphing functions

freq.histogram <- function(
  clinical.data, variable,
  colour = "darkblue", label.colour = colour
){
  #' Mutation histograms from clinical frames
  #'
  #' Produce a histogram with the distribution of the mutational frequencies.
  #' @param clinical.data Clinical data from a project
  #' @param variable The variable to plot. Either `"driver.freq"` or `"passenger.freq"`
  #' @param project.name The name of the project
  #'
  #' Made to work with project-like clinical data with a `project.id` attribute.
  library(ggplot2)

  stopifnot(variable %in% c("driver.freq", "passenger.freq"))

  data <- clinical.data
  title <- paste("Number of driver mutations per patient -", get.name(data))

  tabulated_freqs <- table(data[[variable]])

  p <- ggplot(data, aes(x=.data[[variable]])) +
    theme_minimal() +
    geom_bar(stat="count", show.legend=FALSE, fill=colour, width=1) +
    scale_x_continuous(
      breaks=seq(from=0, to=max(data[[variable]]), by=1), minor_breaks = NULL
    ) +
    xlab("Number of Driver Mutations") + ylab("Number of Patients") +
    ggtitle(title)

  # I need to change dynamically the colours of the labels.
  pdata <- ggplot_build(p)
  set_vjust_val <- function(x, t=10){if(x > t){return(2)} else {return(-2)}}
  set_vjust_col <- function(x, t=10){if(x > t){return("white")} else {return(label.colour)}}
  unnested <- pdata[[1]][[1]]
  vjust_values <- unlist(lapply(unnested$count, set_vjust_val, t=max(unnested$count)/10))
  vjust_colours <- unlist(lapply(unnested$count, set_vjust_col, t=max(unnested$count)/10))

  p <- p + geom_text(stat='count', aes(label=..count..), vjust=vjust_values, col=vjust_colours)
  print(p)
}
