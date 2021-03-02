library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Barplots
# Histograms showing the number of passenger and driver mutations
MutationHist <- function(data, title="Number of mutations per patient"){
  library(ggplot2)
  tabulated_freqs <- table(data$frequency)
  set_vjust_col <- function(count){
    if(count > 15){
      return("white")
    }else{
      return("blue")
    }
  }
  plot <- ggplot(data, aes(x=frequency)) +
    theme_minimal() +
    geom_bar(stat="count", show.legend=FALSE, fill="darkblue", width=1) +
    scale_x_continuous(breaks=seq(from=0, to=max(data$frequency), by=1), minor_breaks = NULL) +
    xlab("Number of Driver Mutations") + ylab("Number of Patients") +
    ggtitle(title)

  pdata <- ggplot_build(plot)
  set_vjust_val <- function(x, t=10){if(x > t){return(2)} else {return(-2)}}
  set_vjust_col <- function(x, t=10){if(x > t){return("white")} else {return("darkblue")}}
  unnested <- pdata[[1]][[1]]
  vjust_values <- unlist(lapply(unnested$count, set_vjust_val, t=max(unnested$count)/10))
  vjust_colours <- unlist(lapply(unnested$count, set_vjust_col, t=max(unnested$count)/10))

  plot <- plot + geom_text(stat='count', aes(label=..count..), vjust=vjust_values, col=vjust_colours)
  print(plot)
}

# BRCA preview
#MutationHist(projects$`TCGA-BRCA`$clinical)


## Heat maps
give.colors <- function(number, spacing="equal"){
  palette <- c("#e3e3e3", "#a6d7ff", "#fa48c5", "#e6193b", "#f28d6b", "#ffd940", "#bbe329", "#6af23d")
  if(number > length(palette)){
    stop("Number of asked colors too large for palette")
  }
  if(spacing=="equal"){
    sample.diff <- floor(length(palette)/number)
    return(palette[seq(from=0, to=number, by=sample.diff)])
  }
  if(spacing=="order"){
    return(palette[seq(from=0, to=number)])
  }
}

# TODO: This needs some tweaking... The colours are all wrong, and the number of mutations are static in the final analysis.

HeatMutPlot = function(data, nr = 10, title="", color.palette="RdBu", max_col_breaks=4){
  mut <- sort(table(data$Hugo), decreasing = TRUE)
  mut <- mut[1:nr]
  mut.names <- names(mut)
  # Filter data to only retain top nr mut lines
  filtered <- subset(data, as.character(data$Hugo) %in% mut.names)
  # Remove unused factor levels
  filtered[] <- lapply(filtered, function(x) if(is.factor(x)) factor(x) else x)
  # Make table
  filtered <- group_by(filtered, ID)
  frequency <- table(filtered$ID, filtered$Hugo)
  # Warn if not enough colours are selected
  if(max_col_breaks < max(frequency)){
    warning("The number of colours selected is smaller than the maximum number of frequencies. This may cause several frequencies to have the same colour.")
  }
  color.breaks <- seq(-0.5, max_col_breaks, by=1)
  color.palette <- suppressWarnings(brewer.pal(n=length(color.breaks)-1, name=color.palette))
  if(length(color.breaks) == 3){
    color.palette <- c(color.palette[1], color.palette[2])
  }
  heatmap.2(frequency,
            trace = "none",
            dendrogram = "row",
            margins = c(7, 2.2),
            density.info = "none",
            main = title,
            xlab = "Mutated genes",
            ylab = "Patient samples",
            col = color.palette,
            srtCol = 45,  # Rotate the column labels
            labRow = NA,
            key.xlab = "Mutation Nr.",
            srtRow = 90,  # Rotate the row labels
            hclustfun = function(x)hclust(x,method = 'ward.D2'),
            breaks = color.breaks
  )

}

HeatMutPlot2 <- function(data, project.id, use_renamed = FALSE, gene_nr = 10, max_freq = NA){
  data <- data$crystal
  if (use_renamed){
    top_mut_genes <- sort(table(data$Renamed.HUGO), decreasing = TRUE)
  }else{
    top_mut_genes <- sort(table(data$Hugo), decreasing = TRUE)
  }
  top_n_mut_genes <- top_mut_genes[1:gene_nr]
  genes_to_include <- names(top_n_mut_genes)
  # Filter data to only retain top nr mut lines
  if (use_renamed){
    filtered <- subset(data, as.character(data$Renamed.HUGO) %in% genes_to_include)
  }else{
    filtered <- subset(data, as.character(data$Hugo) %in% genes_to_include)
  }
  # The previous step leaves a lot of unused factor levels. We remove them here
  filtered[] <- lapply(filtered, function(x) if(is.factor(x)) factor(x) else x)
  filtered <- group_by(filtered, ID)
  if (use_renamed){
    gene_frequencies <- table(filtered$ID, filtered$Renamed.HUGO)
  }else{
    gene_frequencies <- table(filtered$ID, filtered$Hugo)
  }
  # No idea but as.matrix doesn't work, so I need to do this badness
  table_rownames <- rownames(gene_frequencies)
  table_colnames <- colnames(gene_frequencies)
  gene_frequencies <- matrix(gene_frequencies, ncol = length(table_colnames), nrow = length(table_rownames))
  #rownames(gene_frequencies) <- table_rownames
  colnames(gene_frequencies) <- table_colnames

  if (is.na(max_freq)) {
    max_freq <- max(gene_frequencies)
  }

  colours <- colorRamp2(breaks = c(0, 0.0001, max_freq), colors = c("gray", "white", "red"))

  print(Heatmap(
    gene_frequencies,
    col = colours,
    row_title = "Patients",
    column_title = paste("Mutational frequency of top", gene_nr, "most-mutated genes per patient -", project.id),
    name = "Freq",
    show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    column_names_rot = 45
  ))
}

# Preview BRCA
#HeatMutPlot2(projects$`TCGA-BRCA`, "TCGA-BRCA", use_renamed = FALSE)

## Driver vs Passenger scatterplots
plotjitter <- function(
    data,
    response,
    response_label = response,
    na.rm=TRUE,
    title=paste("Driver vs", response_label),
    distribution.breaks = c("red" = 10, "purple" = 25, "blue" = 50, "purple" = 75, "red" = 90)
  ){
  name <- deparse(substitute(data))
  if(na.rm){
    data <- subset(data, !is.na(data[[response]]))
  }
  p <- ggplot(data, aes(y=frequency, x=.data[[response]])) +
    geom_jitter(height=0.1, alpha = 0.6, na.rm = TRUE) +
    xlab(response_label) +
    ylab("Number of Driver Mutations") +
    theme_minimal() +
    scale_x_log10() +
    ggtitle(title)

  if (!is.null(dist)) {
    for (i in seq(from = 1, to = length(distribution.breaks))){
      p <- p + geom_vline(xintercept = quantile(data[[response]], distribution.breaks[i]/100), color=names(distribution.breaks[i]), alpha=0.5)
    }
  }
  print(p)
}

# BRCA preview
#plotjitter(projects$`TCGA-BRCA`$clinical, "passenger_freq", "Number of Passenger mutations", title = "Driver vs Passenger - TCGA-BRCA")
