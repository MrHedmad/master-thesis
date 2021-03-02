load(file = "F:/Data/University/Thesis/Data/CHASMplus/preprocessed_CHASM_data.RData")
projects.names <- names(projects)

generate_graphs <- function(project, save_path, project.id){
  pdf(file.path(save_path, paste("Frequency Histogram - ", project.id, ".pdf", sep="")))
  MutationHist(project$clinical, title = paste("Number of mutations per patient -", project.id))
  dev.off()

  pdf(file.path(save_path, paste("Mutation Heatmap - ", project.id, ".pdf", sep="")))
  HeatMutPlot2(project, project.id, use_renamed = FALSE)
  dev.off()

  pdf(file.path(save_path, paste("Jitterplot - ", project.id, ".pdf", sep="")))
  plotjitter(project$clinical, "passenger_freq", "Number of Passenger mutations", title = paste("Driver vs Passenger -", project.id))
  dev.off()
}

############################################### Runs ###########################################

setwd("F:/Data/University/Thesis/Data/CHASMplusImages/")

for (name in projects.names) {
  generate_graphs(projects[[name]], "F:/Data/University/Thesis/Data/CHASMplusImages/", name)
}
