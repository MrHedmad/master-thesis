## This file contains the data preprocessing steps.
# The end result is a list containing all the tumor types and their relative data
# Saved it for reusage as a R object

library(plyr)
library(tidyverse)
library(reshape2)

projects.base.path <- "F:/Data/University/Thesis/Data/CHASMplus/Projects"

projects.names <- list.files(projects.base.path)

projects.paths <- as.list(rep(NA, length(projects.names)))
names(projects.paths) <- projects.names

for (name in projects.names) {
  projects.paths[[name]] <- list(
    clinical = file.path(projects.base.path, name, paste(name, "_clinical.csv", sep="")),
    crystal = file.path(projects.base.path, name, paste(".proj_", name, "_temp", sep=""), paste(name, "_crystal.csv", sep=""))
  )
}

### Preprocessing functions

preprocess_clinical <- function(data){
  # I need to clean some of the columns as they hold no meaning.
  data$X <- NULL
  data$X0 <- NULL
  # TCGA IDs that are useless
  data$demographic_id <- NULL
  data$diagnosis_id <- NULL
  data$exposure_id <- NULL
  # The following columns hold metadata relative to the TCGA project, and are therefore not relevant
  # Moreover, pandas duplicates them when fusing the clinical data, making them even more meaningless
  # I could compress this but I'm lazy
  data$updated_datetime <- NULL
  data$updated_datetime_x <- NULL # The x and y at the end here are the pandas aberrations
  data$updated_datetime_y <- NULL
  data$created_datetime <- NULL
  data$created_datetime_x <- NULL # The x and y at the end here are the pandas aberrations
  data$created_datetime_y <- NULL
  data$state <- NULL
  data$state_x <- NULL
  data$state_y <- NULL
  # There are many factors, I'll make them into characters so we can filter out the not reported:
  library(dplyr)
  data %>% mutate_if(is.factor, as.character) -> data
  # Change the "not reported" with R-interpretable NAs
  data[data == "not reported"] <- NA
  data[data == "Not Reported"] <- NA
  data[data == ""] <- NA
  data <- subset(data, !is.na(data$passenger_freq))
  # I'm guessing these fields where empty in the pandas dataframe and thus result empty here, so I'm setting them to NA.
  # Back to factors
  data %>% mutate_if(is.character, as.factor) -> data
  return(data)
}

clean_crystal_id <- function(data){
  return(unlist(lapply(data$Identifier, substr, start=1, stop=12)))
}

clean_id_labels <- function(data){
  labs <- names(data)
  labs <- replace(labs, labs %in% c("Identifier", "ParticipantBarcode", "submitter_id"), "ID")
  names(data) <- labs
  return(data)
}

combined_preprocess <- function(dls){
  dls$clinical <- preprocess_clinical(dls$clinical)
  # Crop IDs
  dls$crystal$Identifier <- clean_crystal_id(dls$crystal)
  # Standardize IDs
  dls <- mapply(clean_id_labels, dls)
  return(dls)
}

run_preprocess <- function(projects.paths, name){
  project <- list(clinical = read.csv(projects.paths[[name]]$clinical), crystal = read.csv(projects.paths[[name]]$crystal))
  project <- combined_preprocess(project)
  return(project)
}

projects <- as.list(rep(NA, length(projects.names)))
names(projects) <- projects.names

for (name in projects.names) {
  projects[[name]] <- as.list(run_preprocess(projects.paths, name))
}

# Dirname goes up one folder
save(projects, file = file.path(dirname(projects.base.path), "preprocessed_CHASM_data.RData"))
