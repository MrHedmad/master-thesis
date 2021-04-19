## This file contains the data preprocessing steps.
# The end result is a list containing all the tumor types and their relative data
# Saved it for usage as a R object

library(plyr)
library(tidyverse)
library(reshape2)

preprocess_CHASM_tumors <- function(projects.base.path, refactor.vars) {
  #' Preprocess CHASMplus projects coming from the Python scripts.
  #'
  #' The data is divided into a crystal and a clinical .csv files. This
  #' function reads them all and combines them into a single list, from
  #' now on called a "project", for easy access.
  #'
  #' It also restructures the data, parsing missing data as NA, reshaping the
  #' levels of some variables to be used for modelling purposes.
  #'
  #' @param projects.path The path to the folder containing each project folder.
  #' Each folder is expected to be named with the name of the TCGA project of
  #' reference, like `TCGA-BRCA`.
  #' @param refactor.vars Sent to `refactor_clinical_variables`.

  # Get all of the projects in the folder and read them into project lists
  # First, we generate all the paths to the files that we need for each project
  projects.names <- list.files(projects.base.path)
  projects.paths <- as.list(rep(NA, length(projects.names)))
  names(projects.paths) <- projects.names

  for (name in projects.names) {
    projects.paths[[name]] <- list(
      clinical = file.path(
        projects.base.path, name, paste(name, "_clinical.csv", sep="")
      ),
      crystal = file.path(
        projects.base.path, name, paste(".proj_", name, "_temp", sep=""),
        paste(name, "_crystal.csv", sep="")
      )
    )
  }

  projects <- as.list(rep(NA, length(projects.names)))
  names(projects) <- projects.names

  for (name in projects.names) {
    # Read the data
    project <- list(
      clinical = read.csv(projects.paths[[name]]$clinical),
      crystal = read.csv(projects.paths[[name]]$crystal)
      )
    # Run the preprocessing functions
    project$clinical %>%
      purge_useless_columns() %>%
      reshape_to_na() %>%
      plyr::rename(
        replace = c(frequency = "driver.freq", passenger_freq = "passenger.freq"),
        warn_missing = FALSE
        ) %>%
      refactor_clinical_variables(refactor.vars, name) ->
      project$clinical

    project$crystal %>% clean_crystal_id() -> project$crystal

    project <- mapply(rename_id_labels, project)

    project %>% reshape_crystal() -> project$gene_labels

    projects[[name]] <- project
  }

  return(projects)
}

### Preprocessing functions

purge_useless_columns <- function(clinical.data) {
  #' Remove useless columns from a clinical data frame
  #'
  #' These columns come from python merging, and they hold no meaning for the
  #' analysis. They mostly hold metadata for the entries, or are just merging
  #' artefacts.

  clinical.data %>% select(
    !c(X, X0, demographic_id, diagnosis_id, exposure_id, updated_datetime,
    created_datetime, state, created_datetime_x, created_datetime_y,
    state_x, state_y, updated_datetime_x, updated_datetime_y),
  ) -> clinical.data

  return(clinical.data)
}

reshape_to_na <- function(clinical.data) {
  #' Replace NA-meaning values to actual NAs.
  #'
  #' The TCGA consortium uses `"not reported"` as its NA, as well as
  #' just empty strings as NAs. This replaces them to actual NAs.

  # As I do a global rename, we need character vectors, and not factors
  clinical.data %>% mutate_if(is.factor, as.character) -> clinical.data

  clinical.data[clinical.data == "not reported"] <- NA
  clinical.data[clinical.data == "Not Reported"] <- NA
  clinical.data[clinical.data == ""] <- NA

  clinical.data %>% mutate_if(is.character, as.factor) -> clinical.data

  return(clinical.data)
}

clean_crystal_id <- function(crystal.data){
  #' Clean the IDs of a crystal, which might not be standardized.
  #'
  #' We expect that the ID in the crystal begins with the TCGA patient ID.
  #' This function crops the names to just the ID (the first 12 characters).
  crystal.data$Identifier <-
    unlist(lapply(crystal.data$Identifier, substr, start=1, stop=12))
  return(crystal.data)
}

rename_id_labels <- function(dframe){
  #' Rename columns that contain patient IDs to just `"ID"`
  #'
  #' The TCGA data uses different names to mean the same thing, so I rename
  #' all of them the just the one, uniform name `"ID"`

  # Cannot use rename() here as it fails if it doesn't find the column to rename.
  names(dframe) %>%
    replace(
      names(dframe) %in% c("Identifier", "ParticipantBarcode", "submitter_id"),
      "ID"
      ) -> names(dframe)
  return(dframe)
}

reshape_crystal <- function(project){
  #' Reshape a project to a gene label data frame
  #'
  #' Moves from a project to a large data frame containing the presence or absence
  #' status of mutations encoded into 0 for not present and 1 for present.
  #' Each row is a different
  relevant <- project$crystal[,c("ID", "Renamed.HUGO", "Hugo")]
  # I stole this from stack overflow
  gene_labels <- relevant %>%
    pivot_longer(ends_with("hugo"), names_to = "genes") %>%
    complete(ID, value) %>%
    dplyr::mutate(genes = as.numeric(!is.na(genes))) %>%
    pivot_wider(names_from = value, values_from = genes, values_fn = max)
  # Putting `all.x=TRUE, all.y=TRUE` makes some NAs appear, which is weird....
  # TODO: Check if the length of the gene labels matrix is the same as the crystal,
  # but nothing should have been lost by the trasposition
  temp <- data.frame(ID = project$clinical[,c("ID")])
  gene_labels <- merge(gene_labels, temp, by = "ID", all = TRUE)
  # See disclaimer above. I think that the NAs come from 0-driver patients, that have no ID in the
  # pivoted gene_labels table, so the real labels should be zero
  gene_labels[is.na(gene_labels)] <- 0
  gene_labels <- merge(gene_labels, project$clinical[, c("ID", "passenger.freq", "driver.freq")], by = "ID")
  # These are very long vectors but with just two levels. I could cast them
  # to factors but some functions later on fail.
  return(gene_labels)
}

refactor_clinical_variables <- function(clinical.data, considered.vars, tumor.type) {
  #' Simplify the clinical variable factors
  #'
  #' There are a lot - a lot - of levels in the clinical variables.
  #' This is bad for modelling, so I simplify them here to have less levels.
  #' The comments in the function itself  explain my thought process.
  #' @param clinical.data The clinical data to refactor
  #' @param considered.vars The variables to consider for refactoring. They will
  #' all be cast to factors, then evaluated for refactoring.
  #' @param tumor.type The complete project name (ex. `TCGA-BRCA`). Some refactors
  #' only apply to certain tumor types, this specifies it.
  #'
  #' This function is not particularly efficient as it has to loop through
  #' all columns and test if it needs to reshape them. There is probably
  #' a more efficient method (I'm looking at you, `apply`), but I trade
  #' execution speed for development speed. Sorry, future me.
  #' I also change "unknown" and "not reported" to NAs, even though this
  #' should have already been done before.

  for (colname in colnames(clinical.data)) {
    if (colname %in% considered.vars) {
      # Convert the factorial or factorial-like data to factors for reshaping
      clinical.data[[colname]] <- as.factor(clinical.data[[colname]])
    }
    # Clinical/pathological ajcc scores considerations:
    # The scores are too precise to work well as model factors, so I simplify
    # them here into generic "low" and generic "high"
    if (colname == "ajcc_clinical_m") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "M0" = "M0",
          "M1" = "M1",
          "M1a" = "M1",
          "M1b" = "M1",
          "M1c" = "M1",
          "MX" = NA, # TODO: Is this an NA?
          "cM0 (i+)" = "M0",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_clinical_n") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "N0" = "N0",
          "N0 (i+)" = "N0",
          "N0 (i-)" = "N0",
          "N0 (mol+)" = "N0",
          "N0 (mol-)" = "N0",
          "N1" = "N1",
          "N1a" = "N1",
          "N1b" = "N1",
          "N1bI" = "N1",
          "N1bII" = "N1",
          "N1bIII" = "N1",
          "N1bIV" = "N1",
          "N1c" = "N1",
          "N1mi" = "N1",
          "N2" = "N2+",
          "N2a" = "N2+",
          "N2b" = "N2+",
          "N2c" = "N2+",
          "N3" = "N2+",
          "N3a" = "N2+",
          "N3b" = "N2+",
          "N3c" = "N2+",
          "N4" = "N2+",
          "NX" = NA, # TODO: Is this an NA?
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_clinical_stage") {
      clinical.data[[colname]] <- revalue(
        clinical.data[[colname]], c(
          "Stage 0" = "Stage0",
          "Stage 0a" = "Stage0",
          "Stage 0is" = "Stage0",
          "Stage I" = "Stage1",
          "Stage IA" = "Stage1",
          "Stage IA1" = "Stage1",
          "Stage IA2" = "Stage1",
          "Stage IA3" = "Stage1",
          "Stage IB" = "Stage1",
          "Stage IB1" = "Stage1",
          "Stage IB2" = "Stage1",
          "Stage IC" = "Stage1",
          "Stage II" = "Stage2",
          "Stage IIA" = "Stage2",
          "Stage IIA1" = "Stage2",
          "Stage IIA2" = "Stage2",
          "Stage IIB" = "Stage2",
          "Stage IIC" = "Stage2",
          "Stage IIC1" = "Stage2",
          "Stage III" = "Stage3+",
          "Stage IIIA" = "Stage3+",
          "Stage IIIB" = "Stage3+",
          "Stage IIIC" = "Stage3+",
          "Stage IIIC1" = "Stage3+",
          "Stage IIIC2" = "Stage3+",
          "Stage IS" = "Stage1",
          "Stage IV" = "Stage1",
          "Stage IVA" = "Stage1",
          "Stage IVB" = "Stage1",
          "Stage IVC" = "Stage1",
          # 'Tis' are small mutated cells - some do not consider this as cancer
          "Stage Tis" = "Stage0",
          "Stage X" = NA, # Is this an NA?
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_clinical_t") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "T0" = "T0", # TODO: is this an NA?
          "T1" = "T1",
          "T1a" = "T1",
          "T1a1" = "T1",
          "T1a2" = "T1",
          "T1b" = "T1",
          "T1b1" = "T1",
          "T1b2" = "T1",
          "T1c" = "T1",
          "T1mi" = "T1",
          "T2" = "T2",
          "T2a" = "T2",
          "T2a1" = "T2",
          "T2a2" = "T2",
          "T2b" = "T2",
          "T2c" = "T2",
          "T2d" = "T2",
          "T3" = "T3+",
          "T3a" = "T3+",
          "T3b" = "T3+",
          "T3c" = "T3+",
          "T3d" = "T3+",
          "T4" = "T3+",
          "T4a" = "T3+",
          "T4b" = "T3+",
          "T4c" = "T3+",
          "T4d" = "T3+",
          "T4e" = "T3+",
          "TX" = NA,
          "Ta" = "T1",
          "Tis" = "T0",
          "Tis (DCIS)" = "T0",
          "Tis (LCIS)" = "T0",
          "Tis (Paget's)" = "T0",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_pathologic_m") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "M0" = "M0",
          "M1" = "M1",
          "M1a" = "M1",
          "M1b" = "M1",
          "M1c" = "M1",
          "M1d" = "M1",
          "M2" = "M2",
          "MX" = NA, # TODO: Is this an NA?
          "cM0 (i+)" = "M0",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_pathologic_n") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "N0" = "N0",
          "N0 (i+)" = "N0",
          "N0 (i-)" = "N0",
          "N0 (mol+)" = "N0",
          "N0 (mol-)" = "N0",
          "N1" = "N1",
          "N1a" = "N1",
          "N1b" = "N1",
          "N1bI" = "N1",
          "N1bII" = "N1",
          "N1bIII" = "N1",
          "N1bIV" = "N1",
          "N1c" = "N1",
          "N1mi" = "N1",
          "N2" = "N2+",
          "N2a" = "N2+",
          "N2b" = "N2+",
          "N2c" = "N2+",
          "N3" = "N2+",
          "N3a" = "N2+",
          "N3b" = "N2+",
          "N3c" = "N2+",
          "N4" = "N2+",
          "NX" = NA, # TODO: Is this an NA?
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    # Note for stages. This is a bit arbitrary, but I simplified them to
    # stage 0, 1, 2 and 3 or more (3+)
    if (colname == "ajcc_pathologic_stage") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "Stage 0" = "Stage0",
          "Stage 0a" = "Stage0",
          "Stage 0is" = "Stage0",
          "Stage I" = "Stage1",
          "Stage IA" = "Stage1",
          "Stage IA1" = "Stage1",
          "Stage IA2" = "Stage1",
          "Stage IA3" = "Stage1",
          "Stage IB" = "Stage1",
          "Stage IB1" = "Stage1",
          "Stage IB2" = "Stage1",
          "Stage IC" = "Stage1",
          "Stage II" = "Stage2",
          "Stage IIA" = "Stage2",
          "Stage IIA1" = "Stage2",
          "Stage IIA2" = "Stage2",
          "Stage IIB" = "Stage2",
          "Stage IIC" = "Stage2",
          "Stage IIC1" = "Stage2",
          "Stage III" = "Stage3+",
          "Stage IIIA" = "Stage3+",
          "Stage IIIB" = "Stage3+",
          "Stage IIIC" = "Stage3+",
          "Stage IIIC1" = "Stage3+",
          "Stage IIIC2" = "Stage3+",
          "Stage IS" = "Stage1",
          "Stage IV" = "Stage1",
          "Stage IVA" = "Stage1",
          "Stage IVB" = "Stage1",
          "Stage IVC" = "Stage1",
          # 'Tis' are small mutated cells - some do not consider this as cancer
          "Stage Tis" = "Stage0",
          "Stage X" = NA,
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ajcc_pathologic_t") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "T0" = "T0",
          "T1" = "T1",
          "T1a" = "T1",
          "T1a1" = "T1",
          "T1a2" = "T1",
          "T1b" = "T1",
          "T1b1" = "T1",
          "T1b2" = "T1",
          "T1c" = "T1",
          "T1mi" = "T1",
          "T2" = "T2",
          "T2a" = "T2",
          "T2a1" = "T2",
          "T2a2" = "T2",
          "T2b" = "T2",
          "T2c" = "T2",
          "T2d" = "T2",
          "T3" = "T3+",
          "T3a" = "T3+",
          "T3b" = "T3+",
          "T3c" = "T3+",
          "T3d" = "T3+",
          "T4" = "T3+",
          "T4a" = "T3+",
          "T4b" = "T3+",
          "T4c" = "T3+",
          "T4d" = "T3+",
          "T4e" = "T3+",
          "TX" = NA,  # TODO: Is this an NA?
          "Ta" = "T1",
          "Tis" = "T0",
          "Tis (DCIS)" = "T0",
          "Tis (LCIS)" = "T0",
          "Tis (Paget's)" = "T0",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "gender") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "female" = "female",
          "male" = "male",
          "unknown" = NA,
          "unspecified" = NA,
          "not reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "vital_status") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "Alive" = "Alive",
          "Dead" = "Dead",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "prior_malignancy") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "yes" = "yes",
          "no" = "no",
          "unknown" = NA,
          "not reported" = NA,
          "Not Allowed To Collect" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "alcohol_history") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "Yes" = "Yes",
          "No" = "No",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "ann_arbor_b_symptoms") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "yes" = "yes",
          "no" = "no",
          "unknown" = NA,
          "not reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "figo_stage") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "Stage 0" = "Stage0",
          "Stage I" = "Stage1",
          "Stage IA" = "Stage1",
          "Stage IA1" = "Stage1",
          "Stage IA2" = "Stage1",
          "Stage IB" = "Stage1",
          "Stage IB1" = "Stage1",
          "Stage IB2" = "Stage1",
          "Stage IC" = "Stage1",
          "Stage IC1" = "Stage1",
          "Stage IC2" = "Stage1",
          "Stage IC3" = "Stage1",
          "Stage II" = "Stage2",
          "Stage IIA" = "Stage2",
          "Stage IIA1" = "Stage2",
          "Stage IIA2" = "Stage2",
          "Stage IIB" = "Stage2",
          "Stage IIC" = "Stage2",
          "Stage III" = "Stage3+",
          "Stage IIIA" = "Stage3+",
          "Stage IIIA1" = "Stage3+",
          "Stage IIIA2" = "Stage3+",
          "Stage IIIAi" = "Stage3+",
          "Stage IIIAii" = "Stage3+",
          "Stage IIIB" = "Stage3+",
          "Stage IIIC" = "Stage3+",
          "Stage IIIC1" = "Stage3+",
          "Stage IIIC2" = "Stage3+",
          "Stage IV" = "Stage3+",
          "Stage IVA" = "Stage3+",
          "Stage IVB" = "Stage3+",
          "Unknown" = NA,
          "Not Reported" = NA),
        warn_missing = FALSE
      )
    }
    if (colname == "masaoka_stage") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          "Stage I" = "Stage1",
          "Stage IIa" = "Stage2",
          "Stage IIb" = "Stage2",
          "Stage III" = "Stage3+",
          "Stage IVa" = "Stage3+",
          "Stage IVb" = "Stage3+"),
        warn_missing = FALSE
      )
    }
    # The ICD10 codes are simplified to their more generic counterpats.
    # The source for the meaning of the codes is https://icd.codes/
    # The renaming is cancer-specific as some codes can reappear in different
    # types of cancer and they mean slightly different things.
    if (colname == "icd_10_code" & tumor.type == "TCGA-COAD") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
          # Colon
          "C18.0" = "initialcolon", # Cecum
          "C18.2" = "initialcolon", # Ascending colon
          "C18.3" = "initialcolon", # Hepatic Flexure
          "C18.4" = "medialcolon", # transverse colon
          "C18.5" = "medialcolon", # splenic flexure
          "C18.6" = "terminalcolon", # descending colon
          "C18.7" = "terminalcolon", # sigmoid colon
          "C18.9" = "unspecifiedcolon",
          "C19" = "terminalcolon" # Rectosigmoid Junction
        ),
        warn_missing = FALSE
      )
    }
    if (colname == "icd_10_code" & tumor.type == "TCGA-SKCM") {
      clinical.data[[colname]] <- revalue(
        clinical.data[[colname]], c(
          # Melanoma
          "C07" = "gland", # Parotid Gland
          "C17.9" = "mucosae", # Small intestine
          "C18.9" = "mucosae", # Unspecified colon
          "C22.0" = "gland", # Liver cell carcinoma
          "C34.1" = "mucosae", # Upper lobe of bronchi or lung
          "C34.3" = "mucosae", # Lower lobe of bronchi or lung
          "C34.9" = "mucosae", # Bronchi or lung, unspecified
          "C41.0" = "connective", # Bone of skull and face
          "C43.51" = "skin", # Malignant melanoma (M.m.) of anal skin
          # I wanted to keep the head and neck skin patients separated, but
          # there are not many of them, so I keep them all together
          "C44.2" = "skin", # M.m. of ear
          "C44.3" = "skin", # M.m. of face
          "C44.31" = "skin", # M.m. of eyelid
          "C44.4" = "skin", # M.m. of neck
          "C44.5" = "skin", # M.m. of trunk
          "C44.50" = "skin", # M.m. of trunk
          "C44.6" = "connective", # Malignant neoplasm of penis
          "C44.601" = "connective", # Malignant neoplasm of penis glans
          "C44.7" = "skin", # Lower limb skin
          "C44.701" = "skin", # Lower limb skin
          "C44.9" = "skin", # Generic code for skin malignant neoplasm
          "C48.2" = "mucosae", # Peritoneum, unspecified
          "C49.0" = "connective", # Other connective or soft tissue of head, face or neck
          "C49.1" = "connective", # Other connective or soft tissue of upper limb
          "C49.2" = "connective", # Other connective or soft tissue of lower limb
          "C49.20" = "connective", # See above
          "C49.3" = "connective", # Other connective or soft tissue of thorax
          "C49.4" = "connective", # Other connective or soft tissue of abdomen
          "C49.5" = "connective", # Other connective or soft tissue of pelvis
          "C49.6" = "connective", # Other connective or soft tissue of trunk
          "C49.9" = "connective", # Other connective or soft tissue, unspecified
          "C50.9" = "connective", # Unspecified breast
          "C51.9" = "mucosae", # Malignant neoplasm of Vulva
          "C52" = "mucosae", # Vagina
          "C54.1" = "mucosae", # Corpus uteri, edometrium
          "C71.1" = "nervous", # Brain, except lobes and ventricles
          "C71.3" = "nervous", # Parietal lobe of brain
          "C71.9" = "nervous", # Brain, unspecified
          "C72.0" = "nervous", # CNS, spinal cord
          "C74.9" = "gland", # Adrenal gland, unspecified
          "C76.1" = "generic", # Other site, thorax
          "C76.2" = "generic", # Other site, abdomen
          "C76.3" = "generic", # Other site, pelvis
          "C77.0" = "metastasis", # Secondary neoplasm in lymph nodes - head, face or neck
          "C77.2" = "metastasis", # Secondary neoplasm in lymph nodes - abdominal
          "C77.3" = "metastasis", # Secondary neoplasm in lymph nodes - axilla or upper limb
          "C77.4" = "metastasis", # Secondary neoplasm in lymph nodes - inguinal or lower limb
          "C77.5" = "metastasis", # Secondary neoplasm in lymph nodes - pelvis
          "C77.9" = "metastasis" # Secondary neoplasm in lymph nodes - unspecified
        ),
        warn_missing = FALSE
      )
    }
    if (colname == "icd_10_code" & tumor.type == "TCGA-BRCA") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
        # Breast
        "C50.2" = "breast", # Upper inner quadrant
        "C50.3" = "breast", # Lower inner quadrant
        "C50.4" = "breast", # Upper outer quadrant
        "C50.5" = "breast", # Lower outer quadrant
        "C50.8" = "breast", # Overlapping sites
        "C50.9" = "breast", # Unspecified
        "C50.919" = "breast" # Unspecified, Female
      ),
      warn_missing = FALSE
      )
    }
    if (colname == "icd_10_code" & tumor.type == "TCGA-LUAD") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
        # Lung
        "C34.0" = "mainbronchus", # Main bronchus
        "C34.1" = "lobe", # Upper lobe, bronchus or lung
        "C34.2" = "lobe", # Middle lobe, bronchus or lung
        "C34.3" = "lobe", # Lower lobe, bronchus or lung
        "C34.30" = "lobe", # Lower lobe, unspecified
        "C34.8" = "diffuse", # Overlapping sites
        "C34.9"  = "generic" # Unspecified
      ),
      warn_missing = FALSE
      )
    }
    if (colname == "icd_10_code" & tumor.type == "TCGA-GBM") {
      clinical.data[[colname]] <- revalue(clinical.data[[colname]], c(
        # Colon
        "C71.1" = "cerebrum", # Cerebrum, except lobes and ventricles
        "C71.2" = "temporal", # Temporal Lobe
        "C71.4" = "occipital", # Occipital Lobe
        "C71.8" = "overlapping", # Overlapping sites
        "C71.9" = "generic" # Unspecified
      ),
      warn_missing = FALSE
      )
    }
    # primary and secondary gleason grade are weird...
    # I am not sure how to rename their values, but since it only applies to a single tumor,
    # I just disregard it for now.
  }
  return(clinical.data)
}

# This is the folder containing all the projects' data, each in a separated folder.
projects.base.path <- "F:/Data/University/Thesis/Data/CHASMplus/Projects"

refactor.vars <- c(
  "gender", "vital_status", "icd_10_code",
  "prior_malignancy", "alcohol_history", "ajcc_clinical_m",
  "ajcc_clinical_n", "ajcc_clinical_t", "ajcc_clinical_stage",
  "ajcc_pathologic_m", "ajcc_pathologic_n", "ajcc_pathologic_t",
  "ajcc_pathologic_stage", "ann_arbor_b_symptoms", "figo_stage",
  "masaoka_stage", "primary_gleason_grade", "secondary_gleason_grade"
)

projects <- preprocess_CHASM_tumors(projects.base.path, refactor.vars)
# Dirname goes up one folder
save(projects, file = file.path(dirname(projects.base.path), "preprocessed_CHASM_data.RData"))
# Subset the data to just the tumors of interest and save it again
reduced.projects <- projects[c("TCGA-COAD", "TCGA-BRCA", "TCGA-LUAD", "TCGA-SKCM", "TCGA-GBM")]
save(
  reduced.projects,
  file = file.path(dirname(projects.base.path), "reduced_preprocessed_CHASM_data.RData")
  )
