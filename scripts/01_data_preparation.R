############################################################
# 01_data_preparation.R
#
# Prepare the original continuous symptom data used for
# exploratory factor analysis and Gaussian graphical
# network analysis.
#
# Input:
#   data/raw_data.xlsx
#
# Output:
#   results/01_data_preparation/
#     continuous_symptoms.rds
#     continuous_symptoms.csv
#     data_preparation_summary.csv
#
# Note:
# Participant-level clinical data are not included in the
# public repository because of ethical and privacy restrictions.
############################################################


##############################
# 1. Load package
##############################

library(readxl)


##############################
# 2. File paths
##############################

input_file <- file.path(
  "data",
  "raw_data.xlsx"
)

output_dir <- file.path(
  "results",
  "01_data_preparation"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 3. Check input file
##############################

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Data file not found: ",
      input_file,
      "\nParticipant-level clinical data are not distributed publicly."
    )
  )
}


##############################
# 4. Import original dataset
##############################

data <- read_excel(
  input_file
)


##############################
# 5. Define symptom variables
##############################

symptom_vars <- paste0(
  "Q",
  1:22
)

if (!all(
  symptom_vars %in% names(data)
)) {
  stop(
    "One or more required symptom variables (Q1-Q22) are missing."
  )
}


##############################
# 6. Extract Q1-Q22
##############################

continuous_symptoms <-
  data[
    ,
    symptom_vars
  ]

continuous_symptoms <-
  as.data.frame(
    continuous_symptoms
  )


##############################
# 7. Convert to numeric
##############################

continuous_symptoms[] <-
  lapply(
    continuous_symptoms,
    as.numeric
  )


##############################
# 8. Check dimensions
##############################

stopifnot(
  ncol(continuous_symptoms) == 22,
  identical(
    names(continuous_symptoms),
    symptom_vars
  )
)

cat(
  "Sample size:",
  nrow(continuous_symptoms),
  "\n"
)

cat(
  "Number of symptom variables:",
  ncol(continuous_symptoms),
  "\n"
)


##############################
# 9. Check missing values
##############################

missing_values <-
  colSums(
    is.na(
      continuous_symptoms
    )
  )

if (any(
  missing_values > 0
)) {
  
  print(
    missing_values
  )
  
  stop(
    "Missing symptom values were detected."
  )
}


##############################
# 10. Check score range
##############################

all_values <-
  unlist(
    continuous_symptoms,
    use.names = FALSE
  )

if (any(
  all_values < 0 |
  all_values > 10
)) {
  
  stop(
    "Symptom scores outside the expected 0-10 range were detected."
  )
}


##############################
# 11. Create preparation summary
##############################

data_summary <- data.frame(
  Sample_size =
    nrow(continuous_symptoms),
  
  Number_of_symptoms =
    ncol(continuous_symptoms),
  
  Missing_values =
    sum(
      is.na(
        continuous_symptoms
      )
    ),
  
  Minimum_score =
    min(
      all_values
    ),
  
  Maximum_score =
    max(
      all_values
    )
)


##############################
# 12. Save continuous data
##############################

saveRDS(
  continuous_symptoms,
  file.path(
    output_dir,
    "continuous_symptoms.rds"
  )
)

write.csv(
  continuous_symptoms,
  file.path(
    output_dir,
    "continuous_symptoms.csv"
  ),
  row.names = FALSE
)

write.csv(
  data_summary,
  file.path(
    output_dir,
    "data_preparation_summary.csv"
  ),
  row.names = FALSE
)


##############################
# 13. Completion message
##############################

cat(
  "\nContinuous symptom data preparation completed successfully.\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)