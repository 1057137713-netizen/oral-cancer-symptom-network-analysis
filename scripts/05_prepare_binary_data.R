############################################################
# 05_prepare_binary_data.R
#
# Prepare binary symptom data for Ising model estimation
#
# Coding rule:
#   0   = symptom absent
#   >0  = symptom present
#
# Input:
#   data/raw_data.xlsx
#
# Output:
#   results/05_binary_data/
#     binary_symptoms.csv
#     binary_symptoms.rds
#     binary_symptom_prevalence.csv
############################################################


##############################
# 1. Load package
##############################

library(readxl)


##############################
# 2. File paths
##############################

input_file <- "data/raw_data.xlsx"

output_dir <- "results/05_binary_data"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


##############################
# 3. Import original data
##############################

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Data file not found: ", input_file,
      "\nParticipant-level data are not distributed publicly."
    )
  )
}

data <- read_excel(
  input_file,
  sheet = "整理-ALL"
)


##############################
# 4. Extract symptom variables
##############################

symptom_vars <- paste0("Q", 1:22)

if (!all(symptom_vars %in% names(data))) {
  stop("One or more required symptom variables (Q1-Q22) are missing.")
}

symptoms <- data[, symptom_vars]
symptoms <- as.data.frame(symptoms)

symptoms[] <- lapply(
  symptoms,
  as.numeric
)


##############################
# 5. Basic data checks
##############################

stopifnot(
  ncol(symptoms) == 22,
  identical(
    names(symptoms),
    symptom_vars
  )
)

if (any(is.na(symptoms))) {
  stop(
    "Missing symptom values were detected. The study dataset used for analysis contained no missing symptom data."
  )
}

if (any(
  unlist(symptoms) < 0 |
  unlist(symptoms) > 10
)) {
  stop("Symptom scores outside the expected 0-10 range were detected.")
}

cat("Sample size:", nrow(symptoms), "\n")
cat("Number of symptom variables:", ncol(symptoms), "\n")


##############################
# 6. Dichotomize symptom scores
##############################

# Coding:
# 0 = symptom absent
# >0 = symptom present

binary_symptoms <- as.data.frame(
  lapply(
    symptoms,
    function(x) {
      ifelse(x == 0, 0L, 1L)
    }
  )
)

names(binary_symptoms) <- symptom_vars


##############################
# 7. Verify binary coding
##############################

binary_check <- vapply(
  binary_symptoms,
  function(x) {
    all(x %in% c(0L, 1L))
  },
  logical(1)
)

if (!all(binary_check)) {
  stop("Non-binary values were detected after dichotomization.")
}


##############################
# 8. Calculate symptom prevalence
##############################

symptom_prevalence <- data.frame(
  Symptom = symptom_vars,
  Prevalence = colMeans(binary_symptoms)
)

print(symptom_prevalence)


##############################
# 9. Save outputs
##############################

write.csv(
  binary_symptoms,
  file.path(
    output_dir,
    "binary_symptoms.csv"
  ),
  row.names = FALSE
)

saveRDS(
  binary_symptoms,
  file.path(
    output_dir,
    "binary_symptoms.rds"
  )
)

write.csv(
  symptom_prevalence,
  file.path(
    output_dir,
    "binary_symptom_prevalence.csv"
  ),
  row.names = FALSE
)


##############################
# 10. Completion message
##############################

cat("\nBinary symptom data preparation completed successfully.\n")
cat("Output directory:", output_dir, "\n")