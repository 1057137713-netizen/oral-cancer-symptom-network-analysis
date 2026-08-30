############################################################
# 06_estimate_ising_model.R
#
# Estimate binary Ising model
#
# Input:
#   results/05_binary_data/binary_symptoms.rds
#
# Model:
#   IsingFit
#   family = "binomial"
#   gamma = 0.5
#
# Outputs:
#   results/06_Ising_model/
#     ising_model.rds
#     ising_edge_weight_matrix.csv
#     ising_threshold_parameters.csv
############################################################


##############################
# 1. Load package
##############################

library(IsingFit)


##############################
# 2. File paths
##############################

input_file <- file.path(
  "results",
  "05_binary_data",
  "binary_symptoms.rds"
)

output_dir <- file.path(
  "results",
  "06_Ising_model"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 3. Load binary symptom data
##############################

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ",
      input_file,
      "\nPlease run scripts/05_prepare_binary_data.R first."
    )
  )
}

binary_symptoms <- readRDS(
  input_file
)

binary_symptoms <- as.data.frame(
  binary_symptoms
)


##############################
# 4. Basic checks
##############################

symptom_vars <- paste0(
  "Q",
  1:22
)

stopifnot(
  ncol(binary_symptoms) == 22,
  identical(
    names(binary_symptoms),
    symptom_vars
  )
)

binary_check <- vapply(
  binary_symptoms,
  function(x) {
    all(x %in% c(0, 1))
  },
  logical(1)
)

if (!all(binary_check)) {
  stop(
    "The input dataset contains values other than 0 and 1."
  )
}

if (any(is.na(binary_symptoms))) {
  stop(
    "Missing values were detected in the binary symptom data."
  )
}

cat(
  "Sample size:",
  nrow(binary_symptoms),
  "\n"
)

cat(
  "Number of symptoms:",
  ncol(binary_symptoms),
  "\n"
)


##############################
# 5. Estimate Ising model
##############################

# Prespecified seed retained for reproducibility

SEED_ISING <- 2026

set.seed(
  SEED_ISING
)

ising_model <- IsingFit::IsingFit(
  binary_symptoms,
  family = "binomial",
  plot = FALSE,
  gamma = 0.5
)


##############################
# 6. Extract edge-weight matrix
##############################

edge_weight_matrix <- as.matrix(
  ising_model$weiadj
)

rownames(edge_weight_matrix) <- symptom_vars
colnames(edge_weight_matrix) <- symptom_vars


##############################
# 7. Extract threshold parameters
##############################

threshold_vector <- as.numeric(
  ising_model$thresholds
)

names(threshold_vector) <- symptom_vars

threshold_parameters <- data.frame(
  Symptom = symptom_vars,
  Threshold = threshold_vector
)


##############################
# 8. Parameter checks
##############################

stopifnot(
  nrow(edge_weight_matrix) == 22,
  ncol(edge_weight_matrix) == 22,
  length(threshold_vector) == 22,
  all(is.finite(edge_weight_matrix)),
  all(is.finite(threshold_vector))
)

stopifnot(
  isTRUE(
    all.equal(
      edge_weight_matrix,
      t(edge_weight_matrix)
    )
  )
)


##############################
# 9. Save model and parameters
##############################

saveRDS(
  ising_model,
  file.path(
    output_dir,
    "ising_model.rds"
  )
)

write.csv(
  edge_weight_matrix,
  file.path(
    output_dir,
    "ising_edge_weight_matrix.csv"
  ),
  row.names = TRUE
)

write.csv(
  threshold_parameters,
  file.path(
    output_dir,
    "ising_threshold_parameters.csv"
  ),
  row.names = FALSE
)


##############################
# 10. Save seed information
##############################

writeLines(
  paste0(
    "Ising model seed: ",
    SEED_ISING
  ),
  file.path(
    output_dir,
    "random_seed.txt"
  )
)


##############################
# 11. Summary
##############################

cat(
  "\nIsing model estimation completed successfully.\n"
)

cat(
  "Number of non-zero edges:",
  sum(edge_weight_matrix != 0) / 2,
  "\n"
)

cat(
  "Threshold SD:",
  sd(threshold_vector),
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)