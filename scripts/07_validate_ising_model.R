############################################################
# 07_validate_ising_model.R
#
# Evaluate the binary Ising model
#
# Method:
# - Compare observed symptom activation probabilities
#   with probabilities reproduced by binary data simulated
#   from the estimated Ising-model parameters.
# - Agreement is summarized using the mean absolute
#   difference (MAD) across the 22 symptoms.
#
# Inputs:
#   results/05_binary_data/binary_symptoms.rds
#   results/06_Ising_model/ising_model.rds
#
# Outputs:
#   results/07_ising_validation/
############################################################


##############################
# 1. Load package
##############################

library(IsingFit)


##############################
# 2. File paths
##############################

binary_file <- file.path(
  "results",
  "05_binary_data",
  "binary_symptoms.rds"
)

model_file <- file.path(
  "results",
  "06_Ising_model",
  "ising_model.rds"
)

output_dir <- file.path(
  "results",
  "07_ising_validation"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 3. Check input files
##############################

if (!file.exists(binary_file)) {
  stop(
    paste0(
      "Binary symptom data not found: ",
      binary_file,
      "\nPlease run scripts/05_prepare_binary_data.R first."
    )
  )
}

if (!file.exists(model_file)) {
  stop(
    paste0(
      "Ising model not found: ",
      model_file,
      "\nPlease run scripts/06_estimate_ising_model.R first."
    )
  )
}


##############################
# 4. Load data and model
##############################

binary_symptoms <- readRDS(
  binary_file
)

ising_model <- readRDS(
  model_file
)

binary_symptoms <- as.data.frame(
  binary_symptoms
)

symptom_vars <- paste0(
  "Q",
  1:22
)


##############################
# 5. Basic checks
##############################

stopifnot(
  ncol(binary_symptoms) == 22,
  identical(
    names(binary_symptoms),
    symptom_vars
  )
)

stopifnot(
  nrow(ising_model$weiadj) == 22,
  ncol(ising_model$weiadj) == 22,
  length(ising_model$thresholds) == 22
)

if (any(is.na(binary_symptoms))) {
  stop(
    "Missing values were detected in binary symptom data."
  )
}

binary_check <- vapply(
  binary_symptoms,
  function(x) {
    all(x %in% c(0, 1))
  },
  logical(1)
)

if (!all(binary_check)) {
  stop(
    "Input symptom data contain values other than 0 and 1."
  )
}


##############################
# 6. Extract Ising parameters
##############################

edge_weights <- as.matrix(
  ising_model$weiadj
)

thresholds <- as.numeric(
  ising_model$thresholds
)

rownames(edge_weights) <- symptom_vars
colnames(edge_weights) <- symptom_vars
names(thresholds) <- symptom_vars


##############################
# 7. Define simulation function
##############################

simulate_Ising <- function(
    n,
    W,
    thresholds,
    iterations = 200
) {
  
  p <- length(thresholds)
  
  # Random initial binary states
  X <- matrix(
    rbinom(
      n * p,
      size = 1,
      prob = 0.5
    ),
    nrow = n,
    ncol = p
  )
  
  # Gibbs updates
  for (iter in seq_len(iterations)) {
    
    for (i in seq_len(p)) {
      
      eta <-
        thresholds[i] +
        X %*% W[, i]
      
      prob <-
        1 /
        (1 + exp(-eta))
      
      X[, i] <- rbinom(
        n,
        size = 1,
        prob = prob
      )
    }
  }
  
  colnames(X) <- paste0(
    "Q",
    seq_len(p)
  )
  
  X
}


##############################
# 8. Observed activation probabilities
##############################

observed_prevalence <- colMeans(
  binary_symptoms
)


##############################
# 9. Simulate data from Ising model
##############################

N_SIM <- 5000
N_ITER <- 200
SEED_VALIDATION <- 2026

set.seed(
  SEED_VALIDATION
)

simulated_binary <- simulate_Ising(
  n = N_SIM,
  W = edge_weights,
  thresholds = thresholds,
  iterations = N_ITER
)

simulated_prevalence <- colMeans(
  simulated_binary
)


##############################
# 10. Compare observed and simulated
##############################

validation_results <- data.frame(
  Symptom = symptom_vars,
  Observed_Prevalence =
    as.numeric(observed_prevalence),
  Simulated_Prevalence =
    as.numeric(simulated_prevalence)
)

validation_results$Absolute_difference <-
  abs(
    validation_results$Observed_Prevalence -
      validation_results$Simulated_Prevalence
  )


##############################
# 11. Mean absolute difference
##############################

MAD_activation <- mean(
  validation_results$Absolute_difference
)

validation_summary <- data.frame(
  Number_of_symptoms = 22,
  Number_of_simulated_profiles = N_SIM,
  Gibbs_iterations = N_ITER,
  Mean_absolute_difference = MAD_activation
)


##############################
# 12. Save results
##############################

write.csv(
  validation_results,
  file.path(
    output_dir,
    "observed_vs_simulated_prevalence.csv"
  ),
  row.names = FALSE
)

write.csv(
  validation_summary,
  file.path(
    output_dir,
    "ising_validation_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  validation_results,
  file.path(
    output_dir,
    "ising_validation_results.rds"
  )
)


##############################
# 13. Save seed information
##############################

writeLines(
  c(
    paste0(
      "Validation simulation seed: ",
      SEED_VALIDATION
    ),
    paste0(
      "Number of simulated profiles: ",
      N_SIM
    ),
    paste0(
      "Gibbs iterations: ",
      N_ITER
    )
  ),
  file.path(
    output_dir,
    "simulation_settings.txt"
  )
)


##############################
# 14. Display summary
##############################

print(
  validation_results
)

cat(
  "\nMean absolute difference (MAD):",
  MAD_activation,
  "\n"
)

cat(
  "Rounded MAD:",
  round(MAD_activation, 3),
  "\n"
)

cat(
  "Ising model validation completed successfully.\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)