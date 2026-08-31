############################################################
# 08_NIRA_primary_2SD.R
#
# Primary NIRA analysis using ±2 SD threshold perturbations
#
# Key principles:
# 1. The original binary Ising model is estimated only once.
# 2. The edge-weight matrix is held fixed.
# 3. One symptom threshold is perturbed at a time.
# 4. All other thresholds remain unchanged.
# 5. No post-intervention networks are re-estimated.
# 6. A common simulated baseline is used for all interventions.
#
#  Input:
#   Preferred local input:
#     results/06_Ising_model/ising_model.rds
#   Public fallback inputs:
#     intermediate_parameters/ising_edge_weight_matrix.csv
#     intermediate_parameters/ising_threshold_parameters.csv
#
# Output:
#   results/08_NIRA_primary_2SD/
############################################################


##############################
# 1. Load package
##############################

library(IsingSampler)


##############################
# 2. File paths
##############################

model_file <- file.path(
  "results",
  "06_Ising_model",
  "ising_model.rds"
)

edge_file <- file.path(
  "intermediate_parameters",
  "ising_edge_weight_matrix.csv"
)

threshold_file <- file.path(
  "intermediate_parameters",
  "ising_threshold_parameters.csv"
)

output_dir <- file.path(
  "results",
  "08_NIRA_primary_2SD"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 3. Load fixed Ising parameters
##############################

# If the locally fitted Ising model is available, use it.
# Otherwise, use the public intermediate parameter files
# distributed with the repository.

if (file.exists(model_file)) {
  
  ising_model <- readRDS(
    model_file
  )
  
  edge_weights <- as.matrix(
    ising_model$weiadj
  )
  
  thresholds <- as.numeric(
    ising_model$thresholds
  )
  
  cat(
    "Using locally fitted Ising model.\n"
  )
  
} else {
  
  if (!file.exists(edge_file) ||
      !file.exists(threshold_file)) {
    
    stop(
      paste0(
        "Neither the local fitted Ising model nor the ",
        "public intermediate parameter files were found."
      )
    )
  }
  
  edge_weights <- as.matrix(
    read.csv(
      edge_file,
      row.names = 1,
      check.names = FALSE
    )
  )
  
  threshold_data <- read.csv(
    threshold_file,
    check.names = FALSE
  )
  
  thresholds <- as.numeric(
    threshold_data[[ncol(threshold_data)]]
  )
  
  cat(
    "Using public intermediate Ising parameter files.\n"
  )
}


##############################
# 4. Define parameter names
##############################

symptom_codes <- paste0(
  "Q",
  seq_along(thresholds)
)

rownames(edge_weights) <- symptom_codes
colnames(edge_weights) <- symptom_codes
names(thresholds) <- symptom_codes



##############################
# 5. Basic checks
##############################

stopifnot(
  nrow(edge_weights) == 22,
  ncol(edge_weights) == 22,
  length(thresholds) == 22,
  all(is.finite(edge_weights)),
  all(is.finite(thresholds)),
  all(diag(edge_weights) == 0)
)

stopifnot(
  isTRUE(
    all.equal(
      edge_weights,
      t(edge_weights)
    )
  )
)


##############################
# 6. Define primary perturbation
##############################
# Threshold perturbations represent modeled changes in the
# activation propensity of the target symptom.

threshold_sd <- sd(
  thresholds
)

perturbation_magnitude <-
  2 * threshold_sd

cat(
  "Threshold SD:",
  threshold_sd,
  "\n"
)

cat(
  "Primary perturbation magnitude:",
  perturbation_magnitude,
  "\n"
)


##############################
# 7. Simulation settings
##############################

N_SIM <- 5000

SEED_PRIMARY_NIRA <- 2026

set.seed(
  SEED_PRIMARY_NIRA
)


##############################
# 8. Generate one common baseline
##############################

baseline_sample <-
  IsingSampler::IsingSampler(
    n = N_SIM,
    graph = edge_weights,
    thresholds = thresholds,
    responses = c(0L, 1L)
  )

baseline_symptom_count <-
  rowSums(
    baseline_sample
  )

baseline_mean <-
  mean(
    baseline_symptom_count
  )

cat(
  "Common baseline mean modeled symptom count:",
  baseline_mean,
  "\n"
)


##############################
# 9. Define intervention function
##############################

simulate_intervention <- function(
    node,
    direction,
    edge_weights,
    thresholds,
    perturbation_magnitude,
    n_sim
) {
  
  perturbed_thresholds <-
    thresholds
  
  if (direction == "alleviating") {
    
    perturbed_thresholds[node] <-
      thresholds[node] -
      perturbation_magnitude
    
  } else if (direction == "aggravating") {
    
    perturbed_thresholds[node] <-
      thresholds[node] +
      perturbation_magnitude
    
  } else {
    
    stop(
      "direction must be 'alleviating' or 'aggravating'"
    )
  }
  
  simulated_sample <-
    IsingSampler::IsingSampler(
      n = n_sim,
      graph = edge_weights,
      thresholds = perturbed_thresholds,
      responses = c(0L, 1L)
    )
  
  symptom_count <-
    rowSums(
      simulated_sample
    )
  
  mean(
    symptom_count
  )
}


##############################
# 10. Run alleviating interventions
##############################

alleviating_mean <-
  numeric(
    length(thresholds)
  )

for (i in seq_along(thresholds)) {
  
  alleviating_mean[i] <-
    simulate_intervention(
      node = i,
      direction = "alleviating",
      edge_weights = edge_weights,
      thresholds = thresholds,
      perturbation_magnitude =
        perturbation_magnitude,
      n_sim = N_SIM
    )
}


##############################
# 11. Run aggravating interventions
##############################

aggravating_mean <-
  numeric(
    length(thresholds)
  )

for (i in seq_along(thresholds)) {
  
  aggravating_mean[i] <-
    simulate_intervention(
      node = i,
      direction = "aggravating",
      edge_weights = edge_weights,
      thresholds = thresholds,
      perturbation_magnitude =
        perturbation_magnitude,
      n_sim = N_SIM
    )
}


##############################
# 12. Calculate intervention effects
##############################

# Delta < 0:
# decrease in mean modeled number of symptoms present
#
# Delta > 0:
# increase in mean modeled number of symptoms present

alleviating_delta <-
  alleviating_mean -
  baseline_mean

aggravating_delta <-
  aggravating_mean -
  baseline_mean


##############################
# 13. Create result table
##############################

nira_main_results <- data.frame(
  
  Symptom = symptom_codes,
  
  Baseline_Mean =
    rep(
      baseline_mean,
      length(thresholds)
    ),
  
  Alleviating_Mean =
    alleviating_mean,
  
  Alleviating_Delta =
    alleviating_delta,
  
  Aggravating_Mean =
    aggravating_mean,
  
  Aggravating_Delta =
    aggravating_delta
)


##############################
# 14. Rank symptoms
##############################

# More negative = larger alleviating effect

nira_main_results$Alleviating_Rank <-
  rank(
    nira_main_results$Alleviating_Delta,
    ties.method = "min"
  )

# More positive = larger aggravating effect

nira_main_results$Aggravating_Rank <-
  rank(
    -nira_main_results$Aggravating_Delta,
    ties.method = "min"
  )


##############################
# 15. Sorted result tables
##############################

alleviating_results <-
  nira_main_results[
    order(
      nira_main_results$Alleviating_Delta
    ),
  ]

aggravating_results <-
  nira_main_results[
    order(
      -nira_main_results$Aggravating_Delta
    ),
  ]


##############################
# 16. Save results
##############################

write.csv(
  nira_main_results,
  file.path(
    output_dir,
    "NIRA_primary_results_2SD.csv"
  ),
  row.names = FALSE
)

write.csv(
  alleviating_results,
  file.path(
    output_dir,
    "NIRA_alleviating_ranked_2SD.csv"
  ),
  row.names = FALSE
)

write.csv(
  aggravating_results,
  file.path(
    output_dir,
    "NIRA_aggravating_ranked_2SD.csv"
  ),
  row.names = FALSE
)


##############################
# 17. Save simulation settings
##############################

simulation_settings <- data.frame(
  Seed = SEED_PRIMARY_NIRA,
  Number_of_simulated_profiles = N_SIM,
  Threshold_SD = threshold_sd,
  Perturbation_SD_multiplier = 2,
  Perturbation_magnitude =
    perturbation_magnitude,
  Baseline_mean_symptom_count =
    baseline_mean
)

write.csv(
  simulation_settings,
  file.path(
    output_dir,
    "NIRA_primary_simulation_settings.csv"
  ),
  row.names = FALSE
)


##############################
# 18. Save objects
##############################

saveRDS(
  nira_main_results,
  file.path(
    output_dir,
    "NIRA_primary_results_2SD.rds"
  )
)


##############################
# 19. Display top results
##############################

cat(
  "\nTop five alleviating effects:\n"
)

print(
  head(
    alleviating_results[
      ,
      c(
        "Symptom",
        "Alleviating_Delta",
        "Alleviating_Rank"
      )
    ],
    5
  )
)

cat(
  "\nTop five aggravating effects:\n"
)

print(
  head(
    aggravating_results[
      ,
      c(
        "Symptom",
        "Aggravating_Delta",
        "Aggravating_Rank"
      )
    ],
    5
  )
)

cat(
  "\nPrimary ±2 SD NIRA analysis completed successfully.\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)