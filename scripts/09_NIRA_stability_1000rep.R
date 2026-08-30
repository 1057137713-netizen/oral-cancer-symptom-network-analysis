############################################################
# 09_NIRA_stability_1000rep.R
#
# Repeated stability analysis of the primary ±2 SD NIRA
#
# Purpose:
# Evaluate simulation stability of symptom-level intervention
# effects and rankings across 1,000 complete repetitions.
#
# Key principles:
# 1. Original Ising edge-weight matrix remains fixed.
# 2. One symptom threshold is perturbed at a time.
# 3. Perturbation magnitude = ±2 SD of threshold vector.
# 4. ONE common baseline is generated per repetition.
# 5. The same baseline is used for all 44 intervention
#    comparisons within that repetition.
# 6. Each condition contains 5,000 binary profiles.
# 7. No post-intervention network is re-estimated.
#
# Input:
#   results/06_Ising_model/ising_model.rds
#
# Output:
#   results/09_NIRA_stability/
############################################################


##############################
# 1. Load packages
##############################

library(IsingSampler)
library(parallel)


##############################
# 2. File paths
##############################

model_file <- file.path(
  "results",
  "06_Ising_model",
  "ising_model.rds"
)

output_dir <- file.path(
  "results",
  "09_NIRA_stability"
)

if (!file.exists(model_file)) {
  stop(
    paste0(
      "Input Ising model not found: ",
      model_file,
      "\nPlease run scripts/06_estimate_ising_model.R first."
    )
  )
}

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 3. Load original Ising model
##############################

ising_model <- readRDS(
  model_file
)

edge_weights <- as.matrix(
  ising_model$weiadj
)

thresholds <- as.numeric(
  ising_model$thresholds
)

symptom_codes <- paste0(
  "Q",
  seq_along(thresholds)
)

rownames(edge_weights) <- symptom_codes
colnames(edge_weights) <- symptom_codes
names(thresholds) <- symptom_codes


##############################
# 4. Basic parameter checks
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
# 5. Simulation settings
##############################

N_SIM <- 5000
N_REP <- 1000

MASTER_SEED <- 2026

threshold_sd <- sd(
  thresholds
)

perturbation_magnitude <-
  2 * threshold_sd


##############################
# 6. Generate fixed repetition seeds
##############################

# These seeds make every repetition reproducible,
# independently of the number of parallel workers.

set.seed(
  MASTER_SEED
)

rep_seeds <- sample.int(
  .Machine$integer.max,
  N_REP
)

seed_table <- data.frame(
  Repetition = seq_len(N_REP),
  Seed = rep_seeds
)

write.csv(
  seed_table,
  file.path(
    output_dir,
    "NIRA_repetition_seeds.csv"
  ),
  row.names = FALSE
)


##############################
# 7. Simulate one intervention condition
##############################

simulate_condition <- function(
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
  
  mean(
    rowSums(
      simulated_sample
    )
  )
}


##############################
# 8. Run one complete NIRA repetition
##############################

run_one_nira_rep <- function(
    rep_id,
    seed,
    edge_weights,
    thresholds,
    perturbation_magnitude,
    n_sim
) {
  
  # Repetition-specific deterministic seed
  set.seed(
    seed
  )
  
  
  # ---------------------------
  # Common baseline
  # ---------------------------
  
  baseline_sample <-
    IsingSampler::IsingSampler(
      n = n_sim,
      graph = edge_weights,
      thresholds = thresholds,
      responses = c(0L, 1L)
    )
  
  baseline_mean <-
    mean(
      rowSums(
        baseline_sample
      )
    )
  
  
  # ---------------------------
  # Alleviating simulations
  # ---------------------------
  
  alleviating_mean <-
    numeric(
      length(thresholds)
    )
  
  for (i in seq_along(thresholds)) {
    
    alleviating_mean[i] <-
      simulate_condition(
        node = i,
        direction = "alleviating",
        edge_weights = edge_weights,
        thresholds = thresholds,
        perturbation_magnitude =
          perturbation_magnitude,
        n_sim = n_sim
      )
  }
  
  
  # ---------------------------
  # Aggravating simulations
  # ---------------------------
  
  aggravating_mean <-
    numeric(
      length(thresholds)
    )
  
  for (i in seq_along(thresholds)) {
    
    aggravating_mean[i] <-
      simulate_condition(
        node = i,
        direction = "aggravating",
        edge_weights = edge_weights,
        thresholds = thresholds,
        perturbation_magnitude =
          perturbation_magnitude,
        n_sim = n_sim
      )
  }
  
  
  # ---------------------------
  # Effects vs SAME baseline
  # ---------------------------
  
  alleviating_delta <-
    alleviating_mean -
    baseline_mean
  
  aggravating_delta <-
    aggravating_mean -
    baseline_mean
  
  
  # ---------------------------
  # Rankings
  # ---------------------------
  
  # More negative = larger alleviating effect
  
  alleviating_rank <-
    rank(
      alleviating_delta,
      ties.method = "min"
    )
  
  # More positive = larger aggravating effect
  
  aggravating_rank <-
    rank(
      -aggravating_delta,
      ties.method = "min"
    )
  
  
  # ---------------------------
  # Long-format output
  # ---------------------------
  
  alleviating_results <- data.frame(
    Repetition = rep_id,
    Seed = seed,
    Symptom = symptom_codes,
    Direction = "Alleviating",
    Baseline_Mean = baseline_mean,
    Intervention_Mean =
      alleviating_mean,
    Delta =
      alleviating_delta,
    Rank =
      alleviating_rank
  )
  
  aggravating_results <- data.frame(
    Repetition = rep_id,
    Seed = seed,
    Symptom = symptom_codes,
    Direction = "Aggravating",
    Baseline_Mean = baseline_mean,
    Intervention_Mean =
      aggravating_mean,
    Delta =
      aggravating_delta,
    Rank =
      aggravating_rank
  )
  
  rbind(
    alleviating_results,
    aggravating_results
  )
}


##############################
# 9. Determine parallel workers
##############################

available_cores <-
  parallel::detectCores()

if (is.na(available_cores)) {
  available_cores <- 1
}

N_CORES <- max(
  1,
  min(
    16,
    available_cores - 2
  )
)

cat(
  "Available logical cores:",
  available_cores,
  "\n"
)

cat(
  "Parallel workers used:",
  N_CORES,
  "\n"
)


##############################
# 10. Create parallel cluster
##############################

cl <- parallel::makeCluster(
  N_CORES
)

parallel::clusterEvalQ(
  cl,
  library(IsingSampler)
)

parallel::clusterExport(
  cl,
  varlist = c(
    "edge_weights",
    "thresholds",
    "symptom_codes",
    "perturbation_magnitude",
    "N_SIM",
    "rep_seeds",
    "simulate_condition",
    "run_one_nira_rep"
  ),
  envir = environment()
)


##############################
# 11. Run 1,000 repetitions
##############################

cat(
  "\nStarting 1,000 NIRA repetitions...\n"
)

stability_list <-
  parallel::parLapply(
    cl,
    seq_len(N_REP),
    function(r) {
      
      run_one_nira_rep(
        rep_id = r,
        seed = rep_seeds[r],
        edge_weights = edge_weights,
        thresholds = thresholds,
        perturbation_magnitude =
          perturbation_magnitude,
        n_sim = N_SIM
      )
    }
  )

parallel::stopCluster(
  cl
)

rm(cl)

cat(
  "Parallel simulations completed.\n"
)


##############################
# 12. Combine repetitions
##############################

stability_results <-
  do.call(
    rbind,
    stability_list
  )

rownames(stability_results) <-
  NULL


##############################
# 13. Reproducibility checks
##############################

# Expected:
# 1000 repetitions × 22 symptoms × 2 directions

stopifnot(
  nrow(stability_results) ==
    N_REP * 22 * 2
)

stopifnot(
  length(
    unique(
      stability_results$Repetition
    )
  ) == N_REP
)


# Every repetition must contain exactly 44 rows

rows_per_rep <- table(
  stability_results$Repetition
)

stopifnot(
  all(
    rows_per_rep == 44
  )
)


# Every repetition must use exactly ONE common baseline

baseline_check <- aggregate(
  Baseline_Mean ~ Repetition,
  data = stability_results,
  FUN = function(x) {
    length(
      unique(x)
    )
  }
)

stopifnot(
  all(
    baseline_check$Baseline_Mean == 1
  )
)


# Each repetition must contain one seed only

seed_check <- aggregate(
  Seed ~ Repetition,
  data = stability_results,
  FUN = function(x) {
    length(
      unique(x)
    )
  }
)

stopifnot(
  all(
    seed_check$Seed == 1
  )
)

cat(
  "All reproducibility checks passed.\n"
)


##############################
# 14. Summarize stability
##############################

summary_list <- split(
  stability_results,
  list(
    stability_results$Direction,
    stability_results$Symptom
  ),
  drop = TRUE
)

stability_summary <-
  do.call(
    rbind,
    lapply(
      summary_list,
      function(df) {
        
        data.frame(
          Direction =
            df$Direction[1],
          
          Symptom =
            df$Symptom[1],
          
          Mean_Delta =
            mean(df$Delta),
          
          Median_Delta =
            median(df$Delta),
          
          Delta_2.5 =
            unname(
              quantile(
                df$Delta,
                0.025
              )
            ),
          
          Delta_97.5 =
            unname(
              quantile(
                df$Delta,
                0.975
              )
            ),
          
          Mean_Rank =
            mean(df$Rank),
          
          Median_Rank =
            median(df$Rank),
          
          Rank_First_Frequency =
            mean(
              df$Rank == 1
            ),
          
          Top3_Frequency =
            mean(
              df$Rank <= 3
            )
        )
      }
    )
  )

rownames(stability_summary) <-
  NULL


##############################
# 15. Percentage frequencies
##############################

stability_summary$Rank_First_Percent <-
  100 *
  stability_summary$Rank_First_Frequency

stability_summary$Top3_Percent <-
  100 *
  stability_summary$Top3_Frequency


##############################
# 16. Separate and rank directions
##############################

alleviating_stability <-
  stability_summary[
    stability_summary$Direction ==
      "Alleviating",
  ]

alleviating_stability <-
  alleviating_stability[
    order(
      alleviating_stability$Mean_Delta
    ),
  ]

aggravating_stability <-
  stability_summary[
    stability_summary$Direction ==
      "Aggravating",
  ]

aggravating_stability <-
  aggravating_stability[
    order(
      -aggravating_stability$Mean_Delta
    ),
  ]


##############################
# 17. Focus analysis: Q5, Q9, Q11
##############################

focus_symptoms <- c(
  "Q5",
  "Q9",
  "Q11"
)

focus_results <-
  stability_results[
    stability_results$Symptom %in%
      focus_symptoms,
  ]


focus_stability <-
  stability_summary[
    stability_summary$Symptom %in%
      focus_symptoms,
  ]


##############################
# 18. Ordering of Q5/Q9/Q11
##############################

calculate_focus_order <- function(
    data,
    direction_name
) {
  
  direction_data <-
    data[
      data$Direction ==
        direction_name,
    ]
  
  repetition_ids <-
    sort(
      unique(
        direction_data$Repetition
      )
    )
  
  order_strings <-
    character(
      length(repetition_ids)
    )
  
  for (j in seq_along(repetition_ids)) {
    
    tmp <-
      direction_data[
        direction_data$Repetition ==
          repetition_ids[j],
      ]
    
    if (direction_name ==
        "Alleviating") {
      
      # More negative Delta = stronger effect
      tmp <-
        tmp[
          order(
            tmp$Delta
          ),
        ]
      
    } else {
      
      # More positive Delta = stronger effect
      tmp <-
        tmp[
          order(
            -tmp$Delta
          ),
        ]
    }
    
    order_strings[j] <-
      paste(
        tmp$Symptom,
        collapse = " > "
      )
  }
  
  result <-
    as.data.frame(
      prop.table(
        table(
          order_strings
        )
      )
    )
  
  names(result) <- c(
    "Ordering",
    "Proportion"
  )
  
  result$Percent <-
    100 *
    result$Proportion
  
  result
}


focus_order_alleviating <-
  calculate_focus_order(
    focus_results,
    "Alleviating"
  )

focus_order_aggravating <-
  calculate_focus_order(
    focus_results,
    "Aggravating"
  )


##############################
# 19. Save full results
##############################

write.csv(
  stability_results,
  file.path(
    output_dir,
    "NIRA_stability_all_repetitions.csv"
  ),
  row.names = FALSE
)

write.csv(
  stability_summary,
  file.path(
    output_dir,
    "NIRA_stability_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  alleviating_stability,
  file.path(
    output_dir,
    "NIRA_alleviating_stability_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  aggravating_stability,
  file.path(
    output_dir,
    "NIRA_aggravating_stability_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  focus_stability,
  file.path(
    output_dir,
    "NIRA_Q5_Q9_Q11_stability.csv"
  ),
  row.names = FALSE
)

write.csv(
  focus_order_alleviating,
  file.path(
    output_dir,
    "NIRA_Q5_Q9_Q11_order_alleviating.csv"
  ),
  row.names = FALSE
)

write.csv(
  focus_order_aggravating,
  file.path(
    output_dir,
    "NIRA_Q5_Q9_Q11_order_aggravating.csv"
  ),
  row.names = FALSE
)


##############################
# 20. Save simulation settings
##############################

simulation_settings <- data.frame(
  Master_seed = MASTER_SEED,
  Number_of_repetitions = N_REP,
  Profiles_per_condition = N_SIM,
  Threshold_SD = threshold_sd,
  Perturbation_SD_multiplier = 2,
  Perturbation_magnitude =
    perturbation_magnitude,
  Parallel_workers = N_CORES
)

write.csv(
  simulation_settings,
  file.path(
    output_dir,
    "NIRA_stability_simulation_settings.csv"
  ),
  row.names = FALSE
)


##############################
# 21. Save R object
##############################

saveRDS(
  stability_results,
  file.path(
    output_dir,
    "NIRA_stability_all_repetitions.rds"
  )
)


##############################
# 22. Display key results
##############################

cat(
  "\nTop five alleviating effects:\n"
)

print(
  head(
    alleviating_stability,
    5
  )
)

cat(
  "\nTop five aggravating effects:\n"
)

print(
  head(
    aggravating_stability,
    5
  )
)

cat(
  "\nNIRA stability analysis completed successfully.\n"
)

cat(
  "Results saved in:",
  output_dir,
  "\n"
)