############################################################
# 10_NIRA_sensitivity_analysis.R
#
# Sensitivity analysis of NIRA perturbation magnitude
#
# Perturbation settings:
#   ±1 SD
#   ±1.5 SD
#   ±2 SD
#
# Strategy:
# 1. Re-run the full repeated NIRA procedure for ±1 SD
#    and ±1.5 SD.
# 2. Reuse the completed ±2 SD stability results.
# 3. Use the exact same 1,000 repetition seeds as the
#    primary ±2 SD stability analysis.
# 4. Compare effect directions and mean symptom rankings
#    across perturbation magnitudes.
#
# Inputs:
#   results/06_Ising_model/ising_model.rds
#   results/09_NIRA_stability/NIRA_stability_all_repetitions.csv
#   results/09_NIRA_stability/NIRA_repetition_seeds.csv
#
# Output:
#   results/10_NIRA_sensitivity/
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

stability_2SD_file <- file.path(
  "results",
  "09_NIRA_stability",
  "NIRA_stability_all_repetitions.csv"
)

seed_file <- file.path(
  "results",
  "09_NIRA_stability",
  "NIRA_repetition_seeds.csv"
)

output_dir <- file.path(
  "results",
  "10_NIRA_sensitivity"
)

required_files <- c(
  model_file,
  stability_2SD_file,
  seed_file
)

if (!all(file.exists(required_files))) {
  stop(
    paste0(
      "One or more required input files are missing.\n",
      "Please run scripts/06_estimate_ising_model.R and ",
      "scripts/09_NIRA_stability_1000rep.R first."
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
# 3. Load Ising model
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
# 4. Basic checks
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

threshold_sd <- sd(
  thresholds
)


##############################
# 6. Load exact repetition seeds
##############################

seed_table <- read.csv(
  seed_file,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(seed_table) == N_REP,
  all(
    seed_table$Repetition ==
      seq_len(N_REP)
  )
)

rep_seeds <- seed_table$Seed


##############################
# 7. Simulate one intervention
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
# 8. Run one sensitivity repetition
##############################

run_one_sensitivity_rep <- function(
    rep_id,
    seed,
    sd_multiplier,
    magnitude_label,
    edge_weights,
    thresholds,
    threshold_sd,
    n_sim
) {
  
  set.seed(
    seed
  )
  
  perturbation_magnitude <-
    sd_multiplier *
    threshold_sd
  
  
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
  # Alleviating conditions
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
  # Aggravating conditions
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
  # Effects
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
  
  alleviating_rank <-
    rank(
      alleviating_delta,
      ties.method = "min"
    )
  
  aggravating_rank <-
    rank(
      -aggravating_delta,
      ties.method = "min"
    )
  
  
  alleviating_results <- data.frame(
    Repetition = rep_id,
    Seed = seed,
    Perturbation = magnitude_label,
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
    Perturbation = magnitude_label,
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
# 9. Parallel runner
##############################

run_parallel_magnitude <- function(
    sd_multiplier,
    magnitude_label
) {
  
  available_cores <-
    parallel::detectCores()
  
  if (is.na(available_cores)) {
    available_cores <- 1
  }
  
  n_cores <- max(
    1,
    min(
      16,
      available_cores - 2
    )
  )
  
  cat(
    "\nRunning ",
    magnitude_label,
    " using ",
    n_cores,
    " workers.\n",
    sep = ""
  )
  
  cl <- parallel::makeCluster(
    n_cores
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
      "threshold_sd",
      "N_SIM",
      "rep_seeds",
      "simulate_condition",
      "run_one_sensitivity_rep"
    ),
    envir = environment()
  )
  
  result_list <-
    parallel::parLapply(
      cl,
      seq_len(N_REP),
      function(r) {
        
        run_one_sensitivity_rep(
          rep_id = r,
          seed = rep_seeds[r],
          sd_multiplier =
            sd_multiplier,
          magnitude_label =
            magnitude_label,
          edge_weights =
            edge_weights,
          thresholds =
            thresholds,
          threshold_sd =
            threshold_sd,
          n_sim =
            N_SIM
        )
      }
    )
  
  parallel::stopCluster(
    cl
  )
  
  results <- do.call(
    rbind,
    result_list
  )
  
  rownames(results) <- NULL
  
  
  ############################
  # Checks
  ############################
  
  stopifnot(
    nrow(results) ==
      N_REP * 22 * 2
  )
  
  rows_per_rep <- table(
    results$Repetition
  )
  
  stopifnot(
    all(
      rows_per_rep == 44
    )
  )
  
  baseline_check <- aggregate(
    Baseline_Mean ~ Repetition,
    data = results,
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
  
  seed_check <- aggregate(
    Seed ~ Repetition,
    data = results,
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
    magnitude_label,
    " completed successfully.\n"
  )
  
  results
}


##############################
# 10. Run ±1 SD
##############################

sensitivity_1SD <-
  run_parallel_magnitude(
    sd_multiplier = 1,
    magnitude_label = "1SD"
  )

write.csv(
  sensitivity_1SD,
  file.path(
    output_dir,
    "NIRA_sensitivity_1SD_all_repetitions.csv"
  ),
  row.names = FALSE
)


##############################
# 11. Run ±1.5 SD
##############################

sensitivity_1_5SD <-
  run_parallel_magnitude(
    sd_multiplier = 1.5,
    magnitude_label = "1.5SD"
  )

write.csv(
  sensitivity_1_5SD,
  file.path(
    output_dir,
    "NIRA_sensitivity_1_5SD_all_repetitions.csv"
  ),
  row.names = FALSE
)


##############################
# 12. Load ±2 SD results
##############################

sensitivity_2SD <- read.csv(
  stability_2SD_file,
  stringsAsFactors = FALSE
)

sensitivity_2SD$Perturbation <-
  "2SD"

# Ensure seed column exists
if (!"Seed" %in% names(sensitivity_2SD)) {
  sensitivity_2SD$Seed <-
    rep_seeds[
      sensitivity_2SD$Repetition
    ]
}

sensitivity_2SD <-
  sensitivity_2SD[
    ,
    c(
      "Repetition",
      "Seed",
      "Perturbation",
      "Symptom",
      "Direction",
      "Baseline_Mean",
      "Intervention_Mean",
      "Delta",
      "Rank"
    )
  ]

stopifnot(
  nrow(sensitivity_2SD) ==
    N_REP * 22 * 2
)


##############################
# 13. Verify matched baselines
##############################

extract_baseline <- function(df) {
  
  tmp <- unique(
    df[
      ,
      c(
        "Repetition",
        "Baseline_Mean"
      )
    ]
  )
  
  tmp <- tmp[
    order(
      tmp$Repetition
    ),
  ]
  
  tmp
}

baseline_1SD <-
  extract_baseline(
    sensitivity_1SD
  )

baseline_1_5SD <-
  extract_baseline(
    sensitivity_1_5SD
  )

baseline_2SD <-
  extract_baseline(
    sensitivity_2SD
  )

stopifnot(
  isTRUE(
    all.equal(
      baseline_1SD$Baseline_Mean,
      baseline_1_5SD$Baseline_Mean
    )
  ),
  isTRUE(
    all.equal(
      baseline_1SD$Baseline_Mean,
      baseline_2SD$Baseline_Mean
    )
  )
)

cat(
  "Matched baseline check passed across all perturbation magnitudes.\n"
)


##############################
# 14. Combine all magnitudes
##############################

sensitivity_all <-
  rbind(
    sensitivity_1SD,
    sensitivity_1_5SD,
    sensitivity_2SD
  )

rownames(
  sensitivity_all
) <- NULL


##############################
# 15. Summarize each symptom
##############################

summary_list <- split(
  sensitivity_all,
  list(
    sensitivity_all$Perturbation,
    sensitivity_all$Direction,
    sensitivity_all$Symptom
  ),
  drop = TRUE
)

sensitivity_summary <-
  do.call(
    rbind,
    lapply(
      summary_list,
      function(df) {
        
        data.frame(
          Perturbation =
            df$Perturbation[1],
          
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
          
          Rank_First_Percent =
            100 *
            mean(
              df$Rank == 1
            ),
          
          Top3_Percent =
            100 *
            mean(
              df$Rank <= 3
            )
        )
      }
    )
  )

rownames(
  sensitivity_summary
) <- NULL


##############################
# 16. Direction consistency
##############################

sensitivity_summary$Effect_Sign <-
  ifelse(
    sensitivity_summary$Mean_Delta < 0,
    "Negative",
    ifelse(
      sensitivity_summary$Mean_Delta > 0,
      "Positive",
      "Zero"
    )
  )

direction_consistency <-
  aggregate(
    Effect_Sign ~ Direction + Symptom,
    data = sensitivity_summary,
    FUN = function(x) {
      length(
        unique(x)
      ) == 1
    }
  )

names(
  direction_consistency
)[3] <- "Direction_Consistent"


##############################
# 17. Overall direction summary
##############################

direction_summary <- data.frame(
  Total_symptom_direction_combinations =
    nrow(direction_consistency),
  
  Direction_consistent =
    sum(
      direction_consistency$Direction_Consistent
    ),
  
  Direction_consistency_percent =
    100 *
    mean(
      direction_consistency$Direction_Consistent
    )
)


##############################
# 18. Spearman correlations
##############################

calculate_rank_correlations <- function(
    summary_df,
    direction_name
) {
  
  sub <- summary_df[
    summary_df$Direction ==
      direction_name,
  ]
  
  ref <- sub[
    sub$Perturbation ==
      "1SD",
    c(
      "Symptom",
      "Mean_Rank"
    )
  ]
  
  symptom_order <-
    ref$Symptom
  
  get_rank <- function(label) {
    
    tmp <- sub[
      sub$Perturbation ==
        label,
      c(
        "Symptom",
        "Mean_Rank"
      )
    ]
    
    tmp$Mean_Rank[
      match(
        symptom_order,
        tmp$Symptom
      )
    ]
  }
  
  rank_1SD <-
    get_rank(
      "1SD"
    )
  
  rank_1_5SD <-
    get_rank(
      "1.5SD"
    )
  
  rank_2SD <-
    get_rank(
      "2SD"
    )
  
  data.frame(
    Direction =
      direction_name,
    
    Comparison = c(
      "1SD vs 1.5SD",
      "1SD vs 2SD",
      "1.5SD vs 2SD"
    ),
    
    Spearman_Rho = c(
      cor(
        rank_1SD,
        rank_1_5SD,
        method = "spearman"
      ),
      cor(
        rank_1SD,
        rank_2SD,
        method = "spearman"
      ),
      cor(
        rank_1_5SD,
        rank_2SD,
        method = "spearman"
      )
    )
  )
}


rank_correlations <-
  rbind(
    calculate_rank_correlations(
      sensitivity_summary,
      "Alleviating"
    ),
    calculate_rank_correlations(
      sensitivity_summary,
      "Aggravating"
    )
  )


##############################
# 19. Top five per setting
##############################

get_top5 <- function(
    summary_df,
    magnitude,
    direction_name
) {
  
  tmp <- summary_df[
    summary_df$Perturbation ==
      magnitude &
      summary_df$Direction ==
      direction_name,
  ]
  
  tmp <- tmp[
    order(
      tmp$Mean_Rank
    ),
  ]
  
  head(
    tmp,
    5
  )
}

top5_results <- do.call(
  rbind,
  lapply(
    c(
      "1SD",
      "1.5SD",
      "2SD"
    ),
    function(mag) {
      
      rbind(
        get_top5(
          sensitivity_summary,
          mag,
          "Alleviating"
        ),
        get_top5(
          sensitivity_summary,
          mag,
          "Aggravating"
        )
      )
    }
  )
)

rownames(top5_results) <- NULL


##############################
# 20. Save outputs
##############################

write.csv(
  sensitivity_all,
  file.path(
    output_dir,
    "NIRA_sensitivity_all_repetitions.csv"
  ),
  row.names = FALSE
)

write.csv(
  sensitivity_summary,
  file.path(
    output_dir,
    "NIRA_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  direction_consistency,
  file.path(
    output_dir,
    "NIRA_direction_consistency.csv"
  ),
  row.names = FALSE
)

write.csv(
  direction_summary,
  file.path(
    output_dir,
    "NIRA_direction_consistency_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  rank_correlations,
  file.path(
    output_dir,
    "NIRA_rank_correlations_across_magnitudes.csv"
  ),
  row.names = FALSE
)

write.csv(
  top5_results,
  file.path(
    output_dir,
    "NIRA_top5_across_magnitudes.csv"
  ),
  row.names = FALSE
)


##############################
# 21. Save settings
##############################

sensitivity_settings <- data.frame(
  Number_of_repetitions = N_REP,
  Profiles_per_condition = N_SIM,
  Threshold_SD = threshold_sd,
  Perturbation_levels =
    "1SD; 1.5SD; 2SD",
  Seed_source =
    "results/09_NIRA_stability/NIRA_repetition_seeds.csv"
)

write.csv(
  sensitivity_settings,
  file.path(
    output_dir,
    "NIRA_sensitivity_settings.csv"
  ),
  row.names = FALSE
)


##############################
# 22. Display key results
##############################

cat(
  "\nDirection consistency:\n"
)

print(
  direction_summary
)

cat(
  "\nSpearman correlations of mean ranks:\n"
)

print(
  rank_correlations
)

cat(
  "\nNIRA sensitivity analysis completed successfully.\n"
)

cat(
  "Results saved in:",
  output_dir,
  "\n"
)