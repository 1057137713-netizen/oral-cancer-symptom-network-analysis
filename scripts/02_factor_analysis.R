############################################################
# 02_factor_analysis.R
#
# Exploratory factor analysis of 22 MDASI-H&N symptoms
#
# Methods:
# - Kaiser-Meyer-Olkin (KMO) test
# - Bartlett's test of sphericity
# - Parallel analysis
# - Principal axis factoring
# - Promax rotation
# - Comparison of 2-, 3-, and 4-factor solutions
# - Factor loading threshold: |loading| >= 0.30
# - Cronbach's alpha for retained symptom clusters
############################################################


##############################
# 1. Load packages
##############################

library(readxl)
library(psych)
library(GPArotation)
library(openxlsx)


##############################
# 2. File paths
##############################

# Participant-level data are not included in the public repository
# because of ethical and privacy restrictions.
#
# Authorized users may place the original data file at:
# data/raw_data.xlsx

data_path <- "data/raw_data.xlsx"

output_path <- "results/02_factor_analysis"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}


##############################
# 3. Random seed
##############################

# Prespecified seed for parallel analysis
SEED_PARALLEL <- 20260830
set.seed(SEED_PARALLEL)


##############################
# 4. Import data
##############################

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Data file not found: ", data_path,
      "\nParticipant-level data are not distributed publicly."
    )
  )
}

data <- read_excel(data_path)


##############################
# 5. Select symptom variables
##############################

symptom_vars <- paste0("Q", 1:22)

if (!all(symptom_vars %in% names(data))) {
  stop("One or more required symptom variables (Q1-Q22) are missing.")
}

symptom_data <- data[, symptom_vars]
symptom_data <- as.data.frame(symptom_data)

symptom_data[] <- lapply(
  symptom_data,
  as.numeric
)


##############################
# 6. Basic data checks
##############################

cat("Sample size:", nrow(symptom_data), "\n")
cat("Number of symptoms:", ncol(symptom_data), "\n")

cat("\nMissing values by symptom:\n")
print(colSums(is.na(symptom_data)))

if (any(is.na(symptom_data))) {
  stop(
    "Missing values were detected. The study dataset used for the reported analyses contained no missing symptom data."
  )
}

if (any(
  unlist(symptom_data) < 0 |
  unlist(symptom_data) > 10
)) {
  stop("Symptom scores outside the expected 0-10 range were detected.")
}


##############################
# 7. KMO and Bartlett tests
##############################

KMO_result <- psych::KMO(symptom_data)

bartlett_result <- psych::cortest.bartlett(
  cor(symptom_data),
  n = nrow(symptom_data)
)

KMO_Bartlett <- data.frame(
  KMO = round(KMO_result$MSA, 3),
  Chi_square = round(bartlett_result$chisq, 2),
  df = bartlett_result$df,
  p_value = bartlett_result$p.value
)

write.xlsx(
  KMO_Bartlett,
  file.path(
    output_path,
    "KMO_Bartlett_results.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 8. Parallel analysis
##############################

set.seed(SEED_PARALLEL)

png(
  filename = file.path(
    output_path,
    "Parallel_analysis.png"
  ),
  width = 1800,
  height = 1400,
  res = 300
)

parallel_result <- psych::fa.parallel(
  symptom_data,
  fa = "fa",
  fm = "pa",
  n.iter = 1000,
  main = "Parallel analysis for exploratory factor analysis"
)

dev.off()


##############################
# 9. Traditional scree plot
##############################

eigenvalues <- eigen(
  cor(symptom_data)
)$values

png(
  filename = file.path(
    output_path,
    "Scree_plot.png"
  ),
  width = 1800,
  height = 1400,
  res = 300
)

plot(
  eigenvalues,
  type = "b",
  xlab = "Factor number",
  ylab = "Eigenvalue",
  main = "Scree plot"
)

abline(
  h = 1,
  lty = 2
)

dev.off()


##############################
# 10. EFA function
##############################

run_EFA <- function(n_factors) {
  
  psych::fa(
    symptom_data,
    nfactors = n_factors,
    fm = "pa",
    rotate = "promax",
    scores = "regression"
  )
  
}


##############################
# 11. Estimate 2-, 3-, and 4-factor solutions
##############################

efa_2 <- run_EFA(2)
efa_3 <- run_EFA(3)
efa_4 <- run_EFA(4)


##############################
# 12. Export complete loading matrices
##############################

export_loadings <- function(model, model_name) {
  
  loading_matrix <- as.data.frame(
    unclass(model$loadings)
  )
  
  loading_matrix$Symptom <- rownames(loading_matrix)
  
  loading_matrix <- loading_matrix[
    ,
    c(
      "Symptom",
      setdiff(names(loading_matrix), "Symptom")
    )
  ]
  
  write.xlsx(
    loading_matrix,
    file.path(
      output_path,
      paste0(
        model_name,
        "_factor_loadings.xlsx"
      )
    ),
    rowNames = FALSE
  )
}


export_loadings(efa_2, "2_factor")
export_loadings(efa_3, "3_factor")
export_loadings(efa_4, "4_factor")


##############################
# 13. Assign symptoms to factors
##############################

# Rule used in the manuscript:
# 1. Retain factor loadings with absolute value >= 0.30.
# 2. For symptoms loading on more than one factor,
#    assign the symptom to the factor with the largest
#    absolute loading.
# 3. Symptoms with no loading >= 0.30 are left unassigned.

assign_factors <- function(model, model_name) {
  
  L <- as.matrix(
    unclass(model$loadings)
  )
  
  result <- data.frame(
    Model = model_name,
    Symptom = rownames(L),
    Assigned_factor = NA_integer_,
    Assigned_loading = NA_real_
  )
  
  for (i in seq_len(nrow(L))) {
    
    eligible <- which(
      abs(L[i, ]) >= 0.30
    )
    
    if (length(eligible) > 0) {
      
      best <- eligible[
        which.max(abs(L[i, eligible]))
      ]
      
      result$Assigned_factor[i] <- best
      result$Assigned_loading[i] <- L[i, best]
    }
  }
  
  result
}


assignment_2 <- assign_factors(
  efa_2,
  "2-factor"
)

assignment_3 <- assign_factors(
  efa_3,
  "3-factor"
)

assignment_4 <- assign_factors(
  efa_4,
  "4-factor"
)

factor_assignments <- rbind(
  assignment_2,
  assignment_3,
  assignment_4
)

write.xlsx(
  factor_assignments,
  file.path(
    output_path,
    "Factor_solution_comparison.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 14. Variance explained
##############################

extract_variance <- function(model, model_name) {
  
  variance <- as.data.frame(
    t(model$Vaccounted)
  )
  
  variance$Model <- model_name
  
  variance <- variance[
    ,
    c(
      "Model",
      setdiff(names(variance), "Model")
    )
  ]
  
  variance
}


variance_table <- rbind(
  extract_variance(
    efa_2,
    "2-factor"
  ),
  extract_variance(
    efa_3,
    "3-factor"
  ),
  extract_variance(
    efa_4,
    "4-factor"
  )
)

write.xlsx(
  variance_table,
  file.path(
    output_path,
    "Variance_explained_comparison.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 15. Cronbach's alpha
##############################

calculate_alpha <- function(
    assignment_table,
    n_factors,
    model_name
) {
  
  result <- data.frame()
  
  for (i in seq_len(n_factors)) {
    
    items <- assignment_table$Symptom[
      assignment_table$Assigned_factor == i
    ]
    
    items <- items[!is.na(items)]
    
    if (length(items) >= 2) {
      
      alpha_value <- psych::alpha(
        symptom_data[, items, drop = FALSE],
        warnings = FALSE
      )$total$raw_alpha
      
      temp <- data.frame(
        Model = model_name,
        Factor = paste0("Factor ", i),
        Items = paste(items, collapse = ", "),
        Number_of_items = length(items),
        Cronbach_alpha = round(alpha_value, 3)
      )
      
      result <- rbind(
        result,
        temp
      )
      
    } else {
      
      warning(
        paste(
          model_name,
          "- Factor",
          i,
          "contains fewer than two retained symptoms."
        )
      )
    }
  }
  
  result
}


alpha_results <- rbind(
  
  calculate_alpha(
    assignment_2,
    2,
    "2-factor"
  ),
  
  calculate_alpha(
    assignment_3,
    3,
    "3-factor"
  ),
  
  calculate_alpha(
    assignment_4,
    4,
    "4-factor"
  )
)

write.xlsx(
  alpha_results,
  file.path(
    output_path,
    "Cronbach_alpha_results.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 16. Save analysis objects
##############################

save(
  efa_2,
  efa_3,
  efa_4,
  parallel_result,
  factor_assignments,
  variance_table,
  alpha_results,
  file = file.path(
    output_path,
    "EFA_results_objects.RData"
  )
)


##############################
# 17. Save seed information
##############################

writeLines(
  paste0(
    "Parallel analysis random seed: ",
    SEED_PARALLEL
  ),
  file.path(
    output_path,
    "random_seed.txt"
  )
)


##############################
# 18. Completion message
##############################

cat("\nFactor analysis completed successfully.\n")
cat("Output directory:", output_path, "\n")