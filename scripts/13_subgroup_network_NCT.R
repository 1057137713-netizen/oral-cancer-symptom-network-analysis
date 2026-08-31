############################################################
# 13_subgroup_network_NCT.R
#
# Subgroup symptom network estimation and
# Network Comparison Tests (NCTs)
#
# Subgroup comparisons:
# 1. Sex:
#    Male vs Female
#
# 2. Anatomical site:
#    Oral cavity proper vs Oral vestibule
#
# 3. Treatment modality:
#    Radiotherapy alone vs
#    Radiotherapy plus other treatments
#
# Methods:
# - Original continuous 0-10 symptom scores
# - EBICglasso
# - Pearson correlations
# - EBIC tuning parameter gamma = 0.5
# - NCT with 1,000 permutations
#
# Only subgroup network figures and NCT results are
# generated because subgroup centrality and stability
# analyses are not reported in the final manuscript.
#
# Input:
#   data/raw_data.xlsx
#
# Output:
#   results/13_subgroup_network_NCT/
############################################################


##############################
# 1. Load packages
##############################

library(readxl)
library(bootnet)
library(qgraph)
library(NetworkComparisonTest)


##############################
# 2. File paths
##############################

input_file <- file.path(
  "data",
  "raw_data.xlsx"
)

output_dir <- file.path(
  "results",
  "13_subgroup_network_NCT"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ",
      input_file,
      "\nParticipant-level clinical data are not distributed publicly."
    )
  )
}


##############################
# 3. Import original dataset
##############################

# The subgroup variables are stored in the
# "整理-ALL" worksheet.

data <- read_excel(
  input_file,
  sheet = "整理-ALL"
)


##############################
# 4. Required variables
##############################

symptom_vars <- paste0(
  "Q",
  1:22
)

required_vars <- c(
  "sex",
  "site",
  "way",
  symptom_vars
)

if (!all(
  required_vars %in% names(data)
)) {
  
  missing_vars <- required_vars[
    !required_vars %in% names(data)
  ]
  
  stop(
    paste0(
      "Required variable(s) missing: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  )
}


##############################
# 5. Prepare symptom variables
##############################

data[
  ,
  symptom_vars
] <- lapply(
  data[
    ,
    symptom_vars
  ],
  as.numeric
)

if (any(
  is.na(
    data[
      ,
      symptom_vars
    ]
  )
)) {
  
  stop(
    "Missing values were detected in Q1-Q22."
  )
}

if (any(
  unlist(
    data[
      ,
      symptom_vars
    ]
  ) < 0 |
  unlist(
    data[
      ,
      symptom_vars
    ]
  ) > 10
)) {
  
  stop(
    "Symptom scores outside the expected 0-10 range were detected."
  )
}


##############################
# 6. Define subgroup coding
##############################

# Coding used in the study dataset:
#
# sex:
#   1 = Male
#   2 = Female
#
# site:
#   1 = Oral cavity proper
#   2 = Oral vestibule
#
# way:
#   1 = Radiotherapy alone
#   2 = Radiotherapy plus other treatments


##############################
# 7. Create subgroup datasets
##############################

# Sex

male_data <- data[
  data$sex == 1,
  symptom_vars
]

female_data <- data[
  data$sex == 2,
  symptom_vars
]


# Anatomical site

oral_cavity_data <- data[
  data$site == 1,
  symptom_vars
]

oral_vestibule_data <- data[
  data$site == 2,
  symptom_vars
]


# Treatment modality

RT_alone_data <- data[
  data$way == 1,
  symptom_vars
]

RT_plus_data <- data[
  data$way == 2,
  symptom_vars
]


##############################
# 8. Sample-size checks
##############################

sample_sizes <- data.frame(
  Subgroup = c(
    "Male",
    "Female",
    "Oral cavity proper",
    "Oral vestibule",
    "Radiotherapy alone",
    "Radiotherapy plus other treatments"
  ),
  
  N = c(
    nrow(male_data),
    nrow(female_data),
    nrow(oral_cavity_data),
    nrow(oral_vestibule_data),
    nrow(RT_alone_data),
    nrow(RT_plus_data)
  )
)

print(
  sample_sizes
)

write.csv(
  sample_sizes,
  file.path(
    output_dir,
    "subgroup_sample_sizes.csv"
  ),
  row.names = FALSE
)


##############################
# 9. Estimate subgroup networks
##############################

estimate_subgroup_network <- function(
    subgroup_data
) {
  
  bootnet::estimateNetwork(
    subgroup_data,
    default = "EBICglasso",
    corMethod = "cor",
    tuning = 0.5
  )
}


# Sex networks

male_network <-
  estimate_subgroup_network(
    male_data
  )

female_network <-
  estimate_subgroup_network(
    female_data
  )


# Anatomical-site networks

oral_cavity_network <-
  estimate_subgroup_network(
    oral_cavity_data
  )

oral_vestibule_network <-
  estimate_subgroup_network(
    oral_vestibule_data
  )


# Treatment-modality networks

RT_alone_network <-
  estimate_subgroup_network(
    RT_alone_data
  )

RT_plus_network <-
  estimate_subgroup_network(
    RT_plus_data
  )


##############################
# 10. Obtain common network layout
##############################

# Use the full-sample GGM network, when available,
# to define one common node layout for all subgroup
# figures.

main_network_file <- file.path(
  "results",
  "03_GGM_network_analysis",
  "GGM_network_objects.RData"
)

if (file.exists(main_network_file)) {
  
  load(
    main_network_file
  )
  
  if (!exists("adj_matrix")) {
    stop(
      "adj_matrix was not found in the main GGM results file."
    )
  }
  
  fixed_layout <- qgraph::qgraph(
    adj_matrix,
    layout = "spring",
    DoNotPlot = TRUE
  )$layout
  
} else {
  
  warning(
    paste0(
      "Main GGM object not found. ",
      "A layout will be derived from the full dataset."
    )
  )
  
  overall_network <-
    estimate_subgroup_network(
      data[
        ,
        symptom_vars
      ]
    )
  
  fixed_layout <- qgraph::qgraph(
    overall_network$graph,
    layout = "spring",
    DoNotPlot = TRUE
  )$layout
}


##############################
# 11. Node colors
##############################

# Four-factor symptom structure used in the manuscript

node_colors <- c(
  "#8BC17A",  # Q1
  "#A6CEE3",  # Q2
  "#E9A3A3",  # Q3
  "#A6CEE3",  # Q4
  "#A6CEE3",  # Q5
  "#A6CEE3",  # Q6
  "#A6CEE3",  # Q7
  "#8BC17A",  # Q8
  "#A6CEE3",  # Q9
  "#8BC17A",  # Q10
  "#A6CEE3",  # Q11
  "#E9A3A3",  # Q12
  "#8BC17A",  # Q13
  "#8BC17A",  # Q14
  "#F4A259",  # Q15
  "#F4A259",  # Q16
  "#8BC17A",  # Q17
  "#8BC17A",  # Q18
  "#A6CEE3",  # Q19
  "#8BC17A",  # Q20
  "#8BC17A",  # Q21
  "#8BC17A"   # Q22
)


##############################
# 12. Function to draw subgroup pair
##############################

plot_subgroup_pair <- function(
    network_A,
    network_B,
    title_A,
    title_B,
    output_file
) {
  
  png(
    filename = output_file,
    width = 3200,
    height = 2100,
    res = 300
  )
  
  par(
    mfrow = c(1, 2),
    mar = c(
      1.5,
      1.5,
      3.5,
      1.5
    )
  )
  
  qgraph::qgraph(
    network_A$graph,
    layout = fixed_layout,
    labels = symptom_vars,
    color = node_colors,
    theme = "classic",
    vsize = 7,
    label.cex = 0.95,
    borders = TRUE,
    border.width = 1.2,
    edge.color = ifelse(
      network_A$graph > 0,
      "#2C5EFF",
      "#E64B35"
    ),
    legend = FALSE,
    details = FALSE,
    minimum = 0,
    title = title_A
  )
  
  qgraph::qgraph(
    network_B$graph,
    layout = fixed_layout,
    labels = symptom_vars,
    color = node_colors,
    theme = "classic",
    vsize = 7,
    label.cex = 0.95,
    borders = TRUE,
    border.width = 1.2,
    edge.color = ifelse(
      network_B$graph > 0,
      "#2C5EFF",
      "#E64B35"
    ),
    legend = FALSE,
    details = FALSE,
    minimum = 0,
    title = title_B
  )
  
  dev.off()
}


##############################
# 13. Supplementary Figure S1
#     Sex
##############################

plot_subgroup_pair(
  network_A =
    male_network,
  
  network_B =
    female_network,
  
  title_A =
    paste0(
      "(A) Male (n=",
      nrow(male_data),
      ")"
    ),
  
  title_B =
    paste0(
      "(B) Female (n=",
      nrow(female_data),
      ")"
    ),
  
  output_file = file.path(
    output_dir,
    "Supplementary_Figure_S1_Sex_network.png"
  )
)


##############################
# 14. Supplementary Figure S2
#     Anatomical site
##############################

plot_subgroup_pair(
  network_A =
    oral_cavity_network,
  
  network_B =
    oral_vestibule_network,
  
  title_A =
    paste0(
      "(A) Oral cavity proper (n=",
      nrow(oral_cavity_data),
      ")"
    ),
  
  title_B =
    paste0(
      "(B) Oral vestibule (n=",
      nrow(oral_vestibule_data),
      ")"
    ),
  
  output_file = file.path(
    output_dir,
    "Supplementary_Figure_S2_Anatomical_site_network.png"
  )
)


##############################
# 15. Supplementary Figure S3
#     Treatment modality
##############################

plot_subgroup_pair(
  network_A =
    RT_alone_network,
  
  network_B =
    RT_plus_network,
  
  title_A =
    paste0(
      "(A) Radiotherapy alone (n=",
      nrow(RT_alone_data),
      ")"
    ),
  
  title_B =
    paste0(
      "(B) Radiotherapy plus other treatments (n=",
      nrow(RT_plus_data),
      ")"
    ),
  
  output_file = file.path(
    output_dir,
    "Supplementary_Figure_S3_Treatment_network.png"
  )
)


##############################
# 16. NCT settings
##############################
# Reproducibility seed retained for the revised subgroup workflow
NCT_SEED <- 1234
NCT_PERMUTATIONS <- 1000


##############################
# 17. Sex NCT
##############################

set.seed(
  NCT_SEED
)

sex_NCT <- NetworkComparisonTest::NCT(
  male_data,
  female_data,
  gamma = 0.5,
  it = NCT_PERMUTATIONS,
  binary.data = FALSE,
  paired = FALSE,
  test.edges = FALSE,
  progressbar = TRUE
)


##############################
# 18. Anatomical-site NCT
##############################

set.seed(
  NCT_SEED
)

site_NCT <- NetworkComparisonTest::NCT(
  oral_cavity_data,
  oral_vestibule_data,
  gamma = 0.5,
  it = NCT_PERMUTATIONS,
  binary.data = FALSE,
  paired = FALSE,
  test.edges = FALSE,
  progressbar = TRUE
)


##############################
# 19. Treatment-modality NCT
##############################

set.seed(
  NCT_SEED
)

treatment_NCT <- NetworkComparisonTest::NCT(
  RT_alone_data,
  RT_plus_data,
  gamma = 0.5,
  it = NCT_PERMUTATIONS,
  binary.data = FALSE,
  paired = FALSE,
  test.edges = FALSE,
  progressbar = TRUE
)


##############################
# 20. Extract reported NCT results
##############################

NCT_results <- data.frame(
  
  Comparison = c(
    "Male vs Female",
    "Oral cavity proper vs Oral vestibule",
    paste0(
      "Radiotherapy alone vs ",
      "Radiotherapy plus other treatments"
    )
  ),
  
  Network_invariance_M = c(
    sex_NCT$nwinv.real,
    site_NCT$nwinv.real,
    treatment_NCT$nwinv.real
  ),
  
  Network_invariance_p = c(
    sex_NCT$nwinv.pval,
    site_NCT$nwinv.pval,
    treatment_NCT$nwinv.pval
  ),
  
  Global_strength_S = c(
    sex_NCT$glstrinv.real,
    site_NCT$glstrinv.real,
    treatment_NCT$glstrinv.real
  ),
  
  Global_strength_p = c(
    sex_NCT$glstrinv.pval,
    site_NCT$glstrinv.pval,
    treatment_NCT$glstrinv.pval
  )
)


##############################
# 21. Save Supplementary Table S1
##############################

write.csv(
  NCT_results,
  file.path(
    output_dir,
    "Supplementary_Table_S1_NCT_results.csv"
  ),
  row.names = FALSE
)


##############################
# 22. Save full NCT summaries
##############################

capture.output(
  summary(
    sex_NCT
  ),
  file = file.path(
    output_dir,
    "Sex_NCT_full_summary.txt"
  )
)

capture.output(
  summary(
    site_NCT
  ),
  file = file.path(
    output_dir,
    "Anatomical_site_NCT_full_summary.txt"
  )
)

capture.output(
  summary(
    treatment_NCT
  ),
  file = file.path(
    output_dir,
    "Treatment_NCT_full_summary.txt"
  )
)


##############################
# 23. Save NCT objects
##############################

saveRDS(
  sex_NCT,
  file.path(
    output_dir,
    "Sex_NCT_object.rds"
  )
)

saveRDS(
  site_NCT,
  file.path(
    output_dir,
    "Anatomical_site_NCT_object.rds"
  )
)

saveRDS(
  treatment_NCT,
  file.path(
    output_dir,
    "Treatment_NCT_object.rds"
  )
)


##############################
# 24. Save subgroup networks
##############################

subgroup_networks <- list(
  Male =
    male_network,
  
  Female =
    female_network,
  
  Oral_cavity_proper =
    oral_cavity_network,
  
  Oral_vestibule =
    oral_vestibule_network,
  
  Radiotherapy_alone =
    RT_alone_network,
  
  Radiotherapy_plus_other_treatments =
    RT_plus_network
)

saveRDS(
  subgroup_networks,
  file.path(
    output_dir,
    "subgroup_network_objects.rds"
  )
)


##############################
# 25. Save analysis settings
##############################

analysis_settings <- data.frame(
  Setting = c(
    "Correlation method",
    "Network estimator",
    "EBIC gamma",
    "NCT permutations",
    "NCT random seed",
    "Binary data",
    "Paired comparison"
  ),
  
  Value = c(
    "Pearson",
    "EBICglasso",
    "0.5",
    as.character(
      NCT_PERMUTATIONS
    ),
    as.character(
      NCT_SEED
    ),
    "FALSE",
    "FALSE"
  )
)

write.csv(
  analysis_settings,
  file.path(
    output_dir,
    "subgroup_analysis_settings.csv"
  ),
  row.names = FALSE
)


##############################
# 26. Display results
##############################

cat(
  "\nSubgroup sample sizes:\n"
)

print(
  sample_sizes
)

cat(
  "\nNetwork Comparison Test results:\n"
)

print(
  NCT_results
)

cat(
  "\nSubgroup network and NCT analyses completed successfully.\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)