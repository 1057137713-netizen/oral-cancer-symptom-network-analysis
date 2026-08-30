############################################################
# 11_generate_figure7.R
#
# Generate final Figure 7 from the 1,000-repetition
# primary ±2 SD NIRA stability analysis.
#
# Inputs:
#   results/09_NIRA_stability/
#     NIRA_aggravating_stability_summary.csv
#     NIRA_alleviating_stability_summary.csv
#
# Outputs:
#   results/11_figure7/
#     Figure7_simulated_intervention_effects.tiff
#     Figure7_simulated_intervention_effects.pdf
#     Figure7_simulated_intervention_effects.png
#     Figure7_plotting_data.csv
############################################################


##############################
# 1. Required packages
##############################

required_packages <- c(
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(
        missing_packages,
        collapse = ", "
      )
    )
  )
}

library(ggplot2)
library(patchwork)


##############################
# 2. File paths
##############################

input_dir <- file.path(
  "results",
  "09_NIRA_stability"
)

output_dir <- file.path(
  "results",
  "11_figure7"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}

aggravating_file <- file.path(
  input_dir,
  "NIRA_aggravating_stability_summary.csv"
)

alleviating_file <- file.path(
  input_dir,
  "NIRA_alleviating_stability_summary.csv"
)

if (!file.exists(aggravating_file)) {
  stop(
    paste0(
      "Input file not found: ",
      aggravating_file,
      "\nPlease run scripts/09_NIRA_stability_1000rep.R first."
    )
  )
}

if (!file.exists(alleviating_file)) {
  stop(
    paste0(
      "Input file not found: ",
      alleviating_file,
      "\nPlease run scripts/09_NIRA_stability_1000rep.R first."
    )
  )
}


##############################
# 3. Read stability summaries
##############################

aggravating <- read.csv(
  aggravating_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

alleviating <- read.csv(
  alleviating_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


##############################
# 4. Check required columns
##############################

required_cols <- c(
  "Symptom",
  "Mean_Delta",
  "Delta_2.5",
  "Delta_97.5"
)

stopifnot(
  all(
    required_cols %in%
      names(aggravating)
  ),
  all(
    required_cols %in%
      names(alleviating)
  ),
  nrow(aggravating) == 22,
  nrow(alleviating) == 22
)


##############################
# 5. Symptom labels
##############################

symptom_labels <- c(
  Q1  = "Pain",
  Q2  = "Fatigue",
  Q3  = "Nausea",
  Q4  = "Disturbed sleep",
  Q5  = "Distress",
  Q6  = "Shortness of breath",
  Q7  = "Forgetfulness",
  Q8  = "Lack of appetite",
  Q9  = "Drowsiness",
  Q10 = "Dry mouth",
  Q11 = "Sadness",
  Q12 = "Vomiting",
  Q13 = "Numbness or tingling",
  Q14 = "Mouth or throat mucus",
  Q15 = "Difficulty swallowing or chewing",
  Q16 = "Choking",
  Q17 = "Difficulty speaking",
  Q18 = "Skin pain",
  Q19 = "Constipation",
  Q20 = "Taste changes",
  Q21 = "Mouth or throat soreness",
  Q22 = "Teeth or gum problems"
)

stopifnot(
  all(
    aggravating$Symptom %in%
      names(symptom_labels)
  ),
  all(
    alleviating$Symptom %in%
      names(symptom_labels)
  )
)

aggravating$Symptom_Name <-
  unname(
    symptom_labels[
      aggravating$Symptom
    ]
  )

alleviating$Symptom_Name <-
  unname(
    symptom_labels[
      alleviating$Symptom
    ]
  )


##############################
# 6. Check expected directions
##############################

stopifnot(
  all(
    aggravating$Mean_Delta > 0
  ),
  all(
    alleviating$Mean_Delta < 0
  )
)


##############################
# 7. Sort separately
##############################

# Aggravating:
# more positive Delta = larger modeled increase

aggravating <- aggravating[
  order(
    -aggravating$Mean_Delta
  ),
]

# Alleviating:
# more negative Delta = larger modeled reduction

alleviating <- alleviating[
  order(
    alleviating$Mean_Delta
  ),
]


##############################
# 8. Factor levels
##############################

aggravating$Symptom_Name <- factor(
  aggravating$Symptom_Name,
  levels = rev(
    aggravating$Symptom_Name
  )
)

alleviating$Symptom_Name <- factor(
  alleviating$Symptom_Name,
  levels = rev(
    alleviating$Symptom_Name
  )
)


##############################
# 9. Plot colors
##############################

aggravating_fill <- "#D9793A"
aggravating_edge <- "#B85E28"

alleviating_fill <- "#3D918C"
alleviating_edge <- "#2C6F6B"

reference_color <- "grey45"


##############################
# 10. Shared figure theme
##############################

figure_theme <-
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11.5,
      hjust = 0
    ),
    
    axis.text.y = element_text(
      size = 8.4,
      colour = "black"
    ),
    
    axis.text.x = element_text(
      size = 8.5,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 9.5,
      colour = "black",
      margin = margin(
        t = 7
      )
    ),
    
    axis.title.y =
      element_blank(),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    
    panel.grid =
      element_blank(),
    
    plot.margin = margin(
      t = 7,
      r = 10,
      b = 7,
      l = 7
    )
  )


##############################
# 11. Panel A: Aggravating
##############################

pA <- ggplot(
  aggravating,
  aes(
    x = Mean_Delta,
    y = Symptom_Name
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = reference_color
  ) +
  
  geom_errorbar(
    aes(
      xmin = `Delta_2.5`,
      xmax = `Delta_97.5`
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.60,
    colour = aggravating_fill
  ) +
  
  geom_point(
    shape = 21,
    size = 2.7,
    stroke = 0.7,
    fill = aggravating_fill,
    colour = aggravating_edge
  ) +
  
  scale_x_continuous(
    limits = c(
      -0.02,
      0.65
    ),
    breaks = c(
      0.0,
      0.2,
      0.4,
      0.6
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  
  labs(
    title =
      "(A) Aggravating interventions",
    x =
      "Mean \u0394 in modeled symptom count"
  ) +
  
  figure_theme


##############################
# 12. Panel B: Alleviating
##############################

pB <- ggplot(
  alleviating,
  aes(
    x = Mean_Delta,
    y = Symptom_Name
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    colour = reference_color
  ) +
  
  geom_errorbar(
    aes(
      xmin = `Delta_2.5`,
      xmax = `Delta_97.5`
    ),
    orientation = "y",
    width = 0.14,
    linewidth = 0.60,
    colour = alleviating_fill
  ) +
  
  geom_point(
    shape = 21,
    size = 2.7,
    stroke = 0.7,
    fill = alleviating_fill,
    colour = alleviating_edge
  ) +
  
  scale_x_continuous(
    limits = c(
      -0.85,
      0.02
    ),
    breaks = c(
      -0.8,
      -0.6,
      -0.4,
      -0.2,
      0.0
    ),
    expand = expansion(
      mult = c(
        0.01,
        0
      )
    )
  ) +
  
  labs(
    title =
      "(B) Alleviating interventions",
    x =
      "Mean \u0394 in modeled symptom count"
  ) +
  
  figure_theme


##############################
# 13. Combine Figure 7
##############################

figure7 <- pA + pB +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(
      1,
      1
    )
  )

print(
  figure7
)


##############################
# 14. Save Figure 7
##############################

ggsave(
  filename = file.path(
    output_dir,
    "Figure7_simulated_intervention_effects.tiff"
  ),
  plot = figure7,
  width = 13.2,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7_simulated_intervention_effects.pdf"
  ),
  plot = figure7,
  width = 13.2,
  height = 7.2,
  units = "in",
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure7_simulated_intervention_effects.png"
  ),
  plot = figure7,
  width = 13.2,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)


##############################
# 15. Export plotting data
##############################

plot_data_A <- data.frame(
  Panel = "Aggravating",
  
  Symptom = as.character(
    aggravating$Symptom
  ),
  
  Symptom_Name = as.character(
    aggravating$Symptom_Name
  ),
  
  Mean_Delta =
    aggravating$Mean_Delta,
  
  Lower_95_Simulation_Interval =
    aggravating$`Delta_2.5`,
  
  Upper_95_Simulation_Interval =
    aggravating$`Delta_97.5`
)

plot_data_B <- data.frame(
  Panel = "Alleviating",
  
  Symptom = as.character(
    alleviating$Symptom
  ),
  
  Symptom_Name = as.character(
    alleviating$Symptom_Name
  ),
  
  Mean_Delta =
    alleviating$Mean_Delta,
  
  Lower_95_Simulation_Interval =
    alleviating$`Delta_2.5`,
  
  Upper_95_Simulation_Interval =
    alleviating$`Delta_97.5`
)

figure7_plotting_data <- rbind(
  plot_data_A,
  plot_data_B
)

write.csv(
  figure7_plotting_data,
  file.path(
    output_dir,
    "Figure7_plotting_data.csv"
  ),
  row.names = FALSE
)


##############################
# 16. Final checks
##############################

stopifnot(
  nrow(
    figure7_plotting_data
  ) == 44
)

stopifnot(
  aggravating$Symptom[1] ==
    "Q12"
)

stopifnot(
  alleviating$Symptom[1] ==
    "Q1"
)

cat(
  "\nFigure 7 generated successfully.\n"
)

cat(
  "Top aggravating symptom:",
  aggravating$Symptom[1],
  "\n"
)

cat(
  "Top alleviating symptom:",
  alleviating$Symptom[1],
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)