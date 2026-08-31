############################################################
# 03_GGM_network_analysis.R
#
# Continuous symptom network analysis
#
# Methods:
# - Original 0-10 symptom severity scores
# - Gaussian graphical model (GGM)
# - EBICglasso regularized partial-correlation network
# - Pearson correlations
# - EBIC tuning parameter gamma = 0.5
# - Centrality analysis
# - Bridge strength based on the retained four-factor solution
# - Nonparametric edge bootstrap (1,000)
# - Case-dropping centrality stability bootstrap (1,000)
############################################################


##############################
# 1. Load packages
##############################

library(readxl)
library(bootnet)
library(qgraph)
library(networktools)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)


##############################
# 2. File paths
##############################

data_path <- "data/raw_data.xlsx"

output_path <- "results/03_GGM_network_analysis"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}


##############################
# 3. Random seeds
##############################

SEED_EDGE_BOOTSTRAP <- 20260831
SEED_CENTRALITY_BOOTSTRAP <- 20260901


##############################
# 4. Import and check data
##############################

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Data file not found: ", data_path,
      "\nParticipant-level data are not distributed publicly."
    )
  )
}

data <- read_excel(
  data_path,
  sheet = "整理-ALL"
)

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

if (any(is.na(symptom_data))) {
  stop(
    "Missing values were detected. The study dataset used for the reported analysis contained no missing symptom data."
  )
}

if (any(
  unlist(symptom_data) < 0 |
  unlist(symptom_data) > 10
)) {
  stop("Symptom scores outside the expected 0-10 range were detected.")
}

cat("Sample size:", nrow(symptom_data), "\n")
cat("Number of symptoms:", ncol(symptom_data), "\n")


##############################
# 5. Estimate GGM
##############################

network <- estimateNetwork(
  symptom_data,
  default = "EBICglasso",
  corMethod = "cor",
  tuning = 0.5
)

adj_matrix <- as.matrix(network$graph)
storage.mode(adj_matrix) <- "numeric"

rownames(adj_matrix) <- symptom_vars
colnames(adj_matrix) <- symptom_vars


##############################
# 6. Four-factor symptom clusters
##############################

# Final four-factor solution retained in the manuscript:
#
# Cluster 1: Oral local toxicity
# Q15, Q16
#
# Cluster 2: Chronic psychoneurological
# Q2, Q4, Q5, Q6, Q7, Q9, Q11, Q19
#
# Cluster 3: Head and neck mucosal-physical
# Q1, Q8, Q10, Q13, Q14, Q17, Q18, Q20, Q21, Q22
#
# Cluster 4: Gastrointestinal
# Q3, Q12

groups_list <- list(
  "Cluster 1" = c(15, 16),
  "Cluster 2" = c(2, 4, 5, 6, 7, 9, 11, 19),
  "Cluster 3" = c(1, 8, 10, 13, 14, 17, 18, 20, 21, 22),
  "Cluster 4" = c(3, 12)
)

communities <- c(
  "Cluster3", # Q1
  "Cluster2", # Q2
  "Cluster4", # Q3
  "Cluster2", # Q4
  "Cluster2", # Q5
  "Cluster2", # Q6
  "Cluster2", # Q7
  "Cluster3", # Q8
  "Cluster2", # Q9
  "Cluster3", # Q10
  "Cluster2", # Q11
  "Cluster4", # Q12
  "Cluster3", # Q13
  "Cluster3", # Q14
  "Cluster1", # Q15
  "Cluster1", # Q16
  "Cluster3", # Q17
  "Cluster3", # Q18
  "Cluster2", # Q19
  "Cluster3", # Q20
  "Cluster3", # Q21
  "Cluster3"  # Q22
)

names(communities) <- symptom_vars


##############################
# 7. Draw GGM network figure
##############################

cluster_colors <- c(
  "#F4A261",
  "#8EC7E8",
  "#8BBF6A",
  "#E59AA0"
)

png(
  filename = file.path(
    output_path,
    "Figure3_GGM_network.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

qgraph(
  adj_matrix,
  layout = "spring",
  labels = symptom_vars,
  groups = groups_list,
  color = cluster_colors,
  theme = "classic",
  vsize = 7,
  label.cex = 0.95,
  borders = TRUE,
  border.width = 1.2,
  edge.color = ifelse(
    adj_matrix > 0,
    "#2C5EFF",
    "#E64B35"
  ),
  legend = TRUE,
  details = FALSE,
  minimum = 0
)

dev.off()


##############################
# 8. Export edge-weight matrix
##############################

edge_matrix <- round(
  adj_matrix,
  3
)

diag(edge_matrix) <- 0

write.xlsx(
  as.data.frame(edge_matrix),
  file.path(
    output_path,
    "Supplementary_GGM_edge_weight_matrix.xlsx"
  ),
  rowNames = TRUE
)


##############################
# 9. Export ranked edge table
##############################

edge_table <- data.frame()

for (i in 1:21) {
  
  for (j in (i + 1):22) {
    
    edge_table <- rbind(
      edge_table,
      data.frame(
        Node1 = symptom_vars[i],
        Node2 = symptom_vars[j],
        Weight = adj_matrix[i, j],
        Absolute_weight = abs(adj_matrix[i, j])
      )
    )
  }
}

edge_table <- edge_table[
  order(-edge_table$Absolute_weight),
]

write.xlsx(
  edge_table,
  file.path(
    output_path,
    "GGM_edge_ranking.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 10. Centrality indices
##############################

centrality_results <- qgraph::centralityTable(
  network$graph
)

centrality_long <- centrality_results %>%
  filter(
    measure %in%
      c(
        "Strength",
        "Closeness",
        "Betweenness"
      )
  )

centrality_table <- centrality_long %>%
  select(
    node,
    measure,
    value
  ) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  )

colnames(centrality_table)[1] <- "Symptom"

centrality_table$Strength_rank <- rank(
  -centrality_table$Strength,
  ties.method = "min"
)

centrality_table$Closeness_rank <- rank(
  -centrality_table$Closeness,
  ties.method = "min"
)

centrality_table$Betweenness_rank <- rank(
  -centrality_table$Betweenness,
  ties.method = "min"
)

write.xlsx(
  centrality_table,
  file.path(
    output_path,
    "Centrality_results.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 11. Centrality figure
##############################

centrality_plot <- centrality_results %>%
  filter(
    measure %in%
      c(
        "Strength",
        "Closeness",
        "Betweenness"
      )
  ) %>%
  select(
    node,
    measure,
    value
  )

centrality_plot$measure <- factor(
  centrality_plot$measure,
  levels = c(
    "Strength",
    "Closeness",
    "Betweenness"
  )
)

strength_order <- centrality_plot %>%
  filter(measure == "Strength") %>%
  arrange(desc(value)) %>%
  pull(node)

centrality_plot$node <- factor(
  centrality_plot$node,
  levels = rev(strength_order)
)

p_centrality <- ggplot(
  centrality_plot,
  aes(
    x = value,
    y = node
  )
) +
  geom_point(size = 2.8) +
  geom_line(
    aes(group = 1),
    linewidth = 0.8
  ) +
  facet_wrap(
    ~measure,
    scales = "free_x",
    nrow = 1
  ) +
  theme_bw() +
  labs(
    x = "Raw centrality value",
    y = NULL
  )

ggsave(
  file.path(
    output_path,
    "Figure4_Centrality_raw.png"
  ),
  p_centrality,
  width = 12,
  height = 6,
  dpi = 300
)


##############################
# 12. Bridge strength
##############################

bridge_result <- networktools::bridge(
  adj_matrix,
  communities = communities
)

bridge_table <- data.frame(
  Symptom = names(
    bridge_result$`Bridge Strength`
  ),
  Bridge_strength = as.numeric(
    bridge_result$`Bridge Strength`
  )
)

bridge_table <- bridge_table %>%
  arrange(
    desc(Bridge_strength)
  )

bridge_table$Rank <- seq_len(
  nrow(bridge_table)
)

write.xlsx(
  bridge_table,
  file.path(
    output_path,
    "Bridge_strength_results.xlsx"
  ),
  rowNames = FALSE
)


##############################
# 13. Bridge strength figure
##############################

p_bridge <- ggplot(
  bridge_table,
  aes(
    x = reorder(
      Symptom,
      Bridge_strength
    ),
    y = Bridge_strength
  )
) +
  geom_point(size = 3) +
  geom_line(
    aes(group = 1)
  ) +
  coord_flip() +
  theme_bw() +
  labs(
    x = NULL,
    y = "Bridge strength"
  )

ggsave(
  file.path(
    output_path,
    "Figure5_Bridge_strength.png"
  ),
  p_bridge,
  width = 5,
  height = 7,
  dpi = 300
)


##############################
# 14. Edge-weight accuracy bootstrap
##############################

set.seed(SEED_EDGE_BOOTSTRAP)

boot_edges <- bootnet(
  network,
  nBoots = 1000,
  type = "nonparametric",
  statistics = "edge"
)

png(
  file.path(
    output_path,
    "Figure6A_Edge_bootstrap.png"
  ),
  width = 2600,
  height = 2000,
  res = 300
)

plot(
  boot_edges,
  statistics = "edge",
  plot = "interval",
  labels = TRUE,
  order = "sample",
  meanColor = "black",
  sampleColor = "darkred",
  onlyNonZero = FALSE
)

dev.off()


##############################
# 15. Centrality stability bootstrap
##############################

set.seed(SEED_CENTRALITY_BOOTSTRAP)

boot_centrality <- bootnet(
  network,
  nBoots = 1000,
  type = "case",
  statistics = c(
    "Strength",
    "Closeness",
    "Betweenness"
  )
)

png(
  file.path(
    output_path,
    "Figure6B_Centrality_stability.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

plot(
  boot_centrality,
  statistics = c(
    "Strength",
    "Closeness",
    "Betweenness"
  )
)

dev.off()


##############################
# 16. CS coefficient
##############################

CS_result <- corStability(
  boot_centrality
)

sink(
  file.path(
    output_path,
    "CS_coefficient.txt"
  )
)

print(CS_result)

sink()


##############################
# 17. Save random seeds
##############################

seed_info <- c(
  paste0(
    "Edge bootstrap seed: ",
    SEED_EDGE_BOOTSTRAP
  ),
  paste0(
    "Centrality bootstrap seed: ",
    SEED_CENTRALITY_BOOTSTRAP
  )
)

writeLines(
  seed_info,
  file.path(
    output_path,
    "random_seeds.txt"
  )
)


##############################
# 18. Save analysis objects
##############################

save(
  network,
  adj_matrix,
  edge_matrix,
  edge_table,
  centrality_results,
  centrality_table,
  bridge_result,
  bridge_table,
  boot_edges,
  boot_centrality,
  CS_result,
  file = file.path(
    output_path,
    "GGM_network_objects.RData"
  )
)


##############################
# 19. Completion message
##############################

cat("\nGGM network analysis completed successfully.\n")
cat("Output directory:", output_path, "\n")