############################################################
# 12_session_info.R
#
# Record R version, operating system, package versions,
# and full session information for reproducibility.
#
# Output:
#   results/12_session_info/
#     sessionInfo.txt
#     package_versions.csv
#     R_environment_summary.csv
############################################################


##############################
# 1. Output directory
##############################

output_dir <- file.path(
  "results",
  "12_session_info"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


##############################
# 2. Packages used in this study
##############################

analysis_packages <- c(
  "readxl",
  "psych",
  "GPArotation",
  "openxlsx",
  "bootnet",
  "qgraph",
  "networktools",
  "dplyr",
  "tidyr",
  "ggplot2",
  "IsingFit",
  "IsingSampler",
  "nodeIdentifyR",
  "parallel",
  "patchwork",
  "NetworkComparisonTest"
)


##############################
# 3. Check installed packages
##############################

installed <- installed.packages()

package_versions <- data.frame(
  Package = analysis_packages,
  Installed = analysis_packages %in% rownames(installed),
  Version = NA_character_
)

for (i in seq_along(analysis_packages)) {
  
  pkg <- analysis_packages[i]
  
  if (pkg %in% rownames(installed)) {
    
    package_versions$Version[i] <-
      as.character(
        packageVersion(pkg)
      )
  }
}


##############################
# 4. Save package versions
##############################

write.csv(
  package_versions,
  file.path(
    output_dir,
    "package_versions.csv"
  ),
  row.names = FALSE
)


##############################
# 5. Record R environment
##############################

environment_summary <- data.frame(
  Item = c(
    "R_version",
    "Platform",
    "Operating_system",
    "R_architecture",
    "Locale"
  ),
  
  Value = c(
    R.version.string,
    R.version$platform,
    Sys.info()[["sysname"]],
    R.version$arch,
    Sys.getlocale()
  )
)

write.csv(
  environment_summary,
  file.path(
    output_dir,
    "R_environment_summary.csv"
  ),
  row.names = FALSE
)


##############################
# 6. Save full sessionInfo()
##############################

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


##############################
# 7. Save detailed package session info
##############################

capture.output(
  {
    cat(
      "R SESSION INFORMATION\n"
    )
    
    cat(
      "=====================\n\n"
    )
    
    print(
      sessionInfo()
    )
    
    cat(
      "\n\nPACKAGE VERSIONS USED IN THIS STUDY\n"
    )
    
    cat(
      "===================================\n\n"
    )
    
    print(
      package_versions,
      row.names = FALSE
    )
  },
  
  file = file.path(
    output_dir,
    "sessionInfo_with_package_versions.txt"
  )
)


##############################
# 8. Display summary
##############################

cat(
  "\nSession information saved successfully.\n"
)

cat(
  "R version:",
  R.version.string,
  "\n"
)

cat(
  "Operating system:",
  Sys.info()[["sysname"]],
  "\n"
)

cat(
  "Output directory:",
  output_dir,
  "\n"
)

cat(
  "\nPackage versions:\n"
)

print(
  package_versions,
  row.names = FALSE
)