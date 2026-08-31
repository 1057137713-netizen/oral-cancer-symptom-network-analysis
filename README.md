# oral-cancer-symptom-network-analysis
Study-specific R scripts and reproducibility materials for symptom network analysis and simulation-based perturbation analysis in patients with oral cancer.
## Reproducible analysis workflow

This repository contains the study-specific R scripts used to reproduce the analyses reported in the manuscript.

The workflow includes:

- data preparation
- exploratory factor analysis
- Gaussian graphical network estimation
- bridge-strength sensitivity analysis
- symptom dichotomization
- binary Ising-model estimation
- Ising-model evaluation
- primary ±2 SD NIRA simulations
- 1,000-repetition stability analysis
- perturbation-magnitude sensitivity analysis
- Figure 7 generation
- subgroup network estimation and Network Comparison Tests
- package-version and session-information recording


## Data availability

Participant-level clinical data are not publicly available because of ethical and privacy restrictions.

The original dataset should be placed locally at:

```text
data/raw_data.xlsx
```
The participant-level dataset is not included in this repository.
The analysis scripts expect the following symptom variables:
```text
Q1-Q22
```

For subgroup analyses, the dataset additionally requires:
```text
sex
site
way
```

A detailed variable description is provided in:
```text
data_dictionary.csv
```

To support verification of the simulation workflow without access to the participant-level dataset, the repository also provides intermediate Ising-model parameter files, including the estimated edge-weight matrix and symptom-specific threshold parameters.

Scripts 08-10 automatically use these public intermediate parameter files when the locally fitted Ising-model object is unavailable. Therefore, the NIRA perturbation, repeated stability, and perturbation-magnitude sensitivity workflows can be verified without access to the participant-level dataset.
- the network is not re-estimated after threshold perturbation

## Symptom coding
The 22 MDASI-H&N symptom items are scored from 0 to 10.
For the continuous Gaussian graphical network analyses, the original 0-10 symptom scores are used.
For the binary Ising-model and simulation analyses, symptom scores are dichotomized as:
```text
0  = symptom absent
>0 = symptom present
```

## Analysis workflow
The scripts are intended to be run from the project root directory.
### 01. Data preparation
```text
scripts/01_data_preparation.R
```
Imports the original dataset, extracts Q1-Q22, checks variable names, missing values, and the expected 0-10 score range.
Outputs:
```text
results/01_data_preparation/
```


### 02. Exploratory factor analysis

```text
scripts/02_factor_analysis.R
```

Performs exploratory factor analysis of the 22 continuous symptom items.

The analysis includes:

- Kaiser-Meyer-Olkin (KMO) test
- Bartlett's test of sphericity
- parallel analysis
- principal axis factoring
- Promax rotation
- comparison of 2-, 3-, and 4-factor solutions
- factor-loading assignment using |loading| >= 0.30
- Cronbach's alpha calculation

Outputs:

```text
results/02_factor_analysis/
```


### 03. Gaussian graphical network analysis

```text
scripts/03_GGM_network_analysis.R
```

Estimates the continuous symptom network using the original 0-10 symptom scores.

Network estimation settings:

```text
Estimator: EBICglasso
Correlation method: Pearson
EBIC tuning parameter gamma = 0.5
```

The script also produces:

- edge-weight matrix
- raw strength centrality
- raw bridge strength
- nonparametric edge-weight bootstrap
- case-dropping centrality stability analysis
- correlation-stability (CS) coefficients

Outputs:

```text
results/03_GGM_network_analysis/
```


### 04. Bridge-strength sensitivity analysis

```text
scripts/04_bridge_strength_sensitivity.R
```

Evaluates the sensitivity of bridge-strength estimates to the definition of symptom communities.

Bridge strength is compared using:

```text
Four-factor solution
Two-factor solution
```

Outputs:

```text
results/04_bridge_strength_sensitivity/
```


### 05. Binary symptom-data preparation

```text
scripts/05_prepare_binary_data.R
```

Converts the original 0-10 symptom scores into binary symptom indicators for Ising-model estimation and simulation analyses.

Dichotomization rule:

```text
0  = symptom absent
>0 = symptom present
```

Outputs:

```text
results/05_binary_data/
```


### 06. Binary Ising-model estimation

```text
scripts/06_estimate_ising_model.R
```

Estimates the binary Ising model using the dichotomized Q1-Q22 symptom indicators.

The script exports the fitted model parameters required for subsequent simulation analyses, including:

```text
Ising edge-weight matrix
Symptom-specific threshold parameters
```

Outputs:

```text
results/06_Ising_model/
```


### 07. Ising-model evaluation

```text
scripts/07_validate_ising_model.R
```

Evaluates whether the estimated Ising model reproduces the observed marginal symptom activation probabilities.

The evaluation compares observed and simulated activation probabilities across Q1-Q22.

Agreement is summarized using:

```text
Mean absolute difference (MAD)
```

Outputs:

```text
results/07_ising_validation/
```


### 08. Primary ±2 SD NIRA simulation

```text
scripts/08_NIRA_primary_2SD.R
```

Performs the primary simulation-based perturbation analysis.

For each symptom, the target threshold is perturbed individually by:

```text
-2 SD = alleviating simulation
+2 SD = aggravating simulation
```

During each modeled perturbation:

- the Ising edge-weight matrix remains fixed
- thresholds of all other symptoms remain unchanged
- no post-intervention network is re-estimated
- 5,000 binary symptom profiles are generated

The modeled outcome is the mean number of symptoms present.

The modeled perturbation effect is defined as:

```text
Delta = perturbed-condition mean symptom count
        - common baseline mean symptom count
```

Outputs:

```text
results/08_NIRA_primary_2SD/
```


### 09. Repeated NIRA stability analysis

```text
scripts/09_NIRA_stability_1000rep.R
```

Repeats the complete primary ±2 SD simulation procedure 1,000 times.

Within each repetition, one common simulated baseline is used for all 44 modeled perturbation comparisons.

The stability analysis reports:

- mean Delta
- median Delta
- empirical 2.5th-97.5th percentile simulation interval
- mean rank
- median rank
- rank-first frequency
- top-three frequency

The exact repetition-specific random seeds are also exported for reproducibility.

Outputs:

```text
results/09_NIRA_stability/
```


### 10. Perturbation-magnitude sensitivity analysis

```text
scripts/10_NIRA_sensitivity_analysis.R
```

Evaluates whether the direction and overall ranking pattern of the simulated perturbation effects remain similar...

Perturbation settings:

```text
±1 SD
±1.5 SD
±2 SD
```

The same repetition-specific seeds are used across perturbation magnitudes.

The analysis evaluates:

- effect-direction consistency
- Spearman correlations of mean symptom ranks
- top-ranked symptoms across perturbation magnitudes

Outputs:

```text
results/10_NIRA_sensitivity/
```


### 11. Figure 7 generation

```text
scripts/11_generate_figure7.R
```

Generates the final Figure 7 using the 1,000-repetition primary ±2 SD stability results.

In Figure 7:

```text
Points = mean Delta
Horizontal error bars = empirical 2.5th-97.5th percentile simulation interval
Dashed vertical line = Delta = 0
```

Outputs:

```text
results/11_figure7/
```


### 12. Session information

```text
scripts/12_session_info.R
```

Records the software environment used for the final analyses.

The script documents:

- R version
- operating system
- package versions
- complete sessionInfo()

Outputs:

```text
results/12_session_info/
```


### 13. Subgroup networks and Network Comparison Tests

```text
scripts/13_subgroup_network_NCT.R
```

Estimates separate continuous Gaussian graphical networks for the following subgroup comparisons:

```text
Male vs Female
Oral cavity proper vs Oral vestibule
Radiotherapy alone vs Radiotherapy plus other treatments
```

All subgroup networks use the same settings as the main continuous network:

```text
Estimator: EBICglasso
Correlation method: Pearson
EBIC tuning parameter gamma = 0.5
```

Network Comparison Tests are performed using:

```text
1,000 permutations
```

The reported subgroup comparisons include:

- network structure invariance
- global strength invariance

Outputs:

```text
results/13_subgroup_network_NCT/
```




## Recommended execution order
Run the scripts in numerical order:
```text
01_data_preparation.R
02_factor_analysis.R
03_GGM_network_analysis.R
04_bridge_strength_sensitivity.R
05_prepare_binary_data.R
06_estimate_ising_model.R
07_validate_ising_model.R
08_NIRA_primary_2SD.R
09_NIRA_stability_1000rep.R
10_NIRA_sensitivity_analysis.R
11_generate_figure7.R
12_session_info.R
13_subgroup_network_NCT.R
```


## Random seeds
Random seeds are explicitly defined in the scripts involving simulation, permutation, parallel analysis, or bootstrap procedures.
For the repeated NIRA analyses, repetition-specific seeds are additionally exported so that the exact simulation sequence can be reproduced independently of the number of parallel workers.

## Software environment
R and package versions used for the final analyses are recorded in:
```text
session_info/sessionInfo.txt
```
and
```text
session_info/package_versions.csv
```


## Reproducibility note
Because the participant-level dataset cannot be publicly shared, exact reproduction of analyses that require the original patient-level observations requires authorized access to the source data.
However, the repository provides:
- complete study-specific analysis scripts
- variable definitions and coding rules
- random seeds
- package-version information
- full session information
- estimated Ising edge-weight parameters
- estimated Ising threshold parameters
- scripts for generating simulation outputs and plotting data

These materials allow the analytical workflow and the simulation procedures to be independently inspected and verified.

















