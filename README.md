# Beyond Overall Diagnosis: A Longitudinal Multi-Domain Latent Structure Model for Domain-Specific Parkinson's Disease Progression
### Zijian Ye, Yuanjia Wang, Yu Gu and Kai Kang* ###


## Overview

This repository contains all code and synthetic data needed to reproduce the simulation results in the paper.


## System Requirements

- **Operating System**: Linux (tested on Ubuntu 24.04 LTS)
- **R**: 4.3.3
- **C++ Compiler**: g++ 13.3.0
- **Rscript** must be available on `PATH` (standard with any R installation; verify with `which Rscript`)
> **Note on reproducibility across platforms**: The simulation results were produced on Linux. Running on Windows or macOS may yield slightly different numerical results due to differences in floating-point arithmetic, thread scheduling, and random number generation across platforms.

### Required R Packages

| Package       | Tested Version | Purpose                            |
| ------------- | -------------- | ---------------------------------- |
| Rcpp          | 1.1.1          | C++ integration                    |
| RcppArmadillo | 15.2.3.1       | Linear algebra in C++              |
| RcppParallel  | 5.1.11.1       | Thread-level parallelism in C++    |
| MCMCpack      | 1.7.1          | Dirichlet/Inverse-Gamma sampling   |
| MASS          | 7.3.60.0.1     | Multivariate normal                |
| dqrng         | 0.4.1          | Fast reproducible RNG              |
| doParallel    | 1.0.17         | Replication-level parallelism      |
| doRNG         | 1.8.6.2        | Reproducible parallel RNG          |
| mclust        | 6.1.2          | Adjusted Rand Index                |
| clue          | 0.3.66         | Label-switching (LSAP)             |
| coda          | 0.19.4.1       | MCMC diagnostics (Gelman-Rubin)    |
| ggplot2       | 4.0.2          | Diagnostic plots                   |
| patchwork     | 1.3.2          | Plot composition                   |
| gridExtra     | 2.3            | Plot arrangement                   |
| dplyr         | 1.1.4          | Data manipulation                  |
| tidyr         | 1.3.2          | Data reshaping                     |
| reshape2      | 1.4.5          | Data reshaping for output figures  |
| vcd           | 1.4.13         | Cramer's V (Suggestion simulation) |

Install all at once:

```r
install.packages(c("Rcpp", "RcppArmadillo", "RcppParallel", "MCMCpack", "MASS",
                    "dqrng", "doParallel", "doRNG", "mclust", "clue", "coda",
                    "ggplot2", "patchwork", "gridExtra", "dplyr", "tidyr", "reshape2", "vcd"))
```

## Directory Structure

- `README.md`  # This file
- `code/`  # Simulation & analysis scripts
   - `run_all.sh`  # Master script: run all simulations + tables/figures
   - **C++ Source**
      - `f_MCMC.cpp`  # Core MCMC sampler
      - `f_MCMC_parallel.cpp`  # Within-chain parallel MCMC (for single runs)
      - `f_MCMC_ModelSelection.cpp`  # MCMC sampler with model selection (AIC/BIC/WAIC)
      - `armspp.hpp`  # Adaptive rejection sampling (header-only)
      - `armspp_parallel.hpp`  # Thread-safe ARS (for f_MCMC_parallel.cpp)
   - **Simulation Scripts**
      - `Simulation_Main.R`  # Main simulation, Settings 1-4 (Table S1 and Figure S13)
      - `Simulation_Addition.R`  # Additional theta settings (Tables S4, S5)
      - `Simulation_Suggestion.R`  # Clustering-prior + multiple MCMC initialization (Table S5)
      - `Simulation_Diagnosis.R`  # MCMC diagnostics (Figures S1, S2)
      - `Simulation_ModelSelection.R`  # Model selection simulation (Table S2)
      - `Simulation_MimicPPMI.R`  # Mimic-PPMI simulation with varying G (Table S3)
      - `Simulation_MimicPPMI_ModelSelection.R`  # Mimic-PPMI model selection summary reported in the text
   - **Table/Figure Scripts**
      - `TableS1_MainSimulation.R`  # Table S1 + Figure S13 (requires all 8 settings)
      - `TableS1_MainSimulation_fast.R`  # Table S1 for available settings only (quick check)
      - `TableS2_ModelSelection.R`  # Table S2
      - `TableS3_MimicPPMI.R`  # Table S3
      - `TableS4_AdditionalRobustness.R`  # Table S4
      - `TableS5_EstimationStrategies.R`  # Table S5
      - `MimicPPMI_ModelSelection_Summary.R`  # Unnumbered summary reported in the text
   - **PPMI Data Files**
      - `Mimic_PPMI_data.RData`  # Processed PPMI data
      - `Fitting_PPMI.RData`  # Fitted PPMI model
- `data/`  # Data & demo scripts
   - `Model_Synthetic.R`  # Demo: run MCMC on CSV data (Table 1, Figure 3, Figure 4)
   - `pre-process.R`  # PPMI data preprocessing (raw PPMI → model input format → MCMC)
   - `synthetic_Y.csv`  # Synthetic longitudinal markers
   - `synthetic_X.csv`  # Synthetic covariates
- `output/`  # All outputs (.Rdata, .pdf, .csv, logs)


## Variable and Parameter Notation

The following variables and parameters are used across the codebase:

| Symbol   | Code variable  | Description                                                       |
| -------- | -------------- | ----------------------------------------------------------------- |
| DS       | `DS`           | Number of replications (Monte Carlo repetitions)                  |
| IR       | `IR`           | Total number of MCMC iterations                                   |
| BI       | `BI`           | Burn-in period (number of initial iterations discarded)           |
| N        | `N`            | Sample size (number of subjects)                                  |
| G        | `G`            | Pre-specified number of latent domains                            |
| K        | `J`            | Number of observed variables (items); denoted K in the paper      |
| B        | `B`            | Number of covariates                                              |
| order_s  | `order_s`      | Order of the B-spline basis                                       |
| Lknots   | `Lknots`       | Number of spline basis functions                                  |
| r_i      | `ri`           | Vector of length N; number of visits (time points) per subject    |
| d_j      | `dj`           | Vector of length J; number of response categories per variable    |
| \Beta    | `beta`         | Regression coefficient for covariate b on domain g (B × G matrix) |
| \sigma^2 | `SigmaG`       | Random-effect variance                                            |
| \Theta   | `Theta`        | Emission probabilities (max(dj) × 2J matrix)                      |
| S        | `s` / `mode_S` | Domain assignment vector of length J                              |
| L        | `L_est`        | Domain structure indicator matrix (J × G)                         |


## Quick Start

### Step 0: Adjust CPU cores

**Important**: Before running any script, set the number of CPU cores to match your machine. The default is 102.

**Set working directory** — simulation and table scripts should be run from the `code/` directory:

```bash
cd /path/to/this/repository/code
```

- In `run_all.sh`: edit `MAX_CORES=102` at the top
- For individual runs: set the environment variable before running:
  ```bash
  export N_CORES=4    # set to the number of cores on your machine
  ```

### Step 1: Quick verification (~1-2 hours on 4 cores)

Run a single simulation setting and check the results:

```bash
export N_CORES=4                         # adjust to your machine!
Rscript Simulation_Main.R 1              # Setting 1, N=200 only
Rscript TableS1_MainSimulation_fast.R    # summarise whatever results are available
```

`TableS1_MainSimulation_fast.R` is a variant of `TableS1_MainSimulation.R` that processes whichever result files exist in `output/` and skips missing ones, so it works even if only one setting has been run.

### Step 2 (optional): Even faster with reduced iterations

For a very quick sanity check, you can reduce the MCMC iterations and number of replications. Edit the top of `code/Simulation_Main.R`:

```r
IR <- 2000   # default: 10000
BI <- 1000   # default: 5000
DS <- 50     # default: 100
```

then run as in Step 1. Note that reducing `DS` will noticeably affect IQR estimates (wider intervals); reducing `IR`/`BI` may affect point estimates slightly.

### Step 3: Run the synthetic data demo

`data/Model_Synthetic.R` runs a single MCMC chain on the provided synthetic CSV data and produces:
- **output/Table1.csv**: Regression coefficient table 
- **output/Figure_3.pdf**: Domain structure 
- **output/Figure_4.pdf**: Emission probability plots
- **output/Fitting_Synthetic.RData**: Compact posterior summary

```bash
cd /path/to/this/repository/data
Rscript Model_Synthetic.R
```

### Step 4: Full reproduction (~1.5 hours on 512 cores)

To reproduce all tables and figures:

```bash
chmod +x code/run_all.sh
cd code
nohup ./run_all.sh > ../output/run_all_master.log 2>&1 &
```

Monitor progress:
```bash
tail -f ../output/run_all_master.log
```

`run_all.sh` schedules all 45 simulation jobs (up to `MAX_JOBS=5` concurrent, each using up to `MAX_CORES` workers) and then runs all table/figure scripts sequentially. Progress is logged to `output/logs/`.

### Step 5: Run individual simulations and their tables

Each simulation script accepts a parameter index. After running, call the
corresponding table script to view results.

**Table S1 and Figure S13** (main simulation):
```bash
for i in $(seq 1 8); do Rscript Simulation_Main.R $i; done
Rscript TableS1_MainSimulation.R
```

**Table S2** (model selection):
```bash
for i in $(seq 1 6); do Rscript Simulation_ModelSelection.R $i; done
Rscript TableS2_ModelSelection.R
```

**Table S3** (Mimic-PPMI, varying G):
```bash
for i in 3 4 5 6 7; do Rscript Simulation_MimicPPMI.R $i; done
Rscript TableS3_MimicPPMI.R
```

**Table S4** (additional theta settings):
```bash
for i in $(seq 1 16); do Rscript Simulation_Addition.R $i; done
Rscript TableS4_AdditionalRobustness.R
```

**Table S5** (estimation strategies):
```bash
for i in $(seq 1 16); do Rscript Simulation_Addition.R $i; done
for i in $(seq 1 8); do Rscript Simulation_Suggestion.R $i; done
Rscript TableS5_EstimationStrategies.R
```

**Mimic-PPMI model-selection summary** (reported in the text; not a numbered table):
```bash
Rscript Simulation_MimicPPMI_ModelSelection.R
Rscript MimicPPMI_ModelSelection_Summary.R
```

**Figures S1, S2** (MCMC diagnostics, output saved as PDF):
```bash
Rscript Simulation_Diagnosis.R
```

> All table scripts print results to the console **and** save to `output/`.
> Remember to set `export N_CORES=<your cores>` before each run.

### Numbered outputs reproduced by this package

The current manuscript numbering is:

| Manuscript item | Output file | Generating script |
| --- | --- | --- |
| Table 1 | `output/Table1.csv` | `data/Model_Synthetic.R` |
| Figure 3 | `output/Figure_3.pdf` | `data/Model_Synthetic.R` |
| Figure 4 | `output/Figure_4.pdf` | `data/Model_Synthetic.R` |
| Table S1 | `output/TableS1.csv` | `code/TableS1_MainSimulation.R` |
| Table S2 | `output/TableS2.csv` | `code/TableS2_ModelSelection.R` |
| Table S3 | `output/TableS3.csv` | `code/TableS3_MimicPPMI.R` |
| Table S4 | `output/TableS4.csv` | `code/TableS4_AdditionalRobustness.R` |
| Table S5 | `output/TableS5.csv` | `code/TableS5_EstimationStrategies.R` |
| Figure S1 | `output/Figure_S1.pdf` | `code/Simulation_Diagnosis.R` |
| Figure S2 | `output/Figure_S2.pdf` | `code/Simulation_Diagnosis.R` |
| Figure S13 | `output/Figure_S13.pdf` | `code/TableS1_MainSimulation.R` |

`output/MimicPPMI_ModelSelection_Summary.csv` reproduces the model-selection result reported in the supplementary text, but it is not a numbered table. Tables S6–S9 and Figures S3–S12 and S14–S16 are based on analyses of the original PPMI data and are not generated by the synthetic-data reproducibility workflow.


### Parameter index mapping

| Script                        | Index | Corresponds to                                                               |
| ----------------------------- | ----- | ---------------------------------------------------------------------------- |
| `Simulation_Main.R`           | 1–8   | Setting × Censoring × N: `{1,2,3,4}` × `{20%,40%}` × `{200,500}` (row-major) |
| `Simulation_Addition.R`       | 1–16  | `{1,4}` × `{20%,40%}` × `{200,500}` × `{theta=1,2}`                          |
| `Simulation_Suggestion.R`     | 1–8   | `{1,4}` × `{20%,40%}` × `{200,500}` (theta=2 only)                           |
| `Simulation_ModelSelection.R` | 1–6   | Setting `{1,2,3}` × N `{100,200}`                                            |
| `Simulation_MimicPPMI.R`      | 3–7   | Number of latent groups G = 3, 4, 5, 6, 7                                    |


## Computational Notes

- **CPU core control**: `run_all.sh` exports an `N_CORES` environment variable (default: `MAX_CORES=102`). Each simulation R script reads this variable to set the number of parallel workers. Set to `1` for fully sequential execution.
- **Parallelism**: Each simulation script uses `doParallel` for replication-level parallelism.
- **Within-chain parallelism**: `f_MCMC_parallel.cpp` and `armspp_parallel.hpp` provide thread-level parallelism within a single MCMC run. Use these in `Model_Synthetic.R` (change `sourceCpp("../code/f_MCMC.cpp")` to `sourceCpp("../code/f_MCMC_parallel.cpp")`) for faster single-run execution.


## Using Your Own Data (PPMI or Other Datasets)

Because the PPMI repository undergoes continuous updates and direct data sharing is prohibited under their Data Use Agreement and Publication Policy, we have generated and provided a synthetic dataset. This ensures the long-term computational reproducibility of our code.

In practice, you can:

1. Run the full pipeline directly with the provided synthetic dataset for reproducible computation.
2. Download PPMI data yourself, preprocess it using our script, and then run the same model workflow.

### Running the Model

`data/Model_Synthetic.R` demonstrates how to run the MCMC model on clinical, longitudinal dataset like PPMI.

**Required input files:**

1. **Response data** (`synthetic_Y.csv` format):
   - Columns: `PATNO` (subject ID), `T_AGE` (visit time), `V1`, `V2`, ..., `VJ` (ordinal responses)
   - Each row is one visit for one subject
   - Missing values coded as `0`

2. **Covariate data** (`synthetic_X.csv` format):
   - Columns: `PATNO` (subject ID), then covariate columns (`SEX`, `EDUCYRS`, `PUTAMEN_R`)
   - One row per subject

**To run:**

```bash
cd /path/to/this/repository/data
Rscript Model_Synthetic.R
```

Adjust `G` (number of clusters) and `IR`/`BI` (iterations/burn-in) at the top of the script.

### Obtaining and Running PPMI Data

The Parkinson's Progression Markers Initiative (PPMI) data can be obtained from:

**https://www.ppmi-info.org/access-data-specimens/download-data**

1. Register for an account and request data access.
2. Download the relevant data files from the following categories:
   - **Motor Assessments**: MDS-UPDRS Part I, II, III items
   - **Non-motor Assessments**: Validation variables
   - **Subject Characteristics**: Demographics, DAT scan results (e.g., SEX, EDUCYRS, PUTAMEN_R)
3. Merge the downloaded tables by `PATNO` (subject ID) and visit identifier to create:
   - A longitudinal markers file (Y) containing `PATNO`, `T_AGE`, and MDS-UPDRS items
   - A subject-level covariate file (X) containing `PATNO` and covariates
4. Use `pre-process.R` to clean and convert the data into model input format, and to run MCMC.

> **PPMI Data Disclaimer**: The PPMI database is continuously evolving and updated over time. Changes to the database content may cause code outputs to differ from the published results.

### Preparing PPMI Data for the Model

`data/pre-process.R` demonstrates the full pipeline for converting raw PPMI data into model input format. The key steps are:

1. Read the merged response file (`Y_data.csv`) and covariate file (`X_data.csv`)
2. Remove subjects with missing covariates
3. Filter to subjects with ≥ 5 visits
4. Handle abnormal/missing response values (code as 0)
5. Select MDS-UPDRS variables of interest (30 items from Parts I–III)
6. Merge sparse response categories (categories with < 5% prevalence)
7. Standardize continuous covariates
8. Build the model input matrices: `Y`, `X`, `tig`, `Yijtc`, `ri`, `dj`
9. Run the MCMC algorithm to get the posterior results
