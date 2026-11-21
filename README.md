# Penalized Estimation in Segmented Regression with Covariate-Dependent Change Points

Sample R simulation code accompanying the paper **Penalized Estimation in Segmented Regression with Covariate-Dependent Change Points**.

## Contents

- `sample_simulation.R`: Simulation code comparing SCAD, MCP, lasso, L0, confidence-interval-based selection, and an oracle estimator. It includes three configurable scenarios, variable-selection summaries, and SCAD/MCP inference summaries.
- `SampleCode_Scenario3.R`: Previously included standalone PANR simulation example for Scenario 3.

## Running the sample simulation

Install the required packages once:

```r
install.packages(c("MASS", "ncvreg", "L0Learn"))
```

From the directory containing the script, run:

```r
source("sample_simulation.R")
```

The supplied sample uses Scenario I, sample size 250, 2 replications, and seed 20262027. These settings provide a small example run; they do not reproduce the full Monte Carlo study.

To expand the simulation, edit the settings at the top of `sample_simulation.R`. The commented full-simulation settings select Scenarios I–III, sample sizes 250, 500, 1000, and 2500, and 1,000 replications.

## Output

Results and checkpoints are saved as CSV and RDS files in `simulation_output/` under the current working directory. Outputs include coefficient estimates, method comparisons, change-point variable-selection summaries, and parameter-level and averaged inference summaries for SCAD and MCP.

## Scope

This repository provides sample simulation code. No external data files are required to run `sample_simulation.R`.
