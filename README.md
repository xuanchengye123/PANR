# Estimation and Inference for a Covariate-Dependent Hinge Model: Application to the Ferritin–Hemoglobin Relationship

Sample R code implementing the partitioned alternating Newton–Raphson (PANR) algorithm for Simulation Scenario 3.

## Contents

- `SampleCode_Scenario3.R`: Standalone simulation code for **Scenario 3**, including the data-generating mechanism, the proposed partitioned alternating Newton–Raphson (PANR) estimator, and a Monte Carlo summary table.

## Reproducing Simulation Scenario 3

The script uses **base R only**; on non-Windows systems it can use R's standard `parallel` package, with a serial fallback on Windows. No external data files are required.

Run from R in the directory containing the file:

```r
source("SampleCode_Scenario3.R")
```

Default settings in the supplied script are 1,000 Monte Carlo replications, sample size 1,000, and seed 20260914. The script requests up to 255 parallel workers on non-Windows systems; this may be resource-intensive. Before running on a personal computer, set `N_CORES` to a suitable number, or set it to `1L` for serial execution. Reducing `MC_REPS` is useful for testing, but does not reproduce the original Monte Carlo settings.

The summary reports empirical bias, Monte Carlo standard deviation (MCSD), root mean squared error (RMSE), average estimated standard error (AVESE), and coverage probability (CP) for nominal 95% Wald confidence intervals. The script prints its output but does not save result files.

## Scope

This is **sample code for Simulation Scenario 3 only**. It is not a complete code repository for all simulations or the NHANES data analysis.
