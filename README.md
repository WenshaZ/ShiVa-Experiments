# ShiVa-Experiments

This repository contains simulation code and case study analyses for the paper  **_Detection of Evolutionary Shifts in Variance under an Ornstein–Uhlenbeck Model_**.

---

## File Overview
### Main Files
- **ShiVa.R**  
  Main implementation of the ShiVa method. It supports joint detection of shifts in trait optima and evolutionary variance with optional cross-validation and model comparison via BIC.

- **case-study.Rmd**  
  R Markdown file for case studies. It applies the implemented methods to real-world datasets (e.g., floral diameter in Rafflesiaceae).

### method_helpers/
- **DefineParameterLimits.R** – Defines parameter limits for use in PCMFit.
- **l1ou_fixed.R** – Fixes an outdated issue in the l1ou package that can cause runtime errors.

### simulation_codes/
- **experiment.R** – Main simulation script for evaluating detection methods across settings.
- **experiment_misestimation.R** – Simulations under misestimation of selection strength α.
- **real_data_simulations.R** – Simulations based on models estimated from real datasets.
- **runfile.csv** – Spreadsheet of parameter settings used as input to `experiment.R`.
- **runfile_misestimation.csv** – Settings for `experiment_misestimation.R`.

---

## Required R Packages

Please install the following R packages before running the code:

```r
install.packages(c(
  "phylolm",
  "PhylogeneticEM",
  "PCMBase",
  "PCMBaseCpp",
  "PCMFit",
  "data.table",
  "phytools",
  "devtools"
))
```

For GitHub packages:

```r
library(devtools)
install_github("glmgen/genlasso")
install_github("khabbazian/l1ou")  # Version of l1ou with variance shift support
```

---

## Usage

- The **case-study.Rmd** script reproduces figures and comparisons for empirical datasets like the flower and sunfish data.
- Run **experiment.R** using `runfile.csv` for main experiments.
- Use **experiment_misestimation.R** with `runfile_misestimation.csv` for sensitivity to α.
- Use **real_data_simulations.R** to simulate data from fitted empirical models.

---

