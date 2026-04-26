# ShiVa-Experiments

This repository contains simulation code and case study analyses for the paper **_Detection of Evolutionary Shifts in Variance under an Ornstein–Uhlenbeck Model_**.

The ShiVa method itself is distributed as an R package at <https://github.com/WenshaZ/ShiVa>; the scripts here load it via `library(ShiVa)`. Install with `remotes::install_github("WenshaZ/ShiVa")`.

---

## File Overview

### Top level
- **case-study.Rmd** – R Markdown file applying ShiVa and competing methods to the empirical datasets used in the paper.

### method_helpers/
- **DefineParameterLimits.R** – Parameter limits used by PCMFit.
- **l1ou_fixed.R** – A lightly patched copy of the `l1ou` source. The upstream package is no longer actively maintained on CRAN and a few internal calls no longer match the current versions of its dependencies; the patched copy makes only the minimum modifications needed to keep the simulations runnable. It additionally provides `estimate_shift_configuration_known_alpha`, an alpha-fixed variant of `l1ou::estimate_shift_configuration` that the upstream package does not expose.

### simulation_codes/
- **experiment.R** – Main simulation script. For each row of `runfile.csv`, runs ShiVa, l1ou (pBIC and BIC), phyloEM, and PCMFit on the lizard phylogeny, and records detection rates, predictive log-likelihood, computation time, and Mean Integrated Squared Error (MISE) for both the diffusion variance σ² and the optimal value θ along the tree.
- **experiment_misestimation.R** – Same protocol as `experiment.R`, but ShiVa, l1ou and phyloEM are fit with a misestimated `alpha_hat` (specified in the runfile) and re-fit with the true α for the loglik / MISE evaluation; PCMFit estimates its own α.
- **real_data_simulations.R** – Simulations based on models estimated from the empirical datasets.
- **runfile.csv**, **runfile_misestimation.csv** – Parameter settings consumed by the corresponding scripts.

---

## Required R Packages

```r
install.packages(c(
  "phylolm",
  "PhylogeneticEM",
  "PCMBase",
  "PCMBaseCpp",
  "PCMFit",
  "data.table",
  "phytools",
  "MASS",
  "ape",
  "remotes"
))
```

For GitHub-installed packages:

```r
remotes::install_github("WenshaZ/ShiVa")          # the ShiVa method
remotes::install_github("glmgen/genlasso")
remotes::install_github("khabbazian/l1ou")        # upstream l1ou; see method_helpers/l1ou_fixed.R for compatibility patches
```

---

## Usage

Each simulation script expects an integer index `k` as a single command-line argument; row `ceiling(k/20)` of the corresponding runfile defines the parameter setting, and the script writes results into a `results/` subdirectory of the working directory. Run with, e.g.:

```bash
cd simulation_codes
mkdir -p results
Rscript experiment.R 1
Rscript experiment_misestimation.R 1
```

For the case studies, knit `case-study.Rmd`.

---

## Outputs

Each `experiment*.R` saves an `.rds` file with the following components:

- `ShiVa`, `l1ou`, `phyloEM`, `PCMFit` – per-replicate indicator vectors of detected shift positions
- `loglik` – predictive log-likelihood of each method against the true model on independent test data
- `compute_time` – wall-clock per method
- `MISE` – Mean Integrated Squared Error of the per-branch diffusion variance σ² (ShiVa, l1ou, phyloEM, PCMFit)
- `MISE_mu` – Mean Integrated Squared Error of the per-branch optimal value θ (ShiVa, l1ou, phyloEM, PCMFit; PCMFit's column will be `NA` whenever PCMFit selects a model type without an explicit Theta parameter)
