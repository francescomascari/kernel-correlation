# `kernel-correlation`

**R** code to compute kernel correlation with examples and applications.

**Authors:** Marta Catalano, Hugo Lavenant, Francesco Mascari

**Corresponding author:** Francesco Mascari, [francesco.mascari@phd.unibocconi.it](mailto:francesco.mascari@phd.unibocconi.it)


## Overview

This repository accompanies the following paper:
> Catalano, M., Lavenant, H., Mascari, F. (2025+). **Measuring Exchangeability with Reproducing Kernel Hilbert Spaces.** *Preprint*

The abstract of the manuscript is reported below:
> In Bayesian multilevel models, the data are structured in interconnected groups, and their posteriors borrow information from one another due to prior dependence between latent parameters. However, little is known about the behaviour of the dependence a posteriori. In this work, we develop a general framework for measuring partial exchangeability for parametric and nonparametric models, both a priori and a posteriori. We define an index that detects exchangeability, is invariant by reparametrization, can be estimated through samples, and, crucially, is well-suited for posteriors. We achieve these properties through the use of Reproducing Kernel Hilbert Spaces, which map any random probability to a random object on a Hilbert space. This leads to many convenient properties and tractable expressions, especially a priori and under mixing. We apply our general framework to investigate the dependence a posteriori for the hierarchical Dirichlet process, retrieving a parametric convergence rate under very mild assumptions on the data.


## Repository Structure

```text
`kernel-correlation/`
├── `R/`                        # Core functions and reusable components
├── `scripts/`                  # Main simulation and analysis scripts
├── `output/`                   # All generated results
│   ├── `plots/`                # Plots generated from the simulations
│   └── `results/`              # Simulation output data (CSV format)
├── `.Rprofile`                 # Custom startup settings
├── `.gitignore`                # Git exclusions
├── `LICENSE`                   # Licensing info
├── `README.md`                 # This file
└── `Kernel-correlation.Rproj`  # RStudio project file
```

## Helper Functions

`R/` contains all **R** files that define helper functions used across scripts.

| File                      | Purpose                                                                     |
| ------------------------- | --------------------------------------------------------------------------- |
| `do_outer_mat.R`          | Builds the matrix given by kernel evaluation                                |
| `gibbs_tabs.R`            | Does a Gibbs step to update the table allocation in an hDP                  |
| `hdp_corr_anlys.R`        | Computes the kernel correlation for hDP (analytics-based method)            |
| `hdp_corr_smpl.R`         | Computes the kernel correlation for hDP (analytics-based method)            |
| `hdp_mat_sampler_help.R`  | Helper functions to sample from hDP (output in matrix format)               |
| `hdp_mat_sampler.R`       | Samples from hDP (output in matrix format)                                  |
| `hdp_XT_sampler_help.R`   | Helper functions to sample from hDP (output as vector of values)            |
| `hdp_XT_sampler.R`        | Sample from hDP (output as vector of values)                                |
| `hdp_XT2mat.R`            | Converts an hDP as vector of values in matrix format                        |
| `update_q_probs.R`        | Updates the probabilities for the Gibbs step for table allocation in an hDP |


## Scripts for simulations

`scripts/` contains all **R** scripts to perform simulations.
| File                    | Purpose                                                                                                     |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `convergence_rate.R`    | Studies the converge rate of the hDP a posteriori as the size of data increases                             |
| `kernel_stability.R`    | Studies the dependence of the value of kernel correlation on the tuning of kernel parameters                |
| `par_vs_npar.R`         | Compares borrowing of information between a parametric Gaussian model and an HDP                            |
| `smpl_vs_anlys.R`       | Compares the sampling-based method and the analytics-based method to compute the kernel correlation for hDP |


## Requirements

- **R** (≥ 4.1.0)
- **R packages** : `tidyverse`, `doParallel`, `emulator`, `matrixStats`, `extrafont`
- **System fonts** : [Latin Modern (LM) Roman 10](https://www.ctan.org/tex-archive/fonts/lm)
    > **NOTE:** To register fonts (only once per R installation)
                ```
                extrafont::font_import()
                extrafont::loadfonts()
                ```


##  Automatic Setup

When opening the project `Kernel-correlation.Rproj` (e.g., via RStudio), the `.Rprofile` automatically:

- Loads essential packages:
  - `tidyverse` (data manipulation)
  - `doParallel` (parallel processing)
  - `emulator::quad.form`, `quad.3form` (quadratic forms)
  - `matrixStats::logSumExp` (log-space normalization)
  - `extrafont` (LaTeX-style font support)
- Imports the **Latin Modern font** for plots (ensure it's installed!)
- Sources the scripts in `R/`.

This ensures all functions are loaded and the files in `scripts/` are ready to use.


## Reproducibility

All outputs are fully reproducible from the scripts in `scripts/`, and results are regenerated programmatically.


## License

This project is licensed under the [MIT License](LICENSE).