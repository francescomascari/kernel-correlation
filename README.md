# `kernel-correlation`

**R** code to compute kernel correlation with examples and applications.

**Authors:** Marta Catalano, Hugo Lavenant, Francesco Mascari

**Corresponding author:** Francesco Mascari, [francesco.mascari@phd.unibocconi.it](mailto:francesco.mascari@phd.unibocconi.it)


## Overview

This repository accompanies the following paper:
> Catalano, M., Lavenant, H., Mascari, F. (2025+). **Measuring Exchangeability with Reproducing Kernel Hilbert Spaces.** *Preprint*

The abstract of the manuscript is reported below:
> In Bayesian multilevel models, the data are structured in interconnected groups, and their posteriors borrow information from one another due to prior dependence between latent parameters. However, little is known about the behaviour of the dependence a posteriori. In this work, we develop a general framework for measuring partial exchangeability for parametric and nonparametric models, both a priori and a posteriori. We define an index that detects exchangeability for common models, is invariant by reparametrization, can be estimated through samples, and, crucially, is well-suited for posteriors. We achieve these properties through the use of Reproducing Kernel Hilbert Spaces, which map any random probability to a random object on a Hilbert space. This leads to many convenient properties and tractable expressions, especially a priori and under mixing. We apply our general framework to i) investigate the dependence a posteriori for the hierarchical Dirichlet process, retrieving a parametric convergence rate under very mild assumptions on the data; ii) eliciting the dependence structure of a parametric model for a principled comparison with a nonparametric alternative.


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
| `fancy_png.R`             | Renders a plot to a high-quality PNG using LaTeX/TikZ for consistent fonts  |
| `gibbs_tabs.R`            | Does a Gibbs step to update the table allocation in an hDP                  |
| `hdp_corr_anlys.R`        | Computes the kernel correlation for hDP (analytics-based method)            |
| `hdp_corr_smpl.R`         | Computes the kernel correlation for hDP (analytics-based method)            |
| `hdp_mat_sampler_help.R`  | Helper functions to sample from hDP (output in matrix format)               |
| `hdp_mat_sampler.R`       | Samples from hDP (output in matrix format)                                  |
| `hdp_XT_sampler_help.R`   | Helper functions to sample from hDP (output as vector of values)            |
| `hdp_XT_sampler.R`        | Sample from hDP (output as vector of values)                                |
| `hdp_XT2mat.R`            | Converts an hDP as a vector of values in matrix format                        |
| `update_q_probs.R`        | Updates the probabilities for the Gibbs step for table allocation in an hDP |


## Scripts for simulations

`scripts/` contains all **R** scripts to perform simulations.
| File                    | Purpose                                                                                                     |
| ----------------------- | ----------------------------------------------------------------------------------------------------------- |
| `convergence_rate.R`    | Studies the convergence rate of the hDP a posteriori as the size of data increases                            |
| `kernel_stability.R`    | Studies the dependence of the value of kernel correlation on the tuning of kernel parameters                    |
| `other_indx.R`          | Computes the Pearson correlation coefficient and some RKHS-based indices for the Gaussian example                 |
| `par_vs_npar.R`         | Compares borrowing of information between a parametric Gaussian model and an HDP                                  |
| `penguins.R`         | Compares borrowing of information between a parametric Gaussian model and an HDP for the Palmer Penguins data set |
| `smpl_vs_anlys.R`       | Compares the sampling-based method and the analytics-based method to compute the kernel correlation for hDP       |


## Requirements

- **R** (≥ 4.1.0)
- **R packages**
  - Core workflow: `tidyverse`, `doParallel`, `mvtnorm`
  - Math helpers: `emulator`, `matrixStats`
  - Plot export: `tikzDevice`, `pdftools`, `magick`
- **LaTeX toolchain**
  - A working LaTeX installation with `pdftex`.
    > If you do not already have it, the easiest path is the **R** package `tinytex`. Install with `tinytex::install_tinytex()`.


##  Automatic Setup

When opening the project `Kernel-correlation.Rproj` (e.g., via RStudio), the `.Rprofile` automatically:

- Loads essential packages:
  - `tidyverse` (data manipulation)
  - `doParallel` (parallel processing)
  - `mvtnorm` (multivariate normal distribution)
  - `emulator::quad.form`, `quad.3form` (quadratic forms)
  - `matrixStats::logSumExp` (log-space normalization)
  - `tikzDevice`, `pdftools`, `magick` (plot export via LaTeX/TikZ → PDF → PNG)
- Sets up LaTeX/TikZ options for consistent plot rendering.
- Sources all helper scripts in `R/`.

### Note on LaTeX
The `fancy_png()` function requires a working LaTeX installation.  
If LaTeX is not installed, you can quickly install TinyTeX by running in **R**:
```r
  tinytex::install_tinytex()
```

## Reproducibility

All outputs are fully reproducible from the scripts in `scripts/`, and results are regenerated programmatically.


## License

This project is licensed under the [MIT License](LICENSE).
