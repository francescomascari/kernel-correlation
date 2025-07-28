# Kernel-correlation

**R** code to compute kernel correlation with examples and applications.

**Authors:** Marta Catalano, Hugo Lavenant, Francesco Mascari

**Corresponding author:** Francesco Mascari, [francesco.mascari@phd.unibocconi.it](mailto:francesco.mascari@phd.unibocconi.it)


## Overview

This repository accompanies the following paper:
> Catalano, M., Lavenant, H., Mascari, F. (2025+). **Measuring Exchangeability with Reproducing Kernel Hilbert Spaces.** *Preprint*

The abstract of the manuscript is reported below:
> In Bayesian multilevel models, the data are structured in interconnected groups, and their posteriors borrow information from one another due to prior dependence between latent parameters. However, little is known about the behaviour of the dependence a posteriori. In this work, we develop a general framework for measuring partial exchangeability for parametric and nonparametric models, both a priori and a posteriori. We define an index that detects exchangeability, is invariant by reparametrization, can be estimated through samples, and, crucially, is well-suited for posteriors. We achieve these properties through the use of Reproducing Kernel Hilbert Spaces, which map any random probability to a random object on a Hilbert space. This leads to many convenient properties and tractable expressions, especially a priori and under mixing. We apply our general framework to investigate the dependence a posteriori for the hierarchical Dirichlet process, retrieving a parametric convergence rate under very mild assumptions on the data.


## Helper Functions
`R/` contains all **R** files that define helper functions used across scripts.

| File                    | Purpose                              |
| ----------------------- | ------------------------------------ |
| `scripts/01_load.R`     | Loads and inspects raw data          |
| `scripts/02_clean.R`    | Cleans and processes data            |
| `scripts/03_analysis.R` | Performs main statistical analysis   |
| `scripts/04_plot.R`     | Creates visualizations               |
| `R/utils.R`             | Helper functions used across scripts |


## Scripts for simulations
`script/` contains all **R** scripts to perform simulations.
| File                    | Purpose                              |
| ----------------------- | ------------------------------------ |
| `scripts/01_load.R`     | Loads and inspects raw data          |
| `scripts/02_clean.R`    | Cleans and processes data            |
| `scripts/03_analysis.R` | Performs main statistical analysis   |
| `scripts/04_plot.R`     | Creates visualizations               |
| `R/utils.R`             | Helper functions used across scripts |