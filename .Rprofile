## PACKAGES

# data manipulation
require("tidyverse")

# parallelization
require("doParallel")

# multivariate normals
require("mvtnorm")

# quadratic forms
require("emulator", include.only = c("quad.form", "quad.3form"))

# log-space normalization
require("matrixStats", include.only = "logSumExp")

# plot export stack (LaTeX/TikZ -> PDF -> PNG)
require("tikzDevice")
require("pdftools")
require("magick")


## LATEX SETUP

# if no LaTeX installation with pdftex is present, install TinyTeX with:
# tinytex::install_tinytex()

# TikZ/LaTeX options 
options(
  tikzDefaultEngine = "pdftex",
  tikzLatexPackages = c(
    getOption("tikzLatexPackages"),
    "\\usepackage[T1]{fontenc}",
    "\\usepackage[utf8]{inputenc}",
    "\\usepackage{amsmath,amssymb}"
  )
)


## PROJECT HELPERS

# fancy plot export
source("R/fancy_png.R")

# hDP sampling helpers
source("R/hdp_mat_sampler_help.R")
source("R/hdp_mat_sampler.R")
source("R/hdp_XT_sampler_help.R")
source("R/hdp_XT_sampler.R")
source("R/hdp_XT2mat.R")

# Computation helpers
source("R/do_outer_mat.R")
source("R/gibbs_tabs.R")
source("R/update_q_probs.R")

# Correlation methods
source("R/hdp_corr_anlys.R")
source("R/hdp_corr_smpl.R")