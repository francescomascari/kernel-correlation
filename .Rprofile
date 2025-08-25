#source("renv/activate.R")

#require("jsonlite")

# package for data manipolation
require("tidyverse")

# functions for computation of quadratic forms
require("emulator", include.only = c("quad.form", "quad.3form"))

# function for normalization of probabilities
# in log space
require("matrixStats", include.only = "logSumExp")

# package for LaTeX-style caption in plots
require("extrafont")
# REMINDER: make sure Latin Modern (LM) Roman 10
#           is installed in your system (source: https://www.ctan.org/tex-archive/fonts/lm)
# REMINDER: use `font_import()` the first time and after every R installation
loadfonts() # load fonts

## HELPER FUNCTIONS FOR SAMPLING
## FROM HIERARCHICAL DIRICHLET PROCESS
source("R/hdp_mat_sampler_help.R")
source("R/hdp_mat_sampler.R")
source("R/hdp_XT_sampler_help.R")
source("R/hdp_XT_sampler.R")
source("R/hdp_XT2mat.R")


## HELPER FUNCTIONS FOR COMPUTATIONS
source("R/do_outer_mat.R")
source("R/gibbs_tabs.R")
source("R/update_q_probs.R")

## METHODS TO COMPUTE CORRELATION
source("R/hdp_corr_anlys.R")
source("R/hdp_corr_smpl.R")