#### Aim of script: Compare our approche with refitted Vallet 2006 in Bayesian
## Comments:
# We compare our approach fitted in 01 and 02, and selected in 03 with the Bayesian
# 	version of Vallet 2006 fitted in file 04 for 7 species

## Packages needed to reproduce the study
renv::restore()

library(data.table)
library(cmdstanr)
library(stringi)

## Load data
# Tool functions
source("./tool_functions.R")

# Global variables (paths and others)
source("./global_variables.R")
