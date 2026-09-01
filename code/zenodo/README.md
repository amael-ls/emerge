---
output: html_document
bibliography: references.bib  
---

# Readme

## Introduction

This zenodo repo is to reproduce the results of Anonymous et al. Please contact the main author ANONYMOUS if you are interested in recreating the dataset from raw data, as this is not provided in this repo. You can also find the scripts that were used to create data of Anonymous et al on github ANONYMOUS LINK.

## Description of the scripts

In order to reproduce the results, you just need to run the scripts in numerical order, where:
- `00_pre-run.R` verifies directories, data, and join [the three datasets](#description-of-the-datasets) to create the data used to paramtrise the 14 species. It also check Stan language and compiles the models
- `01_run-model.R` run the species-specific models and two generic models (broadleaf and conifer, not used in the study, but could be useful one day)
- `02_run-model_groups.R` run the pooled species models

## Description of the datasets

- data/emerge_2009-2010.rds: contains data collected in 2009--2010 during the EMERGE project [@Deleuze2013]
- data/inra.rds: contains data collected with Oudin's protocol in France [@Vallet2006;@Oudin1930]
- data/switzerland.rds: contains data collected in Switzerland following the EFM protocol [@Didion2024]

## Packages version

The version of the packages are recorded in the file *renv.lock* and are managed with the R package [`renv`](https://rstudio.github.io/renv/index.html).

The only software version not recorded is Stan language itself. I used the version 2.39.0 for this project. For R, I used *R version 4.6.1 (2026-06-24) -- "Happy Hop"*