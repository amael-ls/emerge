#### Aim of script: Make directories, check data
## Comments
# For cmdstanr, you might need to run (see https://mc-stan.org/cmdstanr/articles/cmdstanr.html):
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

## Packages needed to reproduce the study
renv::restore()
# renv::init() # TO REMOVE
# renv::deactivate(clean = TRUE)

library(data.table)
library(cmdstanr)
library(stringi)


# # To remove after
# library(MetBrewer)
# library(bayesplot)
# library(bayesboot)
# library(posterior)
# library(terra)
# library(loo)

# renv::snapshot() # Select 2

## Check/make directories
source("./global_variables.R")

if (!dir.exists(path_data))
	stop(paste0("The folder <", path_data, "> is missing. Check path and working directory"))

if (!dir.exists(path_models))
	stop(paste0("The folder <", path_models, "> is missing. Check path and working directory"))

if (!dir.exists(path_output))
	dir.create(path_output)

if (!dir.exists(path_pgfplots))
	dir.create(path_pgfplots)

## Check data
tree_file = paste0(path_data, "tree_dt.rds")
if (!file.exists(tree_file))
{
	# Loading Inra data (Oudin's protocol)
	filename = paste0(path_data, "inra.rds")
	if (!file.exists(filename))
		stop(paste0("The file <", filename, "> does not exist! Check path and working directory"))
	inra = readRDS(filename)
	inra = inra[, .(speciesName_sci, tree_id = unique_id, plot_id, fct_type, year, circumference_m, height,
		taper_height = taper_height_flo, bole_volume_m3, total_volume_m3)]
	setkey(inra, speciesName_sci)

	# Loading Swiss data (EFM protocol)
	filename = paste0(path_data, "switzerland.rds")
	if (!file.exists(filename))
		stop(paste0("The file <", filename, "> does not exist! Check path and working directory"))
	swiss = readRDS(filename)
	swiss = unique(swiss[, .(speciesName_sci, tree_id, plot_id, fct_type, year, circumference_m, height,
		taper_height = hdec, bole_volume_m3, total_volume_m3)])
	setkey(swiss, speciesName_sci)

	# Modern Emerge data (Oudin's protocol, data collected in 2 campaigns 2009-2010)
	filename = paste0(path_data, "emerge_2009-2010.rds")
	if (!file.exists(filename))
		stop(paste0("The file <", filename, "> does not exist! Check path and working directory"))
	emerge = readRDS(filename)
	emerge[, year := as.integer(stringi::stri_replace_all(str = dataset, replacement = "", regex = "emerge_"))]
	emerge = emerge[, .(speciesName_sci, tree_id, fct_type, year, circumference_m, height, taper_height,
		bole_volume_m3 = bole_volume_conic_m3, total_volume_m3)]

	# Species groups, determined during EMERGE project
	filename = paste0(path_data, "ls-groups.rds")
	if (!file.exists(filename))
		stop(paste0("The file <", filename, "> does not exist! Check path and working directory"))
	ls_groups = readRDS(filename)

	# Bind everything
	tree_dt = rbindlist(list(inra = inra, swiss = swiss, emerge = emerge), idcol = "origin", fill = TRUE)
	tree_dt[, r := bole_volume_m3/total_volume_m3]
	tree_dt = merge.data.table(x = tree_dt, y = ls_groups, by = "speciesName_sci", all.x = TRUE)
	tree_dt[, any(is.na(group))]
	saveRDS(tree_dt, tree_file)
} else {
	tree_dt = readRDS(tree_file)
}

## Keep only species with more than 150 individual records
ls_species = tree_dt[, .(n_indiv = .N), by = speciesName_sci][order(-n_indiv)]

min_indiv = 150
ls_species = ls_species[n_indiv > min_indiv]

## Join Quercus sp. (609) with other Quercus but petraea, i.e., with ilex, pubescens, robur, and rubra
tree_dt[speciesName_sci %in% c("Quercus ilex", "Quercus pubescens", "Quercus robur", "Quercus rubra"),
	speciesName_sci := "Quercus sp."]

## Subset tree_dt and save the 14 species dataset
tree_dt = tree_dt[speciesName_sci %in% ls_species[, speciesName_sci]]
setkey(tree_dt, speciesName_sci, plot_id)

filename = paste0(path_data, "tree_dt_14species.rds")
if (!file.exists(filename))
	saveRDS(tree_dt, filename)

## Check stan installation, help can be found at https://mc-stan.org/cmdstanr/index.html
cmdstanr::cmdstan_path()
if (cmdstanr::cmdstan_path() == "")
	warning("You might have to install Stan language or to set the correct path")
# cmdstanr::check_cmdstan_toolchain()
# cmdstanr::install_cmdstan(cores = 2)
# cmdstanr::cmdstan_version() # I used the version 2.39.0

## Compile models
# Full model
filename = paste0(path_models, "fullmodel.stan")
if (!file.exists(filename))
	stop(paste0("The model <", filename, "> could not be found"))
fullmodel = cmdstan_model(filename)

# Submodel
filename = paste0(path_models, "submodel.stan")
if (!file.exists(filename))
	stop(paste0("The model <", filename, "> could not be found"))
submodel = cmdstan_model(filename)
