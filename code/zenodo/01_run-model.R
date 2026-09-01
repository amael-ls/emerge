#### Aim of script: Run the models for the 14 species
## Comments
# This file is to run the species-specific models only. For the group models, see 02_run-model_groups.R
# On a not too old computer (say from 2020), this script should take less than an hour to run. The longest
# species to paramtrise are, without surprise, the most abundant.
# For a quick test, use sp = "Fraxinus excelsior", it is the fastest to run

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

# Tree data (14 species)
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
ls_species = tree_dt[, unique(speciesName_sci)]

# Seeds that were used to run the models (full and submodel)
seed_dt = readRDS(paste0(path_data, "ls_seeds.rds"))

# Load stan models
fullmodel = cmdstan_model(paste0(path_models, "fullmodel.stan"))
submodel = cmdstan_model(paste0(path_models, "submodel.stan"))



# --------------------------------------------------------------------------------------
# --------------------    Run the full model, for the 14 species    --------------------
# --------------------------------------------------------------------------------------

for (sp in ls_species)
{
	print(paste("Running", sp))

	filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_fullmodel")

	if (sp == "Quercus sp.")
	{
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")
		sp = ls_species[stri_detect(str = ls_species, regex = "Quercus")]
		warning("For quercus sp., I also add all the other Quercus")
	}
	
	# Subset to targeted species for stanData
	stanData = list(
		N = tree_dt[.(sp)][, .N],
		bole_volume_m3 = tree_dt[.(sp)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(sp)][, total_volume_m3]
	)

	# Run full model
	if (file.exists(paste0(path_output, filename, ".rds")))
	{
		fit = readRDS(paste0(path_output, filename, ".rds"))
	} else {
		fit = fullmodel$sample(data = stanData, chains = n_chains,
			parallel_chains = min(n_chains, 4), seed = seed_dt[.(sp), full],
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		# Save results
		fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
		saveRDS(fit, paste0(path_output, filename, ".rds"))
	}

	div = plot_sp(fit, sp, simplif = FALSE, n_bins = 4, pal = "Hiroshige",
		selected_variable = "height", print_plot = TRUE)

	if (div[["any_div"]])
		warning(paste("Divergences for", sp))
}



# ---------------------------------------------------------------------------------------
# ---------------------    Run the sub model, for the 14 species    ---------------------
# ---------------------------------------------------------------------------------------

for (sp in ls_species)
{
	print(paste("Running", sp))

	filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_submodel")

	if (sp == "Quercus sp.")
	{
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")
		sp = ls_species[stri_detect(str = ls_species, regex = "Quercus")]
		warning("For quercus sp., I also add all the other Quercus")
	}
	
	# Subset to targeted species for stanData
	stanData = list(
		N = tree_dt[.(sp)][, .N],
		bole_volume_m3 = tree_dt[.(sp)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(sp)][, total_volume_m3]
	)

	# Run sub model
	if (file.exists(paste0(path_output, filename, ".rds")))
	{
		fit = readRDS(paste0(path_output, filename, ".rds"))
	} else {
		fit = submodel$sample(data = stanData, chains = n_chains,
			parallel_chains = min(n_chains, 4), seed = seed_dt[.(sp), submodel],
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		# Save results
		fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
		saveRDS(fit, paste0(path_output, filename, ".rds"))
	}

	source("./tool_functions.R")
	div = plot_sp(fit, sp, simplif = TRUE, n_bins = 4, pal = "Hiroshige",
		selected_variable = "height", print_plot = TRUE)

	if (div[["any_div"]])
		warning(paste("Divergences for", sp))
}



# ----------------------------------------------------------------------------------------
# ------------------    Run the full model, for broadleaves/conifers    ------------------
# ----------------------------------------------------------------------------------------

## Run model for broadleaves
filename = "broadleaf"

# Subset to broadleaves
broadleaves = c("Fagus sylvatica", "Fraxinus excelsior", "Quercus petraea", "Quercus sp.")
all(broadleaves %in% tree_dt[, unique(speciesName_sci)])

stanData = list(
	N = tree_dt[.(broadleaves)][, .N],
	bole_volume_m3 = tree_dt[.(broadleaves)][, bole_volume_m3],
	total_volume_m3 = tree_dt[.(broadleaves)][, total_volume_m3]
)

# Run model
if (!file.exists(paste0(path_output, filename, ".rds")))
{
	fit = fullmodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
		iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10, seed = seed_dt["broadleaf", full])

	## Save results
	fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
	saveRDS(fit, paste0(path_output, filename, ".rds"))
}



## Run model for coniefers
filename = "conifer"

# Subset to conifers
conifers = tree_dt[, unique(speciesName_sci)]
conifers = conifers[!(conifers %in% broadleaves)]

length(conifers) + length(broadleaves) == length(tree_dt[, unique(speciesName_sci)])

stanData = list(
	N = tree_dt[.(conifers)][, .N],
	bole_volume_m3 = tree_dt[.(conifers)][, bole_volume_m3],
	total_volume_m3 = tree_dt[.(conifers)][, total_volume_m3]
)

# Run model
if (!file.exists(paste0(path_output, filename, ".rds")))
{
	fit = fullmodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
		iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10, seed = seed_dt["broadleaf", full])

	## Save results
	fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
	saveRDS(fit, paste0(path_output, filename, ".rds"))
}
