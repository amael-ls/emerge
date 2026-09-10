#### Aim of script: Run the models for pooled species
## Comments
# This file is to run the pooled species models

#### Load packages
library(data.table)
library(cmdstanr)
library(stringi)

#### Load data
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Tree data (all species)
tree_dt = readRDS(paste0(path_data, "tree_dt.rds"))
setkey(tree_dt, group)

## Seeds that were used to run the models (full and submodel)
seed_dt = readRDS(paste0(path_data, "ls_seeds.rds"))

## Load stan models
fullmodel = cmdstan_model(paste0(path_models, "fullmodel.stan"))
submodel = cmdstan_model(paste0(path_models, "submodel.stan"))



# --------------------------------------------------------------------------------------
# --------------------    Run the full model, for pooled species    --------------------
# --------------------------------------------------------------------------------------

#### Data table to record problems
success_dt = tree_dt[, .(N_indiv = .N), by = group]
success_dt[, rhat_ok := NA]
success_dt[, no_div := NA]

setkey(success_dt, group)

ls_groups = success_dt[, group]

#### Fit model
for (gr in ls_groups)
{
	print(paste("Running", gr))

	## Subset to targeted group for stanData
	stanData = list(
		N = tree_dt[.(gr)][, .N],
		bole_volume_m3 = tree_dt[.(gr)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(gr)][, total_volume_m3]
	)

	## Run model
	if (file.exists(paste0(path_output, gr, "_fullmodel_theta.rds")))
	{
		fit = readRDS(paste0(path_output, gr, "_fullmodel_theta.rds"))
	} else {
		fit = fullmodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		## Save results
		fit$save_output_files(dir = path_output, basename = gr, random = FALSE)
		saveRDS(fit, paste0(path_output, gr, "_fullmodel_theta.rds"))
	}

	div = plot_gr(fit, gr, simplif = FALSE, n_bins = 4, pal = "Hiroshige", selected_variable = "height",
		print_plot = FALSE)

	rhat_ok = max(div[["rhats"]]) <= 1.01
	no_div = !div[["any_div"]]

	success_dt[.(gr), rhat_ok := ..rhat_ok]
	success_dt[.(gr), no_div := ..no_div]

	if (div[["any_div"]])
		warning(paste("Divergences for", gr))
}

success_dt[, success_full := rhat_ok & no_div]



# ---------------------------------------------------------------------------------------
# ---------------------    Run the sub model, for pooled species    ---------------------
# ---------------------------------------------------------------------------------------

#### Completing data table recording problems
success_dt[, rhat_ok_sub := NA]
success_dt[, no_div_sub := NA]

#### Fit model
for (gr in ls_groups)
{
	print(paste("Running", gr))

	## Subset to targeted group for stanData
	stanData = list(
		N = tree_dt[.(gr)][, .N],
		bole_volume_m3 = tree_dt[.(gr)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(gr)][, total_volume_m3]
	)

	## Run model
	if (file.exists(paste0(path_output, gr, "_submodel.rds")))
	{
		fit = readRDS(paste0(path_output, gr, "_submodel.rds"))
	} else {
		fit = submodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		## Save results
		fit$save_output_files(dir = path_output, basename = paste0(gr, "_submodel"), random = FALSE)
		saveRDS(fit, paste0(path_output, gr, "_submodel.rds"))
	}

	## Check-up
	div = plot_gr(fit, gr, simplif = TRUE, n_bins = 4, pal = "Hiroshige", selected_variable = "height",
		print_plot = FALSE)

	rhat_ok = max(div[["rhats"]]) <= 1.01
	no_div = !div[["any_div"]]

	success_dt[.(gr), rhat_ok_sub := ..rhat_ok]
	success_dt[.(gr), no_div_sub := ..no_div]

	if (div[["any_div"]])
		warning(paste("Divergences for", gr))
}

success_dt[, success_sub := rhat_ok_sub & no_div_sub]

success_dt[, any_success := success_full | success_sub]

if (!file.exists(paste0(path_output, "group-success.rds")))
	saveRDS(success_dt, paste0(path_output, "group-success.rds"))



# ----------------------------------------------------------------------------------------
# ------------------    Run the full model, for broadleaves/conifers    ------------------
# ----------------------------------------------------------------------------------------

#### Fit model for broadleaves
## Common variables
filename = "broadleaf"
setkey(tree_dt, fct_type) # Reorganise tree_dt for fast subset

## Subset to broadleaves
stanData = list(
	N = tree_dt[.("broadleaf")][, .N],
	bole_volume_m3 = tree_dt[.("broadleaf")][, bole_volume_m3],
	total_volume_m3 = tree_dt[.("broadleaf")][, total_volume_m3]
)

## Run model
if (!file.exists(paste0(path_output, filename, ".rds")))
{
	fit = fullmodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
		iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10, seed = seed_dt["broadleaf", full])

	## Save results
	fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
	saveRDS(fit, paste0(path_output, filename, ".rds"))
}



#### Fit model for conifers
filename = "conifer"

## Subset to conifers
stanData = list(
	N = tree_dt[.("conifer")][, .N],
	bole_volume_m3 = tree_dt[.("conifer")][, bole_volume_m3],
	total_volume_m3 = tree_dt[.("conifer")][, total_volume_m3]
)

## Run model
if (!file.exists(paste0(path_output, filename, ".rds")))
{
	fit = fullmodel$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
		iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10, seed = seed_dt["broadleaf", full])

	## Save results
	fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
	saveRDS(fit, paste0(path_output, filename, ".rds"))
}
