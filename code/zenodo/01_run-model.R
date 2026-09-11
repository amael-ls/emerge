#### Aim of script: Run the models for the 14 species
## Comments
# This file is to run the species-specific models only. For the group models, see 02_run-model_groups.R
# On a not too old computer (say from 2020), this script should take less than an hour to run. The longest
# species to paramtrise are, without surprise, the most abundant.
# For a quick test, use sp = "Fraxinus excelsior", it is the fastest to run

#### Load packages
library(data.table)
library(cmdstanr)
library(stringi)

options(max.print = 500)

#### Load data
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Tree data (14 species)
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
ls_species = tree_dt[, unique(speciesName_sci)]

## Seeds that were used to run the models (full and submodel)
seed_dt = readRDS(paste0(path_data, "ls_seeds.rds"))

## Load stan models
fullmodel = cmdstan_model(paste0(path_models, "fullmodel.stan"))
submodel =  cmdstan_model(paste0(path_models, "submodel.stan"))



# --------------------------------------------------------------------------------------
# --------------------    Run the full model, for the 14 species    --------------------
# --------------------------------------------------------------------------------------

for (sp in ls_species)
{
	print(paste("Running", sp))

	filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_fullmodel_theta")
	le_cid = seed_dt[.(sp), full]

	if (sp == "Quercus sp.")
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")

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
			parallel_chains = min(n_chains, 4), seed = le_cid,
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		# Save results
		fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
		saveRDS(fit, paste0(path_output, filename, ".rds"))
	}

	div = plot_sp(fit, sp, simplif = FALSE, n_bins = 4, pal = "Hiroshige",
		selected_variable = "height", print_plot = TRUE)
	plot_divergences(fit, div$loc)

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
	le_cid = seed_dt[.(sp), submodel]

	if (sp == "Quercus sp.")
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")

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
			parallel_chains = min(n_chains, 4), seed = le_cid,
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		# Save results
		fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
		saveRDS(fit, paste0(path_output, filename, ".rds"))
	}

	div = plot_sp(fit, sp, simplif = TRUE, n_bins = 4, pal = "Hiroshige",
		selected_variable = "height", print_plot = TRUE)
	plot_joint_dirty(fit, div$loc)

	if (div[["any_div"]])
		warning(paste("Divergences for", sp))
}



# ----------------------------------------------------------------------------------------
# ---------------------    Run a tailored model for Pinus uncinata   ---------------------
# ----------------------------------------------------------------------------------------

sp = "Pinus uncinata"
pinus_uncinata = cmdstan_model(paste0(path_models, "pinus_uncinata.stan"))

filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_logit")
le_cid = seed_dt[.(sp), submodel]

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
	fit = pinus_uncinata$sample(data = stanData, chains = n_chains,
		parallel_chains = min(n_chains, 4), seed = le_cid,
		iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

	# Save results
	fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
	saveRDS(fit, paste0(path_output, filename, ".rds"))
}

div = posterior::subset_draws(fit$sampler_diagnostics(), variable = "divergent__")
if (any(div != 0))
	stop("There are divergences for Pinus uncinata")

mu_logit_fct = function(x, pars)
	return (inv_logit(pars["logit_alpha"] + exp(-pars["beta_"]*x) *
		(pars["gamma"]*x + pars["delta"])))

vtot_logit_fct = function(x, pars)
	return (x/mu_logit_fct(x, pars))

pars = getParams(model_cmdstan = fit,
	params_names = c("logit_alpha", "beta_", "gamma", "delta"))

## Plot fit on observed data with ratio
plot(stanData$bole_volume_m3, stanData$bole_volume_m3/stanData$total_volume_m3, pch = 19,
	cex = 0.75, axes = FALSE, xlab = "Observed bole", ylab = "Ratio")
curve(mu_logit_fct(x, pars), add = TRUE, col = "#CD212A", lwd = 3)
axis(1)
axis(2, las = 1)

## Plot fit on observed data with total volume
plot(stanData$bole_volume_m3, stanData$total_volume_m3, pch = 19,
	cex = 0.75, axes = FALSE, xlab = "Observed bole volume", ylab = "Observed total volume")
curve(vtot_logit_fct(x, pars), add = TRUE, col = "#CD212A", lwd = 3)
abline(a = 0, b = 1, lty = "dashed", lwd = 0.75, col = "#9A9A9A")
axis(1)
axis(2, las = 1)

## Plot obs vs pred total volume
pred = fit$draws("v_gen_mean") |> apply(MARGIN = 3, FUN = mean)
plot(pred, stanData$total_volume_m3, pch = 19,
	cex = 0.75, axes = FALSE, xlab = "Predicted total volume", ylab = "Observed total volume")
abline(a = 0, b = 1, lty = "dashed", lwd = 0.75, col = "#9A9A9A")
axis(1)
axis(2, las = 1)
