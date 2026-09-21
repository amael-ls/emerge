#### Aim of script: Generate the files for pgfplots table Latex
## Comments:
# This file is not necessary for the study itself. It only creates csv files
#	to be read by latex using pgfplotstable

#### Load packages
library(data.table)
library(stringi)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Load data
tree_dt = readRDS(paste0(path_data, "tree_dt.rds"))
tree_dt_14 = readRDS(paste0(path_data, "tree_dt_14species.rds"))

ls_species = tree_dt_14[, unique(speciesName_sci)]

#### Create table comparison models 14 species
## Load data full vs sub model
comp = data.table(species = ls_species, best = "", elpd_diff = -Inf, se_diff = -Inf,
	warning = FALSE, key = "species")

weights_dt = data.table(species = ls_species, best = "", W_full = -Inf, W_sub = -Inf, key = "species")

rhat_dt = data.table(species = ls_species, Rhat_full = -Inf, Rhat_sub = -Inf, key = "species")

R2D2 = data.table(species = ls_species,
	R2_full = -Inf, R2_sub = -Inf, R2_Vtot_full = -Inf, R2_Vtot_sub = -Inf,
	R2_loo_full = -Inf, R2_loo_sub = -Inf, R2_loo_Vtot_full = -Inf, R2_loo_Vtot_sub = -Inf,
	key = "species")

if (file.exists(paste0(path_output, "comparison_full-sub.rds")))
{
	save_ls = readRDS(paste0(path_output, "comparison_full-sub.rds"))
	rebuilt = rebuild_comp(save_ls)
	weights_dt = rebuilt$weights_dt
	R2D2 = rebuilt$R2D2
	rhat_dt = rebuilt$rhat

	comp = readRDS(paste0(path_output, "comparison_dt.rds"))

	rmse_mape_summary = readRDS(paste0(path_output, "rmse.rds"))
} else {
	stop("You need to run 03_compare-models.R")
}

if (!file.exists(paste0(path_pgfplotstable, "comparison.csv")))
	fwrite(comp, paste0(path_pgfplotstable, "comparison.csv"), na = "NaN")

#### Work on R2
## Load R squared for bole volume
r2_bole = readRDS(paste0(path_data, "rsquared-bole_volume.rds"))

## Merge both R2
R2D2 = merge.data.table(R2D2, r2_bole, by.x = "species", by.y = "speciesName_sci", all.x = TRUE)
R2D2 = merge.data.table(R2D2, comp[, .(species, best)], by = "species")

## Modify manually the column "selected" for Pinus laricio and strobus, as non-significant diff.
R2D2[, selected := best]
if (any(comp[c("Pinus laricio", "Pinus strobus"), elpd_diff > 4]))
	stop("The ELPD was below three for Pinus laricio and strobus but that does not seem the case anymore")

R2D2[c("Pinus laricio", "Pinus strobus"), selected := "full"]

## Add column R2 selected model and remove best
R2D2[, R2_loo_selected := ifelse(selected == "full", R2_loo_full, R2_loo_sub)]
R2D2[, R2_loo_selected_Vtot := ifelse(selected == "full", R2_loo_Vtot_full, R2_loo_Vtot_sub)]

R2D2[, R2_selected := ifelse(selected == "full", R2_full, R2_sub)]
R2D2[, R2_selected_Vtot := ifelse(selected == "full", R2_Vtot_full, R2_Vtot_sub)]

R2D2[, best := NULL]

## Add column RMSE
R2D2 = merge.data.table(R2D2, rmse_mape_summary[, .(species, rmse_med = rmse_q50)], by = "species")

mean_totvol = tree_dt_14[, .(meanV = mean(total_volume_m3), q25 = quantile(total_volume_m3, 0.25),
	q50 = quantile(total_volume_m3, 0.50), q75 = quantile(total_volume_m3, 0.75),
	m = min(total_volume_m3), M = max(total_volume_m3)), by = speciesName_sci]

R2D2 = merge.data.table(R2D2, mean_totvol, by.x = "species", by.y = "speciesName_sci")
R2D2[, rmse_percent := rmse_med/meanV*100]

if (!file.exists(paste0(path_pgfplotstable, "rsquared.csv")))
{
	saveRDS(R2D2, paste0(path_output, "rsquared.rds"))
	fwrite(R2D2, paste0(path_pgfplotstable, "rsquared.csv"), na = "NaN")
}

if (file.exists(paste0(path_output, "longuetaud_VEF.rds")))
{
	longuetaud_pars = readRDS(paste0(path_output, "longuetaud_VEF.rds"))
	if (!file.exists(paste0(path_pgfplotstable, "longuetaud_VEF.csv")))
		fwrite(longuetaud_pars, paste0(path_pgfplotstable, "longuetaud_VEF.csv"))
} else {
	stop("You must run 06_VEF-longuetaud.R before")
}



# -----------------------------------------------------------------------------
# ------------    Extract parameters for all species and groups    ------------
# -----------------------------------------------------------------------------

ls_species = tree_dt[, unique(speciesName_sci)]

success_dt = readRDS(paste0(path_output, "group-success.rds"))
sp_specific_models = readRDS(paste0(path_output, "rsquared.rds"))[, .(species, model = selected)]

group_models = readRDS(paste0(path_output, "comparison_dt_group.rds"))

# Modify manually for groups with non-significant differences
group_models[(elpd_diff < 4) & (best == "sub"), best := full]

# Add best model information to
success_dt = merge.data.table(x = success_dt, y = group_models[, .(group, best)], by = "group")

#### Assign to each species a model (either sp-specific or group or generic model)
## Species-specific models
tree_dt = merge.data.table(tree_dt, sp_specific_models, by.x = "speciesName_sci",
	by.y = "species", all.x = TRUE)

# Link to files
tree_dt[model == "full",
	model := paste0(path_output, stri_replace(speciesName_sci, regex = " ", replacement = "-"), "_fullmodel_theta.rds")]

tree_dt[!is.na(model),
	model := stri_replace(model, regex = "._fullmodel", replacement = "_fullmodel")]

tree_dt[model == "sub",
	model := paste0(path_output, stri_replace(speciesName_sci, regex = " ", replacement = "-"), "_submodel.rds")]

tree_dt[("Pinus uncinata"),
	model := paste0(path_output, stri_replace(speciesName_sci, regex = " ", replacement = "-"), "_logit.rds")]

## Group models (I assume that only the full model was selected, the case in my study)
if (success_dt[, any(best != "full")])
	stop("I assumed that only the full model was selected")

ls_full_gr = success_dt[(success_full), group]
tree_dt[is.na(model) & group %in% ls_full_gr, model := paste0(path_output, group, "_fullmodel_theta.rds")]

## Assign generic models for groups that were not parametrised (A2), and for groups with too little indiv
ls_pb = c(tree_dt[(is.na(model)), unique(group)],
	success_dt[N_indiv < 150, group],
	success_dt[!(any_success), group]) |> unique()
tree_dt[(group %in% ls_pb) & (fct_type == "broadleaf"), model := paste0(path_output, "broadleaf.rds")]
tree_dt[(group %in% ls_pb) & (fct_type == "conifer"), model := paste0(path_output, "conifer.rds")]

## Check that all the species have a group
tree_dt[, any(is.na(model))]

if (!file.exists(paste0(path_output, "species-model.rds")))
	saveRDS(unique(tree_dt[, .(speciesName_sci, group, model)]), paste0(path_output, "species-model.rds"))

#### Compute the averaged parameters for all species
if (!file.exists(paste0(path_output, "avg_params_full.rds")))
{
	params_dt_full = unique(tree_dt[!("Pinus uncinata"), .(speciesName_sci, group)]) # Remove P. uncinata
	pars_names = c("c", "j", "k", "m", "n", "s", "tau")
	params_dt_full[, c(pars_names) :=
		.(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)]

	## Species-specific params
	for (sp in sp_specific_models[, unique(species)])
	{
		if (sp == "Pinus uncinata")
			next;

		print(paste("Species:", sp))
		sp_filename = stri_replace(str = sp, replacement = "-", regex = " ") |>
			stri_replace(replacement = "", regex = "\\.")

		load_sub = FALSE
		if (sp_specific_models[sp, model] == "full")
			filename = paste0(path_output, sp_filename, "_fullmodel_theta", ".rds")

		if (sp_specific_models[sp, model] == "sub")
		{
			filename = paste0(path_output, sp_filename, "_fullmodel_theta", ".rds")
			load_sub = TRUE
		}

		fit = readRDS(filename)

		if (load_sub)
		{
			paramsVec_simplif = getParams(model_cmdstan = fit,
				params_names = c("alpha", "beta_", "gamma", "delta"), type = "mean")
			paramsVec = c(
				c = unname(paramsVec_simplif["alpha"]),
				j = 1,
				k = unname(paramsVec_simplif["beta_"]),
				m = unname(exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] +
					paramsVec_simplif["alpha"]),
				n = unname(paramsVec_simplif["delta"] + paramsVec_simplif["alpha"]),
				s = unname(paramsVec_simplif["beta_"])
			)
		} else {
			paramsVec = getParams(fit, params_names = pars_names, type = "mean")
		}

		rm(fit)
		params_dt_full[.(sp), c(pars_names) := as.list(paramsVec)]
	}

	## Group params for species WITHOUT species-specific models
	for (gp in params_dt_full[, unique(group)])
	{
		print(paste("Group:", gp))
		ls_species = params_dt_full[group == gp & is.na(c), unique(speciesName_sci)] # Non-parametrised species
		filename = tree_dt[ls_species, unique(model)]
		if (length(filename) != 1)
			stop("There should be a unique model to load!")

		fit = readRDS(filename)

		load_sub = TRUE
		generic = (stri_detect(filename, regex = "broadleaf") || stri_detect(filename, regex = "conifer"))
		if ((success_dt[gp, success_full]) || (generic))
		{
			paramsVec = getParams(fit, params_names = pars_names, type = "mean")
			load_sub = FALSE
		}

		if (load_sub)
		{
			warning("I thought there were only full or generic models... Check that out!")
			paramsVec_simplif = getParams(model_cmdstan = fit,
				params_names = c("alpha", "beta_", "gamma", "delta"), type = "mean")
			paramsVec = c(
				c = unname(paramsVec_simplif["alpha"]),
				j = 1,
				k = unname(paramsVec_simplif["beta_"]),
				m = unname(exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] +
					paramsVec_simplif["alpha"]),
				n = unname(paramsVec_simplif["delta"] + paramsVec_simplif["alpha"]),
				s = unname(paramsVec_simplif["beta_"])
			)
		}

		rm(fit)
		params_dt_full[ls_species, c(pars_names) := as.list(paramsVec)]
	}

	## Generic models params
	# Broadleaves
	fit = readRDS(paste0(path_output, "broadleaf.rds"))
	broadleaf = getParams(fit, params_names = pars_names, type = "mean")
	rm(fit)

	# Conifers
	fit = readRDS(paste0(path_output, "conifer.rds"))
	conifer = getParams(fit, params_names = pars_names, type = "mean")
	rm(fit)

	temp_dt = rbindlist(l = list(broadleaf = as.list(broadleaf), conifer = as.list(conifer)),
		idcol = "speciesName_sci")
	temp_dt[, group := NA_character_]
	setcolorder(temp_dt, neworder = names(params_dt_full))

	params_dt_full = rbindlist(l = list(params_dt_full, temp_dt))

	saveRDS(params_dt_full, paste0(path_output, "avg_params_full.rds"))
}
