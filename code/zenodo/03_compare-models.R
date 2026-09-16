#### Aim of script: Run the models for pooled species
## Comments
# This file is to run the pooled species models
#
# Bibliography:
#	Sivula, Tuomas, Måns Magnusson, Asael Alonzo Matamoros, and Aki Vehtari. 2025.
#		Uncertainty in Bayesian Leave-One-Out Cross-Validation Based Model Comparison.
#		Bayesian Analysis 1 (1). https://doi.org/10.1214/25-ba1569.
#
#	McLatchie, Yann, and Aki Vehtari. 2024.
#		Efficient Estimation and Correction of Selection-Induced Bias with Order Statistics.
#		Statistics and Computing 34 (4): 132. https://doi.org/10.1007/s11222-024-10442-4.
#
# Number of individuals per group:
# A1   532
# A3 12236
# B1  1874
# B2  5302
# C1    32
# C2    15
# C3  6465
# D1    19
# D2   116
# E1   174
# E2  2237
#

#### Load packages
library(data.table)
library(cmdstanr)
library(stringi)

#### Load data
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Tree data (14 species)
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
tree_dt_full = readRDS(paste0(path_data, "tree_dt.rds"))
setkey(tree_dt_full, group)

# ------------------------------------------------------------------------------------------
# --------------------    Compare full and sub models for 14 species    --------------------
# ------------------------------------------------------------------------------------------

#### Compare models (ELPD), compute R2, check rhat
## Load/create necessary data
ls_species = tree_dt[, unique(speciesName_sci)]
comp = data.table(species = ls_species, best = "", elpd_diff = -Inf, se_diff = -Inf,
	warning = FALSE, key = "species")

weights_dt = data.table(species = ls_species, best = "", W_full = -Inf, W_sub = -Inf, key = "species")

rhat_dt = data.table(species = ls_species, Rhat_full = -Inf, Rhat_sub = -Inf, key = "species")

R2D2 = data.table(species = ls_species,
	R2_full = -Inf, R2_sub = -Inf, R2_Vtot_full = -Inf, R2_Vtot_sub = -Inf,
	R2_loo_full = -Inf, R2_loo_sub = -Inf, R2_loo_Vtot_full = -Inf, R2_loo_Vtot_sub = -Inf,
	key = "species")

## Run comparison
if (file.exists(paste0(path_output, "comparison_full-sub.rds")))
{
	save_ls = readRDS(paste0(path_output, "comparison_full-sub.rds"))
	rebuilt = rebuild_comp(save_ls)
	comp = rebuilt$comp
	weights_dt = rebuilt$weights_dt
	R2D2 = rebuilt$R2D2
	rhat_dt = rebuilt$rhat

	comp = merge.data.table(comp, rhat_dt, by = "species")
	comp = merge.data.table(comp, tree_dt[, .N, by = speciesName_sci],
		by.x = "species", by.y = "speciesName_sci")

} else {
	save_ls = vector(mode = "list", length = length(ls_species))
	names(save_ls) = ls_species

	for (sp in ls_species)
	{
		temp = comparison_full_sub(sp = sp, tree_dt = tree_dt, path_models = path_models, path_output = path_output)
		save_ls[[sp]] = temp
		comp[.(sp), c("best", "elpd_diff", "se_diff", "warning") :=
			.(temp$comploo[1, "model"], temp$comploo[2, "elpd_diff"], temp$comploo[2, "se_diff"], temp$warning)]

		weights_dt[.(sp), c("best", "W_full", "W_sub") :=
			.(temp$best, temp$weights["full"], temp$weights["sub"])]

		R2D2[.(sp), c("R2_full", "R2_sub") := .(
			median(temp[["rsq_distrib"]][["full"]]),
			median(temp[["rsq_distrib"]][["sub"]])
		)]

		R2D2[.(sp), c("R2_Vtot_full", "R2_Vtot_sub") := .(
			median(temp[["rsq_vtot_distrib"]][["full"]][["rsq_vtot"]]),
			median(temp[["rsq_vtot_distrib"]][["sub"]][["rsq_vtot"]])
		)]

		R2D2[.(sp), c("R2_loo_full", "R2_loo_sub") := .(
			median(temp[["rsq_loo_distrib_r"]][["full"]]),
			median(temp[["rsq_loo_distrib_r"]][["sub"]])
		)]

		R2D2[.(sp), c("R2_loo_Vtot_full", "R2_loo_Vtot_sub") := .(
			median(temp[["rsq_loo_distrib_v"]][["full"]]),
			median(temp[["rsq_loo_distrib_v"]][["sub"]])
		)]

		rhat_dt[sp, c("Rhat_full", "Rhat_sub") := .(temp[["rhat_full"]], temp[["rhat_sub"]])]
	}

	comp = merge.data.table(comp, rhat_dt, by = "species")
	comp = merge.data.table(comp, tree_dt[, .N, by = speciesName_sci], by.x = "species", by.y = "speciesName_sci")
	saveRDS(save_ls, paste0(path_output, "comparison_full-sub.rds"))
}

# Compute proba model A better than B, based on Sivula et al. 2025
comp[, p_better := pnorm(0, mean = elpd_diff, sd = se_diff)]

# Based on the appendix A of McLatchie2024, Equation A2 and A4
# I use pnorm as I noticed the integrate function is unstable for large values of diff_ELPD for Eq A2
comp[!is.na(elpd_diff),
	eq_A2 := pnorm(0, mean = abs(elpd_diff), sd = se_diff, lower.tail = FALSE), by = species]

# Eq A4, pseudp Bayesian model average
f = function(x, diff_elpd, se_elpd)
	return(dnorm(x, mean = 0, sd = se_elpd)/(1 + exp(-diff_elpd - x)))

comp[!is.na(elpd_diff),
	pseudo_bma := integrate(f, lower = -Inf, upper = +Inf, diff_elpd = abs(elpd_diff), se_elpd = se_diff)$value,
	by = species] # Eq A4

if (!file.exists(paste0(path_output, "comparison_dt.rds")))
	saveRDS(comp, paste0(path_output, "comparison_dt.rds"))

## Compute RMSE
rmse_ls = vector(mode = "list", length = length(ls_species))
names(rmse_ls) = ls_species

if (!file.exists(paste0(path_output, "rmse.rds")))
{
	for (sp in ls_species)
	{
		print(sp)
		is_simplif = FALSE
		if (sp %in% c("Fraxinus excelsior", "Picea abies", "Pinus laricio"))
			is_simplif = TRUE

		rmse_ls[[sp]] = RMSE_bayes(sp, tree_dt, path_output, path_models, is_simplif)
	}

	rmse_mape_dt = rbindlist(rmse_ls, idcol = "species")

	rmse_mape_summary = rmse_mape_dt[, .(
		rmse_min = min(rmse), rmse_q025 = quantile(rmse, 0.025), rmse_q50 = median(rmse),
		rmse_q975 = quantile(rmse, 0.975), rmse_max = max(rmse),

		mape_min = min(mape), mape_q025 = quantile(mape, 0.025), mape_q50 = median(mape),
		mape_q975 = quantile(mape, 0.975), mape_max = max(mape)),
	by = species]
	saveRDS(rmse_mape_dt, paste0(path_output, "rmse_ls.rds"))
	saveRDS(rmse_mape_summary, paste0(path_output, "rmse.rds"))
} else {
	rmse_mape_summary = readRDS(paste0(path_output, "rmse.rds"))
}



# ------------------------------------------------------------------------------------------
# -------------------    Compare full and sub models for group models    -------------------
# ------------------------------------------------------------------------------------------

#### Compare models (ELPD), compute R2, check rhat
## Load/create necessary data
ls_groups = tree_dt_full[, unique(group)]
comp = data.table(group = ls_groups, best = "", elpd_diff = -Inf, se_diff = -Inf,
	warning = FALSE, key = "group")

weights_dt = data.table(group = ls_groups, best = "", W_full = -Inf, W_sub = -Inf, key = "group")

rhat_dt = data.table(group = ls_groups, Rhat_full = -Inf, Rhat_sub = -Inf, key = "group")

R2D2 = data.table(group = ls_groups,
	R2_full = -Inf, R2_sub = -Inf, R2_Vtot_full = -Inf, R2_Vtot_sub = -Inf,
	R2_loo_full = -Inf, R2_loo_sub = -Inf, R2_loo_Vtot_full = -Inf, R2_loo_Vtot_sub = -Inf,
	key = "group")

## Run comparison
if (file.exists(paste0(path_output, "comparison_full-sub_groups.rds")))
{

} else {
	save_ls = vector(mode = "list", length = length(ls_groups))
	names(save_ls) = ls_groups

	for (gp in ls_groups)
	{
		print(paste("Running", gp, "N_indiv =", tree_dt_full[.(gp), .N]))
		temp = comparison_full_sub(sp = gp, tree_dt = tree_dt_full, path_models = path_models, path_output = path_output)
		save_ls[[gp]] = temp
		comp[.(gp), c("best", "elpd_diff", "se_diff", "warning") :=
			.(temp$comploo[1, "model"], temp$comploo[2, "elpd_diff"], temp$comploo[2, "se_diff"], temp$warning)]

		weights_dt[.(gp), c("best", "W_full", "W_sub") :=
			.(temp$best, temp$weights["full"], temp$weights["sub"])]

		R2D2[.(gp), c("R2_full", "R2_sub") := .(
			median(temp[["rsq_distrib"]][["full"]]),
			median(temp[["rsq_distrib"]][["sub"]])
		)]

		R2D2[.(gp), c("R2_Vtot_full", "R2_Vtot_sub") := .(
			median(temp[["rsq_vtot_distrib"]][["full"]][["rsq_vtot"]]),
			median(temp[["rsq_vtot_distrib"]][["sub"]][["rsq_vtot"]])
		)]

		R2D2[.(gp), c("R2_loo_full", "R2_loo_sub") := .(
			median(temp[["rsq_loo_distrib_r"]][["full"]]),
			median(temp[["rsq_loo_distrib_r"]][["sub"]])
		)]

		R2D2[.(gp), c("R2_loo_Vtot_full", "R2_loo_Vtot_sub") := .(
			median(temp[["rsq_loo_distrib_v"]][["full"]]),
			median(temp[["rsq_loo_distrib_v"]][["sub"]])
		)]

		rhat_dt[gp, c("Rhat_full", "Rhat_sub") := .(temp[["rhat_full"]], temp[["rhat_sub"]])]
	}

	comp = comp |> merge.data.table(rhat_dt, by = "group") |>
		merge.data.table(tree_dt_full[, .N, by = group], by = "group")

	# Save the comparison files
	saveRDS(comp, paste0(path_output, "comparison_dt_group.rds"))
	saveRDS(R2D2, paste0(path_output, "rsquared_group.rds"))
	saveRDS(weights_dt, paste0(path_output, "weights_dt_group.rds"))

	saveRDS(save_ls, paste0(path_output, "comparison_full-sub_groups.rds"))
}
