#### Aim of script: Run the models for pooled species
## Comments
# This file is to run the pooled species models

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

# --------------------------------------------------------------------------------------
# --------------------    Run the full model, for pooled species    --------------------
# --------------------------------------------------------------------------------------

ls_species = tree_dt[, unique(speciesName_sci)]
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
