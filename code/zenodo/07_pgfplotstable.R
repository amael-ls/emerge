#### Aim of script: Generate the files for pgfplots table Latex
## Comments:
# This file is not necessary for the study itself. It only creates csv files
#	to be read by latex using pgfplotstable

#### Load packages
library(data.table)

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
	comp = rebuilt$comp
	weights_dt = rebuilt$weights_dt
	R2D2 = rebuilt$R2D2
	rhat_dt = rebuilt$rhat

	comp = merge.data.table(comp, rhat_dt, by = "species")
	comp = merge.data.table(comp, tree_dt[, .N, by = speciesName_sci],
		by.x = "species", by.y = "speciesName_sci")

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
	fwrite(R2D2, paste0(path_pgfplotstable, "rsquared.csv"), na = "NaN")

if (file.exists(paste0(path_output, "longuetaud_VEF.rds")))
{
	longuetaud_pars = readRDS(paste0(path_output, "longuetaud_VEF.rds"))
	if (!file.exists(paste0(path_pgfplotstable, "longuetaud_VEF.csv")))
		fwrite(longuetaud_pars, paste0(path_pgfplotstable, "longuetaud_VEF.csv"))
} else {
	stop("You must run 06_VEF-longuetaud.R before")
}
