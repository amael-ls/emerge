#### Aim of script: Prediction on coppices, coppice-with-standards and "outliers"
## Comment:
# This script is used for the appendix "Statistical breakdown per forest structure"

#### Load packages
library(data.table)
library(cmdstanr)
library(stringi)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Loading tree data
tree_dt = readRDS(paste0(path_data, "tree_dt.rds"))
tree_dt_14sp = readRDS(paste0(path_data, "tree_dt_14species.rds"))
ls_14sp = tree_dt_14sp[, unique(speciesName_sci)]

R2D2 = readRDS(paste0(path_output, "rsquared.rds"))[, .(species, selected)]

emerge = readRDS(paste0(path_data, "emerge_2009-2010.rds")) |>
	merge.data.table(unique(tree_dt[, .(speciesName_sci, group)]), by = "speciesName_sci", all.x = TRUE)

emerge[, any(is.na(group))]

ls_species = emerge[, unique(speciesName_sci)]
ls_groups = emerge[, sort(unique(group))]

if ("Pinus uncinata" %in% ls_species)
	stop("I assumed that Pinus uncinata is not present! The code is not adapted to this particular case")

success_dt = readRDS(paste0(path_output, "group-success.rds"))[ls_groups]

woodstock_seed = 1969 - 08 - 18
n_chains = 4

v_gen_list = vector(mode = "list", length = length(ls_species))
names(v_gen_list) = ls_species

v_res_list = vector(mode = "list", length = length(ls_species))
names(v_res_list) = ls_species

dt_sim = emerge[, .(speciesName_sci, dataset, tree_id, str_name)]
sp_model = readRDS(paste0(path_output, "species-model.rds"))

if (!file.exists(paste0(path_output, "v_gen_list.rds")) && !file.exists(paste0(path_output, "v_res_list.rds")))
{
	for (sp in ls_species)
	{
		filename = sp_model[.(sp), model]
		full = stri_detect(str = filename, regex = "fullmodel_theta") ||
			stri_detect(str = filename, regex = "broadleaf") ||
			stri_detect(str = filename, regex = "conifer")

		sp_specific = sp %in% ls_14sp
		fit = readRDS(filename)

		if (sp_specific)
		{
			stanData = list(
				N = tree_dt_14sp[sp, .N],
				bole_volume_m3 = tree_dt_14sp[sp, bole_volume_m3],
				total_volume_m3 = tree_dt_14sp[sp, total_volume_m3],
				N_new = emerge[sp, .N],
				bole_volume_m3_new = emerge[sp, bole_volume_conic_m3],
				total_volume_m3_new = emerge[sp, total_volume_m3]
			)
		} else {
			stanData = list(
				N = tree_dt[sp, .N],
				bole_volume_m3 = tree_dt[sp, bole_volume_m3],
				total_volume_m3 = tree_dt[sp, total_volume_m3],
				N_new = emerge[sp, .N],
				bole_volume_m3_new = emerge[sp, bole_volume_conic_m3],
				total_volume_m3_new = emerge[sp, total_volume_m3]
			)
		}

		if (full)
		{
			genQ = cmdstan_model(paste0(path_models, "fullmodel-genQ.stan"))
		} else {
			genQ = cmdstan_model(paste0(path_models, "submodel-genQ.stan"))
		}

		str_name = emerge[sp, str_name]

		sim = genQ$generate_quantities(fit, data = stanData,
			seed = woodstock_seed, parallel_chains = min(n_chains, 4))

		v_gen = t(posterior::as_draws_matrix(sim$draws("v_gen")))
		opposite_res = v_gen - emerge[sp, total_volume_m3]

		v_gen_list[[sp]] = cbind(dt_sim[sp], as.data.table(v_gen))
		v_res_list[[sp]] = cbind(dt_sim[sp], as.data.table(opposite_res))
	}
	saveRDS(v_gen_list, paste0(path_output, "v_gen_list.rds"))
	saveRDS(v_res_list, paste0(path_output, "v_res_list.rds"))
} else {
	v_gen_list = readRDS(paste0(path_output, "v_gen_list.rds"))
	v_res_list = readRDS(paste0(path_output, "v_res_list.rds"))
}
