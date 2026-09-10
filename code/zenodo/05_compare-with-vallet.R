#### Aim of script: Compare our approche with refitted Vallet 2006 in Bayesian
## Comments:
# We compare our approach fitted in 01 and 02, and selected in 03 with the Bayesian
# 	version of Vallet 2006 fitted in file 04 for 7 species

#### Load packages
library(data.table)
library(cmdstanr)
library(stringi)

#### Load data
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Loading training dataset
tree_dt = readRDS(paste0(path_output, "pred_vallet.rds"))
training_dt = readRDS(paste0(path_data, "tree_dt_14species.rds")) # Dataset used for the full and submodel

## Species parametrised in Vallet et al. 2006
ls_species = c("Abies alba", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris",
	"Pseudotsuga menziesii", "Quercus petraea")
training_dt = training_dt[ls_species]
setkey(training_dt, speciesName_sci)

#### Run comparison
## Compile models
vallet_genQ = cmdstan_model(paste0(path_models, "vallet-genQ.stan"))
full_genQ = cmdstan_model(paste0(path_models, "fullmodel-genQ.stan"))
sub_genQ = cmdstan_model(paste0(path_models, "submodel-genQ.stan"))

comp_dt = data.table(speciesName_sci = ls_species, best = "", elpd_diff = -Inf, se_diff = -Inf,
	wrn_mine = "", wrn_vallet = "", key = "speciesName_sci")

R2D2 = data.table(species = ls_species, R2_vallet = -Inf, key = "species")

rsq_distrib = vector(mode = "list", length = R2D2[, .N])
names(rsq_distrib) = ls_species

for (sp in ls_species)
{
	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")
	filename = paste0(path_output, sp_filename, "_vallet.rds")
	is_douglas = ifelse(sp == "Pseudotsuga menziesii", 1, 0) # Boolean style compatible with Stan language

	stanData_gen = list(
		N = tree_dt[sp, .N],
		N_params = vallet_dt[sp, n_params],
		is_douglas = is_douglas,

		bole_volume_m3 = tree_dt[sp, bole_volume_m3],
		circumference_cm = tree_dt[sp, 100*circumference_m], # Was in cm in Vallet2006!
		height = tree_dt[sp, height],

		total_volume_m3 = tree_dt[sp, total_volume_m3],

		# New data, which are the same...
		N_new = tree_dt[sp, .N],

		bole_volume_m3_new = tree_dt[sp, bole_volume_m3],
		circumference_cm_new = tree_dt[sp, 100*circumference_m], # Was in cm in Vallet2006!
		height_new = tree_dt[sp, height],

		total_volume_m3_new = tree_dt[sp, total_volume_m3]
	)

	# Generate data Bayesian vallet
	fit_vallet = readRDS(filename)
	sim_vallet = vallet_genQ$generate_quantities(fitted_params = fit_vallet, data = stanData_gen,
		parallel_chains = min(4, n_chains))

	# Generate data my model
	if (sp %in% c("Fraxinus excelsior", "Pinus uncinata"))
	{
		sp_filename = paste0(sp_filename, "_submodel.rds")
		fit = readRDS(filename)
		sim = sub_genQ$generate_quantities(fitted_params = fit, data = stanData_gen,
			parallel_chains = min(4, n_chains))
	} else {
		filename = paste0(path_output, sp_filename, "_fullmodel_theta.rds")
		fit = readRDS(filename)
		sim = full_genQ$generate_quantities(fitted_params = fit, data = stanData_gen,
			parallel_chains = min(4, n_chains))
	}

	if (!file.exists(paste0(path_output, "comparison_vallet-mine.rds")))
	{
		## Compute PSIS-LOO...
		# ... for Vallet
		r_eff = loo::relative_eff(exp(sim_vallet$draws("log_lik")), cores = 8)
		loo_vallet = loo::loo(x = sim_vallet$draws("log_lik"), r_eff = r_eff, cores = 8)

		warning_vallet = "none"
		if (any(loo_vallet$diagnostics$pareto_k >= 0.7))
		{
			warning_vallet = "bad"
			n_bad = length(loo_vallet$diagnostics$pareto_k[(loo_vallet$diagnostics$pareto_k >= 0.7) &
				(loo_vallet$diagnostics$pareto_k < 1)])
			if (any(loo_vallet$diagnostics$pareto_k > 1))
			{
				warning_vallet = "very bad"
				n_verybad = sum(loo_vallet$diagnostics$pareto_k > 1)
			}
		}

		# ... for my model/submodel
		r_eff = loo::relative_eff(exp(sim$draws("log_lik")), cores = 8)
		loo_mine = loo::loo(x = sim$draws("log_lik"), r_eff = r_eff, cores = 8)

		warning_mine = "none"
		if (any(loo_mine$diagnostics$pareto_k >= 0.7))
		{
			warning_mine = "bad"
			if (any(loo_mine$diagnostics$pareto_k > 1))
				warning_mine = "very bad"
		}

		comp = loo::loo_compare(list(vallet = loo_vallet, mine = loo_mine))

		## Record comparison output
		comp_dt[sp, c("best", "elpd_diff", "se_diff") :=
			.(comp[1, "model"], comp[2, "elpd_diff"], comp[2, "se_diff"])]
		comp_dt[sp, c("wrn_mine", "wrn_vallet") :=
			.(warning_mine, warning_vallet)]

		## Compute R squared for Vallet (already done in 13_compare_models.qmd for my sub/model)
		# R squared for the volume, based on Gelman 2019 as there are Pareto warnings (disqualify R2-loo)
		var_fit_vallet = apply(X = posterior::as_draws_matrix(sim_vallet$draws("v_gen_mean")),
			MARGIN = 1, FUN = var) # The var contains the correction 1/(n - 1) already!
		var_res_vallet = apply(X = posterior::as_draws_matrix(fit_vallet$draws("sigma")),
			MARGIN = 1, FUN = mean)

		rsq_distrib[[sp]] = var_fit_vallet / (var_fit_vallet + var_res_vallet)
		R2D2[sp, R2_vallet := median(rsq_distrib[[sp]])]
	}

	## Check pred vs obs
	pred_vallet = apply(X = sim_vallet$draws("v_gen_mean"), MARGIN = 3, FUN = mean)
	pred_mine = apply(X = sim$draws("v_gen_mean"), MARGIN = 3, FUN = mean)

	rm(sim, sim_vallet, fit, fit_vallet)

	ind = pred_vallet < stanData_gen$bole_volume_m3

	pgfplots_dt = data.table(obs = tree_dt[.(sp), total_volume_m3], vallet_2006 = tree_dt[.(sp), freq_vallet],
		vallet_2026 = pred_vallet, mine = pred_mine)
	pgfplots_dt[, pb_vallet := as.integer(ind)] # 0 = no problem, 1 = tot < bole
	pgf_file = paste0(path_pgfplotsfig, "vallet-mine_", sp_filename, ".csv")

	if (!file.exists(pgf_file))
		fwrite(pgfplots_dt, pgf_file)
}

if (!file.exists(paste0(path_output, "comparison_vallet-mine.rds")))
{
	saveRDS(comp_dt, paste0(path_output, "comparison_vallet-mine.rds"))
	saveRDS(R2D2, paste0(path_output, "rsquared_vallet.rds"))

	fwrite(comp_dt, paste0(path_pgfplotstable, "comp_dt_vallet.csv"))
	fwrite(R2D2, paste0(path_pgfplotstable, "rsquared_vallet.csv"))
}
