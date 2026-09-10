#### Aim of script: Generate the files for pgfplots Latex
## Comments:
# This file is not necessary for the study itself. It only creates csv files
#	to be read by latex using pgfplots

rm(list = ls())

library(data.table)
library(stringi)
library(mgcv)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Other functions
pred_fct_Vtot = function(sp, stanData_gen, path_models, is_simplif = FALSE, n_points = 5,
	woodstock_seed = 1969 - 08 - 18)
{
	## Generate quantities for full/submodel
	genQ = cmdstanr::cmdstan_model(paste0(path_models, "fullmodel-genQ.stan"))
	filename = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"),
		"_fullmodel_theta.rds")

	if (is_simplif)
	{
		genQ = cmdstanr::cmdstan_model(paste0(path_models, "submodel-genQ.stan"))
		filename = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"),
			"_submodel.rds")
	}

	fit = readRDS(filename)
	n_chains = fit$num_chains()
	n_iter = fit$metadata()$iter_sampling

	sim = genQ$generate_quantities(fit, data = stanData_gen,
		seed = woodstock_seed, parallel_chains = min(n_chains, 4))

	rm(fit)

	## Predictions on volume
	pred = apply(X = sim$draws("v_gen_mean"), MARGIN = 3, FUN = mean)
	pred_quant_mean = apply(X = sim$draws("v_gen_mean"), MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))
	pred_points = sim$draws("v_gen")
	pred_quant = data.table(Vbole = stanData_gen$bole_volume_m3_new,
		q05 = rep(-Inf, stanData_gen$N_new), q95 = rep(-Inf, stanData_gen$N_new))

	## Predictions on ratio
	pred_r = apply(X = sim$draws("r_gen_mean"), MARGIN = 3, FUN = mean)
	pred_r_quant_mean = apply(X = sim$draws("r_gen_mean"), MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))
	pred_r_points = sim$draws("r_gen")
	pred_r_quant = data.table(Vbole = stanData_gen$bole_volume_m3_new,
		q05 = rep(-Inf, stanData_gen$N_new), q95 = rep(-Inf, stanData_gen$N_new))

	rm(sim)

	dt_v = matrix(data = -Inf, nrow = stanData_gen$N_new, ncol = n_points)
	dt_r = matrix(data = -Inf, nrow = stanData_gen$N_new, ncol = n_points)

	set.seed(woodstock_seed)

	for (i in seq_len(stanData_gen$N_new))
	{
		rand_ind = sample(x = 1:(n_iter*n_chains), size = n_points, replace = FALSE)
		dt_v[i,] = pred_points[, , i][rand_ind]
		dt_r[i,] = pred_points[, , i][rand_ind]
		pred_quant[i, c("q05", "q95") := as.list(quantile(pred_points[, , i], prob = c(0.05, 0.95)))]
		pred_r_quant[i, c("q05", "q95") := as.list(quantile(pred_r_points[, , i], prob = c(0.05, 0.95)))]
	}

	return (list(pred = pred, quant_mean = pred_quant_mean, quant = pred_quant,
		pred_r = pred_r, quant_mean_r = pred_r_quant_mean, quant_r = pred_r_quant,
		randomSample_v = dt_v, randomSample_r = dt_r, n_points = n_points))
}

## Load data
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
ls_species = tree_dt[, unique(speciesName_sci)]

range_dt = readRDS(paste0(path_data, "./range_nfi-pred_Vbole.rds"))
range_dt["Quercus sp.", q995_flopp_u := q995_nfi] # Because there is nothing, and I will use q995_flopp_u after
range_dt["Quercus sp.", min_nfi_flopp := min_nfi]

# Species need an extension if 110% max_t is below q995_flopp
range_dt[, need_extension := 1.1*max_t < q995_flopp_u]

#### Predict total volumes and compute quantiles
## Number of new data for predict
n_indiv = 500

for (sp in ls_species[12:14])
{
	print(paste("Running", sp))
	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")

	if (sp == "Quercus sp.")
		sp_filename = stri_replace(str = sp_filename, replacement = "", regex = "\\.")

	if (file.exists(paste0(path_pgfplotsfig, sp_filename, ".csv")))
		next

	## Simulate total above ground volume on nfi or training range (depends which is largest)
	stanData = list(
		N = tree_dt[sp, .N],
		bole_volume_m3 = tree_dt[sp, bole_volume_m3],
		total_volume_m3 = tree_dt[sp, total_volume_m3],
		N_new = n_indiv,
		bole_volume_m3_new = seq(range_dt[sp, min(min_t, min_nfi_flopp)],
			range_dt[sp, max(max_t, q995_flopp_u)], length.out = n_indiv)
	)

	stanData$total_volume_m3_new = 1.3*stanData$bole_volume_m3_new # Useless here, but must be provided for genQ to run

	fwrite(data.table(obs_bole_volume = stanData$bole_volume_m3, obs_total_volume = stanData$total_volume_m3,
		obs_ratio = stanData$bole_volume_m3/stanData$total_volume_m3),
		paste0(path_pgfplotsfig, sp_filename, "_obs-data.csv"))

	simplif = FALSE
	if ((sp == "Fraxinus excelsior") || (sp == "Pinus uncinata"))
		simplif = TRUE # Remember that laricio and strobus are fullmodel (ELPD_diff < 4)

	temp = pred_fct_Vtot(sp = sp_filename, stanData_gen = stanData, path_models = path_models,
		is_simplif = simplif)

	dt = data.table(bole_volume_m3 = stanData$bole_volume_m3_new, pred_avg = temp$pred, pred_avg_r = temp$pred_r,
		q05_v_u = temp$quant[, q05], q95_v_u = temp$quant[, q95],
		q05_v = temp$quant_mean["5%", ], q95_v = temp$quant_mean["95%", ],
		q05_r_u = temp$quant_r[, q05], q95_r_u = temp$quant_r[, q95],
		q05_r = temp$quant_mean_r["5%", ], q95_r = temp$quant_mean_r["95%", ]
	)
	dt[, beyond_rg := bole_volume_m3 > range_dt[sp, max_t]]
	dt[, beyond_rg_int := as.integer(beyond_rg)]

	## Plot
	# Compute boundaries
	max_x = dt[, max(bole_volume_m3)]
	max_y = max(temp$quant[, q95])
	min_y_r = min(temp$quant_r[, q05])
	max_y_r = max(temp$quant_r[, q95])
	if (max_y_r > 1)
		warning("Uncertainty leads beyond 1! How come with a Beta distribution??")

	# Smooth out the uncertainties, only to rmove small wiggles in the quantiles...
	m05_v = mgcv::gam(q05_v_u ~ s(bole_volume_m3, k = 30), data = dt)
	m95_v = mgcv::gam(q95_v_u ~ s(bole_volume_m3, k = 30), data = dt)
	m05_r = mgcv::gam(q05_r_u ~ s(bole_volume_m3, k = 30), data = dt)
	m95_r = mgcv::gam(q95_r_u ~ s(bole_volume_m3, k = 30), data = dt)

	# ... of volumes
	dt[, q05_v_u_s := predict(m05_v)]
	if (dt[, any(q05_v_u_s < 0)])
	{
		warning("Negative values due to smoothing out. Replacing by original value")
		dt[q05_v_u_s < 0, q05_v_u_s := q05_v_u]
	}

	dt[, q95_v_u_s := predict(m95_v)]
	if (dt[, any(q95_v_u_s < 0)])
	{
		warning("Negative values due to smoothing out. Replacing by original value")
		dt[q95_v_u_s < 0, q95_v_u_s := q95_v_u]
	}

	# ... of ratios
	dt[, q05_r_u_s := predict(m05_r)]
	if (dt[, any(q05_r_u_s < 0)])
	{
		warning("Negative values due to smoothing out. Replacing by original value")
		dt[q05_r_u_s < 0, q05_r_u_s := q05_r_u]
	}

	dt[, q95_r_u_s := predict(m95_r)]
	if (dt[, any(q95_r_u_s < 0)])
	{
		warning("Negative values due to smoothing out. Replacing by original value")
		dt[q95_r_u_s < 0, q95_r_u_s := q95_r_u]
	}

	## Save data for pgfplots
	fwrite(dt, paste0(path_pgfplotsfig, sp_filename, ".csv"))
}



#### DRAFT ZONE
## This is for sp = "Pinus uncinata"
params = getParams(fit, c("alpha", "beta_", "gamma", "delta"))

r_func_simplif = function(x, pars)
	return(pars["alpha"] + exp(-pars["beta_"]*x) * (pars["gamma"]*x + pars["delta"]))


curve(r_func_simplif(x, params), to = 40)

params_full = setDT(posterior::as_draws_df(getParams(fit, c("alpha", "beta_", "gamma", "delta", "phi"), "all")))
setnames(params_full, old = c(".chain", ".iteration", ".draw"), new = c("chain", "iteration", "draw"))

shape2 = function(x, pars)
{
	r_func_simplif = function(x, pars)
		return(pars[, alpha] + exp(-pars[, beta_]*x) * (pars[, gamma]*x + pars[, delta]))
	
	return (pars[, phi]*(1 - r_func_simplif(x, pars)))
}


r_func_simplif = function(x, pars)
	return(pars[, alpha] + exp(-pars[, beta_]*x) * (pars[, gamma]*x + pars[, delta]))

curve(r_func_simplif(x, params_full[1, .(alpha, beta_, gamma, delta)]), to = 3, ylim = c(0, 1.2))
for (i in 1:4000)
{
	curve(r_func_simplif(x, params_full[i, .(alpha, beta_, gamma, delta)]), add = TRUE)
	if (i %% 200 == 0)
		print(i)
}
abline(h = 1, col = "#CD212A")


curve(shape2(x, params_full[1, .(alpha, beta_, gamma, delta, phi)]), to = 3, ylim = c(-0.5, 2))
for (i in 1:4000)
{
	curve(shape2(x, params_full[i, .(alpha, beta_, gamma, delta, phi)]), add = TRUE)
	print(i)
}
abline(h = 0, col = "#CD212A")
