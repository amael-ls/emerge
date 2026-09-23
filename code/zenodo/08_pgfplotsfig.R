#### Aim of script: Generate the files for pgfplots Latex
## Comments:
# This file is not necessary for the study itself. It only creates csv files
#	to be read by latex using pgfplots
#
# When I am using the simplified model, and despite it is a submodel of the general case,
#	it seems that I need some cautiousness on the computation with Lambert function, especially
#	with \varepsilon... Check p.12--13 Notebook 4.
#

#### Load packages
library(data.table)
library(stringi)
library(mgcv)

#### Load common data and tool functions
## Tool functions
source("./tool_functions.R")

## Global variables (paths and others)
source("./global_variables.R")

## Other functions
# Function to pred Vtot, with uncertainty from the parameters + residuals
pred_fct_Vtot = function(sp, stanData_gen, path_models, path_output, is_simplif = FALSE,
	n_points = 5, woodstock_seed = 1969 - 08 - 18)
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

	if (sp == "Pinus-uncinata")
	{
		genQ = cmdstanr::cmdstan_model(paste0(path_models, "pinus_uncinata-genQ.stan"))
		filename = paste0(path_output, stringi::stri_replace(str = sp, regex = " ", replacement = "-"),
			"_logit.rds")
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

# Second derivative of mu_1, evaluated at x = tau = j/k (full model)
mu_2nd = function(k, j, m, c, s, n)
	return(-k^2/j*(m - c) - s^2*(c - n)*exp(-s*j/k))

# Second derivative of mu_1, evaluated at x = tau = j/k (submodel equivalent)
mu_2nd_submodel = function(k, m, c, n)
{
	gamma = k*(m - c)*exp(1)
	delta = n - c
	return(-k*gamma*exp(-k*(1/k - delta/gamma)))
}

# Second central moment of mu_1, i.e., variance computed around the mean (not tau) as it is a gamma fct
var_mu = function(tau, j)
	return(tau^2*(j + 1)/j^2)

# Third central moment of mu_1, i.e., skewness
skewness = function(j)
	return(2/sqrt(j + 1))

# Function for Pinus uncinata
mu_logit_fct = function(x, pars)
	return (inv_logit(pars["logit_alpha"] + exp(-pars["beta_"]*x) * (pars["gamma"]*x + pars["delta"])))

vtot_logit_fct = function(x, pars)
	return (x/mu_logit_fct(x, pars))

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

for (sp in ls_species)
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

	if (sp == "Pinus uncinata")
	{
		warning("Pinus uncinata is treated separately. I still think there is not enough info in the data used to parametrise Pinus uncinata")
	}

	temp = pred_fct_Vtot(sp = sp_filename, stanData_gen = stanData, path_models = path_models,
		path_output = path_output, is_simplif = simplif)

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

## Second derivative of mu for Pinus uncinata, evaluated at xM
pinus_uncinata_mu_2nd = function(pars)
{
	logit_alpha <- pars["logit_alpha"]
	beta_ <- pars["beta_"]
	gamma <- pars["gamma"]
	delta <- pars["delta"]
	zstar <- logit_alpha + (gamma / beta_) * exp(-1 + beta_ * delta / gamma)
	fstar <- inv_logit(zstar)
	return (-beta_ * gamma * exp(-1 + beta_ * delta / gamma) * fstar * (1 - fstar))
}

#### Compute threshold with Lambert function, plot ratio on one Fig.
## Common variables
n_points = 750
vbole = seq(0, 7, length.out = n_points) # 7 is a good compromise for the 6 species I am interestied in

threshold_dt = data.table(speciesName_sci = tree_dt[, unique(speciesName_sci)], simplif = FALSE,
	x1 = NA_real_, x2 = NA_real_, wp = NA_real_, c = NA_real_, eps = NA_real_, div_by = NA_integer_,
	delta_mc = NA_real_, diff_eps = NA_real_, key = "speciesName_sci")

params_dt = data.table(speciesName_sci = tree_dt[, unique(speciesName_sci)], c = NA_real_,
	j = NA_real_, k = NA_real_, m = NA_real_, n = NA_real_, s = NA_real_, tau = NA_real_, key = "speciesName_sci")

fct_output = matrix(data = NA_real_, nrow = threshold_dt[, .N] + 1, ncol = n_points)
rownames(fct_output) = c("vbole", threshold_dt[, speciesName_sci])
fct_output["vbole", ] = vbole

loaded = FALSE

if (!file.exists(paste0(path_output, "lambert_calculus.rds")))
{
	for (sp in tree_dt[, unique(speciesName_sci)])
	{
		print(paste("Doing species", sp))

		## Get params
		if ((sp == "Fraxinus excelsior") || (sp == "Pinus uncinata"))
			next # Simplified species, done after

		sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")
		if (sp == "Quercus sp.")
			sp_filename = stri_replace(str = sp_filename, replacement = "", regex = "\\.")

		filename = paste0(path_output, sp_filename, "_fullmodel_theta", ".rds")
		fit = readRDS(filename)

		paramsVec = getParams(model_cmdstan = fit, params_names = c("c", "j", "k", "m", "n", "s", "tau"),
			type = "mean")
		rm(fit)

		m = paramsVec["m"]
		j = paramsVec["j"]
		c = paramsVec["c"]
		k = paramsVec["k"]

		epsilon = 2/100*c # 2% of the asymptotic rate
		div_by = 1

		if (epsilon > (m - c))
		{
			warning(paste("The condition epsilon < m - c is not respeced for", sp))
			count = 0
			while (epsilon > (m - c) && count < 10)
			{
				count = count + 1
				epsilon = epsilon/2
			}

			if (epsilon > (m - c))
			{
				warning(paste("Could not solve the problem for", sp, "by setting epsilon = epsilon/1024"))
				next
			}
			warning(paste("Divided epsilon by", 2^count, "for", sp))
			div_by = 2^count
		}

		xi = -exp(-1) * (epsilon/(m - c))^(1/j)
		wp = -j/k * lamW::lambertW0(xi)
		wm = -j/k * lamW::lambertWm1(xi)

		if (wm < j/k)
			warning("My analytical calculus showed that this should not be possible...")

		s = paramsVec["s"]
		n = paramsVec["n"]
		tau = paramsVec["tau"]

		params_dt[sp, c("c", "j", "k", "m", "n", "s", "tau") := as.list(paramsVec)]

		x2 = 1/s*(log(c - n) - log(epsilon))

		diff = abs(pred_ratio(wm, paramsVec) - c) - epsilon

		threshold_dt[sp, c("x1", "x2","wp", "c", "eps", "div_by", "delta_mc", "diff_eps") :=
			.(wm, ..x2, ..wp, ..c, epsilon, ..div_by, (..m - ..c), diff)]

		fct_output[sp, ] = pred_ratio(vbole, paramsVec)
	}
} else {
	threshold_dt = readRDS(paste0(path_output, "lambert_calculus.rds"))
	fct_output = readRDS(paste0(path_output, "ratio_output.rds"))
	params_dt = readRDS(paste0(path_output, "avg_params.rds"))
	loaded = TRUE
}

## Compute threshold with the Lambert function for species using submodel (Fraxinus excelsior)
if (!loaded)
{
	# Fraxinus excelsior
	sp = "Fraxinus excelsior"
	print(paste("Doing species", sp))

	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")

	filename = paste0(path_output, sp_filename, "_submodel", ".rds")
	fit = readRDS(filename)

	paramsVec_simplif = getParams(model_cmdstan = fit,
		params_names = c("alpha", "beta_", "gamma", "delta"), type = "mean")

	rm(fit)

	alpha = paramsVec_simplif["alpha"] # Frax: 0.7159374
	beta = paramsVec_simplif["beta_"]  # Frax: 5.184622
	gamma = paramsVec_simplif["gamma"] # Frax: 2.829555
	delta = paramsVec_simplif["delta"] # Frax: -0.02979472

	epsilon = 2*alpha/100
	epsilon < gamma/beta*exp(beta*delta/gamma - 1)

	xi = -beta*epsilon/gamma*exp(-beta*delta/gamma)

	wp = -1/beta * (lamW::lambertW0(xi) + beta*delta/gamma)
	wm = -1/beta * (lamW::lambertWm1(xi) + beta*delta/gamma)

	paramsVec = c(
		c = unname(paramsVec_simplif["alpha"]),
		j = 1,
		k = unname(paramsVec_simplif["beta_"]),
		m = unname(exp(-1)*paramsVec_simplif["gamma"]/paramsVec_simplif["beta_"] +
			paramsVec_simplif["alpha"]),
		n = unname(paramsVec_simplif["delta"] + paramsVec_simplif["alpha"]),
		s = unname(paramsVec_simplif["beta_"])
	)
	fct_output[sp, ] = pred_ratio(vbole, paramsVec)

	params_dt[sp, c("c", "j", "k", "m", "n", "s") := as.list(paramsVec)]

	diff = abs(pred_ratio(wm, paramsVec) - alpha) - epsilon

	threshold_dt[sp, c("x1", "wp", "c", "eps", "div_by", "delta_mc", "diff_eps") :=
		.(wm, ..wp, ..alpha, epsilon, 1, gamma/beta*exp(-1), diff)]

	# ------------------------------------------------------------------------------
	# Pinus uncinata treated separatly
	sp = "Pinus uncinata"
	print(paste("Doing species", sp))

	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")

	filename = paste0(path_output, sp_filename, "_logit", ".rds")
	fit = readRDS(filename)

	pars = getParams(model_cmdstan = fit,
		params_names = c("alpha", "beta_", "gamma", "delta", "logit_alpha"), type = "mean")

	rm(fit)

	alpha = pars["alpha"] # Pinus uncinata: 0.8848285
	logit_alpha = pars["logit_alpha"] # Pinus uncinata: 2.104958
	beta = pars["beta_"]  # Pinus uncinata: 1.154964
	gamma = pars["gamma"] # Pinus uncinata: 0.6170549
	delta = pars["delta"] # Pinus uncinata: -0.6233957

	epsilon = 1.25*alpha/100
	(epsilon_logit = logit(epsilon + alpha) - logit_alpha)
	epsilon_logit < gamma/beta*exp(beta*delta/gamma - 1)

	xi = -beta*epsilon_logit/gamma*exp(-beta*delta/gamma)

	wp = -1/beta * (lamW::lambertW0(xi) + beta*delta/gamma)
	wm = -1/beta * (lamW::lambertWm1(xi) + beta*delta/gamma)

	fct_output[sp, ] = mu_logit_fct(vbole, pars)

	diff = abs(mu_logit_fct(wm, pars) - alpha) - epsilon

	threshold_dt[sp, c("x1", "wp", "c", "eps", "div_by", "delta_mc", "diff_eps") :=
		.(wm, ..wp, ..alpha, epsilon, 1, gamma/beta*exp(-1), diff)]

	saveRDS(fct_output, paste0(path_output, "ratio_output.rds"))
	saveRDS(params_dt, paste0(path_output, "avg_params.rds"))
}

## Add curvature at x_M (i.e., second derivative) at maximum x
if (!loaded)
{
	params_dt[, curve := mu_2nd(k, j, m, c, s, n), by = speciesName_sci]
	params_dt["Fraxinus excelsior", curve := mu_2nd_submodel(k, m, c, n)]
}

if (loaded)
{
	sp = "Pinus uncinata"
	print(paste("Doing species", sp))

	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")

	filename = paste0(path_output, sp_filename, "_logit", ".rds")
	fit = readRDS(filename)

	pars = getParams(model_cmdstan = fit,
		params_names = c("alpha", "beta_", "gamma", "delta", "logit_alpha"), type = "mean")

	rm(fit)
	params_dt["Pinus uncinata", curve := pinus_uncinata_mu_2nd(pars)]
}

if (!loaded)
{
	## Add ratio (m - c)/epsilon to know how many times peak above epsilon
	threshold_dt[, eps_mc := delta_mc/eps]
	threshold_dt[, .(speciesName_sci, x1, wp, eps_mc)]

	threshold_dt = merge.data.table(x = threshold_dt,
		y = params_dt[, .(speciesName_sci, curve, percent_c = 100*(m - c)/c, x3 = j/k)], by = "speciesName_sci")
	
	## Add how much percentage of c lies in m - c, and x3 the location of the max
	# Modify manually for Fraxinus excelsior (see notebook 4, p. 55 03 August 2026)
	p = params_dt["Fraxinus excelsior", .(k, n, m, c)]
	threshold_dt["Fraxinus excelsior", x3 := -1/p[, k] * ((p[, n] - p[, c])*exp(1)/(p[, m] - p[, c]) - 1)]

	# Modify manually for Pinus uncinata
	threshold_dt["Pinus uncinata", x3 := (pars["gamma"] - pars["beta_"] * pars["delta"]) / (pars["beta_"] * pars["gamma"])]
}

# Add the computation of the second central moment
params_dt[, var_mu := var_mu(tau, j), by = speciesName_sci]
params_dt["Fraxinus excelsior", var_mu := var_mu(j/k, j)] # Equals to 2*tau^2 = 2/beta^2
params_dt[, .(speciesName_sci, var_mu)]

# Add the computation of the third central moment
params_dt[, skewness := skewness(j), by = speciesName_sci]

if (!loaded)
	threshold_dt = merge.data.table(x = threshold_dt,
		y = params_dt[, .(speciesName_sci, var_mu, skewness)], by = "speciesName_sci")

## Export the data for pgfplots
# Write fct_ouput
filename = paste0(path_pgfplotsfig, "ratio-fct_all-sp.csv")
if (!file.exists(filename))
{
	temp = as.data.table(t(fct_output))
	fwrite(temp, filename)
}

## Write treshold_dt
filename = paste0(path_pgfplotstable, "thresholds.csv")
if (!file.exists(filename))
	fwrite(threshold_dt, filename, na = "NaN")

threshold_dt[percent_c < 4, .(speciesName_sci, percent_c, x1, x3)]
threshold_dt[, .(speciesName_sci, curve, var_mu)]

if (!file.exists(paste0(path_output, "lambert_calculus.rds")))
	saveRDS(threshold_dt, paste0(path_output, "lambert_calculus.rds"))

#### Draft zone for Pinus uncinata (to check location parameters)
plot(fct_output["vbole", ], fct_output["Pinus uncinata", ], type = "l")
abline(v = threshold_dt["Pinus uncinata", x3])
abline(h = mu_logit_fct(threshold_dt["Pinus uncinata", x3], pars))
points(tree_dt["Pinus uncinata", bole_volume_m3], tree_dt["Pinus uncinata", r], pch = 19, cex = 0.75, col = "#FAB255")

broad = c("Fagus sylvatica", "Faxinus excelsior", "Quercus petraea", "Quercus sp.")
plot(0, pch = "", xlim = c(0, 4), ylim = c(0.65, 0.95), axes = FALSE,
	xlab = "Bole volume", ylab = "Ratio")
axis(1)
axis(2, las = 1)
for (sp in ls_species)
{
	if (sp %in% broad)
		curve(pred_ratio(x, params_dt[.(sp)]), add = TRUE, col = "#FAB255", lwd = 4)

	if (!(sp %in% broad))
		curve(pred_ratio(x, params_dt[.(sp)]), add = TRUE, col = "#0F7BA2", lwd = 2)
}
