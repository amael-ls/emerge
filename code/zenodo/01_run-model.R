#### Aim of script: Run the models for the 14 species
## Comments
# This file is to run the species-specific models only. For the group models, see 02_run-model_groups.R
# On a not too old computer (say from 2020), this script should take less than an hour to run. The longest
# species to paramtrise are, without surprise, the most abundant.
# For a quick test, use sp = "Fraxinus excelsior", it is the fastest to run

## Packages needed to reproduce the study
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
ls_species = tree_dt[, unique(speciesName_sci)]

# Seeds that were used to run the models (full and submodel)
seed_dt = readRDS(paste0(path_data, "ls_seeds.rds"))

# Load stan models
fullmodel = cmdstan_model(paste0(path_models, "fullmodel.stan"))
submodel = cmdstan_model(paste0(path_models, "submodel.stan"))



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



#### DRAFT ZONE --------------------------------------------------------------------------------
# Extract posterior draws
draws_df = setDT(posterior::as_draws_df(fit$draws(variables = c("c", "j", "k", "m", "n", "s", "tau", "theta"))))
setnames(draws_df, old = c(".chain", ".iteration", ".draw"), new = c("chain", "iteration", "draw"))

colours = MetBrewer::met.brewer("Egypt")[1:4]

plot(1:1000, draws_df[1:1000, tau], col = colours[1], type = "l", ylim = c(0, 2))
lines(1001:2000 - 1000, draws_df[1001:2000, tau], col = colours[2])
lines(2001:3000 - 2000, draws_df[2001:3000, tau], col = colours[3])
lines(3001:4000 - 3000, draws_df[3001:4000, tau], col = colours[4])
legend("topright", legend = paste("Chain", 1:4), fill = colours)

plot(1:1000, draws_df[1:1000, j], col = colours[1], type = "l")
lines(1001:2000 - 1000, draws_df[1001:2000, j], col = colours[2])
lines(2001:3000 - 2000, draws_df[2001:3000, j], col = colours[3])
lines(3001:4000 - 3000, draws_df[3001:4000, j], col = colours[4])
legend("topright", legend = paste("Chain", 1:4), fill = colours)

mean_params = draws_df[, lapply(.SD, FUN = mean), .SDcols = c("c", "j", "k", "m", "n", "s", "tau"), by = chain]

r_func = function(x, pars)
{
	return((pars[, m] - pars[, c]) * exp(pars[, j] - pars[, k]*x) * (pars[, k]*x/pars[, j])^pars[, j] +
		pars[, c] - (pars[, c] - pars[, n])*exp(-pars[, s]*x));
}

curve(r_func(x, mean_params[chain == 1, .(c, j, k, m, n, s)]), col = colours[1], to = 8)
curve(r_func(x, mean_params[chain == 2, .(c, j, k, m, n, s)]), col = colours[2], add = TRUE)
curve(r_func(x, mean_params[chain == 3, .(c, j, k, m, n, s)]), col = colours[3], add = TRUE)
curve(r_func(x, mean_params[chain == 4, .(c, j, k, m, n, s)]), col = colours[4], add = TRUE)
legend("topright", legend = paste("Chain", 1:4), fill = colours)

plot(draws_df$c, draws_df$j, xlab = "c", ylab = "j", pch = 16, cex = 0.4)
f = function(x, linMod)
{
	qq = coefficients(linMod)
	return(qq[1] + qq[2]*x + qq[3]*x^2)
}
lin = lm(formula = j ~ 1 + c + I(c^2), data = draws_df)
curve(f(x, lin), add = TRUE, lwd = 2, col = "#CD212A")

# 1. Take draws at both ends of the c-j banana
lim_c = draws_df[, quantile(c, prob = c(0.025, 0.975))]
dt = draws_df[c < lim_c["2.5%"] | c > lim_c["97.5%"]]

plot(dt$c, dt$j, xlab = "c", ylab = "j", pch = 16, cex = 0.4)
curve(f(x, lin), add = TRUE, lwd = 2, col = "#CD212A")

curve(r_func(x, dt[1, .(c, j, k, m, n, s)]), to = 8)
for (i in seq_len(dt[, .N]))
{
	if (dt[i, c] < lim_c["2.5%"])
		curve(r_func(x, dt[i, .(c, j, k, m, n, s)]), add = TRUE, col = "#FAB255")
	if (dt[i, c] > lim_c["97.5%"])
		curve(r_func(x, dt[i, .(c, j, k, m, n, s)]), add = TRUE, col = "#0F7BA2")
}

## Check mu_1
mu_1 = function(x, pars)
	return((pars[, m] - pars[, c]) * exp(pars[, j] - pars[, k]*x) * (pars[, k]*x/pars[, j])^pars[, j])

mu_1_alt = function(x, pars)
	return((pars[, m] - pars[, c]) * exp(pars[, j] * (1 - x/pars[, tau])) * (x/pars[, tau])^pars[, j])

curve(mu_1_alt(x, dt[1, .(c, j, k, m, n, s, tau)]), to = 8, ylim = c(0, 0.2))
for (i in seq_len(dt[, .N]))
{
	if (dt[i, c] < lim_c["2.5%"])
		curve(mu_1_alt(x, dt[i, .(c, j, k, m, n, s, tau)]), add = TRUE, col = "#FAB255")
	if (dt[i, c] > lim_c["97.5%"])
		curve(mu_1_alt(x, dt[i, .(c, j, k, m, n, s, tau)]), add = TRUE, col = "#0F7BA2")
}
legend("topright", legend = c("low", "high"), fill = c("#FAB255", "#0F7BA2"))


fun = function(x, j)
	return(exp(j*(1 - x)))

DL_0 = function(x, j)
	return(exp(j) - j*exp(j)*x + j^2*exp(j)/2*x^2)

curve(fun(x, 0.75), from = -0.75, to = 0.57)
curve(DL_0(x, 0.75), add = TRUE, col = "#CD212A", lty = "dashed")

#### END DRAFT ZONE --------------------------------------------------------------------------------



################################################
#### Stuff to find a prior
Ex = 1.2
Vx = 3

sigma = sqrt(log(Vx/Ex^2 + 1))
mu = log(Ex) - 0.5*sigma^2

sim = rlnorm(1e5, mu, sigma)
mean(sim)
var(sim)

curve(dlnorm(x, mu, sigma), to = 5)
