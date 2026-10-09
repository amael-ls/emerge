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

if (!isTRUE(all.equal(key(tree_dt), c("speciesName_sci", "plot_id"))))
	stop("The keys in 00_run were set to species and plot, respectively")

## Seeds that were used to run the models (full and submodel)
seed_dt = readRDS(paste0(path_data, "ls_seeds.rds"))

## Load stan models
hierarchical = cmdstan_model(paste0(path_models, "hierarchical.stan"))

# options(max.print = 500)
# N_test = 1e6

# mu_c = rbeta(N_test, 2, 2);
# kappa_c = rgamma(N_test, shape = 2, rate = 0.1);

# c_beta = rbeta(N_test, mu_c * kappa_c, (1 - mu_c) * kappa_c);

# c = 0.6 + 0.4*c_beta

# hist(c)

# --------------------------------------------------------------------------------------
# --------------------    Run the full model, for the 14 species    --------------------
# --------------------------------------------------------------------------------------
(sp = ls_species[9])
# for (sp in ls_species[1:2])
{
	print(paste("Running", sp))

	filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_hierarchical")
	le_cid = seed_dt[.(sp), full]

	if (sp == "Quercus sp.")
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")

	indices = tree_dt[.(sp)][, .(start_index = min(.I), end_index = max(.I)), by = plot_id]

	# Subset to targeted species for stanData
	stanData = list(
		N = tree_dt[.(sp)][, .N],
		N_plot = tree_dt[.(sp)][, length(unique(plot_id))],
		start_index = indices[, start_index],
		end_index = indices[, end_index],
		bole_volume_m3 = tree_dt[.(sp)][, bole_volume_m3],
		total_volume_m3 = tree_dt[.(sp)][, total_volume_m3]
	)

	# Run full model
	if (file.exists(paste0(path_output, filename, ".rds")))
	{
		fit = readRDS(paste0(path_output, filename, ".rds"))
	} else {
		fit = hierarchical$sample(data = stanData, chains = n_chains,
			parallel_chains = min(n_chains, 4), seed = le_cid,
			iter_warmup = 1500, iter_sampling = 1000, max_treedepth = 10)

		# Save results
		fit$save_output_files(dir = path_output, basename = filename, random = FALSE)
		saveRDS(fit, paste0(path_output, filename, ".rds"))
	}

	ll = fit$diagnostic_summary()

	# Posterior of c (does not really makes sense when there are little plots...)
	pred_c = apply(X = fit$draws("c"), MARGIN = 3, FUN = mean) |>
		as.data.table() |>
		setnames(new = "c")
	hist(pred_c[, c])

	# d_c = density(pred_c[, c], n = 512)
	# density_dt = data.table(x_c = d_c$x, y_c = d_c$y)

	# filename = paste0(path_pgfplotsfig, "hierarchical_c_", sp, ".csv") |>
	# 	stri_replace(regex = " ", replacement = "-")
	# if (!file.exists(filename))
	# 	fwrite(density_dt, filename)

	# Compare the hierarchical and non-hierarchical predictions
	pred_dt = data.table(obs = stanData$total_volume_m3,
		hierarchical = apply(X = fit$draws("v_gen_mean"), MARGIN = 3, FUN = mean))

	filename = paste0(stri_replace(str = sp, regex = " ", replacement = "-"), "_fullmodel_theta")
	if (sp == "Quercus sp.")
		filename = stri_replace(str = filename, regex = "\\.", replacement = "")
	fit = readRDS(paste0(path_output, filename, ".rds"))

	pred_dt[, fullmodel := apply(X = fit$draws("v_gen_mean"), MARGIN = 3, FUN = mean)]


	plot(pred_dt[, hierarchical], pred_dt[, obs], pch = 19, cex = 0.75, xlab = "Predicted",
		ylab = "Observed", axes = FALSE, main = sp)
	points(pred_dt[, fullmodel], pred_dt[, obs], pch = 19, cex = 0.5, col = "#CD212A")
	abline(a = 0, b = 1, lty = "dashed", col = "#9A9A9A", lwd = 0.8)
	axis(1)
	axis(2, las = 1)

	# filename = paste0(path_pgfplotsfig, "hierarchical-vs-fullmodel_", sp, ".csv") |>
	# 	stri_replace(regex = " ", replacement = "-")
	# if (!file.exists(filename))
	# 	fwrite(pred_dt, filename)
}
ll

#### Other information
## Number of individuals per plot etc, by species
probs = c(0, 0.25, 0.5, 0.75, 1)
plot_counts = tree_dt[, .(n_trees = .N), by = .(speciesName_sci, plot_id)]

dt = plot_counts[, c(
	list(N = sum(n_trees), n_plots = .N),
	setNames(as.list(quantile(n_trees, probs, type = 1)), paste0("q", probs * 100))), by = speciesName_sci]

if (!file.exists(paste0(path_pgfplotstable, "n_plots.csv")))
	fwrite(dt, paste0(path_pgfplotstable, "n_plots.csv"))
