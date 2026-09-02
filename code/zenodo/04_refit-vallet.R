#### Aim of script: Refit Vallet 2006 in Bayesian with our data

## Comments:
# The original publication is:
# Vallet, P., Dhôte, J.-F., Moguédec, G. L., Ravart, M., & Pignard, G. (2006).
# 	Development of total aboveground volume equations for seven important forest tree species in France.
# 	Forest Ecology and Management, 229(1-3), 98-110. https://doi.org/10.1016/j.foreco.2006.03.013
# 
# We found that the predictions are quite similar, except for Abies alba. This might originate from
# 	our dataset which is well extended by Swiss data for that species.

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

# Table 5 from Vallet2006
vallet_dt = data.table(speciesName_sci = c("Abies alba", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris", "Pseudotsuga menziesii", "Quercus petraea"),
	alpha = c(0.550, 0.395, 0.631, 0.235, 0.297, 0.534, 0.471),
	sd_alpha = c(0.015, 0.01, 0.007, 0.041, 0.025, 0.01, 0.014),
	beta = c(7.49e-4, 2.66e-4, -9.46e-4, 9.7e-4, 3.18e-4, -5.3e-4, -3.45e-4),
	sd_beta = c(3.9e-5, 4.9e-5, 7.2e-5, 4.06e-4, 1.4e-4, 9.7e-5, 1.3e-5),
	gamma = c(0.277, 0.421, NA, 0.396, 0.384, NA, 0.377),
	sd_gamma = c(0.034, 0.025, NA, 0.057, 0.058, NA, 0.031),
	delta = c(NA, 45.4, NA, 198.8, 204.0, 56.6, NA),
	sd_delta = c(NA, 4, NA, 40, 26.6, 13.8, NA),
	var_res = c(0.0031, 0.0036, 0.0023, 0.0079, 0.0028, 0.0041 , 0.004),
	n_params = c(3, 4, 2, 4, 4, 3, 3), key = "speciesName_sci")

# Loading training dataset
tree_dt = readRDS(paste0(path_data, "tree_dt_14species.rds"))
tree_dt = tree_dt[(!is.na(circumference_m)) & (!is.na(height))]

# Species parametrised in Vallet et al. 2006
ls_species = c("Abies alba", "Fagus sylvatica", "Picea abies", "Pinus pinaster", "Pinus sylvestris",
	"Pseudotsuga menziesii", "Quercus petraea")
tree_dt = tree_dt[ls_species]
setkey(tree_dt, speciesName_sci)

## Refit Vallet 2006 in Bayesian with our dataset
# Compile
vallet_bayesian = cmdstan_model(paste0(path_models, "vallet.stan"))

# Run Bayesian fit and compare with original Freq fit
for (sp in ls_species)
{
	sp_filename = stri_replace(str = sp, replacement = "-", regex = " ")
	filename = paste0(path_output, sp_filename, "_vallet.rds")
	is_douglas = ifelse(sp == "Pseudotsuga menziesii", 1, 0) # Boolean style compatible with Stan language

	if (!file.exists(filename))
	{
		stanData = list(
			N = tree_dt[sp, .N],
			N_params = vallet_dt[sp, n_params],
			is_douglas = is_douglas,

			circumference_cm = tree_dt[sp, 100*circumference_m], # Was in cm in Vallet2006!
			height = tree_dt[sp, height],

			total_volume_m3 = tree_dt[sp, total_volume_m3]
		)

		results = vallet_bayesian$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
			max_treedepth = 12)

		results$save_output_files(dir = path_output, basename = paste0(sp_filename, "_vallet"), random = FALSE)
		saveRDS(results, filename)
	} else {
		results = readRDS(filename)
	}

	bayes_params = getParams(results, paste0("vec_params[", 1:vallet_dt[sp, n_params], "]"), type = "mean")
	
	if (vallet_dt[sp, n_params] == 2)
	{
		names(bayes_params) = c("alpha", "beta")
		freq_params = vallet_dt[sp, c(alpha, beta)]
		names(freq_params) = names(bayes_params)
	}
	
	if (vallet_dt[sp, n_params] == 3)
	{
		names(bayes_params) = c("alpha", "beta", "gamma")
		freq_params = vallet_dt[sp, c(alpha, beta, gamma)]
		names(freq_params) = names(bayes_params)

		if (is_douglas)
		{
			names(bayes_params) = c("alpha", "beta", "delta")
			freq_params = vallet_dt[sp, c(alpha, beta, delta)]
			names(freq_params) = names(bayes_params)
		}
	}
	
	if (vallet_dt[sp, n_params] == 4)
	{
		names(bayes_params) = c("alpha", "beta", "gamma", "delta")
		freq_params = vallet_dt[sp, c(alpha, beta, gamma, delta)]
		names(freq_params) = names(bayes_params)
	}
	
	tree_dt[sp, bayes_vallet :=
		form_vallet(100*circumference_m, height, bayes_params, vallet_dt[sp, n_params], is_douglas)]
	tree_dt[sp, freq_vallet :=
		form_vallet(100*circumference_m, height, freq_params, vallet_dt[sp, n_params], is_douglas)]

	## Plot Bayes vs Freq
	plot(tree_dt[sp, freq_vallet], tree_dt[sp, bayes_vallet], pch = 19, cex = 0.55,
		axes = FALSE, xlab = "Frequentist", ylab = "Bayesian")
	abline(a = 0, b = 1, lwd = 1.5, col = "#CD212A")
	axis(1)
	axis(2, las = 1)
}


