#### Aim of scripts: Contains global variables used in other scripts
## Paths
path_data = "./data/"
path_output = "./results/"
path_models = "./stan_models/"
path_pgfplotsfig = "./pgfplots_figs/"
path_pgfplotstable = "./pgfplots_tables/"

## Variables for running models
n_chains = 4

## Table 5 from Vallet2006
vallet_dt = data.table(speciesName_sci = c("Abies alba", "Fagus sylvatica", "Picea abies",
		"Pinus pinaster", "Pinus sylvestris", "Pseudotsuga menziesii", "Quercus petraea"),
	alpha = c(0.550, 0.395, 0.631, 0.235, 0.297, 0.534, 0.471),
	sd_alpha = c(0.015, 0.01, 0.007, 0.041, 0.025, 0.01, 0.014),
	beta = c(-7.49e-4, 2.66e-4, -9.46e-4, 9.7e-4, 3.18e-4, -5.3e-4, -3.45e-4),
	sd_beta = c(3.9e-5, 4.9e-5, 7.2e-5, 4.06e-4, 1.4e-4, 9.7e-5, 1.3e-5),
	gamma = c(0.277, 0.421, NA, 0.396, 0.384, NA, 0.377),
	sd_gamma = c(0.034, 0.025, NA, 0.057, 0.058, NA, 0.031),
	delta = c(NA, 45.4, NA, 198.8, 204.0, 56.6, NA),
	sd_delta = c(NA, 4, NA, 40, 26.6, 13.8, NA),
	var_res = c(0.0031, 0.0036, 0.0023, 0.0079, 0.0028, 0.0041 , 0.004),
	n_params = c(3, 4, 2, 4, 4, 3, 3), key = "speciesName_sci")

## Parameters and range, verified from Longuetaud
longuetaud_pars = data.table(
	genus = c("betula", "carpinus", "fagus", "fraxinus", "quercus",
		"abies", "larix", "picea", "pinus", "pseudotsuga"),
	b1 = rep(6.83, 10),
	b2 = c(1.016, 1.133, 0.842, 1.112, 1.320,
		1.496, 0.966, 0.953, 0.566, 0.744),
	b3 = rep(1.009, 10),
	b4 = c(0.448, 0.822, 0.627, 0.470, 0.245,
		0.566, 0.454, 0.473, 0.336, 0.783),
	min_d = c(12, 8, 7, 8, 7,
		7, 9, 8, 7, 7),
	max_d = c(31, 38, 79, 43, 102,
		90, 42, 67, 71, 34),
	min_h = c(14.7, 11.3, 5.5, 11.5, 6,
		6.9, 9, 7.7, 5.1, 7),
	max_h = c(28.5, 30, 42.5, 26.7, 40,
		37.8, 25.6, 41.3, 33.4, 29.5),
	is_broadleaf = rep(c(TRUE, FALSE), each = 5), key = "genus")
