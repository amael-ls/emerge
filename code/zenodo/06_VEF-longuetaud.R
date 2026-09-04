#### Aim of script: Compute averaged VEF from Longuetaud et al.
## Comments:
# Bibliography
#	Longuetaud, Fleur, Philippe Santenoise, Frédéric Mothe, Tristan Senga Kiessé, Michaël Rivoire,
#		Laurent Saint-André, Nina Ognouabi, and Christine Deleuze. 2013.
#		Modeling Volume Expansion Factors for Temperate Tree Species in France.
#		Forest Ecology and Management 292: 111–21. https://doi.org/10.1016/j.foreco.2012.12.023.

#### Packages needed to reproduce the study
renv::restore()

library(data.table)
library(cmdstanr)
library(stringi)

#### Load data
## Tool functions
source("./tool_functions.R")

## Other local functions
vef_fct = function(dbh, height, pars)
{
	if (!is.data.table(pars))
		stop("Please provide a data table!")

	b1 = pars[, b1]
	b2 = pars[, b2]
	b3 = pars[, b3]
	b4 = pars[, b4]

	int_broadleaf = ifelse(pars[, is_broadleaf], 1, 0)

	vef = 1 + exp(b2 * (b1 - dbh)) + exp(b3*int_broadleaf + b4)*dbh/height^2

	return(vef)
}

vef_mean = function(pars)
{
	f = function(dbh, height)
		vef_fct(dbh, height, pars)

	res = pracma::integral2(
		fun = f,
		xmin = pars$dbh_cm_q025,
		xmax = pars$dbh_cm_q975,
		ymin = pars$height_q025,
		ymax = pars$height_q975
	)

	areaDomain = (pars$dbh_cm_q975 - pars$dbh_cm_q025) * (pars$height_q975 - pars$height_q025)

	return(res$Q/areaDomain)
}

## Global variables (paths and others)
source("./global_variables.R")

## Tree data (14 species)
tree_dt = readRDS(paste0(path_data, "tree_dt.rds"))[origin %in% c("inra", "emerge")][
	!is.na(circumference_m) & !is.na(height)]

#### Prepare data
## Compute genus-specific quantiles
tree_dt[, genus := tolower(
	stri_sub(speciesName_sci, to = stri_locate_first(speciesName_sci, regex = " ")[, "start"] - 1))]
tree_dt[, dbh_cm := 100*circumference_m/pi]
genus = tree_dt[, as.list(lapply(.SD, quantile, probs = c(0.025, 0.975))), by = "genus",
	.SDcols = c("dbh_cm", "height")]
genus[, vars := rep(c("q025", "q975"), .N/2)]

genus = dcast.data.table(data = genus, formula = genus ~ vars, value.var = c("dbh_cm", "height"))

## Parameters and range, verified
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
	is_broadleaf = rep(c(TRUE, FALSE), each = 5),
	key = "genus"
)

longuetaud_pars = merge.data.table(longuetaud_pars, genus, by = "genus")

#### Compute average VEF from Longuetaud2013
longuetaud_pars[, mean_vef := vef_mean(.SD), by = genus]

if (!file.exists(paste0(path_output, "longuetaud_VEF.rds")))
	saveRDS(longuetaud_pars, paste0(path_output, "longuetaud_VEF.rds"))
