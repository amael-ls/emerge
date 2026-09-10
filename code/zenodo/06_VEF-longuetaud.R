#### Aim of script: Compute averaged VEF from Longuetaud et al.
## Comments:
# Bibliography
#	Longuetaud, Fleur, Philippe Santenoise, Frédéric Mothe, Tristan Senga Kiessé, Michaël Rivoire,
#		Laurent Saint-André, Nina Ognouabi, and Christine Deleuze. 2013.
#		Modeling Volume Expansion Factors for Temperate Tree Species in France.
#		Forest Ecology and Management 292: 111–21. https://doi.org/10.1016/j.foreco.2012.12.023.

#### Packages needed to reproduce the study
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

longuetaud_pars = merge.data.table(longuetaud_pars, genus, by = "genus")

#### Compute average VEF from Longuetaud2013
longuetaud_pars[, mean_vef := vef_mean(.SD), by = genus]

if (!file.exists(paste0(path_output, "longuetaud_VEF.rds")))
	saveRDS(longuetaud_pars, paste0(path_output, "longuetaud_VEF.rds"))
