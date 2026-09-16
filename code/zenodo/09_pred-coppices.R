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

#### Predict volumes on Emerge data, and keep the residual for each draw (4000 per individual)
if (!file.exists(paste0(path_output, "v_gen_list.rds")) && !file.exists(paste0(path_output, "v_res_list.rds")))
{
	for (sp in ls_species[15:17])
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

## Reshape volumes and residuals
v_gen = rbindlist(l = v_gen_list)
opposite_res = rbindlist(l = v_res_list)

setnames(v_gen, old = as.character(1:4000), new = paste0("draws_", 1:4000))
setnames(opposite_res, old = as.character(1:4000), new = paste0("draws_", 1:4000))

## Posterior distribution of the bias
by_struct = v_gen[, lapply(.SD, mean), by = str_name, .SDcols = patterns("^draws_")]
by_struct_res = opposite_res[, lapply(.SD, mean), by = str_name, .SDcols = patterns("^draws_")]

## Quantiles residuals and average individual residuals
mean_res = by_struct_res[, {
	qs = apply(.SD, 1, quantile, probs = c(0, 0.025, 0.5, 0.975, 1)) # matrix: rows = quantiles, cols = str_name
	qs_dt = as.data.table(t(qs)) # Transpose
	setnames(qs_dt, c("q0", "q025", "q50", "q975", "q100"))
	cbind(str_name = str_name, res_mean = rowMeans(.SD), qs_dt)
}, .SDcols = patterns("^draws_")]

individual_res = opposite_res[, .(speciesName_sci, dataset, tree_id, str_name,
	opposite_res = rowMeans(.SD)), .SDcols = patterns("^draws_")]

## Save outputs
if (!file.exists(paste0(path_pgfplotstable, "residuals_structure.csv")))
	fwrite(mean_res, paste0(path_pgfplotstable, "residuals_structure.csv"))

if (!file.exists(paste0(path_pgfplotsfig, "individual_residuals_structure.csv")))
	fwrite(individual_res[, .(str_name, opposite_res)],
		paste0(path_pgfplotsfig, "individual_residuals_structure.csv"))

#### Compute kernel densities of 'biases'
by_struct_res = transpose(l = by_struct_res, make.names = "str_name")
setnames(by_struct_res, new = stri_replace_all(str = names(by_struct_res), regex = "-", replacement = "_"))

d_cws = density(by_struct_res[, coppice_with_standards], n = 512)
d_sl = density(by_struct_res[, single_layered], n = 512)
d_c = density(by_struct_res[, coppice], n = 512)

densities_dt = data.table(x_cws = d_cws$x, y_cws = d_cws$y,
	x_sl = d_sl$x, y_sl = d_sl$y,
	x_c = d_c$x, y_c = d_c$y)

## Save results
if (!file.exists(paste0(path_pgfplotsfig, "posterior_structure.csv")))
	fwrite(densities_dt, paste0(path_pgfplotsfig, "posterior_structure.csv"))

#### Boxplot of opposite res, this way, below zero means underestimated, and above zero means overestimated
aa = boxplot(opposite_res ~ str_name, data = individual_res,
	xlab = "Structure type",
	ylab = "Residuals")
abline(h = 0, lwd = 0.85, lty = "dashed", col = "#595859")

#### Check the type of trees that are rated outliers by boxplot function
## List individuals
individual_res[, is_outlier := opposite_res %in% boxplot.stats(opposite_res)$out, by = str_name]
ls_pb = individual_res[(is_outlier), .(speciesName_sci, dataset, tree_id, opposite_res)]
if (!all.equal(unique(ls_pb), ls_pb))
	stop("ls_pb is not uniquely defined")

## Get all information on these outliers
ls_pb = emerge[ls_pb, on = .(speciesName_sci, tree_id, dataset)]
ls_pb = ls_pb[, .(dataset, tree_id, str_name, speciesName_sci, circumference_m, height, taper_height,
	bole_volume_conic_m3, branch_volume, twig_volume, total_volume_m3, opposite_res, underestimated = opposite_res < 0)]
setkey(ls_pb, str_name, speciesName_sci)



# ------------------------------------------------------------------
# ----------------------    Try some stuff    ----------------------
# ------------------------------------------------------------------
#### Add soe information from the DB. This part of the script will not work for reviewers

os = Sys.info()[['sysname']]
mnt_point = "/mnt/local_share/"
if (os == "Linux" || os == "Darwin")
{
	if (!dir.exists(mnt_point))
		stop(paste0("The mounting point <", mnt_point, "> does not exist"))
	path_data = paste0(mnt_point, "data_orig/")
	if (!dir.exists(path_data))
		stop(paste0("Folder <", path_data, "> does not exist! Check mounting point and mounted folder"))
} else if (os == "Windows") {
	path_data = "//Del1509n015/2024_FairCarbon/data_orig/"
	if (!dir.exists(path_data))
		stop(paste0("Folder <", path_data, "> does not exist!"))
} else {
	stop(paste("Unknown Operating System:", os))
}

# Emerge data
filename_emerge = paste0(path_data, "EMERGE.RData")
if (!file.exists(filename_emerge))
	stop(paste0("The file <", filename_emerge, "> does not exist"))

load(filename_emerge)

setDT(emerge_2010_arbres)
setDT(emerge_bure09_arbres)

emerge_2010_arbres = emerge_2010_arbres[, .(id, insersion_1ere_branche_vivante_m, h_fourche1_m)]
emerge_bure09_arbres = emerge_bure09_arbres[, .(tree, hauteurbranchevivante_m, hauteurfourche1_m)]

setnames(emerge_2010_arbres, new = c("tree_id", "height_1st_alive_branch", "height_fork"))
setnames(emerge_bure09_arbres, new = c("tree_id", "height_1st_alive_branch", "height_fork"))

ll = rbindlist(l = list(emerge_2010 = emerge_2010_arbres, emerge_2009 = emerge_bure09_arbres), idcol = "dataset")
ls_pb = merge.data.table(ls_pb, ll, by = c("dataset", "tree_id"), all.x = TRUE)

hist(ls_pb[(underestimated), total_volume_m3/bole_volume_conic_m3])
hist(ls_pb[(underestimated), bole_volume_conic_m3/total_volume_m3])
hist(ls_pb[!(underestimated), total_volume_m3/bole_volume_conic_m3])


ls_pb[, fork_ratio := height_fork / height]

fork_only = ls_pb[!is.na(height_fork)]
fork_only[, .(mean_ratio = mean(fork_ratio), n = .N), by = underestimated]
# Higher ratio increases proba underestimated. Against my intuition...

wilcox.test(fork_ratio ~ underestimated, data = fork_only)

## Test a GLM
mod = glm(underestimated ~ 1 + height_fork, data = ls_pb[!is.na(height_fork)], family = binomial)

summary(mod)


a = coef(mod)[1]
b = coef(mod)[2]

# jitter the boolean underestimated so that overlapping points are visible (0 = FALSE, 1 = TRUE)
set.seed(woodstock_seed)
ls_pb[, bool_jit := as.integer(underestimated) + runif(.N, -0.03, 0.03)]

# Plot the stuff!
plot(ls_pb[!is.na(height_fork), height_fork], ls_pb[!is.na(height_fork), bool_jit],
	pch = 19, col = "#050505AA",
	xlab = "Height fork", ylab = "P(underestimated)",
	main = "Fitted logiit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)

## Same stuff, but I try a proxy for NAs in height, and then the ratio of fork/height (0 = at the foot, 1 = at the top)
ls_pb[, fork_proxy := fifelse(!is.na(height_fork), height_fork,
	fifelse(!is.na(height_1st_alive_branch), height_1st_alive_branch,
	taper_height))]

mod = glm(underestimated ~ 1 + I(fork_proxy/height),
	data = ls_pb, family = binomial)

summary(mod)


a = coef(mod)[1]
b = coef(mod)[2]

# Plot with the proxy
plot(ls_pb[, fork_proxy/height], ls_pb[, bool_jit],
	pch = 19, col = "#050505AA",
	xlab = "Ratio Hfork/Htot", ylab = "P(underestimated)",
	main = "Fitted logiit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)

# Same signal in both cases, the higher the fork, the more likely I underestimate.
#	It goes against my intuition! I thought larger crown would occur for lower forks!
# CHeck it

plot(ls_pb[, fork_proxy/height], ls_pb[, bole_volume_conic_m3/total_volume_m3],
	pch = 19, col = "#050505AA",
	xlab = "Ratio Hfork/Htot", ylab = "r") # That looks quite random and covers the whole space!

# So there is no rule of lower fork => larger crown => lower ratio bole/tot

## Check in abs residuals
ls_pb[abs(opposite_res) > 1.5] # These trees have high taper height!

## Try a model with taper height
mod = glm(underestimated ~ 1 + taper_height,
	data = ls_pb, family = binomial)

summary(mod)


a = coef(mod)[1]
b = coef(mod)[2]

# Plot with the proxy
plot(ls_pb[, taper_height], ls_pb[, bool_jit],
	pch = 19, col = "#050505AA",
	xlab = "Taper height", ylab = "P(underestimated)",
	main = "Fitted logiit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)

# This time it goes with my intuition! Lower taper => higher proba underestimation.
#	That mean that their are obvious discrepancies between taper height and fork height!
#	Are these stuff reliable?

plot(ls_pb[, taper_height], ls_pb[, height_fork],
	pch = 19, col = "#050505AA", xlab = "Taper height", ylab = "Fork height")

ls_pb[taper_height < height_fork] # There is none! And the relationship fork ~ taper is positive!

# Conclusion: it seems that nothing separates the underestimated from the overestimated...

#### Question: are the outliers different from the non-outliers?
# Comments: I do the same stuff, but on the whole dataset

## Add the information on trees as in ls_pb
full_data = emerge[individual_res, on = .(speciesName_sci, tree_id, dataset, str_name), nomatch = NA]
full_data = ll[full_data, on = .(dataset, tree_id), nomatch = NA]

## Correction for unrealistic values
full_data[height_1st_alive_branch == 0, height_1st_alive_branch := NA]

## Compute explanatoru variables
full_data[, fork_proxy := fifelse(!is.na(height_fork), height_fork, height_1st_alive_branch)]
full_data[, fork_ratio := fork_proxy / height]

## Run model
mod = glm(is_outlier ~ 1 + fork_ratio, data = full_data, family = binomial)
summary(mod)

## Plot stuff
set.seed(woodstock_seed)
full_data[, jit_outlier := as.integer(is_outlier) + runif(.N, -0.03, 0.03)]

a = coef(mod)[1]
b = coef(mod)[2]

# Plot with the proxy
plot(full_data[, fork_ratio], full_data[, jit_outlier],
	pch = 19, col = "#050505AA",
	xlab = "Fork ratio", ylab = "P(outlier)",
	main = "Fitted logit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)

# Does not discriminate at all...

## Run model with taper height
mod = glm(is_outlier ~ 1 + taper_height, data = full_data, family = binomial)
summary(mod)

## Plot stuff
a = coef(mod)[1]
b = coef(mod)[2]

# Plot with the proxy
plot(full_data[, taper_height], full_data[, jit_outlier],
	pch = 19, col = "#050505AA",
	xlab = "Fork ratio", ylab = "P(outlier)",
	main = "Fitted logit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)


## Run model with fork proxy
mod = glm(is_outlier ~ 1 + fork_proxy, data = full_data, family = binomial)
summary(mod)

## Plot stuff
a = coef(mod)[1]
b = coef(mod)[2]

# Plot with the proxy
plot(full_data[, fork_proxy], full_data[, jit_outlier],
	pch = 19, col = "#050505AA",
	xlab = "Fork ratio", ylab = "P(outlier)",
	main = "Fitted logit curve",
	ylim = c(-0.05, 1.05))

curve(inv_logit(a + b * x),
	add = TRUE, col = "#CD212A", lwd = 2)

# OK, whatever I try, there is NO stuff that seems to distinguish outliers from
#	non-outliers. In a way, it is a good sign! Maybe these trees are only
#	exceptional trees!
