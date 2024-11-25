
#### Aim of prog: Estimate the parameters for Vol tot Emerge
## Comments
#	1. The values of the estimated parameters by C. Deleuze can be found in inventR::ListTarEmerge
#	2. For the columns, useless to have v_br1 ... 6 as v_br - (v_br1 + v_br2 + v_br3 + v_br4 + v_br5 + v_br6) always equals 0
#		(provided NA are replaced by 0)
#
## Reference
#	Deleuze, C., Morneau, F., [...], Vallet, P., 2014.
#	Estimer le volume total d’un arbre, quelles que soient l’essence, la taille, la sylviculture, la station.
#	Rendez-vous techniques de l’ONF 44, 22–32.
#

rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(MetBrewer)
library(cmdstanr)
library(stringi)
library(nlme)

#### Prepare data
## Loading
load("../data/EMERGE.RData")
rm(list = setdiff(ls(), "inra_arbres"))
setDT(inra_arbres)

## Keep only column of interest and rename them
inra_arbres[, tot_vol := v_tronc_verif + v_fourche_verif +v_fourche2_verif +v_br_verif +v_menu_verif]
inra_arbres = unique(inra_arbres[, .(nom_fichier, essence, c130, h_tot, tot_vol, genre)])
inra_arbres = na.omit(inra_arbres)

setnames(inra_arbres, new = c("unique_id", "speciesName_sci", "circumference_cm", "height", "tot_vol", "genus"))

inra_arbres[, is_broadleaf := TRUE]
inra_arbres[genus %in% c("Abies", "Cedrus", "Larix", "Picea", "Pinus", "Pseudotsuga", "Thuya", "Tsuga"),
	is_broadleaf := FALSE]

## Compute total volume
setkey(inra_arbres, is_broadleaf, speciesName_sci, unique_id)
inra_arbres[, .N, by = speciesName_sci]

## Indices
# Conifer and broadleaf
lim_broadleaf = range(which(inra_arbres[, is_broadleaf]))
if (lim_broadleaf[1] == 1)
	lim_broadleaf = lim_broadleaf[2] + 1 # Last broadleaf + 1 = first conifer

if ((length(lim_broadleaf) == 2) && (lim_broadleaf[2] == inra_arbres[, .N])) # First broadleaf
	lim_broadleaf = lim_broadleaf[1]

if (length(lim_broadleaf) != 1)
	stop("Error in lim_broadleaf")

# Find start and end indices for each species
ind_species = inra_arbres[, .(start = .I[1], end = .I[.N]), by = .(speciesName_sci)]
ind_species = merge.data.table(ind_species, unique(inra_arbres[, .(speciesName_sci, is_broadleaf)]), by = "speciesName_sci")
setorder(ind_species, start)

## Stan data
stanData = list(
	N = inra_arbres[, .N],
	S = ind_species[, .N],
	ind_start_sp = ind_species[, start],
	ind_end_sp = ind_species[, end],
	lim_broadleaf = lim_broadleaf,
	height = inra_arbres[, height],
	circumference_m = inra_arbres[, circumference_cm/100],
	volume_m3 = inra_arbres[, tot_vol/1000]
)

#### Tool functions
source("./dummy/toolFunctions.R")

## Function to compute tree volume according to model
vol_fct = function(params_vec, predictor_mat, corrected_cyl_vol)
	return((predictor_mat %*% params_vec) * corrected_cyl_vol)

#### Run model
# ## lm
# inra_arbres[, factForm := tot_vol/(1000*corrected_cyl_vol)]
inra_arbres[, circumference_m := circumference_cm/100]
inra_arbres[, hdn := sqrt(circumference_m)/height]
inra_arbres[, slenderness := height/circumference_m]

inra_arbres[, cor(hdn, slenderness)]
# plot(inra_arbres[, hdn], inra_arbres[, slenderness], pch = 19, col = "#FFAF3766")

## Common variables
n_chains = 4
woodstock_seed = 1969 - 08 - 18

## Compile
model = cmdstan_model("./emerge_vtot.stan")

## Fit
if (!file.exists("fit_test.rds"))
{
	fit = model$sample(data = stanData, chains = n_chains, parallel_chains = ifelse(n_chains < 4, n_chains, 4),
		seed = woodstock_seed, refresh = 200, max_treedepth = 12, save_warmup = TRUE)
	fit$save_output_files(dir = "./", basename = paste0("groupEffect_gamma"), random = FALSE)
	saveRDS(fit, "fit_test_groupEffect_gamma.rds")
} else {
	fit = readRDS("./fit_test.rds")
}

fit$cmdstan_diagnose()
fit$summary()

#### Explore results
## Check draws sd
lazyTrace(fit$draws("sigma"))
lazyTrace(fit$draws("b0[1]", inc_warmup = TRUE))

## Check prediction (cmdstanr, Deleuze) vs measured
params = as.data.table(fit$summary())
setkey(params, variable)

inra_arbres[, c("intercept", "hdn", "slenderness") := .(rep(1, .N), sqrt(circumference_cm/100)/height, 100*height/circumference_cm)]
inra_arbres[, corrected_cyl_vol := height*circumference_cm^2/(4e4*pi*(1 - 1.3/height)^2)]

pred_vol_cmdstanr = vol_fct(params_vec = params[c("beta0", "beta1", "beta2"), mean],
	predictor_mat = as.matrix(inra_arbres[, .(intercept, hdn, slenderness)]),
	corrected_cyl_vol = inra_arbres[, corrected_cyl_vol])

params_deleuze_fag = c(beta0 = 0.542, beta1 = 0.661, beta2 = -0.002)

# params_deleuze_abies_alba = c(beta0 = 0.398, beta1 = 1.756, beta2 = 0.002)

pred_vol_deleuze = vol_fct(params_vec = params_deleuze_fag,
	predictor_mat = as.matrix(inra_arbres[, .(intercept, hdn, slenderness)]),
	corrected_cyl_vol = inra_arbres[, corrected_cyl_vol])

pdf("test_deleuze_vs_stan-Fag-syl.pdf")
par(xpd = TRUE)
plot(pred_vol_cmdstanr, inra_arbres[, tot_vol/1000], pch = 20, cex = 0.9, col = "#FFAF37",
	axes = FALSE, xlab = "Predicted volume (m³)", ylab = "Measured volume (m³)", ylim = c(0, max(pred_vol_cmdstanr, pred_vol_deleuze)))
points(pred_vol_deleuze, inra_arbres[, tot_vol/1000], pch = 20, cex = 0.9, col = "#007BA544")
# points(toto*inra_arbres[, corrected_cyl_vol], inra_arbres[, tot_vol/1000], pch = 20, cex = 0.9)
abline(a = 0, b = 1, col = "#F24000", lwd = 2)
axis(1)
axis(2, las = 1)
legend(x = "topleft", inset = c(0, -0.1), legend = c("CmdStanR", "Deleuze"), pch = 15, col = c("#FFAF37", "#007BA5"))
dev.off()


## Residuals
res_cmdstanr = pred_vol_cmdstanr - inra_arbres[, tot_vol/1000]
res_deleuze = pred_vol_deleuze - inra_arbres[, tot_vol/1000]

par(xpd = TRUE)
plot(res_cmdstanr, pch = 20, cex = 0.9, col = "#FFAF37", xlab = "", ylab = "",
	axes = FALSE, ylim = c(min(res_cmdstanr, res_deleuze), max(res_cmdstanr, res_deleuze)))
points(res_cmdstanr, pch = 20, cex = 0.9, col = "#007BA5")
abline(h = 0, col = "#F24000", lwd = 2)
legend(x = "topleft", inset = c(0, -0.1), legend = c("CmdStanR", "Deleuze"), pch = 15, col = c("#FFAF37", "#007BA5"))

## Ratio Deleuze vs cmdstanr
params_deleuze_fag/params[c("beta0", "beta1", "beta2"), mean]
round(params[c("beta0", "beta1", "beta2"), mean], 3)
params_deleuze_fag

cor(fit$draws("beta0"), fit$draws("beta1"))
cor(fit$draws("beta0"), fit$draws("beta2"))
cor(fit$draws("beta1"), fit$draws("beta2"))

lazyPosterior(draws = fit$draws("beta0"), val1 = params_deleuze_fag["beta0"])
lazyPosterior(draws = fit$draws("beta1"), val1 = params_deleuze_fag["beta1"])
lazyPosterior(draws = fit$draws("beta2"), val1 = params_deleuze_fag["beta2"])

#### Heatmap of the posterior marginal distribution of beta1 and beta2
extract_mat = function(fit, ls_params = c("beta1", "beta2"))
{
	n_chains = fit$num_chains()
	n_sampling = fit$metadata()$iter_sampling
	params = matrix(data = NA, ncol = length(ls_params), n_chains*n_sampling)
	draws = fit$draws(ls_params)
	for (i in seq_len(n_sampling))
		for (j in seq_len(ncol(params)))
			params[((i - 1)*n_chains + 1):(i*n_chains), j] = c(draws[i, , j])

	return(params)
}

params_mat = extract_mat(fit)
dim(params_mat)

kern = MASS::kde2d(x = params_mat[, 1], y = params_mat[, 2], n = 100)

image(kern, col = met.brewer(name = "Hokusai3", n = 100, type = "continuous"))
contour(kern, add = TRUE, lwd = 2, col = "#CD212A", labcex = 1.25)

#### Model that was used by C. Deleuze
## In the original model, there are the variables formTot, and formTotNew. I need to rebuild them
inra_arbres[, formTot := 4*pi*tot_vol/(height*circumference_m^2)]
inra_arbres[, formTotNew := formTot * (1 - 1.3/height)^2]
inra_arbres[, feuil.res := "feuillus"]
inra_arbres[!(is_broadleaf), feuil.res := "conifer"]
inra_arbres[, feuil.res := as.factor(feuil.res)]

# Original model, I just change the names of the variables in what follows
# Modele_mai2 = nlme(formTotNew ~ a + b*hdn + d*hsurd, data = Grdata.PV, start = c(a=0.4, 0, b = 1.5, 0, d = 0.0005, 0),
# 	fixed = list(a + b + d ~ feuil.res), random = a + d ~ 1|nomessence2)

Modele_mai2 = nlme(formTotNew ~ a + b*hdn + d*slenderness, data = inra_arbres,
	start = c(a = 0.4, 0, b = 1.5, 0, d = 0.0005, 0),
	fixed = list(a + b + d ~ feuil.res), random = a + d ~ 1|speciesName_sci)

test_lme4 = lme4::lmer(formTotNew ~ feuil.res * (hdn + slenderness) + (1 | speciesName_sci), data = inra_arbres)

params_deleuze_fag = c(beta0 = 0.542, beta1 = 0.661, beta2 = -0.002)


#### Is my model correctly set-up?
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(stringi)

source("./dummy/toolFunctions.R")

set.seed(1969 - 08 - 18) # Woodstock seed

b0 = c(conif = 0.2, broad = 5.3)
b1 = c(conif = 0.27, broad = 1.4)
b2 = c(conif = 4, broad = 0.21)

sigma_beta0 = 4.2
sigma_beta2 = 6.98
sigma = 1.2

n = 3e3
S = 30 # Number of species

generate_random_partition = function(n, k)
{
	if (n < k)
		stop("n should be larger than k")
	parts = c(0, sort(sample(1:(n - 1), k - 1)), n)
	return(diff(parts))
}

rep_species = generate_random_partition(n, S)

lim_broadleaf = sum(rep_species[1:15]) + 1

n_conif = lim_broadleaf - 1
n_broad = n - n_conif

fake_dt = data.table(
	species = rep(paste0("sp_", 1:S), times = rep_species),
	type = c(rep("conif", n_conif), rep("broad", n_broad)),
	fake_hdn = runif(n = n, min = -10, max = 50),
	fake_slenderness = rnorm(n = n, mean = 0, sd = 20),
	b0 = c(rep(b0["conif"], n_conif), rep(b0["broad"], n_broad)),
	b1 = c(rep(b1["conif"], n_conif), rep(b1["broad"], n_broad)),
	b2 = c(rep(b2["conif"], n_conif), rep(b2["broad"], n_broad)))

fake_dt[, beta0 := rnorm(1, unique(b0), sigma_beta0), by = species]
fake_dt[, beta2 := rnorm(1, unique(b2), sigma_beta2), by = species]

# fake_dt[, table(round(beta0, 2))]
# fake_dt[, table(round(beta2, 2))]

# fake_dt[, table(b0)]
# fake_dt[, table(b1)]
# fake_dt[, table(b2)]

# fake_dt[, fake_mu := beta0 + b1*fake_hdn + beta2*fake_slenderness]
fake_dt[, fake_mu := b0 + b1*fake_hdn + beta2*fake_slenderness]
fake_dt[, cor(fake_hdn, fake_slenderness)]
# fake_dt[, fake_vol := rnorm(.N, fake_mu, sigma)]

for (i in 1:n)
	fake_dt[i, fake_vol := rnorm(1, fake_mu, sigma)]

fake_dt[, sd(fake_mu), by = species]

ind_species = fake_dt[, .(start = .I[1], end = .I[.N]), by = .(species)]
ind_species[, n_indiv := end - start + 1, by = species]
ind_species[, sum(n_indiv)] == n

stanData = list(
	N = fake_dt[, .N],
	S = S,
	ind_start_sp = ind_species[, start],
	ind_end_sp = ind_species[, end],
	lim_broadleaf = lim_broadleaf,
	fake_hdn = fake_dt[, fake_hdn],
	fake_slenderness = fake_dt[, fake_slenderness],
	volume_m3 = fake_dt[, fake_vol]
)

# plot(1:n_conif, fake_dt[type == "conif", fake_vol], pch = 19,
# 	col = "#FFAF3766", xlim = c(0, n), xlab = "", ylab = "")
# points((n_conif + 1):n, fake_dt[type == "broad", fake_vol], pch = 19, col = "#AA228866")
# abline(v = c(ind_species[, start], ind_species[.N, end]), lwd = 2, col = "#CD212A")

# fake_dt[, colour := ifelse(type == "conif", "#FFAF3799", "#AA228808")]
# plot(fake_dt[, fake_mu], fake_dt[, fake_vol], pch = 19,
# 	col = fake_dt[, colour], xlab = "", ylab = "")

tt = lm(fake_dt[, fake_vol] ~ 0 + fake_dt[, fake_mu])
summary(tt)$sigma # Should be close to sigma
summary(tt)$coefficients[, "Estimate"] # Should be close to 1

## nlme test
# nlme_test = nlme(formTotNew ~ a + b*hdn + d*slenderness, data = inra_arbres,
# 	start = c(a = 0.4, 0, b = 1.5, 0, d = 0.0005, 0),
# 	fixed = list(a + b + d ~ feuil.res), random = a + d ~ 1|speciesName_sci)

## Common variables
n_chains = 4

## Compile
# model = cmdstan_model("./emerge_vtot.stan")
model = cmdstan_model("./emerge_vtot.stan")

## Fit
if (!file.exists("fit_test.rds"))
{
	fit = model$sample(data = stanData, chains = n_chains, parallel_chains = ifelse(n_chains < 4, n_chains, 4),
		refresh = 200, max_treedepth = 12, save_warmup = TRUE, iter_warmup = 500, iter_sampling = 1500, adapt_delta = 0.8)
	fit$save_output_files(dir = "./", basename = paste0("groupEffect_gamma"), random = FALSE)
	saveRDS(fit, "fit_test_groupEffect_gamma.rds")
} else {
	fit = readRDS("./fit_test.rds")
}

fit$cmdstan_diagnose()
fit$summary()

lazyTrace(fit$draws("sigma", inc_warmup = FALSE))
lazyTrace(fit$draws("b0[1]", inc_warmup = TRUE))
lazyTrace(fit$draws("b0[2]", inc_warmup = TRUE))

lazyTrace(fit$draws("b1[1]"))

beta0_draws = fit$draws("beta0")
round(apply(X = beta0_draws, MARGIN = 3, FUN = mean), 2)
fake_dt[, unique(round(beta0, 2)), by = species]

beta2_draws = fit$draws("beta2")
round(apply(X = beta2_draws, MARGIN = 3, FUN = mean), 2)

ind_species[, beta2_estimate := round(apply(X = beta2_draws, MARGIN = 3, FUN = mean), 2)]
ind_species[, beta2_real := fake_dt[, unique(round(beta2, 2)), by = species][, V1]]

b2_draws = fit$draws("b2")
round(apply(X = b2_draws, MARGIN = 3, FUN = mean), 2)
b2