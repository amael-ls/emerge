
#### Clear space and load packages
# rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(MetBrewer)
library(cmdstanr)
library(stringi)
library(DHARMa)

## Working directory
wd = "/home/ALe-Squin/work/FairCarbon/emerge/communications/cst/2025-10-30/code"
if (getwd() != wd)
	setwd(wd)

## Tool function
source("./toolFunctions.R")

inv_logit = function(x)
	return (1 / (1 + exp(-x)))

#### Parameters
# Intercept
mu_logit_r = 1.386 # Mean logit(ratio)
mu_log_v = -1.194 # Mean total volume

# Slopes
b_logit_r = 0.87
b_log_v = 1.21

sigma_logit_r = 0.5625 # Std. dev. logit(ratio)
sigma_log_v = 1.25 # Std. dev. total volume. Not ideal but larger makes it too fat-tailed compared to real distrib

rho = 0.22 # Covariance between logit(ratio) and log(total volume)

#### Var-Cov matrix
sigma_mat = diag(c(sigma_logit_r, sigma_log_v))
rho_mat = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)

varCov_mat = sigma_mat %*% rho_mat %*% sigma_mat

#### Other variables
n_sim = 2e3
n_chains = 4

#### Create data following our model (intercepts only, with correlation in the residuals)
set.seed(1969 - 08 - 18) # Woodstock seed

data = data.table(explanatory = runif(n = n_sim, min = -2, max = 5), logit_r = rep(-Inf, n_sim),
	log_v = rep(-Inf, n_sim))

for (i in seq_len(n_sim))
	data[i, c("logit_r", "log_v") := as.list(MASS::mvrnorm(1, c(mu_logit_r + b_logit_r*explanatory,
		mu_log_v + b_log_v*explanatory), varCov_mat))]

data[, Vbole := inv_logit(logit_r)*exp(log_v)]

#### Fit the model
## Stan data
stanData = list(
	# Dimension
	N = data[, .N],

	# Explanatory variable
	explanatory = data[, explanatory],

	# Observations
	p = data[, inv_logit(logit_r)],
	total_volume_m3 = data[, exp(log_v)]
)

## Compile model
model = cmdstan_model("./temp_joined.stan")

## Run
fit = model$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4), max_treedepth = 12)

#### Fit a second model assuming uncorrelated random variables
model2 = cmdstan_model("./temp.stan")
fit2 = model2$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4), max_treedepth = 12)

#### Check posteriors
## Posterior of intercept parameters
lazyPosterior(fit$draws("mu_p"), val1 = mu_logit_r, filename = "../Figures/alpha_corr")
lazyPosterior(fit2$draws("mu_p"), val1 = mu_logit_r, filename = "../Figures/alpha_uncorr")

lazyPosterior(fit$draws("sigma_v"), val1 = sigma_log_v, filename = "../Figures/sigma_2_corr")
lazyPosterior(fit2$draws("sigma_v"), val1 = sigma_log_v, filename = "../Figures/sigma_2_uncorr")

#### Plot posterior product
pred_1 = fit$draws("pred_Vbole")
quant_1 = apply(X = pred_1, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))
pred_2 = fit2$draws("pred_Vbole")
quant_2 = apply(X = pred_2, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))

## For a 'small' value of Vbole (around 15th percentile)
quantile(data[, Vbole], c(0.15, 0.5, 0.85))
which(data[, Vbole > 0.05] & data[, Vbole < 0.055])
k = 42

lazyPosterior(pred_1[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/15th-percentile_corr")
lazyPosterior(pred_2[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/15th-percentile_uncorr")

## For a 'medium' value of Vbole (around 50th percentile)
which(data[, Vbole > 1.78] & data[, Vbole < 1.79])
k = 177

lazyPosterior(pred_1[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/50th-percentile_corr")
lazyPosterior(pred_2[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/50th-percentile_uncorr")

## For a 'large' value of Vbole (around 85 percentile)
which(data[, Vbole > 39] & data[, Vbole < 40])
k = 695

lazyPosterior(pred_1[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/85th-percentile_corr")
lazyPosterior(pred_2[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")], n = 4096,
	filename = "../Figures/85th-percentile_uncorr")
