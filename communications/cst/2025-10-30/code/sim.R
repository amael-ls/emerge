
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

fit = model$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4), max_treedepth = 12)

## Check posteriors
# lazyPosterior(fit$draws("mu_p"), val1 = mu_logit_r)
# lazyPosterior(fit$draws("mu_v"), val1 = mu_log_v)
# lazyPosterior(fit$draws("b_p"), val1 = b_logit_r)
# lazyPosterior(fit$draws("b_v"), val1 = b_log_v)
# lazyPosterior(fit$draws("rho"), val1 = rho)

# lazyPosterior(fit$draws("sigma_p"), val1 = sigma_logit_r)
# lazyPosterior(fit$draws("sigma_v"), val1 = sigma_log_v)

# varCov_mat_est = matrix(data = apply(X = fit$draws("Sigma"), FUN = mean, MARGIN = 3), nrow = 2, byrow = 2)

# round(varCov_mat, 2)
# round(varCov_mat_est, 2)

## Fit the model assuming uncorrelated random variables
model2 = cmdstan_model("./temp.stan")
fit2 = model2$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4), max_treedepth = 12)

# lazyPosterior(fit2$draws("mu_p"), val1 = mu_logit_r)
# lazyPosterior(fit2$draws("mu_v"), val1 = mu_log_v)
# lazyPosterior(fit2$draws("b_p"), val1 = b_logit_r)
# lazyPosterior(fit2$draws("b_v"), val1 = b_log_v)

# lazyPosterior(fit2$draws("sigma_p"), val1 = sigma_logit_r)
# lazyPosterior(fit2$draws("sigma_v"), val1 = sigma_log_v)

#### Plot posterior product
pred_1 = fit$draws("pred_Vbole")
quant_1 = apply(X = pred_1, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))
pred_2 = fit2$draws("pred_Vbole")
quant_2 = apply(X = pred_2, MARGIN = 3, FUN = quantile, probs = c(0.05, 0.95))

k = 38
if (k > data[, .N])
	stop("k must be smaller than N")

lazyPosterior(pred_1[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")])
lazyPosterior(pred_2[, , paste0("pred_Vbole[", k, "]")], val1 = data[k, Vbole],
	max_x = quant_1["95%", paste0("pred_Vbole[", k, "]")])

data[k, Vbole]
