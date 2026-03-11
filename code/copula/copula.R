
#### Aim of prog: Generate data according to a Clayton's copula and then fit a Clayton's copula written in Stan

#### Clear space and load packages
rm(list = ls())
graphics.off()

options(max.print = 500)

library(data.table)
library(cmdstanr)
library(copula)

#### Create fake data
## Define Clayton copula with theta = 1.33 (tau = 0.4)
clayton = claytonCopula(param = 1.33, dim = 2)

## Simulate uniform pairs (i.e., in I², then I need to transform back to the real scales)
set.seed(1969 - 08 - 18) # Woodstock seed
N_indiv = 200
uv = as.data.table(rCopula(N_indiv, clayton))

dim(uv) # N_indiv x 2
setnames(uv, new = c("Fbft_unif", "Vtot_unif"))

## Transform back to real scale using the marginals
uv[, Fbft := qgamma(Fbft_unif, shape = 3, rate = 6)] # mean = 0.5
uv[, Vtot := qlnorm(Vtot_unif, meanlog = 2, sdlog = 0.5)] # mean = 0.5

## Compute ranks
uv[, Fbft_rank := rank(Fbft)/(.N + 1)]
uv[, Vtot_rank := rank(Vtot)/(.N + 1)]

## Plot the ranls and ckeck Kendall's tau
round(uv[, cor(Fbft, Vtot, method = "kendall")], 2) # Should be around 0.4, as theta = 2*tau/(1 - tau)

plot(uv[, Fbft_rank], uv[, Vtot_rank], xlab = "Fbft rank", ylab = "Total rank",
	axes = FALSE, pch = 19, cex = 0.65, col = "#A1A1A155")
axis(1)
axis(2)
abline(a = 0, b = 1, col = "#CD212A", lty = "dashed", lwd = 4)
# This is more or less my stuff with real trees...

#### Fit a stan model on these data
## Compile model
model = cmdstan_model("./copula.stan")

## Prepare data
stanData = list(
	N = N_indiv,
	Fbft = uv[, Fbft],
	Vtot = uv[, Vtot]
)

## Common variables
n_chains = 4

## Run model
fit = model$sample(data = stanData, chains = n_chains, parallel_chains = min(n_chains, 4),
	max_treedepth = 12) # Ok it seems to work, YAHOO!

#### Generate data

