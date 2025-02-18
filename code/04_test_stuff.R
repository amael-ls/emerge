
#### Aim of prog: Test some stuff from statistical rethinking, McElreath

library(data.table)
library(cmdstanr)
library(ellipse)
library(MASS)

#### Population parameters, p. 437
a = 3.5 # Mean intercept (i.e., mean waiting time)
b = -1 # Mean slope (i.e., mean improvement waiting time in the afternoon. Improvement because < 0)

sigma_a = 1 # Variance intercept
sigma_b = 0.5 # Variance change morning afternoon (improvement if < 0, worsen otherwise)
sigma_cafe = 0.5 # Variance within cafe

rho = -0.7 # Covariance between intercepts and slopes

#### Var-Cov matrix
sigma_mat = diag(c(sigma_a, sigma_b))
rho_mat = matrix(c(1, rho, rho, 1), nrow = 2, ncol = 2)

varCov_mat = sigma_mat %*% rho_mat %*% sigma_mat

#### Simulate data
## Simulate parameters
n_cafes = 20
set.seed(5)

params_dt = as.data.table(mvrnorm(n_cafes, c(a, b), varCov_mat))
setnames(params_dt, new = c("alpha_cafe", "beta_cafe"))

want_plot = FALSE
if (want_plot)
{
	plot(params_dt[, alpha_cafe], params_dt[, beta_cafe], pch = 19, xlab = "alpha", ylab = "beta", axes = FALSE)
	axis(1)
	axis(2, las = 1)
	for (l in c(0.1, 0.3, 0.5, 0.8, 0.99))
		lines(ellipse(x = varCov_mat, centre = c(a, b), level = l), col = "#22442233", lwd = 2)
}

## Simulate waiting times
set.seed(22)
n_visit = 10 # 5 in the morning, 5 in the afternoon

visit_dt = data.table(cafe = rep(1:n_cafes, each = n_visit),
	is_afternoon = rep(c(rep(FALSE, n_visit/2), rep(TRUE, n_visit/2)), n_cafes))
setkey(visit_dt, cafe)

for (i in 1:n_cafes)
	visit_dt[.(i), waiting_time := rnorm(n = n_visit,
		mean = params_dt[i, alpha_cafe] + params_dt[i, beta_cafe]*is_afternoon, sd = sigma_cafe)]

