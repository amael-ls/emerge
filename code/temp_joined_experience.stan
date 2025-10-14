data {
	int <lower = 1> N;

	// Data
	vector <lower = 0, upper = 1> [N] p; // Ratio of V_bft / V_tot
	vector <lower = 0> [N] total_volume_m3;
}

transformed data {
	vector [N] logit_p = logit(p);
	vector [N] log_V = log(total_volume_m3);
}

parameters {
	// Linear regression
	real mu_p; // Mean logit(p)
	real mu_v; // Mean log(V)

	// Variances (in the variance-covariance matrices)
	real <lower = 0> sigma_p; // Unexplained variance of logit(p)
	real <lower = 0> sigma_v; // Unexplained variance of log(V_tot)
	cholesky_factor_corr[2] L;
}

transformed parameters {
	vector [2] mu_vec = to_vector({mu_p, mu_v});
	vector <lower = 0> [2] sigma_vec = to_vector({sigma_p, sigma_v});
}

model {
	// Priors...
	// ... regressions
	target += normal_lpdf(mu_p | 0, 1);
	target += normal_lpdf(mu_v | 0, 2);

	// ... variance
	target +=  lkj_corr_cholesky_lpdf(L | 2); // It contains rho (non-diag)

	// Variance (residuals)
	target += exponential_lpdf(sigma_p | 1);
	target += exponential_lpdf(sigma_v | 1);

	// Likelihood
	for (i in 1:N)
	{
		target += multi_normal_cholesky_lpdf(to_vector({logit_p[i], log_V[i]}) | mu_vec,
			diag_pre_multiply(sigma_vec, L)); // Regression of p and total volume with correlation
	}
}

generated quantities {
	cov_matrix[2] Sigma;
	real rho;
	
	// This is to recover the variance-covariance matrix instead of having its Cholesky factor
	Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_vec, L));

	// This is to recover the correlation parameter
	rho = tcrossprod(L)[2, 1];
}
