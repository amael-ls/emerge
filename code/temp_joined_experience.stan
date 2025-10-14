data {
	int <lower = 1> N;

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;

	// Observations
	vector <lower = 0, upper = 1> [N] p; // Ratio of V_bft / V_tot
	vector <lower = 0> [N] total_volume_m3;
}

transformed data {
	vector [N] log_circumference_m = log(circumference_m);

	vector [N] logit_p = logit(p);
	vector [N] log_V = log(total_volume_m3);
}

parameters {
	// Linear regression...
	// ... intercepts
	real alpha_0; // logit(p)
	real beta_0; // log(V)

	// ... slopes
	real alpha_1;
	real alpha_2;
	real beta_1;

	// Variances (in the variance-covariance matrices)
	real <lower = 0> sigma_p; // Unexplained variance of logit(p)
	real <lower = 0> sigma_v; // Unexplained variance of log(V_tot)
	cholesky_factor_corr[2] L;
}

transformed parameters {
	vector <lower = 0> [2] sigma_vec = to_vector({sigma_p, sigma_v});
}

model {
	// Priors...
	// ... intercepts regressions
	target += normal_lpdf(alpha_0 | 0, 2);
	target += normal_lpdf(beta_0 | 0, 2);

	// ... slopes regressions
	target += normal_lpdf(alpha_1 | 0, 2);
	target += normal_lpdf(alpha_2 | 0, 2);
	target += normal_lpdf(beta_1 | 0, 2);

	// ... variance
	target +=  lkj_corr_cholesky_lpdf(L | 2); // It contains rho (non-diag)

	// Variance (residuals)
	target += exponential_lpdf(sigma_p | 1);
	target += exponential_lpdf(sigma_v | 1);

	// Likelihood, regression of p and total volume with correlation
	for (i in 1:N)
	{
		// print(mu_vec[i, :]);
		target += multi_normal_cholesky_lpdf(to_vector({logit_p[i], log_V[i]}) |
			[alpha_0 + alpha_1*log_circumference_m[i] + alpha_2*log_circumference_m[i]^2,
				beta_0 + beta_1*log_circumference_m[i]], diag_pre_multiply(sigma_vec, L));
	}
}

generated quantities {
	vector [N] pred_Vbole;
	cov_matrix[2] Sigma;
	real rho;
	
	// This is to recover the variance-covariance matrix instead of having its Cholesky factor
	Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_vec, L));

	// This is to recover the correlation parameter
	rho = tcrossprod(L)[2, 1];

	// Prediction of the bole volume
	{
		vector [2] joined_pred;

		for (i in 1:N)
		{
			joined_pred = multi_normal_rng(
				[alpha_0 + alpha_1*log_circumference_m[i] + alpha_2*log_circumference_m[i]^2,
				beta_0 + beta_1*log_circumference_m[i]], Sigma);

			pred_Vbole[i] = inv_logit(joined_pred[1])*exp(joined_pred[2]);
		}
	}
}
