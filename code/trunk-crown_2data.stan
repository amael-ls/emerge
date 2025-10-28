data {
	// Dimensions
	int <lower = 1> N;

	// Explanatory
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] circumference_m_nfi;

	// Observations
	vector <lower = 0> [N] bole_volume_m3;
	vector <lower = 0> [N] bole_volume_m3_nfi;
	vector <lower = 0> [N] crown_volume_m3;
}

transformed data {
	vector [N] log_circumference_m = log(circumference_m);
	vector [N] log_circumference_min = log_circumference_m - min(log_circumference_m);
	vector [N] log_circumference_min_nfi = log(circumference_m_nfi) - min(log(circumference_m_nfi));
	vector [N] log_b = log(bole_volume_m3);
	vector [N] log_b_nfi = log(bole_volume_m3_nfi);
	vector [N] log_c = log(crown_volume_m3);
}

parameters {
	// Linear regression...
	// ... intercepts
	real alpha_0;
	real beta_0;

	// ... slopes
	real alpha_1;
	real beta_1;
	real <lower = 0> alpha_2;
	real beta_2;
	real beta_3;
	real beta_4;

	// Variances (in the variance-covariance matrices)
	real <lower = 0> sigma_b; // Unexplained variance of log(bole volume)
	real <lower = 0> sigma_b_nfi; // Unexplained variance of log(bole volume) for NFI
	real <lower = 0> sigma_c; // Unexplained variance of log(crown volume)
	cholesky_factor_corr[2] L;
}

model {
	// Priors...
	// ... intercepts regressions
	target += normal_lpdf(alpha_0 | 0, 2);
	target += normal_lpdf(beta_0 | 0, 2);

	// ... slopes regressions
	target += normal_lpdf(alpha_1 | 0, 2);
	target += normal_lpdf(beta_1 | 0, 2);
	target += exponential_lpdf(alpha_2 | 1);
	target += normal_lpdf(beta_2 | 0, 2);
	target += normal_lpdf(beta_3 | 0, 2);
	target += normal_lpdf(beta_4 | 0, 2);

	// ... variance
	target +=  lkj_corr_cholesky_lpdf(L | 2); // It contains rho (non-diag)

	// Variance (residuals)
	target += exponential_lpdf(sigma_b | 1);
	target += exponential_lpdf(sigma_c | 1);

	// Likelihood, regression of bole and crown volumes with correlation
	for (i in 1:N)
	{
		// print(mu_vec[i, :]);
		target += multi_normal_cholesky_lpdf(to_vector({log_b[i], log_c[i]}) |
			[alpha_0 + alpha_1*log_circumference_min[i].^alpha_2,
			beta_0 + beta_1*log_circumference_m[i] + beta_2*log_circumference_m[i].^2 +
			beta_3*log_circumference_m[i].^3 + beta_4*log_circumference_m[i].^4],
			diag_pre_multiply([sigma_b, sigma_c], L));
	}
	target += normal_lpdf(log_b_nfi | alpha_0 + alpha_1*log_circumference_min_nfi.^alpha_2, sigma_b_nfi);
}

generated quantities {
	matrix [2, N] pred_Vb_Cr_joined;
	vector [N] pred_total_volume;
	cov_matrix[2] Sigma;

	real rho = tcrossprod(L)[2, 1];

	// Prediction of the total volume
	{
		vector [2] joined_pred;

		for (i in 1:N)
		{
			Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply([sigma_b, sigma_c], L));
			joined_pred = multi_normal_rng(
				[alpha_0 + alpha_1*log_circumference_min[i].^alpha_2,
				beta_0 + beta_1*log_circumference_m[i] + beta_2*log_circumference_m[i].^2 +
				beta_3*log_circumference_m[i].^3 + beta_4*log_circumference_m[i].^4],
				Sigma);

			pred_total_volume[i] = exp(joined_pred[1]) + exp(joined_pred[2]);
			pred_Vb_Cr_joined[, i] = joined_pred;
		}
	}
}

