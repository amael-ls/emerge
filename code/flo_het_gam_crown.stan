/*
	Model that uses the model of FLorence heterosked Gamma, joined with a crown model
*/

data {
	// Dimensions
	int <lower = 1> N_inra;
	int <lower = 1> N;

	// Explanatory
	vector <lower = 0> [N_inra] circumference_m_inra;
	vector <lower = 0> [N_inra] taper_height_flo_inra;
	vector <lower = 0> [N_inra] height_inra;

	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] taper_height_flo;
	vector <lower = 0> [N] height;

	// Observations
	vector <lower = 0> [N_inra] crown_volume_m3_inra;
	vector <lower = 0> [N_inra] bole_volume_m3_inra;
	vector <lower = 0> [N] bole_volume_m3;
}

transformed data {
	// Predictors bole volume
	vector [N_inra] p1_inra = circumference_m_inra; // Predictor 1
	vector [N_inra] p2_inra = sqrt(circumference_m_inra) ./ taper_height_flo_inra; // Predictor 2
	vector [N_inra] p3_inra = sqrt(taper_height_flo_inra) ./ (circumference_m_inra.^2 .* height_inra); // Pred. 3
	vector [N_inra] p4_inra = 1 - taper_height_flo_inra ./ height_inra; // Predictor 4

	vector [N] p1 = circumference_m; // Predictor 1
	vector [N] p2 = sqrt(circumference_m) ./ taper_height_flo; // Predictor 2
	vector [N] p3 = sqrt(taper_height_flo) ./ (circumference_m.^2 .* height); // Predictor 3
	vector [N] p4 = 1 - taper_height_flo ./ height; // Predictor 4

	// Predictor crown volume
	vector [N_inra] log_circumference_m_inra = log(circumference_m_inra);
	vector [N] log_circumference_m = log(circumference_m);

	// Transformed dependent variables...
	// ... Form factor
	vector [N] fnewbft = 4*pi()*bole_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
	vector [N_inra] fnewbft_inra = 4*pi()*bole_volume_m3_inra ./ (circumference_m_inra.^2 .* height_inra) .*
		(1 - 1.3 ./ height_inra).^2;
	
	// ... Log crown volume
	vector [N_inra] log_crown = log(crown_volume_m3_inra);

}

parameters {
	/*
		alphas are for bole volumes
		betas are for crown volumes
	*/
	// Linear regression...
	// ... intercepts
	real alpha0;
	real beta0;

	// ... slopes
	real alpha1;
	real alpha2;
	real alpha3;
	real alpha4;
	
	real beta1;
	real beta2;
	real beta3;
	real beta4;

	// Variances (in the variance-covariance matrices)
	real <lower = 0> sigma_b; // Unexplained variance of log(bole volume)
	real <lower = 0> sigma_b_inra; // Unexplained variance of log(bole volume) for INRA
	real <lower = 0> sigma_c; // Unexplained variance of log(crown volume)
	cholesky_factor_corr[2] L;
}

model {
	// Priors...
	// ... intercepts regressions
	target += normal_lpdf(alpha0 | 0, 1);
	target += normal_lpdf(beta0 | 0, 2);

	// ... slopes regressions
	target += normal_lpdf(alpha1 | 0, 0.1);
	target += normal_lpdf(alpha2 | 0, 1);
	target += normal_lpdf(alpha3 | 0, 0.5);
	target += normal_lpdf(alpha4 | 0, 0.5);

	target += normal_lpdf(beta1 | 0, 2);
	target += normal_lpdf(beta2 | 0, 2);
	target += normal_lpdf(beta3 | 0, 2);
	target += normal_lpdf(beta4 | 0, 2);

	// ... variance
	target +=  lkj_corr_cholesky_lpdf(L | 2); // It contains rho (non-diag)

	// Variance (residuals)
	target += exponential_lpdf(sigma_b | 1);
	target += exponential_lpdf(sigma_b_inra | 1);
	target += exponential_lpdf(sigma_c | 1);

	// Likelihood...
	// ... on bole volume, regression on the NFI data
	target += normal_lpdf(fnewbft | alpha0 + alpha1*p1 + alpha2*p2 + alpha3*p3 + alpha4*p4, sigma_b);
	
	// ... on bole and crown volumes, regression with correlation on Emerge data (INRA)
	for (i in 1:N_inra)
	{
		// print(mu_vec[i, :]);
		target += multi_normal_cholesky_lpdf(to_vector({fnewbft_inra[i], log_crown[i]}) |
			[alpha0 + alpha1*p1_inra[i] + alpha2*p2_inra[i] + alpha3*p3_inra[i] + alpha4*p4_inra[i],
			beta0 + beta1*log_circumference_m_inra[i] + beta2*log_circumference_m_inra[i].^2 +
			beta3*log_circumference_m_inra[i].^3 + beta4*log_circumference_m_inra[i].^4],
			diag_pre_multiply([sigma_b_inra, sigma_c], L));
	}
}

generated quantities {
	matrix [2, N] pred_Vb_Cr_joined;
	matrix [2, N] pred_Vb_Cr_joined_inra;
	cov_matrix[2] Sigma;
	cov_matrix[2] Sigma_inra;

	real rho = tcrossprod(L)[2, 1];

	// Prediction of the total volume
	{
		vector [2] joined_pred;

		for (i in 1:N)
		{
			Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply([sigma_b, sigma_c], L));
			joined_pred = multi_normal_rng(
				[alpha0 + alpha1*p1[i] + alpha2*p2[i] + alpha3*p3[i] + alpha4*p4[i],
				beta0 + beta1*log_circumference_m[i] + beta2*log_circumference_m[i].^2 +
				beta3*log_circumference_m[i].^3 + beta4*log_circumference_m[i].^4],
				Sigma);

			pred_Vb_Cr_joined[, i] = joined_pred;
		}

		for (i in 1:N_inra)
		{
			Sigma_inra = multiply_lower_tri_self_transpose(diag_pre_multiply([sigma_b_inra, sigma_c], L));
			joined_pred = multi_normal_rng(
				[alpha0 + alpha1*p1_inra[i] + alpha2*p2_inra[i] + alpha3*p3_inra[i] + alpha4*p4_inra[i],
				beta0 + beta1*log_circumference_m_inra[i] + beta2*log_circumference_m_inra[i].^2 +
				beta3*log_circumference_m_inra[i].^3 + beta4*log_circumference_m_inra[i].^4],
				Sigma_inra);

			pred_Vb_Cr_joined_inra[, i] = joined_pred;
		}
	}
}

