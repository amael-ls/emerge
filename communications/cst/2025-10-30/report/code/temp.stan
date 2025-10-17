data {
	// Dimension
	int <lower = 1> N;

	// Explanatory variable
	row_vector [N] explanatory;

	// Observations
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

	real b_p;
	real b_v;

	// Variances (in the variance-covariance matrices)
	real<lower = 0> sigma_p; // Unexplained variance of logit(p)
	real<lower = 0> sigma_v; // Unexplained variance of log(V_tot)
}

model {
	// Priors...
	// ... regressions
	target += normal_lpdf(mu_p | 0, 1);
	target += normal_lpdf(mu_v | 0, 2);

	target += normal_lpdf(b_p | 0, 1);
	target += normal_lpdf(b_v | 0, 1);

	// Variance (residuals)
	target += exponential_lpdf(sigma_p | 1);
	target += exponential_lpdf(sigma_v | 1);

	// Likelihoods
	target += normal_lpdf(logit_p | mu_p + b_p*explanatory, sigma_p);
	target += normal_lpdf(log_V | mu_v + b_v*explanatory, sigma_v);
}

generated quantities {
	array [N] real pred_Vbole;

	for (i in 1:N)
		pred_Vbole[i] = inv_logit(normal_rng(mu_p + b_p*explanatory[i], sigma_p)) *
			exp(normal_rng(mu_v + b_v*explanatory[i], sigma_p));
}
