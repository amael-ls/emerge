data {
	// Dimensions and indices
	int N; // Number of individuals
	int S; // Number of species
	int<lower = 0, upper = S> n_sp_conif; // number of conifer species
	int<lower = S - n_sp_conif, upper = S - n_sp_conif> n_sp_broad; // number of broadleaf species
	array[n_sp_conif] int ind_start_conif; // Conifer species index start
	array[n_sp_broad] int ind_start_broad; // Broadleaf species index start
	array[n_sp_conif] int ind_end_conif; // Conifer species index end
	array[n_sp_broad] int ind_end_broad; // Broadleaf species index end

	// Predictors
	vector [N] fake_hdn;
	vector [N] fake_slenderness;

	// Response variable
	vector[N] volume_m3;
}

parameters {
	// Fixed effects (population parameters)
	vector[2] b0;
	vector[2] b1;
	vector[2] b2;
	
	// Random effects (group parameters)
	vector[S] beta0;
	vector[S] beta2;

	// Variances
	real<lower = 0> sigma; // sd residuals
	real<lower = 0> sigma_beta0; // sd random effect beta0
	real<lower = 0> sigma_beta2; // sd random effect beta2
}

model {
	// Priors
	// --- Population parameters
	target += normal_lpdf(b0 | 0, 10);
	target += normal_lpdf(b1 | 0, 10);
	target += normal_lpdf(b2 | 0, 10);

	// --- Residual variance and population variance
	target += inv_gamma_lpdf(sigma | 1, 1); // Uses shape and scale
	target += inv_gamma_lpdf(sigma_beta0 | 1, 1);
	target += inv_gamma_lpdf(sigma_beta2 | 1, 1);

	// Hierarchy
	target += normal_lpdf(beta0[1:n_sp_conif] | b0[1], sigma_beta0);
	target += normal_lpdf(beta2[1:n_sp_conif] | b2[1], sigma_beta2);
	target += normal_lpdf(beta0[(n_sp_conif + 1):S] | b0[2], sigma_beta0);
	target += normal_lpdf(beta2[(n_sp_conif + 1):S] | b2[2], sigma_beta2);

	for (i in 1:n_sp_conif)
	{	
		// Likelihood conifers
		target += normal_lpdf(volume_m3[ind_start_conif[i]:ind_end_conif[i]] | beta0[i] +
			b1[1]*fake_hdn[ind_start_conif[i]:ind_end_conif[i]] +
			beta2[i]*fake_slenderness[ind_start_conif[i]:ind_end_conif[i]], sigma);
	}

	for (i in 1:n_sp_broad)
	{
		// Likelihood broadleaves
		target += normal_lpdf(volume_m3[ind_start_broad[i]:ind_end_broad[i]] | beta0[n_sp_conif + i] +
			b1[2]*fake_hdn[ind_start_broad[i]:ind_end_broad[i]] +
			beta2[n_sp_conif + i]*fake_slenderness[ind_start_broad[i]:ind_end_broad[i]], sigma);
	}
}


