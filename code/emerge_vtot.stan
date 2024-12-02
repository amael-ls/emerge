data {
	// Dimensions and indices
	int N; // Number of individuals
	int S; // Number of species
	int<lower = 0, upper = S> n_sp_broad; // number of broadleaf species
	int<lower = S - n_sp_broad, upper = S - n_sp_broad> n_sp_conif; // number of conifer species
	array[n_sp_broad] int ind_start_broad; // Broadleaf species index start
	array[n_sp_conif] int ind_start_conif; // Conifer species index start
	array[n_sp_broad] int ind_end_broad; // Broadleaf species index end
	array[n_sp_conif] int ind_end_conif; // Conifer species index end

	// Predictors
	vector<lower = 0> [N] height;
	vector<lower = 0> [N] circumference_m;

	// Response variable
	vector [N] volume_m3;
}

transformed data {
	vector [N] hdn = sqrt(circumference_m) ./ height;
	vector [N] slenderness = height ./ circumference_m;
	vector [N] corrected_cyl_vol = height .* circumference_m^2 ./ (4*pi()*(1 - 1.3/height)^2);
}

parameters {
	// Fixed effects (population parameters) for broadleaf and conifer
	vector[2] b0;
	vector[2] b1;
	vector[2] b2;
	
	// Random effects (group parameters)
	vector[S] eta0;
	vector[S] eta2;
	// vector[S] beta0;
	// vector[S] beta2;

	// Variances
	real<lower = 0> sigma; // sd residuals
	real<lower = 0> sigma_beta0; // sd random effect beta0
	real<lower = 0> sigma_beta2; // sd random effect beta2
}

transformed parameters {
	vector[S] beta0;
	vector[S] beta2;
	
	beta0[1:n_sp_broad] = b0[1] + eta0[1:n_sp_broad]*sigma_beta0;
	beta0[(n_sp_broad + 1):S] = b0[2] + eta0[(n_sp_broad + 1):S]*sigma_beta0;
	
	beta2[1:n_sp_broad] = b2[1] + eta2[1:n_sp_broad]*sigma_beta2;
	beta2[(n_sp_broad + 1):S] = b2[2] + eta2[(n_sp_broad + 1):S]*sigma_beta2;
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

	// target += gamma_lpdf(sigma | 0.1, 0.1); // Uses shape and rate
	// target += gamma_lpdf(sigma_beta0 | 4^2/0.2, 4/0.2);
	// target += gamma_lpdf(sigma_beta2 | 0.1, 0.1);

	// Hierarchy
	// target += normal_lpdf(beta0[1:n_sp_broad] | b0[1], sigma_beta0);
	// target += normal_lpdf(beta0[(n_sp_broad + 1):S] | b0[2], sigma_beta0);
	// target += normal_lpdf(beta2[1:n_sp_broad] | b2[1], sigma_beta2);
	// target += normal_lpdf(beta2[(n_sp_broad + 1):S] | b2[2], sigma_beta2);
	target += normal_lpdf(eta0 | 0, 1);
	target += normal_lpdf(eta2 | 0, 1);

	// Likelihood broadleaves, i = species
	for (i in 1:n_sp_broad)
	{
		target += normal_lpdf(volume_m3[ind_start_broad[i]:ind_end_broad[i]] |
			(beta0[i] +
			b1[1]*hdn[ind_start_broad[i]:ind_end_broad[i]] +
			beta2[i]*slenderness[ind_start_broad[i]:ind_end_broad[i]]) .*
			corrected_cyl_vol[ind_start_broad[i]:ind_end_broad[i]],
			sigma);
	}

	// Likelihood conifers, i = species
	for (i in 1:n_sp_conif)
	{
		target += normal_lpdf(volume_m3[ind_start_conif[i]:ind_end_conif[i]] |
			(beta0[n_sp_broad + i] +
			b1[2]*hdn[ind_start_conif[i]:ind_end_conif[i]] +
			beta2[n_sp_broad + i]*slenderness[ind_start_conif[i]:ind_end_conif[i]]) .*
			corrected_cyl_vol[ind_start_conif[i]:ind_end_conif[i]],
			sigma);
	}
}
