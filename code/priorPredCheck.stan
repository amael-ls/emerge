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
}

transformed data {
	vector [N] hdn = sqrt(circumference_m) ./ height;
	vector [N] slenderness = height ./ circumference_m;
	vector [N] corrected_cyl_vol = height .* circumference_m^2 ./ (4*pi()*(1 - 1.3/height)^2);
}

generated quantities {
	vector[2] b0;
	vector[2] b1;
	vector[2] b2;
	
	vector[S] beta0;
	vector[S] beta2;

	vector[N] volume;

	// Simulate the variance
	real<lower = 0> sigma = inv_gamma_rng(1, 1);
	real<lower = 0> sigma_beta0 = inv_gamma_rng(1, 1);
	real<lower = 0> sigma_beta2 = inv_gamma_rng(1, 1);

	// Simulate the common group intercepts and hdn/slenderness slopes
	for (g in 1:2) // g = functional type, either broadleaf or conifer
	{
		b0[g] = normal_rng(0, 10);
		b1[g] = normal_rng(0, 10);
		b2[g] = normal_rng(0, 10);
	}

	// Simulate group-specific intercepts and slenderness slopes, and individual volume data
	// --- Broadleaves
	for (s in 1:n_sp_broad)
	{
		beta0[s] = normal_rng(b0[1], sigma_beta0);
		beta2[s] = normal_rng(b2[1], sigma_beta2);

		for (i in ind_start_broad[s]:ind_end_broad[s])
			volume[i] = normal_rng((beta0[s] + b1[1]*hdn[i] + beta2[s]*slenderness[i]) .* corrected_cyl_vol[i], sigma);

	}

	// --- Conifers
	for (s in (n_sp_broad + 1):S)
	{
		beta0[s] = normal_rng(b0[2], sigma_beta0);
		beta2[s] = normal_rng(b2[2], sigma_beta2);

		for (i in ind_start_conif[s - n_sp_broad]:ind_end_conif[s - n_sp_broad])
			volume[i] = normal_rng((beta0[s] + b1[2]*hdn[i] + beta2[s]*slenderness[i]) .* corrected_cyl_vol[i], sigma);
	}
}