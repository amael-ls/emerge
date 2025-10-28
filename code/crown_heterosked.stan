data {
	// Dimensions
	int <lower = 1> N;
	int <lower = 1> N_group;
	array [N_group] int <lower = 1, upper = N> start_ind; // Index group starts
	array [N_group] int <lower = 1, upper = N> end_ind; // Index group ends

	// Explanatory
	vector <lower = 0> [N] circumference_m;

	// Observations
	vector <lower = 0> [N] crown_volume_m3;
}

transformed data {
	vector [N] log_circumference_m = log(circumference_m);
	vector [N] log_c = log(crown_volume_m3);
}

parameters {
	// Linear regression
	vector [N_group] alphas;
	vector [N_group] betas;

	// Variances
	vector <lower = 0> [N_group] sigmas;
}

model {
	// Priors
	target += normal_lpdf(alphas | 0, 2);
	target += normal_lpdf(betas | 0, 2);
	target += gamma_lpdf(sigmas | 0.1, 0.1);

	// Likelihood
	for (g in 1:N_group)
	{
		target += normal_lpdf(log_c[start_ind[g]:end_ind[g]] | alphas[g] +
			betas[g]*log_circumference_m[start_ind[g]:end_ind[g]], sigmas[g]);
	}
}

generated quantities {
	array [N] real pred_crown;
	for (g in 1:N_group)
		pred_crown[start_ind[g]:end_ind[g]] = exp(normal_rng(alphas[g] +
			betas[g]*log_circumference_m[start_ind[g]:end_ind[g]], sigmas[g]));
}

