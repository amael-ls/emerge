
data {
	// Dimensions and indices
	int N; // Number of individuals
	int lim_broadleaf; // Number of individuals

	// Predictors
	array[N] real <lower = 0> fake_hdn;
	array[N] real <lower = 0> fake_slenderness;

	// Response variable
	array[N] real volume_m3;
}

parameters {
	// Fixed effects
	array[2] real b0;
	array[2] real b1;
	array[2] real b2;

	// Variance
	real<lower = 0> sigma; // sd residuals
}

model {

	// Priors
	target += normal_lpdf(b0 | 0, 10);
	target += normal_lpdf(b1 | 0, 10);
	target += normal_lpdf(b2 | 0, 10);

	target += inv_gamma_lpdf(sigma | 1, 1);

	for (i in 1:(lim_broadleaf - 1))
		target += normal_lpdf(volume_m3[i] | b0[1] + b1[1]*fake_hdn[i] + b2[1]*fake_slenderness[i], sigma);

	for (i in lim_broadleaf:N)
		target += normal_lpdf(volume_m3[i] | b0[2] + b1[2]*fake_hdn[i] + b2[2]*fake_slenderness[i], sigma);
}
