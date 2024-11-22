
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
	real b0;
	real c0;
	
	real b1;
	real c1;
	
	real b2;
	real c2;

	// Variance
	real<lower = 0> sigma; // sd residuals
}

model {

	// Priors
	target += normal_lpdf(b0 | 0, 10);
	target += normal_lpdf(b1 | 0, 10);
	target += normal_lpdf(b2 | 0, 10);
	target += normal_lpdf(c0 | 0, 10);
	target += normal_lpdf(c1 | 0, 10);
	target += normal_lpdf(c2 | 0, 10);

	target += inv_gamma_lpdf(sigma | 1, 1);

	for (i in 1:(lim_broadleaf - 1))
		target += normal_lpdf(volume_m3[i] | b0 + b1*fake_hdn[i] + b2*fake_slenderness[i], sigma);

	for (i in lim_broadleaf:N)
		target += normal_lpdf(volume_m3[i] | c0 + c1*fake_hdn[i] + c2*fake_slenderness[i], sigma);
}
