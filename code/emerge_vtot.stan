
data {
	// Dimensions and indices
	int N; // Number of individuals
	int S; // Number of species
	array[S] int ind_start_sp; // Species index start
	array[S] int ind_end_sp; // Species index end
	int <lower = 2> lim_broadleaf; // Separation broadleaf/conifer

	// Predictors
	// array[N] real <lower = 0> height;
	// array[N] real <lower = 0> circumference_m;
	array[N] real fake_hdn;
	array[N] real fake_slenderness;

	// Response variable
	array[N] real volume_m3;
}

parameters {
	// Fixed effects for broadleaf/conifer
	array[2] real b0;
	array[2] real b1;
	array[2] real b2;

	// Group effect
	array[S] real beta0;
	// array[S] real beta2;

	// Variances
	real<lower = 0> sigma; // sd residuals
	// real<lower = 0> sigma_beta0; // sd beta0
	// real<lower = 0> sigma_beta2; // sd beta2
}

model {
	// Define local variables
	// real corrected_cyl_vol;
	real mu;
	int k = 1;

	// Priors
	target += normal_lpdf(b0 | 0, 10);
	target += normal_lpdf(b1 | 0, 10);
	target += normal_lpdf(b2 | 0, 10);

	// target += inv_gamma_lpdf(sigma | 1, 1); // Uses shape and scale
	// target += inv_gamma_lpdf(sigma_beta0 | 1, 1);
	// target += inv_gamma_lpdf(sigma_beta2 | 1, 1);

	target += gamma_lpdf(sigma | 0.1, 0.1); // Uses shape and rate
	// target += gamma_lpdf(sigma_beta0 | 0.1, 0.1);
	// target += gamma_lpdf(sigma_beta2 | 0.1, 0.1);

	// Likelihood, i = individual, j = species, k = genus
	for (j in 1:S)
	{
		for (i in ind_start_sp[j]:ind_end_sp[j])
		{
			if (i >= lim_broadleaf)
				k = 2;
			
			// Hierarchy
			target += normal_lpdf(beta0[j] | b0[k], 4.2);
			// target += normal_lpdf(beta2[j] | b2[k], sigma_beta2);

			// corrected_cyl_vol = height[i]*circumference_m[i]^2/(4*pi()*(1 - 1.3/height[i])^2);
			// mu = corrected_cyl_vol*(beta0[j] + b1[k]*sqrt(circumference_m[i])/height[i] +	
			// 	beta2[j]*height[i]/circumference_m[i]);
			
			// mu = beta0[j] + b1[k]*fake_hdn[i] + beta2[j]*fake_slenderness[i];
			mu = beta0[j] + b1[k]*fake_hdn[i] + b2[k]*fake_slenderness[i];
			
			
			target += normal_lpdf(volume_m3[i] | mu, sigma);
		}
	}
}
