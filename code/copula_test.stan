/*
	This file is to test a Clayton Copula on (Fbft, Vtot) for Fagus sylvatica, where:
		- Fbft is the form factor for bole volume
		- Vtot is the total volume
	
	Here, Fbft follows a gamma distribution that was already studied in the report Gohon - Le Squin, while
	Vtot follows a lognormal distribution, with the circumference as a predictor only.

	For Fbft, I use the posteriors of the parameters obtained from NFI's data as a prior.In other
	words, I update these posteriors with the Emerge and Swiss data. Hopefully, that does not modify
	that much the value of the parameters (NFI's data are much larger)

	My source of information for this script are:
		1. Nelsen 2006, p. 116 -- 120
		2. Stan's user guide https://mc-stan.org/docs/stan-users-guide/copulas.html#what-are-copulas
		3. Stan's helpful functions https://spinkney.github.io/helpful_stan_functions/group__clayton.html
*/

functions {
	real clayton_copula_lpdf(real u, real v, real theta)
	{
		if (theta <= 0.0) 
			reject("clayton_copula: theta must > 0");
		
		real result = log1p(theta) - (theta + 1) * (log(u) + log(v)) -
			(2*theta + 1)/theta * log(pow(u, - theta) + pow(v, - theta) - 1) // Check PDF, ð²C(u, v)/ðuðv
		return result;
	}
}

data {
	// Dimensions
	int <lower = 0> N; // Number of observations

	// Parameters priors Fbft
	vector[7] mu_priors; // beta0 -> 4, sigma0, sigma_pow
	vector[7] <lower = 0> sd_priors; // beta0 -> 4, sigma0, sigma_pow

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] taper_height_flo;
	vector <lower = 0> [N] height;

	// Observations
	vector[N] bole_volume_m3;
	vector[N] total_volume_m3;
}

transformed data {
	// Predictors
	vector [N] p1 = circumference_m; // Predictor 1
	vector [N] p2 = sqrt(circumference_m) ./ taper_height_flo; // Predictor 2
	vector [N] p3 = sqrt(taper_height_flo) ./ (circumference_m.^2 .* height); // Predictor 3
	vector [N] p4 = 1 - taper_height_flo ./ height; // Predictor 4
	
	vector [N] ptot1 = log(circumference_m); // Predictor 1 total volume

	// Form factor
	vector [N] fnewbft = 4*pi()*bole_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
}

parameters {
	// Regression parameters Fbft
	real beta0;
	real beta1;
	real beta2;
	real beta3;
	real beta4;

	// Variance Fbft
	real <lower = 0> sigma0;
	real sigma_pow;

	// Regression parameters total volume
	real alpha0;
	real alpha1;

	// Variance total volume
	real <lower = 0> sigma_tot_vol;

	// Clayton copula's parameter
	real <lower = 0> theta;
}

transformed parameters {
	vector <lower = 0> [N] mu = exp(beta0 + beta1*p1 + beta2*p2 + beta3*p3 + beta4*p4);
	vector <lower = 0> [N] sigma = sigma0 * (circumference_m.^2 .* height) .^ sigma_pow;
	
	vector [N] mu_tot = alpha0 + alpha1*ptot;
}

model {
	// Priors...
	// ... Concerning Fbft
	target += normal_lpdf(beta0 | mu_priors[1], sd_priors[1]);
	target += normal_lpdf(beta1 | mu_priors[2], sd_priors[2]);
	target += normal_lpdf(beta2 | mu_priors[3], sd_priors[3]);
	target += normal_lpdf(beta3 | mu_priors[4], sd_priors[4]);
	target += normal_lpdf(beta4 | mu_priors[5], sd_priors[5]);
	
	target += gamma_lpdf(beta3 | mu_priors[6]^2/sd_priors[6]^2, mu_priors[6]^2/sd_priors[6]^2);
	target += normal_lpdf(beta4 | mu_priors[7], sd_priors[7]);

	// ... Concerning total volume
	target += normal_lpdf(alpha0 | 0, 5);
	target += normal_lpdf(alpha1 | 0, 5);
	
	target += normal_lpdf(sigma_tot_vol | 0.25/0.2^2, 0.5/0.2^2); // avg of 0.5, std dev of 0.2

	// ... Concerning copula
	target += gamma_lpdf(theta | 1.5^2, 1.5);

	// Log-likelihood...
	// ... Marginals
	target += gamma_lpdf(fnewbft | mu.^2 ./ sigma.^2, mu ./ sigma.^2);
	target += lognormal_lpdf(total_volume_m3 | mu_tot, sigma_tot_vol);

	// ... Copula
	for (n in 1:N)
		target += clayton_copula(gamma_cdf(fnewbft[i] | mu[i]^2/sigma[i]^2, mu[i]/sigma[i]^2),
			lognormal_cdf(total_volume_m3[i] | mu_tot[i], sigma_tot_vol), theta);
}


