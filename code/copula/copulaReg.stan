/*
	This file is to test a Clayton Copula on (Fbft, Vtot) on simulated data, where:
		- Fbft is the form factor for bole volume
		- Vtot is the total volume
	
	Here, Fbft follows a gamma distribution with shape and rate determined by lin reg beta0 and 1, while
	Vtot follows a lognormal distribution, with meanlog and sdlog determined by lin reg alpha0 and 1.

	My source of information for this script are:
		1. Nelsen 2006, p. 116 -- 120
		2. Stan's user guide https://mc-stan.org/docs/stan-users-guide/copulas.html#what-are-copulas
		3. Stan's helpful functions https://spinkney.github.io/helpful_stan_functions/group__clayton.html
*/

functions {
	real clayton_copula_lpdf(row_vector uv, real theta)
	{
		real result = log1p(theta) - (theta + 1) * (log(uv[1]) + log(uv[2])) -
			(2*theta + 1)/theta * log(pow(uv[1], - theta) + pow(uv[2], - theta) - 1); // Check PDF, ð²C(u, v)/ðuðv
		return result; // See page 132, notebook 3 IGN
	}
}

data {
	// Dimensions
	int <lower = 0> N;

	// Provided parameter
	// real <lower = 0> theta;

	// Predictor
	vector <lower = 0> [N] circumference_m;

	// Observation
	vector <lower = 0> [N] Fbft;
	vector <lower = 0> [N] Vtot;
}

parameters {
	// Marginal of Fbft (gamma distrib)
	real beta0;
	real beta1;
	
	real <lower = 0> sdF;

	// Marginal of Vtot (lognormal distrib)
	real alpha0;
	real alpha1;

	real <lower = 0> sdlog;
	
	// Copula parameter
	real <lower = 0> theta; // Technically can go up to -1 but then non-strict generator...
}

transformed parameters {
	vector [N] shape = exp(sdF*circumference_m);
	vector [N] rate = exp(beta0 + beta1*circumference_m);
	
	vector [N] meanlog = alpha0 + alpha1*log(circumference_m);
}

model {
	// Priors
	target += normal_lpdf(beta0 | 0, 1);
	target += normal_lpdf(beta1 | 0, 1);

	target += normal_lpdf(alpha0 | 0, 1);
	target += normal_lpdf(alpha1 | 0, 1);

	target += gamma_lpdf(sdF | 1.2, 1);
	target += gamma_lpdf(sdlog | 0.15^2/0.05^2, 0.15/0.05^2);
	target += gamma_lpdf(theta | 1.5^2, 1.5); // mean of 1.5, var of 1

	// Log-likelihood...
	// ... Marginals
	target += gamma_lpdf(Fbft | shape, rate);
	target += lognormal_lpdf(Vtot | meanlog, sdlog);

	// ... Copula
	for (i in 1:N)
	{
		real u = gamma_cdf(Fbft[i] | shape[i], rate[i]);
		real v = lognormal_cdf(Vtot[i] | meanlog[i], sdlog);

		target += clayton_copula_lpdf([u, v] | theta);
	}
}

