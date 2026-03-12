/*
	This file is to test a Clayton Copula on (Fbft, Vtot) on simulated data, where:
		- Fbft is the form factor for bole volume
		- Vtot is the total volume
	
	Here, Fbft follows a gamma distribution with fixed shape and rate, while
	Vtot follows a lognormal distribution, with fixed meanlog and sdlog.

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

	// Observation
	vector[N] Fbft;
	vector[N] Vtot;
}

parameters {
	// Marginal of Fbft (gamma distrib)
	real <lower = 0> shape;
	real <lower = 0> rate;

	// Marginal of Vtot (lognormal distrib)
	real meanlog;
	real <lower = 0> sdlog;
	
	// Copula parameter
	real <lower = 0> theta; // Technically can go up to -1 but then non-strict generator...
}

model {
	// Priors
	target += gamma_lpdf(shape | 0.1, 0.1);
	target += gamma_lpdf(rate | 0.1, 0.1);
	target += normal_lpdf(meanlog | 0, 5);
	target += gamma_lpdf(sdlog | 0.1, 0.1);
	target += gamma_lpdf(theta | 1.5^2, 1.5); // mean of 1.5, var of 1

	// Log-likelihood...
	// ... Marginals
	target += gamma_lpdf(Fbft | shape, rate);
	target += lognormal_lpdf(Vtot | meanlog, sdlog);

	// ... Copula
	for (i in 1:N)
		target += clayton_copula_lpdf([gamma_cdf(Fbft[i] | shape, rate),
			lognormal_cdf(Vtot[i] | meanlog, sdlog)] | theta);
}
