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
