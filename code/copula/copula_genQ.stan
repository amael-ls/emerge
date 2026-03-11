/*
	Script to simulate new data based on the Clayton's copula fitted in copula.R
*/

functions {
	real clayton_copula_lpdf(row_vector uv, real theta)
	{
		real result = log1p(theta) - (theta + 1) * (log(uv[1]) + log(uv[2])) -
			(2*theta + 1)/theta * log(pow(uv[1], - theta) + pow(uv[2], - theta) - 1); // Check PDF, ð²C(u, v)/ðuðv
		return result; // See page 132, notebook 3 IGN
	}

	real C_u(real u, real t, real theta) // Conditional distrib of v given u, given Clayton copula
	{
		return u*pow(pow(u, theta) + pow(t, -theta/(theta + 1)) - 1, -1/theta);
	}

	real lognormal_icdf(real p, real mu, real sigma) // inverse CDF lognormal
	{
		return exp(mu + sigma * inv_Phi(p));
	}
	
	real gamma_icdf(real u, real shape, real rate)
	{
		real x = shape / rate;
		real err = 1.0;
		int max_iter = 1000;
		int i = 0;

		while (abs(err) > 1e-6 && i < max_iter)
		{
			err = gamma_p(shape, x * rate) - u;
			x -= err / exp(gamma_lpdf(x | shape, rate));
			x = fmax(x, 1e-10); // prevent x from going negative
			i += 1;
		}
		return x;
	}
}

data {
	// Dimensions
	int <lower = 0> N_new;

	// Observation
	vector[N_new] Fbft_new;
	vector[N_new] Vtot_new;
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

generated quantities {
	vector [N_new] sim_Fbft;
	vector [N_new] sim_Vtot;

	{
		// Clayton's copula variables
		real u;
		real t;
		real v;

		for (i in 1:N_new)
		{
			// Sim Clayton's copula
			u = uniform_rng(0, 1);
			t = uniform_rng(0, 1);
			v = C_u(u, t, theta);

			// Transform back (u, v)...
			// ...For Fbft
			sim_Fbft[i] = gamma_icdf(u, shape, rate);

			// ...For Vtot
			sim_Vtot[i] = lognormal_icdf(v, meanlog, sdlog);
		}
	}
}
