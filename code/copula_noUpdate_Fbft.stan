/*
	This file is to test a Clayton Copula on (Fbft, Vtot), where:
		- Fbft is the form factor for bole volume
		- Vtot is the total volume
	
	Here, Fbft follows a gamma distribution that was already studied in the report Gohon - Le Squin, while
	Vtot follows a lognormal distribution, with the circumference as a predictor only.

	For Fbft, I do not update the posteriors of the parameters obtained from NFI's data. Indeed,
		the Emerge dataset is non-representative of the French forest, and it would not makes sense
		to update parameters on a non-representative dataset. Secondly, the current parametrisation
		of the model (gamma distribution) prevents a good mixing of the parameters as the mean and
		variance are both involved in shape and rate. According to some experiment I made, it would
		have been better to use directly shape and rate rather than the method of moments. However,
		I also failed with shape and rate on non-simulated data...
	
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

	matrix shape_rate_mean(vector c, vector h, vector hdec, row_vector pars)
	{
		/*
			pars[1] -> beta0
			pars[2] -> beta1
			pars[3] -> beta2
			pars[4] -> beta3
			pars[5] -> beta4
			pars[6] -> sigma0
			pars[7] -> sigma_pow

			c is the circumference at breast height, h it total height, and hdec is the
			break taper height (height at which there is at least 10 % decrease within a meter)
		*/
		int N = rows(c);
		// print(N);
		vector [N] mu = exp(pars[1] + pars[2]*c + pars[3]*sqrt(c) ./ hdec +
			pars[4]*sqrt(hdec) ./ (c.^2 .* h) + pars[5]*(1 - hdec ./ h));
		vector [N] sigma = pars[6] * (c.^2 .* h) .^ pars[7];

		matrix [N, 2] res;
		res[:, 1] = mu.^2 ./ sigma.^2;
		res[:, 2] = mu ./ sigma.^2;

		return res;
	}
	
	array[] real F_bft_rng(vector c, vector h, vector hdec, row_vector pars)
	{
		/*
			pars[1] -> beta0
			pars[2] -> beta1
			pars[3] -> beta2
			pars[4] -> beta3
			pars[5] -> beta4
			pars[6] -> sigma0
			pars[7] -> sigma_pow

			c is the circumference at breast height, h it total height, and hdec is the
			break taper height (height at which there is at least 10 % decrease within a meter)
		*/
		int N = rows(c);
		vector [N] mu = exp(pars[1] + pars[2]*c + pars[3]*sqrt(c) ./ hdec +
			pars[4]*sqrt(hdec) ./ (c.^2 .* h) + pars[5]*(1 - hdec ./ h));
		vector [N] sigma = pars[6] * (c.^2 .* h) .^ pars[7];
		return gamma_rng(mu .^2 ./ sigma.^2, mu ./ sigma.^2);
	}
	
	array[] real F_bft_full_rng(vector c, vector h, vector hdec, array[, , ] real pars)
	{
		/*
			pars is of dimensions 1000 x 4 x 7, i.e.,
			n_iter x n_chains x n_params. I draw a random
			number between 1 and 4000, and convert it to
			a couple (row, col). I think it is better than
			drawing a row and col separately

			pars[, , 1] -> beta0
			pars[, , 2] -> beta1
			pars[, , 3] -> beta2
			pars[, , 4] -> beta3
			pars[, , 5] -> beta4
			pars[, , 6] -> sigma0
			pars[, , 7] -> sigma_pow

			c is the circumference at breast height, h it total height, and hdec is the
			break taper height (height at which there is at least 10 % decrease within a meter)
		*/
		int N = rows(c);
		int ind = to_int(uniform_rng(1, 4000.999)); // Do not use ceil or round which are of real type
		int rr = dims(pars)[1];
		int cc = dims(pars)[2];

		int i = ind % rr; // row
		int j = (ind - i) %/% rr + 1; // column, %/% fir integer division but that changes nothing

		vector [N] mu = exp(pars[i, j, 1] + pars[i, j, 2]*c + pars[i, j, 3]*sqrt(c) ./ hdec +
			pars[i, j, 4]*sqrt(hdec) ./ (c.^2 .* h) + pars[i, j, 5]*(1 - hdec ./ h));
		vector [N] sigma = pars[i, j, 6] * (c.^2 .* h) .^ pars[i, j, 7];
		return gamma_rng(mu .^2 ./ sigma.^2, mu ./ sigma.^2);
	}

	real C_u(real u, real t, real theta) // Conditional distrib of v given u, given Clayton copula
	{
		return u*pow(pow(u, theta) + pow(t, -theta/(theta + 1)) - 1, -1/theta);
	}

	real gamma_icdf2(real u, real shape, real rate) // The first version is unstable for 'small' N! uncertainty??
	{
		real lo = 0.0;
		real hi = 2; // Unlikely that Fbft > 2 in real! But maybe a bit dangerous in the real scenario?
		real mid;
		
		for (i in 1:20) // Dichotomy
		{
			mid = (lo + hi) / 2.0;
			if (gamma_cdf(mid | shape, rate) < u)
				lo = mid;
			else
				hi = mid;
		}
		return (lo + hi) / 2.0;
	}

	real lognormal_icdf(real p, real mu, real sigma) // inverse CDF lognormal
	{
		// Code from https://spinkney.github.io/helpful_stan_functions/group__qf.html
		/*
		if (is_nan(p) || p < 0 || p > 1) 
			reject("lognormal_icdf: p must be between 0 and 1; ", "found p = ", p);
		if (is_nan(mu) || is_inf(mu)) 
			reject("lognormal_icdf: mu must be finite; ", "found mu = ", mu);
		if (is_nan(sigma) || is_inf(sigma) || sigma <= 0) 
			reject("lognormal_icdf: sigma must be finite and > 0; ", "found sigma = ", sigma);
		*/
		return exp(mu + sigma * inv_Phi(p));
	}
}

data {
	// Dimensions
	int <lower = 0> N; // Number of observations
	// array [2] int <lower = 0, upper = N> ind_start; // Starting index for dataset random effect
	// array [2] int <lower = 0, upper = N> ind_end; // Ending index for dataset random effect

	// Regression parameters Fbft
	row_vector [7] params_Fbft; // In this order: beta 0 -> 4, sigma0, and lastly sigma_pow, mean values

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] taper_height_flo;
	vector <lower = 0> [N] height;

	// Observations
	vector <lower = 0> [N] bole_volume_m3; // Only used in the likelihood of generated F_bft
	vector <lower = 0> [N] total_volume_m3;
}

transformed data {
	// Predictors
	vector [N] ptot1 = log(circumference_m); // Predictor 1 total volume

	// Form factor (to compute likelihood of generated F bft)
	vector [N] fnewbft = 4*pi()*bole_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
}

parameters {
	// Regression parameters total volume
	real alpha0;
	real alpha1;

	// Variance total volume
	real <lower = 0> sigma_tot_vol;

	// Clayton copula's parameter
	real <lower = 0> theta;
}

transformed parameters {
	vector [N] mu_tot = alpha0 + alpha1*ptot1;
}

model {
	// Generate volumes
	matrix[N, 2] shape_rate = shape_rate_mean(circumference_m, height, taper_height_flo, params_Fbft);

	// print(N);
	// print(rows(circumference_m));
	// print(dims(circumference_m));
	// print(cols(circumference_m));
	// print(shape_rate[1:5, 1:2]);

	// Priors...
	// ... Concerning total volume
	target += normal_lpdf(alpha0 | 0, 1);
	target += normal_lpdf(alpha1 | 0, 1);
	
	target += gamma_lpdf(sigma_tot_vol | 0.25/0.2^2, 0.5/0.2^2); // avg of 0.5, std dev of 0.2

	// ... Concerning copula
	target += gamma_lpdf(theta | 1.5^2, 1.5); // avg 1.5, std dev of 1

	// Log-likelihood...
	// ... Marginals
	target += gamma_lpdf(fnewbft | shape_rate[, 1], shape_rate[, 2]);
	target += lognormal_lpdf(total_volume_m3 | mu_tot, sigma_tot_vol);

	// ... Copula
	for (i in 1:N)
		target += clayton_copula_lpdf([gamma_cdf(fnewbft[i] | shape_rate[i, 1], shape_rate[i, 2]),
			lognormal_cdf(total_volume_m3[i] | mu_tot[i], sigma_tot_vol)] | theta);
}

generated quantities {
	vector [N] sim_bole_volume_m3;
	vector [N] sim_total_volume_m3;

	{
		// Average simulated bole
		matrix[N, 2] shape_rate = shape_rate_mean(circumference_m, height, taper_height_flo, params_Fbft);

		// Clayton's copula variables
		real u;
		real t;
		real v;

		// Intermediate variables		
		real Fbft;

		for (i in 1:N)
		{
			// Sim Clayton's copula
			u = uniform_rng(0, 1);
			t = uniform_rng(0, 1);
			v = C_u(u, t, theta);

			// Transform back (u, v)...
			// ...For bole volume
			Fbft = gamma_icdf2(u, shape_rate[i, 1], shape_rate[i, 2]);
			sim_bole_volume_m3[i] = Fbft * circumference_m[i]^2 * height[i] / 
				(4*pi()) / (1 - 1.3/height[i])^2;

			// ...For total volume
			sim_total_volume_m3[i] = lognormal_icdf(v, mu_tot[i], sigma_tot_vol);
		}
	}
}
