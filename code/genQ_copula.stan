/*
	This file is to generate quantities according to a Clayton Copula on (Fbft, Vtot), where:
		- Fbft is the form factor for bole volume
		- Vtot is the total volume
	
	Here, Fbft follows a gamma distribution that was already studied in the report Gohon - Le Squin, while
	Vtot follows a lognormal distribution, with the circumference as a predictor only.

	For Fbft, I use the posteriors of the parameters obtained from NFI's data as a prior.In other
	words, I update these posteriors with the Emerge and Swiss data. Hopefully, that does not modify
	that much the value of the parameters (NFI's data are much larger)

	My source of information for this script are:
		1. Nelsen 2006, section 2.9 Random variate generation
		2. Stan's user guide https://mc-stan.org/docs/stan-users-guide/copulas.html#what-are-copulas
		3. Stan's helpful functions https://spinkney.github.io/helpful_stan_functions/group__clayton.html
*/

functions {
	real clayton_copula_lpdf(row_vector uv, real theta)
	{
		real result = log1p(theta) - (theta + 1) * (log(uv[1]) + log(uv[2])) -
			(2*theta + 1)/theta * log(pow(uv[1], - theta) + pow(uv[2], - theta) - 1); // Check PDF, ð²C(u, v)/ðuðv
		return result;
	}

	real C_u(real u, real t, real theta) // Conditional distrib of v given u, given Clayton copula
	{
		return u*pow(pow(u, theta) + pow(t, -theta/(theta + 1)) - 1, -1/theta);
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

	// real gamma_icdf(real u, real shape, real rate)
	// {
	// 	real x = shape / rate;
	// 	real err = 1.0;
	// 	int max_iter = 1000;
	// 	int i = 0;

	// 	// Root solver using Newton's method
	// 	while (abs(err) > 1e-6 && i < max_iter) {
	// 		err = gamma_p(shape, x * rate) - u;
	// 		x -= err / exp(gamma_lpdf(x | shape, rate));
	// 		i += 1;
	// 	}
	// 	return x;
	// }

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
	int <lower = 0> N; // Number of observations

	// Parameters priors Fbft
	vector[7] mu_priors; // beta0 -> 4, sigma0, sigma_pow
	vector <lower = 0> [7]sd_priors; // beta0 -> 4, sigma0, sigma_pow

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] taper_height_flo;
	vector <lower = 0> [N] height;

	// Observations
	vector <lower = 0> [N] bole_volume_m3;
	vector <lower = 0> [N] total_volume_m3;

	// New data
	// ...Dimensions
	int <lower = 0> N_new; // Number of new observations

	// ... Explanatory variables
	vector <lower = 0> [N_new] circumference_m_new;
	vector <lower = 0> [N_new] taper_height_flo_new;
	vector <lower = 0> [N_new] height_new;

	// ... Observations --> Uncomment when I will compute likelihood
	// vector <lower = 0> [N_new] bole_volume_m3_new;
	// vector <lower = 0> [N_new] total_volume_m3_new;
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

	// Predictors new data
	vector [N_new] p1_new = circumference_m_new; // Predictor 1
	vector [N_new] p2_new = sqrt(circumference_m_new) ./ taper_height_flo_new; // Predictor 2
	vector [N_new] p3_new = sqrt(taper_height_flo_new) ./ (circumference_m_new.^2 .* height_new); // Predictor 3
	vector [N_new] p4_new = 1 - taper_height_flo_new ./ height_new; // Predictor 4
	
	vector [N_new] ptot1_new = log(circumference_m_new); // Predictor 1 total volume
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
	
	vector [N] mu_tot = alpha0 + alpha1*ptot1;
}

generated quantities {
	vector [N_new] sim_bole_volume_m3;
	vector [N_new] sim_total_volume_m3;
	vector [N_new] shape_temp;
	vector [N_new] rate_temp;

	{
		// Clayton's copula variables
		real u;
		real t;
		real v;

		// Gamma distribution parameters
		real shape;
		real rate;
		
		// Transformed parameters
		vector [N_new] mu_temp = exp(beta0 + beta1*p1_new + beta2*p2_new +
			beta3*p3_new + beta4*p4_new);
		vector [N_new] sigma_temp = sigma0 *
			(circumference_m_new.^2 .* height_new) .^ sigma_pow;
		vector [N_new] mu_tot_temp = alpha0 + alpha1*ptot1_new;

		// Intermediate variables		
		real Fbft;

		for (i in 1:N_new)
		{
			// Sim Clayton's copula
			u = uniform_rng(0, 1);
			t = uniform_rng(0, 1);
			v = C_u(u, t, theta);

			// Transform back (u, v)...
			// ...For bole volume
			shape = mu_temp[i]^2 / sigma_temp[i]^2;
			rate = mu_temp[i] / sigma_temp[i]^2;
			Fbft = gamma_icdf(u, shape, rate);
			sim_bole_volume_m3[i] = Fbft * circumference_m_new[i]^2 * height_new[i] / 
				(4*pi()) / (1 - 1.3/height_new[i])^2;

			// ...For total volume
			sim_total_volume_m3[i] = lognormal_icdf(v, mu_tot_temp[i], sigma_tot_vol);

			// ... Temporary debugging stuff
			shape_temp[i] = shape;
			rate_temp[i] = rate;
		}
	}
}
