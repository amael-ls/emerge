
/*
	Stan code for hierarchical model
*/

functions {
	vector r_6_params(vector x, row_vector pars) // Ratio 6 parameters model
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> c,
			2 -> j,
			3 -> tau,
			4 -> m,
			5 -> n,
			6 -> s
		*/
		return (pars[4] - pars[1]) * exp(pars[2] *(1 - x/pars[3])) .* (x/pars[3]).^pars[2] +
			pars[1] - (pars[1] - pars[5])*exp(-pars[6]*x);
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees
	int <lower = 1> N_plot; // Number of plots

	array[N_plot] int <lower = 1, upper = N> start_index; // Starting index for each plot
	array[N_plot] int <lower = 1, upper = N> end_index; // End index for each plot

	// Predictors
	vector[N] bole_volume_m3;

	// Data
	vector[N] total_volume_m3;
}

transformed data {
	vector[N] ratio = bole_volume_m3 ./ total_volume_m3;
}

parameters {
	// Hyper parameters for c
	real <lower = 0, upper = 1> mu_c; // population mean
	real <lower = 0> kappa_c; // Precision (well, kind of...)

	// Parameters of the 'bumpy' function r_6_params
	vector <lower = 0, upper = 1> [N_plot] c_beta;
	real <lower = 0> j; // Model still defined for j > 1, but it adds an unrealistic inflexion point
	real <lower = 0> tau; // bump location, replaces k
	real <lower = 0, upper = 1> m_beta;
	real <lower = 0, upper = 1> n;
	real <lower = 0> theta; // parameter to force mu_2 negligible at x = tau (see def of s)

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	vector[N_plot] c = 0.6 + 0.4*c_beta; // Forces c to be between 0.6 and 1
	vector[N_plot] m = c + (1 - c)*m_beta; // Forces m to be between c and 1
	real s = (5 + theta)/tau; // Force s to be at least 5*tau, i.e., at least m + exp[-5] for x = tau
	real k = j / tau; // deterministic transofrm, no Jacobian needed

	// Parameters of the beta distribution
	vector [N] shape1;
	vector [N] shape2;
	for (i in 1:N_plot)
	{
		shape1[start_index[i]:end_index[i]] = phi*r_6_params(
			bole_volume_m3[start_index[i]:end_index[i]], [c[i], j, tau, m[i], n, s]);
		shape2[start_index[i]:end_index[i]] = phi*(1 - r_6_params(
			bole_volume_m3[start_index[i]:end_index[i]], [c[i], j, tau, m[i], n, s]));
	}
}

model {
	// Hyper priors
	target += beta_lpdf(mu_c | 2, 2);
	target += gamma_lpdf(kappa_c | 2, 0.1);

	// Priors
	target += beta_lpdf(c_beta | mu_c * kappa_c, (1 - mu_c) * kappa_c); // Hierarchy
	target += lognormal_lpdf(j | -0.38, 1.06); // mean = 1.2, var = 3, 95% interval: 0.085 -- 5.47
	target += lognormal_lpdf(tau | -0.38, 1.06); // mean = 1.2, var = 3, 95% interval: 0.085 -- 5.47
	target += beta_lpdf(m_beta | 1, 8);
	target += beta_lpdf(n | 1, 8);
	target += gamma_lpdf(theta | 3, 0.5); // 95% interval: 1.24 -- 14.45
	target += gamma_lpdf(phi | 3, 0.5);

	target += beta_lpdf(ratio | shape1, shape2);
}

generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N_plot)
	{
		v_gen_mean[start_index[i]:end_index[i]] = bole_volume_m3[start_index[i]:end_index[i]] ./
			r_6_params(bole_volume_m3[start_index[i]:end_index[i]], [c[i], j, tau, m[i], n, s]);
	}

	for (i in 1:N)
		v_gen[i] = bole_volume_m3[i] ./ r_gen[i];
}
