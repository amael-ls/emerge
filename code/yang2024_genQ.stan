
/*
	Stan code to test the model from Yang 2024
*/

functions {
	vector r_yang_6(vector x, row_vector pars) // Function for the 6 parameters model
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> c,
			2 -> j,
			3 -> k,
			4 -> m,
			5 -> n,
			6 -> s
		*/
		return (pars[4] - pars[1]) * exp(pars[2] - pars[3]*x) .* (pars[3]*x/pars[2]).^pars[2] +
			pars[1] - (pars[1] - pars[5])*exp(-pars[6]*x);
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees

	// Predictors
	vector[N] bole_volume_m3;

	// Data
	vector[N] total_volume_m3;
}

transformed data {
	vector[N] ratio = bole_volume_m3 ./ total_volume_m3;
}

parameters {
	// Parameters of the 'bumpy' function r_yang_6
	real <lower = 0, upper = 1> c_beta;
	real <lower = 0> j; // Model still defined for j > 1, but it adds an unrealistic inflexion point
	real <lower = 0> k;
	real <lower = 0, upper = 1> m_beta;
	real <lower = 0, upper = 1> n;
	real <lower = 0> s_multiplier;

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	real c = 0.6 + 0.4*c_beta; // Forces c to be between 0.6 and 1
	real m = 0.8 + 0.2*m_beta; // Forces m to be between 0.8 and 1
	real s = (5 + s_multiplier)*k/j; // Force s to be at least 5
	vector [N] shape1 = phi*r_yang_6(bole_volume_m3, [c, j, k, m, n, s]);
	vector [N] shape2 = phi*(1 - r_yang_6(bole_volume_m3, [c, j, k, m, n, s]));
}

model{
	// Prior linear regression
	target += beta_lpdf(c_beta | 3, 3); // Centred
	target += normal_lpdf(j | 1, 0.1);
	target += gamma_lpdf(k | 2, 10); // Right skewed
	target += beta_lpdf(m_beta | 3, 3); // Centred
	target += beta_lpdf(n | 1, 8); // Right skewed
	target += gamma_lpdf(s_multiplier | 1.5, 0.5); // Right skewed

	target += gamma_lpdf(phi | 3, 0.5); // Right skewed
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}


generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N)
		v_gen[i] = 1/c * bole_volume_m3[i]^( 1 - (log(r_gen[i]) - log(c)) / log(bole_volume_m3[i]) );
	v_gen_mean = 1/c * bole_volume_m3 .^
		( 1 - (log(r_yang_6(bole_volume_m3, [c, j, k, m, n, s])) - log(c)) ./ log(bole_volume_m3) );
}
