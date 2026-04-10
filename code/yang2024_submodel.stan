
/*
	Stan code to test the model from Yang 2024
*/

functions {
	vector r_yang_4(vector x, row_vector pars) // Function for the sub-model (4 parameters)
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> alpha,
			2 -> beta,
			3 -> gamma
			4 -> delta
		*/
		return pars[1] + exp(-pars[2]*x) .* (pars[3]*x + pars[4]);
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
	// Parameters of the 'bumpy' function r_yang_4
	real <lower = 0, upper = 1> c_alpha;
	real <lower = 0> beta_;
	real <lower = 0> gamma;
	real <upper = 0> delta; // Expected to be negative since n = alpha + delta and alpha > 0.6

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	real alpha = 0.6 + 0.4*c_alpha; // Forces alpha (=c the asymptote) to be between 0.6 and 1
	vector [N] shape1 = phi*r_yang_4(bole_volume_m3, [alpha, beta_, gamma, delta]);
	vector [N] shape2 = phi*(1 - r_yang_4(bole_volume_m3, [alpha, beta_, gamma, delta]));
}

model{
	// Prior linear regression
	target += beta_lpdf(c_alpha | 3, 3); // Centred
	target += gamma_lpdf(beta_ | 2, 5); // Right skewed
	target += gamma_lpdf(gamma | 2, 1); // Right skewed
	target += normal_lpdf(delta | -0.3, 0.075); // centred, mostly between -0.5 and -0.1

	target += gamma_lpdf(phi | 3, 0.5); // Right skewed
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}

generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N)
		v_gen[i] = 1/alpha * bole_volume_m3[i]^( 1 - (log(r_gen[i]) - log(alpha)) / log(bole_volume_m3[i]) );
	v_gen_mean = 1/alpha * bole_volume_m3 .^ ( 1 - (log(r_yang_4(bole_volume_m3, [alpha, beta_, gamma, delta])) -
		log(alpha)) ./ log(bole_volume_m3) );
}
