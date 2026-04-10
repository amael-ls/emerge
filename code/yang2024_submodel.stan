
/*
	Stan code to test the model from Yang 2024
*/

functions {
	vector r_yang_3(vector x, row_vector pars) // Function for the simplified parameters model (3 params)
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> c,
			2 -> k,
			3 -> m,
		*/
		return (pars[3] - pars[1]) * pars[2] * x .* exp(1 - pars[2]*x) + pars[1];
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees
	// int <lower = 1, upper = N> G; // Number of genus
	// array[G] int <lower = 1, upper = N> start; // Genus index start
	// array[G] int <lower = 1, upper = N> end; // Genus index end

	// Predictors
	vector[N] SB;

	// Data
	vector[N] AGB;
}

transformed data {
	vector[N] ratio = SB ./ AGB;
}

parameters {
	// Parameters of the 'bumpy' function r_yang_3
	real <lower = 0, upper = 1> c;
	real <lower = 0> k;
	real <lower = 0, upper = 1> m;

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	vector [N] shape1 = phi*r_yang_3(SB, [c, k, m]);
	vector [N] shape2 = phi*(1 - r_yang_3(SB, [c, k, m]));
}

model{
	// Prior linear regression
	target += beta_lpdf(c | 8, 1); // Left skewed
	// target += normal_lpdf(j | 1, 0.1);
	target += beta_lpdf(k | 1, 8); // Right skewed
	target += beta_lpdf(m | 8, 1); // Left skewed
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}


generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N)
		v_gen[i] = 1/c * SB[i]^( 1 - (log(r_gen[i]) - log(c)) / log(SB[i]) );
	v_gen_mean = 1/c * SB .^ ( 1 - (log(r_yang_3(SB, [c, k, m])) - log(c)) ./ log(SB) );
}
