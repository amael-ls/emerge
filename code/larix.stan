
/*
	Stan code to test the model from Yang 2024
*/

functions {
	vector amael(vector x, row_vector pars) // Function for the simplified parameters model (3 params)
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> a, asymptote
			2 -> b,
			3 -> c, two params for small trees (fast dynamics)
			4 -> d, I expect it to be negative
			5 -> h, two params to set up the asymptotic 'landing'
		*/
		return pars[1] + exp(-pars[2]*x).*(pars[3]*x + pars[4]);
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees

	// Predictors
	vector[N] SB;

	// Data
	vector[N] AGB;
}

transformed data {
	vector[N] ratio = SB ./ AGB;
}

parameters {
	// Parameters of the amael function
	real <lower = 0, upper = 1> a; // Asymptote modifier
	real <lower = 0> b;
	real <lower = 0> c;
	real d;
	
	real <lower = 0> phi;
}

transformed parameters {
	real asym = 0.7 + 0.3*a; // Forces asymptote to be between 0.7 and 1
	vector [N] shape1 = phi*amael(SB, [asym, b, c, d]);
	vector [N] shape2 = phi*(1 - amael(SB, [asym, b, c, d]));
}

model{
	// Prior linear regression
	target += beta_lpdf(a | 3, 3); // Centred
	target += gamma_lpdf(b | 3, 0.5); // Right skewed
	target += gamma_lpdf(c | 3, 0.5); // Right skewed
	target += normal_lpdf(d | 0, 1); // Centred
	
	target += gamma_lpdf(phi | 3, 0.5); // Right skewed
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}

generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N)
		v_gen[i] = 1/asym * SB[i]^( 1 - (log(r_gen[i]) - log(asym)) / log(SB[i]) );
	v_gen_mean = 1/asym * SB .^ ( 1 - (log(amael(SB, [asym, b, c, d])) - log(asym)) ./ log(SB) );
}
