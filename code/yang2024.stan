
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
	// Parameters of the 'bumpy' function r_yang_6
	real <lower = 0, upper = 1> c;
	real <lower = 0, upper = 1> j; // Model still defined for j > 1, but it adds an unrealistic inflexion point
	real <lower = 0> k;
	real <lower = 0, upper = 1> m;
	real <lower = 0, upper = 1> n;
	real <lower = 0> s;

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	vector [N] shape1 = phi*r_yang_6(SB, [c, j, k, m, n, s]);
	vector [N] shape2 = phi*(1 - r_yang_6(SB, [c, j, k, m, n, s]));
}

model{
	// Prior linear regression
	target += beta_lpdf(c | 8, 1); // Left skewed
	target += beta_lpdf(j | 8, 1); // Left skewed
	target += beta_lpdf(k | 1, 8); // Right skewed
	target += beta_lpdf(m | 8, 1); // Left skewed
	target += beta_lpdf(n | 1, 8); // Right skewed
	target += gamma_lpdf(s | 3, 0.5); // Right skewed
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}


generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
	vector[N] v_gen;
	vector[N] v_gen_mean;

	for (i in 1:N)
		v_gen[i] = 1/c * SB[i]^( 1 - (log(r_gen[i]) - log(c)) / log(SB[i]) );
	v_gen_mean = 1/c * SB .^ ( 1 - (log(r_yang_6(SB, [c, j, k, m, n, s])) - log(c)) ./ log(SB) );
}
