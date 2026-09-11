
/*
	Stan code for submodel
*/

functions {
	vector r_4_params(vector x, row_vector pars) // Ratio 4 parameters model
	{
		/*
			The vector of parameters, pars, is in this order:
			1 -> logit_alpha, around 1.4 to be around 0.8 on the real scale
			2 -> beta,
			3 -> gamma
			4 -> delta
		*/
		return inv_logit(pars[1] + exp(-pars[2]*x) .* (pars[3]*x + pars[4]));
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
	// Parameters of the 'bumpy' function r_4_params
	real <lower = 0, upper = 1> c_alpha;
	real <lower = 0> beta_;
	real <lower = 0> gamma;
	real delta;

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	real alpha = 0.6 + 0.4*c_alpha; // Forces alpha (=c the asymptote) to be between 0.6 and 1
	real logit_alpha = logit(alpha); // pars[1] in r_4_params
	vector [N] shape1 = phi*r_4_params(bole_volume_m3, [logit_alpha, beta_, gamma, delta]);
	vector [N] shape2 = phi*(1 - r_4_params(bole_volume_m3, [logit_alpha, beta_, gamma, delta]));
}

model{
	// Prior linear regression
	target += beta_lpdf(c_alpha | 3, 3); // Centred
	target += gamma_lpdf(beta_ | 1, 0.5); // Gives a mean of 2 and sd of 2
	target += gamma_lpdf(gamma | 1, 0.5); // Gives a mean of 2 and sd of 2
	target += normal_lpdf(delta | -2.8, 1);

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
	v_gen_mean = 1/alpha *
		bole_volume_m3 .^ ( 1 - (log(r_4_params(bole_volume_m3, [logit_alpha, beta_, gamma, delta])) -
		log(alpha)) ./ log(bole_volume_m3) );
}
