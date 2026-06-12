
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
	int <lower = 1> N_new; // Number of trees

	// Predictors
	vector[N] bole_volume_m3;
	vector[N_new] bole_volume_m3_new;

	// Data
	vector[N] total_volume_m3;
	vector[N_new] total_volume_m3_new;
}

transformed data {
	vector[N] ratio = bole_volume_m3 ./ total_volume_m3;
	vector[N_new] ratio_new = bole_volume_m3_new ./ total_volume_m3_new;
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

generated quantities {
	array[N_new] real r_gen;
	vector [N_new] r_gen_mean = r_yang_4(bole_volume_m3_new, [alpha, beta_, gamma, delta]); // Average ratio
	vector[N_new] v_gen;
	vector[N_new] v_gen_mean;
	vector [N_new] log_lik; // Log likelihood of newly observed volumes given fitted params on (other) data
	vector [N_new] sigma_var; // Variance

	{
		vector [N_new] shape1_new = phi*r_yang_4(bole_volume_m3_new, [alpha, beta_, gamma, delta]);
		vector [N_new] shape2_new = phi*(1 - r_yang_4(bole_volume_m3_new, [alpha, beta_, gamma, delta]));
		r_gen = beta_rng(shape1_new, shape2_new);

		for (i in 1:N_new)
		{
			v_gen[i] = 1/alpha * bole_volume_m3_new[i]^
				( 1 - (log(r_gen[i]) - log(alpha)) / log(bole_volume_m3_new[i]) );
			log_lik[i] = beta_lpdf(ratio_new[i] | shape1_new[i], shape2_new[i]);
		}
		v_gen_mean = 1/alpha * bole_volume_m3_new .^
			( 1 - (log(r_yang_4(bole_volume_m3_new, [alpha, beta_, gamma, delta])) - log(alpha)) ./
			log(bole_volume_m3_new) );

		sigma_var = shape1_new .* shape2_new ./ ((shape1_new + shape2_new).^2 .* (shape1_new + shape2_new + 1));
	}
}
