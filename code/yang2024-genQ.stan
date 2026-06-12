
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
	real m = c + (1 - c)*m_beta; // Forces m to be between c and 1
	real s = 5 + s_multiplier; // Force s to be at least 5
	vector [N] shape1 = phi*r_yang_6(bole_volume_m3, [c, j, k, m, n, s]);
	vector [N] shape2 = phi*(1 - r_yang_6(bole_volume_m3, [c, j, k, m, n, s]));
}

generated quantities {
	array[N_new] real r_gen;
	vector [N_new] r_gen_mean = r_yang_6(bole_volume_m3_new, [c, j, k, m, n, s]); // Average ratio
	vector[N_new] v_gen;
	vector[N_new] v_gen_mean;
	vector [N_new] log_lik; // Log likelihood of newly observed volumes given fitted params on (other) data
	vector [N_new] sigma_var; // Variance

	{
		vector [N_new] shape1_new = phi*r_yang_6(bole_volume_m3_new, [c, j, k, m, n, s]);
		vector [N_new] shape2_new = phi*(1 - r_yang_6(bole_volume_m3_new, [c, j, k, m, n, s]));
		r_gen = beta_rng(shape1_new, shape2_new);

		for (i in 1:N_new)
		{
			v_gen[i] = 1/c * bole_volume_m3_new[i]^( 1 - (log(r_gen[i]) - log(c)) / log(bole_volume_m3_new[i]) );
			log_lik[i] = beta_lpdf(ratio_new[i] | shape1_new[i], shape2_new[i]);
		}
		v_gen_mean = 1/c * bole_volume_m3_new .^
			( 1 - (log(r_yang_6(bole_volume_m3_new, [c, j, k, m, n, s])) - log(c)) ./ log(bole_volume_m3_new) );

		sigma_var = shape1_new .* shape2_new ./ ((shape1_new + shape2_new).^2 .* (shape1_new + shape2_new + 1));
	}
}


