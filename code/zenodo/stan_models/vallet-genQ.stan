functions {
	vector vallet(vector circumference, vector hdn, vector params, int N_params, int is_douglas)
	{
		if (is_douglas == 1)
			return (params[1] + params[2]*circumference).*(1 + params[3]./circumference.^2);
		
		if (N_params == 4)
			return (params[1] + params[2]*circumference + params[3]*hdn).*(1 + params[4]./circumference.^2);
		
		if (N_params == 3)
			return params[1] + params[2]*circumference + params[3]*hdn;
		
		return params[1] + params[2]*circumference; // 2 parameters, Norway spruce
	}
}

data{
	// Dimensions
	int <lower = 1> N; // Number of trees
	int <lower = 2, upper = 4> N_params; // Number of parameters
	int <lower = 0, upper = 1> is_douglas; // 0 if not, and 1 otherwise

	// Predictors
	vector[N] circumference_cm; // In m, while in cm in Vallet2006
	vector[N] height;

	// Observations
	vector[N] total_volume_m3;
	
	// --------------   NEW OBSERVATIONS   --------------
	// Dimensions
	int <lower = 1> N_new;

	// Predictors
	vector[N_new] circumference_cm_new;
	vector[N_new] height_new;

	// Observations
	vector[N_new] total_volume_m3_new;
}

transformed data {
	vector[N] hdn = sqrt(circumference_cm) ./ height;
	vector[N_new] hdn_new = sqrt(circumference_cm_new) ./ height_new;
	vector[N] cylindre_vol = circumference_cm.^2 .* height/(40000*pi());
	vector[N] cylindre_vol_new = circumference_cm_new.^2 .* height_new/(40000*pi());
}

parameters {
	vector[N_params] vec_params;
	real<lower = 0> sigma;
}

generated quantities {
	vector[N_new] v_gen_mean = vallet(circumference_cm_new, hdn_new, vec_params, N_params, is_douglas) .*
		cylindre_vol_new;
	array[N_new] real v_gen = normal_rng(v_gen_mean, sigma);
	
	// Log likelihood of newly observed volumes given fitted params on (other) data
	vector [N_new] log_lik;
	
	for (i in 1:N_new)
		log_lik[i] = normal_lpdf(total_volume_m3_new[i] | v_gen_mean[i], sigma);
}

