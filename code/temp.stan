functions {
	vector surf_Vtaper(vector circumference, vector height, vector params)
	{
		/*
			Description of the parameters, stored in the vector params:
				- params[1] = Max potential taper volume
				- params[2] and [3] = location parameter of max potential in the circum/height space
				- params[4] and [5] = dispersion around location parameters
				- params[6] = correlation between the location/dispersion parameters
		*/
		return params[1] * exp(
			-(circumference - params[2]).^2/params[4]^2 -
			(height - params[3]).^2/params[5]^2 +
			params[6]/(1 - params[6]^2) * (circumference - params[2])/params[4] .* (height - params[3])/params[5]);
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees

	// Predictors
	vector[N] circumference_m;
	vector[N] height;
	vector[N] taper_height;

	// Observations
	vector[N] taper_volume_m3;
}


transformed data {
	vector[N] log_circumference_m = log(circumference_m);
	vector[N] log_height = log(height);
	vector[N] log_taper_height = log(taper_height);
	vector[N] log_taper_volume_m3 = log(taper_volume_m3);
}


parameters {
	real beta0;
	real beta1;
	real beta2;
	real beta3;
	
	// real<lower = 0> alpha;
	// real<lower = 0> beta_; // Increasing variance with volume expected
	real<lower = 0> sigma;
}

// transformed parameters {
// 	vector[N] sigma = alpha*(circumference_m.^2 .* height).^beta_;
// }

model {
	// vector[6] vec_params = [max_potential_v_taper, loc_circum, loc_height, disp_circum, disp_height, rho]';

	target += normal_lpdf(beta0 | 0, 10);
	target += normal_lpdf(beta1 | 0, 10);
	target += normal_lpdf(beta2 | 0, 10);
	target += normal_lpdf(beta3 | 0, 10);

	target += gamma_lpdf(sigma | 1.0^2/10.0, 1.0/10.0); // Mean of 5, var of 10, for var around surface
	
	target += normal_lpdf(log_taper_volume_m3 | beta0 + beta1*log_circumference_m + beta2*log_height +
		beta3*log_taper_height, sigma);
}

generated quantities {
	array[N] real taper_vol_gen = lognormal_rng(beta0 + beta1*log_circumference_m + beta2*log_height +
		beta3*log_taper_height, sigma);
}


