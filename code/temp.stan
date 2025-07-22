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

	// Observations
	vector[N] taper_volume_m3;
}

/*
transformed data {
	vector[N] norm_circumference_m = circumference_m ./ sd(circumference_m);
	vector[N] norm_height = height ./ sd(height);
}
*/

parameters {
	real<lower = 5> max_potential_v_taper;
	real<lower = 0> loc_circum;
	real<lower = 0> loc_height;
	real<lower = 0> disp_circum;
	real<lower = 0> disp_height;
	real<lower = -1, upper = 1> rho; // Positive relation expected (for circum/height, negative for taper c/h?)
	
	real<lower = 0> sigma;
}

model {
	vector[6] vec_params = [max_potential_v_taper, loc_circum, loc_height, disp_circum, disp_height, rho]';

	target += gamma_lpdf(max_potential_v_taper | 10.0^2/20.0, 10.0/20.0); // Mean of 5, var of 10, loc circumference
	target += gamma_lpdf(loc_circum | 5.0^2/10.0, 5.0/10.0); // Mean of 5, var of 10, loc circumference
	target += gamma_lpdf(loc_height | 20.0^2/40.0, 20.0/40.0); // Mean of 20, var of 40, loc height
	target += gamma_lpdf(disp_circum | 0.1, 0.1); // Disp around circum
	target += gamma_lpdf(disp_height | 0.05, 0.05); // Disp around circum
	target += uniform_lpdf(rho | -1, 1); // Correlation circum/height
	
	target += gamma_lpdf(sigma | 5.0^2/10.0, 5.0/10.0); // Mean of 5, var of 10, for var around surface
	
	target += gamma_lpdf(taper_volume_m3 |
		surf_Vtaper(circumference_m, height, vec_params).^2/sigma^2,
		surf_Vtaper(circumference_m, height, vec_params)/sigma^2);
}

generated quantities {
	array[N] real taper_vol_gen;
	{
		vector[6] vec_params = [max_potential_v_taper, loc_circum, loc_height, disp_circum, disp_height, rho]';
		taper_vol_gen = gamma_rng(surf_Vtaper(circumference_m, height, vec_params).^2/sigma^2,
			surf_Vtaper(circumference_m, height, vec_params)/sigma^2);
	}
}


