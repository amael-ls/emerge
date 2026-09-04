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
}

transformed data {
	vector[N] hdn = sqrt(circumference_cm) ./ height;
	vector[N] cylindre_vol = circumference_cm.^2 .* height/(40000*pi());
}

parameters {
	vector[N_params] vec_params;
	real<lower = 0> sigma;
}

model {
	target += normal_lpdf(vec_params[1] | 0, 1);
	target += normal_lpdf(vec_params[2] | 0, 0.001); // Circumference in cm, so not very large!

	if (N_params > 2 && is_douglas != 1)
		target += normal_lpdf(vec_params[3] | 0, 1);
	
	if (is_douglas == 1)
		target += normal_lpdf(vec_params[3] | 45, 10);

	if (N_params == 4)
		target += normal_lpdf(vec_params[4] | 45, 10);
	
	target += gamma_lpdf(sigma | 0.3^2/0.05, 0.3/0.05); // mean = 0.3, sd = 0.2
	
	target += normal_lpdf(total_volume_m3 | vallet(circumference_cm, hdn, vec_params, N_params, is_douglas) .*
		cylindre_vol, sigma);
}

