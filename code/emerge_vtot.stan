
data {
	int N; // Number of individuals

	// Predictors
	array[N] real <lower = 0> height;
	array[N] real <lower = 0> circumference_m;

	// Response variable
	array[N] real <lower = 0> volume_m3;
}

parameters {
	real beta0;
	real beta1;
	real beta2;
	real<lower = 0> sigma;
}

model {
	// Define local variables
	real corrected_cyl_vol;
	real mu;

	// Priors
	target += normal_lpdf(beta0 | 0, 10);
	target += normal_lpdf(beta1 | 0, 10);
	target += normal_lpdf(beta2 | 0, 10);

	target += gamma_lpdf(sigma | 0.1, 0.1);

	// Likelihood
	for (i in 1:N)
	{
		corrected_cyl_vol = height[i]*circumference_m[i]^2/(4*pi()*(1 - 1.3/height[i])^2);
		mu = corrected_cyl_vol*(beta0 + beta1*sqrt(circumference_m[i])/height[i] + beta2*height[i]/circumference_m[i]);
		target += normal_lpdf(volume_m3[i] | mu, sigma);
	}
}
