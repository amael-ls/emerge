/*
	First test: Model Vtot directly, with a gamma distrib.
	I force the function to be 0 at (c, h) = (0, 0). This is physically correct,
	but is that the best?
*/

data {
	// Dimensions
	int <lower = 1> N;

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] height;

	// Observations
	vector <lower = 0> [N] total_volume_m3;
}

transformed data {
	// Std. pred
	vector [N] p1 = (circumference_m - mean(circumference_m))/sd(circumference_m);
	vector [N] p2 = (height - mean(height))/sd(height);

	// Form factor
	vector [N] formTot = 4*pi()*total_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
}

parameters {
	real beta0;
	real beta1;
	real beta2;

	real <lower = 0> sigma0;
	real sigma_pow;
}

transformed parameters {
	vector <lower = 0> [N] mu = exp(beta0 + beta1*p1 + beta2*p2);
	vector <lower = 0> [N] sigma = sigma0 * (circumference_m.^2 .* height) .^ sigma_pow;
}

model {
	// Priors
	target += normal_lpdf(beta0 | 0, 5);
	target += normal_lpdf(beta1 | 0, 5);
	target += normal_lpdf(beta2 | 0, 5);

	target += gamma_lpdf(sigma0 | 0.1, 0.1);
	target += normal_lpdf(sigma_pow | 0, 1);

	// Log-likelihood
	target += gamma_lpdf(formTot | mu.^2 ./ sigma.^2, mu ./ sigma.^2);
}

generated quantities {
	vector [N] Vpred;
	{
		real temp;
		for (i in 1:N)
		{
			temp = gamma_rng(mu[i]^2/sigma[i]^2, mu[i]/sigma[i]^2);
			Vpred[i] = temp*circumference_m[i]^2*height[i] / (4*pi()*(1 - 1.3/height[i])^2);
		}
	}
}
