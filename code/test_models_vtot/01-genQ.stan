/*
	Generator script for model 01
*/

data {
	// Dimensions
	int <lower = 1> N;
	int <lower = 1> N_new;

	// Explanatory variables
	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N_new] circumference_m_new;
	vector <lower = 0> [N] height;
	vector <lower = 0> [N_new] height_new;

	// Observations
	vector <lower = 0> [N] total_volume_m3;
	vector <lower = 0> [N_new] total_volume_m3_new;
}

transformed data {
	// Std. pred
	real m1 = mean(circumference_m);
	real m2 = mean(height);
	real s1 = sd(circumference_m);
	real s2 = sd(height);
	vector [N] p1 = (circumference_m - m1)/s1;
	vector [N] p2 = (height - m2)/s2;

	vector [N_new] p1_new = (circumference_m_new - m1)/s1;
	vector [N_new] p2_new = (height_new - m2)/s2;

	// Form factor
	vector [N] formTot = 4*pi()*total_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
	vector [N_new] formTot_new = 4*pi()*total_volume_m3_new ./ (circumference_m_new.^2 .* height_new) .*
		(1 - 1.3 ./ height_new).^2;
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
	vector [N_new] Vpred;
	vector [N_new] log_lik; // Log likelihood of newly observed volumes given fitted params on other data
	{
		real sigma_new = 0;
		real mu_new = 0;
		real temp = 0;
		real xi = 0;

		for (i in 1:N_new)
		{
			sigma_new = sigma0 * (circumference_m_new[i]^2 * height_new[i])^sigma_pow;
			mu_new = exp(beta0 + beta1*p1_new[i] + beta2*p2_new[i]);
			xi = circumference_m_new[i]^2*height_new[i]/(4*pi()*(1 - 1.3/height_new[i])^2);

			temp = gamma_rng(mu_new^2/sigma_new^2, mu_new/sigma_new^2);
			log_lik[i] = gamma_lpdf(total_volume_m3_new[i] | mu_new^2/sigma_new^2, mu_new/(xi*sigma_new^2));
			Vpred[i] = xi*temp;
		}
	}
}
