/*
	This model is based on:
	Zhou X., Yang M., Liu Z., Li P., Xie B. & Peng C, 2021
	Dynamic allometric scaling of tree biomass and size
	Nature Plants 7
	https://doi.org/10.1038/s41477-020-00815-8
*/

data {
	// Dimensions
	int <lower = 1> N; // Number of trees
	int <lower = 1, upper = N> S; // Number of species
	array[S] int <lower = 1, upper = N> start; // Species index start
	array[S] int <lower = 1, upper = N> end; // Species index end

	// Predictors
	vector[N] circumference;
	vector[N] height;
	vector[N] bole_volume;

	// Data
	vector[N] tot_volume;
}

transformed data {
	vector[N] log_tot_volume = log(tot_volume);
	
	vector[N] DH = circumference.^2 .* height; // Proportional to a cylindre volume
}

parameters {
	vector <lower = 0, upper = 1> [S] m;
	vector <lower = 0, upper = 1> [S] d;
	vector <lower = 0> [S] k;

	real <lower = 0> sigma;
}

model{
	// Priors
	target += normal_lpdf(m | 0.9, 0.05);
	target += normal_lpdf(d | 0.3, 0.05);
	target += gamma_lpdf(k | 0.1, 0.1); // I put a minus in the exponential, hence k > 0
	
	target += gamma_lpdf(sigma | 0.25^2/0.03, 0.25/0.03);

	// Likelihood
	for (i in 1:S)
	{
		target += normal_lpdf(log_tot_volume[start[i]:end[i]] |
			log(bole_volume[start[i]:end[i]] ./ (m[i] - d[i]*exp(-k[i]*DH[start[i]:end[i]]))), sigma);
	}
}
