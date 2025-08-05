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
	
	vector[N] D2H = circumference.^2 .* height; // Divide by pi^2 to keep prior from zhou: he works with dbh
}

parameters {
	// Parameters of Zhou2021 (Dynamic allometric scaling of tree biomass and size)
	real m;
	real d;
	real k;

	// Regression parameters around Zhou's volume
	real beta0;
	real beta1;
	real beta2;
	real beta3;
	real beta4;

	real <lower = 0> sigma;
}

model{
	// Define variables
	vector [N] zhou_ratio = bole_volume ./ exp(m - d*exp(k*D2H));

	// Prior linear regression
	target += normal_lpdf(beta0 | 0, 1);
	target += normal_lpdf(beta1 | 0, 1);
	target += normal_lpdf(beta2 | 0, 1);
	target += normal_lpdf(beta3 | 0, 1);
	target += normal_lpdf(beta4 | 0, 1);

	// Priors Zhou
	target += normal_lpdf(m | 0, 1);
	target += normal_lpdf(d | 0, 1);
	target += normal_lpdf(k | -1, 1);
	
	target += gamma_lpdf(sigma | 0.25^2/0.03, 0.25/0.03);

	// Likelihood
	for (i in 1:S)
	{
		target += normal_lpdf(log_tot_volume | beta0 +
			beta1*log(zhou_ratio) +
			beta2*log(zhou_ratio.^2) +
			beta3*log(zhou_ratio.^3) +
			beta4*log(zhou_ratio.^4), sigma);
	}
}

/*

target += lognormal_lpdf(tot_volume[start[i]:end[i]] | ) // Should be the same parameters, but does it work as well?

print("m = ", m[i]);
print("d = ", d[i]);
print("k = ", k[i]);
print("d^2h = ", D2H[i]);
print("z = ", zhou_ratio[i]);

*/
