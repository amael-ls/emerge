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
	vector [S] m;
	vector [S] d;
	vector [S] k;

	// Regression parameters around Zhou's volume
	vector[S] beta0;
	vector[S] beta1;
	vector[S] beta2;
	vector[S] beta3;
	vector[S] beta4;

	real <lower = 0> sigma;
}

model{
	// Define variables
	vector [N] zhou_ratio;
	for (i in 1:S)
		zhou_ratio[start[i]:end[i]] = bole_volume[start[i]:end[i]] ./ exp(m[i] - d[i]*exp(k[i]*D2H[start[i]:end[i]]));

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
		target += normal_lpdf(log_tot_volume[start[i]:end[i]] | beta0[i] +
			beta1[i]*log(zhou_ratio[start[i]:end[i]]) +
			beta2[i]*log(zhou_ratio[start[i]:end[i]].^2) +
			beta3[i]*log(zhou_ratio[start[i]:end[i]].^3) +
			beta4[i]*log(zhou_ratio[start[i]:end[i]].^4), sigma);
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
