/*
	This model is based on:
	Zhou X., Yang M., Liu Z., Li P., Xie B. & Peng C, 2021
	Dynamic allometric scaling of tree biomass and size
	Nature Plants 7
	https://doi.org/10.1038/s41477-020-00815-8
*/

data {
	// Parameters required for model
	int <lower = 1, upper = N> S; // Number of species

	//  Parameters for simulations...
	// ... Dimensions
	int <lower = 1> N_new; // Number of trees
	int <lower = 1, upper = N> S_new; // Number of species
	array[S] int <lower = 1, upper = N> start_new; // Species index start
	array[S] int <lower = 1, upper = N> end_new; // Species index end

	// ... Predictors
	vector[N] circumference_new;
	vector[N] height_new;
	vector[N] bole_volume_new;
}

transformed data {
	vector[N] DH_new = circumference_new.^2 .* height_new; // Proportional to a cylindre volume
}

parameters {
	vector <lower = 0, upper = 1> [S] m;
	vector <lower = 0, upper = 1> [S] d;
	vector <lower = 0> [S] k;

	real <lower = 0> sigma;
}
