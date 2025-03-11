/*
	This model is based on:
	Zhou X., Yang M., Liu Z., Li P., Xie B. & Peng C, 2021
	Dynamic allometric scaling of tree biomass and size
	Nature Plants 7
	https://doi.org/10.1038/s41477-020-00815-8
*/

data {
	// Parameters required for model
	int <lower = 1> N; // Number of trees
	int <lower = 1, upper = N> S; // Number of species

	//  Parameters for simulations...
	// ... Dimensions
	int <lower = 1> N_new; // Number of trees
	int <lower = 1, upper = S> S_new; // Number of selected species, cannot create a new species
	array[S_new] int <lower = 1, upper = S> selected_species; // Selected species among 1:S fitted species
	array[S_new] int <lower = 1, upper = N_new> start_new; // Species index start
	array[S_new] int <lower = 1, upper = N_new> end_new; // Species index end

	// ... Predictors
	vector[N_new] circumference_new;
	vector[N_new] height_new;
	vector[N_new] bole_volume_new;
}

transformed data {
	vector[N_new] DH_new = circumference_new.^2 .* height_new; // Proportional to a cylindre volume
}

parameters {
	vector <lower = 0, upper = 1> [S] m;
	vector <lower = 0, upper = 1> [S] d;
	vector <lower = 0> [S] k;

	real <lower = 0> sigma;
}

generated quantities {
	array[N_new] real <lower = 0> pred_tot_volume;
	vector[N_new] mean_pred;

	{
		int j = 0;
		for (i in selected_species)
		{
			j = j + 1;
			
			mean_pred[start_new[j]:end_new[j]] = log(bole_volume_new[start_new[j]:end_new[j]] ./
				(m[j] - d[j]*exp(-k[j]*DH_new[start_new[j]:end_new[j]])));
		}
	}
	pred_tot_volume = lognormal_rng(mean_pred, sigma);
}
