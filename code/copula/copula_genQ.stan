/*
	Script to simulate new data based on the Clayton's copula fitted in copula.R
*/

data {
	// Dimensions
	int <lower = 0> N_new;

	// Observation
	vector[N_new] Fbft_new;
	vector[N_new] Vtot_new;
}
