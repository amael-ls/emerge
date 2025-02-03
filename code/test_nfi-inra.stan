data {
	// Dimensions and indices
	int <lower = 1> N; // Total number of individuals
	int <lower = 1, upper = N> N_inra; // Number of individuals with a measured crown
	int <lower = 1, upper = N> S; // Number of species
	int <lower = 0, upper = S> n_sp_broad_inra; // Number of broadleaf species
	int <lower = 0, upper = S - n_sp_broad_inra> n_sp_conif_inra; // Number of conifer species

	int <lower = 1, upper = N - 1> ind_start_broad; // Broadleaf index start
	int <lower = 1, upper = N> ind_end_broad; // Broadleaf index end
	int <lower = 1, upper = N - 1> ind_start_conif; // Conifer index start
	int <lower = 1, upper = N> ind_end_conif; // Conifer index end
	
	array[n_sp_broad_inra] int ind_start_broad_inra; // Broadleaf crown species index start
	array[n_sp_broad_inra] int ind_end_broad_inra; // Broadleaf crown species index end
	array[n_sp_conif_inra] int ind_start_conif_inra; // Conifer crown species index start
	array[n_sp_conif_inra] int ind_end_conif_inra; // Conifer crown species index end

	// Predictors
	vector <lower = 0> [N] circumference_m;

	// Response variables
	vector <lower = 0> [N] bole_volume_m3;
	vector <lower = 0> [N] crown_volume_m3;
}

transformed data {
	// Predictors
	vector [N] log_circumference_m = log(circumference_m);

	// Response variables
	vector [N] log_bole_volume_m3 = log(bole_volume_m3);
	vector [N] log_crown_volume_m3 = log(crown_volume_m3);
}

parameters {
	// Intercepts
	vector [2] a_bole;
	vector [2] a_crown;

	// Slopes
	vector [2] b_bole;
	vector [2] b_crown;
	vector [2] c_bole;
	vector [2] c_crown;

	// Variances
	vector <lower = 0> [2] sigma_bole;
	vector <lower = 0> [2] sigma_crown;
}

model {
	// Priors
	// --- Intercepts
	target += normal_lpdf(a_bole | 0, 2);
	target += normal_lpdf(a_crown | 0, 2);

	// --- Slopes
	target += normal_lpdf(b_bole | 0, 2);
	target += normal_lpdf(b_crown | 0, 2);

	target += normal_lpdf(c_bole | 0, 2);
	target += normal_lpdf(c_crown | 0, 2);

	// --- Variances
	target += gamma_lpdf(sigma_bole | 1.0^2/5.0, 1.0/5.0);
	target += gamma_lpdf(sigma_crown | 1.0^2/5.0, 1.0/5.0);

	// Likelihood
	// --- Broadleaves, bole
	target += normal_lpdf(log_bole_volume_m3[ind_start_broad:ind_end_broad] | a_bole[1] + 
		b_bole[1]*log_circumference_m[ind_start_broad:ind_end_broad] +
		c_bole[1]*log_circumference_m[ind_start_broad:ind_end_broad].^4, sigma_bole[1]);
	
	// --- Conifers, bole
	target += normal_lpdf(log_bole_volume_m3[ind_start_conif:ind_end_conif] | a_bole[2] + 
		b_bole[2]*log_circumference_m[ind_start_conif:ind_end_conif] +
		c_bole[2]*log_circumference_m[ind_start_conif:ind_end_conif].^4, sigma_bole[2]);

	// --- Broadleaves, crown
	for (i in 1:n_sp_broad_inra)
		target += normal_lpdf(log_crown_volume_m3[ind_start_broad_inra[i]:ind_end_broad_inra[i]] | a_crown[1] + 
			b_crown[1]*log_circumference_m[ind_start_broad_inra[i]:ind_end_broad_inra[i]] +
			c_crown[1]*log_circumference_m[ind_start_broad_inra[i]:ind_end_broad_inra[i]].^4, sigma_crown[1]);
	
	// --- Conifers, crown
	for (i in 1:n_sp_conif_inra)
		target += normal_lpdf(log_crown_volume_m3[ind_start_conif_inra[i]:ind_end_conif_inra[i]] | a_crown[2] + 
			b_crown[2]*log_circumference_m[ind_start_conif_inra[i]:ind_end_conif_inra[i]] +
			c_crown[2]*log_circumference_m[ind_start_conif_inra[i]:ind_end_conif_inra[i]].^4, sigma_crown[2]);
}