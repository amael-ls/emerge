/*
	Model that uses the model of FLorence heterosked Gamma, joined with a crown model
*/

data {
	// Dimensions
	int <lower = 1> N_inra;
	int <lower = 1> N;

	// Explanatory
	vector <lower = 0> [N_inra] circumference_m_inra;
	vector <lower = 0> [N_inra] taper_height_inra;
	vector <lower = 0> [N_inra] height_inra;

	vector <lower = 0> [N] circumference_m;
	vector <lower = 0> [N] taper_height;
	vector <lower = 0> [N] height;

	// Observations
	vector <lower = 0> [N_inra] crown_volume_m3_inra;
	vector <lower = 0> [N_inra] bole_volume_m3_inra;
	vector <lower = 0> [N] bole_volume_m3;
}

transformed data {
	// Predictors bole volume
	vector [N] p1_inra = circumference_m_inra; // Predictor 1
	vector [N] p2_inra = sqrt(circumference_m_inra) ./ taper_height_flo_inra; // Predictor 2
	vector [N] p3_inra = sqrt(taper_height_flo_inra) ./ (circumference_m_inra.^2 .* height_inra); // Predictor 3
	vector [N] p4_inra = 1 - taper_height_flo_inra ./ height_inra; // Predictor 4

	vector [N] p1 = circumference_m; // Predictor 1
	vector [N] p2 = sqrt(circumference_m) ./ taper_height_flo; // Predictor 2
	vector [N] p3 = sqrt(taper_height_flo) ./ (circumference_m.^2 .* height); // Predictor 3
	vector [N] p4 = 1 - taper_height_flo ./ height; // Predictor 4

	// Form factor
	vector [N] fnewbft = 4*pi()*bole_volume_m3 ./ (circumference_m.^2 .* height) .* (1 - 1.3 ./ height).^2;
	vector [N] fnewbft_inra = 4*pi()*bole_volume_m3_inra ./ (circumference_m_inra.^2 .* height_inra) .*
		(1 - 1.3 ./ height_inra).^2;
}

