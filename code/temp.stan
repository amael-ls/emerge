functions {
	vector r_zhou(vector x, real beta_, real gamma, real b)
	{
		return exp(-beta_*x) .* (x - b) + gamma*x + b; // Third version, linear trend, intercept = 0
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees
	int <lower = 1, upper = N> G; // Number of genus
	array[G] int <lower = 1, upper = N> start; // Genus index start
	array[G] int <lower = 1, upper = N> end; // Genus index end

	// Predictors
	vector[N] SB;

	// Data
	vector[N] AGB;
}

transformed data {
	vector[N] ratio = SB ./ AGB;

	vector[G] SB_max;
	for (i in 1:G)
		SB_max[i] = max(SB[start[i]:end[i]]);
}

parameters {
	// Parameters of the 'bumpy' function r_zhou
	vector <lower = 0> [G] beta_;
	vector <lower = 0, upper = 1> [G]  b;
	vector <lower = -b ./ SB_max, upper = 0> [G] gamma; // I expect d to be negative

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	vector [N] shape1;
	vector [N] shape2;

	for (i in 1:G)
	{
		shape1[start[i]:end[i]] = phi*r_zhou(SB[start[i]:end[i]], beta_[i], gamma[i], b[i]);
		shape2[start[i]:end[i]] = phi*(1 - r_zhou(SB[start[i]:end[i]], beta_[i], gamma[i], b[i]));
	}
}

model{
	// Prior linear regression
	target += normal_lpdf(beta_ | 0, 10);
	target += normal_lpdf(b | 0.5, 0.15);
	target += normal_lpdf(gamma | 0, 10);
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}

generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
}


