functions {
	vector r_zhou(vector x, real k, real lambda)
	{
		return x.^k .* exp(-x / lambda);
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
}

parameters {
	// Parameters of the 'simple' function r_zhou (simple compared to the flexibility I want)
	vector <lower = 0> [G] lambda;
	vector <lower = 0, upper = exp(1) ./ lambda> [G] k;

	real <lower = 0> phi; // Precision (well kind of...)
}

transformed parameters {
	vector [N] shape1;
	vector [N] shape2;

	for (i in 1:G)
	{
		shape1[start[i]:end[i]] = phi*r_zhou(SB[start[i]:end[i]], k[i], lambda[i]);
		shape2[start[i]:end[i]] = phi*(1 - r_zhou(SB[start[i]:end[i]], k[i], lambda[i]));
	}
}

model{
	// Priors
	target += normal_lpdf(k | 0, 10);
	target += normal_lpdf(lambda | 0, 10);
	target += normal_lpdf(phi | 0, 10);
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}

generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
}


