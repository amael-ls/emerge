functions {
	vector r_zhou(vector x, real m, real k, real d)
	{
		return m - d * exp(-k * x); // Ratio from zhou 2021
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees
	int <lower = 1, upper = N> G; // Number of genus
	array[G] int <lower = 1, upper = N> start; // Genus index start
	array[G] int <lower = 1, upper = N> end; // Genus index end

	// Predictors
	vector[N] Vbft;

	// Data
	vector[N] Vtot;
}

transformed data {
	vector[N] ratio = Vbft ./ Vtot;
}

parameters {
	// Parameters m, d, and k from Zhou2021 (Dynamic allometric scaling of tree biomass and size)
	array[G] real <lower = 0.7, upper = 1> m;
	array[G] real <lower = 0.05, upper = m> d;
	array[G] real <lower = 0, upper = 0.3> k;

	// Precision beta-regression Zhou
	real <lower = 0> phi;
}

transformed parameters {
	vector [N] shape1;
	vector [N] shape2;

	for (i in 1:G)
	{
		shape1[start[i]:end[i]] = phi * r_zhou(Vbft[start[i]:end[i]], m[i], k[i], d[i]);
		shape2[start[i]:end[i]] = phi * (1 - r_zhou(Vbft[start[i]:end[i]], m[i], k[i], d[i]));
	}
}

model{
	// Prior Zhou formula
	target += normal_lpdf(k | 0, 0.1);
	target += normal_lpdf(m | 0.85, 0.05);
	target += normal_lpdf(d | 0.5, 0.15);
	target += gamma_lpdf(phi | 0.01, 0.01);

	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);
}


generated quantities {
	array[N] real r_gen = beta_rng(shape1, shape2);
}


