/*
	This model is based on:
	Zhou X., Yang M., Liu Z., Li P., Xie B. & Peng C, 2021
	Dynamic allometric scaling of tree biomass and size
	Nature Plants 7
	https://doi.org/10.1038/s41477-020-00815-8
*/

functions {
	vector r_zhou(vector x, real alpha, real beta_, real b, real gamma)
	{
		// return alpha*x .* exp(-beta_*x) + b*(1 - exp(-gamma*x)); // First version, with flat tail
		// return alpha*x .* exp(-beta_*x) + b*x + gamma; // Second version, with linear trend, intercept != 0
		return alpha*exp(-beta_*x) .* (x - b/alpha) + gamma*x + b; // Third version, linear trend, intercept = 0
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
	real eps = 1e-6;

	vector[G] SB_max;
	for (i in 1:G)
		SB_max[i] = max(SB[start[i]:end[i]]);
}

parameters {
	// Parameters of the 'bumpy' function r_zhou
	array[G] real <lower = 0> alpha;
	array[G] real <lower = 0> beta_;
	vector <lower = 0, upper = 1> [G]  b;
	vector <lower = -b ./ SB_max, upper = 0> [G] gamma; // I expect d to be negative

	// Std. var.
	// real <lower = 0, upper = 0.1> sigma; // I put 0.09 which is 0.9*(1-  0.9), but then I put 0.1
}

transformed parameters {
	vector [N] r_mean;

	for (i in 1:G)
		r_mean[start[i]:end[i]] = r_zhou(SB[start[i]:end[i]], alpha[i], beta_[i], b[i], gamma[i]);

	// r_mean = (r_mean - min(r_mean))/(max(r_mean) - min(r_mean)) * (1 - 2*eps) + eps;
}

model{
	real sigma = 0.05;
	if (r_mean[1] * (1 - r_mean[1])/sigma^2 - 1 < 0)
	{
		print("-----------");
		print("alpha: ", alpha[1]);
		print("beta_: ", beta_[1]);
		print("b: ", b[1]);
		print("gamma: ", gamma[1]);
		print("sigma: ", sigma);
		print("min: ", min(r_mean .* (1 - r_mean)/sigma^2 - 1));
		print("max: ", max(r_mean .* (1 - r_mean)/sigma^2 - 1));
		print("-----------");
	}
	// Prior linear regression
	target += normal_lpdf(alpha | 0, 10);
	target += normal_lpdf(beta_ | 0, 10);
	target += normal_lpdf(b | 0.5, 15);
	target += normal_lpdf(gamma | 0, 10);
	target += uniform_lpdf(sigma | 0, 0.5);

	// Likelihood
	target += beta_lpdf(ratio | (r_mean .* (1 - r_mean)/sigma^2 - 1) .* r_mean,
		(r_mean .* (1 - r_mean)/sigma^2 - 1) .* (1 - r_mean));

	// target += normal_lpdf(ratio[start[i]:end[i]] | r_zhou(SB[start[i]:end[i]],
	// 	alpha[i], beta_[i], b[i], gamma[i]), sigma);
}

/*
generated quantities {
	array[N] real r_gen;

	for (i in 1:G)
		r_gen[start[i]:end[i]] = normal_rng(r_zhou(SB[start[i]:end[i]], alpha[i], beta_[i], b[i], gamma[i]), sigma);
}
*/
