/*
	This model is based on:
	Zhou X., Yang M., Liu Z., Li P., Xie B. & Peng C, 2021
	Dynamic allometric scaling of tree biomass and size
	Nature Plants 7
	https://doi.org/10.1038/s41477-020-00815-8
*/

functions {
	// vector r_zhou(vector x, real alpha, real beta_, real b, real gamma)
	vector r_zhou(vector x, real k, real lambda)
	{
		// return alpha*x .* exp(-beta_*x) + b*(1 - exp(-gamma*x)); // First version, with flat tail
		// return alpha*x .* exp(-beta_*x) + b*x + gamma; // Second version, with linear trend, intercept != 0
		// return alpha*exp(-beta_*x) .* (x - b/alpha) + gamma*x + b; // Third version, linear trend, intercept = 0
		return x.^k .* exp(-x / lambda);
	}

	// Find to which group the min belongs to
	int which_group(vector x, array[] int ind_start)
	{
		int pos = 1;
		real minmin = x[1];
		for (i in 2:rows(x))
		{
			if (x[i] < minmin)
			{
				minmin = x[i];
				pos = i;
			}
		}

		int group = 1;
		while (group <= dims(ind_start)[1] && ind_start[group] < pos)
			group += 1;
		
		return group - 1;
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

	/*
	vector[G] SB_max;
	for (i in 1:G)
		SB_max[i] = max(SB[start[i]:end[i]]);
	*/
}

parameters {
	// Parameters of the 'bumpy' function r_zhou
	/*
	array[G] real <lower = 0> alpha;
	array[G] real <lower = 0> beta_;
	vector <lower = 0, upper = 1> [G]  b;
	vector <lower = -b ./ SB_max, upper = 0> [G] gamma; // I expect d to be negative
	*/

	vector <lower = 0> [G] lambda;
	vector <lower = 0, upper = exp(1) ./ lambda> [G] k;

	real <lower = 0> phi; // Precision (well kind of...)

	// Std. var.
	// real <lower = 0, upper = 0.1> sigma; // I put 0.09 which is 0.9*(1-  0.9), but then I put 0.1
}

transformed parameters {
	vector [N] shape1;
	vector [N] shape2;

	for (i in 1:G)
	{
		shape1[start[i]:end[i]] = phi*r_zhou(SB[start[i]:end[i]], k[i], lambda[i]);
		shape2[start[i]:end[i]] = phi*(1 - r_zhou(SB[start[i]:end[i]], k[i], lambda[i]));
	}
	// shape2[start[i]:end[i]] = 1.0/r_zhou(SB[start[i]:end[i]], alpha[i], beta_[i], b[i], gamma[i]) - 1;

	// r_mean = (r_mean - min(r_mean))/(max(r_mean) - min(r_mean)) * (1 - 2*eps) + eps;
}

model{
	/*
	int group = 0;
	if (min(shape2) < 0)
	{
		group = which_group(shape2, start);
		print("-----------");
		print("group: ", group);
		print("alpha: ", alpha[group]);
		print("beta_: ", beta_[group]);
		print("b: ", b[group]);
		print("gamma: ", gamma[group]);
		print("-----------");
	}
	*/

	// Prior linear regression
	/*
	target += normal_lpdf(alpha | 0, 10);
	target += normal_lpdf(beta_ | 0, 10);
	target += normal_lpdf(b | 0.5, 0.15);
	target += normal_lpdf(gamma | 0, 10);
	*/
	target += normal_lpdf(k | 0, 10);
	target += normal_lpdf(lambda | 0, 10);
	target += normal_lpdf(phi | 0, 10);
	
	// Likelihood
	target += beta_lpdf(ratio | shape1, shape2);

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
