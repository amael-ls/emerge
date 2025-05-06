
functions {
	vector link_fct(vector x1, vector x2, real b0, real b1, real b2)
	{
		return 1 ./ (1 + exp(-(b0 + b1*x1 + b2*x2)));
	}
}

data {
	// Dimensions
	int <lower = 1> N; // Number of trees

	// Predictors
	vector[N] x1;
	vector[N] x2;

	// Data
	vector[N] Y;
}

parameters {
	real b0;
	real b1;
	real b2;

	real<lower = 0> phi; // Precision
}

transformed parameters {
	vector[N] shape1 = phi*link_fct(x1, x2, b0, b1, b2);
	vector[N] shape2 = phi*(1 - link_fct(x1, x2, b0, b1, b2));
}

model{	
	// Prior linear regression
	target += normal_lpdf(b0 | 0, 1);
	target += normal_lpdf(b1 | 0, 1);
	target += normal_lpdf(b2 | 0, 1);

	target += gamma_lpdf(phi | 2.0^2/0.8, 2.0/0.8);

	// Likelihood
	target += beta_lpdf(Y | shape1, shape2);
}

