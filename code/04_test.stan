
data {
	// Dimensions
	int <lower = 1> N_cafes; // Number of cafes
	int <lower = 1> N_visit; // Number of visits per cafe
	int <lower = N_cafes*N_visit, upper = N_cafes*N_visit> N; // Number of data

	// Explanatory variable
	vector <lower = 0, upper = 1> [N] afternoon; // 0 if morning, 1 otherwise

	// Observations
	vector <lower = 0> [N] wait;
}

parameters {
	// Population parameters...
	// ... intercept and slope
	real a; // Intercept
	real b; // Slope
	
	// ... variances
	real sigma_a;
	real sigma_b;
	cholesky_factor_corr[2] L;

	// Variance (residuals)
	real<lower = 0> sigma; // Variance within same cafe

	// Whatever
	vector [N_cafes] a_cafe;
	vector [N_cafes] b_cafe;
}

transformed parameters {
	vector <lower = 0>[2] sigma_vec = to_vector({sigma_a, sigma_b});
}

model {
	int cc = 0;

	// Population priors
	target += normal_lpdf(a | 5, 2);
	target += normal_lpdf(b | -1, 2);
	
	target +=  lkj_corr_cholesky_lpdf(L | 2); // It contains rho (non-diag)

	// Variance (residuals)
	target += exponential_lpdf(sigma | 1);

	// Likelihood and cafe parameters
	for (i in 1:N_cafes)
	{
		target += multi_normal_cholesky_lpdf(to_vector({a_cafe[i], b_cafe[i]}) | to_vector({a, b}),
			diag_pre_multiply(sigma_vec, L)); // Cafe params
		for (j in 1:N_visit)
		{
			cc = (i - 1)*N_visit + j;
			target += normal_lpdf(wait[cc] | a_cafe[i] + b_cafe[i]*afternoon[cc], sigma); // Likelihood
		}
	}
}

generated quantities {
	cov_matrix[2] Sigma;
	real rho;
	
	// This is to recover the variance-covariance matrix instead of having its Cholesky factor
	Sigma = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_vec, L));

	// This is to recover the correlation parameter
	rho = tcrossprod(L)[2, 1];
}
