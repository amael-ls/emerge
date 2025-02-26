
/*
Remarks:
	1. Note that there is no residual variance here! This is because for Bernoulli variables, B(p), the odds are p,
		and the variance is p(1 - p). So only one parameter!
	2. The non-centred parametrisation
*/

data {
	// Dimensions
	int <lower = 1> N; // Number of data
	int <lower = 1> N_chimp; // Number of chimpanzees
	int <lower = 1> N_block; // Number of blocks
	int <lower = 1> N_treatment; // Number of treatments
	int <lower = 2*(N_chimp - 1), upper = 2*(N_chimp - 1)> N_measure; // Number of measure per group (actor, block)

	// Explanatory variable
	array[N] int <lower = 1, upper = N_treatment> treatment; // Which treatment is applied for observation i
	vector <lower = 0, upper = 1> [N] left_prosocial; // Boolean 1 if left lever is prosocial, 0 otherwise
	vector <lower = 0, upper = N_block> [N] block_id; // Block id, from 1 to N_block

	// Observations
	array[N] int <lower = 0, upper = 1> left_pull; // Boolean, 1 if left lever pulled, 0 otherwise
}

parameters {
	// Population parameters...
	// ... intercept and slopes
	vector[N_treatment] gamma; // Intercept
	array[N_chimp] vector[N_treatment] Z_chimp; // Actor effect
	array[N_block] vector[N_treatment] Z_block; // Block effect
	
	// Variances (in the variance-covariance matrices)
	vector<lower = 0> [N_treatment] sigma_diag_actor;
	vector<lower = 0> [N_treatment] sigma_diag_block;
	cholesky_factor_corr[N_treatment] L_actor;
	cholesky_factor_corr[N_treatment] L_block;
}

transformed parameters {
	array[N_chimp] vector[N_treatment] alpha; // Slope actor
	array[N_block] vector[N_treatment] beta_; // Slope block

	for (ch in 1:N_chimp)
		alpha[ch] = diag_pre_multiply(sigma_diag_actor, L_actor) * Z_chimp[ch];
	
	for (bl in 1:N_block)
		beta_[bl] = diag_pre_multiply(sigma_diag_block, L_block) * Z_block[bl];
}

model {
	real odds = 0;
	int count = 1;
	// Priors...
	// ... variance
	target += lkj_corr_cholesky_lpdf(L_actor | 2);
	target += lkj_corr_cholesky_lpdf(L_block | 2);
	target += exponential_lpdf(sigma_diag_actor | 1);
	target += exponential_lpdf(sigma_diag_block | 1);

	// ... slopes
	target += multi_normal_lpdf(Z_chimp | rep_vector(0, N_treatment), identity_matrix(N_treatment));
	target += multi_normal_lpdf(Z_block | rep_vector(0, N_treatment), identity_matrix(N_treatment));

	// ... Intercept
	target += normal_lpdf(gamma | 0, 1); // Vectorised, applied to each element of the vector

	// Likelihood and cafe parameters
	for (ch in 1:N_chimp)
	{
		for (bl in 1:N_block)
		{
			for (i in 1:N_measure)
			{
				odds = gamma[treatment[count]] + alpha[ch][treatment[count]] + beta_[bl][treatment[count]];
				target += bernoulli_logit_lpmf(left_pull[count] | odds);
				count = count + 1;
			}
		}
	}
}

generated quantities {
	cov_matrix[N_treatment] Sigma_actor;
	cov_matrix[N_treatment] Sigma_block;
	vector[N] new_sim;
	
	// This is to recover the variance-covariance matrix instead of having its Cholesky factor
	Sigma_actor = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_diag_actor, L_actor));
	Sigma_block = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_diag_block, L_block));

	{
		int count = 1;
		real odds = 0;

		for (ch in 1:N_chimp)
		{
			for (bl in 1:N_block)
			{
				for (i in 1:N_measure)
				{
					odds = gamma[treatment[count]] + alpha[ch][treatment[count]] + beta_[bl][treatment[count]];
					new_sim[count] = bernoulli_logit_rng(odds);
					count = count + 1;
				}
			}
		}
	}
}


