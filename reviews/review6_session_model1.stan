data {
	int<lower=0> N;
	vector[N] earn;
	vector[N] height;
}
transformed data {
	vector[N] log_earn; // log transformation
}
parameters {
	vector[2] beta;
	real<lower=0> sigma;
}
model {
	log_earn ~ normal(beta[1] + beta[2] * height, sigma);
}
