data {
	int<lower=0> N;
	vector[N] earn;
	vector[N] height;
}
transformed data {
	vector[N] log10_earn; // log transformation
	vector[N] height_c; // mean centered transformation
	log10_earn = log10(earn);
	height_c = height - mean(height);
}
parameters {
	vector[2] beta;
	real<lower=0> sigma;
}
model {
	log10_earn ~ normal(beta[1] + beta[2] * height_c, sigma);
}


