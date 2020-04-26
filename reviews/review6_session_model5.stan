data {
	int<lower=0> N;
	vector<lower=0>[N] earn;
	vector<lower=0>[N] height;
	vector[N] male;
}
transformed data {
	vector[N] log10_earn; // log transformation
	vector[N] height_c; // mean centered transformation
	vector[N] male_inter; // male interaction
	log10_earn = log10(earn);
	height_c = height - mean(height);
	male_inter = height_c .* male;
}
parameters {
	vector[4] beta;
	real<lower=0> sigma;
}
model {
  beta[1] ~ normal(4.5, 1);
	log10_earn ~ normal(beta[1] + beta[2] * height_c + beta[3] * male + beta[4] * male_inter, sigma);
}
