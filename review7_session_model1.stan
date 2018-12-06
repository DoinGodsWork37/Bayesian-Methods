data {
	int<lower=1> N;
	vector[N] salaries;
	vector[N] rank_one;
	vector[N] rank_two;
	vector[N] discipline;
	vector[N] sex;
	vector[N] yrs_service;
}
transformed data {
	vector[N] log10_sal; // log transformation
	vector[N] yrs_service_rank_one;
	vector[N] yrs_service_rank_two;
	vector[N] yrs_service_discipline;
	vector[N] yrs_service_sex;
	log10_sal = log10(salaries);
	yrs_service_rank_one = rank_one .* yrs_service;
	yrs_service_rank_two = rank_two .* yrs_service;
	yrs_service_discipline = discipline .* yrs_service;
	yrs_service_sex = sex .* yrs_service;
}
parameters {
	vector[6] beta;
	real<lower=0> sigma;
}
model {
  beta[1] ~ normal(5, 1.5);
	log10_sal ~ normal(beta[1] + beta[2] * yrs_service + beta[3] * yrs_service_rank_one +beta[4] * yrs_service_rank_two +
	beta[5] * yrs_service_discipline + beta[6] * yrs_service_sex, sigma);
}


