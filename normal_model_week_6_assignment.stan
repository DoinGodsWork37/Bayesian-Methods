data {
  int<lower=1> Ntotal;
  real y[Ntotal];
  real mean_mu;
  real sd_mu;
}

transformed data {
  real unifLo;
  real unifHi;
  real normalSigma;
  unifLo = sd_mu / 1000;
  unifHi = sd_mu * 1000;
  normalSigma = sd_mu * 100;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  sigma ~ uniform(unifLo, unifHi);
  mu ~ normal(mean_mu, normalSigma);
  y ~ normal(mu, sigma);
}




