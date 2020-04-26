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
  real expLambda;         //New: parameter of prior for nu
  unifLo = sd_mu / 1000;
  unifHi = sd_mu * 1000;
  normalSigma = sd_mu * 100;
  expLambda= 1 / 29.0;      //New: setting value for expLambda
}

parameters {
  real<lower=0> nuMinusOne; //New: definition of additional parameter nu
  real mu;
  real<lower=0> sigma;
}

transformed parameters {
  real<lower=0> nu;           //New: new parameter nu
  nu = nuMinusOne + 1;           //New: shifting nu to avoid zero
}

model {
  sigma ~ uniform(unifLo, unifHi);
  mu ~ normal(mean_mu, normalSigma);
  nuMinusOne~exponential(expLambda);      //New: exponential prior for nu
  y ~ student_t(nu, mu, sigma);           //New: student_t distribution for nu
}

