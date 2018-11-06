data {
  int<lower=1> Ntotal;
  int x[Ntotal];
  real y[Ntotal];
  real meanY;
  real sdY;
}

transformed data {
  real unifLo;
  real unifHi;
  real normalSigma;
  unifLo = sdY/1000;
  unifHi = sdY*1000;
  normalSigma = sdY*100;
}

parameters {
    real mu[2];
    real<lower=0> sigma[2];
}

model {
    sigma ~ uniform(unifLo, unifHi);
    mu ~ normal(meanY, normalSigma);
    for (i in 1:Ntotal) {
        y[i] ~ normal(mu[x[i]] , sigma[x[i]]);
    }
}
