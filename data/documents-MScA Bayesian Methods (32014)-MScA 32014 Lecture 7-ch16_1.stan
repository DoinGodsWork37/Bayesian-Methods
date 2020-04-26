data {
    int<lower=1> Ntotal;
    real y[Ntotal];
    real meanY;
    real sdY;
}
transformed data {
    real unifLo;
    real unifHi;
    real normalSigma;
    unifLo = sdY/100;
    unifHi = sdY*100;
    normalSigma = sdY*100;
}
parameters {
    real mu;
    real<lower=0> sigma;
}
model {
    sigma ~ uniform(unifLo, unifHi);
    mu ~ normal(meanY, normalSigma); 
    y ~ normal(mu, sigma);
}