data {
  int<lower=1> N;  // number of observations
  int<lower=0,upper=1> y[N]; // response variable
  int<lower=1> NSex;
  int<lower=1, upper=NSex> sex[N];
  int<lower=1> NRace;
  int<lower=1, upper=NRace> race[N];
  int<lower=1> NAge;
  int<lower=1, upper=NAge> age[N];
  int<lower=1> NEducation;
  int<lower=1, upper=NEducation> education[N];
  int<lower=1> NState;
  int<lower=1, upper=NState> state[N];
}
parameters {
  real<lower=0> nu;
  real<lower=0> sigma;
  // temporary Intercept:
  real a_0;
  real<lower=0> var_0;
  // one-group alphas:
  vector[NSex]  a_sex;
  vector[NRace] a_race;
  vector[NAge]  a_age;
  vector[NEducation] a_education;
  vector[NState] a_state;
  // one-group variances:
  real<lower=0> var_sex;
  real<lower=0> var_race;
  real<lower=0> var_age;
  real<lower=0> var_education;
  real<lower=0> var_state;
  // two-group alphas (interactions):
  matrix[NSex,NRace] a_sex_race;
  matrix[NSex,NAge] a_sex_age;
  matrix[NSex,NEducation] a_sex_education;
  matrix[NSex,NState] a_sex_state;
  matrix[NRace,NAge] a_race_age;
  matrix[NRace,NEducation] a_race_education;
  matrix[NRace,NState] a_race_state;
  matrix[NAge,NEducation] a_age_education;
  matrix[NAge,NState] a_age_state;
  matrix[NEducation,NState] a_education_state;
  // two-group variances:
  real<lower=0> var_sex_race;
  real<lower=0> var_sex_age;
  real<lower=0> var_sex_education;
  real<lower=0> var_sex_state;
  real<lower=0> var_race_age;
  real<lower=0> var_race_education;
  real<lower=0> var_race_state;
  real<lower=0> var_age_education;
  real<lower=0> var_age_state;
  real<lower=0> var_education_state;
}
model {
  // prior specifications:
  // nu & sigma have vague improper uniform priors
  // intercept:
  a_0 ~ normal(0, sqrt(var_0));
  var_0 ~ scaled_inv_chi_square(nu, sigma);
  // one-group alphas:
  a_sex ~ normal(0, sqrt(var_sex));
  a_race ~ normal(0, sqrt(var_race));
  a_age ~ normal(0, sqrt(var_age));
  a_education ~ normal(0, sqrt(var_education));
  a_state ~ normal(0, sqrt(var_state));
  // one-group variances:
  var_sex ~ scaled_inv_chi_square(nu, sigma);
  var_race ~ scaled_inv_chi_square(nu, sigma);
  var_age ~ scaled_inv_chi_square(nu, sigma);
  var_education ~ scaled_inv_chi_square(nu, sigma);
  var_state ~ scaled_inv_chi_square(nu, sigma);
  // two-group alphas (interactions):
  for (j in 1:NSex) {
    a_sex_race[j,] ~ normal(0, sqrt(var_sex_race));
    a_sex_age[j,] ~ normal(0, sqrt(var_sex_age));
    a_sex_education[j,] ~ normal(0, sqrt(var_sex_education));
    a_sex_state[j,] ~ normal(0, sqrt(var_sex_state));
  }
  for (k in 1:NRace) {
    a_race_age[k,] ~ normal(0, sqrt(var_race_age));
    a_race_education[k,] ~ normal(0, sqrt(var_race_education));
    a_race_state[k,] ~ normal(0, sqrt(var_race_state));
  }
  for (l in 1:NAge) {
    a_age_education[l,] ~ normal(0, sqrt(var_age_education));
    a_age_state[l,] ~ normal(0, sqrt(var_age_state));
  }
  for (m in 1:NEducation) {
    a_education_state[m,] ~ normal(0, sqrt(var_education_state));
  }
  // two-group variances:
  var_sex_race ~ scaled_inv_chi_square(nu, sigma);
  var_sex_age ~ scaled_inv_chi_square(nu, sigma);
  var_sex_education ~ scaled_inv_chi_square(nu, sigma);
  var_sex_state ~ scaled_inv_chi_square(nu, sigma);
  var_race_age ~ scaled_inv_chi_square(nu, sigma);
  var_race_education ~ scaled_inv_chi_square(nu, sigma);
  var_race_state ~ scaled_inv_chi_square(nu, sigma);
  var_age_education ~ scaled_inv_chi_square(nu, sigma);
  var_age_state ~ scaled_inv_chi_square(nu, sigma);
  var_education_state ~ scaled_inv_chi_square(nu, sigma);

  // response var:
  for (i in 1:N) {
    y[i] ~ bernoulli_logit(a_0 + a_sex[sex[i]] + a_race[race[i]]
                           + a_age[age[i]] + a_education[education[i]] + a_state[state[i]]
                           + a_sex_race[sex[i],race[i]] + a_sex_age[sex[i],age[i]]
                           + a_sex_education[sex[i],education[i]] + a_sex_state[sex[i],state[i]]
                           + a_race_age[race[i],age[i]] + a_race_education[race[i],education[i]]
                           + a_race_state[race[i],state[i]]
                           + a_age_education[age[i],education[i]] + a_age_state[age[i],state[i]]
                           + a_education_state[education[i],state[i]]
                          );
  }
}
generated quantities {
  real b_0;  // fixed effects intercept
  // centered REs:
  vector[NSex] b_sex;
  vector[NRace] b_race;
  vector[NAge] b_age;
  vector[NEducation] b_education;
  vector[NState] b_state;
  matrix[NSex,NRace] b_sex_race;
  matrix[NSex,NAge] b_sex_age;
  matrix[NSex,NEducation] b_sex_education;
  matrix[NSex,NState] b_sex_state;
  matrix[NRace,NAge] b_race_age;
  matrix[NRace,NEducation] b_race_education;
  matrix[NRace,NState] b_race_state;
  matrix[NAge,NEducation] b_age_education;
  matrix[NAge,NState] b_age_state;
  matrix[NEducation,NState] b_education_state;

  real mtrx[NSex,NRace,NAge,NEducation,NState];

  for (j in 1:NSex) {
    for (k in 1:NRace) {
      for (l in 1:NAge) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              mtrx[j,k,l,m,n] = a_0 + a_sex[j] + a_race[k] + a_age[l] + a_education[m] + a_state[n]
                               + a_sex_race[j,k] + a_sex_age[j,l] + a_sex_education[j,m] + a_sex_state[j,n]
                               + a_race_age[k,l] + a_race_education[k,m] + a_race_state[k,n]
                               + a_age_education[l,m] + a_age_state[l,n]
                               + a_education_state[m,n];
          }
        }
      }
    }
  }

  // b_0
  b_0 = a_0 + sum(a_sex)/NSex + sum(a_race)/NRace + sum(a_age)/NAge
        + sum(a_education)/NEducation + sum(a_state)/NState
        + sum(a_sex_race)/NSex/NRace + sum(a_sex_age)/NSex/NAge
        + sum(a_sex_education)/NSex/NEducation + sum(a_sex_state)/NSex/NState
        + sum(a_race_age)/NRace/NAge + sum(a_race_education)/NRace/NEducation
        + sum(a_race_state)/NRace/NState + sum(a_age_education)/NAge/NEducation
        + sum(a_age_state)/NAge/NState + sum(a_education_state)/NEducation/NState;
  // b_sex
  for (j in 1:NSex) {
    b_sex[j] = 0;
    for (k in 1:NRace) {
      for (l in 1:NAge) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_sex[j] = b_sex[j] + mtrx[j,k,l,m,n];
          }
        }
      }
    }
    b_sex[j] = b_sex[j]/NRace/NAge/NEducation/NState - b_0;
  }
  // b_race
  for (k in 1:NRace) {
    b_race[k] = 0;
    for (j in 1:NSex) {
      for (l in 1:NAge) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_race[k] = b_race[k] + mtrx[j,k,l,m,n];
          }
        }
      }
    }
    b_race[k] = b_race[k]/NSex/NAge/NEducation/NState - b_0;
  }
  // b_age
  for (l in 1:NAge) {
    b_age[l] = 0;
    for (j in 1:NSex) {
      for (k in 1:NRace) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_age[l] = b_age[l] + mtrx[j,k,l,m,n];
          }
        }
      }
    }
    b_age[l] = b_age[l]/NSex/NRace/NEducation/NState - b_0;
  }
  // b_education
  for (m in 1:NEducation) {
    b_education[m] = 0;
    for (j in 1:NSex) {
      for (k in 1:NRace) {
        for (l in 1:NAge) {
          for (n in 1:NState) {
              b_education[m] = b_education[m] + mtrx[j,k,l,m,n];
          }
        }
      }
    }
    b_education[m] = b_education[m]/NSex/NRace/NAge/NState - b_0;
  }
  // b_state
  for (n in 1:NState) {
    b_state[n] = 0;
    for (j in 1:NSex) {
      for (k in 1:NRace) {
        for (l in 1:NAge) {
          for (m in 1:NEducation) {
              b_state[n] = b_state[n] + mtrx[j,k,l,m,n];
          }
        }
      }
    }
    b_state[n] = b_state[n]/NSex/NRace/NAge/NEducation - b_0;
  }

  // interactions
  // b_sex_race
  for (j in 1:NSex) { for (k in 1:NRace) {
      b_sex_race[j,k] = 0;
      for (l in 1:NAge) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_sex_race[j,k] = b_sex_race[j,k] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_sex_race[j,k] = b_sex_race[j,k]/NAge/NEducation/NState - b_0 - b_sex[j] - b_race[k];
    }
  }
  // b_sex_age
  for (j in 1:NSex) { for (l in 1:NAge) {
      b_sex_age[j,l] = 0;
      for (k in 1:NRace) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_sex_age[j,l] = b_sex_age[j,l] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_sex_age[j,l] = b_sex_age[j,l]/NRace/NEducation/NState - b_0 - b_sex[j] - b_age[l];
    }
  }
  // b_sex_education
  for (j in 1:NSex) { for (m in 1:NEducation) {
      b_sex_education[j,m] = 0;
      for (k in 1:NRace) {
        for (l in 1:NAge) {
          for (n in 1:NState) {
              b_sex_education[j,m] = b_sex_education[j,m] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_sex_education[j,m] = b_sex_education[j,m]/NRace/NAge/NState - b_0 - b_sex[j] - b_education[m];
    }
  }
  // b_sex_state
  for (j in 1:NSex) { for (n in 1:NState) {
      b_sex_state[j,n] = 0;
      for (k in 1:NRace) {
        for (l in 1:NAge) {
          for (m in 1:NEducation) {
              b_sex_state[j,n] = b_sex_state[j,n] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_sex_state[j,n] = b_sex_state[j,n]/NRace/NAge/NEducation - b_0 - b_sex[j] - b_state[n];
    }
  }

  // b_race_age
  for (k in 1:NRace) { for (l in 1:NAge) {
      b_race_age[k,l] = 0;
      for (j in 1:NSex) {
        for (m in 1:NEducation) {
          for (n in 1:NState) {
              b_race_age[k,l] = b_race_age[k,l] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_race_age[k,l] = b_race_age[k,l]/NSex/NEducation/NState - b_0 - b_race[k] - b_age[l];
    }
  }
  // b_race_education
  for (k in 1:NRace) { for (m in 1:NEducation) {
      b_race_education[k,m] = 0;
      for (j in 1:NSex) {
        for (l in 1:NAge) {
          for (n in 1:NState) {
              b_race_education[k,m] = b_race_education[k,m] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_race_education[k,m] = b_race_education[k,m]/NSex/NAge/NState - b_0 - b_race[k] - b_education[m];
    }
  }
  // b_race_state
  for (k in 1:NRace) { for (n in 1:NState) {
      b_race_state[k,n] = 0;
      for (j in 1:NSex) {
        for (l in 1:NAge) {
          for (m in 1:NEducation) {
              b_race_state[k,n] = b_race_state[k,n] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_race_state[k,n] = b_race_state[k,n]/NSex/NAge/NEducation - b_0 - b_race[k] - b_state[n];
    }
  }

  // b_age_education
      for (l in 1:NAge) { for (m in 1:NEducation) {
      b_age_education[l,m] = 0;
      for (j in 1:NSex) {
        for (k in 1:NRace) {
          for (n in 1:NState) {
              b_age_education[l,m] = b_age_education[l,m] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_age_education[l,m] = b_age_education[l,m]/NSex/NRace/NState - b_0 - b_age[l] - b_education[m];
    }
  }
  // b_age_state
  for (l in 1:NAge) { for (n in 1:NState) {
      b_age_state[l,n] = 0;
      for (j in 1:NSex) {
        for (k in 1:NRace) {
          for (m in 1:NEducation) {
              b_age_state[l,n] = b_age_state[l,n] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_age_state[l,n] = b_age_state[l,n]/NSex/NRace/NEducation - b_0 - b_age[l] - b_state[n];
    }
  }

  // b_education_state
  for (m in 1:NEducation) { for (n in 1:NState) {
      b_education_state[m,n] = 0;
      for (j in 1:NSex) {
        for (k in 1:NRace) {
          for (l in 1:NAge) {
              b_education_state[m,n] = b_education_state[m,n] + mtrx[j,k,l,m,n];
          }
        }
      }
      b_education_state[m,n] = b_education_state[m,n]/NSex/NRace/NAge - b_0 - b_education[m] - b_state[n];
    }
  }
}

