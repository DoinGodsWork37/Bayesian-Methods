library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit_robust_adjust <- sampling (
  stanDsoRobustRegPanel,
  data = dataList,
  pars = c(
    "nu",
    "sigma",
    "beta0mu",
    "beta1mu",
    "beta0",
    "beta1",
    "zbeta0sigma",
    "zbeta1sigma"
  ),
  iter = 50000,
  chains = 4,
  cores = 4,
  control = list(
    adapt_delta = 0.999,
    stepsize = 0.01,
    max_treedepth = 15
  )
)
