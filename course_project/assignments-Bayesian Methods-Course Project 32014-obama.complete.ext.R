sampl <- read.csv('./data.csv')

library (rstan)
rstan_options(auto_write = TRUE)

t1 <- Sys.time()
message('Started fitting at ', t1)

model <- stan_model(file="./obama.complete.ext.stan")
t2 <- Sys.time()
message('Model compiled at ', t2)

obama_fit <- sampling(model, 
                      data=list(N = length(sampl$y),
                                y = sampl$y,
                                sex = as.integer(sampl$sex),
                                NSex = nlevels(sampl$sex),
                                race = as.integer(sampl$race),
                                NRace = nlevels(sampl$race),
                                age = as.integer(sampl$age),
                                NAge = nlevels(sampl$age),
                                education = as.integer(sampl$education),
                                NEducation = nlevels(sampl$education),
                                state = as.integer(sampl$state),
                                NState = nlevels(sampl$state) 
                      ),
                      pars=c('b_0','b_sex', 'b_race', 'b_age', 'b_education', 'b_state',
                             'b_sex_race', 'b_sex_age', 'b_sex_education', 'b_sex_state',
                             'b_race_age', 'b_race_education', 'b_race_state',
                             'b_age_education', 'b_age_state', 'b_education_state',
                             'var_0', 'var_sex', 'var_race', 'var_age', 'var_education', 'var_state',
                             'var_sex_race', 'var_sex_age', 'var_sex_education', 'var_sex_state',
                             'var_race_age', 'var_race_education', 'var_race_state',
                             'var_age_education', 'var_age_state', 'var_education_state',
                             'nu', 'sigma'
                      ),
                      control=list(adapt_delta=0.99, max_treedepth=12),
                      iter=1000, chains = 4, cores = 4,
                      verbose = F)

t3 <- Sys.time()
message('Finished fitting at ', t3)
message('Time elapsed: ', difftime(t3,t1, units = 'h'), ' hours')

# save fitted model:
file_suffix <- strftime(t2, format='%Y%m%d_%H%M%S')
fn <- paste0('./fit_ext_', file_suffix, '.Rdata')
save(obama_fit, file=fn)
message('Saved model to ', fn, '. Goodbye!')


#print (obama_fit)
#plot (obama_fit, pars=c('b_sex_race'))
#pairs(obama_fit, pars=c('b_0','b_sex', 'b_race'))

#library(shinystan)
#launch_shinystan(obama_fit)

#str(obama_fit)
#stan_plot(obama_fit)





