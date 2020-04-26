# 1 Model with metric output and one nominal predictor

# 1.1 Sex&Death example of section 19.3 in [K]
source("DBDA2Eprograms/DBDA2E-utilities.R")  # we will use function gammaShRaFromModeSD()

# Read the data, the file is available at [K].

mydf <- read.csv("DBDA2Eprograms/FruitflyDataReduced.csv")
head(mydf)

levels(mydf$CompanionNumber)

(grMeans <- aggregate(Longevity ~ CompanionNumber, mydf, mean))

dataList <- list(
  Ntotal = nrow(mydf),
  y = mydf$Longevity,
  x = as.integer(mydf$CompanionNumber),
  NxLvl = nlevels(mydf$CompanionNumber),
  agammaShRa = unlist(gammaShRaFromModeSD(
    mode = sd(mydf$Longevity) / 2,
    sd = 2 * sd(mydf$Longevity)
  ))
)


library(rstan)

# Define the model according to the diagram.

modelString<-"
data {
  int<lower=1> Ntotal;
  real y[Ntotal];
  int<lower=2> NxLvl;
  int<lower=1, upper=NxLvl> x[Ntotal];
  real<lower=0> agammaShRa[2];
}
transformed data {
  real meanY;
  real sdY;
  meanY = mean(y);
  sdY = sd(y);
}
parameters {
  real a0;
  real<lower=0> aSigma;
  vector[NxLvl] a;
  real<lower=0> ySigma;
}
model {
  a0 ~ normal(meanY, 5*sdY); // prior of the grand mean
  aSigma ~ gamma(agammaShRa[1], agammaShRa[2]);
  a ~ normal(0, aSigma);
  ySigma ~ uniform(sdY/100, sdY*10);
  for ( i in 1:Ntotal ) {
    y[i] ~ normal(a0 + a[x[i]], ySigma);
    }
}
generated quantities {
  // Convert a0,a[] to sum-to-zero b0,b[] :
  real b0;
  vector[NxLvl] b;
  b0 = a0 + mean(a); // mean of all the slopes is a
  b = a - mean(a);
}"

# Run MCMC
stanDso <- stan_model(model_code = modelString)

fit <- sampling(
  stanDso,
  data = dataList,
  pars = c('b0', 'b', 'aSigma', 'ySigma'),
  iter = 5000,
  chains = 2,
  cores = 2
)

# Analyze the chains using shinystan() (not shown).

library(shinystan)
launch_shinystan(fit)

summary(fit)$summary[, c(1, 4, 6, 8, 10)]

cbind(
  GroupMeans = grMeans,
  EstimatedMeans = summary(fit)$summary[2:6, 1] + summary(fit)$summary[1, 1]
)

#Note differences between estimated group means and also shrinkage.

stan_ac(fit, separate_chains = T)

pairs(fit)

# Significant difference betweent he groups
plot(fit)

plot(fit, pars = c("b"))

# No shrinkage; densities are smooth and nice
stan_dens(fit)

# Calculate contrasts for betas.

# Extract chains for parameters b.
fit_ext <- rstan::extract(fit)
head(fit_ext$b)

dim(fit_ext$b)
library(HDInterval)

# Contrast between group 1 ('None0') and group 5 ('Virgin8')
# Take chains difference and plot difference; minimum is greater than zero
contrast_1_5 <- fit_ext$b[, 1] - fit_ext$b[, 5]
plot(contrast_1_5)

hist(contrast_1_5)

(hdiContrast_1_5 <- hdi(contrast_1_5))

# See how wide distribution of contrast is:
(sd.contrast_1_5 <- sd(contrast_1_5))

plot(rank(fit_ext$b[, 1]), rank(fit_ext$b[, 5]))

# Contrast between average of groups 1-3 (None0, Pregnant1, Pregnant8) and average of groups 4-5 (Virgin1, Virgin8)
Comb1 <- (fit_ext$b[, 1] + fit_ext$b[, 2] + fit_ext$b[, 3]) / 3
Comb2 <- (fit_ext$b[, 4] + fit_ext$b[, 5]) / 2
contrast_123_45 <- Comb1 - Comb2
plot(contrast_123_45)

# Significantly different from zero; 4-5 are control groups
hist(contrast_123_45)

(hdiContrast_123_45 <- hdi(contrast_123_45))

(sd.contrast_123_45 <- sd(contrast_123_45))

head(cbind(Comb1, Comb2))

# We normalized; so the relationship is perfected since all groups compared will be offsets
plot(Comb1[1:100], Comb2[1:100])

plot(rank((fit_ext$b[, 1] + fit_ext$b[, 2] + fit_ext$b[, 3]) / 3), rank((fit_ext$b[, 4] + fit_ext$b[, 5]) / 2))

# 2 Adding metric predictor

# Create data with both predictors.

dataList2 <- list(
  Ntotal = nrow(mydf),
  y = mydf$Longevity,
  xMet = mydf$Thorax,
  xNom = as.integer(mydf$CompanionNumber),
  NxLvl = nlevels(mydf$CompanionNumber),
  agammaShRa = unlist(gammaShRaFromModeSD(
    mode = sd(mydf$Longevity) / 2,
    sd = 2 * sd(mydf$Longevity)
  ))
)

modelString<-"
data {
  int<lower=1> Ntotal;
  real y[Ntotal];
  int<lower=2> NxLvl;
  int<lower=1, upper=NxLvl> xNom[Ntotal];
  real xMet[Ntotal];
  real<lower=0> agammaShRa[2];
}
transformed data {
  real meanY;
  real sdY;
  real xMetMean;
  real xMetSD;
  meanY = mean(y);
  sdY = sd(y);
  xMetMean = mean(xMet);
  xMetSD = sd(xMet);
}
parameters {
  real a0;
  real<lower=0> aSigma;
  vector[NxLvl] a;
  real aMet;
  real<lower=0> ySigma;
}
model {
  a0 ~ normal(meanY, 5*sdY);
  aSigma ~ gamma(agammaShRa[1], agammaShRa[2]);
  a ~ normal(0, aSigma);
  aMet ~ normal(0, 2*sdY/xMetSD);
  ySigma ~ uniform(sdY/100, sdY*10);
  for ( i in 1:Ntotal ) {
    y[i] ~ normal(a0 + a[xNom[i]] + aMet*(xMet[i] - xMetMean), ySigma);
    }
}
generated quantities {
  // Convert a0,a[] to sum-to-zero b0,b[] :
  real b0;
  vector[NxLvl] b;
  b0 = a0 + mean(a) - aMet * xMetMean;
  b = a - mean(a);
}
"

# fit model
model2 <- stan_model(model_code = modelString)

#Run MCMC.
fit2 <- sampling(
  model2,
  data = dataList2,
  pars = c('b0', 'b', 'aMet', 'aSigma', 'ySigma'),
  iter = 5000,
  chains = 2,
  cores = 2
)

launch_shinystan(fit2)

# Calculate same contrasts for betas:
# Contrasts between group 1 ('None0') and 5 ('Virgin8').

fit_ext2 <- rstan::extract(fit2)
head(fit_ext2$b)

contrast2_1_5 <- fit_ext2$b[, 1] - fit_ext2$b[, 5]
plot(contrast2_1_5)

hist(contrast2_1_5)

(hdiContrast2_1_5 <- hdi(contrast2_1_5))

(sd.contrast2_1_5 <- sd(contrast2_1_5))

plot(rank(fit_ext2$b[, 1]), rank(fit_ext2$b[, 5]))

# Contrast between average of groups 1-3 and average of groups 4-5
Comb1 <- (fit_ext2$b[, 1] + fit_ext2$b[, 2] + fit_ext2$b[, 3]) / 3
Comb2 <- (fit_ext2$b[, 4] + fit_ext2$b[, 5]) / 2
contrast2_123_45 <- Comb1 - Comb2
plot(contrast2_123_45)

hist(contrast2_123_45)

(hdiContrast2_123_45 <- hdi(contrast2_123_45))

(sd.contrast2_123_45<-sd(contrast2_123_45))

plot(rank(Comb1),rank(Comb2))

# Compare the accuracy of contrast comparison.

# Using standard deviation:

rbind(
  WithoutMetricPred = c(sd.contrast_1_5, sd.contrast_123_45),
  WithMetricPred = c(sd.contrast2_1_5, sd.contrast2_123_45)
)

# Using HDI:

rbind(
  WithoutMetricPred = c(hdiContrast_1_5, hdiContrast_123_45),
  WithMetricPred = c(hdiContrast2_1_5, hdiContrast2_123_45)
)

# Using HDI width:
# Accuracy of the distribution of the parameter is narrower
rbind(
  WithoutMetricPred = c(
    Contrast_1_5 = unname(diff(hdiContrast_1_5)),
    Contrast_123_45 = unname(diff(hdiContrast_123_45))
  ),
  WithMetricPred = c(
    Contrast_1_5 = diff(hdiContrast2_1_5),
    Contrast_123_45 = diff(hdiContrast2_123_45)
  )
)

# 3 Heterogeneous variances

# Prepare the data.

mydf <- InsectSprays
head(mydf)

levels(mydf$spray)

dataListSprays <- list(
  Ntotal = nrow(mydf),
  y = mydf$count,
  x = as.integer(mydf$spray),
  NxLvl = nlevels(mydf$spray),
  aGammaShRa = unlist(gammaShRaFromModeSD(
    mode = sd(mydf$count) / 2,
    sd = 2 * sd(mydf$count)
  ))
)

# Plot the data.

plot(mydf$count ~ mydf$spray)

#Apply traditional ANOVA method.

m1 <- lm(count ~ spray, mydf)
summary(m1)

anova(m1)

# The categorical factor is significant, i.e. utility test fails.

# Create contrasts to identify differences between combinations of sprays.

A <- factor(mydf$spray)
contrasts(A)

contrasts(A) <- cbind(
  "GA vs GB" = c(1, -1, 0, 0, 0, 0),
  "GAB vs GF" = c(1, 1, 0, 0, 0, -2),
  "GC vs GD" = c(0, 0, 1, -1, 0, 0),
  "GC vs GE" = c(0, 0, 1, 0, -1, 0),
  "GABF vs GCDE" = c(1, 1, -1, -1, -1, 1)
)

A

aov(mydf$count ~ A)

# Epsilon has to have the same sigma everywhere...so out of all residuals there will be one sigma estimated
summary.lm(aov(mydf$count ~ A))

# However, separate t-test of C vs. D shows significant difference and equivalence of C and E is almost rejected by t-test.

# Linear model same epsilon causes groups to overlap

# Pairwise picks up on it

t.test(mydf$count[mydf$spray == "C"], mydf$count[mydf$spray == "D"])

t.test(mydf$count[mydf$spray == "C"], mydf$count[mydf$spray == "E"])

# Run MCMC.

modelString<-"
data {
  int<lower=1> Ntotal;
  real y[Ntotal];
  int<lower=2> NxLvl;
  int<lower=1, upper=NxLvl> x[Ntotal];
  real<lower=0> aGammaShRa[2];
}
  transformed data {
  real meanY;
  real sdY;
  meanY = mean(y);
  sdY = sd(y);
}
parameters {
  real<lower=0> nu;
  real a0;
  real<lower=0> aSigma;
  vector[NxLvl] a;
  real<lower=0> ySigma[NxLvl];
  real<lower=0> ySigmaMode;
  real<lower=0> ySigmaSD;
}
transformed parameters{
  real<lower=0> ySigmaSh;
  real<lower=0> ySigmaRa;
  ySigmaRa = ( ( ySigmaMode + sqrt( ySigmaMode^2 + 4*ySigmaSD^2 ) ) / ( 2*ySigmaSD^2 ) );
  ySigmaSh = 1 + ySigmaMode * ySigmaRa;
}
model {
  nu ~ exponential(1/30.0);
  a0 ~ normal(meanY, 10*sdY);
  aSigma ~ gamma(aGammaShRa[1], aGammaShRa[2]);
  a ~ normal(0, aSigma);
  ySigma ~ gamma(ySigmaSh, ySigmaRa);
  ySigmaMode ~ gamma(aGammaShRa[1], aGammaShRa[2]);
  ySigmaSD ~ gamma(aGammaShRa[1], aGammaShRa[2]);
  for ( i in 1:Ntotal ) {
    y[i] ~ student_t(nu, a0 + a[x[i]], ySigma[x[i]]);
    }
}
generated quantities {
  // Convert a0,a[] to sum-to-zero b0,b[] :
  real b0;
  vector[NxLvl] b;
  b0 = a0 + mean(a);
  b = a - mean(a);
}
"

model3 <- stan_model(model_code = modelString)

# Run chains.

fit.sprays <- sampling(
  model3,
  data = dataListSprays,
  pars = c('b0', 'b', 'aSigma', 'ySigma', 'nu', 'ySigmaMode', 'ySigmaSD'),
  iter = 5000,
  chains = 2,
  cores = 2
)

launch_shinystan(fit)

# Analyze estimates.

plot(fit.sprays, pars = c("b0", "aSigma"))

plot(fit.sprays, pars = c("b"))

plot(fit.sprays, pars = c("ySigma"))

# Extract chains.

chains_sprays <- rstan::extract(fit.sprays)
head(chains_sprays$b)

#Look at the same contrasts.

#A vs. B

contrast1_2 <- chains_sprays$b[, 1] - chains_sprays$b[, 2]
plot(contrast1_2)

hist(contrast1_2)

(hdiContrast1_2 <- hdi(contrast1_2))

(sd.contrast1_2 <- sd(contrast1_2))

# Practically independent
plot(rank(chains_sprays$b[, 1]), rank(chains_sprays$b[, 5]))

#The contrast is not different from zero.

#A, B vs. F
Comb1 <- chains_sprays$b[, 1] + chains_sprays$b[, 2]
Comb2 <- 2 * chains_sprays$b[, 6]
contrast12_6 <- Comb1  - Comb2
plot(contrast12_6)

hist(contrast12_6)


(hdiContrast12_6 <- hdi(contrast12_6))

(sd.contrast12_6 <- sd(contrast12_6))

plot(rank(Comb1), rank(Comb2))

# The contrast is not different from zero.

#C vs. D
contrast3_4 <- chains_sprays$b[, 3] - chains_sprays$b[, 4]
plot(contrast3_4)

hist(contrast3_4)

(hdiContrast3_4 <- hdi(contrast3_4))

(sd.contrast3_4<-sd(contrast3_4))

plot(rank(chains_sprays$b[,3]),rank(chains_sprays$b[,4]))

#This contrast is different from zero. ANOVA could not detect that.

#C vs. E
contrast3_5 <- chains_sprays$b[, 3] - chains_sprays$b[, 5]
plot(contrast3_5)

hist(contrast3_5)

(hdiContrast3_5 <- hdi(contrast3_5))
(sd.contrast3_5 <- sd(contrast3_5))

plot(rank(chains_sprays$b[, 3]), rank(chains_sprays$b[, 5]))

# The contrast is not different from zero.

# Combined observations A, B and F vs. C, D, E
Comb1 <- chains_sprays$b[, 1] + chains_sprays$b[, 2] + chains_sprays$b[, 6]
Comb2 <- chains_sprays$b[, 3] + chains_sprays$b[, 4] + chains_sprays$b[, 5]
contrast126_345 <- Comb1  - Comb2
plot(contrast126_345)

hist(contrast126_345)

(hdiContrast126_345<-hdi(contrast126_345))
(sd.contrast126_345<-sd(contrast126_345))
plot(rank(Comb1), rank(Comb2))

# The contrast is significantly different from zero.
