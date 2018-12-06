# 1 2-Way ANOVA in Bayesian Setup
# The resulting model has the following diagram.

library(rstan)
# we are going to use func gammaShRaFromModeSD()
source("DBDA2Eprograms/DBDA2E-utilities.R")

# 1.1 Model description
# Define the model according to the diagram.

modelString<-"
data {
  int<lower=1> Ntotal;
  vector[Ntotal] y;
  int<lower=2> Nx1Lvl;
  int<lower=2> Nx2Lvl;
  int<lower=1, upper=Nx1Lvl> x1[Ntotal];
  int<lower=1, upper=Nx2Lvl> x2[Ntotal];
  real<lower=0> agammaShRa[2];
}
  transformed data {
  real meanY;
  real sdY;
  vector[Ntotal] zy;
  meanY = mean(y);
  sdY = sd(y);
  zy = (y - mean(y)) / sdY;  // center & normalize
}
parameters {
  real a0;
  real<lower=0> a1Sigma; // 1st factor
  real<lower=0> a2Sigma; // 2nd factor
  real<lower=0> a1a2Sigma; // interaction
  vector[Nx1Lvl] a1;
  vector[Nx2Lvl] a2;
  matrix[Nx1Lvl,Nx2Lvl] a1a2;
  real<lower=0> zySigma;
}
model {
  a0 ~ normal(0, 1);
  a1Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
  a1 ~ normal(0, a1Sigma);
  a2Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
  a2 ~ normal(0, a2Sigma);
  a1a2Sigma ~ gamma(agammaShRa[1], agammaShRa[2]);
  for (j1 in 1:Nx1Lvl) {
    a1a2[j1,] ~ normal(0, a1a2Sigma);
    }
  zySigma ~ uniform(1.0/10, 10);
  for ( i in 1:Ntotal ) {
  zy[i] ~ normal(a0 + a1[x1[i]] + a2[x2[i]]+ a1a2[x1[i],x2[i]], zySigma);
  }
}
generated quantities {
  // Convert a to sum-to-zero b :
  real b0;
  vector[Nx1Lvl] b1;
  vector[Nx2Lvl] b2;
  matrix[Nx1Lvl,Nx2Lvl] b1b2;
  matrix[Nx1Lvl,Nx2Lvl] m;
  real<lower=0> b1Sigma;
  real<lower=0> b2Sigma;
  real<lower=0> b1b2Sigma;
  real<lower=0> ySigma;
  for ( j1 in 1:Nx1Lvl ) {
    for ( j2 in 1:Nx2Lvl ) {
    m[j1,j2] = a0 + a1[j1] + a2[j2] + a1a2[j1,j2]; // cell means
    }
  }
  b0 = mean(m);
  for ( j1 in 1:Nx1Lvl ) {
    b1[j1] = mean( m[j1,] ) - b0;
  }
  for ( j2 in 1:Nx2Lvl ) {
    b2[j2] = mean( m[,j2] ) - b0;
  }
  for ( j1 in 1:Nx1Lvl ) {
    for ( j2 in 1:Nx2Lvl ) {
      b1b2[j1,j2] = m[j1,j2] - ( b0 + b1[j1] + b2[j2] );
    }
  }
  // transform to original scale:
  b0 = meanY + sdY * b0;
  b1 = sdY * b1;
  b2 = sdY * b2;
  b1b2 = sdY * b1b2;
  b1Sigma = sdY * a1Sigma;
  b2Sigma = sdY * a2Sigma;
  b1b2Sigma = sdY * a1a2Sigma;
  ySigma = sdY * zySigma;
}"

# Create DSO.

stanDsoANOVA2Way <- stan_model(model_code = modelString)

# 1.2 Data
# Read the data.

# The data show faculty salaries of a university.

# The columns selected as predictors are:

# Position (Assistant Professor, Associate Professor, Full Professor, Full Professor with endowment salary and Distinguished Professor)
# Department (total 60 departments)
# 20_1: Metric Predicted Variable with Two Nominal Predictors

# load data from 'Salary.csv' (see Kruschke)
mydf = read.csv("DBDA2Eprograms/Salary.csv")
mean(mydf$Salary)

head(mydf)

dim(mydf)

colnames(mydf)

table(mydf$Pos)

table(mydf$Org)

length(table(mydf$Org))

y <- mydf$Salary;
x1 <- mydf$Pos;
x2 <- mydf$Org;
dataListSalary <- list(
  Ntotal = length(y),
  y = y,
  x1 = as.integer(x1),
  x2 = as.integer(x2),
  Nx1Lvl = nlevels(x1),
  Nx2Lvl = nlevels(x2),
  agammaShRa = unlist(gammaShRaFromModeSD(mode = 1 /
                                            2, sd = 2))
)

# Create names of variables and their interactions for further reference.

(namesPos <- names(table(mydf$Pos)))

(namesOrg <- names(table(mydf$Org)))

# Use outer() to create names for interactions.

as.vector(outer(1:4, 1:2, paste, sep = "-"))

#Interactions names:

(namesInter <- as.vector(outer(namesOrg, namesPos, paste, sep = "-")))

# Finally, all names:

varNames <- c("Intercept", namesPos, namesOrg, namesInter, rep("Var", 5))

# 1.3 MCMC
# Run MCMC.

# fit model
fit <- sampling (
  stanDsoANOVA2Way,
  data = dataListSalary,
  pars = c(
    'b0',
    'b1',
    'b2',
    'b1b2',
    'b1Sigma',
    'b2Sigma',
    'b1b2Sigma',
    'ySigma'
  ),
  iter = 5000,
  chains = 2,
  cores = 2
)

#Check the results in shinystan (not shown).

library(shinystan)
launch_shinystan(fit)

#Create results including mean value, 2.5%, 50% and 97.5% quantiles.
#Add variable names as row names.

SalaryResults <- summary(fit)$summary[, c(1, 4, 6, 8)]
varNames[nrow(SalaryResults) - (4:0)] <- rownames(SalaryResults)[nrow(SalaryResults) - (4:0)]
rownames(SalaryResults) <- varNames
head(SalaryResults)

# Make plots of mean values and HDIs.

plot(fit, pars = c('b1'))

plot(fit, pars = c('b2'))

plot(fit, pars = c('b1b2'))

# Plots show that not all coefficients of the model equal zero.
# This answers the question of utility test.

# 1.4 Working with chains and contrasts
# Extract chains for the position variables.

fit_ext <- rstan::extract(fit)
names(fit_ext)

fit_ext.b0 <- fit_ext$b0
fit_ext.b1 <- fit_ext$b1
colnames(fit_ext.b1) <- namesPos
head(fit_ext.b1)

# Extract chains for the department variables.

fit_ext.b2 <- fit_ext$b2
colnames(fit_ext.b2) <- namesOrg
head(fit_ext.b2)

# Extract chains for interaction variables.

fit_ext.b1.b2 <- fit_ext$b1b2
dim(fit_ext.b1.b2)

# Interaction chains make a cube.

dimnames(fit_ext.b1.b2)[[2]] <- namesPos
dimnames(fit_ext.b1.b2)[[3]] <- namesOrg
dimnames(fit_ext.b1.b2)

fit_ext.b1.b2[1, , ]

library(HDInterval)
# 1.4.1 Contrast for comparison of departments
# Use contrasts to compare salaries at Business and Finance ("BFIN") with Physics ("PHYS") and with Chemistry ("CHEM") departments.

# To do that select columns of MCMC for departments to "BFIN" and "PHYS" and "BFIN" and "CHEM", take their differences and look at the posterior distribution of the differences
contrast_bfin_phys <- fit_ext.b2[, "BFIN"] - fit_ext.b2[, "PHYS"]
contrast_bfin_chem <- fit_ext.b2[, "BFIN"] - fit_ext.b2[, "CHEM"]
contrast_phys_chem <- fit_ext.b2[, "PHYS"] - fit_ext.b2[, "CHEM"]

# What do you conclude about the differences of salaries?
# BFIN makes significantly more than both!

# Plot histograms of salaries for both departments.
hist(contrast_bfin_phys)
hist(contrast_bfin_chem)
hist(contrast_phys_chem)

# Salaries BFINComp of Business and Finance and PHYSComp of Physics are calculate as described above.
BFINComp <- fit_ext.b2[, "BFIN"] + fit_ext.b0
PHYSComp <- fit_ext.b2[, "PHYS"] + fit_ext.b0
CHEMComp <- fit_ext.b2[, "CHEM"] + fit_ext.b0

distrBFIN <- density(BFINComp)
distrPHYS <- density(PHYSComp)
distrCHEM <- density(CHEMComp)

plot(
  distrPHYS,
  xlim = c(110000, 250000),
  lwd = 2,
  main = "Salaries Distributions",
  col = "blue")
lines(distrBFIN$x, distrPHYS$y,lwd = 2, col = "orange")
lines(distrBFIN$x, distrCHEM$y,lwd = 2, col = "green")
legend(
  "top",
  legend = c("PHYS", "BFIN", "CHEM"),
  col = c("blue", "orange", "green"),
  lty = 1,
  lwd = 2
)

# Do the same comparison for finance professors and chemistry professors.


# 1.4.2 Contrast for comparison of positions
# Use contrasts to compare salaries of Endowment Full Professor ("NDW") and Distinguished Full Professor ("DST"). Compare salaries of Full Professor ("FT1") and Endowment Full Professor ("NDW")

contrast_ndw_dst <- fit_ext.b1[, "NDW"] - fit_ext.b1[, "DST"]
contrast_ft1_ndw <- fit_ext.b1[, "FT1"] - fit_ext.b1[, "NDW"]

# What do you conclude about the differences of salaries?

# Plot histograms of salaries for both departments.
hist(contrast_ndw_dst)
hist(contrast_ft1_ndw)

# Salaries are calculate as described above.
DSTComp <- fit_ext.b1[, "DST"] + fit_ext.b0
NDWComp <- fit_ext.b1[, "NDW"] + fit_ext.b0

distrDST <- density(DSTComp)
distrNDW <- density(NDWComp)
plot(
  distrDST,
  xlim = c(140000, 200000),
  lwd = 2,
  main = "Salaries Distributions",
  col = "blue"
)
lines(distrNDW$x, distrPHYS$y, lwd = 2, col = "orange")
legend(
  "top",
  legend = c("DST", "NDW"),
  col = c("blue", "orange"),
  lty = 1,
  lwd = 2
)

# 1.4.3 Contrast for comparison of spreads
# Use contrasts to compare salaries spreads between Full Professor and Assistant Professor at Physics Department and at Chemistry Department.

FT1_PHYS_Comp <- fit_ext.b0 + fit_ext.b1[, "FT1"] + fit_ext.b2[, "PHYS"] + fit_ext.b1.b2[, "FT1", "PHYS"]
FT3_PHYS_Comp <- fit_ext.b0 + fit_ext.b1[, "FT3"] + fit_ext.b2[, "PHYS"] + fit_ext.b1.b2[, "FT3", "PHYS"]
contrast_phys <- FT1_PHYS_Comp - FT3_PHYS_Comp

FT1_CHEM_Comp <- fit_ext.b0 + fit_ext.b1[, "FT1"] + fit_ext.b2[, "CHEM"] + fit_ext.b1.b2[, "FT1", "CHEM"]
FT3_CHEM_Comp <- fit_ext.b0 + fit_ext.b1[, "FT3"] + fit_ext.b2[, "CHEM"] + fit_ext.b1.b2[, "FT3", "CHEM"]
contrast_chem <- FT1_CHEM_Comp - FT3_CHEM_Comp

contrast_interact <- FT1_PHYS_Comp - FT3_CHEM_Comp

spreads <- contrast_phys - contrast_chem

hist(contrast_phys)
hist(contrast_chem)

distrFT1_PHS <- density(FT1_PHS_Comp)
distrFT3_CHEM <- density(FT3_CHEM_Comp)

plot(
  distrFT1_PHS,
  lwd = 2,
  main = "Salaries Distributions",
  col = "blue"
)
lines(distrFT1_PHS$x, distrFT3_CHEM$y, lwd = 2, col = "orange")
legend(
  "right",
  legend = c("FT1_PHS", "FT3_CHEM"),
  col = c("blue", "orange"),
  lty = 1,
  lwd = 2
)

# Find the highest HDI level for which the spread of differences between "FT1" and "FT3" is significant.
hdi(spreads, .80)
hdi(contrast_phys)
hdi(contrast_chem)

# 2 Understanding the effect of scaling and transformations on interactions
# Nonlinear transformations may affect interactions very significantly.

# Illustrate it on a simple simulated example.

mean00 <- 1
mean10 <- 3
mean01 <- 4
mean11 <- 6
y00 <- rnorm(5, mean00, .1)
y10 <- rnorm(5, mean10, .1)
y01 <- rnorm(5, mean01, .1)
y11 <- rnorm(5, mean11, .1)

# Plot the effects. If the lines are parallel the effects are additive.

plot(
  c(0, 1),
  c(mean(y00), mean(y10)),
  type = "b",
  ylim = c(1, 8),
  col = "darkgreen",
  lwd = 3,
  ylab = "Response",
  xlab = "Predictor 1"
)
lines(c(0, 1),
      c(mean(y01), mean(y11)),
      type = "b",
      col = "lightblue",
      lwd = 3)
legend(
  "topleft",
  legend = c("Predictor2 at 0", "Predictor2 at 1"),
  lty = 1,
  lwd = 3,
  col = c("darkgreen", "lightblue")
)

# Taking exponent of the same data introduces significant interaction.

plot(
  c(0, 1),
  c(mean(exp(y00)), mean(exp(y10))),
  type = "b",
  ylim = c(1, 400),
  col = "darkgreen",
  lwd = 3,
  ylab = "Response",
  xlab = "Predictor 1"
)
lines(c(0, 1),
      c(mean(exp(y01)), mean(exp(y11))),
      type = "b",
      col = "lightblue",
      lwd = 3)
legend(
  "topleft",
  legend = c("Predictor2 at 0", "Predictor2 at 1"),
  lty = 1,
  lwd = 3,
  col = c("darkgreen", "lightblue")
)
