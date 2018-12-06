#1 Ordinal Response

#1.1 Single Group
#Consider the simplest case when latent variable is standard normal does not have predictors, threshold parameters have uniform prior.
#There is only one ordinal variable.


modelString<-"
data {
int<lower=2> K;  # num of y classes
int<lower=0> N;  # num of observations
int<lower=1,upper=K> y[N];
}
parameters {
ordered[K-1] c;
}
model {
vector[K] theta;
for (n in 1:N) {
theta[1] = Phi(c[1]);
for (k in 2:(K-1))
theta[k] = Phi(c[k]) - Phi(c[k-1]);
theta[K] = 1 - Phi(c[K-1]);
y[n] ~ categorical(theta);
}
}
"

#Read the data.
#The sample is just one column of ratings between and 7.
#Calculate frequencies of different levels of the rating.

myData = read.csv( file=paste(dataPath,"OrdinalProbitData-1grp-1.csv",sep="/") )
head(myData)

dim(myData)


table(myData)

(frequences<-table(myData)/dim(myData)[1])

#Compile the model DSO and run Markov Chains.

library(rstan)
# fit model
model <- stan_model(model_code=modelString)

fit <- sampling(model,
                data=list(N=nrow(myData),  # num of observations
                          K=max(myData$Y), # num of outcome classes
                          y=myData$Y),
                pars=c('c'),
                iter=5000, chains = 2, cores = 2
)

# analyze fitted model using shinystan
library(shinystan)
launch_shinystan(fit)

#Check the chains.

pairs(fit)

(fitResults<-summary(fit)$summary[,c(1,4,6,8,10)])

#Returned parameters are 6 thresholds separating 7 ordinal categories.

stan_ac(fit, separate_chains = T)

stan_trace(fit)

stan_dens(fit)

#Compare estimated means with the frequencies of the sample

cbind(CumulativeFrequences=head(cumsum(frequences)),EstimatedProbabilities=pnorm(head(fitResults[,1])))

#Plot HDI of the parameters.

plot(fit)

#HDI intervals overlap which may seem counter-intuitive: thresholds have to be ordered.

#Extract the chains.

fit_ext <- rstan::extract(fit)
fit_ext<-fit_ext$c
head(fit_ext)

#Plot first 6 rows of the chains to see that at each step of the chain thresholds are ordered.

plot(fit_ext[1,],rep(1,6),ylim=c(1,6),xlim=c(0,3),col=1,pch=16)
points(fit_ext[2,],rep(2,6),col=2,pch=16)
points(fit_ext[3,],rep(3,6),col=3,pch=16)
points(fit_ext[4,],rep(4,6),col=4,pch=16)
points(fit_ext[5,],rep(5,6),col=5,pch=16)
points(fit_ext[6,],rep(6,6),col=6,pch=16)

#1.2 Two Groups

#Two ordinal variables may need to be compared, for example, in the following situations:
#  1. Two groups of participants asked to answer a questionnaire with the same scale of ordinal categories, for example, "Strongly Disagree", "Disagree", "Undecided", "Agree", "Strongly Agree". But two groups are asked about their agreement/disagreement with two different statements about social issues: "Left-handed people should have equal rights under the law" - for group 1, and "Disabled people should be given equal rights under the law".
#The assumption is that people in both groups have the same latent variable, call it "sense of fairness". But two social issues have different responses formalized in mean value and standard deviation of the latent variable with respect to the issues.
#The goal of the experiment may be analysis of contrast between the two statements, i.e. how significant the difference is between the two statements.

#Two groups of participants are asked to rate two different products and select level of satisfaction from 1 (lowest) to 10 (highest). Again the scale of ordinal responses is the same. Assumption is that the latent variable of evaluation of each product, "satisfaction", is shared by both groups of participants. But different products result in different mean values and variances of the latent variable. The goal of the study is understanding the difference between the products.


modelString<-"
data {
int<lower=2> K;  // num of y classes
int<lower=0> N;  // num of observations
int<lower=1> D;  // num of groups
int<lower=1,upper=K> y[N]; // response
int<lower=1,upper=D> x[N]; // group index: each obeservation has its own group
}
parameters {
ordered[K-3] c_raw;
vector[D] sigma; // group sigma
vector[D] beta;  // group mean
}
transformed parameters{ // renormalize vector c
vector[K-1] c;
c[1] = 1.5;
for (k in 2:(K-2))
c[k] = 1.5 + (K - 2.0)* inv_logit(c_raw[k-1]); # sigmoid
c[K-1] = K - 0.5;
}
model {
vector[K] theta;
real mu;
real sgm;
beta ~ normal((1+K)/2.0, K);
sigma ~ uniform(K/1000.0, K * 10.0);
for (n in 1:N) {
mu = beta[x[n]];
sgm = sigma[x[n]];
theta[1] = Phi((c[1] - mu)/sgm);
for (k in 2:(K-1))
theta[k] = Phi((c[k] - mu)/sgm) - Phi((c[k-1] - mu)/sgm);
theta[K] = 1 - Phi((c[K-1] - mu)/sgm);
y[n] ~ categorical(theta);
}
}"

#Create DSO.

model <- stan_model(model_code=modelString)

#Read the data.

myData = read.csv( file=paste(dataPath,"OrdinalProbitData1.csv",sep="/") )
head(myData)

table(myData)

#The dataset now has 2 groups A and B, showing different distributions on the 5 ordinal categories.

#Run the chains.

fit <- sampling(model,
                data=list(N=nrow(myData),  # num of observations
                          K=max(myData$Y), # num of outcome classes
                          D=nlevels(myData$X),
                          y=myData$Y,
                          x=as.integer(myData$X)),
                pars=c('c', 'beta', 'sigma'),
                iter=5000, chains = 2, cores = 2
)

#Analyze fitted model using shinystan

library(shinystan)
launch_shinystan(fit)

pairs(fit)

stan_ac(fit, separate_chains = T)

stan_trace(fit)


stan_dens(fit)

#Densities show very similar shapes of posterior distributions for beta[1],beta[2] and for sigma[1],sigma[2].

summary(fit)
plot(fit,pars=c("beta"))


#HDIs show that even at 80% the difference between the two means is insignificant.

plot(fit,pars=c("sigma"))

#Check whether t-test is able to distinguish the means.

subsetA<-subset(myData,X=="A")
subsetB<-subset(myData,X=="B")
t.test(subsetA$Y,subsetB$Y)

#With any level greater than 0.032 the equality hypothesis is rejected.

#In the book this is explained by non-normality of the data.
#Confirm this with qq-plot.

qqnorm(subsetA$Y)
qqline(subsetA$Y)


qqnorm(subsetB$Y)
qqline(subsetB$Y)

#Distributions are indeed not normal. But, probably, the main reason for rejecting the null hypothesis is discrete format of the data.

#Compare the samples assuming that they have multinomial distribution.
#Then apply chi-square test.

chisq.test(c(table(subsetA$Y)/length(subsetA$Y),0),table(subsetB$Y)/length(subsetB$Y))

#With proper method selected the hypothesis of equality of means is not rejected.

#Plot row means against simulated thresholds.

draws <- extract(fit)
c1 <- draws$c[,1]
c2 <- draws$c[,2]
c3 <- draws$c[,3]
c4 <- draws$c[,4]
mc <- rowMeans(draws$c)
head(cbind(c1,c4))

head(mc)


plot(c1, mc, xlim = c(1,5),xlab="c")
points(c2,mc)
points(c3,mc)
points(c4,mc)
abline(v=mean(c2))
abline(v=mean(c3))


#2 Probit Model with Metric Predictors
#2.1 Single predictor: Happiness and Money

#Read the data.

# Example #1: one metric predictor (predict happiness with money)
myData = read.csv( file=paste(dataPath,"HappinessAssetsDebt.csv",sep="/") )
head(myData)

table(myData$Happiness)


#Note that majority of responses are "4".

modelStringHappiness<-"
data {
int<lower=2> K;  // num of y classes
int<lower=0> N;  // num of observations
int<lower=1> D;  // num of metric predictors
int<lower=1,upper=K> y[N];
matrix[N, D] x;
}
transformed data {
vector[D] x_m;  // x means
vector[D] x_sd; // x standard deviations
matrix[N, D] zx;    // normalized x
for (j in 1:D) {
x_m[j] = mean(x[,j]);
x_sd[j] = sd(x[,j]);
zx[,j] = (x[,j] - x_m[j]) / x_sd[j];
}
}
parameters {
ordered[K-1] c;
vector[D] zbeta;
}
model {
vector[K] theta;
real eta;
//zbeta ~ normal(0, K); // ?????????
for (n in 1:N) {
eta = dot_product(zbeta, zx[n,]);
theta[1] = Phi(c[1] - eta);
for (k in 2:(K-1))
theta[k] = Phi(c[k] - eta) - Phi(c[k-1] - eta);
theta[K] = 1 - Phi(c[K-1] - eta);
y[n] ~ categorical(theta);
//y[n] ~ ordered_logistic(dot_product(zbeta, zx[n,]), c);
}
}
generated quantities{
vector[D] beta;
real beta0;
beta = zbeta ./ x_sd;  // < ./ > is elementwise division
beta0 = - sum( zbeta .* x_m ./ x_sd );
}
"

#The model returns scaled variables.

#Run the chains.

# compile model:
model <- stan_model(model_code=modelStringHappiness)

fit <- sampling(model,
                data=list(N=nrow(myData),  # num of observations
                          K=max(myData$Happiness), # num of outcome classes
                          D=1,
                          y=myData$Happiness,
                          x=cbind(x=myData$Assets)),
                pars=c('c', 'beta', 'beta0'),
                iter=5000, chains = 2, cores = 2
)

#Check shiny application.

# analyze fitted model using shinystan
library(shinystan)
launch_shinystan(fit)

#Check the chains and parameters.

summary(fit)$summary[,c(1,4,6,8,10)
                     
pairs(fit)

stan_ac(fit, separate_chains = T)

stan_trace(fit)

stan_dens(fit,pars=c("beta0"))

stan_dens(fit,pars=c("beta"))

plot(fit,pars=c("c"))

plot(fit,pars=c("beta0"))

plot(fit,pars=c("beta"))

#Plot the thresholds.

draws <- extract(fit)
c1 <- draws$c[,1]
c2 <- draws$c[,2]
c3 <- draws$c[,3]
c4 <- draws$c[,4]
mc <- rowMeans(draws$c)
head(cbind(c1,c4))

head(mc)

range(c1)

plot(c1, mc, xlim = c(-3,3),xlab="c")
points(c2,mc)
points(c3,mc)
points(c4,mc)
abline(v=mean(c1))
abline(v=mean(c2))
abline(v=mean(c3))
abline(v=mean(c4))


#Note that range for category 4 is wider than other ranges.
#Why?
library(HDInterval)

#Check HDI to see how significant is the effect of wealth on happiness.

hdi(draws$beta)