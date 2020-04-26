Workshop 1
================
Josh Goldberg
October, 02 2018

# Joint Distribution: 3-Dimensional Case

``` r
u <- c("u1", "u2")
v <- c("v1", "v2", "v3")
w <- c("w1", "w2", "w3")
matr.u0 <- paste("u0", outer(v, w, paste, sep = ","), sep = ",")
dim(matr.u0) <- c(3, 3)
matr.u1 <- paste("u1", outer(v, w, paste, sep = ","), sep = ",")
dim(matr.u1) <- c(3, 3)
matr.u0
```

    ##      [,1]       [,2]       [,3]      
    ## [1,] "u0,v1,w1" "u0,v1,w2" "u0,v1,w3"
    ## [2,] "u0,v2,w1" "u0,v2,w2" "u0,v2,w3"
    ## [3,] "u0,v3,w1" "u0,v3,w2" "u0,v3,w3"

``` r
data3way <-
  read_delim("documents-MScA Bayesian Methods (32014)-MScA 32014 Lecture 1-3WayData.csv",
             delim = ",") %>%
  data.frame()
```

``` r
head(data3way)
```

    ##   w v u
    ## 1 2 1 1
    ## 2 2 3 0
    ## 3 3 1 1
    ## 4 1 1 1
    ## 5 2 1 1
    ## 6 1 1 1

``` r
mat.u0 <- table(subset(data3way, u == 0)[, 1], subset(data3way, u == 0)[, 2])
mat.u1 <- table(subset(data3way, u == 1)[, 1], subset(data3way, u == 1)[, 2])
```

``` r
mat.u0
```

    ##    
    ##     1 2 3
    ##   1 9 7 1
    ##   2 8 2 1
    ##   3 4 2 0

``` r
idx.v1 <- data3way$v == 1
idx.w1 <- data3way$w == 1
idx.u1 <- data3way$u == 1
sum(idx.v1 * idx.w1 * idx.u1) # element (1, 1) of mat.u1
```

    ## [1] 23

``` r
idx.v2 <- data3way$v == 2
sum(idx.v2 * idx.w1 * idx.u1) # element (1, 2) of mat.u1
```

    ## [1] 9

``` r
colnames(mat.u1) <- colnames(mat.u0) <- c("v1", "v2", "v3")
rownames(mat.u1) <- rownames(mat.u0) <- c("w1", "w2", "w3")

data3way.array <-
  array(rep(NA, 18),
        dim = c(3, 3, 2),
        dimnames = list(
          paste("w", 1:3, sep = ""),
          paste("v", 1:3, sep =
                  ""),
          paste("u", 0:1, sep =
                  "")
        ))

data3way.array[, , 1] <- mat.u0
data3way.array[, , 2] <- mat.u1
data3way.array
```

    ## , , u0
    ## 
    ##    v1 v2 v3
    ## w1  9  7  1
    ## w2  8  2  1
    ## w3  4  2  0
    ## 
    ## , , u1
    ## 
    ##    v1 v2 v3
    ## w1 23  9  1
    ## w2 12  8  1
    ## w3  4  7  1

``` r
N <- sum(data3way.array)
(data3way.array.p <- data3way.array / N)
```

    ## , , u0
    ## 
    ##      v1   v2   v3
    ## w1 0.09 0.07 0.01
    ## w2 0.08 0.02 0.01
    ## w3 0.04 0.02 0.00
    ## 
    ## , , u1
    ## 
    ##      v1   v2   v3
    ## w1 0.23 0.09 0.01
    ## w2 0.12 0.08 0.01
    ## w3 0.04 0.07 0.01

## 1.2 Marginal Distributions

Create marginal distribution for `u` as vector `uMarginal`.

``` r
(uMarginal <- apply(data3way.array.p, 3, sum))
```

    ##   u0   u1 
    ## 0.34 0.66

Create marginal distribution for `v` as vector `vMarginal`.

``` r
(vMarginal <- apply(data3way.array.p, 2, sum))
```

    ##   v1   v2   v3 
    ## 0.60 0.35 0.05

``` r
vMarginal["v1"]
```

    ##  v1 
    ## 0.6

Create marginal distribution for `w` as vector `wMarginal`.

``` r
(wMarginal <- apply(data3way.array.p, 1, sum))
```

    ##   w1   w2   w3 
    ## 0.50 0.32 0.18

``` r
wMarginal["w1"]
```

    ##  w1 
    ## 0.5

## 1.3 Conditional Distributions

Create conditional distribution \(p(w,v|u=1)\) as matrix
`cond.v.w.given.u1`.

``` r
cond.v.w.given.u1 <- data3way.array.p[, , "u1"] / sum(data3way.array.p[, , "u1"])
cond.v.w.given.u1["w1", ]
```

    ##         v1         v2         v3 
    ## 0.34848485 0.13636364 0.01515152

``` r
sum(cond.v.w.given.u1) # check
```

    ## [1] 1

Create conditional distribution \(p(v|u=1)\) as vector
`cond.v.given.u1`.

``` r
cond.v.given.u1 <- apply(data3way.array.p[, , "u1"], 2, sum) / sum(data3way.array.p[, , "u1"])
cond.v.given.u1["v1"]
```

    ##        v1 
    ## 0.5909091

``` r
sum(cond.v.given.u1) # check
```

    ## [1] 1

Create conditional distribution \(p(w|v=2,u=1)\) as vector
`cond.w.given.u1.v2`.

``` r
cond.w.given.u1.v2 <- data3way.array.p[, "v2", "u1"] / sum(data3way.array.p[, "v2", "u1"])
cond.w.given.u1.v2["w1"]
```

    ##    w1 
    ## 0.375

``` r
sum(cond.w.given.u1.v2) # check
```

    ## [1] 1

Compare the vectors \(p(w|v2,u1)p(v2|u1)p(u1)\) and
\(p(w,v,u)[,v2,u1]\)

``` r
rbind(uMarginal["u1"]*cond.v.given.u1["v2"]*cond.w.given.u1.v2,data3way.array.p[,"v2","u1"])
```

    ##        w1   w2   w3
    ## [1,] 0.09 0.08 0.07
    ## [2,] 0.09 0.08 0.07

# 2 Simulation Using Conditional Distributions

Let the marginal distribution for random variable `u` be Bernoulli with
\(p(u=0)=0.55\), \(p(u=1)=0.45\). Let conditional distributions for
random variables \((v|u=0)\) and \((v|u=1)\), taking values 1, 2, 3 be

``` r
(pCond.v.given.u0 <- c(.7, .2, .1))
```

    ## [1] 0.7 0.2 0.1

``` r
(pCond.v.given.u1 <- c(.1,.2,.7))
```

    ## [1] 0.1 0.2 0.7

Let random variable \((w|v,u)\) take values \(1,2,3\) with probabilities
\(p(w|v,u)\), given by the following:

``` r
p.given.u0.v1 <- c(.3, .3, .4)
p.given.u0.v2 <- c(.5, .3, .2)
p.given.u0.v3 <- c(.6, .2, .2)
p.given.u1.v1 <- c(.2, .3, .5)
p.given.u1.v2 <- c(.2, .2, .6)
p.given.u1.v3 <- c(.1, .7, .2)
```

Simulate joint sample \((w,v,u)\) of lenth \(n=500\). Use `set.seed(11)`
Start with simulation of u. For each simulated value `u` generate `v`
from the corresponding conditional distribution \(p(v|u)\). Finally, for
each pair \(v,u\) simulate `w` from \(p(w|v,u)\).

``` r
v_probs <- function(u) {
  switch(u + 1, 
         pCond.v.given.u0,
         pCond.v.given.u1)
}

w_probs <- function(u, v) {
  switch(u + 1,
         switch(v,
                p.given.u0.v1,
                p.given.u0.v2,
                p.given.u0.v3),
         switch(v,
                p.given.u1.v1,
                p.given.u1.v2,
                p.given.u1.v3)
         )
}

set.seed(11)
sim.u <- rbinom(500, 1, prob = .45)
sim.v <- map_dbl(map(sim.u, ~ v_probs(.x)), ~ sample(1:3, 1, replace = TRUE, prob = .x))
sim.w <- map_dbl(map2(sim.u, sim.v, ~ w_probs(.x, .y)), ~ sample(1:3, 1, replace = TRUE, prob = .x))

head(cbind(sim.w, sim.v, sim.u), 25)
```

    ##       sim.w sim.v sim.u
    ##  [1,]     1     1     0
    ##  [2,]     1     1     0
    ##  [3,]     3     1     0
    ##  [4,]     1     1     0
    ##  [5,]     2     2     0
    ##  [6,]     1     2     1
    ##  [7,]     2     3     0
    ##  [8,]     1     1     0
    ##  [9,]     2     3     1
    ## [10,]     1     1     0
    ## [11,]     1     1     0
    ## [12,]     3     3     0
    ## [13,]     3     2     1
    ## [14,]     1     2     1
    ## [15,]     2     3     1
    ## [16,]     1     1     1
    ## [17,]     2     1     0
    ## [18,]     1     1     0
    ## [19,]     3     1     0
    ## [20,]     2     1     0
    ## [21,]     1     1     0
    ## [22,]     3     3     1
    ## [23,]     1     1     0
    ## [24,]     2     1     0
    ## [25,]     2     1     0

# 3 Conditional Expected Values

Calculate unconditional expected value \(E[v]\). Random variable `v` can
take values \(1,2,3\) with corresponding probabilities.

``` r
vMarginal
```

    ##   v1   v2   v3 
    ## 0.60 0.35 0.05

Then the unconditional mean value is

``` r
c(1, 2, 3) %*% vMarginal
```

    ##      [,1]
    ## [1,] 1.45

Calculate conditional expected value \(E[v|u]\).

First, find conditional mean values \(E[v|u=u0]=E[v|u0]\) and
\(E[v|u=u1]=E[v|u1]\).

The random variable \((v|u0)\) takes values \(1,2,3\) with corresponding
probabilities \(p(v=1|u0)\), \(p(v=2|u0)\), \(p(v=3|u0)\), given by the
vector.

``` r
(cond.v.given.u0 <- apply(data3way.array.p[, , "u0"], 2, sum) / uMarginal["u0"])
```

    ##         v1         v2         v3 
    ## 0.61764706 0.32352941 0.05882353

Taking conditional expected value with respect to this distribution
\(E[v|u0]\) is:

``` r
(exp.v.given.u0 <- c(1, 2, 3) %*% cond.v.given.u0)
```

    ##          [,1]
    ## [1,] 1.441176

The random variable \((v|u1)\) takes the same values \(1,2,3\), but with
different probabilities \(p(v|u1)\):

``` r
cond.v.given.u1
```

    ##         v1         v2         v3 
    ## 0.59090909 0.36363636 0.04545455

Thus conditional expected value \(E[v|u1]\) is

``` r
(exp.v.given.u1 <- c(1, 2, 3) %*% cond.v.given.u1)
```

    ##          [,1]
    ## [1,] 1.454545

Note that the conditional expected value \(E[v|u]\) takes two different
values: \(E[v|u0]\) and \(E[v|u1]\) with probabilities \(p(v|u0)\) and
\(p(v|u1)\), correspondingly.

This means that \(E[v|u]\) is a random variable and its (unconditional)
expected value is \[E[v|u]=E[v|u0]p(u=u0)+E[v|u1]p(u=u1)\]

Calculate this unconditional expected value:

``` r
(uncond.exp.v.given.u <-
   exp.v.given.u0 * uMarginal["u0"] + exp.v.given.u1 * uMarginal["u1"])
```

    ##      [,1]
    ## [1,] 1.45
