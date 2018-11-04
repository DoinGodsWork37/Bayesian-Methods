Week 5 Assignment: Binomial Output with Two Groups
================
Your Name
October, 28 2018

**This assignment uses data from the course project** **This project is
individual**

# 1 Data

The course project is based on the following data: “Reuters/Ipsos Poll
data, {DATE RANGE}”

Literature for the project is uploaded on the course page:

  - [Gelman
    (2007)](http://ilykei.com/api/fileProxy/assignments%2FBayesian%20Methods%2FCourse%20Project%2032014%2FGelman.%20%202007.%20%20Struggles%20with%20Survey%20Weighting%20...%20SS_.pdf)

  - Ghitza, Gelman

  - Burkner, BMRS Talk
(2016)

<!-- end list -->

``` r
data <- read_csv("assignments-Bayesian Methods-Course Project 32014-MScA_32014_BayesianMethods_CourseProjectData.csv")

data %>% dim()
```

    ## [1] 23223     6

``` r
data %>% head()
```

    ## # A tibble: 6 x 6
    ##   sex    race    age   education       state     y
    ##   <chr>  <chr>   <chr> <chr>           <chr> <int>
    ## 1 1.Male 1.White 18-24 1.NoCollege     GA        0
    ## 2 1.Male 1.White 25-34 1.NoCollege     AZ        0
    ## 3 1.Male 1.White 25-34 2.SomeCollege   SD        0
    ## 4 1.Male 1.White 18-24 3.CollegeOrMore SC        0
    ## 5 1.Male 1.White 18-24 3.CollegeOrMore SC        0
    ## 6 1.Male 1.White 18-24 3.CollegeOrMore SC        0

``` r
data %>% distinct(sex, race)
```

    ## # A tibble: 8 x 2
    ##   sex      race      
    ##   <chr>    <chr>     
    ## 1 1.Male   1.White   
    ## 2 2.Female 1.White   
    ## 3 1.Male   2.Black   
    ## 4 2.Female 2.Black   
    ## 5 1.Male   3.Hispanic
    ## 6 2.Female 3.Hispanic
    ## 7 1.Male   4.Other   
    ## 8 2.Female 4.Other

``` r
(data_summary <- data %>% 
  group_by(sex) %>% 
  summarise(N = n(),
            vote = sum(y)))
```

    ## # A tibble: 2 x 3
    ##   sex          N  vote
    ##   <chr>    <int> <int>
    ## 1 1.Male    8796  2995
    ## 2 2.Female 14427  5519

``` r
fit <- stan(file = "model", data = datalist, iter = 2000, warmup = 500, chains = 2, seed = 42, refresh = 2000)
```

``` r
# Bias sample because divergences
# Proper way to fix is to change parameterization
print(fit)
```

``` r
# adapt adjust steps; algorithm can get stuck in crack; as a result, divergences occur; another way to fix this is to change priors
# Changing steps does not gurantee a better model
fit_cp85 <- stan(file = "model", control = list(adapt_delta=.95))
```
