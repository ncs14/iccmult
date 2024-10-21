
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iccmult

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/ncs14/iccmult/graph/badge.svg)](https://app.codecov.io/gh/ncs14/iccmult)
<!-- badges: end -->

The goal of iccmult is to estimate the intracluster correlation
coefficient (ICC) of clustered categorical response data. It provides
two estimation methods, a resampling based estimator and the method of
moments estimator. These are obtained by specifying a method in the
function `iccmulti::iccmult()`. This package also includes a function to
generate simulated clustered categorical response data:
`iccmulti::rccat()`.

## Installation

You can install iccmult from within R or RStudio with:

``` r
install.packages("iccmult")
```

Alternatively, install the package from [GitHub](https://github.com/)
with:

``` r
# install.packages("pak")
pak::pak("ncs14/iccmult")
```

## Example

This is a basic example which shows you how to generate clustered
categorical response data. The response probabilities must sum 1 and the
desired ICC must be a value between 0 and 1. The number of clusters is
set to 20 and each cluster is of size 25. The output of `rccat` is a two
column data frame of a cluster identifier and a categorical response
vector.

``` r
library(iccmult)
clustdat3 <- rccat(rho=0.25, prop=c(0.2,0.3,0.5), noc=20, csize=25)
```

The `iccmulti()` function is called as follows to estimate the ICC on
the resulting data frame. The function expects two variables: a cluster
identifier and the categorical response vector. The call below requests
both the resampling and the moments estimates.

``` r
iccclust <- iccmulti(cid, y, clustdat3, method=c("rm","mom"))
```

The result is a list of length two, each component holding the estimated
ICC, se(ICC), and confidence interval bounds from each estimation
method.

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
