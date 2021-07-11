
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ovbias

<!-- badges: start -->

<!-- badges: end -->

The goal of ovbias is to compute quantiles of the empirical distribution
of the bias-adjusted treatment effect in a linear econometric model with
omitted variable bias. To analyze such models, one can consider three
regression models: (a) a short regression model where the outcome
variable is regressed on the treatment variable, with or without
additional controls; (b) an intermediate regression model where
additional control variables are added to the short regression; (c) a
hypothetical long regression model, where and index of the omitted
variable(s) is added to the intermediate regressions. The selection on
unobservables is captured by the parameter `delta`; and the R-squared in
the hypothetical long regression is captured by the parameter `Rmax`.

The main function in this package runs three regressions: (a) the short
regression, (b) the intermediate regression, and (c) an auxiliary
regression where the treatment variable is regression on all included
control variables. Relevant parameters are collected from these three
regressions and used to construct a cubic equation in the omitted
variable bias. The cubic equation is solved on a NxN grid over a bounded
box of the (`delta`,`Rmax`) plane, and the roots of the cubic are used
to construct the bias-adjusted treatment effect.

The bounded box of the (`delta`,`Rmax`) plane is defined by the
Cartesian products of two intervals: \[`deltalow`, `deltahigh`\], and
\[`Rtilde`, `Rmax`\]. We posit two possibilities for `Rmax`: (a)
`1.3*Rtilde`, (b) `2.47*Rtilde`. Hence, we can define 4 different
regimes: (1) low `delta`, low `Rmax`, where `deltalow < delta < 1` and
`Rtilde < Rmax < 1.30 Rtilde`; (2) low `delta`, high `Rmax`, where
`deltalow < delta < 1` and `Rtilde < Rmax < 2.47 Rtilde`; (3) high
`delta`, low `Rmax`, where `1 < delta < deltahigh` and `Rtilde < Rmax
< 1.30 Rtilde`; (4) high `delta`, high `Rmax`, where `1 < delta <
deltahigh` and `Rtilde < Rmax < 2.47 Rtilde`.

The output of the function is a list of two elements. The first element
is a data frame of the parameters extracted from the three regressions;
the second element of the list is a matrix containing the empirical
distribution of the bias-adjusted treatment effect for the four regimes.

For more details see Basu, D. (2021). “Bounding Sets for Treatment
Effects with Proportional Selection”. Economics Department Working Paper
Series. 307. University of Massachusetts Amhers. URL:
<https://scholarworks.umass.edu/econ_workingpaper/307>

## Installation

You can install ovbias like so:

``` r
# install.packages("devtools")
devtools::install_github("dbasu-umass/ovbias")
```

## Example

Let us create a simulated data set.

``` r
library(ovbias)

# Parameters
a0 <- 0.5
a1 <- 0.25
a2 <- 0.75
a3 <- 1.0
a4 <- -0.75 
a5 <- -0.25

# Regressors
set.seed(999)
x1 <- rnorm(100, mean=0, sd=1)
x2 <- rnorm(100, mean=1, sd=2)
x3 <- rnorm(100, mean=2, sd=3)
x4 <- rnorm(100, mean=3, sd=4)
x5 <- rep(c(0,1),50)

# Outcome variable
y <- a0 + a1*x1 + a2*x2 + a3*x3 + a4*x4 + a5*x5

# Create data frame
d1 <- data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5)
```

The researcher is interested in estimating the effect of the treatment
variable, `x1`, on the outcome variable, `y`. Suppose the researcher is
unable to collect data on the following two variables `x4` and `x5`.
When the researcher estimates a linear model, the treatment effect is
estimated with bias (because of the two omitted variables).

We would like to use the `bset2()` function to find bounding sets of the
bias-adjusted treatment effect. We will use `deltalow=0.01` and
`deltahigh=5.0` in our call of the `bset2()` function.

``` r
# Generate results
bset2(
  data = d1,
  outcome = "y",
  treatment = "x1",
  shortreg = y ~ x1,
  intreg = y ~ x1 + x2 + x3,
  auxreg = x1 ~ x2 + x3, 
  deltalow = 0.01, 
  deltahigh=5.00
  )
#> $bsetpar
#>    beta0 R0 betatilde Rtilde sigmay sigmax  taux
#> 1 -0.055  0      0.84  0.559  4.325  0.943 0.802
#> 
#> $bsetqnt
#>            Area    2.5%      5%     50%    95%  97.5%
#> LL (URR)  0.740   0.842   0.844   0.883  0.988  1.005
#> LL (NURR) 0.250   0.859   0.865   0.969  1.079  1.092
#> LH (URR)  0.734   0.844   0.847   0.945  1.223  1.269
#> LH (NURR) 0.256   0.868   0.884   1.169  1.474  1.507
#> HL (URR)  0.080  -9.319  -9.300  -9.179 -9.116 -9.112
#> HL (NURR) 0.920 -12.544 -12.124  -9.762 -9.183 -9.155
#> HH (URR)  0.344  -9.700  -9.642  -9.265 -9.095 -9.083
#> HH (NURR) 0.656 -12.595 -12.274 -10.059 -9.234 -9.181
```

To see the relevant parameters that emerged from the short, intermediate
and auxiliary regressions, look at the first element of the list,
`$bsetpar`; to see the quantiles of the empirical distribution of the
bias-adjusted treatment effect look at the second element of the list,
`$bsetqnt`.
