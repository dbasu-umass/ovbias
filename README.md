
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
#>   beta0    R0 betatilde Rtilde sigmay sigmax  taux
#> 1 0.625 0.021     0.442  0.714  4.716  1.096 1.198
#> 
#> $bsetqnt
#>           Area   2.5%    5%      50%     95%     97.5%  
#> LL (URR)  "0.54" "0.417" "0.42"  "0.436" "0.441" "0.442"
#> LL (NURR) "0.45" "0.393" "0.397" "0.42"  "0.438" "0.439"
#> LH (URR)  "0.54" "0.409" "0.413" "0.434" "0.441" "0.442"
#> LH (NURR) "0.45" "0.377" "0.382" "0.413" "0.437" "0.439"
#> HL (URR)  "0"    "NA"    "NA"    "NA"    "NA"    "NA"   
#> HL (NURR) "1"    "NA"    "NA"    "NA"    "NA"    "NA"   
#> HH (URR)  "0"    "NA"    "NA"    "NA"    "NA"    "NA"   
#> HH (NURR) "1"    "NA"    "NA"    "NA"    "NA"    "NA"
```

To see the relevant parameters that emerged from the short, intermediate
and auxiliary regressions, look at the first element of the list,
`$bsetpar`; to see the quantiles of the empirical distribution of the
bias-adjusted treatment effect look at the second element of the list,
`$bsetqnt`.