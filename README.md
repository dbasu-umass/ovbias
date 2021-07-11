
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ovbias

<!-- badges: start -->

<!-- badges: end -->

The goal of the `ovbias` package is to present some functions to compute
quantiles of the empirical distribution of the bias-adjusted treatment
effect in a linear econometric model with omitted variable bias. To
analyze such models, one can consider three regression models: (a) a
short regression model where the outcome variable is regressed on the
treatment variable, with or without additional controls; (b) an
intermediate regression model where additional control variables are
added to the short regression; (c) a hypothetical long regression model,
where and index of the omitted variable(s) is added to the intermediate
regressions. The selection on unobservables is captured by the parameter
`delta`; and the R-squared in the hypothetical long regression is
captured by the parameter `Rmax`. In addition, we need to consider an
auxiliary regression where the treatment variable is regressed on all
observed (and included) control variables.

This package provides two functions to compute quantiles of the
empirical distribution of the bias-adjusted treatment effect, `bset1()`
and `bset2()`.

To use the `bset1()` function, the user will need to run the short,
intermediate and auxiliary regressions, collect relevant parameters and
supply them to `bset`()\`. The relevant parameters are:

  - beta0: treatment effect in short regression
  - R0: R-squared in short regression
  - betatilde: treatment effect in intermediate regression
  - Rtilde: R-squared in intermediate regression
  - taux: variance of residual in auxiliary regression
  - sigmay: standard deviation of the outcome variable
  - sigmax: standard deviation of the treatment variable.

To use the `bset2()` function, the user will need to run supply a data
frame (with all relevant variables), name of the outcome variable, name
of the treatment variable, and formulas for the short, intermediate and
auxiliary regressions. The function will run the short, intermediate and
auxiliary regressions, and collect the relevant parameters internally.

In both functions, the relevant parameters are used to construct a cubic
equation in the omitted variable bias. The cubic equation is solved on a
NxN grid over a bounded box of the `(delta,Rmax)` plane, and the roots
of the cubic are used to construct the bias-adjusted treatment effect.
The bias-adjusted treatment effect is, by definition, the treatment
effect estimated from the intermediate regression, `betatilde`, and the
real root of the cubic equation.

The bounded box of the `(delta,Rmax)` plane is defined by the Cartesian
products of two intervals: `[deltalow, deltahigh]`, and `[Rtilde,
Rmax]`. We consider 4 different regimes: (1) low `delta`, low `Rmax`,
where `deltalow < delta < 1` and `Rtilde < Rmax < 1.30 Rtilde`; (2) low
`delta`, high `Rmax`, where `deltalow < delta < 1` and `Rtilde < Rmax
< 2.47 Rtilde`; (3) high `delta`, low `Rmax`, where `1 < delta <
deltahigh` and `Rtilde < Rmax < 1.30 Rtilde`; (4) high `delta`, high
`Rmax`, where `1 < delta < deltahigh` and `Rtilde < Rmax < 2.47 Rtilde`.

For each regime, the bounded box is divided into two parts: (a) URR
region, over which the cubic equation has a unique real root, and (b)
NURR region, over which the cubic equation has three real roots. For the
NURR region, the root whose distribution is closest to the distribution
of the root over the URR region is chosen. If the URR region is empty,
unique real roots of the cubic cannot be computed. Hence, the
bias-adjusted treatment effect cannot be computed.

The output of the `bset1()` function is a list of two elements. The
first element is a data frame of the parameters about the three
regressions supplied by the user; the second element of the list is a
matrix containing the empirical distribution of the bias-adjusted
treatment effect for the four regimes.

The output of the `bset2()` function is a list of two elements. The
first element is a data frame of the parameters extracted from the three
regressions; the second element of the list is a matrix containing the
empirical distribution of the bias-adjusted treatment effect for the
four regimes.

For more details see Basu, D. (2021). “Bounding Sets for Treatment
Effects with Proportional Selection”. Economics Department Working Paper
Series. 307. University of Massachusetts Amhers. URL:
<https://scholarworks.umass.edu/econ_workingpaper/307>. For latest
version of the paper use URL:
<https://drive.google.com/file/d/1v-IGh9_pKqAu7ALucNGTWJXBx0A0STxJ/view?usp=sharing>

## Installation

You can install ovbias like so:

``` r
# install.packages("devtools")
devtools::install_github("dbasu-umass/ovbias")
```

## Example

### Creating a Simulated Data Set

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

### Using the `bset2()` function

We would like to use the `bset2()` function to find bounding sets of the
bias-adjusted treatment effect. We will use `deltalow=0.01` and
`deltahigh=5.0` in our call of the `bset2()` function. We have to supply
the data frame, the names of the outcome and treatment variables, and
formulas for the short, intermediate and auxiliary regressions.

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
`$bsetpar`, and recall:

  - beta0: treatment effect in short regression
  - R0: R-squared in short regression
  - betatilde: treatment effect in intermediate regression
  - Rtilde: R-squared in intermediate regression
  - taux: variance of residual in auxiliary regression
  - sigmay: standard deviation of the outcome variable
  - sigmax: standard deviation of the treatment variable.

To see the quantiles of the empirical distribution of the bias-adjusted
treatment effect look at the second element of the list, `$bsetqnt`. The
rownames of this matrix is informative about the region used for
constructing the bias-adjusted treatment effect. For example, the LL
(URR) row refers to a low `delta`, low `Rmax` regime, where the roots
have been computed over the URR area. Aa another example, note that the
LH (NURR) row refers to a low `delta`, high `Rmax` regime, where the
roots have been computed over the NURR area. Recall that for the NURR
area, there are three real roots and we choose the root whose empirical
distribution is closest to the empirical distribution of the root
computed over the URR area (using the same regime).

### Using the `bset1()` function

We would like to use the `bset1()` function to compute quantiles of the
bias-adjusted treatment effect. We will use `deltalow=0.01` and
`deltahigh=5.0` in our call of the `bset1()` function - to compare our
results with the use of the `bset2()` function above. We have to run the
short, intermediate and auxiliary regressions, collect the relevant
paramters and supply it to `bset1()`.

Let us run the short regression and collect the relevant parameters.

``` r
sreg <- stats::lm(y~x1,data=d1)
beta0 <- sreg$coefficients["x1"]
R0 <- summary(sreg)$r.squared
```

Let us run the intermediate regression and collect the relevant
parameters.

``` r
ireg <- stats::lm(y~x1+x2+x3,data=d1)
betatilde <- ireg$coefficients["x1"]
Rtilde <- summary(ireg)$r.squared
```

Let us run the auxiliary regression and store the variance of the
residuals.

``` r
auxreg <- stats::lm(x1~x2+x3,data=d1)
taux <- var(auxreg$residuals, na.rm = TRUE)
```

Let us store the standard deviation of the outcome variable.

``` r
sigmay <- sd(d1$y, na.rm = TRUE)
```

Let us store the standard deviation of the treatment variable.

``` r
sigmax <- sd(d1$x1, na.rm=TRUE)
```

We have generated all the relevant parameters and can now use these to
call the `bset1()` function.

``` r
bset1(
 beta0 = beta0,
 R0 = R0,
 betatilde = betatilde,
 Rtilde = Rtilde,
 sigmay = sigmay,
 sigmax = sigmax,
 taux = taux,
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

Note that the results are identical to the ones we got above using the
`bset2()` function.
