
#' Bounding set for treatment effect
#'
#' This funcition runs short, intermediate and auxiliary regressions; collect parameters; and, compute quantiles of the empirical distribution of the treatment effect.
#' @param data dataframe containing all the relevant variables
#' @param outcome name of outcome variable in the dataframe, e.g. "y"
#' @param treatment name of treatment variable in the dataframe, e.g. "x1"
#' @param shortreg short regression as formula; e.g. y ~ x1
#' @param intreg intermediate regression as formula; e.g. y ~ x1 + x2
#' @param auxreg auxiliary regression as formula; e.g. x1 ~ x2
#' @param N used to contruct the grid on which the cubic equation is solved; the grid is N*N; default is N=100
#' @param Rlow Rmax in low regime is Rlow*Rtilde; default is Rlow=1.30
#' @param Rhigh Rmax in low regime is Rlow*Rtilde; default is Rlow=2.47
#' @param deltalow lower bound of delta; this real number must be strictly less than 1
#' @param deltahigh lower bound of delta; this real number must be strictly greater than 1
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @return bset2 returns a list with the following two elements:
#'
#' \item{bsetpar}{a dataframe containing the parameters collected from the three regressions}
#' \item{bsetqnt}{a matrix containing the quantiles of the empirical distribution of the bias-adjusted treatment effect}
#'
#' @references Basu, D. (2021). "Bounding Sets for Treatment Effects with Proportional Selection". Economics Department Working Paper Series. 307. University of Massachusetts Amhers. URL: https://scholarworks.umass.edu/econ_workingpaper/307
#'
#' @export
#'
#'
#' @examples
#' # Generate simulated data
#' # Parameters
#' a0 <- 0.5;a1 <- 0.25;a2 <- 0.75;a3 <- 1.0
#' a4 <- -0.75; a5 <- -0.25
#' # Regressors
#' set.seed(999)
#' x1 <- rnorm(100, mean=0, sd=1);x2 <- rnorm(100, mean=1, sd=2)
#' x3 <- rnorm(100, mean=2, sd=3);x4 <- rnorm(100, mean=3, sd=4)
#' x5 <- rep(c(0,1),50)
#' # Dependent variable
#' y <- a0 + a1*x1 + a2*x2 + a3*x3 + a4*x4 + a5*x5
#' # Create data frame
#' d1 <- data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5)
#' # Generate bounding set for treatment effect with x4 and x5 omitted
#' bset2(data = d1,
#' outcome = "y",
#' treatment = "x1",
#' shortreg = y ~ x1,
#' intreg = y ~ x1 + x2 + x3,
#' auxreg = x1 ~ x2 + x3,
#' deltalow = 0.01,
#' deltahigh=5.00)

bset2 <- function(
  data,
  outcome,
  treatment,
  shortreg,
  intreg,
  auxreg,
  N=100,
  Rlow=1.30,
  Rhigh=2.47,
  deltalow,
  deltahigh
){

  # Is the dependent variable in the short regression the
  # specified outcome variable?
  if (sum(all.vars(shortreg[[2]]) != paste(outcome))==1) {
    stop("Note: Dependent variable in short regression incorrect\n")
  }

  # Is the dependent variable in the intermediate regression the
  # specified outcome variable?
  if (sum(all.vars(intreg[[2]]) != paste(outcome))==1) {
    stop("Note: Dependent variable in short regression incorrect\n")
  }

  # Is the dependent variable in the auxiliary regression the
  # specified treatment variable?
  if (sum(all.vars(auxreg[[2]]) != paste(treatment))==1) {
    stop("Note: Dependent variable in short regression incorrect\n")
  }

  # Make sure xlow <1 and xhigh > 1
  if (deltalow>=1) {
    stop("Note: deltalow must be less than 1\n")
  }
  if (deltahigh<=1) {
    stop("Note: deltahigh must be greater than 1\n")
  }

  # -------------------------------------------------------------------- #
  # ---------------------- Create Funcitons ---------------------------- #

  # ---- Function 1
  # This funciton identifies the region on the plane
  # that gives a unique real solution of the cubic equation
  # corresponding to a given "delta" and "Rmax"
  # by evaluating the discriminant

  mydisc <- function(mydelta,Rmax){
    # Coefficients
    a <- (mydelta -1)*(taux)*((sigmax^2) - taux)
    b <- taux*(beta0 - betatilde)*(sigmax^2)*(mydelta-2)
    c <- (mydelta)*(Rmax - Rtilde)*(sigmay^2)*(sigmax^2 - taux) -
      (Rtilde - R0)*(sigmay^2)*taux - ((beta0- betatilde)^2)*taux*(sigmax^2)
    d <-  (mydelta)*(Rmax - Rtilde)*(sigmay^2)*(beta0- betatilde)*(sigmax^2)
    # Define P and Q
    P <- (3*a*c - (b^2))/(3*(a^2))
    Q <- (27*(a^2)*d + 2*(b^3) - 9*a*b*c)/(27*(a^3))
    # Evaluate discriminant
    D <- 4*(P^3)+27*(Q^2)
    # Is D>0?
    dpos <- ifelse(D>0,1,0)
    return(dpos)
  }



  # ---- Function 2
  # This funciton computes the solutions of a cubic
  # equation given its real coefficients: a, b, c, d
  # Notation for cubic: a*x^3 + b*x^2 + c*x + d = 0
  # The coefficients are functions of mydelta and Rmax

  mycubic <- function(mydelta,Rmax){
    # Coefficients
    a <- (mydelta -1)*(taux)*((sigmax^2) - taux)
    b <- taux*(beta0 - betatilde)*(sigmax^2)*(mydelta-2)
    c <- (mydelta)*(Rmax - Rtilde)*(sigmay^2)*(sigmax^2 - taux) -
      (Rtilde - R0)*(sigmay^2)*taux - ((beta0- betatilde)^2)*taux*(sigmax^2)
    d <-  (mydelta)*(Rmax - Rtilde)*(sigmay^2)*(beta0- betatilde)*(sigmax^2)
    # Solve cubic
    x0 <- RConics::cubic(c(a,b,c,d))
    return(x0)
  }




  # ---- Function 3
  # This function computes the empirical distribution of
  # betastar for a given betatilde, Rtilde and a given value of
  # xlow and xhigh by choosing mydelta to lie in the interval [xlow,xhigh].
  # The calculation is done over the URR region

  qnturr <- function(xlow, xhigh, betatilde,
                     Rtilde, RUB, N=100) {

    z <- NULL

    # Delta
    x1 <- seq(xlow,xhigh,length.out = N)
    # Rmax
    y1 <- seq(Rtilde+0.01,RUB,length.out = N)
    # All combinations of delta, Rmax
    mtrx2d <-  expand.grid(x1,y1)
    # Create data frame for plot
    consol_temp <- data.frame(x=mtrx2d[,1],y=mtrx2d[,2])
    consol <- consol_temp %>%
      dplyr::mutate(
        # Evalutate function at grid points
        z = mydisc(.data$x,.data$y)
        ) %>%
      as.data.frame()


    #---- Empirical distribution of bstar
    # Choose observations with z==1
    # z==1 are points where discriminant>0 (URR)
    data1 <- consol %>%
      dplyr::filter(z==1) %>%
      as.data.frame()
    # How many observations?
    N <- nrow(data1)

    # If N>0 computations possible
    if(N>0){
      # Create placeholder for bias
      bias <- rep(0,N)
      # Compute root of cubic (note:root=bias)
      for (i in 1:N) {
        bias[i] <- Re(mycubic(data1$x[i],data1$y[i])[1])
      }
      # Create data frame
      data2 <- as.data.frame(cbind(bias))
      # Create data frame for betastar
      data2 <- as.data.frame(betatilde-cbind(data2$bias))
      colnames(data2) <- "bstar"

      # Distribution of bstar
      rbdset <- c(
        round(
          stats::quantile(data2$bstar,
                          probs = c(0.025,0.05,0.5,0.95,0.975)),
          digits = 3)
      )
      # Return fraction of area, and vector of quantiles
      return(c(round(N/nrow(mtrx2d),digits = 3),
               rbdset))
    }else{
      # If N=0, bias cannot be computed
      return(c(0,"NA","NA","NA","NA","NA"))
    }
  }



  # ---- Function 4
  # This function computes the empirical distribution of
  # betastar for a given betatilde, Rtilde and a given value of
  # xlow and xhigh by choosing delta = [xlow,xhigh]
  # The calculation is done over the NURR region

  qntnurr <- function(xlow, xhigh, betatilde,
                      Rtilde, RUB, N=100) {
    z <- NULL
    # Delta
    x1 <- seq(xlow,xhigh,length.out = N)
    # Rmax
    y1 <- seq(Rtilde+0.01,RUB,length.out = N)
    # All combinations of delta, Rmax
    mtrx2d <-  expand.grid(x1,y1)
    # Create data frame for plot
    consol_temp <- data.frame(x=mtrx2d[,1],y=mtrx2d[,2])
    consol <- consol_temp %>%
      dplyr::mutate(
        # Evalutate function at grid points
        z = mydisc(.data$x,.data$y)
      ) %>%
      as.data.frame()


    #---- Empirical distribution of betastar

    # -- URR: Choose observations with z==1
    # Choose observations with z==1
    # z==1 are points where discriminant>0 (URR)
    data1 <- consol %>%
      dplyr::filter(z==1) %>%
      as.data.frame()
    # How many observations?
    N <- nrow(data1)

    # If N>0 computations possible
    if(N>0){
      # Create placeholder for bias
      biasurr <- rep(0,N)
      # Compute root of cubic (note:root=bias)
      for (i in 1:N) {
        biasurr[i] <- Re(mycubic(data1$x[i],data1$y[i])[1])
      }
    }else{
      # If N=0, bias cannot be computed
      biasurr <- rep(NA,N)
    }
    # If all elements of "biasurr" is NA, no computaton posssible
    if(sum(is.na(biasurr))==N){
      return(c(1,"NA","NA","NA","NA","NA"))
    }else{
      # -- NURR: Choose observations with z==0
      data1 <- consol %>%
        dplyr::filter(z==0) %>%
        as.data.frame()
      # How many observations?
      N <- nrow(data1)

      # If N>0 computations possible
      if(N>0){
        # Create placeholder for three roots
        biasnurr <- matrix(data=NA,nrow = N,ncol = 3)
        # Compute 3 real roots of cubic and store in matrix
        for (i in 1:N) {
          biasnurr[i,] <- mycubic(data1$x[i],data1$y[i])
        }
        # Create placeholder for vector of bias
        bias <- rep(0,N)
        # Populate 'bias' with the root closest to root from URR region
        c1 <- abs(stats::median(biasurr)-stats::median(biasnurr[,1]))
        c2 <- abs(stats::median(biasurr)-stats::median(biasnurr[,2]))
        c3 <- abs(stats::median(biasurr)-stats::median(biasnurr[,3]))
        if(c1==min(c1,c2,c3)){
          bias <- biasnurr[,1]
        }else{
          if(c2==min(c1,c2,c3)){
            bias <- biasnurr[,2]
          }else{
            bias <- biasnurr[,3]
          }
        }

        # Create data frame
        data2 <- as.data.frame(cbind(bias))
        # Create data frame for betastar
        data2 <- as.data.frame(betatilde-cbind(data2$bias))
        colnames(data2) <- "bstar"

        # Distribution of bstar
        rbdset <- c(
          round(
            stats::quantile(data2$bstar,
                            probs = c(0.025,0.05,0.5,0.95,0.975)),
            digits = 3)
        )
        # Return fraction of area and vector of quantiles
        return(c(round(N/nrow(mtrx2d),digits = 3),
                 rbdset))
      }else{
        # If N=0, bias cannot be computed
        return(c(0,"NA","NA","NA","NA","NA"))
      }

    }

  }

  # ------------------------------------------------------------------ #
  # ------- Run the three regressions and collect parameters --------- #
  # Outcome variable
  myd1 <- data[,paste(outcome)]
  # Treatment variable
  myd2 <- data[,paste(treatment)]

  # Short regression
  reg1 <- stats::lm(shortreg,data = data)
  # Intermediate regression
  reg2 <- stats::lm(intreg,data = data)
  # Auxiliary regression
  reg3 <- stats::lm(auxreg,data = data)

  #--- Parameters
  # From Short regression
  beta0 <- stats::coefficients(reg1)[paste(treatment)]
  R0 <- summary(reg1)$r.squared
  # From Intermediate regression
  betatilde <- stats::coefficients(reg2)[paste(treatment)]
  Rtilde <- summary(reg2)$r.squared
  # From Auxiliary regression
  taux <- stats::var(reg3$residuals, na.rm = TRUE)
  # Outcome variable
  sigmay <- stats::sd(myd1, na.rm = TRUE)
  # Treatment variable
  sigmax <- stats::sd(myd2, na.rm = TRUE)

  # ---- Collect parameters as a data frame
  myparm <- data.frame(
    beta0=round(beta0,digits = 3),
    R0=round(R0,digits = 3),
    betatilde=round(betatilde,digits = 3),
    Rtilde=round(Rtilde,digits = 3),
    sigmay=round(sigmay,digits = 3),
    sigmax=round(sigmax,digits = 3),
    taux=round(taux,digits = 3)
  )

  # ---------------------------------------------------------------- #
  # ------- Quantiles of empirical distribution of betastar -------- #

  # --- Over the URR Area
  # Regime: Low delta, Low Rmax
  qu11_1 <- qnturr(
    xlow = deltalow,
    xhigh = 1,
    RUB = min(Rlow*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: Low delta, High Rmax
  qu12_1 <- qnturr(
    xlow = deltalow,
    xhigh = 1,
    RUB = min(Rhigh*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: High delta, Low Rmax
  qu21_1 <- qnturr(
    xlow = 1.01,
    xhigh = deltahigh,
    RUB = min(Rlow*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: High delta, High Rmax
  qu22_1 <- qnturr(
    xlow = 1.01,
    xhigh = deltahigh,
    RUB = min(Rhigh*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )



  # --- Over the NURR area
  # Regime: Low delta, Low Rmax
  qnu11_1 <- qntnurr(
    xlow = deltalow,
    xhigh = 1,
    RUB = min(Rlow*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: Low delta, High Rmax
  qnu12_1 <- qntnurr(
    xlow = deltalow,
    xhigh = 1,
    RUB = min(Rhigh*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: High delta, Low Rmax
  qnu21_1 <- qntnurr(
    xlow = 1.01,
    xhigh = deltahigh,
    RUB = min(Rlow*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Regime: High delta, High Rmax
  qnu22_1 <- qntnurr(
    xlow = 1.01,
    xhigh = deltahigh,
    RUB = min(Rhigh*myparm[,"Rtilde"],1),
    betatilde = myparm[,"betatilde"],
    Rtilde = myparm[,"Rtilde"]
  )

  # Quantiles of betastar as a matrix
  bset1 <- cbind(
    matrix(c(qu11_1, qnu11_1,
             qu12_1, qnu12_1,
             qu21_1, qnu21_1,
             qu22_1, qnu22_1),
           nrow = 8, ncol = 6, byrow = TRUE)
  )

  colnames(bset1) <- c("Area","2.5%","5%","50%","95%","97.5%")

  rownames(bset1) <- c("LL (URR)","LL (NURR)",
                       "LH (URR)","LH (NURR)",
                       "HL (URR)","HL (NURR)",
                       "HH (URR)","HH (NURR)")

  # ----------------------------------------------------- #
  # ------------------- Return results ------------------ #
  # remove rownmaes from dataframe
  rownames(myparm) <- NULL

  return(
    list(
      # bsetpar: data frame of parameters from
      # the three regressions
      bsetpar = myparm,
      # bsetqnt: matrix of quantiles of betastar
      bsetqnt = bset1
    )
  )

}


