## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message = FALSE---------------------------------------------------
library(RTMBdist)

## ----data---------------------------------------------------------------------
data(ChickWeight)

## ----parameters---------------------------------------------------------------
parameters <- list(
  mua=0,          # Mean slope
  log_sda=1,      # log-Std of slopes
  mub=0,          # Mean intercept
  log_sdb=1,      # log-Std of intercepts
  log_sigma=0,    # log-Scale of BCCG distribution
  nu = 0.1,       # Skewness of BCCG distribution
  a=rep(0, 50),   # Random slope by chick
  b=rep(5, 50)    # Random intercept by chick
)

## ----likelihood---------------------------------------------------------------
nll_chick <- function(parms) {
  getAll(ChickWeight, parms, warn=FALSE)
  # Optional (enables extra RTMB features)
  weight <- OBS(weight)
  # Initialise joint negative log likelihood
  nll <- 0
  # Random slopes
  sda <- exp(log_sda); ADREPORT(sda)
  nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
  # Random intercepts
  sdb <- exp(log_sdb); ADREPORT(sdb)
  nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
  # Data
  predWeight <- exp(a[Chick] * Time + b[Chick])
  sigma <- exp(log_sigma); ADREPORT(sigma)
  nll <- nll - sum(dbccg(weight, mu=predWeight, sigma=sigma, nu=nu, log=TRUE))
  # Get predicted weight uncertainties
  ADREPORT(predWeight)
  # Return
  nll
}

## ----model fitting------------------------------------------------------------
obj_chick <- MakeADFun(nll_chick, parameters, random=c("a", "b"), silent = TRUE)
opt_chick <- nlminb(obj_chick$par, obj_chick$fn, obj_chick$gr)

## ----laplace_check------------------------------------------------------------
set.seed(1)
checkConsistency(obj_chick)

## ----residuals----------------------------------------------------------------
osa_chick <- oneStepPredict(obj_chick, discrete=FALSE, trace=FALSE)
qqnorm(osa_chick$res); abline(0,1)

## ----cleanup1, include=FALSE--------------------------------------------------
TMB::FreeADFun(obj_chick)

## -----------------------------------------------------------------------------
data(InsectSprays)

## ----pardat-------------------------------------------------------------------
# Creating the model matrix
X <- model.matrix(~ spray - 1, data = InsectSprays)

par <- list(
  beta0 = log(mean(InsectSprays$count)),
  beta = rep(0, length(levels(InsectSprays$spray))),
  log_phi = log(1),
  log_sigma = log(1)
)
dat <- list(
  count = InsectSprays$count,
  spray = InsectSprays$spray,
  X = X
)

## ----nll_insect---------------------------------------------------------------
nll_insect <- function(par) {
  getAll(par, dat, warn=FALSE)
  count <- OBS(count)
  # Random effect likelihood
  sigma <- exp(log_sigma); ADREPORT(sigma)
  nll <- - sum(dnorm(beta, 0, sigma, log = TRUE))
  # Data likelihood
  lambda <- exp(beta0 + as.numeric(X %*% beta)); ADREPORT(lambda)
  phi <- exp(log_phi); ADREPORT(phi)
  nll <- nll -sum(dgenpois(count, lambda, phi, log = TRUE))
  nll
}

## ----model fit3---------------------------------------------------------------
obj_insect <- MakeADFun(nll_insect, par, random = "beta", silent = TRUE)
opt_insect <- nlminb(obj_insect$par, obj_insect$fn, obj_insect$gr)

# Checking if the Laplace approximation is adequate
checkConsistency(obj_insect)
# Check okay

# Calculating quantile residuals
osa_insect <- oneStepPredict(obj_insect, method = "oneStepGeneric", 
                             discrete=TRUE, trace=FALSE)
qqnorm(osa_insect$res); abline(0,1)

## ----cleanup2, include=FALSE--------------------------------------------------
TMB::FreeADFun(obj_insect)

## ----packages, message = FALSE------------------------------------------------
library(gamlss.data)   # Data
library(LaMa)          # Creating model matrices
library(Matrix)        # Sparse matrices

## ----data_boys----------------------------------------------------------------
data(dbbmi)
# Subset (just for speed here)
set.seed(1)
ind <- sample(1:nrow(dbbmi), 2000)
dbbmi <- dbbmi[ind, ]

## ----creating model matrices--------------------------------------------------
k <- 10 # Basis dimension
modmat <- make_matrices(~ s(age, bs="cs"), data = dbbmi)
X <- modmat$Z                              # Design matrix
S <- Matrix(modmat$S[[1]], sparse = TRUE)  # Sparse penalty matrix

## ----likelihood2--------------------------------------------------------------
nll_dbbmi <- function(par) {
  getAll(par, dat, warn=FALSE)
  bmi <- OBS(bmi)
  # Calculating response parameters
  mu <- exp(X %*% c(beta0_mu, beta_age_mu)); ADREPORT(mu) # Location
  sigma <- exp(X %*% c(beta0_sigma, beta_age_sigma)); ADREPORT(sigma) # Scale
  nu <- X %*% c(beta0_nu, beta_age_nu); ADREPORT(nu) # Skewness
  tau <- exp(log_tau); ADREPORT(tau) # Kurtosis
  # Data likelihood: Box-Cox power exponential distribution
  nll <- - sum(dbcpe(bmi, mu, sigma, nu, tau, log=TRUE))
  # Penalised splines as random effects: log likelihood / penalty
  lambda <- exp(log_lambda); REPORT(lambda)
  nll <- nll - dgmrf(beta_age_mu, 0, lambda[1] * S, log=TRUE)
  nll <- nll - dgmrf(beta_age_sigma, 0, lambda[2] * S, log=TRUE)
  nll <- nll - dgmrf(beta_age_nu, 0, lambda[3] * S, log=TRUE)
  nll
}

## ----parameters and data------------------------------------------------------
par <- list(
  beta0_mu = log(18), beta0_sigma = log(0.15),
  beta0_nu = -1, beta_age_mu = rep(0, k-1),
  beta_age_sigma = rep(0, k-1), beta_age_nu = rep(0, k-1),
  log_tau = log(2),
  log_lambda = log(rep(1e4, 3))
)
dat <- list(
  bmi = dbbmi$bmi,
  age = dbbmi$age,
  X = X,
  S = S
)

## ----REML fit-----------------------------------------------------------------
# Restricted maximum likelihood (REML) - also integrating out fixed effects
random <- names(par)[names(par) != "log_lambda"]
obj_dbbmi <- MakeADFun(nll_dbbmi, par, random = random, silent = TRUE)
opt_dbbmi <- nlminb(obj_dbbmi$par, obj_dbbmi$fn, obj_dbbmi$gr)

## ----results------------------------------------------------------------------
sdr <- sdreport(obj_dbbmi, ignore.parm.uncertainty = TRUE)
par <- as.list(sdr, "Est", report = TRUE)
par_sd <- as.list(sdr, "Std", report = TRUE)

## ----effects_plot-------------------------------------------------------------
age <- dbbmi$age
ord <- order(age)

# Plotting estimated effects
oldpar <- par(mfrow = c(1,3))
plot(age[ord], par$mu[ord], type = "l", lwd = 2, bty = "n", xlab = "Age", ylab = "Mu")
polygon(c(age[ord], rev(age[ord])),
        c(par$mu[ord] + 2*par_sd$mu[ord], rev(par$mu[ord] - 2*par_sd$mu[ord])),
        col = "#00000020", border = "NA")
plot(age[ord], par$sigma[ord], type = "l", lwd = 2, bty = "n", xlab = "Age", ylab = "Sigma")
polygon(c(age[ord], rev(age[ord])),
        c(par$sigma[ord] + 2*par_sd$sigma[ord], rev(par$sigma[ord] - 2*par_sd$sigma[ord])),
        col = "#00000020", border = "NA")
plot(age[ord], par$nu[ord], type = "l", lwd = 2, bty = "n", xlab = "Age", ylab = "Nu")
polygon(c(age[ord], rev(age[ord])),
        c(par$nu[ord] + 2*par_sd$nu[ord], rev(par$nu[ord] - 2*par_sd$nu[ord])),
        col = "#00000020", border = "NA")
par(oldpar)

## ----cond_dist----------------------------------------------------------------
# Plotting conditional distribution
plot(dbbmi$age, dbbmi$bmi, pch = 16, col = "#00000020",
     xlab = "Age", ylab = "BMI", bty = "n")
lines(age[ord], par$mu[ord], lwd = 3, col = "deepskyblue")

# Compute quantiles (point estimates)
par <- lapply(par, as.numeric)
ps <- seq(0, 1, length = 8)
ps[1] <- 0.005 # avoid 0 and 1
ps[length(ps)] <- 0.995 # avoid 0 and 1
for(p in ps) {
  q <- qbcpe(p, par$mu, par$sigma, par$nu, par$tau) # quantiles
  lines(age[ord], q[ord], col = "deepskyblue")
}
legend("topleft", lwd = c(3, 1), col = "deepskyblue", legend = c("Mean", "Quantiles"), bty = "n")

## ----cleanup3, include=FALSE--------------------------------------------------
TMB::FreeADFun(obj_dbbmi)

## ----data4--------------------------------------------------------------------
library(gamlss.data)
head(aep)

## ----binom4-------------------------------------------------------------------
# Defininig the model matrix for the model reported in Gange et al. (1996)
X <- model.matrix(~ age + ward + loglos * year, data = aep)

# (zero-inflated) binomial likelihood
nll_aep <- function(par) {
  getAll(par, dat)
  y <- OBS(y); size <- OBS(size)
  prob <- plogis(X %*% beta); ADREPORT(prob) # linear predictor and link
  zeroprob <- plogis(logit_zeroprob); ADREPORT(zeroprob)
  - sum(dzibinom(y, size, prob, zeroprob, log = TRUE))
}

# Initial parameters
beta_init <- c(-1, rep(0, ncol(X)-1))
names(beta_init) <- colnames(X)
par <- list(beta = beta_init)

dat <- list(
  y = aep$y[,1],
  size = aep$los,
  X = X
)

# Fitting the binomial model (zeroprob fixed at 0)
map <- list(logit_zeroprob = factor(NA)) # fixing at initial value
par$logit_zeroprob <- qlogis(0) # set to zero
obj_aep1 <- MakeADFun(nll_aep, par, silent = TRUE, map = map)
opt_aep1 <- nlminb(obj_aep1$par, obj_aep1$fn, obj_aep1$gr)

# Fitting the zero-inflated binomial model, no parameter restrictions
par$logit_zeroprob <- qlogis(1e-2) # more sensible initial value
obj_aep2 <- MakeADFun(nll_aep, par, silent = TRUE)
opt_aep2 <- nlminb(obj_aep2$par, obj_aep2$fn, obj_aep2$gr)

# Reporting
sdr_aep1 <- sdreport(obj_aep1)
sdr_aep2 <- sdreport(obj_aep2)

beta1 <- as.list(sdr_aep1, "Est")$beta
beta2 <- as.list(sdr_aep2, "Est")$beta
(zeroprob2 <- as.list(sdr_aep2, "Est", report = TRUE)$zeroprob)

round(rbind(beta1, beta2), 3)

## ----betabinom----------------------------------------------------------------
# Beta-binomial likelihood
nll_aep2 <- function(par) {
  getAll(par, dat)
  y <- OBS(y); size <- OBS(size)
  theta <- plogis(X_theta %*% beta_theta); ADREPORT(theta) # overdispersion parameter
  prob <- plogis(X %*% beta); ADREPORT(prob) # linear predictor and link
  - sum(dbetabinom(y, size, prob / theta, (1-prob) / theta, log = TRUE))
}

# Design matrices
X <- model.matrix(~ ward + loglos + year, data = aep)
X_theta <- model.matrix(~ year, data = aep)

# Initial parameters
beta <- c(-1, rep(0, ncol(X)-1)); names(beta) <- colnames(X)
beta_theta <- c(1, 0); names(beta_theta) <- colnames(X_theta)

par <- list(beta = beta, beta_theta = beta_theta)
dat <- list(
  y = aep$y[,1],
  size = aep$los,
  X = X, 
  X_theta = X_theta
)

obj_aep3 <- MakeADFun(nll_aep2, par, silent = TRUE)
opt_aep3 <- nlminb(obj_aep3$par, obj_aep3$fn, obj_aep3$gr)

sdr_aep3 <- sdreport(obj_aep3)

beta3 <- as.list(sdr_aep3, "Est")$beta

round(beta3, 3)

## ----cleanup4, include=FALSE--------------------------------------------------
TMB::FreeADFun(obj_aep1)
TMB::FreeADFun(obj_aep2)
TMB::FreeADFun(obj_aep3)

## ----faithful-----------------------------------------------------------------
data(faithful)

## ----copula likelihood--------------------------------------------------------
nll_copula <- function(par) {
  getAll(par, faithful)
  REPORT(mu1); REPORT(mu2)
  sigma1 <- exp(log_sigma1); REPORT(sigma1) # marginal sds component 1
  sigma2 <- exp(log_sigma2); REPORT(sigma2) # marginal sds component 2
  theta <- exp(log_theta); REPORT(theta) # dependence parameters
  alpha <- exp(log_alpha); REPORT(alpha) # mixture weights
  # Marginal densities
  # Margin 1: Waiting
  d1 <- cbind(dnorm(waiting, mu1[1], sigma1[1], log=TRUE), # Component 1
              dnorm(waiting, mu2[1], sigma2[1], log=TRUE)) # Component 2
  # Margin 2: Eruptions
  d2 <- cbind(dnorm(eruptions, mu1[2], sigma1[2], log=TRUE), # Component 1
              dnorm(eruptions, mu2[2], sigma2[2], log=TRUE)) # Component 2
  # Marginal CDFs
  # Margin 1: Waiting
  p1 <- cbind(pnorm(waiting, mu1[1], sigma1[1]), # Component 1
              pnorm(waiting, mu2[1], sigma2[1])) # Component 2
  # Margin 2: Eruptions
  p2 <- cbind(pnorm(eruptions, mu1[2], sigma1[2]), # component 1
              pnorm(eruptions, mu2[2], sigma2[2])) # component 2
  
  # Computing mixture likelihood:
  ll1 <- dcopula(d1[,1], d2[,1], p1[,1], p2[,1], cclayton(theta[1]), log=TRUE) # f1(x,y)
  ll2 <- dcopula(d1[,2], d2[,2], p1[,2], p2[,2], cclayton(theta[2]), log=TRUE) # f2(x,y)
  # alpha * f1(x,y) + (1-alpha) * f2(x,y) on log scale for each obervation
  ll <- logspace_add(log_alpha + ll1, log1p(-alpha) + ll2) 
  - sum(ll) # returning negative sum
}

## ----copula fit---------------------------------------------------------------
# Initial parameters
par <- list(
  mu1 = c(55, 2), mu2 = c(80, 4),
  log_sigma1 = log(c(10, 1)), log_sigma2 = log(c(10, 1)),
  log_theta = log(c(0.5, 0.5)),
  log_alpha = log(0.5)
)

obj_copula <- MakeADFun(nll_copula, par, silent = TRUE)
opt_copula <- nlminb(obj_copula$par, obj_copula$fn, obj_copula$gr)

mod_copula <- obj_copula$report()

# Extract transformed parameters
mu1    <- mod_copula$mu1
mu2    <- mod_copula$mu2
sigma1 <- mod_copula$sigma1
sigma2 <- mod_copula$sigma2
theta  <- mod_copula$theta
alpha  <- mod_copula$alpha

## ----copula results-----------------------------------------------------------
# Scatterplot
plot(faithful$waiting, faithful$eruptions, pch = 20, bty = "n",
     xlab = "Waiting time", ylab = "Eruption time", col = "#00000070")

# Grid for evaluation
xseq <- seq(min(faithful$waiting), max(faithful$waiting), length.out = 80)
yseq <- seq(min(faithful$eruptions), max(faithful$eruptions), length.out = 80)

# Evaluate component densities on grid
f1 <- outer(xseq, yseq, function(x,y){
  d1c1 <- dnorm(x, mu1[1], sigma1[1])
  d2c1 <- dnorm(y, mu1[2], sigma1[2])
  p1c1 <- pnorm(x, mu1[1], sigma1[1])
  p2c1 <- pnorm(y, mu1[2], sigma1[2])
  dcopula(d1c1, d2c1, p1c1, p2c1, cclayton(theta[1]))
})
f2 <- outer(xseq, yseq, function(x,y){
  d1c2 <- dnorm(x, mu2[1], sigma2[1])
  d2c2 <- dnorm(y, mu2[2], sigma2[2])
  p1c2 <- pnorm(x, mu2[1], sigma2[1])
  p2c2 <- pnorm(y, mu2[2], sigma2[2])
  dcopula(d1c2, d2c2, p1c2, p2c2, cclayton(theta[2]))
})

# Add contours
contour(xseq, yseq, alpha * f1, add = TRUE, nlevels = 8,
        drawlabels = FALSE, col = "orange", lwd = 2)
contour(xseq, yseq, (1-alpha) * f2, add = TRUE, nlevels = 8,
        drawlabels = FALSE, col = "deepskyblue", lwd = 2)

## ----cleanup5, include=FALSE--------------------------------------------------
TMB::FreeADFun(obj_copula)

## ----tape_config, include=FALSE-----------------------------------------------
# generalised to multivariate t distributed responses
TapeConfig(atomic="disable") ## Optional (speeds up this model)

## ----svt data-----------------------------------------------------------------
source("https://raw.githubusercontent.com/kaskr/RTMB/master/tmb_examples/sdv_multi_data.R")

## ----mvt sv-------------------------------------------------------------------
# Multivatiate SV model from Table 5 of Skaug and Yu "A flexible and automated likelihood based 
# framework for inference in stochastic volatility models." Computational Statistics & Data Analysis 76 (2014): 642-654.

## Parameter initial guess
par <- list(
  logit_phi  = qlogis(rep(0.97,p)),       # See eqn (12) in Skaug and Yu (2014)
  log_sigma  = log(rep(0.2,p)),           #       ---------||---------
  mu_y       = rep(-0.5,p),               #       ---------||---------
  off_diag_x = rep(0,p),                  #       ---------||---------
  h          = matrix(0,nrow=n,ncol=p),   #       ---------||---------
  log_df     = log(20)                    # this allows for heavier tails
)

# Negative joint likelihood (nll) of data and parameters
nll_svt <- function(par) {
  getAll(par)
  # Optionally mark the observation object
  y <- OBS(y)
  
  # Parameters on natural scale
  sigma <- exp(log_sigma)
  phi <- plogis(logit_phi)
  sigma_init <- sigma / sqrt(1-phi^2)
  df <- exp(log_df)
  
  nll <- 0  # Start collecting contributions
  
  # Prior on state variable (h)
  nll <- nll - sum(dnorm(h[1,], 0, sigma_init, log=TRUE))   # Start in stationary distribution
  for(j in 1:p) {
    nll <- nll - sum(dnorm(h[-1,j], phi[j]*h[-n,j], sigma[j], log=TRUE))  # AR(1) process for each dimension
  }
  
  # Parameterizes correlation matrix of X in terms of Cholesky factor
  L <- diag(p)
  L[lower.tri(L)] <- off_diag_x   
  row_norms <- apply(L, 1, function(row) sqrt(sum(row^2)))
  L <- L / row_norms
  R <- L %*% t(L)  # Correlation matrix of X (guarantied positive definite)
  
  # Likelihood of data y given h
  for(i in 1:n){
    sigma_y <- exp(0.5 * (mu_y + h[i,])) # variance depending on state
    Cov <- diag(sigma_y) %*% R %*% diag(sigma_y)
    nll <- nll - sum(dmvt(y[i,], 0, Cov, df, log = TRUE))
  }
  nll
}

obj_svt <- MakeADFun(nll_svt, par, random="h", silent = TRUE)
system.time(
  opt_svt <- nlminb(obj_svt$par, obj_svt$fn, obj_svt$gr)
)
rep <- sdreport(obj_svt)
rep

