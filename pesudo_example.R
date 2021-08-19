###
### R code for the sample size calculations using simulation (written by Gary Collins)
###
### Riley RD, Collins GS, Ensor J, Archer L, Booth S, Mozumder S, Rutherford M,
### van Smeden M, Lambert P, Snell K. "Minimum sample size calculations for 
### external validation of a clinical prediction model with a time-to-event outcome"
### (under review)

#here::here()
require(openxlsx)
require(progress)
require(pseudo)
require(simsurv)
require(survival)
require(pec)
require(prodlim)
require(rms)
require(mfp)

dev.off(dev.list()["RStudioGD"]) # remove any plot windows

### (c) work out the rate parameter needed to get an S(3) of 0.83 as 
### expected for validation sample
N       <- 5000
# specify LP distribution from the developed prediction model
moments <- c(mean = 4.60, variance = 0.65, skewness = -0.5, kurtosis = 5)
LP      <- data.frame(id = 1:N, value = PearsonDS::rpearson(N, moments = moments))

rate_parameter <- function(lambda, LP = LP, time.point = 3){
  X.sim <- simsurv(dist = 'exponential', lambdas = lambda, x = LP, betas = c(value = 1))
  X.sim$eventtime <- X.sim$eventtime
  X.sim$dead <- rep(1, length(LP$value))
  X.sim$dead[X.sim$eventtime > time.point] <- 0
  X.sim$eventtime[X.sim$eventtime > time.point] <- time.point
  X.sim$LP <- LP$value
  
  fit.survfit <- survfit(Surv(eventtime, dead)~1, data = X.sim)
  St          <- summary(fit.survfit, times = time.point)$surv
  fit.coxph   <- coxph(Surv(eventtime, dead)~LP, data = X.sim)
  c.index     <- as.numeric(fit.coxph$concordance[6])
  D           <- as.numeric(royston(fit.coxph)[1])
  return(list(St = St, c.index = c.index, D = D))
}
# trial and error finds that lambda~0.00055 for S(3) = 0.83
# leads to a c-index of ~0.69, and a D~1.19

### small simulation to get lambda to account for sampling variability in survsim
N.SIM       <- 100
N           <- 5000
lambda.grid <- seq(0.00025, 0.00075, length = 20)
St          <- matrix(ncol = length(lambda.grid), nrow = N.SIM)
c.index     <- matrix(ncol = length(lambda.grid), nrow = N.SIM)
D           <- matrix(ncol = length(lambda.grid), nrow = N.SIM)

pb <- progress_bar$new(format = "  simulation :what [:bar] :percent eta: :eta",
                       clear = FALSE, 
                       total = N.SIM * length(lambda.grid), 
                       width = 60)

## generate large data with LP values from skew normal distribution (VTE example)
moments <- c(mean = 4.60, variance = 0.65, skewness = -0.5, kurtosis = 5)

for(j in 1:length(lambda.grid)){  ## This is slow
  for(i in 1:N.SIM){
      pb$tick()
      LP           <- data.frame(id=1:N, value = PearsonDS::rpearson(N, moments = moments))
      OUT          <- rate_parameter(lambda = lambda.grid[j], LP = LP, time.point = 3)
      St[i,j]      <- OUT$St
      c.index[i,j] <- OUT$c.index
      D[i,j]       <- OUT$D
    }
}
apply(St, 2, mean)
apply(c.index, 2, mean)
apply(D, 2, mean)

plot(lambda.grid, colMeans(St), type = 'b', pch = 20, xlab = 'lambda_c', ylab = 'S(3)')
grid()
# get lambdas for S(3) = 0.83
approx(x = colMeans(St), y = lambda.grid, xout = c(0.83))
approx(x = lambda.grid, y = colMeans(c.index), xout = approx(x = colMeans(St), y = lambda.grid, xout = c(0.83))$y)
approx(x = lambda.grid, y = colMeans(D), xout = approx(x = colMeans(St), y = lambda.grid, xout = c(0.83))$y)

  
### (d) checking censoring rate in VTE model development dataset
### censoring was high in the development dataset
### 864 (72%) out of 1200 participants were censored before 3 years in the VTE dataset

### small simulation to account for sampling variability in survsim
N.SIM <- 200
N.OBS <- 10000
lambda.grid <- seq(0.3, 0.6, 0.05)
cens <- matrix(ncol = length(lambda.grid), nrow = N.SIM)
pb <- progress_bar$new(format = "  simulation :what [:bar] :percent eta: :eta",
                       clear = FALSE, 
                       total = N.SIM * length(lambda.grid), 
                       width = 60)

for(j in 1: length(lambda.grid)){
  for(i in 1:N.SIM){
    pb$tick()
    X.sim                            <- simsurv(dist = 'exponential', lambdas = lambda.grid[j], x = data.frame(ids = seq(1:N.OBS)))
    X.sim$cens                       <- rep(1, nrow(X.sim))
    X.sim$cens[X.sim$eventtime > 3]  <- 0
    X.sim$eventtime                  <- X.sim$eventtime
    X.sim$eventtime[X.sim$eventtime > 3] <- 3
    fit.survfit                      <- survfit(Surv(eventtime, status)~1, data = X.sim)
    cens[i, j]                       <- summary(fit.survfit, times = 2.99999999)$surv
  }
}
apply(cens, 2, mean)
plot(lambda.grid, colMeans(cens), type = 'b', pch = 20, xlab = 'lambda_c', ylab = '% censored')
grid()
# get lambdas for 18%, 28%, 38% change of not being censored
approx(x = colMeans(cens), y = lambda.grid, xout = c(0.18, 0.28, 0.38)) # 0.572, 0.426, 0.323


### PART B: Sample size calculation #####
### STEPS 2 to 9: simulate data according to a particular size 
### then estimate calibration slope and its standard error
### then repeat over many (e.g. 1000) simulation 
### Run all this code in one go from START till END 
### START 
rm(list = ls())
set.seed(123456)
N          <- 100000
N.VAL      <- 3600
N.SIM      <- 100
time.point <- 3
se.slope   <- vector(mode = 'numeric', length = N.SIM)
S3         <- vector(mode = 'numeric', length = N.SIM)
c.index    <- matrix(ncol = 2, nrow = N.SIM)
moments    <- c(mean = 4.60, variance = 0.65, skewness = -0.5, kurtosis = 5)
LP         <- data.frame(id=1:N, value = PearsonDS::rpearson(N, moments = moments))
X.sim      <- simsurv(dist = 'exponential', lambdas = 0.00051, x = LP, betas = c(value = 1))
X.sim$LP   <- LP$value
rm(LP)

# censor at max follow-up time, say 3 years
X.sim$status <- rep(1, nrow(X.sim))
X.sim$status[X.sim$eventtime > 3]    <- 0
X.sim$eventtime[X.sim$eventtime > 3] <- 3

# generate censoring times from the exp dist identified ealier
X.cens <- simsurv(dist = 'exponential', lambdas = 0.418, x = data.frame(ids=seq(1:nrow(X.sim))))
X.sim$status[X.cens$eventtime < X.sim$eventtime] <- 0
X.sim$eventtime[X.cens$eventtime < X.sim$eventtime] <- X.cens$eventtime[X.cens$eventtime < X.sim$eventtime]

# set up the plot for the calibration curves
plot(c(0,1), c(0,1), type= 'n', xlab = 'predicted risk', ylab = 'observed risk')
  
pb <- progress_bar$new(
  format = "  simulation :what [:bar] :percent eta: :eta",
  clear = FALSE, total = N.SIM, width = 60)

for(i in 1:N.SIM){
  pb$tick()
  X.sim$val <- rep(0, nrow(X.sim))
  index <- sample(1:nrow(X.sim), size = N.VAL, replace = F)
  X.sim$val[index] <- 1
  fit <- coxph(Surv(eventtime, status)~LP, data = X.sim[X.sim$val == 0, ])

  S3[i]      <- summary(survfit(fit), times = time.point)$surv
  pred.val   <- 1 - S3[i]^predict(fit, newdata = X.sim[X.sim$val == 1, ], type = 'risk')
  f.val      <- prodlim(Hist(eventtime, status)~1, data = X.sim[X.sim$val == 1, ])
  pseudo.val <- jackknife(f.val, times = time.point)

  c.index[i,] <- rcorr.cens(1-pred.val, S=Surv(X.sim$eventtime[X.sim$val==1], X.sim$status[X.sim$val==1]))[c(1, 3)]
  c.index[i, 2] <- c.index[i, 2] / 2

  X <- matrix(ncol=2, nrow = length(pred.val))
  X[,1] <- pred.val
  X[,2] <- pseudo.val
  X <- X[order(X[,1]),]
  # Plot calibration curves at each iteration
  lines(X[, 1], 1 - predict(loess(X[,2]~X[,1])), col = 'grey')
  
  pv.val <- 1 - pseudo.val ## pseudo values are generated for survival (so 1 minus for the event)
  mv.val <- pred.val
  xx.val <- log(-log(1 - mv.val))

  # Fit calibration model
  calfit.val.2 <- geese(pv.val~xx.val, 
                        jack      = T, 
                        scale.fix = T, 
                        id        = 1:length(xx.val), 
                        family    = gaussian, 
                        mean.link = 'cloglog', 
                        corstr    = 'independence')

  se.slope[i] <- summary(calfit.val.2)$mean[2,2] # pull out robust se
}
abline(a = 0, b = 1, col = 'black', lwd = 2)

summary(se.slope)
plot(cumsum(se.slope)/seq_along(se.slope), pch = 20, ylab='running mean', xlab='iteration')
colMeans(c.index) 

### based on 200 simulations
# SS of  3600 mean SE is 0.10292 (mean ci width of the c index is 0.01425576*1.96*2=0.056)
# SS of  5000 mean SE is 0.08689
# SS of  7000 mean SE is 0.07330
# SS of  9000 mean SE is 0.06482
# SS of 12000 mean SE is 0.05602
# SS of 13000 mean SE is 0.05365
# SS of 14000 mean SE is 0.05186


############ PART 2: Sample size calculation based on D statistic ############
  
# PART A: find the right distribution to simulate survival times
# using D statistic of 1.115 ***
# then SD of LP is 1.115 / sqrt(8/pi) = 0.699 
# let us assume mean is 0, so centered
# then need to find rate of the baseline to get an S(3) of 0.83

N        <- 100000
# specify LP distribution from the developed prediction model
LP       <- rnorm(N, 0, 0.699)
X.sim    <- simsurv(dist = 'exponential', lambdas = 0.052, x = data.frame(value=LP), betas = c(value = 1))
X.sim$LP <- LP
rm(LP)

X.sim$status <- rep(1, nrow(X.sim))
X.sim$status[X.sim$eventtime > 3]    <- 0
X.sim$eventtime[X.sim$eventtime > 3] <- 3

# generate censoring times from the exp dist identified ealier
X.cens <- simsurv(dist = 'exponential', lambdas = 0.418, x = data.frame(ids=seq(1:nrow(X.sim))))
X.sim$status[X.cens$eventtime < X.sim$eventtime] <- 0
X.sim$eventtime[X.cens$eventtime < X.sim$eventtime] <- X.cens$eventtime[X.cens$eventtime < X.sim$eventtime]

N.VAL <- 3400
N.SIM <- 50
se.slope   <- vector(mode = 'numeric', length = N.SIM)
for(i in 1:N.SIM){
  X.sim$val <- rep(0, nrow(X.sim))
  index <- sample(1:nrow(X.sim), size = N.VAL, replace = F)
  X.sim$val[index] <- 1
  fit <- coxph(Surv(eventtime, status)~LP, data = X.sim[X.sim$val == 0, ])
  
  S3[i]      <- summary(survfit(fit), times = time.point)$surv
  pred.val   <- 1 - S3[i]^predict(fit, newdata = X.sim[X.sim$val == 1, ], type = 'risk')
  f.val      <- prodlim(Hist(eventtime, status)~1, data = X.sim[X.sim$val == 1, ])
  pseudo.val <- jackknife(f.val, times = time.point)
  
  pv.val <- 1 - pseudo.val ## pseudo values are generated for survival (so 1 minus for the event)
  
  mv.val <- pred.val
  xx.val <- log(-log(1 - mv.val))
  
  # Fit calibration model
  calfit.val.2 <- geese(pv.val~xx.val, 
                        jack      = T, 
                        scale.fix = T, 
                        id        = 1:length(xx.val), 
                        family    = gaussian, 
                        mean.link = 'cloglog', 
                        corstr    = 'independence')
  
  se.slope[i] <- summary(calfit.val.2)$mean[2,2] # pull out robust se
}
summary(se.slope)

# SS of 5000 (513 events) gives mean SE is 0.08738
# SS of 4000 (403 events) gives mean SE of 0.09745
# SS of 3400 (366 events) gives mean SE of 0.10790
