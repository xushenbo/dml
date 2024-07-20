library(dplyr)
library(ggplot2)
library(survival)
library(pec)
library(data.table)
library(riskRegression)
library(DescTools)
library(zoo)
library(foreach)
library(doParallel)
library(simsurv)
library(mvtnorm)
library(mstate)
library(MASS)

ntrue <- 1e+06
admin.cens <- 10
taus <- c(1,2,3,4)

dgp4_true <- function(ntrue){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=ntrue, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(ntrue,1,0.5)
  x5 <- rbinom(ntrue,1,0.5)
  x6 <- rbinom(ntrue,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # -1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6
  a <- rbinom(ntrue, 1, ps) # -2+1.2*x1+1.6*x2+1.6*x3+1.6*x4, -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  # ps <- Winsorize(minval = 1e-4, maxval = 1-1e-4, ps)
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  data.df$naive <- 1
  data.df$ow.tilt <- with(data.df, ps*(1-ps))
  data.df$mw.tilt <- with(data.df, pmin(ps, 1-ps))
  data.df$enw.tilt <- with(data.df, ifelse(ps==0 | ps==1, 0, -(ps*log(ps)+(1-ps)*log(1-ps))))
  
  data.df$Tj1a0 <- ceiling(simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime*1000)/1000
  data.df$Tj1a1 <- ceiling(simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime*1000)/1000
  data.df$Tj2a0 <- ceiling(simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime*1000)/1000
  data.df$Tj2a1 <- ceiling(simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime*1000)/1000
  
  data.df$Ta0 <- pmin(data.df$Tj1a0,data.df$Tj2a0) # round(pmin(data.df$Tj1a0,data.df$Tj2a0),digits=3)
  data.df$Ta1 <- pmin(data.df$Tj1a1,data.df$Tj2a1) # round(pmin(data.df$Tj1a1,data.df$Tj2a1),digits=3)
  data.df$Da0 <- with(data.df,ifelse(Tj1a0<Tj2a0,1,2))
  data.df$Da1 <- with(data.df,ifelse(Tj1a1<Tj2a1,1,2))
  return(data.df)
}

start.time <- Sys.time()
large.df <- dgp4_true(ntrue)
end.time <- Sys.time()
end.time - start.time

large.df$Ta0[large.df$Ta0>taus[4]+1] <- taus[4]+1
large.df$Ta1[large.df$Ta1>taus[4]+1] <- taus[4]+1

compute_true <- function(time, event, tilt)
{
  n <- length(time)
  cause <- 1
  s <- sort(unique(time))
  ns <- length(s)
  ds <- diff(c(0, s))
  return(1/mean(tilt)*colMeans(matrix(tilt, ncol=ns, nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) (s[u]-pmin(time, s[u]))*(event == cause)))))
}

start.time <- Sys.time()
true.a0 <- data.frame(s=sort(unique(large.df$Ta0)))
true.a0$rmtlj.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$naive)
true.a0$rmtlj.ow.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$ow.tilt)
true.a0$rmtlj.mw.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$mw.tilt)
true.a0$rmtlj.enw.a0 <- compute_true(time=large.df$Ta0, event=large.df$Da0, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
true.a1 <- data.frame(s=sort(unique(large.df$Ta1)))
true.a1$rmtlj.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$naive)
true.a1$rmtlj.ow.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$ow.tilt)
true.a1$rmtlj.mw.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$mw.tilt)
true.a1$rmtlj.enw.a1 <- compute_true(time=large.df$Ta1, event=large.df$Da1, tilt=large.df$enw.tilt)
end.time <- Sys.time()
end.time - start.time

find_value <- function(true.a0, true.a1, tau)
{
  true.value <- data.frame(tau, rmtlj.a0=NA, rmtlj.ow.a0=NA, rmtlj.mw.a0=NA, rmtlj.enw.a0=NA,
                           rmtlj.a1=NA, rmtlj.ow.a1=NA, rmtlj.mw.a1=NA, rmtlj.enw.a1=NA,
                           rmtlj.diff=NA, rmtlj.ow.diff=NA, rmtlj.mw.diff=NA, rmtlj.enw.diff=NA)
  # a=0
  ind.a0 <- max(which(true.a0$s<=tau))
  true.value$rmtlj.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.ow.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.ow.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.mw.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.mw.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  true.value$rmtlj.enw.a0 <- as.numeric(approx(x=true.a0$s[ind.a0:(ind.a0+1)], y=true.a0$rmtlj.enw.a0[ind.a0:(ind.a0+1)], xout=tau))[2]
  
  # a=1
  ind.a1 <- max(which(true.a1$s<=tau))
  true.value$rmtlj.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.ow.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.ow.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.mw.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.mw.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  true.value$rmtlj.enw.a1 <- as.numeric(approx(x=true.a1$s[ind.a1:(ind.a1+1)], y=true.a1$rmtlj.enw.a1[ind.a1:(ind.a1+1)], xout=tau))[2]
  
  # diff
  true.value$rmtlj.diff <- true.value$rmtlj.a1-true.value$rmtlj.a0
  true.value$rmtlj.ow.diff <- true.value$rmtlj.ow.a1-true.value$rmtlj.ow.a0
  true.value$rmtlj.mw.diff <- true.value$rmtlj.mw.a1-true.value$rmtlj.mw.a0
  true.value$rmtlj.enw.diff <- true.value$rmtlj.enw.a1-true.value$rmtlj.enw.a0
  return(true.value)
}

true.value.tau1 <- find_value(true.a0, true.a1, taus[1])
true.value.tau2 <- find_value(true.a0, true.a1, taus[2])
true.value.tau3 <- find_value(true.a0, true.a1, taus[3])
true.value.tau4 <- find_value(true.a0, true.a1, taus[4])

save(true.value.tau1, true.value.tau2, true.value.tau3, true.value.tau4, file="true_competing_7.RData")
