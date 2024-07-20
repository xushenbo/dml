# library(survRM2)
# library(ggplot2)
# library(ggsci)
# library(fastDummies)
# library(gee)
# library(geepack)
# library(surv2sampleComp)
# library(relsurv)
# library(fastpseudo)
# library(data.table)
# library(tmle) # tmle3, drtmle
# library(pseudo)
# library(eventglm)
# library(riskRegression)
# library(causalCmprsk)
# 
# library(drgee)
# library(PSweight)
# library(xtable)
# library(mstate)

##
library(foreach)
library(doParallel)
library(mvtnorm)
library(simsurv)
library(pec)
library(DescTools)
library(zoo)
library(dplyr)
library(survival)
library(quantmod)
library(tidyr)
library(fdrtool)
library(simecol)
library(MASS)

registerDoParallel(cores=100)

nsim <- 1000
n <- 2000
nboot <- 100
# ntrue <- 1e+06

# nsim <- 10
# n <- 1000
# ntrue <- 1000
# nboot <- 10

time.grid <- 0.01
admin.cens <- 10
taus <- c(1,2,3,4)
cens.prob <- data.frame(matrix(NA, nrow=nsim, ncol=4))
colnames(cens.prob) <- c("tau1","tau2","tau3","tau4") # ,"tau6"
# counting.sl.lib <- c("SL.mean","SL.glm","SL.step","SL.glmnet","SL.rpartPrune","SL.gam","SL.xgboost","SL.ranger") # ,"SL.nnet","SL.kernelKnn","SL.earth" # "SL.xgboost":33s; "SL.ranger":30s; "SL.ksvm":48.19s; "SL.bartMachine":slow; "SL.polymars": not working
# a.sl.lib <- c("SL.mean","SL.nnet","SL.kernelKnn","SL.rpartPrune","SL.xgboost","SL.ranger") # ,"SL.glm","SL.step","SL.gam","SL.glmnet","SL.earth"
# event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.loglogreg","survSL.pchreg","survSL.rfsrc","survSL.coxph","survSL.expreg","survSL.weibreg") # ; "survSL.gam","survSL.pchSL"
confounders <- c("x1","x2","x3","x4","x5","x6")

# outcome regression rmtlj by integrating survival
or_rmtlj_integral_point <- function(id, a, time, event, tilt, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  
  estimate.a0 <- 1/mean(tilt)*cumsum(cummax(colMeans(term2.a0))*ds)
  estimate.a1 <- 1/mean(tilt)*cumsum(cummax(colMeans(term2.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# ipcw rmtlj by integrating survival
ipcw_rmtlj_integral <- function(id, a, time, event, bw, G.a0, G.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Njt <- do.call(cbind, lapply(1:ns, function(u)  (time <= s[u])*(event == cause)))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  estimate.a0 <- 1/mean(bw*(a==0))*cumsum(cummax(colMeans(term1.a0))*ds)
  estimate.a1 <- 1/mean(bw*(a==1))*cumsum(cummax(colMeans(term1.a1))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# dr rmtlj by integrating survival
dr_rmtlj_integral <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, Sj.a0, Sj.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  cause <- 1
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  Sj.a0 <- t(na.locf(t(ifelse(Sj.a0 < 1e-3, 1e-3, Sj.a0))))
  Sj.a1 <- t(na.locf(t(ifelse(Sj.a1 < 1e-3, 1e-3, Sj.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  Fj.dHazard.a0 <- t(apply(cbind(0, -log(Sj.a0)), 1, diff))
  Fj.dHazard.a1 <- t(apply(cbind(0, -log(Sj.a1)), 1, diff))
  Fj.a0 <- t(apply(cbind(1, S.a0[, 1:(ns-1)])*Fj.dHazard.a0, 1, cumsum))
  Fj.a1 <- t(apply(cbind(1, S.a1[, 1:(ns-1)])*Fj.dHazard.a1, 1, cumsum))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  Njt <- do.call(cbind, lapply(1:ns, function(u)  (time <= s[u])*(event == cause)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Fj.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a0)*Fj.a0/(S.a0*G.a0), 1, cumsum))
  estimate.a0 <- cumsum((cummax(1/mean(bw*(a==0))*colMeans(term1.a0+term3.a0-term4.a0)+1/mean(tilt)*colMeans(term2.a0)))*ds)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Njt/matrix(rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*Fj.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Fj.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply((dNct - Yt*G.dHazard.a1)*Fj.a1/(S.a1*G.a1), 1, cumsum))
  estimate.a1 <- cumsum((cummax(1/mean(bw*(a==1))*colMeans(term1.a1+term3.a1-term4.a1)+1/mean(tilt)*colMeans(term2.a1)))*ds)
  return(data.frame(s, estimate.a0, estimate.a1))
}

# lognormal <- function(data.df, betas, sigma, n)
# {return(exp(sigma*qnorm(runif(n))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# loglogistic <- function(data.df, betas, gamma, n)
# {return(exp(gamma*(log(runif(n))-log(1-runif(n)))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# aft <- function(data.df, betas, sigma, n)
# {return(exp(sigma*rnorm(n, mean=0, sd=1)+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

## proportional hazards
dgp4_competing <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(0.3+0.2*x1+0.3*x2+0.3*x3-0.2*x4-0.3*x5-0.2*x6))) # 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  
  data.df$Tj1a0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Tj1a1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Tj1 <- (1-a)*data.df$Tj1a0 + a*data.df$Tj1a1
  
  data.df$Tj2a0 <- simsurv(dist = "exponential", lambdas = 0.1, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.12, x1=-0.1, x2=0.3, x3=0.1, x4=0.2, x5=-0.4, x6=0.5))$eventtime
  data.df$Tj2a1 <- simsurv(dist = "exponential", lambdas = 0.08, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=-0.2, x2=-0.1, x3=0.2, x4=0.3, x5=0.3, x6=-0.3))$eventtime
  data.df$Tj2 <- (1-a)*data.df$Tj2a0 + a*data.df$Tj2a1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.12, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.14, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Tj1, data.df$Tj2)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df, ifelse(pmin(Tj1,Tj2)>C,0, ifelse(Tj1<Tj2,1,2)))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

extract_point <- function(estimate.df, tau){
  estimate.a0 <- approx(estimate.df$s, estimate.df$estimate.a0, xout=tau)$y
  estimate.a1 <- approx(estimate.df$s, estimate.df$estimate.a1, xout=tau)$y
  estimate.diff <- estimate.a1-estimate.a0
  return(data.frame(estimate.a0, estimate.a1, estimate.diff))
}

f.sim.result <- data.frame(matrix(NA,nrow=14,ncol=0))

start.time <- Sys.time()
for(i in 1:nsim) {
  sim.df <- dgp4_competing(2*n, admin.cens, time.grid)
  
  while(max(sim.df$end.time[sim.df$a==0 & sim.df$event==0])<taus[4] | max(sim.df$end.time[sim.df$a==1 & sim.df$event==0])<taus[4] |
        max(sim.df$end.time[sim.df$a==0 & sim.df$event==1])<taus[4] | max(sim.df$end.time[sim.df$a==1 & sim.df$event==1])<taus[4]) {
    sim.df <- dgp4_competing(2*n, admin.cens, time.grid)}
  
  sim.df$naive <- 1
  sim.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict(glm(a~x1+x2+x3+x4+x5+x6, data=sim.df, family=binomial()), type="response", newdata=sim.df)) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
  sim.df$iptw <- with(sim.df,1/(ps*a+(1-ps)*(1-a)))
  sim.df$iptw.trunc <- with(sim.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
  sim.df$ow <- with(sim.df,ifelse(a,1-ps,ps))
  sim.df$mw <- with(sim.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
  sim.df$enw <- with(sim.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
  sim.df$ow.tilt <- with(sim.df, ps*(1-ps))
  sim.df$mw.tilt <- with(sim.df, pmin(ps, 1-ps))
  sim.df$enw.tilt <- with(sim.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
  
  # conditional survival
  S.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event!=0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  S.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event!=0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  
  # conditional censoring
  G.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  G.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  
  # conditional competing
  Sj.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==0,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  Sj.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sim.df[sim.df$a==1,], x=TRUE), newdata=sim.df, times=sort(unique(sim.df$end.time)))
  
  f.sim <- list("or.integral"=NA, "or.integral.ow"=NA, "or.integral.mw"=NA, "or.integral.enw"=NA,
                  "ipcw.integral.iptw"=NA, "ipcw.integral.iptw.trunc"=NA, "ipcw.integral.ow"=NA, "ipcw.integral.mw"=NA, "ipcw.integral.enw"=NA,
                  "dr.integral.driptw"=NA, "dr.integral.driptw.trunc"=NA, "dr.integral.drow"=NA, "dr.integral.drmw"=NA, "dr.integral.drenw"=NA)
  
  f.sim$or.integral <- or_rmtlj_integral_point(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$naive, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$or.integral.ow <- or_rmtlj_integral_point(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$ow.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$or.integral.mw <- or_rmtlj_integral_point(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$mw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$or.integral.enw <- or_rmtlj_integral_point(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, tilt=sim.df$enw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  
  f.sim$ipcw.integral.iptw <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw, G.a0=G.a0, G.a1=G.a1)
  f.sim$ipcw.integral.iptw.trunc <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw.trunc, G.a0=G.a0, G.a1=G.a1)
  f.sim$ipcw.integral.ow <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$ow, G.a0=G.a0, G.a1=G.a1)
  f.sim$ipcw.integral.mw <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$mw, G.a0=G.a0, G.a1=G.a1)
  f.sim$ipcw.integral.enw <- ipcw_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$enw, G.a0=G.a0, G.a1=G.a1)
  
  f.sim$dr.integral.driptw <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw, tilt=sim.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$dr.integral.driptw.trunc <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$iptw.trunc, tilt=sim.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$dr.integral.drow <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$ow, tilt=sim.df$ow.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$dr.integral.drmw <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$mw, tilt=sim.df$mw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  f.sim$dr.integral.drenw <- dr_rmtlj_integral(id=sim.df$id, a=sim.df$a, time=sim.df$end.time, event=sim.df$event, bw=sim.df$enw, tilt=sim.df$enw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
  
  f.sim.point <- data.frame(matrix(NA, nrow=length(f.sim), ncol=3))
  colnames(f.sim.point) <- c("estimate.a0", "estimate.a1", "estimate.diff")
  rownames(f.sim.point) <- names(f.sim)
  
  f.sim.point["or.integral", ] <- extract_point(f.sim$or.integral, taus[4])
  f.sim.point["or.integral.ow", ] <- extract_point(f.sim$or.integral.ow, taus[4])
  f.sim.point["or.integral.mw", ] <- extract_point(f.sim$or.integral.mw, taus[4])
  f.sim.point["or.integral.enw", ] <- extract_point(f.sim$or.integral.enw, taus[4])
  
  f.sim.point["ipcw.integral.iptw", ] <- extract_point(f.sim$ipcw.integral.iptw, taus[4])
  f.sim.point["ipcw.integral.iptw.trunc", ] <- extract_point(f.sim$ipcw.integral.iptw.trunc, taus[4])
  f.sim.point["ipcw.integral.ow", ] <- extract_point(f.sim$ipcw.integral.ow, taus[4])
  f.sim.point["ipcw.integral.mw", ] <- extract_point(f.sim$ipcw.integral.mw, taus[4])
  f.sim.point["ipcw.integral.enw", ] <- extract_point(f.sim$ipcw.integral.enw, taus[4])
  
  f.sim.point["dr.integral.driptw", ] <- extract_point(f.sim$dr.integral.driptw, taus[4])
  f.sim.point["dr.integral.driptw.trunc", ] <- extract_point(f.sim$dr.integral.driptw.trunc, taus[4])
  f.sim.point["dr.integral.drow", ] <- extract_point(f.sim$dr.integral.drow, taus[4])
  f.sim.point["dr.integral.drmw", ] <- extract_point(f.sim$dr.integral.drmw, taus[4])
  f.sim.point["dr.integral.drenw", ] <- extract_point(f.sim$dr.integral.drenw, taus[4])
  
  f.sim.boot <- foreach(j=1:nboot,.combine=cbind) %dopar% {
    sample.df <- sim.df[sample(1:nrow(sim.df), size = nrow(sim.df), replace = TRUE),]
    
    while(max(sample.df$end.time[sample.df$a==0 & sample.df$event==0])<taus[4] | max(sample.df$end.time[sample.df$a==1 & sample.df$event==0])<taus[4] |
          max(sample.df$end.time[sample.df$a==0 & sample.df$event==1])<taus[4] | max(sample.df$end.time[sample.df$a==1 & sample.df$event==1])<taus[4]) {
      sample.df <- sim.df[sample(1:nrow(sim.df), size = nrow(sim.df), replace = TRUE),]}
    
    sample.df$naive <- 1
    sample.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict(glm(a~x1+x2+x3+x4+x5+x6, data=sample.df, family=binomial()), type="response", newdata=sample.df)) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
    sample.df$iptw <- with(sample.df,1/(ps*a+(1-ps)*(1-a)))
    sample.df$iptw.trunc <- with(sample.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
    sample.df$ow <- with(sample.df,ifelse(a,1-ps,ps))
    sample.df$mw <- with(sample.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
    sample.df$enw <- with(sample.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
    sample.df$ow.tilt <- with(sample.df, ps*(1-ps))
    sample.df$mw.tilt <- with(sample.df, pmin(ps, 1-ps))
    sample.df$enw.tilt <- with(sample.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
    
    # conditional survival
    S.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    S.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    
    # conditional censoring
    G.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    G.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==0))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    
    # conditional competing
    Sj.a0 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==0,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    Sj.a1 <- predictSurvProb(coxph(Surv(end.time, 1*(event==1))~x1+x2+x3+x4+x5+x6, data=sample.df[sample.df$a==1,], x=TRUE), newdata=sample.df, times=sort(unique(sample.df$end.time)))
    
    f.sim <- list("or.integral"=NA, "or.integral.ow"=NA, "or.integral.mw"=NA, "or.integral.enw"=NA,
                  "ipcw.integral.iptw"=NA, "ipcw.integral.iptw.trunc"=NA, "ipcw.integral.ow"=NA, "ipcw.integral.mw"=NA, "ipcw.integral.enw"=NA,
                  "dr.integral.driptw"=NA, "dr.integral.driptw.trunc"=NA, "dr.integral.drow"=NA, "dr.integral.drmw"=NA, "dr.integral.drenw"=NA)
    
    f.sim$or.integral <- or_rmtlj_integral_point(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$naive, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$or.integral.ow <- or_rmtlj_integral_point(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$ow.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$or.integral.mw <- or_rmtlj_integral_point(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$mw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$or.integral.enw <- or_rmtlj_integral_point(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, tilt=sample.df$enw.tilt, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    
    f.sim$ipcw.integral.iptw <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw, G.a0=G.a0, G.a1=G.a1)
    f.sim$ipcw.integral.iptw.trunc <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw.trunc, G.a0=G.a0, G.a1=G.a1)
    f.sim$ipcw.integral.ow <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$ow, G.a0=G.a0, G.a1=G.a1)
    f.sim$ipcw.integral.mw <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$mw, G.a0=G.a0, G.a1=G.a1)
    f.sim$ipcw.integral.enw <- ipcw_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$enw, G.a0=G.a0, G.a1=G.a1)
    
    f.sim$dr.integral.driptw <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw, tilt=sample.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$dr.integral.driptw.trunc <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$iptw.trunc, tilt=sample.df$naive, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$dr.integral.drow <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$ow, tilt=sample.df$ow.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$dr.integral.drmw <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$mw, tilt=sample.df$mw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    f.sim$dr.integral.drenw <- dr_rmtlj_integral(id=sample.df$id, a=sample.df$a, time=sample.df$end.time, event=sample.df$event, bw=sample.df$enw, tilt=sample.df$enw.tilt, G.a0=G.a0, G.a1=G.a1, S.a0=S.a0, S.a1=S.a1, Sj.a0=Sj.a0, Sj.a1=Sj.a1)
    
    f.sim.point <- data.frame(matrix(NA, nrow=length(f.sim), ncol=3))
    colnames(f.sim.point) <- c(paste0("estimate.a0.v",j), paste0("estimate.a1.v",j), paste0("estimate.diff.v",j))
    rownames(f.sim.point) <- names(f.sim)
    
    f.sim.point["or.integral", ] <- extract_point(f.sim$or.integral, taus[4])
    f.sim.point["or.integral.ow", ] <- extract_point(f.sim$or.integral.ow, taus[4])
    f.sim.point["or.integral.mw", ] <- extract_point(f.sim$or.integral.mw, taus[4])
    f.sim.point["or.integral.enw", ] <- extract_point(f.sim$or.integral.enw, taus[4])
    
    f.sim.point["ipcw.integral.iptw", ] <- extract_point(f.sim$ipcw.integral.iptw, taus[4])
    f.sim.point["ipcw.integral.iptw.trunc", ] <- extract_point(f.sim$ipcw.integral.iptw.trunc, taus[4])
    f.sim.point["ipcw.integral.ow", ] <- extract_point(f.sim$ipcw.integral.ow, taus[4])
    f.sim.point["ipcw.integral.mw", ] <- extract_point(f.sim$ipcw.integral.mw, taus[4])
    f.sim.point["ipcw.integral.enw", ] <- extract_point(f.sim$ipcw.integral.enw, taus[4])
    
    f.sim.point["dr.integral.driptw", ] <- extract_point(f.sim$dr.integral.driptw, taus[4])
    f.sim.point["dr.integral.driptw.trunc", ] <- extract_point(f.sim$dr.integral.driptw.trunc, taus[4])
    f.sim.point["dr.integral.drow", ] <- extract_point(f.sim$dr.integral.drow, taus[4])
    f.sim.point["dr.integral.drmw", ] <- extract_point(f.sim$dr.integral.drmw, taus[4])
    f.sim.point["dr.integral.drenw", ] <- extract_point(f.sim$dr.integral.drenw, taus[4])
    
    f.sim.point
  }
  f.sim.point$estimate.a0.se <- apply(f.sim.boot[, grep("estimate.a0",names(f.sim.boot))],1, sd)
  f.sim.point$estimate.a1.se <- apply(f.sim.boot[, grep("estimate.a1",names(f.sim.boot))],1, sd)
  f.sim.point$estimate.diff.se <- apply(f.sim.boot[, grep("estimate.diff",names(f.sim.boot))],1, sd)
  
  names(f.sim.point) <- paste0(names(f.sim.point),".v",i)
  f.sim.result <- cbind(f.sim.result, f.sim.point)
  
  # if (i%%10==0){
  print(i)
  # flush.console()
  # }
}

end.time <- Sys.time()
end.time - start.time

save(f.sim.result, file="sim7-baseline.RData") # , f.sim.result.tau6
