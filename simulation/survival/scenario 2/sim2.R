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

library(SuperLearner)
library(survSuperLearner)
library(xgboost)
# library(bartMachine)
library(KernelKnn)
library(earth)
library(fdrtool)
library(simecol)
library(MASS)

registerDoParallel(cores=255)

nsim <- 1000
n <- 1000
# ntrue <- 1e+06
# nboot <- 100

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
a.sl.lib <- c("SL.mean","SL.nnet","SL.kernelKnn","SL.rpartPrune","SL.xgboost","SL.ranger","SL.glm","SL.step","SL.gam","SL.glmnet","SL.earth") # 
event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.loglogreg","survSL.pchreg","survSL.rfsrc","survSL.coxph","survSL.expreg","survSL.weibreg") # ; "survSL.gam","survSL.pchSL"
confounders <- c("x1","x2","x3","x4","x5","x6")

# doubly robust rmst by integrating doubly robust survival with censoring martingale
dr_rmst_integral_mc_ht <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  cens.martingale.integral.a0 <- t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*S.a0*(cens.martingale.integral.a0-1)
  ueif.a0 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a0+term2.a0+term3.a0), 1, cumsum))
  
  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  ueif.a1 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a1+term2.a1+term3.a1), 1, cumsum))
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

dr_rmst_integral_mc_ha <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  cens.martingale.integral.a0 <- t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*S.a0*(cens.martingale.integral.a0-1)
  ueif.a0 <- 1/mean(bw*(a==0))*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a0+term3.a0), 1, cumsum))+
    1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*term2.a0, 1, cumsum))
  
  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  ueif.a1 <- 1/mean(bw*(a==1))*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a1+term3.a1), 1, cumsum))+
    1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*term2.a1, 1, cumsum))
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# doubly robust rmst by integrating doubly robust survival with event martingale
dr_rmst_integral_mt_ht <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  
  event.martingale.integral.a0 <- t(apply((dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term2.a0 <- matrix(ifelse(a, 0, bw),ncol=ns,nrow=n,byrow=FALSE)*S.a0*event.martingale.integral.a0
  ueif.a0 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a0-term2.a0), 1, cumsum))
  
  event.martingale.integral.a1 <- t(apply((dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term2.a1 <- matrix(ifelse(a, bw, 0),ncol=ns,nrow=n,byrow=FALSE)*S.a1*event.martingale.integral.a1
  ueif.a1 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a1-term2.a1), 1, cumsum))
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

dr_rmst_integral_mt_ha <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  
  event.martingale.integral.a0 <- t(apply((dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term2.a0 <- matrix(ifelse(a, 0, bw),ncol=ns,nrow=n,byrow=FALSE)*S.a0*event.martingale.integral.a0
  ueif.a0 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a0), 1, cumsum))-1/mean(bw*(a==0))*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term2.a0), 1, cumsum))
  
  event.martingale.integral.a1 <- t(apply((dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term2.a1 <- matrix(ifelse(a, bw, 0),ncol=ns,nrow=n,byrow=FALSE)*S.a1*event.martingale.integral.a1
  ueif.a1 <- 1/mean(tilt)*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term1.a1), 1, cumsum))-1/mean(bw*(a==1))*t(apply(matrix(ds,ncol=ns,nrow=n,byrow=TRUE)*(term2.a1), 1, cumsum))
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# doubly robust rmst direct with censoring martingale
dr_rmst_direct_mc_ht <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0+term4.a0-term5.a0)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1+term4.a1-term5.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

dr_rmst_direct_mc_ha <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))
  
  Delta.tau.over.G.a0 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a0+matrix(event/rowSums(G.a0*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Delta.tau.over.G.a1 <- do.call(cbind, lapply(1:ns, function(u) ifelse(time > s[u], 1, 0)))/G.a1+matrix(event/rowSums(G.a1*do.call(cbind, lapply(1:ns, function(u) ifelse(time == s[u], 1, 0)))),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u) ifelse(time <= s[u], 1, 0)))
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))
  
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a0*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a0)/G.a0, 1, cumsum))
  term4.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*(t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))-1)
  term5.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*term2.a0+1/mean(bw*(a==0))*(term1.a0+term3.a0+term4.a0-term5.a0)
  
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*Delta.tau.over.G.a1*pmin(matrix(time,ncol=ns,nrow=n,byrow=FALSE), matrix(s,ncol=ns,nrow=n,byrow=TRUE))
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(matrix(s,ncol=ns,nrow=n,byrow=TRUE)*(dNct - Yt*G.dHazard.a1)/G.a1, 1, cumsum))
  term4.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*(t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))-1)
  term5.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*term2.a1+1/mean(bw*(a==1))*(term1.a1+term3.a1+term4.a1-term5.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# doubly robust rmst direct with event martingale
dr_rmst_direct_mt_ht <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*t(apply((dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*(term1.a0-term2.a0+term3.a0)
  
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*t(apply((dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*(term1.a1-term2.a1+term3.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

dr_rmst_direct_mt_ha <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))
  
  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  RMST.a0 <- t(apply(S.a0*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  RMST.a1 <- t(apply(S.a1*matrix(ds,ncol=ns,nrow=n,byrow=TRUE), 1, cumsum))
  S.dHazard.a0 <- t(apply(cbind(0, -log(S.a0)), 1, diff))
  S.dHazard.a1 <- t(apply(cbind(0, -log(S.a1)), 1, diff))
  
  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNt <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 1)))
  
  term1.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a0
  term2.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*RMST.a0*t(apply((dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a0*(dNt - Yt*S.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  ueif.a0 <- 1/mean(tilt)*term1.a0+1/mean(bw*(a==0))*(-term2.a0+term3.a0)
  
  term1.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*RMST.a1
  term2.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*RMST.a1*t(apply((dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*t(apply(RMST.a1*(dNt - Yt*S.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  ueif.a1 <- 1/mean(tilt)*term1.a1+1/mean(bw*(a==1))*(-term2.a1+term3.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

# lognormal <- function(data.df, betas, sigma, n)
# {return(exp(sigma*qnorm(runif(n))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

# loglogistic <- function(data.df, betas, gamma, n)
# {return(exp(gamma*(log(runif(n))-log(1-runif(n)))+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

aft <- function(data.df, betas, sigma, n)
{return(exp(sigma*rnorm(n, mean=0, sd=1)+as.numeric((as.matrix(data.df) %*% as.matrix(betas)))))}

## proportional hazards
dgp4_survival_ph <- function(n, admin.cens, time.grid){
  corr.mat <- matrix(0.5, nrow=3, ncol=3)
  diag(corr.mat) <- 1
  x <- mvrnorm(n=n, mu=rep(0,3), Sigma=corr.mat)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- rbinom(n,1,0.5)
  x5 <- rbinom(n,1,0.5)
  x6 <- rbinom(n,1,0.5)
  ps <- 1/(1+exp(-(-1+x1+1.5*x2+1.5*x3-x4-1.5*x5-x6)))
  a <- rbinom(n, 1, ps) # -1.5+0.9*x1+1.2*x2+1.2*x3+1.2*x4
  data.df <- data.frame(a, x1, x2, x3, x4, x5, x6)
  
  data.df$Ta0 <- simsurv(dist = "exponential", lambdas = 0.12, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.1, x1=0.1, x2=-0.2, x3=0.2, x4=0.1, x5=0.8, x6=-0.2))$eventtime
  data.df$Ta1 <- simsurv(dist = "exponential", lambdas = 0.15, x = data.frame(x0=1, x1, x2, x3, x4, x5, x6), betas = c(x0=0.17, x1=0.2, x2=-0.1, x3=0.4, x4=0.2, x5=0.3, x6=0.4))$eventtime
  data.df$Ta <- (1-a)*data.df$Ta0 + a*data.df$Ta1
  
  data.df$C[data.df$a==0] <- simsurv(dist = "exponential", lambdas = 0.06, x = cbind(x0=1, data.df[data.df$a==0, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0.1, x1=0.4, x2=-0.7, x3=-0.4, x4=-0.5, x5=0.8, x6=-0.6))$eventtime
  data.df$C[data.df$a==1] <- simsurv(dist = "exponential", lambdas = 0.08, x = cbind(x0=1, data.df[data.df$a==1, c("x1", "x2", "x3", "x4", "x5", "x6")]), betas = c(x0=0, x1=0.5, x2=-0.6, x3=0.2, x4=0.6, x5=0.9, x6=-0.5))$eventtime
  data.df$C <- pmin(data.df$C, admin.cens, na.rm=TRUE)
  data.df$end.time <- ceiling(pmin(data.df$C, data.df$Ta)/time.grid)*time.grid # , data.df$Tj2, ceil(pmin(data.df$C,data.df$Tj1,data.df$Tj2)*365.25) # round(pmin(data.df$C,data.df$Tj1,data.df$Tj2),digits=2)
  data.df$event <- with(data.df,ifelse(Ta<=C, 1, 0))
  data.df <- data.df[order(data.df$end.time,-data.df$event),]
  data.df$id <- 1:n
  return(data.df)
}

extract.avarage <- function(train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time, tau){
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])
  
  estimate.a0.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a0))))], y=colMeans(ueif.a0)[1:max(which(!is.na(colMeans(ueif.a0))))], type="lcm")[c("x.knots","y.knots")]
  names(estimate.a0.lcm) <- c("t", "estimate.a0")
  ueif.df <- merge(data.frame(t=sort(unique(c(train.time, estimation.time)))), estimate.a0.lcm, all=TRUE, by="t")
  
  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))))[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])
  estimate.a1.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a1))))], y=colMeans(ueif.a1)[1:max(which(!is.na(colMeans(ueif.a1))))], type="lcm")[c("x.knots","y.knots")]
  names(estimate.a1.lcm) <- c("t", "estimate.a1")
  ueif.df <- merge(ueif.df, estimate.a1.lcm, all=TRUE, by="t")
  
  ueif.df <- approxTime(ueif.df, xout=ueif.df$t)
  ueif.df$estimate.diff <- ueif.df$estimate.a1-ueif.df$estimate.a0
  
  ## standard error
  ueif.df$estimate.a0.se <- cummax(sqrt(colMeans((ueif.a0-matrix(ueif.df$estimate.a0, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  ueif.df$estimate.a1.se <- cummax(sqrt(colMeans((ueif.a1-matrix(ueif.df$estimate.a1, nrow=nrow(ueif.a1), ncol=ncol(ueif.a1), byrow=TRUE))^2))/sqrt(nrow(ueif.a1)))
  ueif.df$estimate.diff.se <- cummax(sqrt(colMeans(((ueif.a1-ueif.a0)-matrix(ueif.df$estimate.diff, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(nrow(ueif.a0)))
  
  tau.ind <- max(which((tau-ueif.df$t)>0))
  ueif.df <- ueif.df[tau.ind,]+(tau-ueif.df$t[tau.ind])*(ueif.df[tau.ind+1,]-ueif.df[tau.ind,])/(ueif.df$t[tau.ind+1]-ueif.df$t[tau.ind])
  return(ueif.df[,-1])
}

start.time <- Sys.time()
f.sim.result <- foreach(i=1:nsim, .combine=function(x,y) mapply(cbind, x, y, SIMPLIFY=FALSE)) %dopar% {
  sim.df <- dgp4_survival_ph(2*n, admin.cens, time.grid)
  # cens.prob[i,"tau2"] <- sum(sim.df$event==0 & sim.df$end.time < taus[1])/(2*n)
  # cens.prob[i,"tau4"] <- sum(sim.df$event==0 & sim.df$end.time < taus[2])/(2*n)
  # cens.prob[i,"tau6"] <- sum(sim.df$event==0 & sim.df$end.time < taus[3])/(2*n)
  
  train.ind <- sample(1:(2*n), n, replace=FALSE) # sample training set
  train.df <- sim.df[train.ind,]
  train.df <- train.df[order(train.df$id),]
  estimation.df <- sim.df[setdiff(1:(2*n), train.ind),]
  estimation.df <- estimation.df[order(estimation.df$id),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]
  
  while(max(train.df$end.time[train.df$a==0 & train.df$event==0])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==0])<taus[4] |
        max(train.df$end.time[train.df$a==0 & train.df$event==1])<taus[4] | max(train.df$end.time[train.df$a==1 & train.df$event==1])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==0])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==0])<taus[4] |
        max(estimation.df$end.time[estimation.df$a==0 & estimation.df$event==1])<taus[4] | max(estimation.df$end.time[estimation.df$a==1 & estimation.df$event==1])<taus[4]) {
    sim.df <- dgp4_survival_ph(2*n, admin.cens, time.grid)
    train.ind <- sample(1:(2*n), n, replace=FALSE) # sample training set
    train.df <- sim.df[train.ind,]
    train.df <- train.df[order(train.df$id),]
    estimation.df <- sim.df[setdiff(1:(2*n), train.ind),]
    estimation.df <- estimation.df[order(estimation.df$id),]}
  
  ## training set
  train.df$naive <- 1
  train.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = estimation.df$a, X = estimation.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                             cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=train.df[, confounders])$pred) # train.df$ps <- predict(glm(a~x1+x2+x3+x4,data=train.df,family=binomial()),type="response")
  train.df$iptw <- with(train.df,1/(ps*a+(1-ps)*(1-a)))
  train.df$iptw.trunc <- with(train.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
  train.df$ow <- with(train.df,ifelse(a,1-ps,ps))
  train.df$mw <- with(train.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
  train.df$enw <- with(train.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
  train.df$ow.tilt <- with(train.df, ps*(1-ps))
  train.df$mw.tilt <- with(train.df, pmin(ps, 1-ps))
  train.df$enw.tilt <- with(train.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
  
  ## estimation set
  estimation.df$naive <- 1
  estimation.df$ps <- Winsorize(minval = 1e-3, maxval = 1-1e-3, predict.SuperLearner(SuperLearner(Y = train.df$a, X = train.df[, confounders], family = binomial(), SL.library = a.sl.lib,
                                                                                                  cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL)), newdata=estimation.df[, confounders])$pred) # estimation.df$ps <- predict(glm(a~x1+x2+x3+x4,data=estimation.df,family=binomial()),type="response")
  estimation.df$iptw <- with(estimation.df,1/(ps*a+(1-ps)*(1-a)))
  estimation.df$iptw.trunc <- with(estimation.df,ifelse(iptw>quantile(iptw,0.99),quantile(iptw,0.99),iptw))
  estimation.df$ow <- with(estimation.df,ifelse(a,1-ps,ps))
  estimation.df$mw <- with(estimation.df,pmin(ps,1-ps)*ifelse(a,1/ps,1/(1-ps)))
  estimation.df$enw <- with(estimation.df,-(ps*log(ps)+(1-ps)*log(1-ps))*ifelse(a,1/ps,1/(1-ps)))
  estimation.df$ow.tilt <- with(estimation.df, ps*(1-ps))
  estimation.df$mw.tilt <- with(estimation.df, pmin(ps, 1-ps))
  estimation.df$enw.tilt <- with(estimation.df,-(ps*log(ps)+(1-ps)*log(1-ps)))
  
  ## train
  train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a0 <- train.all.survival.super.a0$cens.SL.predict
  train.S.a0 <- train.all.survival.super.a0$event.SL.predict
  
  train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                  newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  train.G.a1 <- train.all.survival.super.a1$cens.SL.predict
  train.S.a1 <- train.all.survival.super.a1$event.SL.predict
  
  ## estimation
  estimation.all.survival.super.a0 <- survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a0 <- estimation.all.survival.super.a0$cens.SL.predict
  estimation.S.a0 <- estimation.all.survival.super.a0$event.SL.predict
  
  estimation.all.survival.super.a1 <- survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                       newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
  estimation.G.a1 <- estimation.all.survival.super.a1$cens.SL.predict
  estimation.S.a1 <- estimation.all.survival.super.a1$event.SL.predict
  
  ## train
  f.train <- list("dr.integral.mc.ht.driptw.ueif"=NA, "dr.integral.mc.ht.driptw.trunc.ueif"=NA, "dr.integral.mc.ht.drow.ueif"=NA, "dr.integral.mc.ht.drmw.ueif"=NA, "dr.integral.mc.ht.drenw.ueif"=NA,
                  "dr.integral.mc.ha.driptw.ueif"=NA, "dr.integral.mc.ha.driptw.trunc.ueif"=NA, "dr.integral.mc.ha.drow.ueif"=NA, "dr.integral.mc.ha.drmw.ueif"=NA, "dr.integral.mc.ha.drenw.ueif"=NA,
                  "dr.integral.mt.ht.driptw.ueif"=NA, "dr.integral.mt.ht.driptw.trunc.ueif"=NA, "dr.integral.mt.ht.drow.ueif"=NA, "dr.integral.mt.ht.drmw.ueif"=NA, "dr.integral.mt.ht.drenw.ueif"=NA,
                  "dr.integral.mt.ha.driptw.ueif"=NA, "dr.integral.mt.ha.driptw.trunc.ueif"=NA, "dr.integral.mt.ha.drow.ueif"=NA, "dr.integral.mt.ha.drmw.ueif"=NA, "dr.integral.mt.ha.drenw.ueif"=NA,
                  "dr.direct.mc.ht.driptw.ueif"=NA, "dr.direct.mc.ht.driptw.trunc.ueif"=NA, "dr.direct.mc.ht.drow.ueif"=NA, "dr.direct.mc.ht.drmw.ueif"=NA, "dr.direct.mc.ht.drenw.ueif"=NA,
                  "dr.direct.mc.ha.driptw.ueif"=NA, "dr.direct.mc.ha.driptw.trunc.ueif"=NA, "dr.direct.mc.ha.drow.ueif"=NA, "dr.direct.mc.ha.drmw.ueif"=NA, "dr.direct.mc.ha.drenw.ueif"=NA,
                  "dr.direct.mt.ht.driptw.ueif"=NA, "dr.direct.mt.ht.driptw.trunc.ueif"=NA, "dr.direct.mt.ht.drow.ueif"=NA, "dr.direct.mt.ht.drmw.ueif"=NA, "dr.direct.mt.ht.drenw.ueif"=NA,
                  "dr.direct.mt.ha.driptw.ueif"=NA, "dr.direct.mt.ha.driptw.trunc.ueif"=NA, "dr.direct.mt.ha.drow.ueif"=NA, "dr.direct.mt.ha.drmw.ueif"=NA, "dr.direct.mt.ha.drenw.ueif"=NA)
  #
  f.train$dr.integral.mc.ht.driptw.ueif <- dr_rmst_integral_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ht.driptw.trunc.ueif <- dr_rmst_integral_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ht.drow.ueif <- dr_rmst_integral_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ht.drmw.ueif <- dr_rmst_integral_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ht.drenw.ueif <- dr_rmst_integral_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)

  f.train$dr.integral.mc.ha.driptw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ha.driptw.trunc.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ha.drow.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ha.drmw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mc.ha.drenw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)

  f.train$dr.integral.mt.ht.driptw.ueif <- dr_rmst_integral_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ht.driptw.trunc.ueif <- dr_rmst_integral_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ht.drow.ueif <- dr_rmst_integral_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ht.drmw.ueif <- dr_rmst_integral_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ht.drenw.ueif <- dr_rmst_integral_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)

  f.train$dr.integral.mt.ha.driptw.ueif <- dr_rmst_integral_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ha.driptw.trunc.ueif <- dr_rmst_integral_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ha.drow.ueif <- dr_rmst_integral_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ha.drmw.ueif <- dr_rmst_integral_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.integral.mt.ha.drenw.ueif <- dr_rmst_integral_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
 
  f.train$dr.direct.mc.ht.driptw.ueif <- dr_rmst_direct_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ht.driptw.trunc.ueif <- dr_rmst_direct_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ht.drow.ueif <- dr_rmst_direct_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ht.drmw.ueif <- dr_rmst_direct_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ht.drenw.ueif <- dr_rmst_direct_mc_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  
  f.train$dr.direct.mc.ha.driptw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ha.driptw.trunc.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ha.drow.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ha.drmw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mc.ha.drenw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  
  f.train$dr.direct.mt.ht.driptw.ueif <- dr_rmst_direct_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ht.driptw.trunc.ueif <- dr_rmst_direct_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ht.drow.ueif <- dr_rmst_direct_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ht.drmw.ueif <- dr_rmst_direct_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ht.drenw.ueif <- dr_rmst_direct_mt_ht(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  
  f.train$dr.direct.mt.ha.driptw.ueif <- dr_rmst_direct_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ha.driptw.trunc.ueif <- dr_rmst_direct_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ha.drow.ueif <- dr_rmst_direct_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ha.drmw.ueif <- dr_rmst_direct_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  f.train$dr.direct.mt.ha.drenw.ueif <- dr_rmst_direct_mt_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
  
  ## estimation
  f.estimation <- list("dr.integral.mc.ht.driptw.ueif"=NA, "dr.integral.mc.ht.driptw.trunc.ueif"=NA, "dr.integral.mc.ht.drow.ueif"=NA, "dr.integral.mc.ht.drmw.ueif"=NA, "dr.integral.mc.ht.drenw.ueif"=NA,
                       "dr.integral.mc.ha.driptw.ueif"=NA, "dr.integral.mc.ha.driptw.trunc.ueif"=NA, "dr.integral.mc.ha.drow.ueif"=NA, "dr.integral.mc.ha.drmw.ueif"=NA, "dr.integral.mc.ha.drenw.ueif"=NA,
                       "dr.integral.mt.ht.driptw.ueif"=NA, "dr.integral.mt.ht.driptw.trunc.ueif"=NA, "dr.integral.mt.ht.drow.ueif"=NA, "dr.integral.mt.ht.drmw.ueif"=NA, "dr.integral.mt.ht.drenw.ueif"=NA,
                       "dr.integral.mt.ha.driptw.ueif"=NA, "dr.integral.mt.ha.driptw.trunc.ueif"=NA, "dr.integral.mt.ha.drow.ueif"=NA, "dr.integral.mt.ha.drmw.ueif"=NA, "dr.integral.mt.ha.drenw.ueif"=NA,
                       "dr.direct.mc.ht.driptw.ueif"=NA, "dr.direct.mc.ht.driptw.trunc.ueif"=NA, "dr.direct.mc.ht.drow.ueif"=NA, "dr.direct.mc.ht.drmw.ueif"=NA, "dr.direct.mc.ht.drenw.ueif"=NA,
                       "dr.direct.mc.ha.driptw.ueif"=NA, "dr.direct.mc.ha.driptw.trunc.ueif"=NA, "dr.direct.mc.ha.drow.ueif"=NA, "dr.direct.mc.ha.drmw.ueif"=NA, "dr.direct.mc.ha.drenw.ueif"=NA,
                       "dr.direct.mt.ht.driptw.ueif"=NA, "dr.direct.mt.ht.driptw.trunc.ueif"=NA, "dr.direct.mt.ht.drow.ueif"=NA, "dr.direct.mt.ht.drmw.ueif"=NA, "dr.direct.mt.ht.drenw.ueif"=NA,
                       "dr.direct.mt.ha.driptw.ueif"=NA, "dr.direct.mt.ha.driptw.trunc.ueif"=NA, "dr.direct.mt.ha.drow.ueif"=NA, "dr.direct.mt.ha.drmw.ueif"=NA, "dr.direct.mt.ha.drenw.ueif"=NA)
  #
  f.estimation$dr.integral.mc.ht.driptw.ueif <- dr_rmst_integral_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ht.driptw.trunc.ueif <- dr_rmst_integral_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ht.drow.ueif <- dr_rmst_integral_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ht.drmw.ueif <- dr_rmst_integral_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ht.drenw.ueif <- dr_rmst_integral_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.integral.mc.ha.driptw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ha.driptw.trunc.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ha.drow.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ha.drmw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mc.ha.drenw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.integral.mt.ht.driptw.ueif <- dr_rmst_integral_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ht.driptw.trunc.ueif <- dr_rmst_integral_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ht.drow.ueif <- dr_rmst_integral_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ht.drmw.ueif <- dr_rmst_integral_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ht.drenw.ueif <- dr_rmst_integral_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.integral.mt.ha.driptw.ueif <- dr_rmst_integral_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ha.driptw.trunc.ueif <- dr_rmst_integral_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ha.drow.ueif <- dr_rmst_integral_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ha.drmw.ueif <- dr_rmst_integral_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.integral.mt.ha.drenw.ueif <- dr_rmst_integral_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.direct.mc.ht.driptw.ueif <- dr_rmst_direct_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ht.driptw.trunc.ueif <- dr_rmst_direct_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ht.drow.ueif <- dr_rmst_direct_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ht.drmw.ueif <- dr_rmst_direct_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ht.drenw.ueif <- dr_rmst_direct_mc_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.direct.mc.ha.driptw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ha.driptw.trunc.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ha.drow.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ha.drmw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mc.ha.drenw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.direct.mt.ht.driptw.ueif <- dr_rmst_direct_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ht.driptw.trunc.ueif <- dr_rmst_direct_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ht.drow.ueif <- dr_rmst_direct_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ht.drmw.ueif <- dr_rmst_direct_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ht.drenw.ueif <- dr_rmst_direct_mt_ht(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.estimation$dr.direct.mt.ha.driptw.ueif <- dr_rmst_direct_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ha.driptw.trunc.ueif <- dr_rmst_direct_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ha.drow.ueif <- dr_rmst_direct_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ha.drmw.ueif <- dr_rmst_direct_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  f.estimation$dr.direct.mt.ha.drenw.ueif <- dr_rmst_direct_mt_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
  
  f.sim <- rep(list(data.frame(estimate.a0=rep(NA, 40), estimate.a1=rep(NA, 40), estimate.diff=rep(NA, 40), estimate.a0.se=rep(NA, 40), estimate.a1.se=rep(NA, 40), estimate.diff.se=rep(NA, 40))), length(taus))
  names(f.sim) <- paste0("tau", taus)
  f.sim <- lapply(f.sim, "colnames<-", paste0(colnames(f.sim[[1]]), ".v", i))
  
  f.sim <- lapply(f.sim, "rownames<-", c("dr.integral.mc.ht.driptw", "dr.integral.mc.ht.driptw.trunc", "dr.integral.mc.ht.drow", "dr.integral.mc.ht.drmw", "dr.integral.mc.ht.drenw",
                                         "dr.integral.mc.ha.driptw", "dr.integral.mc.ha.driptw.trunc", "dr.integral.mc.ha.drow", "dr.integral.mc.ha.drmw", "dr.integral.mc.ha.drenw",
                                         "dr.integral.mt.ht.driptw", "dr.integral.mt.ht.driptw.trunc", "dr.integral.mt.ht.drow", "dr.integral.mt.ht.drmw", "dr.integral.mt.ht.drenw",
                                         "dr.integral.mt.ha.driptw", "dr.integral.mt.ha.driptw.trunc", "dr.integral.mt.ha.drow", "dr.integral.mt.ha.drmw", "dr.integral.mt.ha.drenw",
                                         "dr.direct.mc.ht.driptw", "dr.direct.mc.ht.driptw.trunc", "dr.direct.mc.ht.drow", "dr.direct.mc.ht.drmw", "dr.direct.mc.ht.drenw",
                                         "dr.direct.mc.ha.driptw", "dr.direct.mc.ha.driptw.trunc", "dr.direct.mc.ha.drow", "dr.direct.mc.ha.drmw", "dr.direct.mc.ha.drenw",
                                         "dr.direct.mt.ht.driptw", "dr.direct.mt.ht.driptw.trunc", "dr.direct.mt.ht.drow", "dr.direct.mt.ht.drmw", "dr.direct.mt.ht.drenw",
                                         "dr.direct.mt.ha.driptw", "dr.direct.mt.ha.driptw.trunc", "dr.direct.mt.ha.drow", "dr.direct.mt.ha.drmw", "dr.direct.mt.ha.drenw"))
  for (j in 1:length(taus))
  {
    #
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ht.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ht.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ht.driptw.ueif$ueif.a0,
                                                                                    estimation.ueif.a1=f.estimation$dr.integral.mc.ht.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ht.driptw.ueif$ueif.a0,
                                                                                    train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ht.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                          estimation.ueif.a1=f.estimation$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                          train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ht.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ht.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ht.drow.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mc.ht.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ht.drow.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ht.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ht.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ht.drmw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mc.ht.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ht.drmw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ht.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ht.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ht.drenw.ueif$ueif.a0,
                                                                                   estimation.ueif.a1=f.estimation$dr.integral.mc.ht.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ht.drenw.ueif$ueif.a0,
                                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ha.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ha.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ha.driptw.ueif$ueif.a0,
                                                                                    estimation.ueif.a1=f.estimation$dr.integral.mc.ha.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ha.driptw.ueif$ueif.a0,
                                                                                    train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ha.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ha.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ha.driptw.trunc.ueif$ueif.a0,
                                                                                          estimation.ueif.a1=f.estimation$dr.integral.mc.ha.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ha.driptw.trunc.ueif$ueif.a0,
                                                                                          train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ha.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ha.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ha.drow.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mc.ha.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ha.drow.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ha.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ha.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ha.drmw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mc.ha.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ha.drmw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mc.ha.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mc.ha.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mc.ha.drenw.ueif$ueif.a0,
                                                                                   estimation.ueif.a1=f.estimation$dr.integral.mc.ha.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mc.ha.drenw.ueif$ueif.a0,
                                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ht.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ht.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ht.driptw.ueif$ueif.a0,
                                                                                    estimation.ueif.a1=f.estimation$dr.integral.mt.ht.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ht.driptw.ueif$ueif.a0,
                                                                                    train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ht.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ht.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ht.driptw.trunc.ueif$ueif.a0,
                                                                                          estimation.ueif.a1=f.estimation$dr.integral.mt.ht.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ht.driptw.trunc.ueif$ueif.a0,
                                                                                          train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ht.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ht.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ht.drow.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mt.ht.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ht.drow.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ht.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ht.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ht.drmw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mt.ht.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ht.drmw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ht.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ht.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ht.drenw.ueif$ueif.a0,
                                                                                   estimation.ueif.a1=f.estimation$dr.integral.mt.ht.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ht.drenw.ueif$ueif.a0,
                                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ha.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ha.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ha.driptw.ueif$ueif.a0,
                                                                                    estimation.ueif.a1=f.estimation$dr.integral.mt.ha.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ha.driptw.ueif$ueif.a0,
                                                                                    train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ha.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ha.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ha.driptw.trunc.ueif$ueif.a0,
                                                                                          estimation.ueif.a1=f.estimation$dr.integral.mt.ha.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ha.driptw.trunc.ueif$ueif.a0,
                                                                                          train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ha.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ha.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ha.drow.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mt.ha.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ha.drow.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ha.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ha.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ha.drmw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.integral.mt.ha.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ha.drmw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.integral.mt.ha.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.integral.mt.ha.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.integral.mt.ha.drenw.ueif$ueif.a0,
                                                                                   estimation.ueif.a1=f.estimation$dr.integral.mt.ha.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.integral.mt.ha.drenw.ueif$ueif.a0,
                                                                                   train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ht.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ht.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ht.driptw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.direct.mc.ht.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ht.driptw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ht.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ht.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                        estimation.ueif.a1=f.estimation$dr.direct.mc.ht.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ht.driptw.trunc.ueif$ueif.a0,
                                                                                        train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ht.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ht.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ht.drow.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mc.ht.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ht.drow.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ht.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ht.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ht.drmw.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mc.ht.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ht.drmw.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ht.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ht.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ht.drenw.ueif$ueif.a0,
                                                                                 estimation.ueif.a1=f.estimation$dr.direct.mc.ht.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ht.drenw.ueif$ueif.a0,
                                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ha.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ha.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.driptw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.direct.mc.ha.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.driptw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ha.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ha.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.driptw.trunc.ueif$ueif.a0,
                                                                                        estimation.ueif.a1=f.estimation$dr.direct.mc.ha.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.driptw.trunc.ueif$ueif.a0,
                                                                                        train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ha.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ha.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.drow.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mc.ha.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.drow.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ha.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ha.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.drmw.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mc.ha.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.drmw.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mc.ha.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mc.ha.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.drenw.ueif$ueif.a0,
                                                                                 estimation.ueif.a1=f.estimation$dr.direct.mc.ha.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.drenw.ueif$ueif.a0,
                                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ht.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ht.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ht.driptw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.direct.mt.ht.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ht.driptw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ht.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ht.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ht.driptw.trunc.ueif$ueif.a0,
                                                                                        estimation.ueif.a1=f.estimation$dr.direct.mt.ht.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ht.driptw.trunc.ueif$ueif.a0,
                                                                                        train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ht.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ht.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ht.drow.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mt.ht.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ht.drow.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ht.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ht.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ht.drmw.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mt.ht.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ht.drmw.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ht.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ht.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ht.drenw.ueif$ueif.a0,
                                                                                 estimation.ueif.a1=f.estimation$dr.direct.mt.ht.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ht.drenw.ueif$ueif.a0,
                                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ha.driptw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ha.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ha.driptw.ueif$ueif.a0,
                                                                                  estimation.ueif.a1=f.estimation$dr.direct.mt.ha.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ha.driptw.ueif$ueif.a0,
                                                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ha.driptw.trunc",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ha.driptw.trunc.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ha.driptw.trunc.ueif$ueif.a0,
                                                                                        estimation.ueif.a1=f.estimation$dr.direct.mt.ha.driptw.trunc.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ha.driptw.trunc.ueif$ueif.a0,
                                                                                        train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ha.drow",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ha.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ha.drow.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mt.ha.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ha.drow.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ha.drmw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ha.drmw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ha.drmw.ueif$ueif.a0,
                                                                                estimation.ueif.a1=f.estimation$dr.direct.mt.ha.drmw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ha.drmw.ueif$ueif.a0,
                                                                                train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    f.sim[[paste0("tau", taus[j])]]["dr.direct.mt.ha.drenw",] <- extract.avarage(train.ueif.a1=f.train$dr.direct.mt.ha.drenw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mt.ha.drenw.ueif$ueif.a0,
                                                                                 estimation.ueif.a1=f.estimation$dr.direct.mt.ha.drenw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mt.ha.drenw.ueif$ueif.a0,
                                                                                 train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)), tau=taus[j])
    }
  # if (i%%10==0){
  print(i)
  # flush.console()
  # }
  f.sim
}

end.time <- Sys.time()
end.time - start.time

save(cens.prob, f.sim.result, file="sim2.RData") # , f.sim.result.tau6
