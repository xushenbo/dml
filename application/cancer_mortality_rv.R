## RAM: 8 GB; core: 104; 10 min
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
# library(KernelKnn)
library(earth)
library(fdrtool)
library(simecol)
library(MASS)
library(data.table)
library(ggplot2)
# library(mvnfast)
# library(lqmm)
library(Rfast)

registerDoParallel(cores=100)

npath <- 1000
a.sl.lib <- c("SL.mean","SL.glm","SL.nnet","SL.rpartPrune","SL.xgboost","SL.ranger","SL.step","SL.glmnet","SL.earth") # a.sl.lib <- c("SL.mean","SL.glm") # ,"SL.kernelKnn","SL.gam"
event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.loglogreg","survSL.rfsrc","survSL.coxph","survSL.expreg","survSL.weibreg") #  event.sl.lib <- cens.sl.lib <- c("survSL.km","survSL.coxph") "survSL.pchreg": too slow; # ; "survSL.gam","survSL.pchSL"

####################################################################################
load("antidiabetes_aje.RData")
data.df <- data.df %>% data.table()
data.df[prescription_date<af_eventdate,"af"] <- 0
data.df[prescription_date<hp_combined_eventdate,"hypertension_combined"] <- 0
data.df[prescription_date<chd_combined_eventdate,"chd_combined"] <- 0
data.df[prescription_date<pvd_eventdate,"pvd"] <- 0
data.df[prescription_date<hf_eventdate,"hf"] <- 0
data.df[prescription_date<copd_eventdate,"copd"] <- 0
data.df[prescription_date<stroke_combined_eventdate,"stroke_combined"] <- 0
data.df$diabetes_after_prescription <- with(data.df, pmax(as.numeric(difftime(prescription_date, diabetes_combined_eventdate, units="days")), 0, na.rm=TRUE))
data.df$tod_n[which(with(data.df, tod_n<cancer_eventdate|tod_n<deathdate_n|tod_n<prescription_date))] <- NA # with(data.df, pmin(cancer_eventdate, deathdate_n, na.rm=TRUE)[which(data.df$tod_n<data.df$cancer_eventdate|data.df$tod_n<data.df$deathdate_n)])

data.df %<>%
  # rowwise() %>%
  filter(prescription_date > cancer_eventdate) %>%
  mutate(event1 = death) %>%
  mutate(end.time1 = as.numeric(difftime(pmin(deathdate_n, tod_n, lcd_n, na.rm=TRUE), prescription_date, units="days"))) %>% # /365.242
  data.frame()

confounders <- c("gender","prescription_age","prescription_year","region","imd_imputed","bmi","hba1c_converge_1year_latest",
                 "smoke","hf","chd_combined","af","stroke_combined","hypertension_combined","pvd","copd","diabetes_after_prescription")
data.df <- data.df[complete.cases(data.df[, confounders]), ]
# data.df <- data.df[1:2000,]
time.grid <- 1
admin.cens <- max(data.df$end.time1)
n <- nrow(data.df)
data.df$a <- data.df$ref_sulf
data.df$end.time <- data.df$end.time1
data.df$event <- data.df$event1
data.df$id <- data.df$patid

####################################################################################
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

survival_inference <- function(npath, train.ueif.a1, train.ueif.a0, estimation.ueif.a1, estimation.ueif.a0, train.time, estimation.time){
  n <- nrow(train.ueif.a1)+nrow(estimation.ueif.a1)
  ## point estimate
  # a=0
  ueif.a0 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a0)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"), all=TRUE, by="t")
  ueif.a0 <- approxTime(ueif.a0, xout=ueif.a0$t, rule=2)
  if(sum(is.na(ueif.a0))>0)
  {ueif.a0 <- approxTime(rbind(0, ueif.a0), xout=c(0, sort(unique(c(train.time, estimation.time)))), rule=2)[-1,]}
  ueif.a0 <- t(ueif.a0[,-1])
  
  estimate.a0.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a0))))], y=colMeans(ueif.a0)[1:max(which(!is.na(colMeans(ueif.a0))))], type="lcm")[c("x.knots","y.knots")]
  names(estimate.a0.lcm) <- c("t", "estimate.a0")
  ueif.df <- merge(data.frame(t=sort(unique(c(train.time, estimation.time)))), estimate.a0.lcm, all=TRUE, by="t")
  
  # a=1
  ueif.a1 <- merge(merge(data.frame(t=matrix(train.time), t(train.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t"),
                   merge(data.frame(t=matrix(estimation.time), t(estimation.ueif.a1)), data.frame(t=sort(unique(c(train.time, estimation.time)))), all=TRUE, by="t") , all=TRUE, by="t")
  ueif.a1 <- approxTime(ueif.a1, xout=ueif.a1$t, rule=2)
  if(sum(is.na(ueif.a1))>0)
  {ueif.a1 <- approxTime(rbind(0, ueif.a1), xout=c(0, sort(unique(c(train.time, estimation.time)))), rule=2)[-1,]}
  ueif.a1 <- t(ueif.a1[,-1])
  
  estimate.a1.lcm <- gcmlcm(x=sort(unique(c(train.time, estimation.time)))[1:max(which(!is.na(colMeans(ueif.a1))))], y=colMeans(ueif.a1)[1:max(which(!is.na(colMeans(ueif.a1))))], type="lcm")[c("x.knots","y.knots")]
  names(estimate.a1.lcm) <- c("t", "estimate.a1")
  ueif.df <- merge(ueif.df, estimate.a1.lcm, all=TRUE, by="t")
  
  # diff
  ueif.df <- approxTime(ueif.df, xout=ueif.df$t, rule=2)
  ueif.df$estimate.diff <- ueif.df$estimate.a1-ueif.df$estimate.a0
  
  ## standard error
  ueif.df$estimate.a0.se <- cummax(sqrt(colMeans((ueif.a0-matrix(ueif.df$estimate.a0, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(n))
  ueif.df$estimate.a1.se <- cummax(sqrt(colMeans((ueif.a1-matrix(ueif.df$estimate.a1, nrow=nrow(ueif.a1), ncol=ncol(ueif.a1), byrow=TRUE))^2))/sqrt(n))
  ueif.df$estimate.diff.se <- cummax(sqrt(colMeans(((ueif.a1-ueif.a0)-matrix(ueif.df$estimate.diff, nrow=nrow(ueif.a0), ncol=ncol(ueif.a0), byrow=TRUE))^2))/sqrt(n))
  
  ueif.df$estimate.a0.lower.ci <- ueif.df$estimate.a0+qnorm(0.025)*ueif.df$estimate.a0.se
  ueif.df$estimate.a0.upper.ci <- ueif.df$estimate.a0+qnorm(0.975)*ueif.df$estimate.a0.se
  ueif.df$estimate.a1.lower.ci <- ueif.df$estimate.a1+qnorm(0.025)*ueif.df$estimate.a1.se
  ueif.df$estimate.a1.upper.ci <- ueif.df$estimate.a1+qnorm(0.975)*ueif.df$estimate.a1.se
  ueif.df$estimate.diff.lower.ci <- ueif.df$estimate.diff+qnorm(0.025)*ueif.df$estimate.diff.se
  ueif.df$estimate.diff.upper.ci <- ueif.df$estimate.diff+qnorm(0.975)*ueif.df$estimate.diff.se
  
  ## confidence bands
  # a=0
  tl <- ifelse(sum(ueif.df$estimate.a0.se==0), max(which(ueif.df$estimate.a0.se==0))+1, 1)
  ceif.a0 <- (ueif.a0[, tl:nrow(ueif.df)]-matrix(ueif.df$estimate.a0[tl:nrow(ueif.df)], nrow=nrow(ueif.a0), ncol=(ncol(ueif.a0)-tl+1), byrow=TRUE))/matrix(ueif.df$estimate.a0.se[tl:nrow(ueif.df)], nrow=nrow(ueif.a0), ncol=(ncol(ueif.a0)-tl+1), byrow=TRUE)
  c.a0 <- foreach(i=1:npath, .combine=c) %dopar% {
    max(abs(rbind(Rfast::Rnorm(n)/sqrt(n)) %*% ceif.a0))
  } %>% quantile(., probs=0.95)
  ueif.df$estimate.a0.lower.cb <- ueif.df$estimate.a0-ueif.df$estimate.a0.se*c.a0/sqrt(n)
  ueif.df$estimate.a0.upper.cb <- ueif.df$estimate.a0+ueif.df$estimate.a0.se*c.a0/sqrt(n)
  
  # a=1
  tl <- ifelse(sum(ueif.df$estimate.a1.se==0), max(which(ueif.df$estimate.a1.se==0))+1, 1)
  ceif.a1 <- (ueif.a1[, tl:nrow(ueif.df)]-matrix(ueif.df$estimate.a1[tl:nrow(ueif.df)], nrow=nrow(ueif.a1), ncol=(ncol(ueif.a1)-tl+1), byrow=TRUE))/matrix(ueif.df$estimate.a1.se[tl:nrow(ueif.df)], nrow=nrow(ueif.a1), ncol=(ncol(ueif.a1)-tl+1), byrow=TRUE)
  c.a1 <- foreach(i=1:npath, .combine=c) %dopar% {
    max(abs(rbind(Rfast::Rnorm(n)/sqrt(n)) %*% ceif.a1))
  } %>% quantile(., probs=0.95)
  ueif.df$estimate.a1.lower.cb <- ueif.df$estimate.a1-ueif.df$estimate.a1.se*c.a1/sqrt(n)
  ueif.df$estimate.a1.upper.cb <- ueif.df$estimate.a1+ueif.df$estimate.a1.se*c.a1/sqrt(n)
  
  # diff
  tl <- ifelse(sum(ueif.df$estimate.diff.se==0), max(which(ueif.df$estimate.diff.se==0))+1, 1)
  ueif.diff <- ueif.a1-ueif.a0
  ceif.diff <- (ueif.diff[, tl:nrow(ueif.df)]-matrix(ueif.df$estimate.diff[tl:nrow(ueif.df)], nrow=nrow(ueif.diff), ncol=(ncol(ueif.diff)-tl+1), byrow=TRUE))/matrix(ueif.df$estimate.diff.se[tl:nrow(ueif.df)], nrow=nrow(ueif.diff), ncol=(ncol(ueif.diff)-tl+1), byrow=TRUE)
  c.diff <- foreach(i=1:npath, .combine=c) %dopar% {
    max(abs(rbind(Rfast::Rnorm(n)/sqrt(n)) %*% ceif.diff))
  } %>% quantile(., probs=0.95)
  ueif.df$estimate.diff.lower.cb <- ueif.df$estimate.diff-ueif.df$estimate.diff.se*c.diff/sqrt(n)
  ueif.df$estimate.diff.upper.cb <- ueif.df$estimate.diff+ueif.df$estimate.diff.se*c.diff/sqrt(n)
  return(ueif.df)
}

####################################################################################
start.time <- Sys.time()

train.ind <- sample(1:n, ceiling(n/2), replace=FALSE) # sample training set
train.df <- data.df[train.ind,]
train.df <- train.df[order(train.df$patid),]
estimation.df <- data.df[setdiff(1:n, train.ind),]
estimation.df <- estimation.df[order(estimation.df$patid),] # estimation.df <- estimation.df[order(estimation.df$end.time,-estimation.df$event),]

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
# conditional survival, censoring
train.all.survival.super.a0 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==0], event = 1*(estimation.df$event[estimation.df$a==0]!=0), X = estimation.df[estimation.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a0 <- train.all.survival.super.a0$cens.SL.predict
train.S.a0 <- train.all.survival.super.a0$event.SL.predict

train.all.survival.super.a1 <- survSuperLearner(time = estimation.df$end.time[estimation.df$a==1], event = 1*(estimation.df$event[estimation.df$a==1]!=0), X = estimation.df[estimation.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                newX = train.df[, confounders], new.times = sort(unique(train.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
train.G.a1 <- train.all.survival.super.a1$cens.SL.predict
train.S.a1 <- train.all.survival.super.a1$event.SL.predict

## estimation
# conditional survival, censoring
estimation.all.survival.super.a0 <- survSuperLearner(time = train.df$end.time[train.df$a==0], event = 1*(train.df$event[train.df$a==0]!=0), X = train.df[train.df$a==0, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
estimation.G.a0 <- estimation.all.survival.super.a0$cens.SL.predict
estimation.S.a0 <- estimation.all.survival.super.a0$event.SL.predict

estimation.all.survival.super.a1 <- survSuperLearner(time = train.df$end.time[train.df$a==1], event = 1*(train.df$event[train.df$a==1]!=0), X = train.df[train.df$a==1, confounders], cvControl = list(V = 2L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL),
                                                     newX = estimation.df[, confounders], new.times = sort(unique(estimation.df$end.time)), event.SL.library = event.sl.lib, cens.SL.library = cens.sl.lib, verbose = FALSE)
estimation.G.a1 <- estimation.all.survival.super.a1$cens.SL.predict
estimation.S.a1 <- estimation.all.survival.super.a1$event.SL.predict

end.time <- Sys.time()
end.time - start.time

start.time <- Sys.time()
## train
f.train <- list("dr.integral.mc.ha.driptw.ueif"=NA, "dr.integral.mc.ha.driptw.trunc.ueif"=NA, "dr.integral.mc.ha.drow.ueif"=NA, "dr.integral.mc.ha.drmw.ueif"=NA, "dr.integral.mc.ha.drenw.ueif"=NA,
                "dr.direct.mc.ha.driptw.ueif"=NA, "dr.direct.mc.ha.driptw.trunc.ueif"=NA, "dr.direct.mc.ha.drow.ueif"=NA, "dr.direct.mc.ha.drmw.ueif"=NA, "dr.direct.mc.ha.drenw.ueif"=NA)
#
# f.train$dr.integral.mc.ha.driptw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.integral.mc.ha.driptw.trunc.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.integral.mc.ha.drow.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.integral.mc.ha.drmw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.integral.mc.ha.drenw.ueif <- dr_rmst_integral_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)

f.train$dr.direct.mc.ha.driptw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.direct.mc.ha.driptw.trunc.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$iptw.trunc, tilt=train.df$naive, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
f.train$dr.direct.mc.ha.drow.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$ow, tilt=train.df$ow.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.direct.mc.ha.drmw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$mw, tilt=train.df$mw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)
# f.train$dr.direct.mc.ha.drenw.ueif <- dr_rmst_direct_mc_ha(id=train.df$id, a=train.df$a, time=train.df$end.time, event=train.df$event, bw=train.df$enw, tilt=train.df$enw.tilt, G.a0=train.G.a0, G.a1=train.G.a1, S.a0=train.S.a0, S.a1=train.S.a1)

## estimation
f.estimation <- list("dr.integral.mc.ha.driptw.ueif"=NA, "dr.integral.mc.ha.driptw.trunc.ueif"=NA, "dr.integral.mc.ha.drow.ueif"=NA, "dr.integral.mc.ha.drmw.ueif"=NA, "dr.integral.mc.ha.drenw.ueif"=NA,
                     "dr.direct.mc.ha.driptw.ueif"=NA, "dr.direct.mc.ha.driptw.trunc.ueif"=NA, "dr.direct.mc.ha.drow.ueif"=NA, "dr.direct.mc.ha.drmw.ueif"=NA, "dr.direct.mc.ha.drenw.ueif"=NA)
#
# f.estimation$dr.integral.mc.ha.driptw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.integral.mc.ha.driptw.trunc.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.integral.mc.ha.drow.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.integral.mc.ha.drmw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.integral.mc.ha.drenw.ueif <- dr_rmst_integral_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)

f.estimation$dr.direct.mc.ha.driptw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.direct.mc.ha.driptw.trunc.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$iptw.trunc, tilt=estimation.df$naive, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
f.estimation$dr.direct.mc.ha.drow.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$ow, tilt=estimation.df$ow.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.direct.mc.ha.drmw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$mw, tilt=estimation.df$mw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)
# f.estimation$dr.direct.mc.ha.drenw.ueif <- dr_rmst_direct_mc_ha(id=estimation.df$id, a=estimation.df$a, time=estimation.df$end.time, event=estimation.df$event, bw=estimation.df$enw, tilt=estimation.df$enw.tilt, G.a0=estimation.G.a0, G.a1=estimation.G.a1, S.a0=estimation.S.a0, S.a1=estimation.S.a1)

end.time <- Sys.time()
end.time - start.time

f.sim <- list("dr.integral.mc.ha.driptw"=NA, "dr.integral.mc.ha.driptw.trunc"=NA, "dr.integral.mc.ha.drow"=NA, "dr.integral.mc.ha.drmw"=NA, "dr.integral.mc.ha.drenw"=NA,
     "dr.direct.mc.ha.driptw"=NA, "dr.direct.mc.ha.driptw.trunc"=NA, "dr.direct.mc.ha.drow"=NA, "dr.direct.mc.ha.drmw"=NA, "dr.direct.mc.ha.drenw"=NA)

start.time <- Sys.time()
f.sim$dr.direct.mc.ha.driptw <- survival_inference(npath=npath, train.ueif.a1=f.train$dr.direct.mc.ha.driptw.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.driptw.ueif$ueif.a0,
                                                    estimation.ueif.a1=f.estimation$dr.direct.mc.ha.driptw.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.driptw.ueif$ueif.a0,
                                                    train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
f.sim$dr.direct.mc.ha.drow <- survival_inference(npath=npath, train.ueif.a1=f.train$dr.direct.mc.ha.drow.ueif$ueif.a1, train.ueif.a0=f.train$dr.direct.mc.ha.drow.ueif$ueif.a0,
                                                  estimation.ueif.a1=f.estimation$dr.direct.mc.ha.drow.ueif$ueif.a1, estimation.ueif.a0=f.estimation$dr.direct.mc.ha.drow.ueif$ueif.a0,
                                                  train.time=sort(unique(train.df$end.time)), estimation.time=sort(unique(estimation.df$end.time)))
end.time <- Sys.time()
end.time - start.time

save(f.sim, file="cancer_mortality_rv.RData")
