
#DGM2: Conditional normality with unequal variances

######################################################################################

library("readstata13"); library("pROC"); library("robustbase");library("scoring");
library("rms"); library("EnvStats"); library("safeBinaryRegression");library("MASS");
library("Matrix");library("speedglm"); library("BaylorEdPsych");library("sn");
library("pracma");require(ggplot2);require(gridExtra); library("openxlsx")
library(plyr);library(ggpubr);
if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
#source("Z:\\clustered data\\analysis_functions.R")
#source("Z:\\smmr_hvs\\roc.R")
#source("Z:\\clustered data\\analysis_functions.R")
#source("Z:\\smmr_hvs\\roc.R")

invlogit    <- function(x){exp(x)/(1+exp(x))}

######################################################################################

meas_true_unequal<- function(nevents,p,c,f0=0.67,f1=1.33,cmu=-0.1) {
  
  n1 <- nevents
  n0 <- round(round(nevents/p)-nevents)
  n<-n0+n1
    
  sigma <- sqrt(2)*qnorm(c)
  sigma1 <- 1.33*sigma
  sigma0 <- 0.67*sigma
  
  mu1 <- (sigma^2)/2 + log(n1/n0)
  mu0 <- mu1-sigma^2-cmu
  
  
  eta0 <- rnorm(n0,mu0,sigma1) 
  eta1 <- rnorm(n1,mu1,sigma0) 
  
  eta    <- c(eta0,eta1)
  prob <- invlogit(eta) 
  
  #Retrospective (gets C right)
  
  #if (prospective==FALSE) y_val  <- c(rep(0,n0), rep(1,n1))
  #n <-round(n_events/p)
  #pa_val           <- invlogit(eta)  
  
  
  #Prospective (gets calibration slope right)
  
  y <- as.vector(rbinom(n,1,prob=prob))  
  n1<-sum(y)
  
  c <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c_est <- as.vector(c$auc)

  cs    <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs_est <- as.numeric(coef(cs))[2]
  out<-c(mean(y),c_est,cs_est,mu0,mu1,sigma0,sigma1)
  round(out,3)
}


nevents=50000;p=0.05
meas_true_unequal(nevents,p,c=0.609,f0=0.67,f1=1.33,cmu=0.06)
meas_true_unequal(nevents,p,c=0.677,f0=0.67,f1=1.33,cmu=0.1)
meas_true_unequal(nevents,p,c=0.752,f0=0.67,f1=1.33,cmu=0.253)
meas_true_unequal(nevents,p,c=0.803,f0=0.67,f1=1.33,cmu=0.4)
meas_true_unequal(nevents,p,c=0.865,f0=0.67,f1=1.33,cmu=0.6)


nevents=50000;p=0.1
meas_true_unequal(nevents,p,c=0.612,f0=0.67,f1=1.33,cmu=0.06)
meas_true_unequal(nevents,p,c=0.682,f0=0.67,f1=1.33,cmu=0.1)
meas_true_unequal(nevents,p,c=0.76,f0=0.67,f1=1.33,cmu=0.2)
meas_true_unequal(nevents,p,c=0.815,f0=0.67,f1=1.33,cmu=0.3)
meas_true_unequal(nevents,p,c=0.875,f0=0.67,f1=1.33,cmu=0.5)

nevents=50000;p=0.3
meas_true_unequal(nevents,p,c=0.624,f0=0.67,f1=1.33,cmu=0.02)
meas_true_unequal(nevents,p,c=0.702,f0=0.67,f1=1.33,cmu=0.05)
meas_true_unequal(nevents,p,c=0.785,f0=0.67,f1=1.33,cmu=0.16)
meas_true_unequal(nevents,p,c=0.837,f0=0.67,f1=1.33,cmu=0.30)
meas_true_unequal(nevents,p,c=0.887,f0=0.67,f1=1.33,cmu=0.5)

meas_unequal<- function(nevents,p,c,f0=0.67,f1=1.33,cmu) {
  
  n1 <- nevents
  n0 <- round(round(nevents/p)-nevents)
  n<-n0+n1
  
  sigma <- sqrt(2)*qnorm(c)
  sigma1 <- 1.33*sigma
  sigma0 <- 0.67*sigma
  
  mu1 <- (sigma^2)/2 + log(n1/n0)
  mu0 <- mu1-sigma^2-cmu
  
  
  eta0 <- rnorm(n0,mu0,sigma1) 
  eta1 <- rnorm(n1,mu1,sigma0) 
  
  eta    <- c(eta0,eta1)
  prob <- invlogit(eta) 
  
  #Retrospective (gets C right)
  
  #if (prospective==FALSE) y_val  <- c(rep(0,n0), rep(1,n1))
  #n <-round(n_events/p)
  #pa_val           <- invlogit(eta)  
  
  
  #Prospective (gets calibration slope right)
  
  y <- as.vector(rbinom(n,1,prob=prob))  
  n1<-sum(y)
  
  cstat <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c     <- as.vector(cstat$auc)
  cs    <- speedglm(y ~ eta, family=binomial(link='logit')) 
  #cs    <- withWarnings(speedglm(Y_val ~ eta, family=binomial(link='logit')) )
  # if (length(calibration$error)==0) {
  #  calibration<-calibration$value
  # out<-summary(calibration)
  cs <- as.numeric(coef(cs))[2]
  out<-c(p,c,cs)
  out
  
}

meas_unequal(nevents,p,c=0.612,f0=0.67,f1=1.33,cmu=0.06)

Nsim <- 10000

res<-NULL

for (p in c(0.05,0.1,0.3)){
  
  for (c in c(0.64,0.72,0.8,0.85,0.9)){
    
    if (p==0.05 & c==0.64) {c1=0.609;cmu=0.06}
    if (p==0.05 & c==0.72) {c1=0.677;cmu=0.1}
    if (p==0.05 & c==0.80) {c1=0.752;cmu=0.253}
    if (p==0.05 & c==0.85) {c1=0.803;cmu=0.4}
    if (p==0.05 & c==0.90) {c1=0.865;cmu=0.6}
    
    if (p==0.1 & c==0.64) {c1=0.612;cmu=0.06}
    if (p==0.1 & c==0.72) {c1=0.682;cmu=0.1}
    if (p==0.1 & c==0.80) {c1=0.760;cmu=0.2}
    if (p==0.1 & c==0.85) {c1=0.815;cmu=0.3}
    if (p==0.1 & c==0.90) {c1=0.875;cmu=0.5}
    
    if (p==0.3 & c==0.64) {c1=0.624;cmu=0.02}
    if (p==0.3 & c==0.72) {c1=0.702;cmu=0.05}
    if (p==0.3 & c==0.80) {c1=0.785;cmu=0.16}
    if (p==0.3 & c==0.85) {c1=0.837;cmu=0.30}
    if (p==0.3 & c==0.90) {c1=0.887;cmu=0.50}
    
    for (nevents in c(50,100,200,400)){
      
    print(c(p,c,nevents))
      
true<-meas_true_unequal(100000,p=p,c=c1,f0=0.67,f1=1.33,cmu=cmu)
p_true <-true[1]; 
c_true <- true[2]

a<-replicate(Nsim, meas_unequal(nevents,p=p,c=c1,f0=0.67,f1=1.33,cmu=cmu) )
a<-t(a)
colMeans(a)

n<-ceiling(nevents/p_true)

var_emp<- apply(a,2,var,na.rm=TRUE)[2:3]

var_emp_c=var_emp[1]
var_emp_cs=var_emp[2]
A <- 2*p_true*(1-p_true)*qnorm(c_true)^2
var_app_cs <-1/(A*n)+2/(n-2)   
var_app_c=((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3)))-c_true^2)/(n*p_true-n*p_true^2) 

a_sum <- c(round(c(p_true,c_true),2),nevents,sqrt(var_emp_c), sqrt(var_app_c),sqrt(var_app_c)/sqrt(var_emp_c),  
           sqrt(var_emp_cs), sqrt(var_app_cs),sqrt(var_app_cs)/sqrt(var_emp_cs),true[6],true[7] )
res<-rbind(res, a_sum)
  }}}

res_se <- res

res_se<-data.frame(res_se)

res_se[,4:5]<-round(res_se[,4:5],4)
res_se[,7:8]<-round(res_se[,7:8],4)
res_se[,6]<-round(res_se[,6],2)
res_se[,9]<-round(res_se[,9],2)

colnames(res_se) <- c("p","c","n_events","se_emp_c","se_app_c","se_app_c/se_emp_c","se_emp_cs","se_app_cs","se_app_cs/se_emp_cs","sigma0","sigma1")
res_se_unequal<-res_se
View(res_se_unequal)

cd("C:\\Users\\Menelaos\\Dropbox\\Current Papers\\Validation\\Tables")

wb <- createWorkbook()
addWorksheet(wb, "se")
writeData(wb, 1, res_se)
addFilter(wb, 1, row = 1, cols = 1:ncol(res_se))
saveWorkbook(wb, file = "dgm3_se.xlsx", overwrite = TRUE)

