
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

################################################################################################

# Simulation Standard error

meas_true<- function(nevents,p,c,fc=1,fp=1) {
  n<-nevents/p
  sigmain <- sqrt(2)*qnorm(c)*fc
  #mu=logit(p)*fp
  
  mu<-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
  sigma<-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))
  
  eta <- rnorm(n,mu,sigma) 
  p_est <- invlogit(eta) 
  y <- rbinom(n,1,p_est)
  
  n1<-sum(y)
  n0<-n-n1
  cstat <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c_est <- as.vector(cstat$auc)
  sigma0<-sd(eta[y==0]); sigma1<-sd(eta[y==1])
  
  sigma0<-round(sigma0,2)
  sigma1<-round(sigma1,2)
  out<-c(mean(p_est),c_est,mu,sigma,round(sigma0,2),round(sigma1,2))
  round(out,3)
}

meas<- function(nevents,p,c,fc=1,fp=1) {
  
  n<-nevents/p
  
  sigmain <- sqrt(2)*qnorm(c)*fc
  #mu=logit(p)*fp
  
  mu<-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
  sigma<-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))
  
  
  eta <- rnorm(n,mu,sigma) 
  p <- invlogit(eta) 
  y <- rbinom(n,1,p)
  n1<-sum(y)
  n0<-n-n1
  
  cstat <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c     <- as.vector(cstat$auc)
  cs    <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs <- as.numeric(coef(cs))[2]
  out<-c(mean(p),c,cs)
  out
  
}

#p<-0.5
#meas(100000,p=0.2,c=0.9,fc=1.055)
#p*100

Nsim <-10000

res<-NULL

for (p in c(0.05,0.1,0.3)){
  
  for (c in c(0.64,0.72,0.8,0.85,0.9)){
    
    fc=1
    if (c==0.85) fc=1.025
    if (c==0.9) fc=1.045
    
    for (nevents in c(50,100,200,400)){
      
    print(c(p,c,nevents))

n<-nevents/p
a<-replicate(Nsim, meas(nevents,p=p,c=c,fc=fc) )
a<-t(a)
colMeans(a)

true<-meas_true(100000,p,c,fc=fc)

p_true <- true[1]; c_true <- true[2]

var_emp<- apply(a,2,var,na.rm=TRUE)[2:3]

var_emp_c=var_emp[1]
var_emp_cs=var_emp[2]
A <- 2*p_true*(1-p_true)*qnorm(c_true)^2
var_app_cs <-1/(A*n)+2/(n-2) 
var_app_cs <-1/(A*n)+2/(n-2) 

sigmain <- sqrt(2)*qnorm(c)*fc
mu<-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
sigma<-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))

nsim<-1000000
eta_sim<-rnorm(n_sim,mu,sigma)
prob_sim         <- invlogit(eta_sim)
omega            <- prob_sim*(1-prob_sim)
mean_omega       <- mean(omega)
mean_omega_eta   <- mean(omega*eta_sim)
mean_omega_etasq <- mean(omega*eta_sim^2)


numer <- N*mean(omega)
denom <- N^2 *mean_omega_eta_sq*mean_omega -N^2*mean_omega_eta^2

var_app_cs2=numer/denom


var_app_c=((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3)))-c_true^2)/(n*p_true-n*p_true^2) 

a_sum <- c(round(c(p_true,c_true),2),nevents,sqrt(var_emp_c), sqrt(var_app_c),sqrt(var_app_c)/sqrt(var_emp_c),  
           sqrt(var_emp_cs), sqrt(var_app_cs),sqrt(var_app_cs)/sqrt(var_emp_cs),sqrt(var_app_cs2),sqrt(var_app_cs2)/sqrt(var_emp_cs),true[5],true[6] )
res<-rbind(res, a_sum)
  }}}

res_se <- res

res_se<-data.frame(res_se)

res_se[,4:5]<-round(res_se[,4:5],4)
res_se[,7:8]<-round(res_se[,7:8],4)
res_se[,6]<-round(res_se[,6],2)
res_se[,9]<-round(res_se[,9],2)
res_se[,10:11]<-round(res_se[,10:11],4)

colnames(res_se) <- c("p","c","n_events","se_emp_c","se_app_c","se_app_c/se_emp_c","se_emp_cs","se_app_cs","se_app_cs/se_emp_cs","se_emp_cs","se_app_cs2","se_app_cs2/se_emp_cs","sigma0","sigma1")
res_se_marg_norm<-res_se
View(res_se_marg_norm)

#cd("C:\\Users\\Menelaos\\Dropbox\\Current Papers\\Validation\\Tables\\Tables_se")

#wb <- createWorkbook()
#addWorksheet(wb, "se")
#writeData(wb, 1, res_se)
#addFilter(wb, 1, row = 1, cols = 1:ncol(res_se))
#saveWorkbook(wb, file = "dgm3_se.xlsx", overwrite = TRUE)

