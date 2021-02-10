library("readstata13"); library("pROC"); library("robustbase");library("scoring");
library("rms"); library("EnvStats"); library("safeBinaryRegression");library("MASS");
library("Matrix");library("speedglm"); library("sn");
library("pracma");require(ggplot2);require(gridExtra); library("openxlsx")
library("BaylorEdPsych");
library(plyr);library(ggpubr);


invlogit    <- function(x){exp(x)/(1+exp(x))}

meas<- function(n,beta) {

  cov              <- cbind(1,rbinom(n,1,prob=0.1),rbinom(n,1,prob=0.2),rbinom(n,1,prob=0.3),rbinom(n,1,prob=0.5))
  eta              <- as.vector(cov%*%beta)
  mu               <- mean(eta)
  sigma            <- sd(eta)
  pa               <- invlogit(eta)                                    #probabilties for generating Y in the validating data 
  y                <- as.vector(rbinom(n,1,prob=pa)) 
  

  cstat  <- pROC::roc(y,eta,quiet=TRUE,ci=FALSE)
  auc    <- as.vector(cstat$auc)
  cs_mod <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs     <- as.numeric(coef(cs_mod))[2]
  csl_mod <- speedglm(y ~ offset(eta), family=binomial(link='logit')) 
  csl     <- as.numeric(coef(csl_mod))[1]
  out     <- c(mean(y),auc,cs,csl,mu,sigma)
  out
}

gerse_dgm4 <- function(n,beta){
  sim              <- replicate(10000, meas(n,beta)) 
  true             <- t(sim) 
  mean             <- apply(true,2,mean)[c(1,2,5,6)]
  p_true           <- mean[1]; c_true <- mean[2]; mu <-mean[3];sigma<-mean[4]
  print(c(c_true,p_true))
  #sample           <- replicate(Nsim, meas(n,beta))

  sigmain     <- sqrt(2)*qnorm(c_true)
  mu1         <-0.5*(2*p_true-1)*(sigmain^2)+log(p_true/(1-p_true))
  sigma1      <-sqrt((sigmain^2)*(1+p_true*(1-p_true)*(sigmain^2)))
    pin1      <- exp(mu1)/(1+exp(mu1))
  term1m1     <- n* pin1*(1-pin1)
  term2m1     <- n* (1/2)*(1-6*pin1+6*pin1^2)*pin1*(1-pin1)*sigma1^2
  var_app_cl1 <- 1/ (term1m1+term2m1)
  
  
  var_emp          <- apply(true,2,var,na.rm=TRUE)[2:4]
  var_emp_c        =  var_emp[1]
  var_emp_cs       =  var_emp[2] 
  var_emp_cl       =  var_emp[3]

  var_app_c        =  ((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3)))-c_true^2)/(n*p_true-n*p_true^2) 
  
  A              <- 2*p_true*(1-p_true)*qnorm(c_true)^2
  var_app_cs     <- 1/(A*n)+2/(n-2) 

  pin            <- exp(mu)/(1+exp(mu))
  term1m         <- n* pin*(1-pin)
  term2m         <- n* (1/2)*(1-6*pin+6*pin^2)*pin*(1-pin)*sigma^2
  var_app_cl     <- 1/ (term1m+term2m)
  
  nevents        <- n*p_true
  n_sim<-100000
  eta_sim           <-  rnorm(n_sim,mu,sigma)
  prob_sim          <-  invlogit(eta_sim)
  omega_sim         <-  prob_sim*(1-prob_sim)
  omega_eta_sim     <-  omega_sim*eta_sim
  omega_eta_sq_sim  <-  omega_sim*(eta_sim^2)
  mean_omega        <- mean(omega_sim)
  mean_omega_eta    <- mean(omega_eta_sim)
  mean_omega_eta_sq <- mean(omega_eta_sq_sim) 
  numer            <-  mean_omega
  denom            <-  n* (mean_omega *mean_omega_eta_sq  - mean_omega_eta^2)
  var_app_cs2      <- numer/denom
 
  #fisher
  var_app_cl2      <- 1/(n*mean_omega)
  
  
  eta_sim1           <-  rnorm(n_sim,mu1,sigma1)
  prob_sim1          <-  invlogit(eta_sim1)
  omega_sim1         <-  prob_sim1*(1-prob_sim1)
  omega_eta_sim1     <-  omega_sim1*eta_sim1
  omega_eta_sq_sim1  <-  omega_sim1*(eta_sim1^2)
  mean_omega1        <- mean(omega_sim1)
  mean_omega_eta1    <- mean(omega_eta_sim1)
  mean_omega_eta_sq1 <- mean(omega_eta_sq_sim1) 
  numer1            <-  mean_omega1
  denom1            <-  n* (mean_omega1 *mean_omega_eta_sq1  - mean_omega_eta1^2)
  var_app_cs3       <- numer1/denom1
  
  #fisher
  var_app_cl3      <- 1/(n*mean_omega1)
  
  
out       <- c(round(c(p_true,c_true),3),nevents,
             sqrt(var_emp_c),  sqrt(var_app_c),(sqrt(var_app_c)/sqrt(var_emp_c)-1)*100,  
             sqrt(var_emp_cs), sqrt(var_app_cs),(sqrt(var_app_cs)/sqrt(var_emp_cs)-1)*100,
             sqrt(var_app_cs2),(sqrt(var_app_cs2)/sqrt(var_emp_cs)-1)*100,
             sqrt(var_app_cs3),(sqrt(var_app_cs3)/sqrt(var_emp_cs)-1)*100,
             sqrt(var_emp_cl), sqrt(var_app_cl),(sqrt(var_app_cl)/sqrt(var_emp_cl)-1)*100,
             sqrt(var_app_cl1),(sqrt(var_app_cl1)/sqrt(var_emp_cl)-1)*100,
             sqrt(var_app_cl2),(sqrt(var_app_cl2)/sqrt(var_emp_cl)-1)*100,
             sqrt(var_app_cl3),(sqrt(var_app_cl3)/sqrt(var_emp_cl)-1)*100)
out
}


#p=5%
beta_c1p1       <- c(-3.6,1.1,0.5,0.5,0.3)
beta_c2p1       <- c(-4,1.3,1.1,1,0.3)
beta_c3p1       <- c(-4.6,2.3,1.4,1.2,0.4)
beta_c4p1       <- c(-4.92,3,1.4,1.2,0.4) #c=0.85,p=0.05
beta_c5p1       <- c(-5.37,4,1.5,1.2,0.4) #c=0.90,p=0.05
#p=10%
beta_c1p2       <- c(-2.8,1.1,0.5,0.6,0.3)
beta_c2p2       <- c(-3.3,1.5,1.1,0.9,0.4)
beta_c3p2       <- c(-3.8,2.4,1.4,1.3,0.4)
beta_c4p2       <- c(-4.2,3.3,1.6,1.5,0.4) #c=0.85,p=0.10
beta_c5p2       <- c(-4.72,4.5,1.8,1.5,0.4)#c=0.90,p=0.10

#P=30%
beta_c1p3  <- c(-1.4,1.2,0.5,0.6,0.3)
beta_c2p3  <- c(-1.85,1.8,1.2,0.9,0.5)
beta_c3p3  <- c(-2.4,2.8,1.6,1.5,0.6)
beta_c4p3  <- c(-2.8,4.2,2,2,0.6)   #c=0.85,p=0.30
beta_c5p3  <- c(-3.3,5.2,3.3,2,0.6) #c=0.90,p=0.30

result_p1 <- rbind(gerse_dgm4(1000,beta_c1p1),gerse_dgm4(2000,beta_c1p1),gerse_dgm4(4000,beta_c1p1),gerse_dgm4(8000,beta_c1p1),
                gerse_dgm4(1000,beta_c2p1),gerse_dgm4(2000,beta_c2p1),gerse_dgm4(4000,beta_c2p1),gerse_dgm4(8000,beta_c2p1),
                gerse_dgm4(1000,beta_c3p1),gerse_dgm4(2000,beta_c3p1),gerse_dgm4(4000,beta_c3p1),gerse_dgm4(8000,beta_c3p1),
                gerse_dgm4(1000,beta_c4p1),gerse_dgm4(2000,beta_c4p1),gerse_dgm4(4000,beta_c4p1),gerse_dgm4(8000,beta_c4p1),
                gerse_dgm4(1000,beta_c5p1),gerse_dgm4(2000,beta_c5p1),gerse_dgm4(4000,beta_c5p1),gerse_dgm4(8000,beta_c5p1))

result_p2 <- rbind(gerse_dgm4(500,beta_c1p2),gerse_dgm4(1000,beta_c1p2),gerse_dgm4(2000,beta_c1p2),gerse_dgm4(4000,beta_c1p2),
                   gerse_dgm4(500,beta_c2p2),gerse_dgm4(1000,beta_c2p2),gerse_dgm4(2000,beta_c2p2),gerse_dgm4(4000,beta_c2p2),
                   gerse_dgm4(500,beta_c3p2),gerse_dgm4(1000,beta_c3p2),gerse_dgm4(2000,beta_c3p2),gerse_dgm4(4000,beta_c3p2),
                   gerse_dgm4(500,beta_c4p2),gerse_dgm4(1000,beta_c4p2),gerse_dgm4(2000,beta_c4p2),gerse_dgm4(4000,beta_c4p2),
                   gerse_dgm4(500,beta_c5p2),gerse_dgm4(1000,beta_c5p2),gerse_dgm4(2000,beta_c5p2),gerse_dgm4(4000,beta_c5p2))

result_p3 <- rbind(gerse_dgm4(167,beta_c1p3),gerse_dgm4(333,beta_c1p3),gerse_dgm4(667,beta_c1p3),gerse_dgm4(1333,beta_c1p3),
                   gerse_dgm4(167,beta_c2p3),gerse_dgm4(333,beta_c2p3),gerse_dgm4(667,beta_c2p3),gerse_dgm4(1333,beta_c2p3),
                   gerse_dgm4(167,beta_c3p3),gerse_dgm4(333,beta_c3p3),gerse_dgm4(667,beta_c3p3),gerse_dgm4(1333,beta_c3p3),
                   gerse_dgm4(167,beta_c4p3),gerse_dgm4(333,beta_c4p3),gerse_dgm4(667,beta_c4p3),gerse_dgm4(1333,beta_c4p3),
                   gerse_dgm4(167,beta_c5p3),gerse_dgm4(333,beta_c5p3),gerse_dgm4(667,beta_c5p3),gerse_dgm4(1333,beta_c5p3))

rbind(result_p1,result_p2,result_p3)
