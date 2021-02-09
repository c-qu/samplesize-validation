#The linear predictor is conditionally Normally distributed given Y with the corresponding variances being equal. Data are generated from a LDA model under Assumption 2. 

library("readstata13");library("pROC");library("robustbase");library("scoring")
library("rms");library("EnvStats");library("safeBinaryRegression");library("MASS")
library("Matrix");library("speedglm");library("BaylorEdPsych");library("sn")
library("pracma");require(ggplot2);require(gridExtra);library("openxlsx")

#Compare the empirical SE and approx. SE of performance measures 

invlogit    <- function(x){exp(x)/(1+exp(x))}

variancetest<- function(nevents,p,truec,truecs) {


  sigma <- sqrt(2)*qnorm(c)
  
  n1 <- nevents
  n0 <- round(round(nevents/p)-nevents)
  n<-n0+n1

  mu1 <- (sigma^2)/2 + log(n1/n0)
  mu0 <- mu1-sigma^2
  
  eta0 <- rnorm(n0,mu0,sigma) 
  eta1 <- rnorm(n1,mu1,sigma) 
  
  eta    <- c(eta0,eta1)
  pa_val <- invlogit(eta) 
  y  <- c(rep(0,n0), rep(1,n1))

  cc <- roc(y,eta,quiet=TRUE,ci=TRUE)
  c <- as.vector(cc$auc)
  
  css    <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs  <- as.numeric(coef(css))[2]
  

  out<-c(round(p,2),truec,truecs,c,cs)
  out
  
}


Nsim <- 10000

res<-NULL


for (p in c(0.05,0.1,0.3)){
  
  for (c in c(0.64,0.72,0.8,0.85,0.9)){
    
    print(c(p,c))
    
    for (nevents in c(50,100,200,400)){
      
a<-matrix(NA,nrow=Nsim, ncol=5)

for (i in 1:Nsim){
  
set.seed(i)  
  
a[i,]<-variancetest(nevents,p,c,1)
  
}

var_emp<- apply(a,2,var,na.rm=TRUE)[4:5]
n<-nevents/p
var_emp_c=var_emp[1]
var_emp_cs=var_emp[2]
A <- 2*p*(1-p)*qnorm(c)^2
var_app_cs <-1/(A*n)+2/(n-2)   
var_app_c=((c-2*T.Owen(-qnorm(c),1/sqrt(3)))-c^2)/(n*p-n*p^2) 

a_sum <- c(p,c,nevents,sqrt(var_emp_c), sqrt(var_app_c),sqrt(var_app_c)/sqrt(var_emp_c),  
           sqrt(var_emp_cs), sqrt(var_app_cs),sqrt(var_app_cs)/sqrt(var_emp_cs) )

res<-rbind(res, a_sum)
  }}}

res_se <- res

res_se<-data.frame(res_se)

colnames(res_se) <- c("p","c","n_events","se_emp_c","se_app_c","se_app_c/se_emp_c","se_emp_cs","se_app_cs","se_app_cs/se_emp_cs")
View(res_se)


#Power and type1 error under dgm1

power_c<- function(p,c0,c1,alpha=0.05,beta=0.1) {
  
  #alpha = type I error, significance level
  #beta  = type-II error, 1-beta=power
  #H0: c=c0
  #H1: c=c1
  
  d=c1-c0
  
  sd_c0=sqrt((c0-2*T.Owen(-qnorm(c0),1/sqrt(3))-c0^2) /(p-p^2))
  sd_c1=sqrt((c1-2*T.Owen(-qnorm(c1),1/sqrt(3))-c1^2) /(p-p^2))
  
  n <- ceiling(((qnorm(1-alpha)*sd_c0 + qnorm(1-beta)*sd_c1))^2/d^2) ;n
  
  #Generarate data under alternative hypothesis
  sigma <- sqrt(2)*qnorm(c1)
  
  n1  <- ceiling(n*p); 
  n0  <- n-n1
  
  mu1 <- (sigma^2)/2 + log(n1/n0)
  mu0 <- mu1-sigma^2
  
  eta0 <- rnorm(n0,mu0,sigma) 
  eta1 <- rnorm(n1,mu1,sigma) 
  
  eta    <- c(eta0,eta1)

  y      <- c(rep(0,n0), rep(1,n1))
  
  c      <- roc(y,eta,quiet=TRUE,ci=TRUE)
  c_est  <- as.vector(c$auc)     
  se_c   <- (c$ci[3]-c_est)/qnorm(0.975)
  se_c
  cov_c  <- ifelse( (c_est-qnorm(1-alpha) *se_c) >= c0 , 1, 0)
 
  out<-c(p,c0,c1,c_est,n1,cov_c)
  out
}


type1_error_c<- function(p,c0,c1,alpha=0.05,beta=0.1) {
  
  #alpha = type I error, significance level
  #beta  = type-II error, 1-beta=power
  #H0: c=c0
  #H1: c=c1
  
  d=c1-c0
  
  se_c0=sqrt((c0-2*T.Owen(-qnorm(c0),1/sqrt(3))-c0^2) /(p-p^2))
  se_c1=sqrt((c1-2*T.Owen(-qnorm(c1),1/sqrt(3))-c1^2) /(p-p^2))
  
  n <- ceiling(((qnorm(1-alpha)*se_c0 + qnorm(1-beta)*se_c1))^2/d^2) ;n
  
  #Generarate data under the NULL hypothesis
  sigma <- sqrt(2)*qnorm(c0)
  n1 <- ceiling(n*p)  
  n0 <- n-n1
  n0;n1
  
  mu1 <- (sigma^2)/2 + log(n1/n0)
  mu0 <- mu1-sigma^2
  
  eta0 <- rnorm(n0,mu0,sigma) 
  eta1 <- rnorm(n1,mu1,sigma) 
  
  eta    <- c(eta0,eta1)
  pa_val <- invlogit(eta) 
  y  <- c(rep(0,n0), rep(1,n1))
  
  c     <- roc(y,eta,quiet=TRUE,ci=TRUE)
  c_est <- as.vector(c$auc)     
  se_c  <- (c$ci[3]-c_est)/qnorm(0.975)
  
  cov_c <- ifelse( (c_est+qnorm(1-alpha) *se_c) <= c0 , 1, 0)
  out <- c(p,c0,c1,c_est,n1,cov_c)
  out
  
}


Nsim <- 10000
power_type1 <-NULL

for (dc in c(0.03,0.05)){
  
  for (p in c(0.05, 0.1,0.3)){
    
    for (c in c(0.64, 0.72, 0.8, 0.85, 0.9)){
      
      print(c(p,c))
      
      a<-replicate(Nsim, type1_error_c(p=p,c0=c,c1=c+dc,alpha=0.05,beta=0.1) ); a<-t(a)
      type1_error_n<-round(colMeans(a),2);type1_error_n
      
      b<-replicate(Nsim, power_c(p=p,c0=c,c1=c+dc,alpha=0.05,beta=0.1) ); b<-t(b)
      power_n<-round(colMeans(b),2);power_n
      
      
      a_sum <- c(dc,p,c,c+dc,type1_error_n[5:6],power_n[6])
      
      power_type1<-rbind(power_type1, a_sum)
    }}}

res_power           <-data.frame(rpower_type1)
colnames(res_power) <- c("Difference","p","c0","c1","nevents", "alpha","power")
View(res_power)




