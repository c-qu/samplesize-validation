
library("readstata13"); library("pROC"); library("robustbase");library("scoring");
library("rms"); library("EnvStats"); library("safeBinaryRegression");library("MASS");
library("Matrix");library("speedglm"); library("sn");
library("pracma");require(ggplot2);require(gridExtra); library("openxlsx")
library("BaylorEdPsych");
library(plyr);library(ggpubr);
if(!require(devtools)) install.packages("devtools")

invlogit    <- function(x){exp(x)/(1+exp(x))}

################################################################################################
############compare empirical SE vs Approx. SE of performance measures##########################
################################################################################################

meas_true<- function(nevents,p,c,fc) {
  n       <- nevents/p
  sigmain <- sqrt(2)*qnorm(c)*fc
  
  mu      <-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
  sigma   <-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))
  
  eta     <- rnorm(n,mu,sigma) 
  pi      <- invlogit(eta) 
  y       <- rbinom(n,1,pi)
  
  n1      <- sum(y)
  n0      <- n-n1
  cstat <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c_est <- as.vector(cstat$auc)
  sigma0<-sd(eta[y==0]); sigma1<-sd(eta[y==1])
  
  sigma0<-round(sigma0,2)
  sigma1<-round(sigma1,2)
  out<-c(mean(y),c_est,mu,sigma,round(sigma0,2),round(sigma1,2))
  round(out,3)
}


meas<- function(nevents,p,c,fc) {
  
  n<-nevents/p
  
  sigmain <- sqrt(2)*qnorm(c)*fc
   mu     <-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
  sigma   <-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))
  
  
  eta <- rnorm(n,mu,sigma) 
  pi <- invlogit(eta) 
  y  <- rbinom(n,1,pi)

  
  cstat  <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c      <- as.vector(cstat$auc)
  cs_mod <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs     <- as.numeric(coef(cs_mod))[2]
  cl_mod <- speedglm(y ~ offset(eta), family=binomial(link='logit')) 
  cl     <- as.numeric(coef(cl_mod))[1]
  out<-c(mean(y),c,cs,cl)
  out
  
}


Nsim <-10000

res<-NULL

for (p in c(0.05,0.1,0.3)){
  for (c in c(0.64, 0.72, 0.8, 0.85, 0.9)){

    fc=1
    if (c==0.8)            fc=1.01
    if (p<0.3  & c==0.85)  fc=1.03
    if (p<0.3  & c==0.9)   fc=1.05
    if (p==0.3 & c==0.85)  fc=1.02
    if (p==0.3 & c==0.9)   fc=1.04
    
    for (nevents in c(50, 100, 200,400)){
      
    print(c(p,c,nevents))

    n<-nevents/p
    a<-replicate(Nsim, meas(nevents,p=p,c=c,fc=fc) )
    a<-t(a)
    

true   <- meas_true(300000,p,c,fc=fc)
p_true <- true[1]; c_true <- true[2]

var_emp<- apply(a,2,var,na.rm=TRUE)[2:4]

var_emp_c=var_emp[1]    #empirical SE of C
var_emp_cs=var_emp[2]   #empirical SE of CS
var_emp_cl=var_emp[3]   #empirical SE of CSL

#approx. SE of C-statistics using equation (6)
var_app_c   <- ((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3)))-c_true^2)/(n*p_true-n*p_true^2) 

#approx. SE of calibration slope using equation (7)
A           <- 2*p_true*(1-p_true)*qnorm(c_true)^2
var_app_cs  <-1/(A*n)+2/(n-2)   

sigmain     <- sqrt(2)*qnorm(c_true)
mu          <-0.5*(2*p_true-1)*(sigmain^2)+log(p_true/(1-p_true))
sigma       <-sqrt((sigmain^2)*(1+p_true*(1-p_true)*(sigmain^2)))
n_sim<-3000000
eta_sim           <-  rnorm(n_sim,mu,sigma)
prob_sim          <-  invlogit(eta_sim)
omega_sim         <-  prob_sim*(1-prob_sim)
omega_eta_sim     <-  omega_sim*eta_sim
omega_eta_sq_sim  <-  omega_sim*(eta_sim^2)
mean_omega        <- mean(omega_sim)
mean_omega_eta    <- mean(omega_eta_sim)
mean_omega_eta_sq <- mean(omega_eta_sq_sim) 

#approx. SE of calibration slope using equation (13)
numer       <-n*mean_omega
denom       <-n* mean_omega * n*mean_omega_eta_sq  - (n^2)*mean_omega_eta^2
var_app_cs2 <- numer/denom

#approx. SE of Calibration in the large using equation (8)
pin    <-exp(mu)/(1+exp(mu))
term1m <- n* pin*(1-pin)
term2m <- n* (1/2)*(1-6*pin+6*pin^2)*pin*(1-pin)*sigma^2
var_app_cl <- 1/ (term1m+term2m)


#approx. SE of Calibration in the large using equation (12)
var_app_cl2 <- 1/(n*mean_omega)


a_sum <- c(round(c(p_true,c_true),3),nevents,
           sqrt(var_emp_c), sqrt(var_app_c),(sqrt(var_app_c)/sqrt(var_emp_c)-1)*100,  
           sqrt(var_emp_cs), sqrt(var_app_cs),(sqrt(var_app_cs)/sqrt(var_emp_cs)-1)*100,
           sqrt(var_app_cs2),(sqrt(var_app_cs2)/sqrt(var_emp_cs)-1)*100,
           sqrt(var_emp_cl), sqrt(var_app_cl),(sqrt(var_app_cl)/sqrt(var_emp_cl)-1)*100,
           sqrt(var_app_cl2),(sqrt(var_app_cl2)/sqrt(var_emp_cl)-1)*100,
           true[5],true[6])
res<-rbind(res, a_sum)
    }}}

res_se<-res

colnames(res_se) <- c("p","c","n_events","se_emp_c","se_app_c","Bias_se_app",
                      "se_emp_cs","se_app_cs","Bias_se_app_cs",
                      "se_app_cs2","Bias_se_app_cs2",
                      "se_emp_cl","se_app_cl","Bias_se_app_cl",
                      "se_app_cl2","Bias_se_app_cl2","sigma0","sigma1")

View(res_se)



#################################################################################################
###########################events required to achieve target SE##################################
#################################################################################################

# approx. No of events vs true No of events required to achieve the SE of C-statistic = 0.0125,0.025,0.05 
cstat<- function(nevent,p,mu,sigma) {
  n     <- nevent/p
  eta   <- rnorm(n,mu,sigma) 
  p_est <- invlogit(eta) 
  y     <- rbinom(n,1,p_est)
  
  cstat <- pROC::roc(y,eta,quiet=TRUE,ci=FALSE)
  c_est <- as.vector(cstat$auc)
}


find_events_req_c <- function(se_true,p,mu,sigma,nstart1,nstart2,nstart_sim, tol=0.00001){
  
  #initial values
  n1<-nstart1
  n2<-nstart2
  nsim<-nstart_sim
  
  #set.seed(0)
  
  se1<-sd(replicate(nsim, cstat(n1,p,mu,sigma)))
  se2<-sd(replicate(nsim, cstat(n2,p,mu,sigma)))
  se1;se2
  
  d1<-(se1-se_true)
  d2<-(se2-se_true)
  d1;d2
  
  i<-1
  
  if (abs(d1)<abs(d2)) out<-c(n1,se1) else out<-c(n2,se2)
  
  miniter<-5
  
  while (min(abs(d1),abs(d2))>tol & i<10){
    
    #set.seed(i)
    
    if (abs(d1)>abs(d2) & sign(d1)!=sign(d2)) n1<-(n1+n2)/2 
    if (abs(d1)<abs(d2) & sign(d1)!=sign(d2)) n2<-(n1+n2)/2
    if (sign(d1)==sign(d2) & n1>n2)  {n1<-n1+5; n2<-n2-5}
    if (sign(d1)==sign(d2) & n1<n2)  {n1<-n1-5; n2<-n2+5}
    
    nsimnew<-nsim
    
    se1<-sd(replicate(nsimnew, cstat(n1,p,mu,sigma) ))
    se2<-sd(replicate(nsimnew, cstat(n2,p,mu,sigma) ))
    se1;se2
    d1=se1-se_true
    d2<-se2-se_true
    i<-i+1
    print(i-1)
    if (abs(d1)<abs(d2)) out<-c(n1,se1) else out<-c(n2,se2)
  }
  
  out
}


req_event_c <-NULL

for (se_true in c(0.0125,0.025,0.05)){
    for (p in c(0.05,0.1,0.3)){
        for (c in c(0.8,0.85,0.9)){
          fc=1
          if (c==0.8)            fc=1.01
          if (p<0.3  & c==0.85)  fc=1.03
          if (p<0.3  & c==0.9)   fc=1.05
          if (p==0.3 & c==0.85)  fc=1.02
          if (p==0.3 & c==0.9)   fc=1.04
          
      true            <- meas_true(1000000,p,c,fc); p_true<- true[1];c_true<-true[2];mu_true<-true[3];sigma_true<-true[4]
      
      events_req_app  <- ceiling((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3))-c_true^2)/((p_true-p_true^2)*se_true^2)*p_true) 
      
      
      if (c==0.64)   {nstart1 <- events_req_app*0.97 ; nstart2<-events_req_app*1.05}  
      if (c==0.72)   {nstart1 <- events_req_app*0.97 ; nstart2<-events_req_app*1.05}
      if (c==0.8)    {nstart1 <- events_req_app*0.97 ; nstart2<-events_req_app*1.05}
      if (c==0.85)   {nstart1 <- events_req_app*0.80 ; nstart2<-events_req_app*1.05}
      if (c==0.9)    {nstart1 <- events_req_app*0.7 ;  nstart2 <-events_req_app*1.5}
      
      
      nstart1<-ceiling(nstart1)
      nstart2<-ceiling(nstart2)
      
      error  <- NULL
      
      a <- tryCatch({find_events_req_c(se_true,p_true,mu_true,sigma_true,nstart1,nstart2,10000,tol=0.0001) 
      }, error=function(e){cat("error:",conditionMessage(e), "\n")})  
      if(inherits(a,"NULL")){
        error          <- error + 1
      }
      else{
        events_req<-ceiling(a[1])
        se_emp<- a[2]
        
        a_sum <- c(p_true,c_true,se_true, events_req, events_req_app, (events_req_app/events_req-1)*100, se_emp)
        
        req_event_c <-rbind(req_event_c, a_sum)
        print(c)
      }}}}


req_event_c <-data.frame(req_event_c) 

colnames(req_event_c) <- c("p","c","se_true", "events_req", "events_req_app", "bias", "se_emp_c")

##########################################################################################################
#approx. No of events vs true No of events required to achieve the SE of Calibration slope = 0.05,0.10,0.15 

cal<- function(nevents,p,mu,sigma) {
  n       <- nevents/p
  eta     <- rnorm(n,mu,sigma) 
  p_est   <- invlogit(eta) 
  y       <- rbinom(n,1,p_est)
  cs      <- coef(speedglm(y ~ eta, family=binomial(link='logit')))[2]
  cs
}



find_events_req_cs <- function(se_true,p,mu,sigma,nstart1,nstart2,nstart_sim, tol=0.0001){
  
  #initial values
  n1<-nstart1
  n2<-nstart2
  nsim<-nstart_sim
  
  #set.seed(0)
  
  se1<-sd(replicate(nsim, cal(n1,p,mu,sigma) ))
  se2<-sd(replicate(nsim, cal(n2,p,mu,sigma) ))
  se1;se2
  
  d1<-(se1-se_true)
  d2<-(se2-se_true)
  d1;d2
  
  i<-1
  
  if (abs(d1)<abs(d2)) out<-c(n1,se1) else out<-c(n2,se2)
  
  miniter<-5
  
  while (min(abs(d1),abs(d2))>tol){
    
    #set.seed(i)
    if (abs(d1)>abs(d2) & sign(d1)!=sign(d2)) n1<-(n1+n2)/2 
    if (abs(d1)<abs(d2) & sign(d1)!=sign(d2)) n2<-(n1+n2)/2
    if (sign(d1)==sign(d2) & n1>n2)  {n1<-n1+5; n2<-n2-5}
    if (sign(d1)==sign(d2) & n1<n2)  {n1<-n1-5; n2<-n2+5}
    
    nsimnew<-nsim
    se1<-sd(replicate(nsimnew, cal(n1,p,mu,sigma)))
    se2<-sd(replicate(nsimnew, cal(n2,p,mu,sigma)))
    se1;se2
    d1=se1-se_true
    d2<-se2-se_true
    i<-i+1
    print(i-1)
    if (abs(d1)<abs(d2)) out<-c(n1,se1) else out<-c(n2,se2)
  }
  out
}


req_event_cs<-NULL
for (se_true in c(0.05,0.1,0.15)){
         for (p in c(0.05, 0.1,0.3)){
             for (c in c(0.64,0.72,0.8,0.85,0.9)){
      fc=1
      if (c==0.8)            fc=1.01
      if (p<0.3  & c==0.85)  fc=1.03
      if (p<0.3  & c==0.9)   fc=1.05
      if (p==0.3 & c==0.85)  fc=1.02
      if (p==0.3 & c==0.9)   fc=1.04
      
      true            <- meas_true(1000000,p,c,fc); p_true<- true[1];c_true<-true[2];mu_true<-true[3];sigma_true<-true[4]
      #No of events required to achieve target SE using equation (6) 
      A              <- 2*p_true*(1-p_true)*qnorm(c_true)^2
      events_req_app <- ceiling(1/se_true^2 * (1/A+2)*p_true)
      
      #No of events required to achieve target SE using equation (13) 
      sigmain           <- sqrt(2)*qnorm(c_true)
      mu_est            <- 0.5*(2*p_true-1)*(sigmain^2)+log(p_true/(1-p_true))
      sigma_est         <- sqrt((sigmain^2)*(1+p_true*(1-p_true)*(sigmain^2)))
      eta_sim           <- rnorm(300000,mu_true,sigma_true)
      prob_sim          <- invlogit(eta_sim)
      omega_sim         <- prob_sim*(1-prob_sim)
      omega_eta_sim     <- omega_sim*eta_sim
      omega_eta_sq_sim  <- omega_sim*(eta_sim^2)
      mean_omega        <- mean(omega_sim)
      mean_omega_eta    <- mean(omega_eta_sim)
      mean_omega_eta_sq <- mean(omega_eta_sq_sim) 
      numer             <-  mean_omega
      denom             <-  mean_omega*mean_omega_eta_sq  - mean_omega_eta^2
      events_req_app_mc <-  (1/se_true^2)*(numer/denom)*p_true
    
      #find the true No of events to achieve target SE  
      if (c==0.64)   {nstart1 <- events_req_app*1.00 ; nstart2<-events_req_app*1.05}  
      if (c==0.72)   {nstart1 <- events_req_app*1.02 ; nstart2<-events_req_app*1.1}
      if (c==0.8)    {nstart1 <- events_req_app*1.1  ; nstart2<-events_req_app*1.25}
      if (c==0.85)   {nstart1 <- events_req_app*1.1  ; nstart2<-events_req_app*1.35}
      if (c==0.9)    {nstart1 <- events_req_app*1.2  ; nstart2<-events_req_app*1.45}
      
      nstart1<-ceiling(nstart1)
      nstart2<-ceiling(nstart2)
      
      a          <-find_events_req_cs(se_true,p_true,mu_true,sigma_true,nstart1,nstart2,15000, tol=0.0001)
      events_req <- ceiling(a[1])
      se_emp<-a[2]
      
      #result
      a_sum <- c(p_true,c_true,se_true, events_req, events_req_app, (events_req_app/events_req-1)*100,events_req_app_mc,(events_req_app_mc/events_req-1)*100,se_emp)
      
      req_event_cs<-rbind(req_event_cs, a_sum)
      print(c)
    }}}

req_event_cs           <-data.frame(req_event_cs)
colnames(req_event_cs) <- c("p","c","se_true", "events_req", "events_req_app","bias","events_req_app_mc","bias","se_emp_cs", "se_app")



#################################################################################################################
# approx. No of events vs true No of events required to achieve the SE of calibration in the large= 0.05,0.1,0.15 

cs_in_l<- function(nevents,p,mu,sigma){
    n    <- nevents/p
  eta    <- rnorm(n,mu,sigma) 
  pa     <- invlogit(eta) 
  y      <- rbinom(n,1,pa)
  csl    <- coef(speedglm(y ~ offset(eta), family=binomial(link='logit')))[1]
  csl
}


req_event_csl<-NULL

for (se_true in c(0.05, 0.1, 0.15)){
  
  for (p in c( 0.1,0.3)){
    
    for (c in c(0.85,0.9)){
      
      true            <- meas_true(1000000,p,c,fc); p_true<- true[1];c_true<-true[2];mu_true<-true[3];sigma_true<-true[4]
      sigmain<- sqrt(2)*qnorm(c_true)
      mu     <- 0.5*(2*p_true-1)*(sigmain^2)+log(p_true/(1-p_true))
      sigma  <- sqrt((sigmain^2)*(1+p_true*(1-p_true)*(sigmain^2)))
      
      #No of events required to achieve target SE using equation (8)
      pin    <-  exp(mu)/(1+exp(mu))
      term1m <-  pin*(1-pin)
      term2m <-  (1/2)*(1-6*pin+6*pin^2)*pin*(1-pin)*sigma^2
      events_req_app <-  ceiling(1/ (se_true^2*(term1m+term2m))*p)
      
      #No of events required to achieve target SE using equation (13)
      eta_sim           <-  rnorm(300000,mu,sigma)
      prob_sim          <-  invlogit(eta_sim)
      omega_sim         <-  prob_sim*(1-prob_sim)
      omega_eta_sim     <-  omega_sim*eta_sim
      omega_eta_sq_sim  <-  omega_sim*(eta_sim^2)
      mean_omega        <- mean(omega_sim)
      mean_omega_eta    <- mean(omega_eta_sim)
      mean_omega_eta_sq <- mean(omega_eta_sq_sim) 
      events_req_app_mc  <-1/(se_true^2*mean_omega)*p_true
      
      n0 <- events_req_app # Set start value to supplied lower bound
      i  <- 1 
      
      if (c==0.64)   {n1 <- ceiling(events_req_app*0.9) ; n2<-ceiling(events_req_app*1.1)}  
      if (c==0.72)   {n1 <- ceiling(events_req_app*0.9) ; n2<-ceiling(events_req_app*1.1)}
      if (c==0.8)    {n1 <- ceiling(events_req_app*0.9) ; n2<-ceiling(events_req_app*1.1)}
      if (c==0.85)   {n1 <- ceiling(events_req_app*1) ;  n2<-ceiling(events_req_app*1.35)}
      if (c==0.9)    {n1 <- ceiling(events_req_app*1) ;  n2<-ceiling(events_req_app*1.5)}
      
      
      f <- function(x) abs(se_true -sd(replicate(15000,cs_in_l(x,p,c))))
      # Check the upper and lower bounds to see if approximations result in 0
      d1 <- abs(se_true -sd(replicate(15000,cs_in_l(n1,p_true,mu_true,sigma_true))))
      d2 <- abs(se_true -sd(replicate(15000,cs_in_l(n2,p_true,mu_true,sigma_true))))
      
      d <- 1
      i <- 1
      while (i < 50 && d > 0.0001) {
        if (d1 > d2) {
          n1   <- ceiling((n1+n2)/2)
          se1  <- sd(replicate(15000,cs_in_l(n1,p_true,mu_true,sigma_true)))
          d1    <- abs(se_true -se1)
          d     <- d1
          events_req      <- n1
          se_emp <- se1
        }
        else{
          n2   <- ceiling((n1+n2)/2)
          se2  <- sd(replicate(15000,cs_in_l(n2,p_true,mu_true,sigma_true)))
          d2   <- abs(se_true -se2) 
          d    <- d2
          events_req    <- n2
          se_emp <- se2
        }
        print(i)
        i = i+1
        events_req
        se_emp
      }
      
      a_sum <- c(p_true,c_true,se_true, events_req, events_req_app,events_req_app_mc,events_req_app/events_req,events_req_app_mc/events_req ,se_emp)
      
      req_event_csl<-data.frame(rbind(req_event_csl, a_sum))
      print(c)
    }}}

colnames(req_event_csl) <- c("p","c","se_true", "events_req", "events_req_app", "events_req_app_mc", "events_req_app/events_req",  "events_req_app_mc/events_req", "se_emp_cs")


data.frame(rbind(req_event_c[,1:6],req_event_cs[,3:8],req_event_csl[,3:8]))


#########################################################################################################################################################################################
########################################################Power and type 1 error###########################################################################################################
#########################################################################################################################################################################################

power_c <- function(p0=for_c0[1],p1=for_c1[1],c0=for_c0[2],c1=for_c1[2],mu=for_c1[3], sigma=for_c1[4],alpha=0.05,beta=0.1) {
  
  #alpha = type I error, significance level
  #beta  = type-II error, 1-beta=power
  #H0: c=c0
  #H1: c=c1
  
  d=c1-c0
  
  sd_c0=sqrt((c0-2*T.Owen(-qnorm(c0),1/sqrt(3))-c0^2) /(p0-p0^2))
  sd_c1=sqrt((c1-2*T.Owen(-qnorm(c1),1/sqrt(3))-c1^2) /(p1-p1^2))
  
  #Formula for ss calculation
  
  n <- ceiling(((qnorm(1-alpha)*sd_c0 + qnorm(1-beta)*sd_c1))^2/d^2) ;n
  
  #Generate data under the alternative hypothesis
  
  eta    <- rnorm(n,mu,sigma)
  y      <- rbinom(n,1,invlogit(eta))
  
  c      <- roc(y,eta,quiet=TRUE,ci=TRUE)
  c_est  <- as.vector(c$auc) 
  se_c   <- (c$ci[3]-c_est)/qnorm(0.975)
  c_est+qnorm(1-alpha) *se_c
  cov_c  <- ifelse( (c_est-qnorm(1-alpha) *se_c) >= c0 , 1, 0)
  n1<-ceiling(n*p0)
  
  out<-c(mean(y),c0,c1,c_est,n1,cov_c)
  out
}

type1_error_c<- function(p0=for_c0[1],p1=for_c1[1],c0=for_c0[2],c1=for_c1[2],mu=for_c0[3], sigma=for_c0[4],alpha=0.05,beta=0.1) {
  
  #alpha = type I error, significance level
  #beta  = type-II error, 1-beta=power
  #H0: c=c0
  #H1: c=c1
  
  d=c1-c0
  
  sd_c0=sqrt((c0-2*T.Owen(-qnorm(c0),1/sqrt(3))-c0^2) /(p0-p0^2))
  sd_c1=sqrt((c1-2*T.Owen(-qnorm(c1),1/sqrt(3))-c1^2) /(p1-p1^2))
  
  n <- ceiling(((qnorm(1-alpha)*sd_c0 + qnorm(1-beta)*sd_c1))^2/d^2) ;n
  
  n1<-ceiling(n*p0)
  
  #Generate data under the NULL hypothesis
  
  eta  <- rnorm(n,mu,sigma)
  y    <- rbinom(n,1,invlogit(eta))
  
  c     <- roc(y,eta,quiet=TRUE,ci=TRUE)
  c_est <- as.vector(c$auc)     
  se_c  <- (c$ci[3]-c_est)/qnorm(0.975)
  
  cov_c <- ifelse( (c_est+qnorm(1-alpha) *se_c) <= c0 , 1, 0)
  out   <- c(mean(y),c0,c1,c_est,n1,cov_c)
  out
  
}


#################

Nsim <- 100000
res<-NULL

for (dc in c(0.03, 0.05)){
  
  for (p in c(0.05,0.1,0.3)){
    
    for (c in c(0.64, 0.72, 0.8, 0.85)){
      
      print(c(p,c))
      
      for_c0<- meas_true(2000000,p,c,fc=1)
      for_c1<- meas_true(2000000,p,c+dc,fc=1)    
      
      a<-replicate(Nsim, type1_error_c(p0=for_c0[1],p1=for_c1[1],c0=for_c0[2],c1=for_c1[2],mu=for_c0[3], sigma=for_c0[4],alpha=0.05,beta=0.1) ); a<-t(a)
      type1_error_n<-round(colMeans(a),3);type1_error_n
      
      b<-replicate(Nsim, power_c(p0=for_c0[1],p1=for_c1[1],c0=for_c0[2],c1=for_c1[2],mu=for_c1[3], sigma=for_c1[4],alpha=0.05,beta=0.1) ); b<-t(b)
      power_n<-round(colMeans(b),3);power_n
      
      a_sum <- c(dc,p,c,c+dc,type1_error_n[5:6],power_n[6])
      
      res<-rbind(res, a_sum)
    }}}
res_power <- res

res_power<-data.frame(res_power)


colnames(res_power) <- c("Difference","p","c0","c1","nevents", "alpha","power")
View(res_power)

