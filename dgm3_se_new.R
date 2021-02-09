
library("readstata13"); library("pROC"); library("robustbase");library("scoring");
library("rms"); library("EnvStats"); library("safeBinaryRegression");library("MASS");
library("Matrix");library("speedglm"); library("sn");
library("pracma");require(ggplot2);require(gridExtra); library("openxlsx")
library("BaylorEdPsych");
library(plyr);library(ggpubr);
if(!require(devtools)) install.packages("devtools")

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
p=0.1
c=0.7
fc=1
true <- meas_true(100000,p,c,fc=fc); true

meas<- function(nevents,p,c,fc=1,fp=1) {
  
  n<-nevents/p
  
  sigmain <- sqrt(2)*qnorm(c)*fc
  #mu=logit(p)*fp
  
  mu<-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
  sigma<-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))
  
  
  eta <- rnorm(n,mu,sigma) 
  p <- invlogit(eta) 
  y <- rbinom(n,1,p)
  #n1<-sum(y)
  #n0<-n-n1
  
  cstat  <- roc(y,eta,quiet=TRUE,ci=FALSE)
  c      <- as.vector(cstat$auc)
  cs_mod <- speedglm(y ~ eta, family=binomial(link='logit')) 
  cs     <- as.numeric(coef(cs_mod))[2]
  cl_mod <- speedglm(y ~ offset(eta), family=binomial(link='logit')) 
  cl     <- as.numeric(coef(cl_mod))[1]
  out<-c(mean(p),c,cs,cl)
  out
  
}

#p<-0.5
#meas(100000,p=0.2,c=0.9,fc=1.055)
#p*100

Nsim <-5000

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
colMeans(a)

true <- meas_true(300000,p,c,fc=fc)


p_true <- true[1]; c_true <- true[2]

var_emp<- apply(a,2,var,na.rm=TRUE)[2:4]

var_emp_c=var_emp[1]
var_emp_cs=var_emp[2]
var_emp_cl=var_emp[3]
A <- 2*p_true*(1-p_true)*qnorm(c_true)^2
var_app_cs <-1/(A*n)+2/(n-2) 

sigmain <- sqrt(2)*qnorm(c)*fc
mu<-0.5*(2*p-1)*(sigmain^2)+log(p/(1-p))
sigma<-sqrt((sigmain^2)*(1+p*(1-p)*(sigmain^2)))

n_sim<-3000000
eta_sim           <-  rnorm(n_sim,mu,sigma)
prob_sim          <-  invlogit(eta_sim)
omega_sim         <-  prob_sim*(1-prob_sim)
omega_eta_sim     <-  omega_sim*eta_sim
omega_eta_sq_sim  <-  omega_sim*(eta_sim^2)
mean_omega        <- mean(omega_sim)
mean_omega_eta    <- mean(omega_eta_sim)
mean_omega_eta_sq <- mean(omega_eta_sq_sim) 


#Calibration in the large 

var_app_cl2 <- 1/(n*mean_omega)

pin<-exp(mu)/(1+exp(mu))
term1m <- n* pin*(1-pin)
term2m <- n* (1/2)*(1-6*pin+6*pin^2)*pin*(1-pin)*sigma^2
term3m<-  n* (1/24)*(1-30*pin+150*pin^2-240*pin^3+120*pin^4)*pin*(1-pin)*3*sigma^4
term1m+term2m+term3m
term1m+term2m
var_app_cl <- 1/ (term1m+term2m)

#Fisher
var_app_cl2 <- 1/(n*mean_omega)


#Gareth

a <- exp(mu)/(exp(mu)+1)^2 

b <- exp(mu)*(exp(2*mu)-4*exp(mu)+1)/(2*(exp(mu)+1)^4)
var_app_cl = 1/(n*(a+b*sigma^2))



#Chen

a_2              <- exp(-mu)/(1+exp(-mu))^2+exp(-mu)*(exp(-mu)-1)*(-mu)/(1+exp(-mu))^3 + (exp(-3*mu)-4*exp(-2*mu)+exp(-mu))*mu^2/(2*(1+exp(-mu))^4) 
b_2              <- exp(-mu)*(exp(-mu)-1)/(1+exp(-mu))^3- mu*(exp(-3*mu)-4*exp(-2*mu)+exp(-mu))/(1+exp(-mu))^4 
c_2              <- (exp(-3*mu)-4*exp(-2*mu)+exp(-mu))/(2*(1+exp(-mu))^4)
var_app_cs2      <- abs((1/n)*(a_2+b_2*mu+c_2*(mu^2+sigma^2))/((a_2+b_2*mu+c_2*(mu^2+sigma^2))*(a_2*(mu^2+sigma^2)+b_2*(mu^3+3*mu*sigma^2)+c_2*(3*sigma^4+mu^4+6*mu^2*sigma^2))-(a_2*mu+b_2*(mu^2+sigma^2)+c_2*(mu^3+3*mu*sigma^2))^2)) 


a_3              <- exp(-mu)/(1+exp(-mu))^2+exp(-mu)*(exp(-mu)-1)*(-mu)/(1+exp(-mu))^3 + (exp(-3*mu)-4*exp(-2*mu)+exp(-mu))*mu^2/(2*(1+exp(-mu))^4)-mu^3*(exp(-4*mu)-11*exp(-3*mu)+11*exp(-2*mu)-exp(-mu))/(6*(1+exp(-mu))^5)
b_3              <- exp(-mu)*(exp(-mu)-1)/(1+exp(-mu))^3- mu*(exp(-3*mu)-4*exp(-2*mu)+exp(-mu))/(1+exp(-mu))^4+3*mu^2*(exp(-4*mu)-11*exp(-3*mu)+11*exp(-2*mu)-exp(-mu))/(6*(1+exp(-mu))^5)
c_3              <- (exp(-3*mu)-4*exp(-2*mu)+exp(-mu))/(2*(1+exp(-mu))^4)-3*mu*(exp(-4*mu)-11*exp(-3*mu)+11*exp(-2*mu)-exp(-mu))/(6*(1+exp(-mu))^5)
d_3              <- (exp(-4*mu)-11*exp(-3*mu)+11*exp(-2*mu)-exp(-mu))/(6*(1+exp(-mu))^5)
var_app_cs2      <- abs((1/n)*(a_3+b_3*mu+c_3*(mu^2+sigma^2)+d_3*(mu^3+3*mu*sigma^2))/((a_3+b_3*mu+c_3*(mu^2+sigma^2)+d_3*(mu^3+3*mu*sigma^2))*(a_3*(mu^2+sigma^2)+b_3*(mu^3+3*mu*sigma^2)+c_3*(3*sigma^4+mu^4+6*mu^2*sigma^2)+d_3*(10*mu^3*sigma^2+15*mu*sigma^4+mu^5))-(a_3*mu+b_3*(mu^2+sigma^2)+c_3*(mu^3+3*mu*sigma^2)+d_3*(3*sigma^4+mu^4+6*mu^2*sigma^2))^2))


#var_app_cl2<- 1/(n*(a_3+b_3*mu+c_3*(mu^2+sigma^2)+d_3*(mu^3+3*mu*sigma^2)))

#var_app_cl2<- 1/(n*(a_2+b_2*mu+c_2*(mu^2+sigma^2)))





#1/ (term1m+term2m)
#1/(n*(a_2+b_2*mu+c_2*(mu^2+sigma^2)))


#s<-sigma
#var_app_cl2 <- 1/(n/(4 + mu^2 + sigma^2 + 1/12*(mu^4 + 6*mu^2*sigma^2 + 3*sigma^4)))
#term<-1/sqrt(n/(4 + mu^2 + s^2 + 1/12*(mu^4 + 6*mu^2*s^2 + 3*s^4)))
#var_app_cl2<- term^2


#a = exp(mu)/(exp(mu)+1)^2 
#b = exp(mu)*(exp(2*mu)-4*(exp(mu)+1)/(2*(exp(mu)+1)^4))

#se = 1/sqrt(n*(a + b*s^2))

#var_app_cl2 <- se^2

#Calibration slope

numer <-  n*mean_omega
denom <-  n* mean_omega * n*mean_omega_eta_sq  - (n^2)*mean_omega_eta^2


#approx
term1<-term1m+term2m+term3m
term2<-p*(1-p)*log(p/1-p)


#fisher
var_app_cs2=numer/denom


var_app_c=((c_true-2*T.Owen(-qnorm(c_true),1/sqrt(3)))-c_true^2)/(n*p_true-n*p_true^2) 

a_sum <- c(round(c(p_true,c_true),3),nevents,
           sqrt(var_emp_c), sqrt(var_app_c),(1-sqrt(var_app_c)/sqrt(var_emp_c))*100,  
           sqrt(var_emp_cs), sqrt(var_app_cs),(1-sqrt(var_app_cs)/sqrt(var_emp_cs)*100,
           sqrt(var_app_cs2),(1-sqrt(var_app_cs2)/sqrt(var_emp_cs)*100,
           sqrt(var_emp_cl), (1-sqrt(var_app_cl),sqrt(var_app_cl)/sqrt(var_emp_cl),
           sqrt(var_app_cl2),sqrt(var_app_cl2)/sqrt(var_emp_cl),
           true[5],true[6])
res<-rbind(res, a_sum)
  }}}

res_se <- res

res_se<-data.frame(res_se)

res_se[,4:16]<-round(res_se[,4:16],4)

res_se[,6]<-round(res_se[,6],2)
res_se[,9]<-round(res_se[,9],2)
res_se[,11]<-round(res_se[,11],2)
res_se[,14]<-round(res_se[,14],2)
res_se[,16]<-round(res_se[,16],2)

colnames(res_se) <- c("p","c","n_events","se_emp_c","se_app_c","se_app_c/se_emp_c",
                      "se_emp_cs","se_app_cs","se_app_cs/se_emp_cs",
                      "se_app_cs2","se_app_cs2/se_emp_cs",
                      "se_emp_cl","se_app_cl","se_app_cl/se_emp_cl",
                      "se_app_cl2","se_app_cl2/se_emp_cl","sigma0","sigma1")
res_se_marg_norm<-res_se
View(res_se_marg_norm)

cd("C:\\Users\\Menelaos\\Dropbox\\Current Papers\\Validation\\Tables\\Tables_se")

wb <- createWorkbook()
addWorksheet(wb, "se")
writeData(wb, 1, res_se)
addFilter(wb, 1, row = 1, cols = 1:ncol(res_se))
#saveWorkbook(wb, file = "dgm3_se_revise.xlsx", overwrite = TRUE)

