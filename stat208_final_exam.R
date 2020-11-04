########################################################################
## 208 final
########################################################################
## Problem 1
########################################################################

par(mfrow=c(2,1))
x<-seq(from=-2,to=2,length.out=100)
plot(x,pnorm(x),type="l",main="N(0,1) C.D.F.")

n=6
j=c(1:n)
qnorm((j-0.5)/n)
sum(qnorm((j-0.5)/n))
qnorm(1-(1/40))
qnorm(0.05-(1/40))

n=6
j=c(1:n)
par(mfrow=c(2,1))
x=seq(from=0,to=1,0.01)
plot(x,qnorm(x),type="l",main="N(0,1) Inverse C.D.F.")
abline(h=0,col="red")
abline(v=0.5,col="blue")
for(i in 1:n){
  abline(v=j[i]/n,col="green")
}
plot(x,qnorm(x),type="l",main="N(0,1) Inverse C.D.F. Location Transformation")
abline(h=0,col="red")
abline(v=0.5,col="blue")
for(i in 1:n){
  abline(v=(j[i]-0.5)/n,col="green")
}

## part (c)
dev.off()
t=19
l=c(1:t)
cos(2*pi*l/t)
cos(2*pi*l/t) %>% sum()
x=seq(from=0,to=1,0.01)
plot(x,cos(2*pi*x),type="l")
abline(h=0,col="red")
for(i in 1:t){
  abline(v=cos(l[i]/t),col="green")
}
for(i in 1:t){
  abline(v=cos(2*pi*l[i]/t),col="orange")
}

########################################################################
## Problem 5
########################################################################
# S12 mauna loa model
{
  setwd("~/Desktop/Stat 208/")
  library(dplyr)
  library(readr)
  library(numbers)
  library(readr)
  library(dplyr)
  library(janitor)
  
  loa <- read_csv("~/Desktop/Stat 208/MaunaLoaCO2.csv") %>% clean_names()
  d=60#2018-1959+1
  T=12
  N=d*T
  p2=14 # updated this
  
  # Define Observations and Times 
  # Initialize time (tax) and observation (Y) arrays
  tax=1:N
  tyr=1959+tax/12
  
  Y=loa$v3
  #Plot the data
  # Seasonal effect matrix
  S=matrix(0,12,11)
  diag(S)=rep(1,11)
  S[12,]=rep(-1,11)
  # Define Design Matrix
  D2=matrix(0,N,p2)
  for(i in 1:N){
    D2[i,1]=1
    D2[i,2]=i
    D2[i,3]=i^2
    if(mod(i-1,12) == 0) {
      D2[i:(i+11),4:14] = S
    }
  }
  
  # Compute the OLS estimators 
  H1.2=t(D2) %*% D2
  H2.2=solve(H1.2)
  H3.2= H2.2 %*% t(D2)
  theta_hat2= H3.2 %*% Y
  
  # Compute Estimated Y
  Yhat2= D2 %*% theta_hat2
  
  # Compute Residuals
  Resid2=Y-Yhat2
  
  # Compute Parameter Standard Errors
  se2=matrix(0,p2,1)
  SSE2=t(Resid2) %*% Resid2
  sighat2=SSE2/(N-p2)
  for(i in 1:p2){
    se2[i]=sighat2^(1/2)*H2.2[i,i]^(1/2)
  }
  
  R.sq2 <- 1-SSE2/sum((Y-mean(Y))^2)
}
# sample variance of the residuals as \hat{\sigma^2}
sigma.sq.hat<-var(Resid2)
# sample 1 lag correlation as the estimate of \hat{\rho}
rho.hat<-acf(Resid2,plot=F)$acf[2]
gamma.hat<-diag(N)
# create the \Gamma_0 matrix 
for(i in 1:nrow(gamma.hat)){
  for(j in 1:ncol(gamma.hat)){
    gamma.hat[i,j]=rho.hat^(abs(i-j))
  }
}
gamma.hat.inv<-solve(gamma.hat)
# GLM Diagnostics Slide 17
beta.hat.gls<-solve(t(D2)%*%gamma.hat.inv%*%D2)%*%t(D2)%*%gamma.hat.inv%*%Y
## report all parameter estimators and R^2. 
## Is the quadratic coefficient significantly positive? Report a p-value
Yhat3=D2%*%beta.hat.gls
Resid3=Y-Yhat3
se3=matrix(0,p2,1)
SSE3=t(Resid3) %*% Resid3
sighat3=SSE3/(N-p2)
for(i in 1:p2){
  se3[i]=sighat3^(1/2)*H2.2[i,i]^(1/2)
}

R.sq3 <- 1-SSE3/sum((Y-mean(Y))^2)
# parameter estimates
beta.hat.gls
# R^2
R.sq3
# confidence interval for quadratic term
alpha <- 0.05
beta_hat_ci_upp <- beta.hat.gls + qt(1-alpha/2, N-p2)*sqrt(se3)
beta_hat_ci_low <- beta.hat.gls - qt(1-alpha/2, N-p2)*sqrt(se3)
beta_hat_ci_low[3];beta_hat_ci_upp[3] # insignificant
# p-value for quadratic term
dt(beta.hat.gls[3],df=N-p2)
# weak evidence to reject the null hypothesis that the quadratic term = 0


########################################################################
## Problem 7
########################################################################
y1<-c(8.42,14.68,21.42,25.45,27.14,30.53,34.51,34.54,33.24,39.63,43.98,47.77)
y2<-c(9.86,9.54,11.96,12.46,11.38,14.69,16.48,20.11)
y3<-c(6.52,5.11,7.75,6.84,7.65,9.49,7.03,9.41,12.01)
x1<-c(1,3,5,6,7,8,9,9,10,11,12,14)
x2<-c(3,3,4,5,6,8,9,12)
x3<-c(2,5,7,8,10,15,16,18,20)

n1=length(y1)
n2=length(y2)
n3=length(y3)

mu1<-matrix(c(rep(1,n1),rep(0,n2),rep(0,n3)),ncol=1)
mu2<-matrix(c(rep(0,n1),rep(1,n2),rep(0,n3)),ncol=1)
mu3<-matrix(c(rep(0,n1),rep(0,n2),rep(1,n3)),ncol=1)
a1<-matrix(c(x1,rep(0,n2),rep(0,n3)),ncol=1)
a2<-matrix(c(rep(0,n1),x2,rep(0,n3)),ncol=1)
a3<-matrix(c(rep(0,n1),rep(0,n2),x3),ncol=1)

D<-cbind(mu1,mu2,mu3,a1,a2,a3)

n<-n1+n2+n3
y<-matrix(c(y1,y2,y3))
r<-3

# complete null model, just an intercept term
x.o<-matrix(rep(1,n))
m.o<-x.o%*%solve(t(x.o)%*%x.o)%*%t(x.o)
m.a<-D%*%solve(t(D)%*%D)%*%t(D)
F.stat<-((t(y)%*%(m.a-m.o)%*%y)/(r-1))*((n-r)/(t(y)%*%(diag(n)-m.a)%*%y))
F.stat
qf(p=0.95,df1=r-1,df2=n-r)
F.stat > qf(p=0.95,df1=r-1,df2=n-r) # If true, reject Ho

# testing all mu = 0
# just switch out null model design matrix
x.o=cbind(a1,a2,a3)
m.o<-x.o%*%solve(t(x.o)%*%x.o)%*%t(x.o)
m.a<-D%*%solve(t(D)%*%D)%*%t(D)
F.stat<-((t(y)%*%(m.a-m.o)%*%y)/(r-1))*((n-r)/(t(y)%*%(diag(n)-m.a)%*%y))
F.stat
qf(p=0.95,df1=r-1,df2=n-r)
F.stat > qf(p=0.95,df1=r-1,df2=n-r) 

# testing all alpha = 0
# just switch out null model design matrix
x.o=cbind(mu1,mu2,mu3)
m.o<-x.o%*%solve(t(x.o)%*%x.o)%*%t(x.o)
m.a<-D%*%solve(t(D)%*%D)%*%t(D)
F.stat<-((t(y)%*%(m.a-m.o)%*%y)/(r-1))*((n-r)/(t(y)%*%(diag(n)-m.a)%*%y))
F.stat
qf(p=0.95,df1=r-1,df2=n-r)
F.stat > qf(p=0.95,df1=r-1,df2=n-r) 
