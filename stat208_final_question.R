## 208 Final
## Last Question
## Analyze data using methods in class
###########################################################################
library(dplyr)
library(readr)
library(janitor)
setwd("~/Desktop/Stat 208/")
dat<-read_csv("final_exam_data.csv")
y<-dat$V1
x<-(1:length(y))
###########################################################################
# EDA
plot(x,y)
par(mfrow=c(2,2))
hist(y,breaks=20) # multimodel distribution
plot(x,log(y)) # the dots fall very closely along this line
plot(log(x),y) # they fall closely along this one as well
plot(log(x),log(y))
# the linear relationship looks more appropriate than log transformations
###########################################################################
# Models
# linear, quadratic, log(x), exp(x), cosine(), NULL
# fit a model on the log transformed data as well?
m1<-lm(y~1+x)
m2<-lm(y~1+x+I(x^2))
m3<-lm(y~1+x+I(x^2)+factor(mod(x-1,12)))
m4<-lm(y~1+x+I(x^2)+cos((2*pi*(x))/(12))+sin((2*pi*(x))/(12)))
m5<-lm(y~1+x+I(x^2)+cos((2*pi*(x))/(6))+sin((2*pi*(x))/(6))) # change the period? didnt really do anything

dev.off()
plot(x,y,main="Data and Models 1-4")
abline(a=m1$coefficients[1],b=m1$coefficients[2],col="red")
points(x,m2$fitted.values,type="l",col="blue")
points(x,m3$fitted.values,type="l",col="orange")
points(x,m4$fitted.values,type="l",col="pink")
legend("topleft", 
       legend = c("Model 1", "Model 2", "Model 3", "Model 4"), 
       col = c("red","blue","orange","green"), 
       pch = 1, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
###########################################################################
# residual plots
par(mfrow=c(2,2))
plot(m1,main="Model 1")
par(mfrow=c(2,2))
plot(m2,main="Model 2")
par(mfrow=c(2,2))
plot(m3,main="Model 3")
par(mfrow=c(2,2))
plot(m4,main="Model 4")
# portmanteau test and chi-squared test of independence
m.step<-step(m4,direction="both")
# f-tests
anova(m1,m2,m3,m4)
###########################################################################
# conditioning on variance as we did in question 5
# just pull that code and update the design matrix
{
  N=nrow(dat)
  p2=3 # updated this
  Y=dat$V1
  tax=1:length(Y)
  # Define Design Matrix
  D2=matrix(0,N,p2)
  for(i in 1:N){
    D2[i,1]=1
    D2[i,2]=i
    D2[i,3]=i^2
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
R.sq3
# summary(m2)$r.squared 
# this lowers R^2 a bit
alpha <- 0.05
beta_hat_ci_upp <- beta.hat.gls + qt(1-alpha/2, N-p2)*sqrt(se3)
beta_hat_ci_low <- beta.hat.gls - qt(1-alpha/2, N-p2)*sqrt(se3)
beta_hat_ci_low[3];beta_hat_ci_upp[3] 
# insignificant quadratic term when we incorporate the error variance
dt(beta.hat.gls[3],df=N-p2)
par(mfrow=c(2,1))
plot(x,Resid3)
plot(x,m2$residuals)
# they look very similar
anova(m1,m2)
###########################################################################
# Questions:
# Linear vs Quadratic vs cosine model / seasonal effects?
# - quadratic w/ intercept term and slope is best. hope its not temp
# pattern in the residual plots?
# - maybe a slight one in model 2
# conditioning on variance as we did in question 5
# - this model says the quadratic term is not significant, but i like it anyway
