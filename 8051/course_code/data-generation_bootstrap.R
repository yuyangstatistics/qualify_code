
===================================
  
  Optimize a function

R commands: optimize, uniroot, nls (for computing nonlinear least squares)

Examples: 
  
  f <- function (x,a) (x-a)^2-0.5
xmin <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)
xmin


xmax <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3, maximum = T)
xmax

uniroot(f,c(0,1),tol=0.0001,a=1/4)

======================================
  
  A demonstrtion of Central Limit Theorem

1. the distribution of the sample mean from uniform distribution

density.estimate<-function(n,N){
  y<-rep(0,N)
  for (i in 1:N){
    x<-runif(n)
    y[i]<-mean(x)
  }
  par(mfrow=c(2,2))
  hist(runif(N))
  hist(y,main="Histogram of the sample mean")
  plot(density(y),main="Kernel density estimate of the sample mean")
  qqnorm(y)
}

density.estimate(1,1000)
density.estimate(2,1000)
density.estimate(10,1000)
density.estimate(100,1000)
density.estimate(200,1000)
density.estimate(500,1000)

2. the distribution of the sample mean from gamma distribution


density.estimate<-function(n,N){
  y<-rep(0,N)
  for (i in 1:N){
    x<-rgamma(n,shape=2)
    y[i]<-mean(x)
  }
  par(mfrow=c(2,2))
  hist(rgamma(N,shape=2))
  hist(y,main="Histogram of the sample mean")
  plot(density(y),main="Kernel density estimate of the sample mean")
  qqnorm(y)
}

density.estimate(1,1000)
density.estimate(2,1000)
density.estimate(10,1000)
density.estimate(100,1000)
density.estimate(200,1000)
density.estimate(500,1000)


3. the distribution of the sample mean from binomial distribution

density.estimate<-function(n,N){
  y<-rep(0,N)
  for (i in 1:N){
    x<-rbinom(n,size=10,prob=0.8)
    y[i]<-mean(x)
  }
  par(mfrow=c(2,2))
  hist(rbinom(N,size=10,prob=0.8))
  hist(y,main="Histogram of the sample mean")
  plot(density(y),main="Kernel density estimate of the sample mean")
  qqnorm(y)
}

density.estimate(1,1000)
density.estimate(2,1000)
density.estimate(10,1000)
density.estimate(100,1000)
density.estimate(200,1000)
density.estimate(500,1000)


==================================
  
  Monte Carlo integration

For one-dimensional integration, can use the command "integrate".

Suppose we want to find the integral of sin((x1+x2+x3)^2) on [0,1]^3.

Step 1: Generate 1000 copies of (X1,X2,X3) uniformlly distributed on [0,1]^3

Step 2: Find the mean and the standard error of sin((x1+x2+x3)^2)

f<-function(x1,x2,x3) sin((x1+x2+x3)^2)

integ<-function(f,N){
  x1<-runif(N)
  x2<-runif(N)
  x3<-runif(N)
  m<-mean(f(x1,x2,x3))
  se<-sqrt(var(f(x1,x2,x3)))/sqrt(N)   # could use sd directly for standard deviation
  CI<-c(m-se*qnorm(0.975),m+se*qnorm(0.975)) # 95% CI
  list(mean=m,standard.error=se,CI=CI)
}  

integ(f,1000)

integ(f,5000)

integ(f,10000)

===================================
  
  # Generate a linear regression data set with independent and uniformly distributed covariates
  
  data.reg<-function(n,d,beta,sig){
    e<-sig*rnorm(n)
    a<-runif(d*n)
    x<-matrix(a,nrow=n,ncol=d,byr=T)
    y<-beta[1]+x%*%beta[-1]+e
    list(y=y,x=x,beta=beta,n=n)
  }

a<-data.reg(100,6,c(1,1,2,3,0,0,0),1)

a<-data.reg(100,1,c(1,0.1),1)

par(mfrow=c(2,2))
a<-data.reg(100,1,c(1,0.3),0.1)
plot(a$x,a$y,main="Sigma = 0.1")
a<-data.reg(100,1,c(1,0.3),0.5)
plot(a$x,a$y,main="Sigma = 0.5")
a<-data.reg(100,1,c(1,0.3),1)
plot(a$x,a$y,main="Sigma = 1")
a<-data.reg(100,1,c(1,0.3),2)
plot(a$x,a$y,main="Sigma = 2")


par(mfrow=c(2,2))
a<-data.reg(100,1,c(1,3),0.1)
plot(a$x,a$y,main="Sigma = 0.1")
a<-data.reg(100,1,c(1,3),0.5)
plot(a$x,a$y,main="Sigma = 0.5")
a<-data.reg(100,1,c(1,3),1)
plot(a$x,a$y,main="Sigma = 1")
a<-data.reg(100,1,c(1,3),2)
plot(a$x,a$y,main="Sigma = 2")

# Generate a linear regression data set with independent and normally distributed covariates

data.reg<-function(n,d,beta,sig){
  e<-sig*rnorm(n)
  a<-rnorm(d*n)
  x<-matrix(a,nrow=n,ncol=d,byr=T)
  y<-beta[1]+x%*%beta[-1]+e
  list(y=y,x=x,beta=beta,n=n)
}

a<-data.reg(100,6,c(1,1,2,3,0,0,0),1)

a<-data.reg(100,1,c(1,0.1),1)

par(mfrow=c(2,2))
a<-data.reg(100,1,c(1,0.3),0.1)
plot(a$x,a$y,main="Sigma = 0.1")
a<-data.reg(100,1,c(1,0.3),0.5)
plot(a$x,a$y,main="Sigma = 0.5")
a<-data.reg(100,1,c(1,0.3),1)
plot(a$x,a$y,main="Sigma = 1")
a<-data.reg(100,1,c(1,0.3),2)
plot(a$x,a$y,main="Sigma = 2")


par(mfrow=c(2,2))
a<-data.reg(100,1,c(1,3),0.1)
plot(a$x,a$y,main="Sigma = 0.1")
a<-data.reg(100,1,c(1,3),0.5)
plot(a$x,a$y,main="Sigma = 0.5")
a<-data.reg(100,1,c(1,3),1)
plot(a$x,a$y,main="Sigma = 1")
a<-data.reg(100,1,c(1,3),2)
plot(a$x,a$y,main="Sigma = 2")

par(mfrow=c(2,2))
a<-data.reg(500,1,c(1,1),0.1)
plot(a$x,a$y,main="Sigma = 0.1")
a<-data.reg(500,1,c(1,1),0.5)
plot(a$x,a$y,main="Sigma = 0.5")
a<-data.reg(500,1,c(1,1),1)
plot(a$x,a$y,main="Sigma = 1")
a<-data.reg(500,1,c(1,1),2)
plot(a$x,a$y,main="Sigma = 2")

====================================
  
  #Bootstrap method
  
  library(boot)
help(boot) 

#Example 1: bootstrap CI for the mean

x<-rchisq(20,2)
m<-mean(x)

se.norm<-sd(x)/sqrt(length(x))  

se.norm

CI.norm<-c(m-qt(0.975,19)*se.norm,m+qt(0.975,19)*se.norm)

# Q: Is the CI reliable?

coverage<-function(n,N) {
  k=0
  for (i in 1:N) {
    x<-rchisq(n,2)
    m<-mean(x)
    se.norm<-sd(x)/sqrt(length(x))  
    CI.norm<-c(m-qt(0.975,n-1)*se.norm,m+qt(0.975,n-1)*se.norm)
    if (m-qt(0.975,n-1)*se.norm<=2&m+qt(0.975,n-1)*se.norm>=2) k=k+1    # mean of chi_2 is 2
  }
  k/N
}

coverage(20,1000) # the coverage is smaller than 95%

stat<-function(x,D) mean(x[D])   # D is a vector of indices  
a<-boot(x,stat,1000)

CI.boot<-c(quantile(a$t,0.025),quantile(a$t,0.975))

CI.norm

CI.boot

se.boot<-sd(a$t)

se.boot

# Q: Is the coverage better?

coverage1<-function(n,N) {
  k=0
  for (i in 1:N) {
    x<-rchisq(n,2)
    a<-boot(x,stat,5000)
    CI.boot<-c(quantile(a$t,0.025),quantile(a$t,0.975))
    if (quantile(a$t,0.025)<=2&quantile(a$t,0.975)>=2) k=k+1
  }
  k/N
}

coverage1(20,1000) # It may not necessarily work better at a small sample size.


# inference on the median
x<-rchisq(20,2)
med<-median(x)
med
qchisq(0.5,2)

smedian <- function(x, D) {
  return(median(x[D]))
}

a<-boot(x,smedian,1000)

CI.boot<-c(quantile(a$t,0.025),quantile(a$t,0.975))

CI.boot

mean(a$t)-med  # bootstrap estimate of the bias of the sample median

# Delta method and residual-based bootsrap for linear regression

datag<-data.reg(100,6,c(1,1,2,3,0,0,0),1)  # generate data; suppose we are interested in estimating theta = beta0 + beta1*beta3

out=lm(datag$y~datag$x)

thetahat=out$coef[1]+out$coef[2]*out$coef[4]  # point estimate of theta
der=c(1,out$coef[4],0,out$coef[2],0,0,0)  # gradient, estimated
se=sqrt(t(der)%*%vcov(out)%*%der)  # a standard error of thetahat
c(thetahat-qnorm(0.025,lower.tail =F),thetahat+qnorm(0.025,lower.tail =F))  # An asymptotic CI for theta based on delta method

res=residuals(out)               # get residuals
rres=res/sqrt(1-hatvalues(out))   # rescaled residuals
yhat=fitted(out)                  # fitted values

stat <- function(dat, D) {
  ystar=yhat+rres[D]
  outstar=lm(ystar~dat$x)
  thetastar=outstar$coef[1]+outstar$coef[2]*outstar$coef[4]
  thetastar
}

library(boot)
bs=boot(datag, statistic=stat, R=1000)
bs      # note that the "original" is not really the same as thetahat

c(quantile(bs$t,0.025),quantile(bs$t,0.975))   # a bootstrap 95% CI for theta

# Ordinary nonparametric bootstrap for linear regression

datf=data.frame(y=datag$y,x1=datag$x[,1],x2=datag$x[,2],x3=datag$x[,3],x4=datag$x[,4],x5=datag$x[,5],x6=datag$x[,6])

stat1 <- function(dat, D) {
  outstar=lm(y~x1+x2+x3+x4+x5+x6,subset=D,data=dat)
  thetastar=outstar$coef[1]+outstar$coef[2]*outstar$coef[4]
  thetastar
}

bs1=boot(datf, statistic=stat1, R=5000)

bs1     # note that the "original" here is the same as thetahat
c(quantile(bs1$t,0.025),quantile(bs1$t,0.975))