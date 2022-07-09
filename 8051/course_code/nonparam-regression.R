library(faraway)

# polynomial regression

data(exa)  # generated with f(x)=(sin(2*pi*x^3))^3
plot(y~x,exa,pch=".")
lines(m~x,exa)

m1=glm(y~poly(x,3,raw=T),data=exa) # raw=T means the usual polynomial basis instead of orthogonal 										polynomials
summary(m1)
lines(m1$fitted~x,exa,type="p")

orders<-function(data,M){       		# order selection by AIC
  crit<-rep(0,M)
  x<-data$x
  y<-data$y
  for (i in 1:M) {
    crit[i]<-glm(y~poly(x,i,raw=T))$aic
  }
  sel=min(which.min(crit))
  list(sel,crit)
}

ord=orders(exa,20)   # had convergence problems

orders<-function(data,M){
  crit<-rep(0,M)
  x<-data$x
  y<-data$y
  for (i in 1:M) {
    crit[i]<-glm(y~poly(x,i))$aic   # use orthogonal polynomials
  }
  sel=min(which.min(crit))
  list(sel,crit)
}

ord=orders(exa,25)   # no convergence issue with the orthogonal polynomials

m2=glm(y~poly(x,21,raw=T),data=exa) 

plot(y~x,exa,pch=".")
lines(m~x,exa)
lines(m2$fitted~x,exa,type="p")

# piecewise linear splines

rhs<-function(x,c) ifelse(x>c,x-c,0)
curve(rhs(x,0.5),0,1)
knots<-0:9/10
dm<-outer(exa$x,knots,rhs)  # outer product of dimension of x by that of knots, entry being the value of function rhs on the corresponding elements of x and knots
dm   # the design matrix
matplot(exa$x,dm,type="l")
pl<-glm(exa$y~dm)
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(pl))

# piecewise quadratic splines
library(splines)
matplot(bs(seq(0,1,length=1000),df=12,degree=2),type="l",ylab="")  # quadratic B-splines with 12 basis functions
pq<-glm(y~bs(x,12,degree=2),data=exa)
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(pq),type="p")

# piecewise cubic splines
matplot(bs(seq(0,1,length=1000),df=12),type="l",ylab="")  # cubic B-splines with 12 basis functions
pc<-glm(y~bs(x,12),data=exa)
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(pc),type="p")

matplot(ns(seq(0,1,length=1000),df=12),type="l",ylab="")  # natural spline with 12 basis functions
pcn<-glm(y~ns(x,12),data=exa)
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(pcn),type="p")

c(m1$aic,m2$aic,pl$aic,pq$aic,pc$aic,pcn$aic) # compare the models based on AIC

pc2<-glm(y~bs(x,knots=seq(0.6,0.9,length=6)),data=exa)  # choosing knots based on inspection
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(pc2),type="p")

c(m1$aic,m2$aic,pl$aic,pq$aic,pc$aic,pcn$aic,pc2$aic)
df<-length(exa$x)-c(m1$df.r,m2$df.r,pl$df.r,pq$df.r,pc$df.r,pcn$df.r,pc2$df.r)
c(m1$aic,m2$aic,pl$aic,pq$aic,pc$aic,pcn$aic,pc2$aic)+(log(length(exa$x))-2)*df  # From BIC 								perspective, the model with manually chosen knots works the best.

# smoothing spline

ss<-smooth.spline(exa$x,exa$y)
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lines(exa$x,fitted(ss),type="p")

predict(ss,x=c(0.5,1.5))   # prediction

data(exb)    # simulated data with mean function 0
ss2<-smooth.spline(exb$x,exb$y)
plot(y~x,exb,pch=".",xlab="x",ylab="y")
lines(m~x,exb)
lines(exb$x,fitted(ss2),lty=2)   # severe overfitting

ss3<-smooth.spline(exb$x,exb$y,cv=T)
plot(y~x,exb,pch=".",xlab="x",ylab="y")
lines(m~x,exb)
lines(exb$x,fitted(ss3),lty=2)   # the use of LOO CV instead of GCV does not help

# lowess

plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lo<-loess(y~x,exa)
lines(exa$x,fitted(lo),lty=2)  # overly smooth

plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(m~x,exa)
lo2<-loess(y~x,exa,span=0.2)
lines(exa$x,fitted(lo2),lty=2) # the reduced span of 0.2 is better than the default choice of 0.75

# kernel smoothing

par(mfrow=c(3,2))
layout(matrix(c(1:6), 3, 2, byrow = FALSE))
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"normal",0.1))   # Gaussian kernel
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"normal",0.5))
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"normal",1))

plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"box",0.1))   # box kernel or moving window averaging
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"box",0.5))
plot(y~x,exa,pch=".",xlab="x",ylab="y")
lines(ksmooth(exa$x,exa$y,"box",1))     # the Gaussian kernel gives smoother estiamtes

ksmooth(exa$x,exa$y,"normal",0.1,x.points=c(0.5,0.9,1.5))  # prediction at three new points; it does not work at x=1.5.

library(sm)
hm<-hcv(exa$x,exa$y,display="lines")  # smoothing parameter selection by LOO CV
sm.regression(exa$x,exa$y,h=hm,xlab="x",ylab="y") # estimation with the selected bandwidth

# Assessing a parametric model via comparison with a nonparametric estimate; a type of informal goodness of fit test

data(faithful)
plot(waiting~eruptions,faithful)
linear<-lm(waiting~eruptions,faithful)
lines(faithful$eruptions,fitted(linear),lty=2) # the fit seems to be fine
hm<-hcv(faithful$eruptions,faithful$waiting)
sm.regression(faithful$eruptions,faithful$waiting,h=hm,xlab="x",ylab="y") # the estimate seems somewhat different from the linear version. Is linear model OK? A formal comparison can be done via CV.

compare<-function(x,y,r,M){
  win=0
  n2=ceiling(length(x)*r)
  dat<-data.frame(x=x,y=y)
  h<-hcv(x,y)      # it is better to reselect the bandwidth each time below, but it 						may fail
  for (i in 1:M) {
    eva<-sort(sample(length(x),n2))
    lin<-lm(y~x,subset=-eva,dat)
    new=data.frame(x=x[eva])
    err.l<-sum((y[eva]-predict(lin,new))^2)
    #		h<-hcv(x[-eva],y[-eva])
    ker<-ksmooth(x[-eva],y[-eva],"normal",h,x.points=x[eva])  # Note that the new x values are ordered in the output and the predictions correspond to that ordering. But x[eva] is not ordered if we use sample(length(x),n2) without ordering. 
    err.k<-sum((y[eva]-ker$y)^2)
    if (err.k>err.l) win=win+1
  }
  win/M
}

compare(faithful$eruptions,faithful$waiting,0.25,100)  # it seems clear that the linear model is better for estimation

compare(exa$x,exa$y,0.25,100)   # the nonparametric method is way better in this case




