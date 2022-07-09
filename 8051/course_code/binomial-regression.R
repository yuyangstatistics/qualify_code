# first install the library "faraway"

library(faraway)

# Challenger O-rings data analysis

data(orings)
orings
attach(orings)
plot(damage/6~temp,xlim=c(25,85),ylim=c(0,1),xlab="Temp", ylab="Prop of damage")

m1<-glm(cbind(damage,6-damage)~temp,family=binomial)
summary(m1)

x <- seq(25,85,1)
lines(x,ilogit(11.6630-0.2162*x)) # estimated ilogit curve


## Prediction
fitted(m1)
predict(m1, type='response')
m1$fitted.values
ilogit(m1$linear.predictors)
# prediction at 31
ilogit(11.6630-0.2162*31)


## Inference and Goodness of fit

# Is temp really needed? Test based on deviance is more reliable than normal approximation
pchisq(38.898-16.912,1,lower=FALSE)

# goodness of fit test: is the model with intercept and temp adequate?
pchisq(16.912,21,lower=FALSE)

# goodness of fit test: is the model with intercept only adequate?
pchisq(38.898,22,lower=FALSE)

# 95% CI for coefficient of temp
c(-0.21623-1.96*0.05318, -0.21623+1.96*0.05318) # based on normal approximation
library(MASS)
confint(m1,level=0.95) # based on profile likelihood

drop1(m1, test = 'Chisq') 


## Consider different link function

m2<-glm(cbind(damage,6-damage)~temp,family=binomial(link=probit))
summary(m2)

m3<-glm(cbind(damage,6-damage)~temp,family=binomial(link=cloglog))
summary(m3)

cbind(fitted(m1),fitted(m2),fitted(m3))
tp<-seq(25,85,0.5)
pl<-ilogit(m1$coef[1]+m1$coef[2]*tp)
pp<-pnorm(m2$coef[1]+m2$coef[2]*tp)
pc<-1-exp(-exp(m3$coef[1]+m3$coef[2]*tp))
plot(tp,pl,type="l",ylab="Probability of Failure",xlab="Temp")
lines(tp,pp,lty="dotted")
lines(tp,pc,lty="longdash")

matplot(tp,cbind(pp/pl,pc/pl),type="l",lty=c("dotted","twodash"),xlab="Temp",ylab="Prob Ratio")

# baby feeding data
data(babyfood)
babyfood
attach(babyfood)
xtabs(disease/(disease+nondisease)~sex+food)
m1<-glm(cbind(disease,nondisease)~sex+food,family=binomial)
summary(m1)

m2<-update(m1,~.+sex:food)
summary(m2)

drop1(m1,test="Chi")

summary(m1)

exp(-0.3126) # girls reduced odds of respiratory disease to 73% of that for the boys

library(MASS)
exp(confint(m1))

# Finding effective "dose"

data(orings)
attach(orings)
m1<-glm(cbind(damage,6-damage)~temp,family=binomial)

a<-predict(m1,newdata=data.frame(temp=30),se=T)
ilogit(c(a$fit-1.96*a$se.fit,a$fit+1.96*a$se.fit)) #CI for failure prob at temp=30


#Q: At what temperature the failure prob of an O-ring is 80%? 
# Effective Dose

that<-(-m1$coef[1]+logit(0.8))/m1$coef[2] # Estimate of ED80
V<-vcov(m1) # covariance matrix of the parameter estimators
dg<-c(-1/m1$coef[2],(m1$coef[1]-logit(0.8))/m1$coef[2]^2) #gradient
se<-sqrt(t(dg)%*%V%*%dg) 
library(MASS)
dose.p(m1,p=c(0.8, 0.5))

# overdispersion

data(troutegg)
attach(troutegg)
troutegg
ftable(xtabs(cbind(survive,total) ~ location+period, troutegg))

m<-glm(cbind(survive,total-survive)~location+period,family=binomial)

summary(m) # goodness of fit not satisfactory

halfnorm(residuals(m)) # Two largest values are labeled. The number can be changed using "nlab=5", for instance. Note also residuals(m,"response") gives response residuals.

plot(location,residuals(m))  

plot(residuals(m)) 

halfnorm(rnorm(20),nlab=1)

halfnorm(rchisq(20,10),nlab=3)

# Is the interaction important?

elogits=log((troutegg$survive+0.5)/(troutegg$total-troutegg$survive+0.5)) # get empirical logits

with(troutegg,interaction.plot(period, location,elogits)) # The with( ) function applys an expression to a dataset.

summary(glm(cbind(survive,total-survive)~location*period,family=binomial))

# The intended conclusison from the textbook is that the ineractions are not important. On the other hand, the AIC of the saturated model is much smaller. So it is hard to say... 

# Note that there are only 20 observations and the model with interaction terms leaves no df

(s2 <- sum(residuals(m,type="pearson")^2)/(20-8)) # estimate dispersion parameter

drop1(m,scale=s2,test="F")

summary(m,dispersion=s2) # re-estimate the se

# when the covariate classes are not roughly equal in size, use other methods
library(dispmod)
dmod <- glm.binomial.disp(m)
summary(dmod)


# Ungrouped data

library(MASS)
library(randomForest)


hosmerlem = function(y, yhat, g=10) {
  cutyhat = cut(yhat,
                breaks = quantile(yhat, probs=seq(0,
                                                  1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq,p.value=P))
}         # Hosmer-Lemeshow test
## Read p147 for Agresti's book for the detail of Hosmer-Lemeshow test

loglikel<-function(a,y,delta){
  n<-length(a)
  l<-0
  a[a>1-delta]<-1-delta
  a[a<delta]<-delta
  for (i in 1:n) {
    l<-l+y[i]*log(a[i])+(1-y[i])*log(1-a[i])
  }
  l
}  # computing predictive log-likelhood, with truncation via delta to avoid 0 or 1 in probability 

delta=0.05

data.g<-function(n,p1,p2,rho,beta0,beta,gamma){
  V<-matrix(0,p2,p2)
  y<-rep(0,n)
  x1<-matrix(rnorm(n*p1),byr=T,ncol=p1)
  for (i in 1:p2) {
    for (j in 1:p2)
      V[i,j]=rho^(abs(i-j))
  }
  x2<-mvrnorm(n, rep(0,p2), V)
  x<-cbind(x1,x2)
  p=p1+p2
  lin<-beta0+x%*%beta+gamma*x[,1]*x[,2]
  prob<-exp(lin)/(1+exp(lin))
  for (k in 1:n) {
    y[k]<-rbinom(1,1,prob[k])
  }
  list(x=x,y=y)
}      # data generation based on a logit model, with gamma controlling the interaction term X1*X2


dat<-data.g(100,4,2,0.5,0.5,c(2,2,1,0,0,1),0) # no interaction in data generation
x<-dat$x
y<-dat$y
dat<-data.frame(y=y,x=x)
summary(glm(y ~ . , data=dat,binomial))

dat<-data.g(100,4,2,0.5,0.5,c(2,2,1,0,0,1),2) # there is an interaction term in the data generation 
x<-dat$x
y<-dat$y
dat<-data.frame(y=y,x=x)
summary(glm(y ~ . , data=dat,binomial))  # the deviance cannot be used to justify the model



goft<-function(x,y,r,m,delta){
  n<-length(y)
  dat<-data.frame(y=y,x=x)	
  mod=glm(y ~ ., data=dat,binomial)
  hl_p=hosmerlem(y, yhat=fitted(mod))$p.value
  sel1<-0
  sel2<-0
  win1=0
  win2=0
  for (j in 1:m) {
    ss<-sample(n,r)
    dat_T<-dat[ss,]
    dat_E<-dat[-ss,]
    m_T<-update(mod,.~.,data=dat_T)
    pred<-predict(m_T,dat_E,type="response")
    l_m<-loglikel(pred,dat_E[,1],delta)
    rf<-randomForest(factor(y) ~ ., data=dat_T, importance=TRUE)
    bg<-randomForest(factor(y) ~ ., mtry=6, data=dat_T, importance=TRUE)
    a<-predict(rf, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE,proximity=FALSE, nodes=FALSE)
    b<-predict(bg, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    l_rf=loglikel(a[,2],dat_E[,1],delta)
    l_bg=loglikel(b[,2],dat_E[,1],delta)
    if (l_m>l_rf) win1=win1+1
    if (l_m>l_bg) win2=win2+1
  }
  win1_f<-win1/m
  win2_f<-win2/m
  if (win1_f>0.5) sel1<-1
  if (win2_f>0.5) sel2<-1
  list(hl_p=hl_p,win_rf=win1_f,win_bg=win2_f,sel_rf=sel1, sel_bg=sel2)
}

dat<-data.g(100,4,2,0.5,0.5,c(2,2,1,0,0,1),0) # no interaction in data generation
x<-dat$x
y<-dat$y
dat<-data.frame(y=y,x=x)

r=50
m=100
delta=0.05

goft(x,y,r,m,delta)

set.seed(1005)

dat<-data.g(100,4,2,0.5,0.5,c(2,2,1,0,0,1),2) #  there is an interaction term in the data generation
x<-dat$x
y<-dat$y
dat<-data.frame(y=y,x=x)

r=50
m=100
delta=0.05

goft(x,y,r,m,delta)

# a data example

data(wcgs, package="faraway")

wcgs$bmi <- with(wcgs, 703*wcgs$weight/(wcgs$height^2))
wcgs$chd<-unclass(wcgs$chd)-1

lmod <- glm(chd ~ age + height + weight +bmi + sdp + dbp + chol + dibep + cigs +arcus, family=binomial, wcgs) 
lmodr <- step(lmod, trace=0) # the selected variables are: age + height +bmi + sdp + chol + dibep + cigs +arcus
sumary(lmodr)

wcgs=wcgs[,-c(11,12)]
wcgs=data.frame(y=wcgs[,10],x=wcgs[,-10])
wcgs=na.omit(wcgs)
attach(wcgs)

n<-length(chd)
r=n/2

goft<-function(m,delta){
  mod=glm(chd ~ age + height +bmi + sdp + chol + dibep + cigs +arcus, family=binomial, wcgs)
  hl_p=hosmerlem(chd, yhat=fitted(mod))$p.value
  sel1<-0
  sel2<-0
  win1=0
  win2=0
  for (j in 1:m) {
    ss<-sample(n,r)
    dat_T<-wcgs[ss,-10]
    dat_E<-wcgs[-ss,-10]
    y_T=wcgs[ss,10]
    y_E=wcgs[-ss,10]
    m_T<-update(mod,.~.,data=wcgs[ss,])
    pred<-predict(m_T,wcgs[-ss,],type="response")
    l_m<-loglikel(pred,y_E,delta)
    rf<-randomForest(factor(y_T) ~ ., data=dat_T, importance=TRUE)
    bg<-randomForest( factor(y_T)~ ., mtry=11, data=dat_T, importance=TRUE)
    a<-predict(rf, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE,proximity=FALSE, nodes=FALSE)
    b<-predict(bg, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    l_rf=loglikel(a[,2],y_E,delta)
    l_bg=loglikel(b[,2],y_E,delta)
    if (l_m>l_rf) win1=win1+1
    if (l_m>l_bg) win2=win2+1
  }
  win1_f<-win1/m
  win2_f<-win2/m
  if (win1_f>0.5) sel1<-1
  if (win2_f>0.5) sel2<-1
  list(hl_p=hl_p,win_rf=win1_f,win_bg=win2_f,sel_rf=sel1, sel_bg=sel2)
}
goft(100,0.05)  # CV strongly prefers the model over random forest or bagging

# check on a poor model

goft<-function(m,delta){
  mod=glm(chd ~ age + height , family=binomial, wcgs) # most variables dropped
  hl_p=hosmerlem(chd, yhat=fitted(mod))$p.value
  sel1<-0
  sel2<-0
  win1=0
  win2=0
  for (j in 1:m) {
    ss<-sample(n,r)
    dat_T<-wcgs[ss,-10]
    dat_E<-wcgs[-ss,-10]
    y_T=wcgs[ss,10]
    y_E=wcgs[-ss,10]
    m_T<-update(mod,.~.,data=wcgs[ss,])
    pred<-predict(m_T,wcgs[-ss,],type="response")
    l_m<-loglikel(pred,y_E,delta)
    rf<-randomForest(factor(y_T) ~ ., data=dat_T, importance=TRUE)
    bg<-randomForest( factor(y_T)~ ., mtry=11, data=dat_T, importance=TRUE)
    a<-predict(rf, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE,proximity=FALSE, nodes=FALSE)
    b<-predict(bg, dat_E, type="prob",norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
    l_rf=loglikel(a[,2],y_E,delta)
    l_bg=loglikel(b[,2],y_E,delta)
    if (l_m>l_rf) win1=win1+1
    if (l_m>l_bg) win2=win2+1
  }
  win1_f<-win1/m
  win2_f<-win2/m
  if (win1_f>0.5) sel1<-1
  if (win2_f>0.5) sel2<-1
  list(hl_p=hl_p,win_rf=win1_f,win_bg=win2_f,sel_rf=sel1, sel_bg=sel2)
}
goft(100,0.05)

# result: The model is not rejected by HL test, but is beaten easily by random forest and bagging.

