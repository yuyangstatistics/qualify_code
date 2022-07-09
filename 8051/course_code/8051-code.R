# import packages
library(MASS)
library(alr4)
library(leaps)
library(faraway)
library(glmnet) # Read http://cran.r-project.org/web/packages/glmnet/glmnet.pdf
library(ncvreg) # for SCAD and MCP
library(glmvsd)
library(randomForest)
library(SOIL)

##### model-selection-updated.R #####
library(alr4)
data(BGSall)

m<-lm(WT18~Sex+WT2+HT2+WT9+HT9+LG9+ST9+HT18+LG18+ST18+Soma,data=BGSall)
m<-lm(WT18~.,data=BGSall) # same as above
summary(m)
step(m) # if scope is missing, then backward deletion

step(m,scope=list(upper=~.,lower=~1)) # default direction is both forward and backward
m0 <- step(m,scope=list(upper=~.,lower=~1), trace = FALSE) # only output the final selected one
extractAIC(m0)

step(lm(WT18~1,data=BGSall),scope=list(
  upper=~Sex+WT2+HT2+WT9+HT9+LG9+ST9+HT18+LG18+ST18+Soma,lower=~1))
# m is the initial model in the search

step(lm(WT18~1,data=BGSall),scope=list(
  upper=~Sex+WT2+HT2+WT9+HT9+LG9+ST9+HT18+LG18+ST18+Soma,lower=~1),direction="forward")


step(m,k=log(length(BGSall$WT18))) # backward deletion based on BIC instead of AIC (but name unchanged)


library(leaps) # a library for all subset selection

attach(BGSall)
x<-cbind(Sex,WT2,HT2,WT9,HT9,LG9,ST9,HT18,LG18,ST18,Soma)
print(outs<-leaps(x,WT18)) # the default gives 10 best (or max # of) models at each dimension
plot(outs$size, outs$Cp, log = "y", xlab = "p", ylab = expression(C[p]))   # Y axis is to be logarithmic

print(leaps(x,WT18,nbest=1))  # give the best model at each size

reg.all<-regsubsets(x,WT18,nvmax=11) # the default gives the best model for each dimension

reg.summary<-summary(reg.all)

names(reg.summary)

reg.summary$rsq

reg.summary$which

coef(reg.all,id=3) # coefficient of the third model on the list (with 3 predictors plus the intecept)

(indi=reg.summary$which[3,])

fit3=lsfit(x[,indi[-1]],WT18) # fitting a least square model; "-1" removes the intercept

(fit3b <- glm(WT18~x[,indi[-1]])) # glm can give us AIC value directly
fit3b$aic


extractAIC(fit3b,k=2)  #Calculate AIC

npar <- fit3b$df.null - fit3b$df.residual + 1
n <- length(WT18)
extractAIC(fit3b,k=2)+2*npar*(npar+1)/(n-npar-1)  #Calculate AICc

extractAIC(fit3b,k=log(n))  #Calculate BIC


indi=reg.summary$which[2,]

fit2=lsfit(x[,indi[-1]],WT18) # fitting the model; "-1" removes the intercept

fit2b=glm(WT18~x[,indi[-1]]) # glm can give us AIC value directly
fit2b$aic


extractAIC(fit2b,k=2)  #Calculate AIC
extractAIC(fit2b,k=log(n))  #Calculate BIC. It seems that the BIC value by leaps is a shifted value of the true BIC.


par(mfrow=c(2,2))
plot(reg.summary$rsq ,xlab="Model Size",ylab="R^2",
     type="l")
plot(reg.summary$adjr2 ,xlab="Model Size ",
     ylab="Adjusted R^2",type="l")
plot(reg.summary$cp ,xlab="Model Size",ylab="Cp",
     type="l")
plot(reg.summary$bic ,xlab="Model Size ",
     ylab="BIC",type="l")
par(mfrow=c(1,1))

plot(reg.all,scale="bic")  # may choose other criteria, e.g, cp.
?plot.regsubsets  # The plot gives visual understanding on which variables are involved in the models at different levels of the criterion chosen

which.min(reg.summary$bic)
which.max(reg.summary$adjr2)

summary(lm(WT18~x))

reg.fwd=regsubsets(x,WT18,nvmax=11, method ="forward") # forward addition
summary(reg.fwd)
coef(reg.fwd,id=3)

reg.bwd=regsubsets(x,WT18,nvmax=11, method ="backward")  # backward deletion
summary(reg.bwd)  # Note the forward and backward selections do not necessarily give the same results or agree with the best subset models from exhaustive search

detach(BGSall)




# Now try Lasso and ridge regression and compare them

library(glmnet) # Read http://cran.r-project.org/web/packages/glmnet/glmnet.pdf
library(ncvreg) # for SCAD and MCP

library(MASS)

data.generation<-function(n,p,intercept=0,beta,rho,sig){
  x<-matrix(0,n,p)
  Sigma <- matrix(rep(rho,p*p),byr=T,ncol=p) # set all predictors to have the same correlation
  diag(Sigma)<-rep(1,p) # compound symmatry correlation
  x<-mvrnorm(n,rep(0,p),Sigma)
  y<-intercept+x%*%beta+sig*rnorm(n)
  list(x=x,y=y)
}



p<-200
n<-50
rho<-0.5
sig<-1
beta<-c(rep(1,7),rep(0,193)) # only the first 7 are non-zero

dat<-data.generation(n,p,intercept=0,beta,rho,sig)
x=dat$x; y=dat$y

grid=10^seq(10,-3,length=100) # specify lambda values to be considered for 

ridge.mod=glmnet(x,y,alpha=0,lambda=grid)  # alpha = 0 does ridge regression, alpha = 1 is lasso.

ridge.mod=glmnet(x,y,alpha=0)  # defaut choice for lambda      
# plot(ridge.mod,type.coef="2norm")  #If type.coef="2norm" then a single curve per variable, else if type.coef="coef",
#   a coefficient plot per response. The command does not seem to work, with a few warnings.
plot(ridge.mod,xvar="norm") # the xlab should be L2 norm for ridge regression, but due to the package design, all of the 
# xlab is L1 norm. Should pay attention to this.
# Ridge regression does no variable selection, so we could notice the number of variables remains 200 in the plot.
plot(ridge.mod,xvar="lambda",label=TRUE)   # choice for xvar: "norm", "lambda", "dev" (the percent deviance explained)
plot(ridge.mod)
cv.out=cv.glmnet(x,y,alpha=0)  # default is 10 fold. Can change using, e.g., nfold=5
cv.glmnet(x,y,alpha=0,nfold=5)
cv.out$lambda.min   # best lambda, one may also use lambda.1se (largest lambda within 1 sd of the smallest cv value)
cv.out$lambda.1se   # a larger cv would lead to a simpler model
plot(cv.out)    # all the 200 predictors are involved at the different lambda values
ridge.pred=predict(ridge.mod,s=cv.out$lambda.min,newx=x)
predict(ridge.mod,type="coefficients",s=cv.out$lambda.min) # Look at the estimated coefficients
predict(ridge.mod,s=cv.out$lambda.min,newx=x[1,,drop=FALSE]) # predict for case 1

lasso.mod=glmnet(x,y,alpha=1,lambda=grid) 
plot(lasso.mod)    # the numbers on top indicate how many predictos have non-zero coefficients
plot(lasso.mod,xvar="lambda",label=TRUE)
lasso.mod=glmnet(x,y,alpha=1) 
cv.out.l=cv.glmnet(x,y,alpha=1)  
plot(cv.out.l)
cv.out.l$lambda.min  
lasso.pred=predict(lasso.mod,s=cv.out.l$lambda.min,newx=x)
coeff.lasso=predict(lasso.mod,type="coefficients",s=cv.out.l$lambda.min)
which(abs(coeff.lasso)>0)

mcp.mod <- cv.ncvreg(x,y,family='gaussian',penalty='MCP',nfolds=5)
pre.mcp <- predict(mcp.mod,X=x[1,,drop=FALSE],lambda=mcp.mod$lambda.min)
mcp.fit=mcp.mod$fit
coeff.mcp=mcp.fit$beta[,mcp.mod$min]
which(abs(coeff.mcp)>0)

scad.mod <- cv.ncvreg(x,y,family='gaussian',penalty='SCAD',nfolds=5)
pre.scad <- predict(scad.mod,X=x[1,,drop=FALSE],lambda=scad.mod$lambda.min)
scad.fit=scad.mod$fit
coeff.scad=scad.fit$beta[,scad.mod$min]
which(abs(coeff.scad)>0)

cv.compare<-function(x,y,n2,r){
  n<-nrow(x)
  err.r<-rep(0,r)
  err.l<-rep(0,r)
  err.m<-rep(0,r)
  err.s<-rep(0,r)
  for (i in 1:r) {
    train=sample(1:n, n-n2)
    test=(-train)
    y.test=y[test]
    
    ridge.mod=glmnet(x[train,],y[train],alpha=0)
    cv.out=cv.glmnet(x[train,],y[train],alpha=0)  
    ridge.pred=predict(ridge.mod,s=cv.out$lambda.min,newx=x[test,])
    err.r[i]=mean((ridge.pred-y.test)^2)
    
    lasso.mod=glmnet(x[train,],y[train],alpha=1)
    cv.out.l=cv.glmnet(x[train,],y[train],alpha=1)  
    lasso.pred=predict(lasso.mod,s=cv.out.l$lambda.min,newx=x[test,])
    err.l[i]=mean((lasso.pred-y.test)^2)
    
    mcp.mod <- cv.ncvreg(x[train,],y[train],family='gaussian',penalty='MCP',nfolds=5)
    mcp.pred <- predict(mcp.mod,X=x[test,],lambda=mcp.mod$lambda.min)
    err.m[i]=mean((mcp.pred-y.test)^2)
    
    scad.mod <- cv.ncvreg(x[train,],y[train],family='gaussian',penalty='SCAD',nfolds=5)
    scad.pred <- predict(scad.mod,X=x[test,],lambda=scad.mod$lambda.min)
    err.s[i]=mean((scad.pred-y.test)^2)
  }
  par(mfrow=c(2,2))
  plot(err.r,err.l,xlab="Ridge MSPE",ylab="Lasso MSPE")
  abline(0,1)
  plot(err.m,err.l,xlab="MCP MSPE",ylab="Lasso MSPE")
  abline(0,1)
  plot(err.s,err.l,xlab="SCAD MSPE",ylab="Lasso MSPE")
  abline(0,1)
  plot(err.s,err.m,xlab="SCAD MSPE",ylab="MCP MSPE")
  abline(0,1)
  
  
  c(ridge.err=mean(err.r),ridge.sd=sd(err.r),lasso.err=mean(err.l),lasso.sd=sd(err.l), mcp.err=mean(err.m),mcp.sd=sd(err.m), scad.err=mean(err.s),scad.sd=sd(err.s))		
}

cv.compare(x,y,15,50)

beta<-c(1/seq(1:15),rep(0,185))
dat<-data.generation(n,p,intercept=0,beta,rho,sig)
x=dat$x; y=dat$y

cv.compare(x,y,15,50)


#####  Model-selection-Diagnosis.R   #####
# REGRESSION CASE: Low-Dimensional
# generate simulation data
library(MASS)
n <- 50
p <- 8
beta <- c(3,1.5,0,0,2,0,0,0)
sigma <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 1*e
# user provide a model to be checked
model_check <- c(0,1,1,1,0,0,0,1)

library(glmvsd)

glmvsd(x, y, n_train = ceiling(n/2), no_rep = 100,
       n_train_bound = n_train - 2, n_bound = n - 2,
       model_check, psi = 1, family = c("gaussian",
                                        "binomial"), method = c("union", "customize"),
       candidate_models, weight_type = c("BIC", "AIC",
                                         "ARM"), prior = TRUE, reduce_bias = FALSE)

# compute VSD for model_check using ARM with prior4 glmvsd
v_ARM <- glmvsd(x, y, n_train = ceiling(n/2),
                no_rep=50, model_check = model_check, psi=1,
                family = "gaussian", method = "union",
                weight_type = "ARM", prior = TRUE)
# compute VSD for model_check using BIC
v_BIC <- glmvsd(x, y,
                model_check = model_check,
                family = "gaussian", method = "union",
                weight_type = "BIC", prior = TRUE)
# user supplied candidate models
candidate_models = rbind(c(0,0,0,0,0,0,0,1),
                         c(0,1,0,0,0,0,0,1), c(0,1,1,1,0,0,0,1),
                         c(0,1,1,0,0,0,0,1), c(1,1,0,1,1,0,0,0),
                         c(1,1,0,0,1,0,0,0))
v1_BIC <- glmvsd(x, y,
                 model_check = model_check, psi=1,
                 family = "gaussian",
                 method = "customize",
                 candidate_models = candidate_models,
                 weight_type = "BIC", prior = TRUE)

v_ARM
v_BIC
v1_BIC

# REGRESSION CASE: High-Dimensional
# generate simulation data
n <- 50
p <- 100
#p <- 8
beta <- c(3,1.5,0,0,2,0,0,0,rep(0,92))
#beta <- c(3,1.5,0,0,2,0,0,0)
sigma <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 0.4*e
# user provide a model to be checked
# model_check <- c(1,0,1,0,1,0,0,1,rep(0,92))

cv.out=cv.glmnet(x,y,alpha=1,nfold=10)
lambdah=cv.out$lambda.1se
lasso.mod=glmnet(x,y,alpha=1,lambda=lambdah) 
betah=lasso.mod$beta
betah[abs(betah)>0]<-1
model_check <- as.numeric(betah)

# compute VSD for model_check using ARM with prior4 glmvsd
v_ARM <- glmvsd(x, y, n_train = ceiling(n/2),
                no_rep=100, model_check = model_check, psi=1,
                family = "gaussian", method = "union",
                weight_type = "ARM", prior = TRUE)
# compute VSD for model_check using BIC
v_BIC <- glmvsd(x, y,
                model_check = model_check,
                family = "gaussian", method = "union",
                weight_type = "BIC", prior = TRUE)


c(v_ARM$VSD, v_ARM$VSD_minus, v_ARM$VSD_plus) 

c(v_BIC$VSD, v_BIC$VSD_minus, v_BIC$VSD_plus) 

model_check

v_ARM$candidate_models_cleaned

v_ARM$weight


# Measuring Variable Importance

# Random Forest Importance 

library(randomForest)

rf=randomForest(x,y,mtry=8,importance =TRUE)  # mtry: number of variables randomly sampled as candidates at each split.
rf$importance
importance(rf,scale=FALSE) # scale only affects type 1 measure; some researchers suggest scale=FALSE only!
varImpPlot(rf,scale=FALSE)

# SOIL Importance: Based on linear modeling

library(SOIL)

a=SOIL(x, y, n_train = ceiling(n/2), no_rep = 100,
       n_train_bound = ceiling(n/2)- 2, n_bound = n - 2,
       psi = 1, family = "gaussian", method = "union",
       weight_type = "ARM", prior = TRUE)

a$importance
order(a$importance)

# REGRESSION CASE: Moderately High-Dimensional
# generate simulation data
n <- 50
p <- 40
#p <- 8
beta <- c(3,1.5,0,0,2,0,0,0,rep(0,32))
sigma <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 0.9*e

lasso.mod=glmnet(x,y,alpha=1,lambda=lambdah) 
cv.out=cv.glmnet(x,y,alpha=1,nfold=10)
lambdah=cv.out$lambda.1se
betah=lasso.mod$beta
betah[abs(betah)>0]<-1
model_check <- as.numeric(betah)

# compute VSD for model_check using ARM with prior4 glmvsd
v_ARM <- glmvsd(x, y, n_train = ceiling(n/2),
                no_rep=100, model_check = model_check, psi=2,
                family = "gaussian", method = "union",
                weight_type = "ARM", prior = TRUE)
# compute VSD for model_check using BIC
v_BIC <- glmvsd(x, y,
                model_check = model_check,
                family = "gaussian", method = "union",
                weight_type = "BIC", prior = TRUE)


c(v_ARM$VSD, v_ARM$VSD_minus, v_ARM$VSD_plus) 

c(v_BIC$VSD, v_BIC$VSD_minus, v_BIC$VSD_plus) 

model_check

v_ARM$cand

v_ARM$weight


# CLASSIFICATION CASE
# generate simulation data
n = 300
p = 8
b <- c(1,1,1,-3*sqrt(2)/2)
x=matrix(rnorm(n*p, mean=0, sd=1), n, p)
feta=x[, 1:4]%*%b
fprob=exp(feta)/(1+exp(feta))
y=rbinom(n, 1, fprob)
# user provide a model to be checked
model_check <- c(0,1,1,1,0,0,0,1)
# compute VSD for model_check using BIC with prior
b_BIC <- glmvsd(x, y, n_train = ceiling(n/2),
                family = "binomial",
                no_rep=50, model_check = model_check, psi=1,
                method = "union", weight_type = "BIC",
                prior = TRUE)
candidate_models =
  rbind(c(0,0,0,0,0,0,0,1),
        c(0,1,0,0,0,0,0,1),
        c(1,1,1,1,0,0,0,0),
        c(0,1,1,0,0,0,0,1),
        c(1,1,0,1,1,0,0,0),
        c(1,1,0,0,1,0,0,0),
        c(0,0,0,0,0,0,0,0),
        c(1,1,1,1,1,0,0,0))
# compute VSD for model_check using BIC
# user supplied candidate models
b_BIC1 <- glmvsd(x, y,
                 family = "binomial",
                 model_check = model_check, psi=1,
                 method = "customize",
                 candidate_models = candidate_models,
                 weight_type = "BIC")

b_BIC1

# Instability
library(MASS)
# generate simulation data
n <- 50
p <- 8
beta<-c(2.5,1.5,0.5,rep(0,5))
sigma<-matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 0.01*e
ins_seq <- stability.test(x, y, method = "seq",
                          penalty = "SCAD", nrep = 100,
                          remove = 0.1, tau = 0.2, nfolds = 5)

ins_seq

n <- 50
p <- 8
beta<-c(2.5,1.5,0.5,rep(0,5))
sigma<-matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 0.1*e
ins_bs <- stability.test(x, y, method = "bs",
                         penalty = "SCAD", nrep = 100,
                         remove = 0.1, tau = 0.2, nfolds = 5)

ins_bs

n <- 50
p <- 8
beta<-c(2.5,1.5,0.5,rep(0,5))
sigma<-matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 0.1*e
ins_per <- stability.test(x, y, method = "per",
                          penalty = "SCAD", nrep = 100,
                          remove = 0.1, tau = 0.1, nfolds = 5)

ins_per

# High-dimensional case
n <- 50
p <- 100
beta <- c(3,1.5,0,0,2,0,0,0,rep(0,92))
sigma <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p) sigma[i,j] <- 0.5^abs(i-j)
}
x <- mvrnorm(n, rep(0,p), sigma)
e <- rnorm(n)
y <- x %*% beta + 1*e
ins_seq <- stability.test(x, y, method = "seq",
                          penalty = "LASSO", nrep = 100,
                          remove = 0.04, tau = 0.2, nfolds = 5)
ins_seq


#####  nonlinear-reg.R  #####
library(alr4)
data(turk0)
n1 <- nls(Gain ~ th1 + th2*(1 - exp(-th3 * A)), data=turk0,
          start=list(th1=620, th2=200, th3=10))
summary(n1)

plot(Gain ~ A, turk0, xlab="Amount (percent of diet)", ylab="Weight gain, g")
x <- (0:44)/100
lines(x, predict(n1, data.frame(A=x)))

m1 <- lm(Gain ~ as.factor(A), turk0) # the one-way ANOVA model has no bias
anova(n1,m1)    # lack-of-fit test

n2 <- nls(Gain ~ th1 + 150*(1 - exp(-th3 * A)), data=turk0,
          start=list(th1=620, th3=10))

anova(n2,n1)    # compare the two nested models

# May need to create indicator variables to deal with factors

data(turkey)

temp <- model.matrix( ~ as.factor(S) -1, turkey) # "-1" drops intercept
colnames(temp) <- c("S1", "S2", "S3")
turkey <- cbind(turkey, temp)  # the modified data matrix with the indicators

m4 <- nls( Gain ~ (th1 + th2*(1-exp(-th3*A))),
           data=turkey,start=list(th1=620, th2=200, th3=10),weights=m)
# most general
m1 <- nls( Gain ~ S1*(th11 + th21*(1-exp(-th31*A)))+
             S2*(th12 + th22*(1-exp(-th32*A)))+
             S3*(th13 + th23*(1-exp(-th33*A))),data=turkey,
           start= list(th11=620, th12=620, th13=620,
                       th21=200, th22=200, th23=200,th31=10, th32=10, th33=10),weights=m)
# common intercept
m2 <- nls(Gain ~ th1 + S1*(th21*(1-exp(-th31*A)))+
            S2*(th22*(1-exp(-th32*A)))+S3*(th23*(1-exp(-th33*A))),
          data=turkey,start= list(th1=620,th21=200, th22=200, th23=200,
                                  th31=10, th32=10, th33=10),weights=m)
# common intercept and asymptote
m3 <- nls( Gain ~ (th1 + th2 *(S1*(1-exp(-th31*A))+S2*(1-exp(-th32*A))+
                                 S3*(1-exp(-th33*A)))),data=turkey,start= list(th1=620, th2=200,
                                                                               th31=10,th32=10,th33=10),weights=m)

anova(m4, m2 ,m1)
anova(m4,m1)

# boostrap
library(car)
data(segreg)
s1 <- nls(C ~ th0 + th1 * (pmax(0, Temp - gamma)), data=segreg,
          start=list(th0=70, th1=.5, gamma=40))
summary(s1)

s1.boot <- Boot(s1, R=999)
summary(s1.boot)
confint(s1.boot,type="bca")  # The intervals calculated using the adjusted bootstrap percentile (BCa) method. Other choices are "norm", "basic", "stud", and "perc".




#####  binomial-regression.R  #####
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

n<-length(wcgs$chd)
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


#####  Poisson-regression.R  #####

library(faraway)

# gala data (# of species of plants on each of 30 islands)

data(gala)
gala
help(gala)
attach(gala)
gala <- gala[, -2]

m<-glm(Species~Area+Elevation+Nearest+Scruz+Adjacent,family=poisson)
summary(m)

halfnorm(residuals(m))

par(mfrow=c(2,3))
plot(Area,residuals(m))
plot(Elevation,residuals(m))
plot(Nearest,residuals(m))
plot(Scruz,residuals(m))
plot(Adjacent,residuals(m))  # May consider dropping one case.
par(mfrow = c(1, 1))

(summary(m1<-glm(Species~(Area+Elevation+Nearest+Scruz+Adjacent)^2,family=poisson)))
# Interaction terms do not improve too much (the residual deviance still seems to be large). AIC value is much reduced.

(dp<-sum(residuals(m,type="pearson")^2)/m$df.res)
(dp1<-sum(residuals(m1,type="pearson")^2)/m1$df.res)


summary(m,dispersion=dp)
summary(glm(Species ~ ., family = quasipoisson, gala))

drop1(m,test="F")  # preferred
drop1(m1,test="F") 
anova(m,m1,test="F") # It does not give you the desired F test!
(f<-(deviance(m)-deviance(m1))/((m$df.res-m1$df.res)*dp1))
1-pf(f,(m$df.res-m1$df.res),m1$df.res) # The simpler model may be preferred.

step(m1)
step(m1,k=log(30)) # Both AIC and BIC prefer m1, which is different from F test.



## Rate Models

data(dicentric)
help(dicentric)
dicentric
attach(dicentric)
xtabs(ca/cells~doseamt+doserate,dicentric)
interaction.plot(doseamt,doserate,ca/cells)
doseamt<-factor(dicentric$doseamt)
m1<-glm(ca~cells+log(doserate)*doseamt,family=poisson) # the use of cells is problematic
summary(m1)
m2<-glm(ca~log(cells)+log(doserate)*doseamt,family=poisson)
summary(m2)

m3<-glm(ca~offset(log(cells))+log(doserate)*doseamt,family=poisson)
summary(m3)

## Negative Binomial

data(solder)
help(solder)
solder
attach(solder)

m1<-glm(skips~.,family=poisson,solder)
summary(m1)

m2<-glm(skips~(Opening+Solder+Mask+PadType+Panel)^2,family=poisson,solder)
summary(m2)
pchisq(deviance(m2),df.residual(m2),lower=FALSE)

library(MASS)
m3<-glm(skips~.,family=negative.binomial(1),solder) # k=1
summary(m3)

m4<-glm.nb(skips~.,solder) # k being estimated antomatically
summary(m4)
pchisq(deviance(m4),df.residual(m4),lower=FALSE)

#####  contingency-table.R  #####
library(faraway)

# wafer data

y<-c(320,14,80,36)
particle<-gl(2,1,4,labels=c("no","yes")) # gl generates factor levels; the second argument indicates number of repeats
quality<-gl(2,2,4,labels=c("good","bad"))
particle
quality
wafer<-data.frame(y,particle,quality)
wafer

80/320    # odds of "bad" without particle
36/14     # odds of "bad" with particle
(36/14)/(80/320)  # odds ratio

# binomial analysis (400 wafers with no particles, 50 with)


bd<-matrix(y,nrow=2)  # the default is byr=F.
bd

fisher.test(bd) # odds ratio is y11*y22/(y12*y21)

fisher.test(bd,alternative = "greater")


# effects of smoking

data(femsmoke)
femsmoke

ct<-xtabs(y~smoker+dead,femsmoke)
ct
(a <- prop.table(ct,1))

summary(ct)  # gives chisquare test of independence
fisher.test(ct)  

prop.table(xtabs(y~smoker+age,femsmoke),2) # column sums to 1
xtabs(y~smoker+age,femsmoke)
prop.table(xtabs(y~smoker+age,femsmoke),1) # row sums to 1
xtabs(y~smoker+age,femsmoke)

ct1<-xtabs(y~smoker+dead,femsmoke,subset=(age=="18-24"))
ct1
prop.table(ct1,1)

ct2<-xtabs(y~smoker+dead,femsmoke,subset=(age=="25-34"))
ct2
prop.table(ct2,1)

ct3<-xtabs(y~smoker+dead+age,femsmoke)
ct3
apply(ct3,3,function(x) (x[1,1]*x[2,2]/(x[1,2]*x[2,1]))) # observed odds ratio at the different age groups

ct4<-xtabs(y~smoker+dead,femsmoke,subset=(age=="55-64"))
ct4
prop.table(ct4,1)

m<-glm(ct4~factor(c(0,1)),family=binomial)
summary(m)
drop1(m)
summary(ct4) # Pearson chisquare gives similar result on independence 

mantelhaen.test(ct3,exact=TRUE) # test conditional independence in each stratum 
# given age, smoker is correlated with dead

summary(ct3) # mutual independence (not pairwise independence)
# mutual independence: P(ABC) = P(A)P(B)P(C)

modm<-glm(y~smoker+dead+age,femsmoke,family=poisson) # mutual independence model
summary(modm)

modj<-glm(y~smoker*dead+age,femsmoke,family=poisson) # joint independence model
summary(modj)
# joint independece means smoker*dead as a group is independent with age

modc<-glm(y~smoker*age+age*dead,femsmoke,family=poisson) # conditional independence
summary(modc)
# conditional independence means given age, smoker has no interaction with dead,
# namely, given age, the other two factors are independent.

modu<-glm(y~(smoker+age+dead)^2,femsmoke,family=poisson) # uniform association
summary(modu)
# for conditional independence, given age, the odds ratio are all 1; 
# for uniform association, given age, the odds ratio are all the same, but not 
# necessarily 1. In such kind of model, it considers all interaction except the 
# three-factor interaction. With three-factor interaction, the model will become
# the saturated model.

anova(modc,modu,test="Chisq") # modc is not preferred

ctf<-xtabs(fitted(modu)~smoker+dead+age,femsmoke)
apply(ctf,3,function(x) (x[1,1]*x[2,2]/(x[1,2]*x[2,1])))

modsat<-glm(y~smoker*dead*age,femsmoke,family=poisson)
summary(modsat)
step(modsat) # Use AIC to do stepwise selection

# binomial regression when one variable is treated as the response
femsmoke 

ybin<-matrix(femsmoke$y,ncol=2) # create binomial data treating life status as the response

summary(modbin<-glm(ybin~smoker*age,femsmoke[1:14,],family=binomial))
drop1(modbin,test="Chi")
modbinr<-glm(ybin~smoker+age,femsmoke[1:14,],family=binomial) # Note this model corresponds to the Poisson regression with uniform association

c(deviance(modu),deviance(modbinr))

summary(modbina<-glm(ybin~age,femsmoke[1:14,],family=binomial)) # Note this model corresponds to the Poisson regression with conditional independence (modc)
c(deviance(modc),deviance(modbina))

summary(modbins<-glm(ybin~smoker,femsmoke[1:14,],family=binomial)) # Note this model does not correspond to the Poisson regression with joint independence (modj)
c(deviance(modj),deviance(modbins))


summary(modbinn<-glm(ybin~1,femsmoke[1:14,],family=binomial)) # Note this model does not corresponds to the Poisson regression with mutual independence (modc)
c(deviance(modm),deviance(modbinn))  # not equal


# an example of ordinal variables
data(nes96)
help(nes96)
nes96
xtabs(~PID+educ,nes96)
partyed<-as.data.frame.table(xtabs(~PID+educ,nes96))
partyed

modno<-glm(Freq~PID+educ,partyed,family=poisson) # a nominal model
summary(modno)
pchisq(40.743,36,lower=FALSE)  

(partyed$oPID<-unclass(partyed$PID))
(partyed$oeduc<-unclass(partyed$educ))

modod<-glm(Freq~PID+educ+I(oPID*oeduc),partyed,family=poisson) # linear-by-linear association model
summary(modod)

anova(modno,modod,test="Chi")

modod2<-glm(Freq~PID+educ+oPID:educ,partyed,family=poisson) # an ordinal-by-nominal model
summary(modod2)
anova(modod,modod2,test="Chi") # linear-by-linear association model is preferred
#####  simpson-paradox.R   #####
y<-c(360,30,310,390)
n<-c(900,100,450,600)
tr<-c(1,2,1,2) # 1 = 3 mammograms; 2 =  1 mammogram
h<-c(1,1,2,2)  # 1 = With family history; 2 = Without family history
demo<-cbind(Survival=y,Total=n,Treatment=tr,History=h)
xtabs(y~tr,demo)/xtabs(n~tr,demo)
xtabs(y/n~tr+h,demo)
m<-glm(cbind(y,n-y)~tr,family=binomial)
summary(m)


m1<-glm(cbind(y,n-y)~tr+h,family=binomial)
summary(m1) # Missing a lurking variable gives a totally wrong conclusion! The lack of fit is seen from summary(m). But if the history variable is not in the data, no lack of fit can be detected. The new data: y<-c(670,420); n<-c(1350,700); tr<-c(1,2); mm<-glm(cbind(y,n-y)~tr,family=binomial); summary(mm)

mm<-glm(cbind(y,n-y)~tr+h,family=binomial(link=probit))
summary(mm) 


# How did the paradox occur?

a<-xtabs(n~h+tr,demo)

a/(a%*%matrix(1,2,2))     # Reason: Most women with family history went with 3 mammograms; the problme is avoided if women are randomly assigned with 3 mammograms or 1 mammogram. The commond is the same as:  
prop.table(a,1)

m2<-glm(cbind(y,n-y)~tr*h,family=binomial) # no df left.
summary(m2)

# We re-organize the above data and run Poisson regressions

z<-c(360,30,310,390,540,70,140,210)
tr<-c(1,2,1,2,1,2,1,2)
h<-c(1,1,2,2,1,1,2,2)
s<-c(1,1,1,1,2,2,2,2) # 1 = survival; 2 = death

g<-glm(z~s*tr+h*tr,family=poisson)
summary(g) 
summary(m) # Note that the deviances of m and g are the same, though the null deviances are different. Also, the coefficient for s:tr in g is the same as that for tr in m, up to sign change. If we use other link than logit, the deviances are not the same (for cbind(y,n-y)~tr, the deviance does not depend on the choice of the link).

g1<-glm(z~s*tr+s*h+tr*h,family=poisson)
summary(g1)
summary(m1)

g2<-glm(z~s*tr*h,family=poisson)
summary(g2)
summary(m2)


summary(glm(cbind(y,n-y)~1,family=binomial))
summary(glm(z~s+tr*h,family=poisson))




#####  tree.R  #####
# Regression trees

library(MASS)
library(tree)
set.seed(1)       # to make the result repeatable, may fix the seed

train = sample(1:nrow(Boston), nrow(Boston)/2)  # generate the training and test data

tree.boston=tree(medv~., Boston, subset=train)
summary(tree.boston)  # residual mean deviance is just RSS/error df, error df = nrow(Boston)/2-8

plot(tree.boston)
text(tree.boston,pretty=0) # If pretty = 0 then the level names of a factor split attributes are used unchanged. If pretty = NULL, the levels are presented by a, b, . . . z, 0 . . . 5. If pretty is a positive integer, 	abbreviate is applied to the labels with that value for its argument minlength.

cv.boston=cv.tree(tree.boston)  # CV for error estimation in number of nodes, default = 10-fold

plot(cv.boston$size, cv.boston$dev, type='b') # "b" for both "p" (point) and "l" (line)

cv2.boston=cv.tree(tree.boston, K=5)  # 5-fold

plot(cv2.boston$size, cv2.boston$dev, type='b')

prune.boston=prune.tree(tree.boston, best=5)
plot(prune.boston)
text(prune.boston, pretty=0) 

yhat=predict(tree.boston, newdata=Boston[-train, ])
boston.test=Boston[-train, "medv"]
plot(yhat,boston.test)
abline(0,1)
mean((yhat-boston.test)^2)

# classification tree

library(tree)
library(ISLR)    # ISLR: Data for An Introduction to Statistical Learning with Applications in R, by Gareth James, Daniela Witten, Trevor Hastie and Robert Tibshirani

attach(Carseats)
High=ifelse(Sales<=8,"No","Yes")

Carseats =data.frame(Carseats, High)

tree.carseats =tree(High~.-Sales, Carseats )
summary(tree.carseats)  
# the residual mean deviance is -2*maximized binomial loglikelihood diveided by n - tree size 

plot(tree.carseats)
text(tree.carseats,pretty=0)

tree.carseats

set.seed(2)
train=sample(1:nrow(Carseats), 200)
Carseats.test=Carseats[-train ,]
High.test=High[-train]
tree.carseats=tree(High~.-Sales,Carseats, subset=train)
tree.pred=predict(tree.carseats,Carseats.test,type="class")
table(tree.pred,High.test)
(30+27) /200 # error rate

set.seed (3)
cv.carseats =cv.tree(tree.carseats,FUN=prune.misclass) # FUN=prune.misclass uses classification error rate to guide the cross-validation and pruning process, rather than the default for the cv.tree() function, which is deviance

names(cv.carseats)
cv.carseats  # k controls fitness and complexity; deviance is actually the number of errors; size 1 (no splitting at all)

par(mfrow=c(1,2))
plot(cv.carseats$size,cv.carseats$dev,type="b")  
plot(cv.carseats$k,cv.carseats$dev,type="b")

par(mfrow=c(1,1))
prune.carseats=prune.misclass(tree.carseats,best=9)  # prune.misclass is an abbreviation for prune.tree(method = "misclass") for use with cv.tree. Besides "misclass", the other method is "deviance". For regression trees, only the default, deviance, is accepted. For classification trees, the default is deviance and the alternative is misclass (number of misclassifications or total loss).
plot(prune.carseats )
text(prune.carseats,pretty=0)

tree.pred=predict(prune.carseats,Carseats.test,type="class")
table(tree.pred ,High.test)

(24+22) /200 

prune.carseats=prune.misclass(tree.carseats,best=15)
plot(prune.carseats )
text(prune.carseats,pretty=0)
tree.pred=predict(prune.carseats,Carseats.test,type="class")
table(tree.pred ,High.test)
(30+22) /200    # the error rate is higher

# random forests and bagging

library(randomForest)
set.seed(1)

train = sample(1:nrow(Boston), nrow(Boston)/2)  # generate the training and test data
bag.boston=randomForest(medv~.,data=Boston,subset=train,mtry=13,importance =TRUE)    # there are 13 predictors; so mtry=13 does bagging. Default: ntree=500

bag.boston

yhat.bag = predict(bag.boston,newdata=Boston[-train ,])
plot(yhat.bag, boston.test)
abline(0,1)
mean((yhat.bag-boston.test)^2)

bag.boston=randomForest(medv~.,data=Boston,subset=train, mtry=13,ntree=25) # number of trees to grow is 25
yhat.bag = predict(bag.boston ,newdata=Boston[-train ,])
mean((yhat.bag-boston.test)^2)

set.seed(1)
rf.boston=randomForest(medv~.,data=Boston,subset=train,mtry=6,importance =TRUE)  # mtry: number of variables randomly sampled as candidates at each split.
yhat.rf = predict(rf.boston ,newdata=Boston[-train ,])
mean((yhat.rf-boston.test)^2)

importance(rf.boston)	# a variable importance measure. type: 1=mean decrease in accuracy (MSE for regreesion, error rate for classification) based on permutation of each variable values using OOB observations for each tree, averaged over all trees; 2=mean total decrease in node impurity from splitting on the variable (for classification, the node impurity is measured by the Gini index; for regression, it is measured by residual sum of squares)
# scale option:	For permutation based measures, should the measures be divided by their “standard errors”?  

importance(rf.boston,scale=FALSE) # scale only affects type 1 measure; some researchers suggest scale=FALSE only!

varImpPlot(rf.boston,scale=FALSE)

importance(rf.boston,scale=FALSE)

lstat2<-Boston$lstat+0.5*rnorm(nrow(Boston))

Boston2 =data.frame(Boston, lstat2)
rf.boston2=randomForest(medv~.,data=Boston2,subset=train,mtry=14,importance =TRUE)
importance(rf.boston2,scale=FALSE)

dis2<-Boston$dis+0.5*rnorm(nrow(Boston))

Boston3 =data.frame(Boston, dis2)
rf.boston3=randomForest(medv~.,data=Boston3,subset=train,mtry=14,importance =TRUE)

importance(rf.boston3,scale=FALSE)

lstat2<-Boston$lstat+0.5*rnorm(nrow(Boston))

Boston2 =data.frame(Boston, lstat2)
rf.boston2=randomForest(medv~.,data=Boston2,subset=train,mtry=7,importance =TRUE)

importance(rf.boston2,scale=FALSE)

#####  nonparam-regression.R  #####
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
par(mfrow=c(1,1))

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





#####  generalized-additive-modelling.R  #####
# two dimensional kernel estimation

library(faraway)
library(sm)

data(savings)
help(savings)

y<-savings$sr
xx<-cbind(savings$pop15,savings$ddpi)
sm.regression(xx,y,h=c(1,1),xlab="pop15",ylab="growth",zlab="saving rate")
sm.regression(xx,y,h=c(5,5),xlab="pop15",ylab="growth",zlab="saving rate")

# Additive modeling

data(ozone)
help(ozone)
pairs(ozone)

olm<-lm(O3~temp+ibh+ibt,ozone)
summary(olm)
plot(olm)   # it shows the linear model does not fit well

library(gam)
help(gam)

am1<-gam(O3~lo(temp)+lo(ibh)+lo(ibt),data=ozone) # lo: loess smoother
summary(am1)
par(mfrow=c(1,3))
plot(am1,residuals=TRUE,se=TRUE,pch=".")

am2<-gam(O3~temp+lo(ibh)+lo(ibt),data=ozone) # linear in temp 
summary(am2)
par(mfrow=c(1,3))
plot(am2,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am2))
plot(ozone$ibt,residuals(am2))
plot(ozone$ibh,residuals(am2))

anova(am2,am1,test="F")  # the simpler model is rejected

am3<-gam(O3~lo(vh)+lo(wind)+lo(humidity)+lo(ibh)+lo(ibt)+lo(dpg)+lo(vis)+lo(doy),data=ozone)
summary(am3)

am4<-gam(O3~vh+lo(humidity)+lo(ibh)+lo(ibt)+lo(dpg)+lo(vis)+lo(doy),data=ozone) 
summary(am4)

anova(am1,am4,test="F")
par(mfrow=c(3,3))
plot(am4,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am4))
plot(ozone$ibt,residuals(am4))
plot(ozone$ibh,residuals(am4))

# we may also try natural splines to fit additive models
attach(ozone)
am5<-glm(O3~ns(temp,2)+ibt+ns(ibh))
summary(am5)
plot(am5,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am4))
plot(ozone$ibt,residuals(am4))
plot(ozone$ibh,residuals(am4))








