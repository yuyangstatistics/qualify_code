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

# compute VSD for model_check using ARM with prior for glmvsd
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
v_BIC  # seems like BIC has closer results
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
(model_check <- as.numeric(betah))

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
order(a$importance) # the 1st is the least important, the last is the most important

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

cv.out=cv.glmnet(x,y,alpha=1,nfold=10)
lambdah=cv.out$lambda.1se
lasso.mod=glmnet(x,y,alpha=1,lambda=lambdah) 
betah=lasso.mod$beta
betah[abs(betah)>0]<-1
(model_check <- as.numeric(betah))

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
b_BIC

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
# generate simulation data
library(glmvsd)
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

