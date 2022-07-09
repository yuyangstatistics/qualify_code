###### This script saves functions used in 2019Q2.Rmd ######


##### ========== Data Exploration ============ #####

# plot the histograms, qqplots, and boxplots for the covariates
hist.qq.box.plot <- function(data) {
  # plot the histograms and qqplots for the covariates
  
  # assumes the name of the response is "y"
  X <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  # histogram
  par(mfrow = c(2, 3), mar = c(3, 3, 3, 3))
  for (i in 1:ncol(X)) {
    hist(X[, i], main = bquote(Histogram ~ of ~ V * .(i)), 
         xlab = bquote(V * .(i)))
  }
  
  # qqplot
  for (i in 1:ncol(X)) {
    qqnorm(X[, i], main = bquote(QQ ~ Plot ~ of ~ V * .(i)))
    qqline(X[, i])
  }
  
  # boxplot
  for (i in 1:ncol(X)) {
    boxplot(X[, i], main = bquote(Boxplot ~ of ~ V * .(i)))
  }
  par(mfrow = c(1, 1))
  
  # plot y
  par(mfrow = c(2, 2))
  hist(y, main = "Histogram of y")
  qqnorm(y, main = "QQplot of y")
  qqline(y)
  boxplot(y, main = "Boxplot of y")
  par(mfrow = c(1, 1))
}


# check outliers
outlier.check <- function(lmod) {
  par(mfrow = c(2, 2))
  plot(lmod, which = 1)
  plot(lmod, which = 3)
  halfnorm(cooks.distance(lmod))
  plot(lmod, which = 4)
  par(mfrow = c(1, 1))
  print(outlierTest(lmod))
}


# plot scatter plot using the training data
scatter.plot <- function(data) {
  # Note: since scatter plot requires labels, so it should use training set.
  
  # get X and y
  X <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  # check scatter plot
  par(mfrow = c(2, 3))
  for (i in 1:ncol(X)) {
    plot(X[, i], y, main = bquote(Scatter ~ Plot ~ of ~ V * .(i)), 
         xlab = bquote(V * .(i)), ylab = "y")
    # fit a linear line
    abline(lm(y ~ X[, i]))
    # fit a quadratic curve
    lines(sort(X[, i]), fitted(lm(y ~ X[, i] + I(X[, i]^2)))[order(X[, i])], 
          col='red', lwd = 2)
  }
  par(mfrow = c(1, 1))
  
}


##### ========== Print confusion matrices ============ #####

print.confmat <- function(preds, actuals, title=""){
  print(paste(title, "Confusion Matrix:"))
  print(table(preds, actuals))
}

confmat.abic <- function(data){
  # assumes the name of the response column is "y"
  
  n = nrow(data)
  y = data$y
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # selection using AIC, both
  lmod <- glm(y ~ ., family = binomial, data[train, ]) # fit a full model first
  aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                  direction = 'both', trace=FALSE)
  aic.probpred <- predict(aic.mod, newdata = data[test, ], type = "response")
  aic.pred <- as.numeric(aic.probpred > 0.5)
  print.confmat(aic.pred, y.test, "AIC")
  
  # selection using BIC, both
  bic.mod <- step(lmod, k = log(nrow(data)),
                  scope = list(upper=~.,lower=~1),
                  direction = 'both', trace=FALSE)
  bic.probpred <- predict(bic.mod, newdata = data[test, ], type = "response")
  bic.pred <- as.numeric(bic.probpred > 0.5)
  print.confmat(bic.pred, y.test, "BIC")

}

confmat.regularized <- function(data){
  # print out the confusion matrix
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # selection using lasso
  lasso.mod <- glmnet(x[train,], y[train], family = "binomial", alpha=1)
  cv.out.l <- cv.glmnet(x[train,], y[train], family = "binomial", alpha=1)  
  lasso.probpred <- drop(predict(lasso.mod, s = cv.out.l$lambda.1se, 
                                 newx = x[test, ], type = "response"))
  lasso.pred <- as.numeric(lasso.probpred > 0.5)
  print.confmat(lasso.pred, y.test, "LASSO")
  
  # selection using mcp
  mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'binomial',
                       penalty = 'MCP', nfolds = 5)
  mcp.probpred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min, 
                          type = "response")
  mcp.pred <- as.numeric(mcp.probpred > 0.5)
  print.confmat(mcp.pred, y.test, 'MCP')
  
  # selection using scad
  scad.mod <- cv.ncvreg(x[train,], y[train], family = 'binomial',
                        penalty = 'SCAD', nfolds = 5)
  scad.probpred <- predict(scad.mod, X=x[test,], lambda = scad.mod$lambda.min, 
                           type = "response")
  scad.pred <- as.numeric(scad.probpred > 0.5)
  print.confmat(scad.pred, y.test, "SCAD")
  
}

confmat.xgb.rf <- function(data, xgb.params, xgb.nrounds, 
                        seed=99){
  
  require(ggplot2)
  require(ggpubr)
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  
  set.seed(seed) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  X.test <- xgb.DMatrix(x[test, ])
  
  # xgboost model fit
  dtrain <- xgb.DMatrix(data = x[train,], label = y[train])
  xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
  xgb.pred <- predict(xgb.mod, X.test)
  xgb.plot <- plot.pred(xgb.pred, y.test, "XGB")
  err.xgb <- mean((xgb.pred - y.test)^2)
  
  # random forest model fit
  rf.mod <- randomForest(y ~ ., data = data[train, ])
  rf.pred <- predict(rf.mod, newdata = data[test, ])
  rf.plot <- plot.pred(rf.pred, y.test, "RF")
  err.rf <- mean((rf.pred - y.test)^2)
  
  print(c(xgb.err = err.xgb, rf.err=err.rf))
  
  figure <- ggarrange(xgb.plot, rf.plot, 
                      labels = c('A', 'B'), 
                      ncol = 2, nrow = 1)
  figure
  
}



confmat.glm <- function(data, glm.form){
  # assumes the name of the response column is "y"
  n = nrow(data)
  y = data$y
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # fit a glm model
  glm.mod <- glm(glm.form, family=binomial, data[train, ])
  glm.probpred <- predict(glm.mod, newdata = data[test, ], type = "response")
  glm.pred <- as.numeric(glm.probpred > 0.5)
  
  print.confmat(glm.pred, y.test, "GLM")
}


##### ========== Evaluation Metrics ============ #####
accuracy_score <- function(preds, actuals) {
  # preds here should be binary: 0/1
  mean(preds == actuals)
}

error_rate <- function(preds, actuals) {
  # preds here should be binary: 0/1
  mean(preds != actuals)
}

auc_score <- function(probpreds, actuals) {
  # probpreds here should be estimated probabilities
  require(mltools)
  auc_roc(probpreds, actuals)
}

precision_score <- function(preds, actuals) {
  # preds here should be binary: 0/1
  # precision = TP / (TP + FP)
  t <- table(preds, actuals)
  p <- t[2, 2] / (t[2, 1] + t[2, 2])
  p
}

recall_score <- function(preds, actuals) {
  # preds here should be binary: 0/1
  # recall/sensitivity = TP / (TP + FN)
  t <- table(preds, actuals)
  r <- t[2, 2] / (t[1, 2] + t[2, 2])
  r
}

f1_score <- function(preds, actuals) {
  # preds here should be binary: 0/1
  # recall/sensitivity = TP / (TP + FN)
  t <- table(preds, actuals)
  p <- t[2, 2] / (t[2, 1] + t[2, 2])
  r <- t[2, 2] / (t[1, 2] + t[2, 2])
  f1 <- (2 * p * r) / (p + r)
  f1
}

##### ========== Cross Validation ============ #####


paste.mean.sd <- function(a) {
  # paste mean and standard deviation for better visualization
  paste(round(mean(a), 4), "(", round(sd(a), 4), ")", sep = "")
}

# compare regularized model selection methods
# can be used for data with or without interactions, just change data input
cv.regularized.compare <- function(data,n2,r){
  # use repeated random sub-sampling cross validation to compare 
  # regularized model selection methods
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  
  result <- foreach (i = 1:r, .combine = rbind) %dopar% {
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    y.test <- y[test]
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], family = "binomial", alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], family = "binomial", alpha=1)  
    lasso.probpred <- drop(predict(lasso.mod, s = cv.out.l$lambda.1se, 
                                   newx = x[test, ], type = "response"))
    lasso.pred <- as.numeric(lasso.probpred > 0.5)
    acc.l <- accuracy_score(lasso.pred, y.test)
    auc.l <- auc_score(lasso.probpred, y.test)
    
    # selection using mcp
    mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'binomial',
                         penalty = 'MCP', nfolds = 5)
    mcp.probpred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min, 
                        type = "response")
    mcp.pred <- as.numeric(mcp.probpred > 0.5)
    acc.m <- accuracy_score(mcp.pred, y.test)
    auc.m <- auc_score(mcp.probpred, y.test)
    
    # selection using scad
    scad.mod <- cv.ncvreg(x[train,], y[train], family = 'binomial',
                         penalty = 'SCAD', nfolds = 5)
    scad.probpred <- predict(scad.mod, X=x[test,], lambda = scad.mod$lambda.min, 
                            type = "response")
    scad.pred <- as.numeric(scad.probpred > 0.5)
    acc.s <- accuracy_score(scad.pred, y.test)
    auc.s <- auc_score(scad.probpred, y.test)
    
    c(acc.l, acc.m, acc.s, auc.l, auc.m, auc.s)
  }
  
  acc.l <- result[, 1]
  acc.m <- result[, 2]
  acc.s <- result[, 3]
  auc.l <- result[, 4]
  auc.m <- result[, 5]
  auc.s <- result[, 6]
  
  list(
    Accuracy = c(lasso.acc = paste.mean.sd(acc.l), 
                 mcp.acc = paste.mean.sd(acc.m),
                 scad.acc = paste.mean.sd(acc.s)),
    AUC = c(lasso.auc = paste.mean.sd(auc.l), 
                 mcp.auc = paste.mean.sd(auc.m),
                 scad.auc = paste.mean.sd(auc.s))
  )

}



# compare AIC, BIC model selection methods, including forward, backward and both
# the seeds are the same as the ones in cv.regularized.compare
# used for data without interactions
cv.abic.compare<-function(data,n2,r){
  # use cross validation to compare aic and bic model selection methods
  
  # the seeds are the same as those in cv.compare.regularized, so the
  # results from these two functions are comparable.
  
  # the reason I split the comparison into two functions is that for
  # some dataset, AIC and BIC cannot work due to large dimension size
  
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n = nrow(data)
  y = data$y
  result <- foreach (i = 1:r, .combine = rbind) %dopar% {
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    y.test <- y[test]
    
    # selection using AIC, both
    lmod <- glm(y ~ ., family = binomial, data[train, ]) # fit a full model first
    aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    aic.probpred <- predict(aic.mod, newdata = data[test, ], type = "response")
    aic.pred <- as.numeric(aic.probpred > 0.5)
    auc.aic <- auc_score(aic.probpred, y.test)
    acc.aic <- accuracy_score(aic.pred, y.test)
    
    # selection using BIC, both
    bic.mod <- step(lmod, k = log(nrow(data)),
                    scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    bic.probpred <- predict(bic.mod, newdata = data[test, ], type = "response")
    bic.pred <- as.numeric(bic.probpred > 0.5)
    auc.bic <- auc_score(bic.probpred, y.test)
    acc.bic <- accuracy_score(bic.pred, y.test)
    
    c(auc.aic, acc.aic, auc.bic, acc.bic)
  }
  
  auc.aic <- result[, 1]
  acc.aic <- result[, 2]
  auc.bic <- result[, 3]
  acc.bic <- result[, 4]
  
  list(
    Accuracy = c(aic.acc = paste.mean.sd(acc.aic), 
                 bic.acc = paste.mean.sd(acc.bic)),
    AUC = c(aic.auc = paste.mean.sd(auc.aic), 
            bic.auc = paste.mean.sd(auc.bic))
  )		
}

# compare AIC, BIC, ridge, lasso, mcp, and scad
cv.compare<-function(data,n2,r){
  # such a cross-validation is called repeated random sub-sampling validation.
  # reference: https://en.wikipedia.org/wiki/Cross-validation_(statistics)
  
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  result <- foreach (i = 1:r, .combine = rbind) %dopar% {
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    y.test <- y[test]
    
    # selection using AIC
    lmod <- lm(y ~ ., data[train, ]) # fit a full model first
    aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    aic.pred <- predict(aic.mod, newdata = data[test, ])
    err.a <- mean((aic.pred - y.test)^2)
    
    # selection using BIC
    bic.mod <- step(lmod, k = log(nrow(data)),
                    scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    bic.pred <- predict(bic.mod, newdata = data[test, ])
    err.b <- mean((bic.pred - y.test)^2)
    
    # selection using ridge
    ridge.mod <- glmnet(x[train,], y[train], alpha=0)
    cv.out <- cv.glmnet(x[train,], y[train], alpha=0)  
    ridge.pred <- predict(ridge.mod, s = cv.out$lambda.min, newx = x[test,])
    err.r <- mean((ridge.pred - y.test)^2)
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)  
    lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.min, newx=x[test,])
    err.l <- mean((lasso.pred - y.test)^2)
    
    # selection using mcp
    mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                         penalty = 'MCP', nfolds = 5)
    mcp.pred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min)
    err.m <- mean((mcp.pred - y.test)^2)
    
    # selection using scad
    scad.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                          penalty = 'SCAD', nfolds = 5)
    scad.pred <- predict(scad.mod, X = x[test,], lambda = scad.mod$lambda.min)
    err.s <- mean((scad.pred - y.test)^2)
    
    c(err.a, err.b, err.r, err.l, err.m, err.s)
  }
  
  err.a <- result[, 1]
  err.b <- result[, 2]
  err.r <- result[, 3]
  err.l <- result[, 4]
  err.m <- result[, 5]
  err.s <- result[, 6]
  
  par(mfrow = c(2, 3))
  plot(err.a,err.l,xlab="AIC MSPE",ylab="Lasso MSPE")
  abline(0, 1)
  plot(err.b,err.l,xlab="BIC MSPE",ylab="Lasso MSPE")
  abline(0, 1)
  plot(err.r,err.l,xlab="Ridge MSPE",ylab="Lasso MSPE")
  abline(0, 1)
  plot(err.m,err.l,xlab="MCP MSPE",ylab="Lasso MSPE")
  abline(0, 1)
  plot(err.s,err.l,xlab="SCAD MSPE",ylab="Lasso MSPE")
  abline(0, 1)
  plot(err.s,err.m,xlab="SCAD MSPE",ylab="MCP MSPE")
  abline(0, 1)
  par(mfrow = c(1, 1))
  
  c(aic.err = mean(err.a), aic.sd = sd(err.a),
    bic.err = mean(err.b), bic.sd = sd(err.b), 
    ridge.err = mean(err.r), ridge.sd = sd(err.r), 
    lasso.err = mean(err.l), lasso.sd = sd(err.l), 
    mcp.err = mean(err.m), mcp.sd = sd(err.m), 
    scad.err = mean(err.s), scad.sd = sd(err.s))		
}


# compare XGBoost, Random Forest, and SVM
cv.ml.compare<-function(data,n2,r, xgb.params, xgb.nrounds){
  
  # assumes the name of the response column is "y"
  require(xgboost)
  require(randomForest)
  require(e1071)
  require(doParallel)
  registerDoParallel(cores = 4)
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  result <- foreach (i = 1:r, .combine = rbind) %dopar% {
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    X.test <- xgb.DMatrix(as.matrix(data[test, ] %>% select(-y)))
    y.test <- y[test]
    
    # fit an XGBoost model using the tuned parameters
    dtrain <- xgb.DMatrix(data = as.matrix(data[train, ] %>% select(-y)), 
                          label = data[train, ]$y)
    xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
    xgb.probpred <- predict(xgb.mod, X.test)
    xgb.pred <- as.numeric(xgb.probpred > 0.5)
    auc.xgb <- auc_score(xgb.probpred, y.test)
    acc.xgb <- accuracy_score(xgb.pred, y.test)
    
    # fit a Random Forest model
    rf.mod <- randomForest(as.factor(y) ~ ., data = data[train, ])
    rf.probpred <- predict(rf.mod, newdata = data[test, ], type = "prob")[, 2]
    rf.pred <- predict(rf.mod, newdata = data[test, ])
    auc.rf <- auc_score(rf.probpred, y.test)
    acc.rf <- accuracy_score(rf.pred, y.test)
    
    # fit an SVM model
    svm.mod <- svm(formula = as.factor(y) ~ ., data = data[train, ], 
                   type = 'C-classification', kernel = 'linear', 
                   probability = TRUE)
    svm.pred <- predict(svm.mod, newdata = data[test, ], probability = TRUE)
    svm.probpred <- as.data.frame(attr(svm.pred, "probabilities"))$`1`
    auc.svm <- auc_score(svm.probpred, y.test)
    acc.svm <- accuracy_score(svm.pred, y.test)
    
    c(auc.xgb, acc.xgb, auc.rf, acc.rf, auc.svm, acc.svm)
    
  }
  
  auc.xgb <- result[, 1]
  acc.xgb <- result[, 2]
  auc.rf <- result[, 3]
  acc.rf <- result[, 4]
  auc.svm <- result[, 5]
  acc.svm <- result[, 6]
  
  
  list(
    Accuracy = c(xgb.acc = paste.mean.sd(acc.xgb), 
                 rf.acc = paste.mean.sd(acc.rf), 
                 svm.acc = paste.mean.sd(acc.svm)),
    AUC = c(xgb.auc = paste.mean.sd(auc.xgb), 
            rf.auc = paste.mean.sd(auc.rf), 
            svm.auc = paste.mean.sd(auc.svm))
  )				
}

# compare with glm models
cv.glm.compare<-function(data,n2,r, glm.form){
  
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  result <- foreach (i = 1:r, .combine = rbind) %dopar% {
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    y.test <- y[test]
    
    # fit a glm model
    glm.mod <- glm(glm.form, family=binomial, data[train, ])
    glm.probpred <- predict(glm.mod, newdata = data[test, ], type = "response")
    glm.pred <- as.numeric(glm.probpred > 0.5)
    auc.glm <- auc_score(glm.probpred, y.test)
    acc.glm <- accuracy_score(glm.pred, y.test)
    
    c(auc.glm, acc.glm)
  }
  
  auc.glm <- result[, 1]
  acc.glm <- result[, 2]
  
  c(glm.acc = paste.mean.sd(acc.glm), glm.auc = paste.mean.sd(auc.glm))		
}



##### ========== Select Features ============ #####
select.features.regularized <- function(data, seed = 99) {
  # seed is set to guarantee reproducibility
  
  # pass in the data, it will output the features selected by regularized methods.
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  
  # use LASSO to get the selected features
  set.seed(seed)
  cv.out <- cv.glmnet(x,y,family="binomial",alpha=1,nfold=10)
  lambdah <- cv.out$lambda.1se
  lasso.mod <- glmnet(x,y,family="binomial",alpha=1,lambda=lambdah) 
  betah <- lasso.mod$beta
  lasso.features <- betah[which(abs(betah) > 0)]
  names(lasso.features) <- colnames(x)[which(abs(betah) > 0)]
  
  # use MCP to get the selected features
  set.seed(seed)
  mcp.mod <- cv.ncvreg(x,y,family='binomial',penalty='MCP',nfolds=5)
  mcp.fit <- mcp.mod$fit
  coeff.mcp <- mcp.fit$beta[,mcp.mod$min]
  mcp.features <- coeff.mcp[which(abs(coeff.mcp)>0)][-1]
  
  # use SCAD to get the selected features
  set.seed(seed)
  scad.mod <- cv.ncvreg(x,y,family='binomial',penalty='SCAD',nfolds=5)
  scad.fit <- scad.mod$fit
  coeff.scad <- scad.fit$beta[,scad.mod$min]
  scad.features <- coeff.scad[which(abs(coeff.scad)>0)][-1]
  
  list(lasso.features = lasso.features, 
       mcp.features = mcp.features, 
       scad.features = scad.features)
}

select.features.abic <- function(data) {
  
  lmod <- glm(y ~ ., family=binomial, data)
  
  # use AIC to get the selected features
  aic.mod <- step(lmod, k = 2, scope = list(upper=~.,lower=~1), 
                  direction = 'both', trace=FALSE)
  
  # use BIC to get the selected features
  bic.mod <- step(lmod, k = log(nrow(data)), 
                  scope = list(upper=~.,lower=~1), 
                  direction = 'both', trace=FALSE)
  
  list(aic.mod = aic.mod, bic.mod = bic.mod)
}

# check feature importance by SOIL, Random Forest, and XGBoost
check.importance <- function(data, xgb.params, xgb.nrounds) {
  
  # SOIL importance
  # print("SOIL important features")
  # print(soil(data))
  
  # Random Forest importance
  rf.mod <- randomForest(as.factor(y) ~ ., data, importance=TRUE)
  print("Random Forest Importance Measure")
  print(importance(rf.mod))
  print("Check varImpPlot by Random Forest...")
  varImpPlot(rf.mod, scale = FALSE)
  
  # XGBoost importance
  dtrain <- xgb.DMatrix(data = as.matrix(data %>% select(-y)), 
                        label = data$y)
  xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
  mat <- xgb.importance(feature_names = colnames(data %>% select(-y)), 
                        model = xgb.mod)
  print("Check Importance Plot by XGBoost...")
  xgb.plot.importance(importance_matrix = mat)
}


##### ========== Data Manipulation ============ #####

# obtain model matrix and dataframe with order 2 interactions
order2.model.matrix <- function(data) {
  ## assumes the response name is "y"
  require(foreach)
  # get x matrix
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  # get all but quadratic terms
  x.inter <- model.matrix(rep(1, length(y)) ~ .^2, data = data)
  
  # get quadratic terms manually
  x.q <- foreach (i = 1:ncol(x), .combine = cbind) %do% {
    x[, i] * x[, i]
  }
  
  # assign names to the quadratic terms
  x.colnames <- colnames(x)
  x.q.colnames <- vector(length = length(x.colnames))
  for (i in 1:length(x.colnames)) {
    x.q.colnames[i] <- paste(x.colnames[i], "^2", sep = "")
  }
  colnames(x.q) <- x.q.colnames
  
  # combine x.inter and x.q
  x.order2 <- cbind(x.inter[, -1], x.q) # don't include intercept
  
  list(X = x.order2, df = data.frame(cbind(y, x.order2)))
}



##### ========== XGBoost ============ #####
tune.xgb.reg <- function(X_train, y_train) {
  # use cv to tune parameters for xgboost regressor
  # X_train should be xgb.DMatrix
  require(caret)
  require(xgboost)
  
  xgb_trcontrol <- trainControl(
    method = "cv", 
    number = 5, 
    allowParallel = TRUE, 
    verboseIter = FALSE, 
    returnData = FALSE
  )
  
  # may change the grid based on different datasets
  xgbGrid <- expand.grid(nrounds = c(50, 100, 150, 200), 
                         max_depth = c(2, 3, 4, 5), 
                         colsample_bytree = seq(0.5, 0.9, length.out = 5), 
                         eta = c(0.01, 0.1, 0.2), 
                         gamma = 0, 
                         min_child_weight = 1, 
                         subsample = c(0.6, 0.7, 0.8, 1))
  
  set.seed(0)
  xgb_model <- train(
    X_train, y_train,
    trControl = xgb_trcontrol, 
    tuneGrid = xgbGrid, 
    method = "xgbTree"
  )
  
  # the optimal parameter
  xgb_model$bestTune
}


##### ========== Others ============ #####
# get SOIL important features
soil <- function(data) {
  require(SOIL)
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  a=SOIL(x, y, n_train = ceiling(n/2), no_rep = 100,
         n_train_bound = ceiling(n/2)- 2, n_bound = n - 2,
         psi = 1, family = "gaussian", method = "union",
         weight_type = "ARM", prior = TRUE)
  
  out <- a$importance[, order(a$importance, decreasing = TRUE)]
  out[out > 0]
}

# Calculate AIC, BIC, and AICc
AIC.BIC <- function(lmod) {
  
  n <- length(lmod$residuals)
  npar <- n - lmod$df.residual
  aic <- extractAIC(lmod, k = 2)
  bic <- extractAIC(lmod, k = log(n))
  aicc <- extractAIC(lmod,k = 2) + 2*npar*(npar+1)/(n-npar-1)  
  
  list(aic=aic, bic=bic, aicc=aicc)
}


##### ========== Prediction ============ #####

# LOO to estimate the performance of linear models
loo.lm <- function(data, lm.form){
  # use leave-one-out to estimate the accuracy of the model
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  result <- foreach (i = 1:n, .combine = rbind) %dopar% {
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # fit a glm model
    lm.mod <- glm(lm.form, family = binomial, data[train, ])
    lm.pred <- predict(lm.mod, newdata = data[test, ])
    mse.lm <- mean((lm.pred - y.test)^2)
    mae.lm <- mean(abs(lm.pred - y.test))
    
    c(mse.lm, mae.lm)
  }
  
  mse.lm <- result[, 1]
  mae.lm <- result[, 2]
  
  c(lm.mse = mean(mse.lm), lm.mse.sd = sd(mse.lm), 
    lm.mae = mean(mae.lm), lm.mae.sd = sd(mae.lm))		
}


loo.avg <- function(data, glm.form, bic.form, lasso.form, 
                    mcp.form){
  # use leave-one-out to estimate the accuracy of the avg of lm and lasso
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n <- nrow(data)
  y <- data$y
  result <- foreach (i = 1:n, .combine = rbind) %dopar% {
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # prediction by glm
    glm.mod <- glm(glm.form, family=binomial, data[train, ])
    glm.probpred <- predict(glm.mod, newdata = data[test, ], type = "response")
    glm.pred <- as.numeric(glm.probpred > 0.5)
    auc.glm <- auc_score(glm.probpred, y.test)
    acc.glm <- accuracy_score(glm.pred, y.test)
    
    # prediction by BIC
    bic.mod <- glm(bic.form, family=binomial, data[train, ])
    bic.probpred <- predict(bic.mod, newdata = data[test, ], type = "response")
    bic.pred <- as.numeric(bic.probpred > 0.5)
    auc.bic <- auc_score(bic.probpred, y.test)
    acc.bic <- accuracy_score(bic.pred, y.test)
    
    # prediction by lasso
    lasso.mod <- glm(lasso.form, family=binomial, data[train, ])
    lasso.probpred <- predict(lasso.mod, newdata = data[test, ], type = "response")
    lasso.pred <- as.numeric(lasso.probpred > 0.5)
    auc.lasso <- auc_score(lasso.probpred, y.test)
    acc.lasso <- accuracy_score(lasso.pred, y.test)
    
    # prediction by mcp
    mcp.mod <- glm(mcp.form, family=binomial, data[train, ])
    mcp.probpred <- predict(mcp.mod, newdata = data[test, ], type = "response")
    mcp.pred <- as.numeric(mcp.probpred > 0.5)
    auc.mcp <- auc_score(mcp.probpred, y.test)
    acc.mcp <- accuracy_score(mcp.pred, y.test)
    
    # pred <- mean(c(lm.pred, bic.pred, lasso.pred, mcp.pred))
    # use max voting
    # avg.pred <- as.numeric(((glm.pred+bic.pred+lasso.pred+mcp.pred)/4) > 0.5)
    avg.probpred <- (glm.probpred + bic.probpred + lasso.probpred + mcp.probpred) / 4
    avg.pred <- as.numeric(avg.probpred > 0.5)
    auc.avg <- auc_score(avg.probpred, y.test)
    acc.avg <- accuracy_score(avg.pred, y.test)
    
    c(auc.glm, acc.glm, auc.bic, acc.bic, auc.lasso, acc.lasso, 
      auc.mcp, acc.mcp, auc.avg, acc.avg)
  }
  
  auc.glm <- result[, 1]
  acc.glm <- result[, 2]
  auc.bic <- result[, 3]
  acc.bic <- result[, 4]
  auc.lasso <- result[, 5]
  acc.lasso <- result[, 6]
  auc.mcp <- result[, 7]
  acc.mcp <- result[, 8]
  auc.avg <- result[, 9]
  acc.avg <- result[, 10]
  
  list(
    Accuracy = c(glm.acc = paste.mean.sd(acc.glm), 
                 bic.acc = paste.mean.sd(acc.bic), 
                 lasso.acc = paste.mean.sd(acc.lasso), 
                 mcp.acc = paste.mean.sd(acc.mcp),
                 avg.acc = paste.mean.sd(acc.avg)),
    AUC = c(glm.auc = paste.mean.sd(auc.glm), 
            bic.auc = paste.mean.sd(auc.bic), 
            lasso.auc = paste.mean.sd(auc.lasso), 
            mcp.auc = paste.mean.sd(auc.mcp),
            avg.auc = paste.mean.sd(auc.avg))
  )	
}


pred.avg <- function(data, newdata, glm.form, bic.form, lasso.form, 
                     mcp.form){
  # use leave-one-out to estimate the accuracy of the avg of lm and lasso
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  # prediction by glm
  glm.mod <- glm(glm.form, family=binomial, data)
  glm.probpred <- predict(glm.mod, newdata = newdata, type = "response")
  
  # prediction by BIC
  bic.mod <- glm(bic.form, family=binomial, data)
  bic.probpred <- predict(bic.mod, newdata = newdata, type = "response")
  
  # prediction by lasso
  lasso.mod <- glm(lasso.form, family=binomial, data)
  lasso.probpred <- predict(lasso.mod, newdata = newdata, type = "response")
  
  # prediction by mcp
  mcp.mod <- glm(mcp.form, family=binomial, data)
  mcp.probpred <- predict(mcp.mod, newdata = newdata, type = "response")
  
  # pred <- mean(c(lm.pred, bic.pred, lasso.pred, mcp.pred))
  # use max voting
  # avg.pred <- as.numeric(((glm.pred+bic.pred+lasso.pred+mcp.pred)/4) > 0.5)
  avg.probpred <- (glm.probpred + bic.probpred + lasso.probpred + mcp.probpred) / 4
  avg.pred <- as.numeric(avg.probpred > 0.5)
  
  avg.pred  
}


##### ========== Obsolete Functions ============ #####

# in the last stage, we should use regularized methods and AIC/BIC to select 
# linear models, and use LOO to estimate the performance of the linear 
# models. Therefore, it doesn't make much sense to use loo.regularized and
# loo.abic.
loo.regularized <- function(data){
  # use LOO to estimate the performance of regularized methods.
  
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  result <- foreach (i = 1:n, .combine = rbind) %dopar% {
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # selection using ridge
    ridge.mod <- glmnet(x[train,], y[train], alpha=0)
    cv.out <- cv.glmnet(x[train,], y[train], alpha=0)
    ridge.pred <- predict(ridge.mod, s = cv.out$lambda.min, 
                          newx = matrix(x[test,], nrow = 1))
    mse.r <- mean((ridge.pred - y.test)^2)
    mae.r <- mean(abs(ridge.pred - y.test))
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)
    lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.min, 
                          newx = matrix(x[test,], nrow = 1))
    mse.l <- mean((lasso.pred - y.test)^2)
    mae.l <- mean(abs(lasso.pred - y.test))
    
    # selection using mcp
    mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                         penalty = 'MCP', nfolds = 5)
    mcp.pred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min)
    mse.m <- mean((mcp.pred - y.test)^2)
    mae.m <- mean(abs(mcp.pred - y.test))
    
    # selection using scad
    scad.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                          penalty = 'SCAD', nfolds = 5)
    scad.pred <- predict(scad.mod, X = x[test,], lambda = scad.mod$lambda.min)
    mse.s <- mean((scad.pred - y.test)^2)
    mae.s <- mean(abs(scad.pred - y.test))
    
    c(mse.r, mae.r, mse.l, mae.l, mse.m, mae.m, mse.s, mae.s)
  }
  
  mse.r <- result[, 1]
  mae.r <- result[, 2]
  mse.l <- result[, 3]
  mae.l <- result[, 4]
  mse.m <- result[, 5]
  mae.m <- result[, 6]
  mse.s <- result[, 7]
  mae.s <- result[, 8]
  
  list(
    MSE = c(ridge.mse = mean(mse.r), ridge.mse.sd = sd(mse.r),
            lasso.mse = mean(mse.l), lasso.mse.sd = sd(mse.l),
            mcp.mse = mean(mse.m), mcp.mse.sd = sd(mse.m),
            scad.mse = mean(mse.s), scad.mse.sd = sd(mse.s)),
    
    MAE = c(ridge.mae = mean(mae.r), ridge.mae.sd = sd(mae.r),
            lasso.mae = mean(mae.l), lasso.mae.sd = sd(mae.l),
            mcp.mae = mean(mae.m), mcp.mae.sd = sd(mae.m),
            scad.mae = mean(mae.s), scad.mae.sd = sd(mae.s))
  )
  
}

loo.abic<-function(data){
  # use loo to estimate the performance of AIC, BIC, both direction secection
  
  # the reason I split the comparison into two functions is that for
  # some dataset, AIC and BIC cannot work due to large dimension size
  
  # assumes the name of the response column is "y"
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n = nrow(data)
  y = data$y
  result <- foreach (i = 1:n, .combine = rbind) %dopar% {
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # selection using AIC, both
    lmod <- lm(y ~ ., data[train, ]) # fit a full model first
    aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    aic.pred <- predict(aic.mod, newdata = data[test, ])
    mse.a <- mean((aic.pred - y.test)^2)
    mae.a <- mean(abs(aic.pred - y.test))
    
    # selection using BIC, both
    bic.mod <- step(lmod, k = log(nrow(data)),
                    scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    bic.pred <- predict(bic.mod, newdata = data[test, ])
    mse.b <- mean((bic.pred - y.test)^2)
    mae.b <- mean(abs(bic.pred - y.test))
    
    c(mse.a, mae.a, mse.b, mae.b)
  }
  
  mse.a <- result[, 1]
  mae.a <- result[, 2]
  mse.b <- result[, 3]
  mae.b <- result[, 4]
  
  list(
    MSE = c(aic.mse = mean(mse.a), aic.mse.sd = sd(mse.a),
            bic.mse = mean(mse.b), bic.mse.sd = sd(mse.b)),
    
    MAE = c(aic.mae = mean(mae.a), aic.mae.sd = sd(mae.a),
            bic.mae = mean(mae.b), bic.mae.sd = sd(mae.b))
  )
}




# this function is obsolete, just for future reference
loo.avg.old <- function(lm.data, lm.form, bic.data, bic.form, 
                        lasso.data, lasso.form, mcp.data, mcp.form){
  # use leave-one-out to estimate the accuracy of the avg of lm and lasso
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n <- nrow(lm.data)
  y <- lm.data$y
  result <- foreach (i = 1:n, .combine = rbind) %dopar% {
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # prediction by lm
    lm.mod <- lm(lm.form, lm.data[train, ])
    lm.pred <- predict(lm.mod, newdata = lm.data[test, ])
    mse.lm <- mean((lm.pred - y.test)^2)
    mae.lm <- mean(abs(lm.pred - y.test))
    
    # prediction by BIC
    bic.mod <- lm(bic.form, bic.data[train, ])
    bic.pred <- predict(bic.mod, newdata = bic.data[test, ])
    mse.bic <- mean((bic.pred - y.test)^2)
    mae.bic <- mean(abs(bic.pred - y.test))
    
    # prediction by lasso
    lasso.mod <- lm(lasso.form, lasso.data[train, ])
    lasso.pred <- predict(lasso.mod, newdata = lasso.data[test, ])
    mse.lasso <- mean((lasso.pred - y.test)^2)
    mae.lasso <- mean(abs(lasso.pred - y.test))
    
    # prediction by mcp
    mcp.mod <- lm(mcp.form, mcp.data[train, ])
    mcp.pred <- predict(mcp.mod, newdata = mcp.data[test, ])
    mse.mcp <- mean((mcp.pred - y.test)^2)
    mae.mcp <- mean(abs(mcp.pred - y.test))
    
    # pred <- mean(c(lm.pred, bic.pred, lasso.pred, mcp.pred))
    avg.pred <- 0.3*lm.pred + 0.4*lasso.pred + 0.15*bic.pred + 0.15*mcp.pred
    mse.avg <- mean((avg.pred - y.test)^2)
    mae.avg <- mean(abs(avg.pred - y.test))
    
    c(mse.lm, mae.lm, mse.bic, mae.bic, mse.lasso, mae.lasso, 
      mse.mcp, mae.mcp, mse.avg, mae.avg)
  }
  
  mse.lm <- result[, 1]
  mae.lm <- result[, 2]
  mse.bic <- result[, 3]
  mae.bic <- result[, 4]
  mse.lasso <- result[, 5]
  mae.lasso <- result[, 6]
  mse.mcp <- result[, 7]
  mae.mcp <- result[, 8]
  mse.avg <- result[, 9]
  mae.avg <- result[, 10]
  
  list(
    MSE = c(lm.mse = mean(mse.lm), lm.mse.sd = sd(mse.lm),
            bic.mse = mean(mse.bic), bic.mse.sd = sd(mse.bic),
            lasso.mse = mean(mse.lasso), lasso.mse.sd = sd(mse.lasso),
            mcp.mse = mean(mse.mcp), mcp.mse.sd = sd(mse.mcp), 
            avg.mse = mean(mse.avg), avg.mse.sd = sd(mse.avg)),
    
    MAE = c(lm.mae = mean(mae.lm), lm.mae.sd = sd(mae.lm),
            bic.mae = mean(mae.bic), bic.mae.sd = sd(mae.bic),
            lasso.mae = mean(mae.lasso), lasso.mae.sd = sd(mae.lasso),
            mcp.mae = mean(mae.mcp), mcp.mae.sd = sd(mae.mcp), 
            avg.mae = mean(mae.avg), avg.mae.sd = sd(mae.avg))
  )		
}
