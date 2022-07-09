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
  par(mfroo = c(1, 1))
}


# check outliers
outlier.check <- function(lmod) {
  par(mfrow = c(2, 3))
  plot(lmod, which = 1)
  plot(lmod, which = 3)
  plot(lmod, which = 2)
  halfnorm(cooks.distance(lmod))
  plot(lmod, which = 4)
  plot(lmod, which = 5)
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


##### ========== Plot Fit ============ #####
plot.pred <- function(predicted, y_test, title="") {
  
  options(repr.plot.width=8, repr.plot.height=4)
  my_data = as.data.frame(cbind(predicted = predicted,
                                observed = y_test))
  
  # Plot predictions vs test data
  ggplot(my_data,aes(predicted, observed)) + geom_point(color = "darkred", alpha = 0.5) + 
    geom_smooth(method=lm)+ ggtitle('Linear Regression ') + ggtitle(paste(title, "Prediction vs Test Data")) +
    xlab("Predecited Power Output ") + ylab("Observed Power Output") + 
    theme(plot.title = element_text(color="darkgreen",size=16,hjust = 0.5),
          axis.text.y = element_text(size=12), axis.text.x = element_text(size=12,hjust=.5),
          axis.title.x = element_text(size=14), axis.title.y = element_text(size=14))
}


plot.regularized <- function(data){
  
  require(ggplot2)
  require(ggpubr)
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  n <- nrow(x)
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # selection using ridge
  ridge.mod <- glmnet(x[train,], y[train], alpha=0)
  cv.out <- cv.glmnet(x[train,], y[train], alpha=0)  
  ridge.pred <- predict(ridge.mod, s = cv.out$lambda.1se, newx = x[test,])
  ridge.plot <- plot.pred(ridge.pred, y.test, "Ridge")
  err.r <- mean((ridge.pred - y.test)^2)
  
  # selection using lasso
  lasso.mod <- glmnet(x[train,], y[train], alpha=1)
  cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)  
  lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.1se, newx=x[test,])
  lasso.plot <- plot.pred(lasso.pred, y.test, "LASSO")
  err.l <- mean((lasso.pred - y.test)^2)
  
  # selection using mcp
  mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                       penalty = 'MCP', nfolds = 5)
  mcp.pred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min)
  mcp.plot <- plot.pred(mcp.pred, y.test, "MCP")
  err.m <- mean((mcp.pred - y.test)^2)
  
  # selection using scad
  scad.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                        penalty = 'SCAD', nfolds = 5)
  scad.pred <- predict(scad.mod, X = x[test,], lambda = scad.mod$lambda.min)
  scad.plot <- plot.pred(scad.pred, y.test, "SCAD")
  err.s <- mean((scad.pred - y.test)^2)
  
  print(c(ridge.err = err.r, lasso.err=err.l, 
          mcp.err=err.m, scad.err=err.s))
  
  figure <- ggarrange(ridge.plot, lasso.plot, mcp.plot, scad.plot, 
                      labels = c('A', 'B', 'C', 'D'), 
                      ncol = 2, nrow = 2)
  figure
  
}

plot.xgb.rf <- function(data, xgb.params, xgb.nrounds, 
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

plot.abic <- function(data){
  # assumes the name of the response column is "y"
  require(ggplot2)
  require(ggpubr)
  
  n = nrow(data)
  y = data$y
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # selection using AIC, both
  lmod <- lm(y ~ ., data[train, ]) # fit a full model first
  aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                  direction = 'both', trace=FALSE)
  aic.pred <- predict(aic.mod, newdata = data[test, ])
  aic.both.plot <- plot.pred(aic.pred, y.test, "AIC Both")
  err.a.both <- mean((aic.pred - y.test)^2)
  
  # selection using BIC, both
  bic.mod <- step(lmod, k = log(nrow(data)),
                  scope = list(upper=~.,lower=~1),
                  direction = 'both', trace=FALSE)
  bic.pred <- predict(bic.mod, newdata = data[test, ])
  bic.both.plot <- plot.pred(bic.pred, y.test, "BIC Both")
  err.b.both <- mean((bic.pred - y.test)^2)
  
  # selection using AIC, backward
  lmod <- lm(y ~ ., data[train, ]) # fit a full model first
  aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                  direction = 'backward', trace=FALSE)
  aic.pred <- predict(aic.mod, newdata = data[test, ])
  aic.back.plot <- plot.pred(aic.pred, y.test, "AIC Back")
  err.a.back <- mean((aic.pred - y.test)^2)
  
  # selection using BIC, backward
  bic.mod <- step(lmod, k = log(nrow(data)),
                  scope = list(upper=~.,lower=~1),
                  direction = 'backward', trace=FALSE)
  bic.pred <- predict(bic.mod, newdata = data[test, ])
  bic.back.plot <- plot.pred(bic.pred, y.test, "BIC Back")
  err.b.back <- mean((bic.pred - y.test)^2)
  
  print(c(aic.both.err = err.a.both, bic.both.err=err.b.both, 
          aic.back.err=err.a.back, bic.back.err=err.b.back))
  
  figure <- ggarrange(aic.both.plot, bic.both.plot, 
                      aic.back.plot, aic.back.plot, 
                      labels = c('A', 'B', 'C', 'D'), 
                      ncol = 2, nrow = 2)
  figure
}

plot.lmod <- function(data, lm.form){
  # assumes the name of the response column is "y"
  require(ggplot2)
  require(ggpubr)
  
  n = nrow(data)
  y = data$y
  
  set.seed(99) # set seed to ensure reproducibility
  train <- sample(1:n, round(n/10*7))
  test <- (-train)
  y.test <- y[test]
  
  # selection using AIC, both
  lm.mod <- lm(lm.form, data[train, ]) # fit a full model first
  lm.pred <- predict(lm.mod, newdata = data[test, ])
  lm.plot <- plot.pred(lm.pred, y.test, "LM")
  err.lm <- mean((lm.pred - y.test)^2)
  
  print(c(lm.err = err.lm))
  
  lm.plot
}


##### ========== Cross Validation ============ #####

mse <- function(preds, actuals) {
  mean((preds - actuals)^2)
}

mae <- function(preds, actuals) {
  mean(abs(preds - actuals))
}


paste.mean.sd <- function(a) {
  # paste mean and standard deviation for better visualization
  paste(round(mean(a), 4), "(", round(sd(a), 4), ")", sep = "")
}


# do the idxth comparison and add the results to the comp.df
new.compare <- function(comp.df, idx, data, n2, r, 
                        xgb.params, xgb.nrounds, abic=TRUE){
  # add the new comparison result to comp.df
  if (abic) {
    comp.df[c(1:2), idx] <- cv.abic.compare(data, n2, r)
  } else {
    comp.df[c(1:2), idx] <- 0
  }
  comp.df[c(3:6), idx] <- cv.regularized.compare(data, n2, r)
  comp.df[c(7:8), idx] <- cv.ml.compare(data, n2, r, 
                                               xgb.params, xgb.nrounds)
  
  comp.df
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
    
    # selection using ridge
    ridge.mod <- glmnet(x[train,], y[train], alpha=0)
    cv.out <- cv.glmnet(x[train,], y[train], alpha=0)  
    ridge.pred <- predict(ridge.mod, s = cv.out$lambda.1se, newx = x[test,])
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)  
    lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.1se, newx=x[test,])
    
    # selection using mcp
    mcp.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                         penalty = 'MCP', nfolds = 5)
    mcp.pred <- predict(mcp.mod, X=x[test,], lambda = mcp.mod$lambda.min)
    
    # selection using scad
    scad.mod <- cv.ncvreg(x[train,], y[train], family = 'gaussian',
                          penalty = 'SCAD', nfolds = 5)
    scad.pred <- predict(scad.mod, X = x[test,], lambda = scad.mod$lambda.min)
    
    # calculate the mse error of the predictions
    preds <- cbind(ridge.pred, lasso.pred, mcp.pred, scad.pred)
    err <- apply(preds, 2, function(pred) mse(pred, y.test))
    err
  }
  
  model.names <- c('ridge', 'lasso', 'mcp', 'scad')
  colnames(result) <- sapply(model.names, function(x) paste(x, "err", sep="."))
  apply(result, 2, paste.mean.sd)
  
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
    lmod <- lm(y ~ ., data[train, ]) # fit a full model first
    aic.mod <- step(lmod, scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    aic.pred <- predict(aic.mod, newdata = data[test, ])
    
    # selection using BIC, both
    bic.mod <- step(lmod, k = log(nrow(data)),
                    scope = list(upper=~.,lower=~1),
                    direction = 'both', trace=FALSE)
    bic.pred <- predict(bic.mod, newdata = data[test, ])
    
    # calculate the mse error of the predictions
    preds <- cbind(aic.pred, bic.pred)
    err <- apply(preds, 2, function(pred) mse(pred, y.test))
    err
  }
  
  model.names <- c('aic', 'bic')
  colnames(result) <- sapply(model.names, function(x) paste(x, "err", sep="."))
  apply(result, 2, paste.mean.sd)
}

# compare XGBoost and Random Forest
cv.ml.compare<-function(data,n2,r, xgb.params, xgb.nrounds){
  
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
    X.test <- xgb.DMatrix(as.matrix(data[test, ] %>% select(-y)))
    y.test <- y[test]
    
    # fit an XGBoost model using the tuned parameters
    dtrain <- xgb.DMatrix(data = as.matrix(data[train, ] %>% select(-y)), 
                          label = data[train, ]$y)
    xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
    xgb.pred <- predict(xgb.mod, X.test)
    
    # fit a Random Forest model
    rf.mod <- randomForest(y ~ ., data = data[train, ])
    rf.pred <- predict(rf.mod, newdata = data[test, ])
    
    # calculate the mse error of the predictions
    preds <- cbind(xgb.pred, rf.pred)
    err <- apply(preds, 2, function(pred) mse(pred, y.test))
    err
  }
  
  model.names <- c('xgb', 'rf')
  colnames(result) <- sapply(model.names, function(x) paste(x, "err", sep="."))
  apply(result, 2, paste.mean.sd)		
}

# compare with linear models
cv.lm.compare <- function(data,n2,r, lm.forms, lm.names){
  # use leave-one-out to estimate the accuracy of the avg of glm models
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n <- nrow(data)
  y <- data$y
  
  result <- foreach (i = 1:r, .combine = 'rbind', 
                     .packages = c('foreach', 'stats')) %dopar% {
    
    set.seed(i) # set seed to ensure reproducibility
    train <- sample(1:n, n-n2)
    test <- (-train)
    y.test <- y[test]
    
    # prediction by lm models
    preds <- foreach(lm.form=lm.forms, .combine = 'cbind') %do% {
     lm.mod <- lm(lm.form, data[train, ])
     lm.pred <- predict(lm.mod, newdata = data[test, ])
     lm.pred
    }
    
    # calculate the mse error of the predictions
    err <- apply(preds, 2, function(pred) mse(pred, y.test))
    err
  }

  model.names <- lm.names
  colnames(result) <- sapply(model.names, function(x) paste(x, "err", sep="."))
  apply(result, 2, paste.mean.sd)	
  
}




##### ========== Select Features ============ #####
select.features.regularized <- function(data, seed = 99) {
  # seed is set to guarantee reproducibility
  
  # pass in the data, it will output the features selected by regularized methods.
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  
  # use LASSO to get the selected features
  set.seed(seed)
  cv.out <- cv.glmnet(x,y,alpha=1,nfold=10)
  lambdah <- cv.out$lambda.1se
  lasso.mod <- glmnet(x,y,alpha=1,lambda=lambdah) 
  betah <- lasso.mod$beta
  lasso.features <- betah[which(abs(betah) > 0)]
  names(lasso.features) <- colnames(x)[which(abs(betah) > 0)]
  
  # use MCP to get the selected features
  set.seed(seed)
  mcp.mod <- cv.ncvreg(x,y,family='gaussian',penalty='MCP',nfolds=5)
  mcp.fit <- mcp.mod$fit
  coeff.mcp <- mcp.fit$beta[,mcp.mod$min]
  mcp.features <- coeff.mcp[which(abs(coeff.mcp)>0)][-1]
  
  # use SCAD to get the selected features
  set.seed(seed)
  scad.mod <- cv.ncvreg(x,y,family='gaussian',penalty='SCAD',nfolds=5)
  scad.fit <- scad.mod$fit
  coeff.scad <- scad.fit$beta[,scad.mod$min]
  scad.features <- coeff.scad[which(abs(coeff.scad)>0)][-1]
  
  list(lasso.features = lasso.features, 
       mcp.features = mcp.features, 
       scad.features = scad.features)
}

select.features.abic <- function(data) {
  
  lmod <- lm(y ~ ., data)
  
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
  print("SOIL important features")
  print(soil(data))
  
  # Random Forest importance
  rf.mod <- randomForest(y ~ ., data, importance=TRUE)
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
order2.df <- function(data) {
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
  
  data.frame(cbind(y, x.order2))
}


quadratic.df <- function(data) {
  ## assumes the response name is "y"
  require(foreach)
  # get x matrix
  x <- as.matrix(data[, -which(names(data) %in% c("y"))])
  y <- data$y
  
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
  
  data.frame(cbind(y, x, x.q))
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
loo.lm <- function(data, lm.forms, lm.names){
  # use leave-one-out to estimate the accuracy of the avg of glm models
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n <- nrow(data)
  y <- data$y
  
  result <- foreach (i = 1:n, .combine = 'rbind', 
                     .packages = c('foreach', 'stats')) %dopar% {
    
    test <- c(i)
    train <- c(-i)
    y.test <- y[test]
    
    # prediction by lm models
    preds <- foreach(lm.form=lm.forms, .combine = 'cbind') %do% {
     lm.mod <- lm(lm.form, data[train, ])
     lm.pred <- predict(lm.mod, newdata = data[test, ])
     lm.pred
    }
    
    mse <- apply(preds, 2, function(pred) mse(pred, y.test))
    mae <- apply(preds, 2, function(pred) mae(pred, y.test))
    c(mse, mae)
  
  }
  
  model.names <- lm.names
  colnames(result) <- c(sapply(model.names, function(x) paste(x, "mse", sep=".")), 
                        sapply(model.names, function(x) paste(x, "mae", sep=".")))
  apply(result, 2, paste.mean.sd)	
}


loo.avg <- function(data, lm.forms, lm.names, ridge.data,  
                    xgb.params, xgb.nrounds, avg.weights){
  # use leave-one-out to estimate the accuracy of the avg of lm and lasso
  # data here should be the whole training set
  # assumes the name of the response column is "y"
  
  require(doParallel)
  registerDoParallel(cores = 4)
  
  n <- nrow(data)
  y <- data$y
  
  x <- as.matrix(ridge.data[, -which(names(ridge.data) %in% c("y"))])
  
  # tune lambda for ridge regression
  cv.out <- cv.glmnet(x, y, alpha=0, nfold=10)
  ridge.lambda <- cv.out$lambda.1se
  names(ridge.lambda) <- "ridge.lambda"
  print(ridge.lambda)
  
  result <- foreach (i = 1:n, .combine = 'rbind', 
                     .packages = c('foreach', 'stats', 'xgboost', 
                                   'randomForest', 'e1071')) %dopar% {
    test <- c(i)
    train <- c(-i)
    X.test <- xgb.DMatrix(as.matrix(data[test, ] %>% select(-y)))
    y.test <- y[test]
    
    # prediction by lm models
    lm.preds <- foreach(lm.form=lm.forms, .combine = 'cbind') %do% {
     lm.mod <- lm(lm.form, data[train, ])
     lm.pred <- predict(lm.mod, newdata = data[test, ])
     lm.pred
    }
    
    ridge.pred <- rep(0, length(y.test))
    xgb.pred <- rep(0, length(y.test))
    rf.pred <- rep(0, length(y.test))
    svm.pred <- rep(0, length(y.test))
    # prediction by Ridge
    ridge.mod <- glmnet(x[train, ], y[train], alpha=0, lambda=ridge.lambda)
    ridge.pred <- drop(predict(ridge.mod, newx = matrix(x[test, ], 
                                                           nrow = 1)))
    
    # prediction by XGBoost
    dtrain <- xgb.DMatrix(data = as.matrix(data[train, ] %>% select(-y)), 
                         label = data[train, ]$y)
    xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
    xgb.pred <- predict(xgb.mod, X.test)
    
    # prediction by RF
    rf.mod <- randomForest(y ~ ., data = data[train, ])
    rf.pred <- predict(rf.mod, newdata = data[test, ])
    
    preds <- cbind(lm.preds, ridge.pred, xgb.pred, rf.pred)
    avg.pred <- preds %*% avg.weights
    preds <- cbind(preds, avg.pred)
    
    mse <- apply(preds, 2, function(pred) mse(pred, y.test))
    mae <- apply(preds, 2, function(pred) mae(pred, y.test))
    c(mse, mae)
    
  }
  
  model.names <- c(lm.names, 'ridge', 'xgb', 'rf', 'avg')
  colnames(result) <- c(sapply(model.names, function(x) paste(x, "mse", sep=".")), 
                        sapply(model.names, function(x) paste(x, "mae", sep=".")))
  apply(result, 2, paste.mean.sd)	
}


pred.avg <- function(data, newdata, lm.forms, lm.names, 
                     ridge.data, ridge.newdata, ridge.lambda, 
                     xgb.params, xgb.nrounds, avg.weights){
  # assumes the name of the response column is "y"
  require(foreach)
  # prediction by lm models
  lm.preds <- foreach(lm.form=lm.forms, .combine = 'cbind') %do% {
    lm.mod <- lm(lm.form, data)
    lm.pred <- predict(lm.mod, newdata = newdata)
    lm.pred
  }
  
  ridge.pred <- rep(0, nrow(newdata))
  xgb.pred <- rep(0, nrow(newdata))
  rf.pred <- rep(0, nrow(newdata))
  svm.pred <- rep(0, nrow(newdata))
  
  x <- as.matrix(ridge.data[, -which(names(ridge.data) %in% c("y"))])
  y <- ridge.data$y
  x.test <- as.matrix(ridge.newdata[, -which(names(ridge.newdata) %in% c("y"))])
  
  # prediction by Ridge
  ridge.mod <- glmnet(x, y, alpha=0, lambda=ridge.lambda)
  ridge.pred <- drop(predict(ridge.mod, newx = x.test))
  
  # prediction by XGBoost
  X.test <- xgb.DMatrix(as.matrix(newdata %>% select(-y)))
  dtrain <- xgb.DMatrix(data = as.matrix(data %>% select(-y)),
                        label = data$y)
  xgb.mod <- xgb.train(params = xgb.params, dtrain, nrounds = xgb.nrounds)
  xgb.pred <- predict(xgb.mod, X.test)
  
  # prediction by RF
  rf.mod <- randomForest(y ~ ., data = data)
  rf.pred <- predict(rf.mod, newdata = newdata)
  
  preds <- cbind(lm.preds, ridge.pred, xgb.pred, rf.pred)
  avg.pred <- preds %*% avg.weights
  preds <- cbind(preds, avg.pred)
  
  pred.names <- c(lm.names, 'ridge', 'xgb', 'rf', 'avg')
  par(mfrow=c(2, 2))
  for (i in 2:ncol(preds)) {
    plot(preds[, i], preds[, 1], xlab = pred.names[i], ylab = pred.names[1], 
         main = "Prediction Comparison")
  }
  par(mfrow=c(1, 1))
  
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
    ridge.pred <- predict(ridge.mod, s = cv.out$lambda.1se, 
                          newx = matrix(x[test,], nrow = 1))
    mse.r <- mean((ridge.pred - y.test)^2)
    mae.r <- mean(abs(ridge.pred - y.test))
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)
    lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.1se, 
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


######=========== Obsolete Functions ==============########

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
    ridge.pred <- predict(ridge.mod, s = cv.out$lambda.1se, newx = x[test,])
    err.r <- mean((ridge.pred - y.test)^2)
    
    # selection using lasso
    lasso.mod <- glmnet(x[train,], y[train], alpha=1)
    cv.out.l <- cv.glmnet(x[train,], y[train], alpha=1)  
    lasso.pred <- predict(lasso.mod, s = cv.out.l$lambda.1se, newx=x[test,])
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

