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
%IncMSE IncNodePurity
crim     7.8630295    1076.63281
zn       0.2323719      72.70382
indus    3.5644244     996.59033
chas     0.2227144      51.26195
nox      6.0449609    1152.46188
rm      33.3402561    6046.40570
age      3.1135431     532.15045
dis      5.6931708    1297.82972
rad      0.9306152     100.38939
tax      2.3064281     421.78975
ptratio  4.0510269     928.40562
black    1.3841616     337.93872
lstat   72.5821213    7597.27092


lstat2<-Boston$lstat+0.5*rnorm(nrow(Boston))

Boston2 =data.frame(Boston, lstat2)
rf.boston2=randomForest(medv~.,data=Boston2,subset=train,mtry=14,importance =TRUE)
importance(rf.boston2,scale=FALSE)

# the results are not very stable: two runs:

> 		importance(rf.boston2,scale=FALSE)
%IncMSE IncNodePurity
crim     5.587385306     603.78854
zn       0.003929573      17.00770
indus    0.970085508      93.86411
chas     0.054447330      45.30679
nox      1.300204913     137.48697
rm      37.201557612    4547.09689
age      4.168140444     357.35333
dis      9.499668009    1355.48631
rad      0.935039516      89.80899
tax      1.676764505     181.62171
ptratio  2.033051955     270.83076
black    0.812649214     221.16099
lstat   20.667434736    1722.73097
lstat2  68.050157545    7905.25390

strangely, dis is boosted!
  
  > 		importance(rf.boston2,scale=FALSE)
%IncMSE IncNodePurity
crim     5.51403130     562.21399
zn       0.05622546      22.09538
indus    1.39230759     105.28391
chas     0.03953455      29.25108
nox      1.49042026     130.71349
rm      42.04431132    6028.51543
age      3.70976396     407.63156
dis     10.40776736    1512.92567
rad      0.97653973     111.30825
tax      1.78627132     193.02054
ptratio  2.85141586     312.88908
black    0.34706078     226.01347
lstat   40.00104779    4951.40458
lstat2  42.28665536    2941.91117

dis2<-Boston$dis+0.5*rnorm(nrow(Boston))

Boston3 =data.frame(Boston, dis2)
rf.boston3=randomForest(medv~.,data=Boston3,subset=train,mtry=14,importance =TRUE)

importance(rf.boston3,scale=FALSE)

> 		importance(rf.boston3,scale=FALSE)
%IncMSE IncNodePurity
crim     5.90009264     599.67900
zn       0.04770139      22.95781
indus    1.08764288      91.77337
chas     0.05445325      38.78853
nox      1.30969383     109.10227
rm      45.51288239    6193.91954
age      3.42364286     326.43915
dis     10.01860733    1369.45159
rad      0.94910214     113.96019
tax      2.15889104     200.67480
ptratio  2.74186231     302.13950
black    0.57809258     187.17115
lstat   83.35103856    7622.85699
dis2     1.85905683     272.44531

lstat2<-Boston$lstat+0.5*rnorm(nrow(Boston))

Boston2 =data.frame(Boston, lstat2)
rf.boston2=randomForest(medv~.,data=Boston2,subset=train,mtry=7,importance =TRUE)

importance(rf.boston2,scale=FALSE)

> 		importance(rf.boston2,scale=FALSE)
%IncMSE IncNodePurity
crim     5.9698457     817.90516
zn       0.1592810      48.92280
indus    1.6585188     416.05049
chas     0.0843438      33.75989
nox      4.8682153     674.16628
rm      26.6955362    5303.09305
age      2.7957815     385.23329
dis      3.6537742     906.06272
rad      0.8612008      71.36219
tax      1.4269833     269.60268
ptratio  2.6607913     510.74290
black    0.6703289     250.01348
lstat   34.4657442    4910.29898
lstat2  31.1787020    5196.56410		