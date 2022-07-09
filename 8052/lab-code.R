##### Lab 1 #####
set.seed(8052)
library(AlgDesign)
library(gmodels)
library(gtools)

df <- read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/ex3.5", header = TRUE)
df$trt = as.factor(df$trt)
fit1 = lm(deg ~ trt, data = df)
anova(fit1)
