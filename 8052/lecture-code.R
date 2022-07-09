##### Lecture Code 1: ANOVA #####
donuts <- read.delim("/Users/yuyang/Documents/Courses/19Spring/STAT8052/data/donut-fat.txt")
attach(donuts)

fat <- as.factor(fat)
m <- lm(y ~ fat)
anova(m)

### Diagnostic plots
qqnorm(residuals(m))
qqline(residuals(m))

plot(jitter(fitted(m)), residuals(m), xlab = "Fitted", ylab = "Residuals",
     main = "Residuals vs. Fitted")

par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))


### Parametrization of the linear mmodel
summary(m)
(X <- model.matrix(~ fat))
names(m)
options(contrasts = c("contr.sum", "contr.poly"))
m1 <- lm(y ~ fat)
summary(m1)
model.matrix(m1)
# restore the default settings
options(contrasts = c("contr.treatment", "contr.poly"))

detach(donuts)
rm(list = ls())

##### Lecture Code 2: Contrasts #####
library(gmodels)
donuts <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/donut-fat.txt", 
                     header = T)
attach(donuts)

donuts$fat <- as.factor(donuts$fat)
m <- lm(y ~ fat, data = donuts)
options(contrasts=c("contr.treatment", "contr.poly"))
### Contrast for mu1 - mu3
fit.contrast(m, "fat", coeff = c(1, 0, -1, 0))
fit.contrast(m, "fat", coeff = c(1, 0, -1, 0), conf.int = 0.95)

### Test several contrasts at once
coeff.matrix <- rbind("1 vs 3" = c(1, 0, -1, 0),
                      "1 vs ave 2, 3" = c(1, -0.5, -0.5, 0))
fit.contrast(m, "fat", coeff = coeff.matrix, conf.int = 0.95)

### Estimates for mu and alpha
aovmod <- aov(y ~ fat)
summary(aovmod)

# estimates of alpha
model.tables(aovmod, type = "effects")

# estimates of mu
model.tables(aovmod, type = "means")

### using estimable functions
# mu1
estimable(m, cm = c(1, 0, 0, 0), conf.int = 0.95)
# mu2
estimable(m, cm = c(1, 1, 0, 0), conf.int = 0.95)
# all mu
estimable(m, cm = rbind("m1" = c(1, 0, 0, 0),
                        "m2" = c(1, 1, 0, 0),
                        "m3" = c(1, 0, 1, 0),
                        "m4" = c(1, 0, 0, 1)), conf.int = 0.95)

### Change model.matrix to do contrast, calculate SSw

## 1 vs ave 2, 3, 4
contrasts(fat, 1) <- c(-3, 1, 1, 1)
m5 <- lm(y ~ fat)
anova(m5)
model.matrix(m5)

## 2 vs ave 3, 4
contrasts(fat, 1) <- c(0, -2, 1, 1)
m6 <- lm(y ~ fat)
anova(m6)

## 3 vs 4
contrasts(fat, 1) <- c(0, 0, -1, 1)
m7 <- lm(y ~ fat)
anova(m7)

# return contrasts to R default
options(contrasts = c("contr.treatment", "contr.poly"))


##### Lecture Code Chapter 5: Multiple comparisons #####
donuts <- read.delim("/Users/yuyang/Documents/Courses/19Spring/STAT8052/data/donut-fat.txt")
attach(donuts)
donuts$fat <- as.factor(donuts$fat)

### 1. Tukey's HSD
aovmod <- aov(y ~ fat, data = donuts)

(cis = TukeyHSD(aovmod, which = "fat", ordered = TRUE, conf.level = 0.95))
plot(cis)

### 2. Scheffe, Bonferroni and unadjusted t CIs
library(cfcdae)
m <- lm(y ~ fat)

contr.matrix <- cbind("c1" = c(-3, 1, 1, 1), 
                      "c2" = c(0, -2, 1, 1), 
                      "c3" = c(0, 0, -1, 1))
linear.contrast(m, term = fat, contr.coefs = contr.matrix, scheffe = TRUE, 
                confidence = 0.95)

linear.contrast(m, term = fat, contr.coefs = contr.matrix, bonferroni = TRUE, 
                confidence = 0.95)

linear.contrast(m, term = fat, contr.coefs = contr.matrix, confidence = 0.95)

# CIs and p-values based on t-distribution
library(gmodels)
fit.contrast(m, "fat", coeff = t(contr.matrix), conf.int = 0.95)

contr.matrix <- cbind("1 - 2" = c(1, -1, 0, 0), 
                      "1 - 3" = c(1, 0, -1, 0), 
                      "1 - 4" = c(1, 0, 0, -1), 
                      "2 - 3" = c(0, 1, -1, 0), 
                      "2 - 4" = c(0, 1, 0, -1), 
                      "3 - 4" = c(0, 0, 1, -1))
linear.contrast(m, term = fat, contr.coefs = contr.matrix, bonferroni = TRUE, 
                confidence = 0.95)
linear.contrast(m, term = fat, contr.coefs = contr.matrix, scheffe = TRUE, 
                confidence = 0.95)
fit.contrast(m, "fat", coeff = t(contr.matrix)[1:3, ], conf.int = 0.95)

# all pairwise comparisons with t-critical values
linear.contrast(m, term = fat, allpairs = TRUE, confidence = 0.95)
linear.contrast(m, term = fat, allpairs = TRUE, confidence = 0.95, bonferroni = T)
linear.contrast(m, term = fat, allpairs = TRUE, confidence = 0.95, scheffe = T)



### 3. Pairwise comparisons using pairwise function
pairwise(m, term = fat, confidence = 0.95, type = "hsd")
# type could be: hsd, bsd, lsd, regwr, snk

cis <- pairwise(m, term = fat, confidence = 0.95, type = "hsd")
plot(cis)

cis <- pairwise(m, term = fat, confidence = 0.95, type = "regwr")
plot(cis)


### 4. Dunnet's CIs for comparison to one control and comparison to the best
m <- lm(y ~ fat, data = donuts)
compare.to.control(m, term = fat, control.level = 3, confidence = 0.95)
compare.to.best(m, term = fat, lowisbest = F, confidence = 0.95)


### 5. Bonferroni adjustment for comparison with control
linear.contrast(m, term = fat, controlpairs = 3, confidence = 0.95, bonferroni = T)


### 6. use glht{multcomp}
library(multcomp)
m <- lm(y ~ fat, data = donuts)
g <- glht(m, linfct = mcp(fat = contrMat(1:4, type = "Dunnett", base = 3)), 
          alternative = "two.sided")
summary(g)
confint(g)


### 7. Tukey's HSD critical value
# the percentile of studentized range statistic
qtukey(0.95, nmeans = 4, df = 20) / sqrt(2)

abs(qtukey(0.95, nmeans = 2, df = 20) / sqrt(2) - qt(0.975, df = 20)) < 1e-4


##### Lecture Code Chapter 6: Checking assumptions #####
bhh <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/BHH-survival-times.txt",
                  header = TRUE)
names(bhh)
attach(bhh)
trt12 <- as.factor(trt12)

m <- lm(y ~ trt12)
anova(m)

# Diagnostic plots
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))


# Box-cox plot
library(MASS)
boxcox(m)

bc <- boxcox(m)
i <- which.max(bc$y)
bc$x[i]

# 95% CI for lambda
(max.l <- bc$y[i])
n <- nrow(bhh)
(const <- n / 2 * log( 1 + qt(0.025, m$df.residual)^2 / m$df.residual))
lower.l <- max.l - const
i.lower <- min(which(bc$y >= lower.l))
i.upper <- max(which(bc$y >= lower.l))
(ci <- bc$x[c(i.lower, i.upper)])

# ANOVA with transformed data
rate <- 1/y
m2 <- lm(rate ~ trt12)
anova(m2)
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))


##### Lecture Code Chapter 7: Sample Size Selection #####
power.anova.test(groups = 4, n = NULL, between.var = var(c(15/2, -15/2, 0, 0)),
                 within.var = 100, sig.level = 0.01, power = 0.95)
df_E <- 20
(power <- 1 - pf(qf(1-0.05, 1, df_E), 1, df_E, ncp = 3))


##### Lecture Code Chapter 8: Factorial Design #####
bhh <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/BHH-survival-times.txt",
                  header = TRUE)
names(bhh)
bhh$poison <- as.factor(bhh$poison)
bhh$method <- as.factor(bhh$method)
bhh$trt12 <- as.factor(bhh$trt12)
attach(bhh)
rate <- 1/y

options(contrasts = c("contr.treatment", "contr.poly"))

### (1) 2-factor anova model
m <- lm(rate ~ poison * method, data = bhh)
anova(m)
m$contrasts
model.matrix(m)

## 1-factor anova model
anova(lm(rate ~ trt12))

## check the one without interaction
anova(lm(rate ~ poison + method))

### interaction plots
interaction.plot(x.factor = method, trace.factor = poison, response = rate, 
                 fun = mean, legend = T, xpd = NA, type = 'b', 
                 trace.label = 'Poison', xlab = "Method", ylab = "Mean Rate", 
                 lty = c(1, 2, 4), col = c("black", "red", "blue"), 
                 main = "Interaction Plot for rate ~ poison * method")

interaction.plot(x.factor = poison, trace.factor = method, response = rate, 
                 fun = mean, legend = T, xpd = NA, type = "b", 
                 trace.label = 'Method', xlab = "Poison", ylab = "Mean Rate", 
                 lty = 1:4, col = c("black", "red", "blue", "green"), 
                 main = "Interaction Plot for rate ~ poison * method")

# profile plot for method
interaction.plot(x.factor = method, trace.factor = rep(1, length(rate)), 
                 response = rate, fun = mean, legend = F, xpd = NA, 
                 type = "l", xlab = 'Method', ylab = 'Mean Rate', 
                 main = 'Profile Plot for Method')

# pairwise comparisons for method
m.aov <- aov(rate ~ poison * method)
TukeyHSD(m.aov, "method")
sort(model.tables(m.aov, type = "means")$table$method)

# profile plot for poison
interaction.plot(x.factor = poison, trace.factor = rep(1, length(rate)), 
                 response = rate, fun = mean, legend = F, xpd = NA, 
                 type = 'l', xlab = "Poison", ylab = "Mean Rate", 
                 main = "Profile Plot for Poison")

# pairwise comparisons for poison
TukeyHSD(m.aov, "poison")
plot(TukeyHSD(m.aov, "poison"))
sort(model.tables(m.aov, type = "means")$table$poison)


### Checking model assumptions
par(mfrow = c(2, 2))
plot(m)

detach(bhh)

### (2) 3-factor anova model
pwg <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/pig-weight-gain.txt",
                  header = T)
names(pwg)
pwg <- within(pwg, {
  L <- as.factor(L)
  P <- as.factor(P)
  S <- as.factor(S)
  blocks <- as.factor(blocks)
  }
)
attach(pwg)

m <- lm(y ~ blocks + L*P*S, data = pwg)
# check model assumption
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
# anova model 
anova(m)

## interaction plot
interaction.plot(x.factor = L, trace.factor = P, response = y, 
                 fun = mean, type = "b", legend = T, xpd = NA, 
                 trace.label = "Protein", xlab = "Lysine", 
                 ylab = "Mean Weight Gain", lty = 1:2, 
                 col = c('black', 'red'), 
                 main = "Interaction plot for y ~ blocks + L*P*S")

m.aov <- aov(y ~ blocks + L*P*S)
model.tables(m.aov, type = "means")$tables$'L:P'
sort(model.tables(m.aov, type = "means")$tables$'L:P')
TukeyHSD(m.aov, "L:P")

detach(pwg)

### (3) Tukey one-degree-of-freedom test
bc <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/bacteria-counts.txt",
                 header = T)
names(bc)
attach(bc)
day <- as.factor(day)
store <- as.factor(store)

m <- lm(y ~ day * store)
anova(m)  # only one observation in each category

# fit an additive anova model
ma <- lm(y ~ day + store)
anova(ma)
# check model assumptions
m <- lm(y ~ day + store)
par(mfrow = c(2, 2))
plot(m)  # suggesting additive model is not good
# check the interaction plot
par(mfrow = c(1, 1))
interaction.plot(x.factor = day, trace.factor = store, response = y,
                 fun = mean, type = 'b', legend = T, trace.label = 'store', 
                 xpd = NA, xlab = "day", ylab = 'Bacteria count', 
                 lty = 1:5, col = c('red', "black", 'yellow','green','blue'),
                 main = 'Plot of bacteria counts')

# the above procedure suggests some non-additivity, so we consider using
# Tukey's one-degree-of-freedom test for transformable non-additivity

# step 1: fit additive model
ma.aov <- aov(y ~ day + store)

# step 2: calculate covariate for augmented model
ypred <- ma.aov$fitted.values
mu <- model.tables(ma.aov, type = 'means')$table$'Grand mean'
onedf <- ypred * ypred / (2 * mu) # this is a vector

# step 3: fit augmented model
m.onedf <- lm(y ~ day + store + onedf)
anova(m.onedf) # eta is significant

# step 4: If onedf is significant, then do transformation
eta <- m.onedf$coefficients['onedf']
1-eta

# step 5: transform response
logy <- log(y)
mt.aov <- aov(logy ~ day + store)
summary(mt.aov)

# step 6: make sure the new model has no non-additivity
ypred <- mt.aov$fitted.values
mu <- model.tables(mt.aov, type = 'means')$tables$'Grand mean'
onedf <- ypred * ypred / (2 * mu) # this is a vector
mt.onedf <- lm(logy ~ day + store + onedf)
anova(mt.onedf) # eta is not significant, the transformed model works

# step 7: do analysis on the transformed model
mt <- lm(logy ~ day + store)
par(mfrow = c(2, 2))
plot(mt)
par(mfrow = c(1, 1))

TukeyHSD(mt.aov, "store")
sort(model.tables(mt.aov, type = "means")$tables$store)
exp(sort(model.tables(mt.aov, type = "means")$tables$store))


##### Lecture Code Chapter 10: Unbalanced Data #####
rm(list = ls())
unbalanced <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/weight-gain-unbalanced.txt",
                         header = T)
names(unbalanced)
unbalanced <- within(unbalanced, {
  a <- as.factor(a)
  b <- as.factor(b)
})
attach(unbalanced)

### (1) SAS Type I: sequential sum of squares
m0 <- lm(y ~ 1)
anova(m0)

m1 <- lm(y ~ a)
anova(m1)

m2 <- lm(y ~ b)
anova(m2)

m3 <- lm(y ~ a + b)
anova(m3)

m4 <- lm(y ~ b + a)
anova(m4)

m5 <- lm(y ~ a * b)
anova(m5)

m6 <- lm(y ~ b * a)
anova(m6)

### (2) SAS Type II
## using lm and anova
# A
anova(lm(y ~ b + a))
# B
anova(lm(y ~ a + b))
# A*B
anova(lm(y ~ a * b))

## using Anova function in car package
library(car)
m <- lm(y ~ a * b)
Anova(m, type = "II")


### (3) SAS Type III
Anova(m, type = "III") # don't use this one

# check the parametrization now using
options("contrasts")
# change the parametrization and re-fit the model
options(contrasts = c('contr.sum', 'contr.poly'))
m2.aov <- aov(y ~ a * b)
Anova(m2.aov, type = "III") # use this one instead
Anova(lm(y ~ a * b), type = "III") # lm has the same effect as aov here


### (4) Contrasts with unbalance data
m.aov <- aov(y ~ a * b)

# use gmodels package
library(gmodels)
library(cfcdae)
# when using gmodels to do this, ensure the model matrix is right
options(contrasts = c('contr.treatment', 'contr.poly'))
m.aov <- aov(y ~ a * b)
fit.contrast(m.aov, "a", coeff = c(1, -1))
linear.contrast(m.aov, term = a, contr.coefs = c(1, -1))

options(contrasts = c('contr.sum', 'contr.poly'))
fit.contrast(m.aov, "a", coeff = c(1, -1))

# use cfcdae package
# note that cfcdae is not in CRAN, so if we are working for companies or 
# publishing papers, then use gmodels

linear.contrast(m.aov, term = a, contr.coefs = c(1, -1))
# contr.sum or contr.treatment would not influence linear.contrast()

### (5) Estimating marginal means
# use emmeans package
library(emmeans)
m.lm <- lm(y ~ a * b)
emmeans(m.lm, specs = "a")
emmeans(m.lm, specs = ~ a) # these two ways are the same, prefer this one

emmeans(m.lm, specs = "b")
emm <- print(emmeans(m.lm, specs = ~ a + b)) 
# the type of emmeans(m.lm, specs = ~ a + b) is S4, adding print() will 
# change it to a dataframe
emm

mean(emm$emmean[emm$a == 1])
mean(emm$emmean[emm$a == 2])

detach(unbalanced)

### (6) Interaction plots
# use emmeans package
library(emmeans)
rm(list = ls())
pwg <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/pig-weight-gain.txt",
                  header = T)
names(pwg)
pwg <- within(pwg, {
  L <- as.factor(L)
  P <- as.factor(P)
  S <- as.factor(S)
  blocks <- as.factor(blocks)  
})

m <- lm(y ~ blocks + L * P * S, data = pwg)
emmip(m, P ~ L)
emmip(m, P ~ L | S)

detach(pwg)

### tests to check the relationships between three types of SS
dat <- read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr10.5",
                  header = T)
# since the unbalance of this dataset has some structure, modifying some of
# them would help better check the relationship between these three types
dat <- dat[c(-19, -27),] # delete 4.51 and 3.00 to make it more unbalanced

attach(dat)
variety <- as.factor(variety)
location <- as.factor(location)
temperature <- as.factor(temperature)

### Type I SS
m <- lm(size ~ variety * location * temperature)
anova(m)

# A
anova(lm(size ~ variety))

# B
anova(lm(size ~ variety + location))

# C
anova(lm(size ~ variety + location + temperature))

# AB
anova(lm(size ~ variety + location + temperature + variety*location))

# AC
anova(lm(size ~ variety + location + temperature + variety*location + 
           variety * temperature))

# BC
anova(lm(size ~ variety + location + temperature + variety*location + 
           variety * temperature + location * temperature))

### Type II SS
library(car)
Anova(m, type = "II")

# AB
anova(lm(size ~ variety + location + temperature + variety * temperature + 
           location * temperature+ variety*location ))



### Type III SS
Anova(m, type = "III")

detach(dat)

### 2^k factorial design
rm(list = ls())
dat2k <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-single-rep.txt", 
                    header = T)
dat2k <- within(dat2k, {
  a <- as.factor(a)
  b <- as.factor(b)
  c <- as.factor(c)
  d <- as.factor(d)
})
attach(dat2k)
summary(dat2k)

## ANOVA
options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ a * b * c * d)
# dfE is 0, so we cannot get the diagnostic plots
model.matrix(m)
anova(m)
coef(m) # actually estimating the effects
# check the meaning of coefficents
m.aov <- aov(y ~ a * b * c * d)
fit.contrast(m.aov, "a", c(1, -1))
model.tables(m.aov, type = "means")$table$a

## Identify important effects using normal plot
qqnorm(coef(m)[-1], main = "")  # remove the effect of the intercept
title(main = "Normal Probability Plot for Effects \n y~a*b*c*d")

library(faraway)  # with the number on the plot
qqnorml(coef(m)[-1], main = "")
title(main = "Normal Probability Plot for Effects qqnorml() \n y~a*b*c*d")

## (1) half-normal plot
halfnorm(coef(m)[-1], nlab = 3)
title(main = "Half-normal Plot for Effects \n y~a*b*c*d")

halfnorm(coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "Half-normal Plot for Effects \n y~a*b*c*d")

## (2) Lenth's pseudo standard error for "p-values"
# Step 1: Define the "pseudo random error"
z <- 2 * coef(m)[-1]
(s0 <- median(abs(z)) * 1.5)
(z2 <- z[abs(z) < s0 *2.5])
(pse <- median(abs(z2)) * 1.5)

# Step 2: Use the "pseudo random error" for t-tests
(t0 <- z / pse)
2 * pt(-abs(t0), 15/3)

## (3) Pooling higher order interaction effects into the error
m2 <- lm(y ~ (a + b + c + d)^2)
anova(m2)

par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))

## (3) Collapse into 2-factorial design
m3 <- lm(y ~ b * c, data = dat2k)
par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))
library(MASS)
boxcox(m3)

m4 <- lm(y ~ b + c, data = dat2k)
par(mfrow = c(2, 2))
plot(m4)
par(mfrow = c(1, 1))
boxcox(m4)

interaction.plot(x.factor = b, trace.factor = c, response = y, 
                 fun = mean, legend = T, xpd = NA, type = 'b', 
                 trace.label = 'Factor C', xlab = "Factor B", ylab = "Mean Response", 
                 lty = c(1, 2, 4), col = c("black", "red", "blue"), 
                 main = "Interaction Plot for y ~ b * c")

TukeyHSD(aov(m3), "b:c")

detach(dat2k)

##### Lecture Code Chapter 13: Complete Block Designs #####
library(BHH2)
data("penicillin.data")
names(penicillin.data)
attach(penicillin.data)

m <- lm(yield ~ blend + treat)
anova(m)

detach(penicillin.data)


### Latin Square Design
rm(list = ls())
rabbits <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/rabbits1-latin-sq.txt",
                      header = T)
rabbits2 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/rabbits2-latin-sq.txt",
                      header = T)
names(rabbits)
rabbits$treat <- rabbits2$treat[1:16]
attach(rabbits)
day <- as.factor(day)
subject <- as.factor(subject)
treat <- as.factor(treat)

m <- lm(y ~ day + subject + treat)
anova(m)
par(mfrow = c(2, 2))
plot(m)  # not good
par(mfrow = c(1, 1))
boxcox(m)

# Tukey's 1df test
m.aov <- aov(y ~ day + subject + treat)
ypred <- m.aov$fitted.values
(mu <- model.tables(m.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m.onedf <- lm(y ~ day + subject + treat + onedf)
anova(m.onedf)
eta <- m.onedf$coefficients["onedf"]
1 - eta

# try square transformation
y2 <- y*y
m2 <- lm(y2 ~ day + subject + treat)
anova(m2)
par(mfrow = c(2, 2))
plot(m2)  # still not good enough
par(mfrow = c(1, 1))
boxcox(m2)

# Multiple comparisons
m2.aov <- aov(y2 ~ day + subject + treat)
TukeyHSD(m2.aov, "treat")
sort(model.tables(m2.aov, type = "means")$tables$treat)

detach(rabbits)
detach(rabbits2)

### Replicated Latin Squares
rm(list = ls())
rabbits2 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/rabbits2-latin-sq.txt",
                       header = T)
names(rabbits2)
attach(rabbits2)
day <- as.factor(day)
order <- as.factor(order)
subject <- as.factor(subject)
treat <- as.factor(treat)
reps <- as.factor(reps)

# ANOVA
m <- lm(y ~ reps + day + subject + treat)
anova(m)

par(mfrow = c(2, 2))
plot(m)  # still not good enough
par(mfrow = c(1, 1))
boxcox(m)

# Tukey's 1df test
m.aov <- aov(y ~ reps + day + subject + treat)
ypred <- m.aov$fitted.values
(mu <- model.tables(m.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m.onedf <- lm(y ~ day + subject + treat + onedf)
anova(m.onedf)
eta <- m.onedf$coefficients["onedf"]
1 - eta

# do square transformation
y2 <- y*y
m2 <- lm(y2 ~ reps + day + subject + treat)
m2.aov <- aov(m2)
ypred <- m2.aov$fitted.values
(mu <- model.tables(m2.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m2.onedf <- lm(y2 ~ day + subject + treat + onedf)
anova(m2.onedf)
eta <- m2.onedf$coefficients["onedf"]
1 - eta

par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2.aov)

summary(m2.aov)
TukeyHSD(m2.aov, "treat")
# balanced design, so we can use model.tables
sort(model.tables(m2.aov, type="means")$tables$treat)

detach(rabbits2)


##### Lecture Code Chapter 14: Incomplete Balanced Designs #####
bibd <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/BIBD-cloths.txt",
                   header = T)
names(bibd)
bibd <- within(bibd, {
  blocks <- as.factor(blocks)
  cloths <- as.factor(cloths)
})
with(bibd, {table(blocks, cloths)})


m <- lm(y ~ blocks + cloths, data = bibd)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
m.aov <- aov(m)
summary(m.aov)

options(contrasts = c('contr.sum', 'contr.poly'))
library(gmodels)
fit.contrast(m.aov, "cloths", c(-1, 1, 0, 0, 0))

## Don't use TukeyHSD, the result of TukeyHSD is different from the other two.
TukeyHSD(m.aov, "cloths")

## instead, use the following two
library(multcomp)
g <- glht(m.aov, linfct=mcp(cloths="Tukey"), alternative = "two.sided")
summary(g)
plot(g)

# library(cfcdae)
# linear.contrast(m, term = cloths, allpairs = TRUE, confidence = 0.95)


### Interblock recovery 
attach(bibd)
## Step 1: Regressio with response y.j
y.inter <- tapply(y, blocks, sum)
(x <- table(blocks, cloths))
x1 <- x[, 1]
x2 <- x[, 2]
x3 <- x[, 3]
x4 <- x[, 4]
x5 <- x[, 5]

m.inter <- lm(y.inter ~ x1 + x2 + x3 + x4 + x5 - 1)
summary(m.inter)

## Step 2: Inter-block estimate for the mu1 and the contrast
(mu.i.inter <- coef(m.inter))
(w.inter.est <- coef(m.inter)["x1"] - coef(m.inter)["x2"])
names(w.inter.est) <- NULL  # get rid of the x1 name

## Step 3: Inter-block estimate for the variance of the contrast
b <- 10
N <- 30
g <- 5
k <- 3
r <- 6
lambda <- 3

m2 <- lm(y ~ cloths + blocks)
anova(m2)
(mse.m2 <- anova(m2)$'Mean Sq')
(s2b <- (b - 1) / (N - g) * (mse.m2[2] - mse.m2[3]))

# Inter-block estimate of the variance of the contrast w
(w.inter.var <- (k ^ 2 * s2b + k * mse.m2[3]) * 2 / (r - lambda)) 

## Step 4: Combined estimate for the contrast
# Intra-block estimate of contrast and its variance
library(gmodels)
a <- fit.contrast(m.aov, cloths, c(1, -1, 0, 0, 0))
(w.intra.est <- a[1, "Estimate"])
(w.intra.var <- a[1, "Std. Error"]^2)

# combined estimate:
c1 <- 1 / w.intra.var
c2 <- 1 / w.inter.var
(w.comb.est <- (c1 * w.intra.est + c2 * w.inter.est) / (c1 + c2))

# variance of the combined estimate
(w.comb.var <- (c1^2 * w.intra.var + c2^2 * w.inter.var) / (c1+c2)^2)

## Step 5: 95% CI for the contrast
(df.E <- anova(m2)$Df[2])
hw <- qt(0.975, df.E) * sqrt(w.comb.var)
(ci <- c("lower"=w.comb.est - hw, "upper" = w.comb.est + hw))

detach(bibd)

##### Lecture Code Chapter 15: Confounding #####
## 2 blocks
rm(list = ls())
d2k.2bl <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-confounding-2blocks.txt",
                      header = T)
d2k.2bl$a <- as.factor(d2k.2bl$a)
d2k.2bl$b <- as.factor(d2k.2bl$b)
d2k.2bl$c <- as.factor(d2k.2bl$c)
d2k.2bl$blocks <- as.factor(d2k.2bl$blocks)
summary(d2k.2bl)
attach(d2k.2bl)

# ANOVA
options(contrasts = c('contr.sum', 'contr.poly'))
m <- lm(y ~ blocks + a * b * c)
anova(m)
alias(m)

# use half-normal plot to determine the significant effects
library(faraway)
halfnorm(-2*coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "Half-normal plot, effects in y~blocks+a*b*c")

detach(d2k.2bl)

## 4 blocks
rm(list = ls())
d2k.4bl <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-confounding-4blocks.txt",
                      header = T)
d2k.4bl$a <- as.factor(d2k.4bl$a)
d2k.4bl$b <- as.factor(d2k.4bl$b)
d2k.4bl$c <- as.factor(d2k.4bl$c)
d2k.4bl$d <- as.factor(d2k.4bl$d)
d2k.4bl$blocks <- as.factor(d2k.4bl$blocks)
summary(d2k.4bl)

attach(d2k.4bl)
# ANOVA
options(contrasts = c('contr.sum', 'contr.poly'))
m <- lm(y ~ blocks + a*b*c*d)
anova(m)
alias(m)

# determine important effects
library(faraway)
halfnorm(-2*coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "Half-normal plot, effects in y~blocks+a*b*c*d")

detach(d2k.4bl)

## Construct a design
library(conf.design)
G <- rbind(c(1, 1, 0, 1), c(1, 0, 1, 1))
colnames(G) <- LETTERS[1:4]
G
conf.set(G, p=2)
conf.design(G, p = 2, block.name = 'blocks')


## Partial confounding, when we have replicates
rm(list = ls())
d2k.p <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-partial-confounding.txt",
                    header = T)
d2k.p$a <- as.factor(d2k.p$a)
d2k.p$b <- as.factor(d2k.p$b)
d2k.p$c <- as.factor(d2k.p$c)
d2k.p$blocks <- as.factor(d2k.p$blocks)
summary(d2k.p)

attach(d2k.p)
options(contrasts = c('contr.sum', 'contr.poly'))
m <- lm(y ~ blocks + a*b*c)
anova(m)
alias(m)

library(faraway)
halfnorm(-2*coef(m)[-1], labs = names(coef(m)[-1]), nlab=3)
title(main = 'Half-normal plot, effects in y~blocks+a*b*c')

detach(d2k.p)

##### Lecture Code Chapter 18: Fractional factorials #####
rm(list = ls())
BHH <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-BHH-household-liquids.txt", 
                  header = T)
BHH <- within(BHH, {
  a <- as.factor(a)
  b <- as.factor(b)
  c <- as.factor(c)
  d <- as.factor(d)
})

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ a * b * c * d, data = BHH)
anova(m)
library(faraway)
halfnorm(-2 * coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "BHH household liquids\n effects in y ~ a*b*c*d")
coef(m)
alias(m)

m2 <- lm(y ~ an + bn, data = BHH)
coef(m2)


## Example: 2^(5-1)
rm(list = ls())
d2k5.1 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-2-frac-factorial.txt", 
                     header = T)
d2k5.1 <- within(d2k5.1, {
  a <- as.factor(a)
  b <- as.factor(b)
  c <- as.factor(c)
  d <- as.factor(d)
  e <- as.factor(e)
})

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ a * b * c * d * e, data = d2k5.1)
anova(m)

library(faraway)
halfnorm(-2 * coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "Effects in y ~ a*b*c*d*e")
coef(m)
alias(m)


## Example: 2^(7-4)
rm(list = ls())
d2k7.4 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/2k-BHH-filtration-cycle-foldover.txt", 
                     header = T)
d2k7.4 <- within(d2k7.4, {
  a <- as.factor(a)
  b <- as.factor(b)
  c <- as.factor(c)
  d <- as.factor(d)
  e <- as.factor(e)
  f <- as.factor(f)
  g <- as.factor(g)
  blocks <- as.factor(blocks)
})
summary(d2k7.4)

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ blocks + a * b * c * d * e * f * g, data = d2k7.4)
anova(m)
library(faraway)
halfnorm(-2 * coef(m)[-1], labs = names(coef(m)[-1]), nlab = 3)
title(main = "BHH household liquids\n effects in y ~ a*b*c*d*e*f*g")


library(FrF2)
aliases(m)

# obtain the defining relation
dr <- alias(m)$Complete[, "(Intercept)"]
dr[dr %in% c(1, -1)]


## use R to construct the first and second fractions of the design
library(FrF2)
# note that the results can be random, we will get different results for different runs
(design.frf2 <- FrF2(nruns = 8, nfactors = 7, generators = c("AB", "AC", "BC", "ABC")))
FrF2(nruns = 8, nfactors = 7, generators = c("-AB", "-AC", "-BC", "ABC"))

design.info(design.frf2)$aliased$main


## construct 2^(k-q) fractional factorials with blocks
library(FrF2)
FrF2(nruns = 16, nfactors = 7, blocks = 2)


## construct 2^(k-q) fractional factorials with a given resolution
FrF2(nfactors = 7, resolution = 4)


##### Lecture Code Chapter 17: Designs with covariates #####
rm(list = ls())
library(ISwR)
data("hellung")
# factorize glucose and use yes and no as labels
hellung$glucose <- factor(hellung$glucose, labels = c('yes', 'no'))
attach(hellung)
summary(hellung)

# graphic overview
library(gplots)
plot(conc, diameter, pch = as.numeric(glucose), col = ifelse(glucose=="yes", 'red', 'blue'), 
     xlab = "Concentration", ylab = "Diameter")
legend("topright", legend = c('glucose', 'no glucose'), pch = 1:2, col = c('red', 'blue'))
color = c('red', 'blue')
label = c('yes', 'no')
for (i in 1:2){
  lines(lowess(diameter~conc, data = hellung[glucose == label[i], ]), col = color[i])
}
title(main = "Tetrahymena cells, diameter vs cell concentration \n with lowess lines")

# do log10 transformation
plot(log10(conc), log10(diameter), pch = as.numeric(glucose), col = ifelse(glucose=="yes", 'red', 'blue'), 
     xlab = "Concentration", ylab = "Diameter")
legend("topright", legend = c('glucose', 'no glucose'), pch = 1:2, col = c('red', 'blue'))
color = c('red', 'blue')
label = c('yes', 'no')
for (i in 1:2){
  lines(lowess(log10(diameter) ~ log10(conc), data = hellung[glucose == label[i], ]), col = color[i])
}
title(main = "Tetrahymena cells, diameter vs cell concentration \n with lowess lines")

options(contrasts = c("contr.treatment", "contr.poly"))

# ANOVA
m1 <- lm(log10(diameter) ~ log10(conc) * glucose)
library(car)
Anova(m1, type = 2)

par(mfrow = c(2, 2))
plot(m1)
par(mfrow = c(1, 1))
boxcox(m1)
summary(m1)

# Classical ANCOVA
m2 <- lm(log10(diameter) ~ log10(conc) + glucose)
Anova(m2, type = 2)
summary(m2)$coefficients

# CI for the difference between glucose groups
library(gmodels)
fit.contrast(m2, glucose, coeff = c(-1, 1), conf.int = 0.95)

# covariate-adjusted means
logconc.c <- log10(conc) - mean(log10(conc))
m3 <- lm(log10(diameter) ~ logconc.c + glucose)
Anova(m3)
summary(m3)
# If I run this part of code after options, then the result will be glucose1. If not, then it would be 
# glucoseno. So, which one should I use? Should I set options in the first place or not?

library(emmeans)
emmeans(m3, ~glucose)

# analysis without adjusting for the covariates
m4 <- lm(log10(diameter) ~ glucose)
anova(m4)

detach(hellung)


##### Lecture Code Chapter 11: Random Effects #####
rm(list = ls())
leaves <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/calcium-in-leaves-1.txt", 
                     header = T)
names(leaves)
leaves$leaf <- as.factor(leaves$leaf)
attach(leaves)
summary(leaves)

options(contrasts = c('contr.sum', 'contr.poly'))
m <- lm(y ~ leaf)
anova(m)
coef(m)
detach(leaves)


##### Lecture Code Chapter 12: Nested Designs #####
rm(list = ls())
nb <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/nested-batches.txt",
                 header = T)
names(nb)
nb$supplier <- as.factor(nb$supplier)
nb$batch <- as.factor(nb$batch)
summary(nb)

options(contrasts = c('contr.sum', 'contr.poly'))

## ANOVA for nested design, and both A and B are fixed
m <- aov(y ~ supplier/batch, nb)
summary(m)

## ANOVA for cross-over design, and both A and B are fixed
summary(aov(y ~ supplier * batch, nb))

## ANOVA for nested design, and B (batch) is random, A(supplier) if fixed
# F-test for A
m1 <- aov(y ~ supplier + Error(batch %in% supplier), nb)
summary(m1)

# F-test for B
m2 <- aov(y ~ supplier + batch %in% supplier, nb)
summary(m2)

# or use linear mixed model
library(nlme)
nb$sb <- with(nb, as.factor(supplier:batch))
m.lme <- lme(y ~ supplier, data = nb, random = ~1|sb)
summary(m.lme)
# F-test for A
anova(m.lme)

# LRT test for B
m.lm <- lm(y ~ supplier, data = nb)
anova(m.lm)
library(RLRsim)
exactLRT(m.lme, m.lm)

library(gmodels)
fit.contrast(m.lme, "supplier", coeff = c(1, -1, 0), conf.int = 0.95)
intervals(m.lme)
fixef(m.lme)
ranef(m.lme)
fitted(m.lme)
coef(m.lme)
predict(m.lme, level = 1)
predict(m.lme, level = 0)


## when both A and B are random, use linear mixed model
library(nlme)
m.lme <- lme(y ~ 1, data = nb, random = ~1|supplier/batch)
summary(m.lme)

# ANOVA for fixed effect
anova(m.lme) # testing wheter mu is significant, and we don't care about it
fixef(m.lme)

# Inference for random effects
intervals(m.lme)
ranef(m.lme)


### Hasse Diagrams
rm(list = ls())
glucose <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/glucose-crossed-nested.txt",
                      header = T)
glucose$conc <- as.factor(glucose$conc)
glucose$day <- as.factor(glucose$day)
glucose$run <- as.factor(glucose$run)
summary(glucose)

# ANOVA test for AC and C (using MSE as denominator)
summary(aov(y ~ conc*(day/run), glucose))

# F-test for A (AB as denominator)
summary(aov(y ~ conc + day + Error(conc:day), glucose))

# F-test for B (C as denominator)
summary(aov(y ~ day + Error(day:run), glucose))

# F-test for AB (AC as denominator)
summary(aov(y ~ conc*day + day/run + Error(conc:day:run), glucose))

# diagnostics
m <- lm(y ~ conc*(day/run), glucose)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))

##### Lecture Code Chapter R: Repeated measures and longitudinal data #####
library(nlme)
data("Orthodont")
head(Orthodont)
summary(Orthodont)
data.class(Orthodont)
formula(Orthodont)

## Recreate groupedData object
Ortho <- groupedData(distance ~ age | Subject, 
                     data = as.data.frame(Orthodont), 
                     outer = ~ Sex, inner = ~ age, 
                     labels = list(x="Age", y="Distance"), 
                     units = list(x="(in years)", y="(in mm)"))
summary(Ortho)  
data.class(Ortho)
formula(Ortho)
for(i in 1:4) {print(c(names(Ortho)[i], data.class(Ortho[, i])))}

par(mfrow = c(1, 1))
plot(Ortho, main = "Orthonut, by Subject")


## Fitting longitudinal model with lme()
options(contrasts = c("contr.sum", "contr.poly"))
m.lme <- lme(distance ~ Sex + age, data = Ortho, random = ~ 1 | Subject)
anova(m.lme)
fixef(m.lme)
ranef(m.lme)
intervals(m.lme)$fixed

# fitted values
library(lattice)
xyplot(fitted(m.lme) ~ age | Subject, data = Ortho, type = "b")

# observations + fitted
plot(augPred(m.lme), layout = c(5, 6), main = "Predicted, per Subject")


## Fit centered age variable for easier interpretation
m2.lme <- update(m.lme, fixed = distance ~ Sex + I(age - 11))
anova(m2.lme)
fixef(m2.lme)

## Or treat age = 8 to be time = 0
Ortho$time <- Ortho$age - 8
summary(Ortho$time)
m3.lme <- update(m.lme, fixed = distance ~ Sex + time)
anova(m3.lme)
fixef(m3.lme)


## Residual plot
plot(m.lme)
plot(m.lme, id = 0.1)

xyplot(resid(m.lme) ~ age | Subject, data = Ortho, layout = c(5, 6),
       panel = function(x, y, ...){
         panel.xyplot(x, y, ...)
         panel.lines(range(x), c(0, 0))
       })

qqnorm(m.lme)
qqnorm(m.lme, form = ~ranef(.), id = 0.1)


## Should there be different slop for age by Sex?
m4.lme <- update(m.lme, fixed = distance ~ Sex * age)
anova(m4.lme)
fixef(m4.lme)
summary(Ortho$Sex)
plot(augPred(m4.lme), layout = c(5, 6), main = "Predicted, per Subject, m4.lme")


## Should there be different slope for age by Subject?
m5.lme <- lme(distance ~ Sex * age, data = Ortho, 
              random = ~ age | Subject)
anova(m5.lme)

# test for the random slope effect
m6.lme <- update(m5.lme, random = ~ 1 | Subject)
anova(m6.lme, m5.lme)

# comparison of random effects for REML and ML
m5.lmeM <- update(m5.lme, method = 'ML')
anova(m5.lmeM)

plot(compareFits(ranef(m5.lme), ranef(m5.lmeM)), mark = c(0, 0),
     main = 'Compare m5.lme, m5.lmeM')

# Correlation structure for within-group error
m5.lme$call
library(lattice)
xyplot(resid(m5.lme) ~ age | Subject, data = Ortho, layout = c(5, 6),
       panel = function(x, y){
         panel.xyplot(x, y)
         panel.lines(range(x), c(0, 0))
       })

m7.lme <- lme(distance ~ Sex*age, data = Ortho, random = ~ age | Subject, 
              correlation = corAR1())
anova(m5.lme, m7.lme)
plot(augPred(m7.lme), layout = c(5, 6), main = "Predicted, per Subject, m4.lme")

# Investigate the covariance matrices
library(mgcv)
extract.lme.cov2(m7.lme, data = Ortho)$V[[1]]
getVarCov(m5.lme, type ='random.effects')
getVarCov(m5.lme, individual = 'M01', type ='conditional')
getVarCov(m5.lme, individual = 'M01', type ='marginal')
intervals(m5.lme)

# consider quadratic  polynomial fit
mq.lme <- lme(distance ~ Sex * (I(age - 11) + I((age - 11)^2)), 
              data = Ortho, random = ~ 1|Subject)
anova(mq.lme)

m1.lmeM <- lme(distance ~ Sex * I(age - 11), data = Ortho, 
               random = ~ age|Subject, method = "ML")
mq.lmeM <- lme(distance ~ Sex * (I(age - 11) + I((age - 11)^2)), 
               data = Ortho, random = ~ age + I(age^2)|Subject, method = "ML")
anova(m1.lmeM, mq.lmeM)


##### Lecture Code Chapter P: Repeated measures with nested grouping factors #####
library(nlme)
data("Pixel")
head(Pixel)
data.class(Pixel)
formula(Pixel)

library(lattice)
plot(Pixel, main = "Pixel density, by Subject")
xyplot(pixel ~ day | Dog, data = Pixel, groups = ~Side, type = c('g', 'p', 'l'), 
       auto.key = T, main = 'Pixel density, by Subject, by Side', layout = c(5, 2))

m1s.lme <- lme(pixel ~ day + I(day^2) + I(day^3) + Side, data = Pixel, 
               random = ~ 1 | Dog/Side)
anova(m1s.lme)

# Test for systematic effect of Side with LRT 
m1s.lmeM <- update(m1s.lme, method = 'ML')   # (differs in fixed effect)
m1.lmeM <- update(m1s.lmeM, fixed = pixel ~ day + I(day^2) + I(day^3))
anova(m1s.lmeM, m1.lmeM)

# Structure of random effect
m1.lme <- lme(pixel ~ day + I(day^2) + I(day^3), data = Pixel, random = ~1|Dog/Side)
m1nrs.lme <- update(m1.lme, random = ~1|Dog)
anova(m1.lme, m1nrs.lme) # although p-val is small, still keep the nested random effect

library(lattice)
plot(augPred(m1.lme), main = 'Pixel density, m1.lme', layout = c(6, 4))

# Investigate the covariance matrices
library(mgcv)
extract.lme.cov2(m1.lme, data = Pixel)$V[[8]]
index <- with(data = Pixel, {which(Dog==8)})
Pixel[index, ]

# Test for Dog-specific slope for day
m1.lme$call
m2.lme <- lme(pixel ~ day + I(day^2) + I(day^3), data = Pixel, 
          random = list(Dog = ~1, Side = ~1))
m3.lme <- update(m2.lme, random = list(Dog = ~day, Side = ~1))
anova(m1.lme, m2.lme, m3.lme)  # choose m3.lme

m4.lme <- update(m2.lme, random = list(Dog = ~day, Side = ~day))

library(lattice)
plot(augPred(m3.lme), main = "Pixel density, m3.lme", layout = c(6, 4))

# consider AR(1) correlation structure
m5.lme <- update(m3.lme, correlation = corCAR1())
anova(m3.lme, m5.lme)
# always remember to check the graphic result
plot(augPred(m5.lme), main = "Pixel density, m5.lme", layout = c(6, 4))
plot(compareFits(ranef(m3.lme, level = 1), ranef(m5.lme, level = 1)),
     mark = c(0, 0), main = "Random effects, m3.lme and m5.lme")

# variance proportion
intervals(m3.lme)

### Inference using lme4 library
library(lme4)
options(contrasts = c('contr.sum', 'contr.poly'))
data("Orthodont", package = 'nlme')
head(Orthodont)
m.lmer <- lmer(distance ~ Sex * age + (1|Subject), data = Orthodont)
summary(m.lmer)

# ANOVA for fixed effects
anova(m.lmer)
fixef(m.lmer)
ranef(m.lmer)

# Tests for fixed effects
# (1) LRT
m1 <- lmer(distance ~ Sex + age + (1|Subject), data = Orthodont)
m2 <- lmer(distance ~ Sex * age + (1|Subject), data = Orthodont)
anova(m1,m2)

# (2) F-test
library(pbkrtest)
pbkrtest::KRmodcomp(m1, m2)  # Kenward-Rogers

# (3) p-value calculated with anova{lmerTest}
library(lmerTest)
m2.reml <- lmerTest::lmer(distance ~ Sex * age + (1|Subject), data = Orthodont, 
                          REML = T)
anova(m2.reml, type = 2)

# (4) p-value calculated with parametric bootstrap
boot <- pbkrtest::PBmodcomp(m2, m1) # Parametric bootstrap
summary(boot)

# LRT for Subject-specific slope
m3 <- lmer(distance ~ Sex * age + (1|Subject), data = Orthodont)
m4 <- lmer(distance ~ Sex * age + (1|Subject) + (0 + age | Subject), data = Orthodont)
m5 <- lmer(distance ~ Sex * age + (age|Subject), data = Orthodont)
anova(m3, m4)
anova(m3, m5)



##### Lecture Code Chapter 16: Split Plots #####
library(nlme) # data is contained in this library. We don't need to use it
              # for the purpose of analysis
data(Oats)
head(Oats)
Oats$nitro <- as.factor(Oats$nitro)
oats4 <- Oats[Oats$Block %in% list("I", "II", "III", "IV"),]
summary(oats4)

## ANOVA for split-plot design
options(contrasts = c('contr.sum', 'contr.poly'))
m <- aov(yield ~ Block + Variety + Error(Block:Variety) + nitro + 
           nitro:Variety, data = oats4)
summary(m)

## A closer look at the Whole Plot analysis
oats4.gd <- groupedData(yield ~ Variety | Block/nitro, data = oats4, 
                        FUN = mean)
summary(oats4.gd)
# data are collapsed over all levels below "Block", need to specify nlme
oats4.wp <- nlme::collapse(oats4.gd, collapseLevel = "Block") 
head(oats4.wp)
summary(oats4.wp)

# ANOVA for the whole plot analysis
mwp <- lm(yield ~ Block + Variety, data = oats4.wp)
anova(mwp)
library(gmodels)
fit.contrast(mwp, "Variety", coeff = c(1, -1, 0), conf.int = 0.95)

## A closer look at the Split Plot analysis
msp <- lm(yield ~ Block:Variety + nitro + nitro:Variety, data = oats4)
anova(msp)     
fit.contrast(msp, "nitro", coeff = c(1, 0, 0, -1), conf.int = 0.95)
par(mfrow = c(2, 2))
plot(msp)
par(mfrow = c(1, 1))

## Use orthogonal polynomials to check linear, quadratic, or cubic
options(contrasts=c("contr.treatment", "contr.poly"))
oats4$nitro <- as.ordered(oats4$nitro)
    # for ordered factors, R uses orthogonal contrasts by default
m <- lm(yield ~ Block:Variety + nitro + nitro:Variety, data = oats4)
summary(m)

## Use linear regression
oats4$nitro <- as.numeric(oats4$nitro)
m2 <- lm(yield ~ Block:Variety + nitro + I(nitro^2), data = oats4)
m1 <- lm(yield ~ Block:Variety + nitro, data = oats4)
m0 <- lm(yield ~ Block:Variety, data = oats4)
anova(m0, m1, m2)
summary(m2)

## What happens if we ignore the split plot structure?
anova(lm(yield ~ (Block + Variety + nitro)^2, data = oats4))

anova(lm(yield ~ Block + Variety * nitro, data = oats4))

### Example: Oats data(cont.)
library(nlme)
data(Oats)
head(Oats)
Oats$nitro <- as.ordered(Oats$nitro)
oats4 <- Oats[Oats$Block %in% list("I", "II", "III", "IV"),]
summary(oats4)

m1.lme <- lme(yield ~ Variety * nitro, data = oats4, random = ~ 1|Block/Variety)
anova(m1.lme)
summary(m1.lme)

lmer.mod <- lmer(yield ~ Variety * nitro + (1|Block) + (1|Block:Variety), oats4)
sumary(lmer.mod)
rmod <- lmer(yield ~ Variety + nitro + (1|Block) + (1|Block:Variety), oats4)
KRmodcomp(lmer.mod, rmod)
rmod <- lmer(yield ~ nitro + (1|Block) + (1|Block:Variety), oats4)
KRmodcomp(lmer.mod, rmod)


# confidence intervals for split plots
# options(contrasts = c('contr.treatment', 'contr.poly'))
m.lme <- lme(yield ~ Variety * nitro, data = oats4, random = ~1|Block/Variety)
anova(m.lme)

library(gmodels)
fit.contrast(m.lme, varname = "Variety", coeff = c(1, -1, 0))
fit.contrast(m.lme, varname = "nitro", coeff = c(1, 0, 0, -1))


## Use separate lm()
# Step1: collapse the data
oats4.gd <- groupedData(yield ~ Variety | Block/nitro, data = oats4, FUN = mean)
oats4.wp <- collapse(oats4.gd, collapseLevel = "Block")
# Step2: whole plot model
m.wp <- lm(yield ~ Block + Variety, data = oats4.wp)
anova(m.wp)
fit.contrast(m.wp, varname="Variety", coeff = c(1, -1, 0))
# Step3: split plot model
m.sp <- lm(yield ~ Block:Variety + nitro + nitro:Variety, data = oats4)
anova(m.sp)
fit.contrast(m.sp, varname = "nitro", coeff = c(1, 0, 0, -1))  # should under contr.sum

aov.sp <- aov(yield ~ Block:Variety + nitro + nitro:Variety, data = oats4)
model.tables(aov.sp, "means")


options(contrasts = c("contr.treatment", "contr.poly"))  # not correct
msp <- lm(yield ~ Block:Variety + nitro + nitro:Variety, data = oats4)
anova(msp)
library(gmodels)
fit.contrast(msp, varname = "nitro", coeff = c(1, 0, 0, -1))



##### Lecture Code Faraway Chapter10 & 13: GEE & GLMM #####
library(faraway)
data(ctsib)
ctsib$stable <- ifelse(ctsib$CTSIB == 1, 1, 0)
head(ctsib)

### GEE
options(contrasts = c('contr.treatment', 'contr.poly'))
library(geepack)
g1 <- geeglm(stable ~ Sex + Age + Height + Weight + Surface + Vision, 
             id = Subject, family = binomial, data = ctsib, 
             corstr = "exchangeable", scale.fix = TRUE)
summary(g1)

# comparing nested models
g2 <- update(g1, stable ~ Sex + Height + Surface + Vision) # delete Age and Weight
anova(g1, g2)

# linear hypothesis
library(doBy)
cm <- rbind("dome vs open" = c(0, 0, 0, 0, 0, 0, 1, -1), 
            "open vs closed" = c(0, 0, 0, 0, 0, 0, 0, 1))
esticon(g1, L = cm, beta0 = 0, conf.int = TRUE, level = 0.95)


### GLMM
library(lme4)
g <- glmer(stable ~ Sex + Age + Height + Weight + Surface + Vision + 
             (1|Subject), family = binomial, data = ctsib)

library(MASS)
m <- glmmPQL(stable ~ Sex + Age + Height + Weight + Surface + Vision, 
             random = ~ 1 | Subject, family = binomial, data = ctsib)
summary(m)

### GEE vs GLMM
data("toenail")
head(toenail)
m.gee <- geeglm(outcome ~ treatment + month, id = ID, 
                family = binomial, data = toenail, corstr = "exchangeable")
summary(m.gee)

m.glmm <- glmer(outcome ~ treatment + month + (1|ID), 
                family = binomial, data = toenail)
summary(m.glmm)

# Figure for the subject-specific mean
time <- seq(0, 10, by = 0.02)
cm <- coef(m.glmm)$ID
dim(cm)
ID.A <- unique(toenail$ID[toenail$treatment==0])  # subject in A
ind.A <- which(rownames(cm) %in% ID.A)      # rows for Subjects in A
cm <- cm[ind.A,]
p.glmm <- ilogit(cm[, "(Intercept)"] + outer(cm[, "month"], time))

plot(x = time, y = p.glmm[1, ], ylim = c(0, 1), ylab = "Estimated p_i", 
     xlab = "Month", main = "Subject-specific Estimated Means p_i", type = "l")
for(i in 1:30) {lines(x = time, y = p.glmm[i, ])}

# Add population average
p.ave <- colMeans(p.glmm)
lines(x = time, y = p.ave, lwd = 4, col = "red", lty = 2)

# Add X*beta
p.pop1 <- ilogit(fixef(m.glmm)["(Intercept)"] + fixef(m.glmm)["month"] * time)
lines(x = time, y = p.pop1, lwd = 4, col = "green", lty =3)
p.pop2 <- ilogit(mean(cm[, "(Intercept)"] + fixef(m.glmm)["month"]) * time)
lines(x = time, y = p.pop2, lwd = 4, col = "brown", lty = 4)

# GEE model
cm.gee <- coefficients(m.gee)
p.gee <- ilogit(cm.gee["(Intercept)"] + cm.gee["month"] * time)
lines(x = time, y = p.gee, lty = 5, lwd = 4, col = "blue")
legend(x = 7, y = 1, lty = c(2, 5, 4, 3), lwd = 4, 
       col = c("red", "blue", "brown", "green"), 
       legend = c("Ave of GLMM est. p_i", "GEE estimate", 
       "GLMM estimate using mean of coefs", "GLMM est. X*beta"), bg = "white")



##### Lecture Code Chapter19 Response Surface #####
ds1 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/response-surface-1.txt",
                  header = TRUE)
ds1
m1 <- lm(y ~ x1 + x2, ds1)
summary(m1)
anova(m1)

library(rsm)
m1.rsm <- rsm(y ~ FO(x1, x2), ds1)
summary(m1.rsm)
contour(m1.rsm, ~ x1 + x2)

m0 <- lm(y ~ as.factor(x1) + as.factor(x2), ds1)
anova(m0) # Q: doesn't agree with what Prof. taught?


ds2 <- read.table("http://users.stat.umn.edu/~birgit/classes/8052/data/response-surface-2.txt",
                  header = TRUE)

ds2$block <- as.factor(ds2$block)
options(contrasts=c("contr.treatment", "contr.poly"))
m2.rsm <- rsm(y ~ block + SO(a, b, c), ds2)
summary(m2.rsm)

par(mfrow = c(1, 3))
contour(m2.rsm, ~ a + b + c, at = summary(m2.rsm)$canonical$xs)  ## Q?
par(mfrow = c(1, 1))

steepest(m2.rsm, dist = seq(0, 1, by = 0.2))
