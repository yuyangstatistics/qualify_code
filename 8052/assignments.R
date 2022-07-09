#####  Assignment 1  #####

## Randomization Test
y = c(1, 6, 4, 6, 2, 3, 5)
library(AlgDesign)
T = as.matrix(gen.factorial(2, 7)) %*% y
hist(T, main = "Randomization Distribution and Test Statistic")
abline(v = 19, col="red", lty = 2)
(p <- length(T[T >= 19]) / length(T))

## P3.2
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr3.2", header = T)
df$trt = factor(df$trt)
options(contrasts = c("contr.sum", "contr.poly"))
fit = lm(days ~ trt, df)

par(mfrow = c(1, 2))
plot(fit, which = c(1, 2))
par(mfrow = c(1, 1))
anova(fit)
boxplot(days ~ trt, df, main = "Fruit Flies Data", 
        xlab = "Treatment", ylab = "Days")

model.tables(aov(fit), type = "effects")
library(gmodels)
estimable(fit, cm = c(0, 1, 0, 0, 0), conf.int = 0.95)
fit.contrast(aov(fit), "trt", 
             coeff = c(0.8, -0.2, -0.2, -0.2, -0.2), conf.int = 0.95)

estimable(fit, cm = c(0, 1, 0, 0, -1), conf.int = 0.95)
fit.contrast(aov(fit), "trt", 
             coeff = c(1, 0, 0, -1, 0), conf.int = 0.95)
library(cfcdae)
linear.contrast(fit, term = trt, contr.coefs = c(1, 0, 0, -1, 0))

TukeyHSD(aov(fit), "trt")
pairwise(fit, term = trt, confidence = 0.95, type = "hsd")

pairwise(fit, term = trt, confidence = 0.95, type = "snk")

sort(model.tables(aov(fit), type = "means")$tables$trt)


#####  Assignment 2  #####
n <- 1:25
n[qt(1 - 0.05/(2*n), 25) <= sqrt(4 * qf(0.95, 4, 25))]


## E6.3
df <- read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/ex6.3", header = T)
df$drug <- as.factor(df$drug)
attach(df)

m0 <- lm(dose ~ drug)
par(mfrow = c(2, 2))
plot(m0)
par(mfrow = c(1, 1))

library(MASS)
bc <- boxcox(aov(m0))
i <- which.max(bc$y)
bc$x[i]

# 95% CI for lambda
(max.l <- bc$y[i])
n <- nrow(df)
(const <- n / 2 * log( 1 + qt(0.025, m0$df.residual)^2 / m0$df.residual))
lower.l <- max.l - const
i.lower <- min(which(bc$y >= lower.l))
i.upper <- max(which(bc$y >= lower.l))
(ci <- bc$x[c(i.lower, i.upper)])

logdose <- log(df$dose)
mt <- lm(logdose ~ drug)
par(mfrow = c(2, 2))
plot(mt)
par(mfrow = c(1, 1))

anova(mt)

sort(model.tables(aov(mt), type = "effects")$tables$drug)
TukeyHSD(aov(mt), "drug")

compare.to.best(mt, term = drug, confidence = 0.95, lowisbest = T)

detach(df)

## P6.1
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr6.1", header = T)
df$trt = factor(df$trt)

m0 <- lm(y ~ trt, df)
par(mfrow = c(2, 2))
plot(m0)
par(mfrow = c(1, 1))

library(MASS)
boxcox(aov(m0))

anova(m0)

sort(model.tables(aov(m0), type = "means")$tables$trt)
library(gmodels)
compare.to.best(m0, term = trt, confidence = 0.95, lowisbest = F)


#####  Assignment 3  #####

## E7.2
stats::power.anova.test(groups = 3, n = NULL, between.var = var(c(10, 11, 11)), 
                 within.var = 4, sig.level = 0.05, power = 0.9)


## E7.5
w <- c(-1/4, 1/2, -1/4, 1/2, -1/4, -1/4)
MSE <- 0.25
n <- 4
(gamma <- 1 / (MSE * t(w) %*% w / n))
dfE <- 24 - 6
q_null <- qf(1 - 0.01, 1, dfE, ncp = 0)
1 - pf(q_null, 1, dfE, ncp = gamma)

w <- c(1/2, 1/2, -1/2, -1/2, 0, 0)
n <- 2:10
se <- sqrt(MSE /n * c(t(w) %*% w) )
t_critical <- qt(1 - 0.05 / 2, 6*(n-1))
width <- 2 * t_critical * se
n[min(which(width <=2))]


## P7.2
stats::power.anova.test(groups = 3, n = NULL, 
                        between.var = var(c(45, 32, 60)), 
                        within.var = 35, sig.level = 0.01, 
                        power = 0.95)

n <- 10:80
gamma <- 5^2 / (35 * 2 / n)
dfE <- 3 * (n - 1)
q_null <- qf(1 - 0.01 / 2, 1, dfE, ncp = 0)
(power <- 1 - pf(q_null, 1, dfE, ncp = gamma))
n[min(which(power >= 0.95))]

# greedy search
n = 2
while (T) {
  a = c(1, -1, 0)
  delta = 5^2 / (35 / n * t(a) %*% a)
  q_null = qf(1 - .01/2, 1, 3*n - 3, ncp = 0)
  p = 1 - pf(q_null, 1, 3 * n - 3, ncp = delta)
  if (p <= .95) {
    n = n + 1
  } else {
    break
  }
}
print(n)


## P8.5
rm(list = ls())
df <- read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr8.5", header = T)
df <- within(df, {
  shape <- as.factor(shape)
  trt <- as.factor(trt)
})
attach(df)

m <- lm(y ~ shape * trt)
anova(m)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
library(MASS)
boxcox(m)

# Tukey's 1df test
m <- lm(y ~ shape + trt)
m.aov <- aov(m)
ypred <- m.aov$fitted.values
(mu <- model.tables(m.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m.onedf <- update(m, . ~ . + onedf)
anova(m.onedf)
eta <- m.onedf$coefficients["onedf"]
1 - eta

y2 <- sqrt(y)
m <- lm(y2 ~ shape + trt)
m.aov <- aov(m)
ypred <- m.aov$fitted.values
(mu <- model.tables(m.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m.onedf <- update(m, . ~ . + onedf)
anova(m.onedf)
eta <- m.onedf$coefficients["onedf"]
1 - eta
boxcox(m)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))

summary(m.aov)

# interaction plot
interaction.plot(x.factor = shape, trace.factor = trt, response = y2, 
                 fun = mean, legend = T, xpd = NA, type = 'b', 
                 trace.label = 'Trt', xlab = "Shape", ylab = "Mean Response", 
                 lty = c(1, 2, 4), col = c("black", "red", "blue"), 
                 main = "Interaction Plot for y2 ~ shape + trt")

# profile plot
interaction.plot(x.factor = shape, trace.factor = rep(1, length(y2)), 
                 response = y2, fun = mean, legend = F, xpd = NA, 
                 type = "l", xlab = 'Shape', ylab = 'Mean Response', 
                 main = 'Profile Plot for Shape')

interaction.plot(x.factor = trt, trace.factor = rep(1, length(y2)), 
                 response = y2, fun = mean, legend = F, xpd = NA, 
                 type = "l", xlab = 'Trt', ylab = 'Mean Response', 
                 main = 'Profile Plot for Trt')

# pairwise comparison
TukeyHSD(m.aov, "shape")
sort(model.tables(m.aov, type = "means")$tables$shape)

TukeyHSD(m.aov, "trt")
sort(model.tables(m.aov, type = "means")$tables$trt)


## (b)
options(contrasts = c("contr.sum", "contr.poly"))
m2 <- lm(y2 ~trt * shape)
m2.aov <- aov(m2)
summary(m2.aov)
model.matrix(m2)
TukeyHSD(m2.aov, "trt:shape")

detach(df)


## P8.6
rm(list = ls())
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr8.6", header = T)
df <- within(df, {
  ht <- as.factor(ht)
  fert <- as.factor(fert)
  int <- as.factor(int)
})
attach(df)

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ ht * fert * int)
anova(m)
library(faraway)
halfnorm(coef(m)[-1], nlab = 3, labs = names(coef(m))[-1])
title(main = "Half-normal plot, effects in y ~ ht * fert * int")

m <- lm(y ~ (ht + fert + int)^2)
anova(m)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

m.aov <- aov(m)

interaction.plot(x.factor = ht, trace.factor = int, response = y, 
                 fun = mean, legend = T, xpd = NA, type = 'b', 
                 trace.label = 'int', xlab = "ht", ylab = "Mean Response", 
                 lty = c(1, 2, 4), col = c("black", "red", "blue"), 
                 main = "Interaction Plot for y ~ (ht + fert + int)^2")

interaction.plot(x.factor = fert, trace.factor = int, response = y, 
                 fun = mean, legend = T, xpd = NA, type = 'b', 
                 trace.label = 'int', xlab = "fert", ylab = "Mean Response", 
                 lty = c(1, 2, 4), col = c("black", "red", "blue"), 
                 main = "Interaction Plot for y ~ (ht + fert + int)^2")

# (b)
ht1 = as.numeric(as.character(df$ht))
fit2 = lm(y ~ I(ht1)*(fert+int) + I(ht1^2)*(fert+int) + fert*int, data = df)
anova(fit2)

detach(df)

#####  Assignment 4  #####

## P10.3
rm(list = ls())
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr10.3", header = T)
df <- within(df, {
  gum <- as.factor(gum)
  protein <- as.factor(protein)
})
attach(df)
options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(response ~ gum * protein)
library(car)
Anova(m, type = 'II')
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

emmip(m, protein ~ gum)

m2 <- lm(response ~ gum + protein)
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2)
Anova(m2, type = "II")
library(emmeans)
emmeans(m2, specs = "gum")

library(multcomp)
g <- glht(aov(m2), linfct = mcp(gum = "Tukey"), 
          alternative = "two.sided")
summary(g)
plot(g)

TukeyHSD(aov(m2), "gum")
linear.contrast(m2, term = gum, allpairs = T) # the same as glht

detach(df)

## P10.5
rm(list = ls())
library(car)
library(MASS)
library(emmeans)
library(multcomp)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr10.5", header = T)
data <- within(data, {
  variety <- as.factor(variety)
  location <- as.factor(location)
  temperature <- as.factor(temperature)
})
attach(data)

options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(size ~ variety * location * temperature)
Anova(m, type = "II")
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

rate <- 1 / size
m2 <- lm(rate ~ variety * location * temperature)
Anova(m2, type = "II")
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2)

m3 <- lm(rate ~ variety + location + temperature)
Anova(m3, type = "II")
par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))
boxcox(m3)

emmeans(aov(m3), specs = "variety")
summary(glht(aov(m3), linfct = mcp(variety = "Tukey"), 
             alternative = "two.sided"))

emmeans(aov(m3), specs = "location")
summary(glht(aov(m3), linfct = mcp(location = "Tukey"), 
     alternative = "two.sided"))

emmeans(aov(m3), specs = "temperature")
summary(glht(aov(m3), linfct = mcp(temperature = "Tukey"), 
     alternative = "two.sided"))

detach(data)

## E13.5
rm(list = ls())
library(MASS)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/ex13.5", header = T)
names(data)
data <- within(data, {
  operator <- as.factor(operator)
  machine <- as.factor(machine)
  day <- as.factor(day)
  trt <- as.factor(trt)
})
attach(data)

m <- lm(y ~ operator + machine + day + trt)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

# Tukey's 1df test
m.aov <- aov(m)
ypred <- m.aov$fitted.values
(mu <- model.tables(m.aov, type = "means")$table$`Grand mean`)
onedf <- ypred * ypred / (2 * mu)
m.onedf <- update(m, . ~ . + onedf)
anova(m.onedf)
eta <- m.onedf$coefficients["onedf"]
1 - eta

anova(m)

detach(data)

## P13.4
rm(list = ls())
library(gmodels)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr13.4", header = T)
names(data)
data <- within(data, {
  student <- as.factor(student)
  grader <- as.factor(grader)
  exam <- as.factor(exam)
})
attach(data)
options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(score ~ student + grader + exam)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

y <- score^2
m <- lm(y ~ student + grader + exam)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

anova(m)
TukeyHSD(aov(m), "exam", ordered = T)
sort(model.tables(aov(m), type = "means")$tables$exam)

w <- rbind("new - old" = c(1/2, 1/2, -1/3, -1/3, -1/3))
fit.contrast(aov(m), "exam", w)

detach(data)


## P13.5
rm(list = ls())
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr13.5", header = T)
names(data)
data <- within(data, {
  year.loc <- as.factor(year.loc)
  rotation <- as.factor(rotation)
  variety <- as.factor(variety)
})
attach(data)
options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(y ~ year.loc + rotation * variety)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)

anova(m)

m2 <- lm(y ~ year.loc + rotation)
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2)

anova(m2)

TukeyHSD(aov(m2), "rotation", ordered = T)
sort(model.tables(aov(m2), type = "means")$tables$rotation)



#####  Assignment 5  #####

## P14.3
rm(list = ls())
library(MASS)
library(car)
library(multcomp)
library(emmeans)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr14.3", header = T)
names(data)
data <- within(data, {
  exam <- as.factor(exam)
  grader <- as.factor(grader)
})
attach(data)

options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(score ~ exam + grader)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)  # no need of transformation

Anova(m, type = 2)

g <- glht(aov(m), linfct = mcp(grader = "Tukey"), alternative = "two.sided")
# summary(g)
# plot(g)
confint(g)

emmeans(m, spec = "grader")

# or use manual calculation
B <- sapply(31:54, function(i) {coef(m)[1] + coef(m)[i]})
names(B) <- names(coef(m))[31:54]
B["grader25"] <- coef(m)[1] - sum(coef(m)[31:54])
sort(B)

(mse <- anova(m)$`Mean Sq`[3])
g <- 25
k <- 5
b <- 30
(r <- k * b / g)
(n_bibd <- r * g * (k - 1) / (g - 1) / k)
qtukey(0.95, 25, 96) / sqrt(2) * sqrt(mse * 2 / n_bibd)  
# the results conincide with the ones given by glht

detach(data)


## P10.7
rm(list = ls())
library(faraway)
library(FrF2)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr10.7", header = T)
names(data)
(cols <- names(data)[1:4])
data[cols] <- lapply(data[cols], factor)
attach(data)

options(contrasts = c("contr.sum", "contr.poly"))

# half-normal plot
m <- lm(y ~ complaint * age * timeknown * maturity)
anova(m)
halfnorm(-2 * coef(m)[-1], labs = names(coef(m))[-1], nlab = 3)
title(main = "Effects of y ~ complaint * age * timeknown * maturity")
DanielPlot(m, half = T)

# pool high interaction terms
m2 <- lm(y ~ (complaint + age + timeknown + maturity)^2)
anova(m2)
Anova(m2, type = 2)
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2) 

m3 <- lm(y ~ complaint + age + timeknown + maturity)
anova(m3)
Anova(m3, type = 2)
par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))
boxcox(m3) 


# Lenth's pseudo p-values
z <- 2 * coef(m)[-1]
s0 <- median(abs(z)) * 1.5
z2 <- z[abs(z) < s0 *2.5]
pse <- median(abs(z2)) * 1.5
t0 <- z / pse
k <- 4
p <- 2 * pt(-abs(t0), (2^k - 1)/3)
p[which(p < 0.05)]


## E15.2
library(conf.design)
G <- rbind(c(1, 0, 1, 1), c(1, 1, 0, 1))
colnames(G) <- LETTERS[1:4]
conf.set(G, p = 2)
conf.design(G, p = 2, block.name = "blocks")


## 15.3
library(conf.design)
G <- rbind(c(1, 1, 1, 0, 0, 1, 0, 0), 
           c(1, 1, 0, 1, 1, 0, 0, 0), 
           c(1, 0, 1, 1, 1, 0, 0, 0), 
           c(0, 1, 1, 1, 0, 0, 0, 1))
colnames(G) <- LETTERS[1:8]
conf.set(G, p=2)


## P15.6
rm(list = ls())
library(faraway)
library(FrF2)
library(MASS)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr15.6", header = T)
names(data)
cols <- names(data)[1:5]
data[cols] <- lapply(data[cols], factor)

attach(data)

options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(days ~ cooler + A * B * C * D)
anova(m)
alias(m)

halfnorm(-2 * coef(m)[-1], labs = names(coef(m))[-1], nlab = 3)
DanielPlot(m, half = T)

m1 <- lm(days ~ cooler + D)
par(mfrow = c(2, 2))
plot(m1)
par(mfrow = c(1, 1))
boxcox(m1)
anova(m1)
coef(m1)

detach(data)

#####  Assignment 6  #####

## E18.2
library(FrF2)
design.frf2 <- FrF2(nruns = 16, nfactors = 6, 
                    generators = c("-BCD", "ABD"))
design.frf2

# factor levels
(factor_comb <- apply(design.frf2, 1, function(row) paste(names(row[which(row == 1)]), collapse = "")))


## 18.5
library(FrF2)
FrF2(nruns = 64, nfactors = 8, blocks = 4)

FrF2(nfactors = 8, resolution = 5)


## 18.3
# create fake data to get aliases
rm(list = ls())
library(FrF2)
design.frf2 <- FrF2(nruns = 16, nfactors = 8, 
                    generators = c("BCD", "ACD", "ABC", "ABD"))
data <- as.data.frame(design.frf2)
data$y <- 1
data
sapply(data, class)

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(y ~ A*B*C*D*E*F*G*H, data = data)
anova(m)
aliases(m)

# get the defining relations
dr <- alias(m)$Complete[, "(Intercept)"]
dr[dr %in% c(1, -1)]


## 18.6
rm(list = ls())
library(FrF2)
library(MASS)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr18.6", header = T)
names(data)
cols <- names(data)[-length(names(data))]
data[cols] <- lapply(data[cols], factor)

options(contrasts = c("contr.sum", "contr.poly"))

m <- lm(thickness ~ A*B*C*D*E*F*G*H, data = data)
anova(m)
DanielPlot(m, half = T)

m2 <- lm(thickness ~ D, data = data)
anova(m2)

aliases(m)



## P17.1

# (a)
rm(list = ls())
library(car)
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr17.1", header = T)
names(data)
data$pesticide <- as.factor(data$pesticide)
attach(data)

options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(calcium.conc ~ diameter + pesticide)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m)
Anova(m, type = 2)
summary(m)
plot(TukeyHSD(aov(m), "pesticide", ordered = T, conf.level = 0.95))
library(emmeans)
emmeans(m, ~pesticide)
library(multcomp)
g <- glht(aov(m), linfct = mcp(pesticide="Tukey"), alternative = "two.sided")
summary(g)
plot(g)


# (b)
options(contrasts = c("contr.sum", "contr.poly"))
m2 <- lm(calcium.conc ~ diameter * pesticide)
par(mfrow = c(2, 2))
plot(m2)
par(mfrow = c(1, 1))
boxcox(m2)
Anova(m2, type = 2)

# (c)
# Influential Observations
# added variable plots
library(car)
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(data)-length(m2$coefficients)-2)) 
plot(m2, which=4, cook.levels=cutoff)
# Influence Plot 
influencePlot(m2, id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )

outlierTest(m2)
qqPlot(m2, main="QQ Plot") #qq plot for studentized resid 
leveragePlots(m2) # leverage plots

boxplot(data$diameter)

m3 <- lm(calcium.conc ~ diameter * pesticide, data = data[-2, ])
par(mfrow = c(2, 2))
plot(m3)
par(mfrow = c(1, 1))
boxcox(m3)
Anova(m3, type = 2)

# (d)
m4 <- lm(diameter ~ pesticide, data = data[-2, ])
par(mfrow = c(2, 2))
plot(m4)
par(mfrow = c(1, 1))
boxcox(m4)
anova(m4)

detach(data)


#####  Assignment 7  #####

## E11.3
rm(list = ls())
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/ex11.3", header = T)
names(data)
data$brand <- as.factor(data$brand)

m <- lm(y ~ brand, data = data)
anova(m)

MS.Trt <- anova(m)$`Mean Sq`[1]
MS.E <- anova(m)$`Mean Sq`[2]
df.Trt <- anova(m)$Df[1]
df.E <- anova(m)$Df[2]
n <- 6
L <- (MS.Trt / MS.E * qf(0.005, df.E, df.Trt) - 1) / n
U <- (MS.Trt / MS.E * qf(0.995, df.E, df.Trt) - 1) / n
c(L, U)

par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
library(MASS)
boxcox(m)

# (b)
sigma_a <- (MS.Trt - MS.E) / n
sigma_a / (sigma_a + MS.E)

# (c)
L <- (MS.Trt / MS.E * qf(0.025, df.E, df.Trt) - 1) / n
U <- (MS.Trt / MS.E * qf(0.975, df.E, df.Trt) - 1) / n
c(L / (1 + L), U / (1 + U))


## E11.5
(1.288 - 0.0921) / 3
c(1.842 / qchisq(0.975, 20), 1.842 / qchisq(0.025, 20))

## E11.6
a <- 2:50
(power <- 1 - pf(qf(0.99, a - 1, a) / 3, a-1, a))
a[min(which(power > 0.9))]


## P11.1
a <- c(10, 50)
n <- 100 / a
(qf(0.975, a*n-a, a-1) - qf(0.025, a*n-a, a-1))/n


## P11.2
rm(list = ls())
data = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr11.2", header = T)
names(data)
cols <- names(data)[-length(names(data))]
data[cols] <- lapply(data[cols], factor)
sapply(data, class)

options(contrasts = c("contr.sum", "contr.poly"))

m <- aov(bacteria ~ lab * sample, data)
summary(m)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
library(MASS)
boxcox(m)

# 95% CI for lambda
bc <- boxcox(m)
i <- which.max(bc$y)
bc$x[i]
(max.l <- bc$y[i])
n <- nrow(data)
(const <- n / 2 * log( 1 + qt(0.025, m$df.residual)^2 / m$df.residual))
lower.l <- max.l - const
i.lower <- min(which(bc$y >= lower.l))
i.upper <- max(which(bc$y >= lower.l))
(ci <- bc$x[c(i.lower, i.upper)])

data$newy <- (data$bacteria)^(bc$x[i])
m1 <- aov(newy ~ lab * sample, data)
summary(m1)
par(mfrow = c(2, 2))
plot(m1)
par(mfrow = c(1, 1))
boxcox(m1)

m2 <- aov(newy ~ lab + sample + Error(lab:sample), data)
summary(m2)


#####  Assignment 9  #####

### Oxboys problem

## (a)
library(lattice)
library(nlme)
data(Oxboys, package = "nlme")
xyplot(height ~ age | Subject, data = Oxboys, strip = T, grid = T)

m1 <- lme(height ~ age, data = Oxboys, random = ~ age | Subject)
anova(m1)

m2 <- update(m1, random = ~ 1 | Subject)
anova(m1, m2)

m3 <- update(m2, method = "ML")
m4 <- update(m3, fixed = height ~ age + I(age^2))
anova(m3, m4)

m <- lme(height ~ age + I(age^2), data = Oxboys, 
         random = ~ age + I(age^2) | Subject, method = "REML")
summary(m)

## (b)
plot(augPred(m), layout = c(5, 6), main = "Predicted, per Subject", 
     grid = T)
xyplot(resid(m) ~ age | Subject, data = Oxboys, layout = c(5, 6),
       panel = function(x, y, ...){
         panel.xyplot(x, y, ...)
         panel.lines(range(x), c(0, 0))
       })

### nepali problem
## (a)
library(dplyr)
library(lme4)
library(ggplot2)
data(nepali, package = "faraway")
nepali <- nepali %>% select(-c(ht)) %>% filter(!is.na(wt)) %>% 
  mutate(sex = as.factor(ifelse(sex == 1, "male", "female"))) %>%
  mutate(id = as.factor(id), lit = as.factor(lit))
ggplot(nepali, aes(x=age, y=wt, group=id)) + geom_point() +
  facet_wrap(~ sex)

## (b)
library(MASS)
library(car)
m.lm <- lm(wt ~ age + sex + mage + lit + died, data = nepali)
par(mfrow = c(2, 2))
plot(m.lm)
par(mfrow = c(1, 1))
boxcox(m.lm)
Anova(m.lm, type = 2)

## (c) & (d)
library(nlme)
m.lme <- lme(wt ~ age * sex + mage + lit, data = nepali, 
             random = ~ 1 | id)
Anova(m.lme, type = 2)
plot(m.lme)
qqnorm(residuals(m.lme), main = "QQ plot for residuals")
qqline(residuals(m.lme))
qqnorm(ranef(m.lme)$"(Intercept)", main = "QQ plot for the random intercept")
qqline(ranef(m.lme)$"(Intercept)")
fixef(m.lme)

## (e)
m.lme2 <- update(m.lme, method = "ML")
m.lme3 <- update(m.lme2, fixed = ~ age + mage)
anova(m.lme2, m.lme3)


### Oxboys problem 2
## (a)
library(nlme)
library(mgcv)
library(car)
data(Oxboys, package = "nlme")
head(Oxboys)
m <- lme(height ~ age, data = Oxboys, random = ~ 1 | Subject)
summary(m)
intervals(m)


#####  Assignment 10  #####

### P13.5
## (a)
library(nlme)
library(lme4)
library(lmerTest)
library(MASS)
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr13.5", header = T)
df <- within(df, {
  year.loc <- as.factor(year.loc)
  rotation <- as.factor(rotation)
  variety <- as.factor(variety)
})

m <- lm(y ~ year.loc + rotation * variety, data = df)
anova(m)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(m, data = df)

## (b)
m2 <- lme(y ~ rotation * variety, random = ~ 1 | year.loc, data = df)
anova(m2)

m3 <- lmerTest::lmer(y ~ rotation * variety + (1|year.loc), data = df, 
                     REML = T)
anova(m3, type = 2)

## (c)
?corClasses
m4 <- update(m2, correlation = corCompSymm())
anova(m2, m4)


### Faraway 10.6
## (a)
library(ggplot2)
library(lme4)
data(breaking, package = "faraway")
ggplot(breaking, aes(x = supplier, y = y, shape = day, color = operator)) + 
  geom_point()

## (b)
m <- lm(y ~ day + operator + supplier, data = breaking)
anova(m)

## (c)
options(contrasts = c("contr.treatment", "contr.poly"))
m2 <- lmer(y ~ supplier + (1|day) + (1|operator), data = breaking)
anova(m2, type = 2)
fixef(m2)

## (d)
library(RLRsim)
m.day <- lmer(y ~ supplier + (1|day), data = breaking)
m.op <- lmer(y ~ supplier + (1|operator), data = breaking)
# test operator
exactRLRT(lmer(y ~ (1|operator), data = breaking), m2, m.day)
# test days
exactRLRT(lmer(y ~ (1|day), data = breaking), m2, m.op)

## (e)
library(pbkrtest)
m3 <- lmer(y ~ (1|day) + (1|operator), data = breaking)
KRmodcomp(m2, m3)

## (f)
fixef(m2)
summary(m2)
pnorm(1000, mean = 622.5 + 411.25, sd = sqrt(219.1 + 4990.3))


### P16.8
## (a)
rm(list = ls())
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr16.8", header = T)
df <- within(df, {
  sign <- as.factor(sign)
  timing <- as.factor(timing)
  interchange <- as.factor(interchange)
})

m <- aov(speed ~ sign + Error(interchange:sign) + timing + sign:timing, 
         data = df)
m <- aov(speed ~ sign * timing + Error(interchange:sign), data = df)
summary(m)

## (b)
library(lme4)
library(emmeans)
library(multcomp)
m <- lmer(speed ~ sign * timing + (1|interchange:sign), data = df)
anova(m, type = 2)
pairs(emmeans(m, "timing"))
confint(pairs(emmeans(m, "timing")))
g <- glht(m, linfct=mcp(timing="Tukey"), alternative = "two.sided")
summary(g)
plot(g)

## (d)
head(getME(m, "X"))
head(getME(m, "Z"))
getME(m, "sigma")^2 * diag(nrow = 3)
VarCorr(m)


### P16.3
## (a)
rm(list = ls())
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr16.3", header = T)
df <- within(df, {
  replicate <- as.factor(replicate)
  pesticide <- as.factor(pesticide)
  irigation <- as.factor(irigation)
  variety <- as.factor(variety)
})

library(nlme)
m <- lme(yield ~ pesticide * irigation * variety, data = df, 
         random = ~1|replicate/pesticide/irigation)
anova(m)
summary(contrast(emmeans(m, "pesticide"), method="pairwise"),
        infer=c(T,T))
summary(contrast(emmeans(m, "irigation"), method="pairwise"),
        infer=c(T,T))
summary(contrast(emmeans(m, "variety"), method="pairwise"),
        infer=c(T,T))


#####  Assignment 11  #####

### Faraway Problem 13.4
## (a)
library(ggplot2)
data(nitrofen, package = "boot")
nitrofen$conc <- as.factor(nitrofen$conc)
ggplot(nitrofen, aes(x=conc, y=total)) + geom_point()

library(reshape2)
dat <- melt(nitrofen, measure.vars = c("brood1","brood2","brood3"))
names(dat)[3] <- "brood"
dat$brood <- as.factor(rep(1:3, each = 50))
dat$id <- rep(1:50, 3)
head(dat)
ggplot(dat, aes(x=conc, y=value)) + geom_point() + facet_wrap(~brood)
ggplot(dat, aes(x=conc, y=value)) + geom_boxplot() + facet_grid(~brood)

## (b)
library(lme4)
dat$conc <- scale(as.numeric(dat$conc))
options(contrasts = c("contr.treatment", "contr.poly"))
g0 <- glmer(value ~ brood * conc + (1|id), family = poisson, data = dat, 
            control = glmerControl(check.conv.singular = "warning", 
                                   optCtrl = list(maxfun = 10000)))
summary(g0)

## (c)
library(geepack)
g1 <- geeglm(value ~ brood * conc, id = id, family = poisson, data = dat, 
             corstr = "exchangeable", scale.fix = TRUE)
summary(g1)


### Faraway Problem 13.2
rm(list = ls())
library(reshape2)
library(dplyr)
library(ggplot2)
data(potuse, package = "faraway")
head(potuse)
potuseMelt <- melt(potuse, measure.vars = c(colnames(potuse[2:6])), 
                   variable.name = "year", value.name = "category")
levels(potuseMelt$year) <- c(76:80)
potuseMelt <- within(potuseMelt, {
  year <- as.numeric(year) + 75
  category <- as.factor(category)
  sex <- ifelse(sex==1, 'Male', 'Female')
})
CumSum <- potuseMelt %>% group_by(sex, year, category) %>% 
  summarise(total = sum(count))
ggplot(CumSum, aes(x=year, y=total, linetype=category)) + geom_line() + 
  facet_grid(~sex)


## (b)
potuse[, 2:6] <- ifelse(potuse[,2:6]==1, 0, 1)
potuseExpanded <- potuse[rep(row.names(potuse), potuse$count),1:6]
potuseExpanded$id <- c(1:nrow(potuseExpanded))
dat <- melt(potuseExpanded, measure.vars = c(colnames(potuseExpanded)[2:6]), 
            variable.name = "year", value.name = "y")
levels(dat$year) <- c(1:5)
dat <- within(dat, {
  year <- as.numeric(year)
  sex <- as.factor(sex)
  sex <- ifelse(sex==1, 'Male', 'Female')
})

options(contrasts = c("contr.treatment", "contr.poly"))
mod1 <- glmer(y ~ year*sex + (1|id), family = binomial, 
              control = glmerControl(check.conv.singular="warning", 
                                     optCtrl = list(maxfun = 10000)), 
              data = dat)
summary(mod1)
mod2 <- glmer(y ~ sex + year + (1|id), family = binomial, data = dat)
anova(mod1, mod2)
summary(mod2)

## (c)
mod3 <- glmer(y ~ year + (1|id), family = binomial, data = dat)
anova(mod2, mod3)

## (d)
mod4 <- glmer(y ~ sex + as.factor(year) + (1|id), family = binomial, data = dat)
anova(mod2, mod4)
summary(mod4)

## (e)
library(geepack)
mod.gee <- geeglm(y ~ sex + year, id = id, family = binomial, data = dat, 
                  corstr = "exchangeable", scale.fix = TRUE)
summary(mod.gee)


### Faraway Problem 13.5
## (a)
library(dplyr)
library(ggplot2)
data(toenail, package = "faraway")
head(toenail)
toenailProp <- toenail %>% group_by(visit, treatment) %>% 
  summarise(prop = mean(outcome))
toenailProp$treatment <- as.factor(toenailProp$treatment)
ggplot(toenailProp, aes(x=visit, y=prop, linetype = treatment)) + 
  geom_line() + labs(x="Visit", y="Proportion of patients with severe condition")

## (b)
library(lme4)
options(contrasts = c("contr.treatment", "contr.poly"))
toenail <- within(toenail, {
  ID <- as.factor(ID)
  visit <- as.factor(visit)
  fu <- ifelse(visit==1, 0, 1)
})
m1 <- glmer(outcome ~ treatment * fu + month + (1|ID), family = binomial, 
            data = toenail)
summary(m1)

## (c)
m2 <- glmer(outcome ~ fu + month + (1|ID), family = binomial, 
            data = toenail)
anova(m1, m2)




#####  Assignment 12  #####
rm(list = ls())
library(rsm)
library(MASS)
df = read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr19.5", header = T)
names(df)[c(2, 3)] <- c(names(df)[c(3, 2)])
df$day <- as.factor(df$day)
fit1 <- rsm(strength ~ day + FO(time, temperature), data = df)
summary(fit1)
fit2 <- rsm(strength ~ day + SO(time, temperature), data = df)
summary(fit2)
fit <- rsm(strength ~ day + FO(time, temperature) + PQ(temperature), 
           data = df)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
par(mfrow = c(1, 1))
boxcox(aov(fit))

contour(fit, ~ time + temperature)
canonical.path(fit)
steepest(fit)


### P19.3
rm(list = ls())
df <- read.table("http://www.stat.umn.edu/~gary/book/fcdae.data/pr19.3", header=T)
df$block <- as.factor(df$block)
library(rsm)
library(MASS)
fit1 <- rsm(y ~ block + FO(load, speed, advance, ratio, egr), data = df)
summary(fit1)
fit2 <- rsm(y ~ block + SO(load, speed, advance, ratio, egr), data = df)
summary(fit2)
par(mfrow = c(2, 2))
plot(fit2)
par(mfrow = c(1, 1))
boxcox(aov(fit2))

fit3 <- rsm(log(y) ~ block + SO(load, speed, advance, ratio, egr), data = df)
summary(fit3)

fit4 <- rsm(log(y) ~ block + FO(load, speed, advance, ratio) + 
              TWI(load, ratio) + TWI(advance, ratio) + 
              PQ(speed, advance, ratio), data = df)
summary(fit4)

canonical.path(fit4)
