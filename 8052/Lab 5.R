# Lab 5

library(MASS)
library(nlme)
library(lmerTest)

data <- read.table("/Users/yuyang/Documents/Courses/19Spring/STAT8052/data/cd4.txt", header = TRUE)
data$treat <- as.factor(data$treat)
data$gender <- as.factor(data$gender)
data$race <- as.factor(data$race)
data$aids <- as.factor(data$aids)

# a.
options(contrasts = c("contr.sum", "contr.poly"))
m <- lm(cd4.m0 ~ treat, data = data)
par(mfrow = c(2,2))
plot(m) # check assumptions
par(mfrow = c(1,1))
# right skew
anova(m) # p = 1.225e-11
# it appears that some treatments are more likely to be prescribed to patients with lower CD4 cell counts
dummy.coef(m) # coefficient for treatment 3 is much lower 
boxplot(cd4.m0 ~ treat, data = data)
# m2 <- lme(cd4 ~ cd4.m0*treat + vl, data = data, random = ~ 1|pid)

# b
m2 <- lm(cd4 ~ cd4.m0*treat + gender + race + aids + age, data = data)
par(mfrow = c(2,2))
plot(m2)
anova(m2)

# treatments are not randomized, and the treatment assignment may have been related to starting CD4 count
# so one way to deal with this is to center the response by starting CD4 count

data$diff <- data$cd4 - data$cd4.m0
m2 <- lme(diff ~ age + treat + aids + month + age:treat + treat:month, data = data,
          random = ~1|pid) # used model selection to get these specific interactions
summary(m2)

m3 <- lme(cd4 ~ cd4.m0 + age + treat + aids + month + age:treat + treat:month, data = data,
          random = ~1|pid) # used model selection to get these specific interactions
summary(m3)

m4 <- lme(cd4 ~ age + treat + aids + month + age:treat + treat:month, data = data,
          random = ~1|pid) # used model selection to get these specific interactions
summary(m4)

# because we set contr.sum, slopes of each treatment: 1 = 13.5, 2 = -23.8, 3 = 10.4
# treatment 1 is most effective

m4 <- update(m2, fixed = diff ~ .)
anova(m4)

# another method of capturing the info of the treatments from the baseline
data$fu <- 1
data$fu[data$month == 0] <- 0
data$fu <- as.factor(data$fu)

m3 <- lme(cd4 ~ 1 + treat + fu + treat:fu, data = data, random = ~1|pid)
summary(m3)

# c.
# Mixed model assumptions, random effects are independent (uncorrelated)

# d.
# note that viral load suppression is binary; we'll get to a mixed model analysis for such a response later



# gee
library(geepack)
modgeep <- geeglm(vl ~ age + treat + aids + month + age:treat + treat:month, 
                  id=pid, corstr="exchangeable", scale.fix=TRUE,  data=data,
                  family=binomial)
summary(modgeep)



modgeep <- geeglm(vl ~ age + treat + aids + month, 
                  id=pid, corstr="exchangeable", scale.fix=TRUE,  data=data,
                  family=binomial)
summary(modgeep)

modgeep2 <- geeglm(vl ~ age + aids + month, 
                   id=pid, corstr="exchangeable", scale.fix=TRUE,  data=data,
                   family=binomial)
summary(modgeep2)
anova(modgeep, modgeep2)

modgh <- glmer(vl ~ age + treat + aids + month + age:treat + treat:month + (1|pid), nAGQ=25, family=binomial, data=data)
summary(modgh)