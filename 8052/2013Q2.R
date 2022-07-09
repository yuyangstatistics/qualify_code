library(gmodels)
library(dplyr)
library(MASS)
library(car)

##### ============= part (a) ============= ##### 
da <- read.table("http://users.stat.umn.edu/~wangx346/artstudy-widedata.txt", 
                 header = T)
View(da)
da <- within(da, {
  ID <- as.factor(ID)
  gender <- as.factor(gender)
  country <- as.factor(country)
  treat <- as.factor(treat)
  np.chg <- np.8 - 0.5 * (np.0 + np.2)
})

m <- lm(np.chg ~ treat, data = da)
fit.contrast(m, varname="treat", coeff = c(1, -1), 
             conf.int = 0.95)
par(mfrow = c(2, 2))
plot(m)
par(mfrow = c(1, 1))
boxcox(lm(np.chg+300 ~ treat, data = da))

nrow(da %>% filter(treat == 1))

##### ============= part (b) ============= ##### 
with(da, scatter.smooth(x = age, y = np.chg))

m <- lm(np.chg ~ (treat + I(age - 40) + I((age - 40)^2) + gender + 
                    country + art.pre0 + artdur2)^2, data = da)
stepAIC(m, scope = list(lower = ~ treat), direction = "backward", trace = F)

stepAIC(m, scope = list(lower = ~ (treat + art.pre0) * I(age - 40)), direction = "backward", trace = F)

m2 <- lm(np.chg ~ treat * (I(age - 40) + I((age - 40)^2)) + I((age - 40)^3) + 
           gender * country + art.pre0 * (I(age - 40) + I((age - 40)^2)), 
         data = da)
fit.contrast(m2, varname="treat", coeff = c(1, -1), 
             conf.int = 0.95)
summary(m2)


##### ============= part (e) ============= ##### 
require(nlme)
require(lattice)
dc <- read.table("http://users.stat.umn.edu/~wangx346/artstudy.txt", 
                 header = T)
dc <- within(dc, {
  ID <- as.factor(ID)
  gender <- as.factor(gender)
  country <- as.factor(country)
  treat2 <- ifelse(treat ==2 & t== 8, 1, 0)
})
dc <- dc %>% filter(t > 0)

m.lme <- lme(np ~ (t + I(t^2)) * (country + gender + I(age - 40) + I((age - 40)^2)) + 
               country:gender + treat2 * (I(age - 40) + I((age - 40)^2)) + 
               I((age - 40)^3) + art, random = ~ 1|ID, data = dc)
