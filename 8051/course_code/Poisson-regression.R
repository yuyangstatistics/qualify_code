
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
