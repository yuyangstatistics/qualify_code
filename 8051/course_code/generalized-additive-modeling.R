# two dimensional kernel estimation

library(faraway)
library(sm)

data(savings)
help(savings)

y<-savings$sr
xx<-cbind(savings$pop15,savings$ddpi)
sm.regression(xx,y,h=c(1,1),xlab="pop15",ylab="growth",zlab="saving rate")
sm.regression(xx,y,h=c(5,5),xlab="pop15",ylab="growth",zlab="saving rate")

# Additive modeling

data(ozone)
help(ozone)
pairs(ozone)

olm<-lm(O3~temp+ibh+ibt,ozone)
summary(olm)
plot.lm(olm)   # it shows the linear model does not fit well

library(gam)
help(gam)

am1<-gam(O3~lo(temp)+lo(ibh)+lo(ibt),data=ozone) # lo: loess smoother
summary(am1)
par(mfrow=c(1,3))
plot(am1,residuals=TRUE,se=TRUE,pch=".")

am2<-gam(O3~temp+lo(ibh)+lo(ibt),data=ozone) # linear in temp 
summary(am2)
par(mfrow=c(1,3))
plot(am2,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am2))
plot(ozone$ibt,residuals(am2))
plot(ozone$ibh,residuals(am2))

anova(am2,am1,test="F")  # the simpler model is rejected

am3<-gam(O3~lo(vh)+lo(wind)+lo(humidity)+lo(ibh)+lo(ibt)+lo(dpg)+lo(vis)+lo(doy),data=ozone)
summary(am3)

am4<-gam(O3~vh+lo(humidity)+lo(ibh)+lo(ibt)+lo(dpg)+lo(vis)+lo(doy),data=ozone) 
summary(am4)

anova(am1,am4,test="F")
par(mfrow=c(3,3))
plot(am4,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am4))
plot(ozone$ibt,residuals(am4))
plot(ozone$ibh,residuals(am4))

# we may also try natural splines to fit additive models
attach(ozone)
am5<-glm(O3~ns(temp,2)+ibt+ns(ibh))
summary(am5)
plot(am5,residuals=TRUE,se=TRUE,pch=".")
plot(ozone$temp,residuals(am4))
plot(ozone$ibt,residuals(am4))
plot(ozone$ibh,residuals(am4))








