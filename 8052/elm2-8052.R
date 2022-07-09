#####============== Chapter 10: Random Effects ============#####
library(faraway)
data(pulp, package="faraway")
op <- options(contrasts=c("contr.sum", "contr.poly"))
lmod <- aov(bright ~ operator, pulp)
summary(lmod)
library(ggplot2)
ggplot(pulp, aes(x=operator, y=bright))+geom_point(position = position_jitter(width=0.1, height=0.0))
coef(lmod)
options(op)
(0.447-0.106)/5
library(lme4)
mmod <- lmer(bright ~ 1+(1|operator), pulp)
summary(mmod)
sumary(mmod)
smod <- lmer(bright ~ 1+(1|operator), pulp, REML=FALSE)
sumary(smod)

nullmod <- lm(bright ~ 1, pulp)
lrtstat <- as.numeric(2*(logLik(smod)-logLik(nullmod)))
pvalue <- pchisq(lrtstat,1,lower=FALSE)
data.frame(lrtstat, pvalue)
y <- simulate(nullmod)
lrstat <- numeric(1000)
set.seed(123)
for(i in 1:1000){
  y <- unlist(simulate(nullmod))
  bnull <- lm(y ~ 1)
  balt <- lmer(y ~ 1 + (1|operator), pulp, REML=FALSE)
  lrstat[i] <- as.numeric(2*(logLik(balt)-logLik(bnull)))
}
mean(lrstat < 0.00001)
mean(lrstat > 2.5684)
sqrt(0.019*0.981/1000)
library(RLRsim)
exactLRT(smod, nullmod)
exactRLRT(mmod)
VarCorr(mmod)
as.data.frame(VarCorr(mmod))
bsd <- numeric(1000)
for(i in 1:1000){
  y <- unlist(simulate(mmod))
  bmod <- refit(mmod, y)
  bsd[i] <- as.data.frame(VarCorr(bmod))$sdcor[1]
}
quantile(bsd, c(0.025, 0.975))
confint(mmod, method="boot")

#####============== 10.3 Estimating Random Effects ============#####
ranef(mmod)$operator
(cc <- model.tables(lmod))
cc[[1]]$operator/ranef(mmod)$operator
library(lattice)
dotplot(ranef(mmod, condVar=TRUE))

#####============== 10.4 Prediction ============#####
fixef(mmod)+ranef(mmod)$operator
predict(mmod, re.form=~0)[1]
predict(mmod, newdata=data.frame(operator="a"))
group.sd <- as.data.frame(VarCorr(mmod))$sdcor[1]
resid.sd <- as.data.frame(VarCorr(mmod))$sdcor[2]
pv <- numeric(1000)
for(i in 1:1000){
  y <- unlist(simulate(mmod))
  bmod <- refit(mmod, y)
  pv[i] <- predict(bmod, re.form=~0)[1] + rnorm(n=1,sd=group.sd) + rnorm(n=1,sd=resid.sd)
}
quantile(pv, c(0.025, 0.975))
for(i in 1:1000){
  y <- unlist(simulate(mmod, use.u=TRUE))
  bmod <- refit(mmod, y)
  pv[i] <- predict(bmod, newdata=data.frame(operator="a")) + rnorm(n=1,sd=resid.sd)
}
quantile(pv, c(0.025, 0.975))

#####============== 10.5 Diagnostics ============#####
qqnorm(residuals(mmod),main="")
qqline(residuals(mmod))
plot(fitted(mmod),residuals(mmod),xlab="Fitted",ylab="Residuals")
abline(h=0)

#####============== 10.6 Blocks as Random Effects ============#####
library(faraway)
library(ggplot2)
data(penicillin, package="faraway")
summary(penicillin)
penicillin$Blend <- gl(5,4)
ggplot(penicillin, aes(y=yield, x=treat, shape=Blend))+geom_point()+xlab("Treatment")
ggplot(penicillin, aes(y=yield, x=Blend, shape=treat))+geom_point()
op <- options(contrasts=c("contr.sum", "contr.poly"))
lmod <- aov(yield ~ blend + treat, penicillin)
summary(lmod)
coef(lmod)
mmod <- lmer(yield ~ treat + (1|blend), penicillin)
sumary(mmod)
options(op)
ranef(mmod)$blend
anova(mmod)

amod <- aov(yield ~ treat + Error(blend), penicillin)
summary(amod)
# or use the following
summary(aov(yield ~ treat + blend + Error(blend:treat), penicillin))

lme.mod <- lme(yield ~ treat, data = penicillin, random = ~ 1|blend)
anova(lme.mod)


## test fixed effects
library(pbkrtest)
amod <- lmer(yield ~ treat + (1|blend), penicillin, REML=FALSE)
nmod <- lmer(yield ~ 1 + (1|blend), penicillin, REML=FALSE)
KRmodcomp(amod, nmod)
as.numeric(2*(logLik(amod)-logLik(nmod)))
1-pchisq(4.0474,3)
lrstat <- numeric(1000)
for(i in 1:1000){
  ryield <- unlist(simulate(nmod))
  nmodr <- refit(nmod, ryield)
  amodr <- refit(amod, ryield)
  lrstat[i] <- 2*(logLik(amodr)-logLik(nmodr))
}
mean(lrstat > 4.0474)
pmod <- PBmodcomp(amod, nmod)
summary(pmod)

## test random effects
rmod <- lmer(yield ~ treat + (1|blend), penicillin)
nlmod <- lm(yield ~ treat, penicillin)
as.numeric(2*(logLik(rmod)-logLik(nlmod,REML=TRUE)))
lrstatf <- numeric(1000)
for(i in 1:1000){
  ryield <-  unlist(simulate(nlmod))
  nlmodr <- lm(ryield ~ treat, penicillin)
  rmodr <- refit(rmod, ryield)
  lrstatf[i] <- 2*(logLik(rmodr)-logLik(nlmodr,REML=TRUE))
}
mean(lrstatf < 0.00001)
mean(lrstatf > 2.7629)
library(RLRsim)
exactRLRT(rmod)

#####============== 10.7 Split Plots ============#####
## Split-plots example
data(irrigation,package="faraway")
summary(irrigation)
ggplot(irrigation, aes(y=yield, x=field, shape=irrigation, color=variety)) + geom_point()
lmer(yield ~ irrigation * variety + (1|field) +(1|field:variety),  irrigation)
lmod <- lmer(yield ~ irrigation * variety + (1|field), irrigation)
sumary(lmod)
library(pbkrtest)
lmoda <- lmer(yield ~ irrigation + variety + (1|field),data=irrigation)
KRmodcomp(lmod, lmoda)
lmodi <- lmer(yield ~ irrigation + (1|field), irrigation)
KRmodcomp(lmoda, lmodi)
lmodv <- lmer(yield ~  variety + (1|field), irrigation)
KRmodcomp(lmoda, lmodv)
plot(fitted(lmod),residuals(lmod),xlab="Fitted",ylab="Residuals")
qqnorm(residuals(lmod),main="")
qqline(residuals(lmod))

library(RLRsim)
exactRLRT(lmod)
mod <- lm(yield ~ irrigation * variety, data=irrigation)
anova(mod)

#####============== 10.8 Nested Effects ============#####
data(eggs, package="faraway")
summary(eggs)
library(ggplot2)
ggplot(eggs, aes(y=Fat, x=Lab, color=Technician, shape=Sample)) + geom_point(position = position_jitter(width=0.1, height=0.0))+scale_color_grey()
cmod <- lmer(Fat ~ 1 + (1|Lab) + (1|Lab:Technician) +  (1|Lab:Technician:Sample), data=eggs)
sumary(cmod)
cmodr <- lmer(Fat ~ 1 + (1|Lab) + (1|Lab:Technician), data=eggs)
lrstat <- numeric(1000)
for(i in 1:1000){
  rFat <- unlist(simulate(cmodr))
  nmod <- lmer(rFat ~ 1 + (1|Lab) + (1|Lab:Technician), data=eggs)
  amod <- lmer(rFat ~ 1 + (1|Lab) + (1|Lab:Technician) +
                 (1|Lab:Technician:Sample), data=eggs)
  lrstat[i] <- 2*(logLik(amod)-logLik(nmod))
}
mean(lrstat > 2*(logLik(cmod)-logLik(cmodr)))
library(RLRsim)
cmods <- lmer(Fat ~ 1 + (1|Lab:Technician:Sample), data=eggs)
exactRLRT(cmods, cmod, cmodr)
VarCorr(cmodr)
confint(cmod, method="boot")

cmodt <- lmer(Fat ~ 1 + (1|Lab:Technician), data=eggs)
cmodl <- lmer(Fat ~ 1 + (1|Lab), data=eggs)
exactRLRT(cmodt, cmodr, cmodl)

#####============== 10.9 Crossed Effects ============#####
data(abrasion, package="faraway")
matrix(abrasion$material,4,4)
library(ggplot2)
ggplot(abrasion,aes(x=material, y=wear, shape=run, color=position))+geom_point(position = position_jitter(width=0.1, height=0.0))+scale_color_grey()
lmod <- aov(wear ~ material + run + position, abrasion)
summary(lmod)
mmod <- lmer(wear ~ material + (1|run) + (1|position), abrasion)
sumary(mmod)
library(RLRsim)
mmodp <- lmer(wear ~ material + (1|position), abrasion)
mmodr <- lmer(wear ~ material + (1|run), abrasion)
exactRLRT(mmodp, mmod, mmodr)
exactRLRT(mmodr, mmod, mmodp)
library(pbkrtest)
mmod <- lmer(wear ~ material + (1|run) + (1|position), abrasion,REML=FALSE)
nmod <- lmer(wear ~ 1+ (1|run) + (1|position), abrasion,REML=FALSE)
KRmodcomp(mmod, nmod)

#####============== 10.10 Multilevel Models example ============#####
data(jsp, package="faraway")
jspr <- jsp[jsp$year==2,]
ggplot(jspr, aes(x=raven, y=math))+xlab("Raven Score")+ylab("Math Score")+geom_point(position = position_jitter(),alpha=0.3)
ggplot(jspr, aes(x=social, y=math))+xlab("Social Class")+ylab("Math Score")+geom_boxplot()
glin <- lm(math ~ raven*gender*social,jspr)
anova(glin)
glin <- lm(math ~ raven*social,jspr)
anova(glin)
glin <- lm(math ~ raven+social,jspr)
summary(glin)
table(jspr$school)
mmod <- lmer(math ~ raven*social*gender+(1|school)+(1|school:class),  data=jspr, 
             REML = F)
mmodr <- lmer(math ~ raven*social+(1|school)+(1|school:class),  data=jspr, 
              REML = F)
KRmodcomp(mmod, mmodr)

all3 <- lmer(math ~ raven*social*gender+(1|school)+(1|school:class), data=jspr, 
             REML=FALSE)
all2 <- update(all3, . ~ . - raven:social:gender)
notrs <- update(all2, . ~ . -raven:social)
notrg <- update(all2, . ~ . -raven:gender)
notsg <- update(all2, . ~ . -social:gender)
onlyrs <- update(all2, . ~ . -social:gender - raven:gender)
all1 <-  update(all2, . ~ . -social:gender - raven:gender - social:raven)
nogen <- update(all1, . ~ . -gender)
anova(all3, all2, notrs, notrg, notsg, onlyrs, all1, nogen)[,1:4]

jspr$craven <- jspr$raven-mean(jspr$raven)
mmod <- lmer(math ~ craven*social+(1|school)+(1|school:class),jspr)
sumary(mmod)
diagd <- fortify(mmod)
ggplot(diagd,aes(sample=.resid))+stat_qq()
ggplot(diagd,aes(x=.fitted,y=.resid)) +geom_point(alpha=0.3) +geom_hline(yintercept=0) +xlab("Fitted") +ylab("Residuals")
qqnorm(ranef(mmod)$school[[1]],main="School effects")
qqline(ranef(mmod)$school[[1]])
qqnorm(ranef(mmod)$"school:class"[[1]],main="Class effects")
qqline(ranef(mmod)$"school:class"[[1]])
adjscores <- ranef(mmod)$school[[1]]
sort(adjscores)
rawscores <- coef(lm(math ~ school-1,jspr))
rawscores <- rawscores-mean(rawscores)
plot(rawscores,adjscores)
sint <- c(9,14,29)
text(rawscores[sint],adjscores[sint]+0.2,c("9","15","30"))
library(RLRsim)
mmodc <- lmer(math ~ craven*social+(1|school:class),jspr)
mmods <- lmer(math ~ craven*social+(1|school),jspr)
exactRLRT(mmodc, mmod, mmods)
exactRLRT(mmods, mmod, mmodc)
schraven <- lm(raven ~ school, jspr)$fit
mmodc <- lmer(math ~ craven*social+schraven*social+(1|school)+  (1|school:class),jspr)
KRmodcomp(mmod, mmodc)




#####======== Chapter 11 Repeated Measures and Longitudinal Data ========#####

library(faraway)
#####================ 11.1 Longitudinal Data ================#####
data(psid, package="faraway")
head(psid)
library(dplyr)
psid20 <- psid %>% filter(person <=20)
library(ggplot2)
ggplot(psid20, aes(x=year, y=income))+geom_line()+facet_wrap(~ person)
ggplot(psid20, aes(x=year, y=income+100, group=person)) +geom_line() + facet_wrap(~ sex) + scale_y_log10()
lmod <- lm(log(income) ~ I(year-78), subset=(person==1), psid)
coef(lmod)
library(lme4)
ml <- lmList(log(income) ~ I(year-78) | person, psid)
intercepts <- sapply(ml,coef)[1,]
slopes <- sapply(ml,coef)[2,]
plot(intercepts,slopes,xlab="Intercept",ylab="Slope")
psex <- psid$sex[match(1:85,psid$person)]
boxplot(split(slopes,psex))
t.test(slopes[psex=="M"],slopes[psex=="F"])
t.test(intercepts[psex=="M"],intercepts[psex=="F"])
library(lme4)
psid$cyear <- psid$year-78
options(contrasts = c("contr.treatment", "contr.poly"))
mmod <- lmer(log(income) ~ cyear*sex +age+educ+(cyear|person),psid)
sumary(mmod, digits=3)

library(pbkrtest)
mmod <- lmer(log(income) ~ cyear*sex +age+educ+(cyear|person),psid, REML=FALSE)
mmodr <- lmer(log(income) ~ cyear + sex +age+educ+(cyear|person),psid, REML=FALSE)
KRmodcomp(mmod,mmodr)
confint(mmod, method="boot")

diagd <- fortify(mmod)
ggplot(diagd,aes(sample=.resid))+stat_qq()+facet_grid(~sex)
diagd$edulevel <- cut(psid$educ,c(0,8.5,12.5,20), labels=c("lessHS","HS","moreHS"))
ggplot(diagd, aes(x=.fitted,y=.resid)) + geom_point(alpha=0.3) + geom_hline(yintercept=0) + facet_grid(~ edulevel) + xlab("Fitted") + ylab("Residuals")


#####================ 11.2 Repeated Measures ================#####
data(vision, package="faraway")
vision$npower <- rep(1:4,14)
ggplot(vision, aes(y=acuity, x=npower, linetype=eye)) + geom_line() + facet_wrap(~ subject, ncol=4) + scale_x_continuous("Power",breaks=1:4,labels=c("6/6","6/18","6/36","6/60"))
mmod <- lmer(acuity~power + (1|subject) + (1|subject:eye),vision)
sumary(mmod)
4.64^2/(4.64^2+3.21^2+4.07^2)
(4.64^2+3.21^2)/(4.64^2+3.21^2+4.07^2)
library(pbkrtest)
mmod <- lmer(acuity~power+(1|subject)+(1|subject:eye),vision,REML=FALSE)
nmod <- lmer(acuity~1+(1|subject)+(1|subject:eye),vision,REML=FALSE)
KRmodcomp(mmod, nmod)
mmodr <- lmer(acuity~power+(1|subject)+(1|subject:eye),vision,REML=FALSE, subset=-43)
nmodr <- lmer(acuity~1+(1|subject)+(1|subject:eye),vision,REML=FALSE, subset=-43)
KRmodcomp(mmodr, nmodr)
op <- options(contrasts=c("contr.helmert", "contr.poly"))
mmodr <- lmer(acuity~power+(1|subject)+(1|subject:eye),vision,subset=-43)
sumary(mmodr)
options(op)
contr.helmert(4)
plot(resid(mmodr) ~ fitted(mmodr),xlab="Fitted",ylab="Residuals")
abline(h=0)
qqnorm(ranef(mmodr)$"subject:eye"[[1]],main="")


#####============ 11.3 Multiple Response Multilvel Models ============#####

data(jsp, package="faraway")
jspr <- jsp[jsp$year==2,]
mjspr <- data.frame(rbind(jspr[,1:6],jspr[,1:6]),  subject=factor(rep(c("english","math"),c(953,953))),  score=c(jspr$english/100,jspr$math/40))
ggplot(mjspr, aes(x=raven, y=score))+geom_jitter(alpha=0.25)+facet_grid(gender ~ subject)
mjspr$craven <- mjspr$raven-mean(mjspr$raven)
options(contrasts = c("contr.treatment", "contr.poly"))
mmod <- lmer(score ~ subject*gender + craven*subject + social + (1|school) + (1|school:class) + (1|school:class:id),mjspr)
sumary(mmod, digit=3)
library(pbkrtest)
mmod <- lmer(score ~ subject*gender+craven*subject+social+  (1|school)+(1|school:class)+(1|school:class:id),mjspr, REML=FALSE)
mmodr <- lmer(score ~ subject*gender+craven+subject+social+(1|school)+(1|school:class)+(1|school:class:id),mjspr, REML=FALSE)
KRmodcomp(mmod, mmodr)
0.101^2/(0.101^2+0.117^2)
diagd <- fortify(mmod)
ggplot(diagd, aes(x=.fitted,y=.resid)) + geom_point(alpha=0.3) + geom_hline(yintercept=0) + facet_grid(~ subject) + xlab("Fitted") + ylab("Residuals")






#####==== Chapter 13 Mixed Effect Models for Nonnormal Responses =====#####

#####================ 13.3 Binary Response ================#####
library(faraway)
data(ctsib, package="faraway")
ctsib$stable <- ifelse(ctsib$CTSIB==1,1,0)
xtabs(stable ~ Surface + Vision, ctsib)/80
library(dplyr)
subsum <- ctsib %>% group_by(Subject) %>% 
  summarise(Height=Height[1],  Weight=Weight[1], stable=mean(stable), Age=Age[1], Sex=Sex[1])
library(ggplot2)
ggplot(subsum, aes(x=Height,y=stable))+geom_point()
ggplot(subsum, aes(x=Weight,y=stable))+geom_point()
ggplot(subsum, aes(x=Age,y=stable))+geom_point()
ggplot(subsum, aes(x=Sex,y=stable))+geom_boxplot()

# ignore the subject information entirely
gf <- glm(stable ~ Sex+Age+Height+Weight+Surface+Vision,binomial,data=ctsib)
sumary(gf)
# add a fixed subject effect
gfs <- glm(stable ~ Sex + Age + Height + Weight + Surface + Vision + factor(Subject), 
           binomial,data=ctsib)
# PQL
library(MASS)
modpql <- glmmPQL(stable ~ Sex + Age + Height + Weight + Surface + Vision,  
                  random=~1|Subject,  family=binomial,data=ctsib)
summary(modpql)
# numerical integration
library(lme4)
modlap <- glmer(stable ~ Sex + Age + Height + Weight + Surface +  Vision + 
                  (1|Subject), family=binomial, data=ctsib)
modgh <- glmer(stable ~ Sex + Age + Height + Weight + Surface +  Vision + 
                 (1|Subject), nAGQ=25, family=binomial, data=ctsib)
summary(modgh)
modgh2 <- glmer(stable ~  Surface +  Vision + (1|Subject), nAGQ=25, 
                family=binomial, data=ctsib)
anova(modgh, modgh2)
dd <- fortify(modgh2)
ggplot(dd, aes(sample=.resid))+stat_qq() + facet_grid(Surface~Vision)
# Bayesian method
library(INLA)
formula <- stable ~  Surface +  Vision + f(Subject, model="iid")
result <- inla(formula, family="binomial", data=ctsib)
sigmaalpha <- inla.tmarginal(function(x) 1/sqrt(x), result$marginals.hyperpar$"Precision for Subject")
x <- seq(0,7,length.out = 100)
sdf <- data.frame(yield = x, density=inla.dmarginal(x, sigmaalpha))
ggplot(sdf,aes(x=yield,y=density))+geom_line()
restab <- sapply(result$marginals.fixed, function(x) inla.zmarginal(x,silent=TRUE))
restab <- cbind(restab, inla.zmarginal(sigmaalpha,silent=TRUE))
colnames(restab) = c("mu","norm","dome","open","alpha")
data.frame(restab)
x <- seq(-2,11,length.out = 100)
rden <- sapply(result$marginals.fixed,function(y) inla.dmarginal(x, y))[,-1]
ddf <- data.frame(yield=rep(x,3), density=as.vector(rden), treat=gl(3,100, labels=c("norm","dome","open")))
ggplot(ddf, aes(x=yield, y=density, linetype=treat))+geom_line()
2*inla.pmarginal(0,result$marginals.fixed$Visiondome)
xm <- model.matrix(~ Sex + Age + Height + Weight + Surface + Vision, ctsib)
stabledat <- with(ctsib, list(Nobs=nrow(ctsib),
                              Nsubs=length(unique(ctsib$Subject)), Npreds=ncol(xm),
                              y=stable, subject=Subject, x=xm))
library(rstan)
rt <- stanc("glmmbin.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
fit <- sampling(sm, data=stabledat)
traceplot(fit,pars="sigmasubj", inc_warmup=FALSE)
print(fit,pars=c("sigmasubj","beta"))
ipars <- data.frame(extract(fit, pars=c("sigmasubj","beta")))
colnames(ipars)[-1] <- colnames(xm)
library(reshape2)
rdf <- melt(ipars)
ggplot(rdf, aes(x=value))+geom_density() + facet_wrap(~ variable, scales="free")+geom_vline(xintercept=0)
ppars <- data.frame(extract(fit, pars="subeff"))
sort(colMeans(ppars))

#####================ 13.4 Count Response ================#####
library(faraway)
library(ggplot2)
data(epilepsy, package="faraway")
epilepsy$period <- rep(0:4, 59)
epilepsy$drug <- factor(c("placebo","treatment")[epilepsy$treat+1])
epilepsy$phase <- factor(c("baseline","experiment")[epilepsy$expind+1])
epilepsy[epilepsy$id < 2.5,]
library(dplyr)
epilepsy %>%
  group_by(drug, phase) %>%
  summarise(rate=mean(seizures/timeadj)) %>%
  xtabs(formula=rate ~ phase + drug)
ggplot(epilepsy, aes(x=period, y=seizures,  linetype=drug, group=id)) + 
  geom_line() + xlim(1,4) + scale_y_sqrt(breaks=(0:10)^2) + 
  theme(legend.position = "top", legend.direction = "horizontal")
ratesum <- epilepsy %>%
  group_by(id, phase, drug) %>%
  summarise(rate=mean(seizures/timeadj))
library(tidyr)
comsum <- spread(ratesum, phase, rate)
ggplot(comsum, aes(x=baseline, y=experiment, shape=drug)) + geom_point() + 
  scale_x_sqrt() + scale_y_sqrt() + geom_abline(intercept=0, slope=1) + 
  theme(legend.position = "top", legend.direction = "horizontal")
epilo <- epilepsy %>% filter(id != 49)
modglm <- glm(seizures ~offset(log(timeadj)) + expind + treat + 
                I(expind*treat), family=poisson, data=epilo)
# this is a rate model, so we need to use offset(log) term.
# Q: why use I() here?
sumary(modglm)

library(MASS)
modpql <- glmmPQL(seizures ~offset(log(timeadj)) + expind + treat + 
                    I(expind*treat), random = ~1|id, family=poisson, 
                  data=epilo)
summary(modpql)
library(lme4)
modgh <- glmer(seizures ~offset(log(timeadj)) + expind + treat + 
                 I(expind*treat)+ (1|id), nAGQ=25, family=poisson, data=epilo)
summary(modgh)
exp(-0.302)

epilo$id[epilo$id == 59] <- 49
xm <- model.matrix( ~ expind + treat + I(expind*treat), epilo)
epildat <- with(epilo,list(Nobs=nrow(epilo), Nsubs=length(unique(id)),
                           Npreds=ncol(xm),
                           y=seizures,
                           subject=id,
                           x=xm, offset=timeadj))
library(rstan)
rt <- stanc("glmmpois.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
fit <- sampling(sm, data=epildat)
traceplot(fit,pars="sigmasubj", inc_warmup=FALSE)
ipars <- data.frame(rstan::extract(fit, pars=c("sigmasubj","beta")))
colnames(ipars) <- c("subject","intercept","expind","treat","interaction")
ggplot(ipars, aes(x=subject))+geom_density()
ggplot(ipars, aes(x=interaction))+geom_density() +geom_vline(xintercept=0)
bayespval <- function(x) {p <- mean(x > 0); 2*min(p,1-p)}
smat <-  apply(ipars, 2, function(x) c(mean(x), quantile(x,c(0.025, 0.975)), bayespval(x)))
row.names(smat) <- c("mean","LCB","UCB","pvalue")
t(smat)
formula <- seizures ~offset(log(timeadj)) + expind + treat + I(expind*treat) + f(id,model="iid")
result <- inla(formula, family="poisson", data = epilo)
sigmaalpha <- inla.tmarginal(function(x) 1/sqrt(x), result$marginals.hyperpar$"Precision for id")
restab <- sapply(result$marginals.fixed, function(x) inla.zmarginal(x,silent=TRUE))
restab <- cbind(restab, inla.zmarginal(sigmaalpha,silent=TRUE))
colnames(restab) = c("mu","expind","treat","interaction","alpha")
data.frame(restab)

#####============== 13.5 Generalized Estimating Equations ==============#####
data(ctsib, package="faraway")
ctsib$stable <- ifelse(ctsib$CTSIB==1,1,0)
library(geepack)
options(contrasts=c("contr.treatment", "contr.poly"))
modgeep <- geeglm(stable ~ Sex + Age + Height + Weight + Surface +  Vision, 
                  id=Subject, corstr="exchangeable", scale.fix=TRUE,  
                  data=ctsib, family=binomial)
summary(modgeep)
modgeep2 <- geeglm(stable ~ Sex + Age + Height + Weight +  Surface, 
                   id=Subject, corstr="exchangeable", scale.fix=TRUE,  
                   data=ctsib, family=binomial)
anova(modgeep2, modgeep)
data(epilepsy, package="faraway")
modgeep <- geeglm(seizures ~offset(log(timeadj)) + expind + treat + 
                    I(expind*treat), id=id, family=poisson, corstr="ar1", 
                  data=epilepsy, subset=(id!=49))
summary(modgeep)
