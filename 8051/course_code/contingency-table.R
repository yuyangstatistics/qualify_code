library(faraway)

# wafer data

y<-c(320,14,80,36)
particle<-gl(2,1,4,labels=c("no","yes")) # gl generates factor levels; the second argument indicates number of repeats
quality<-gl(2,2,4,labels=c("good","bad"))
particle
quality
wafer<-data.frame(y,particle,quality)
wafer

80/320    # odds of "bad" without particle
36/14     # odds of "bad" with particle
(36/14)/(80/320)  # odds ratio

# binomial analysis (400 wafers with no particles, 50 with)


bd<-matrix(y,nrow=2)  # the default is byr=F.
bd

fisher.test(bd) # odds ratio is y11*y22/(y12*y21)

fisher.test(bd,alternative = "greater")


# effects of smoking

data(femsmoke)
femsmoke

ct<-xtabs(y~smoker+dead,femsmoke)
ct
a <- prop.table(ct,1)

summary(ct)  # gives chisquare test of independence
fisher.test(ct)  

prop.table(xtabs(y~smoker+age,femsmoke),2) # column sums to 1
xtabs(y~smoker+age,femsmoke)
prop.table(xtabs(y~smoker+age,femsmoke),1) # row sums to 1
xtabs(y~smoker+age,femsmoke)

ct1<-xtabs(y~smoker+dead,femsmoke,subset=(age=="18-24"))
ct1
prop.table(ct1,1)

ct2<-xtabs(y~smoker+dead,femsmoke,subset=(age=="25-34"))
ct2
prop.table(ct2,1)

ct3<-xtabs(y~smoker+dead+age,femsmoke)
ct3
apply(ct3,3,function(x) (x[1,1]*x[2,2]/(x[1,2]*x[2,1]))) # observed odds ratio at the different age groups

ct4<-xtabs(y~smoker+dead,femsmoke,subset=(age=="55-64"))
ct4
prop.table(ct4,1)

m<-glm(ct4~factor(c(0,1)),family=binomial)
summary(m)
drop1(m)
summary(ct4) # Pearson chisquare gives similar result on independence 

mantelhaen.test(ct3,exact=TRUE) # test conditional independence in each stratum 
# given age, smoker is correlated with dead

summary(ct3) # mutual independence (not pairwise independence)
# mutual independence: P(ABC) = P(A)P(B)P(C)

modm<-glm(y~smoker+dead+age,femsmoke,family=poisson) # mutual independence model
summary(modm)

modj<-glm(y~smoker*dead+age,femsmoke,family=poisson) # joint independence model
summary(modj)
# joint independece means smoker*dead as a group is independent with age

modc<-glm(y~smoker*age+age*dead,femsmoke,family=poisson) # conditional independence
summary(modc)
# conditional independence means given age, smoker has no interaction with dead,
# namely, given age, the other two factors are independent.

modu<-glm(y~(smoker+age+dead)^2,femsmoke,family=poisson) # uniform association
summary(modu)
# for conditional independence, given age, the odds ratio are all 1; 
# for uniform association, given age, the odds ratio are all the same, but not 
# necessarily 1. In such kind of model, it considers all interaction except the 
# three-factor interaction. With three-factor interaction, the model will become
# the saturated model.

anova(modc,modu,test="Chisq") # modc is not preferred

ctf<-xtabs(fitted(modu)~smoker+dead+age,femsmoke)
apply(ctf,3,function(x) (x[1,1]*x[2,2]/(x[1,2]*x[2,1])))

modsat<-glm(y~smoker*dead*age,femsmoke,family=poisson)
summary(modsat)
step(modsat) # Use AIC to do stepwise selection

# binomial regression when one variable is treated as the response
femsmoke 

ybin<-matrix(femsmoke$y,ncol=2) # create binomial data treating life status as the response

summary(modbin<-glm(ybin~smoker*age,femsmoke[1:14,],family=binomial))
drop1(modbin,test="Chi")
modbinr<-glm(ybin~smoker+age,femsmoke[1:14,],family=binomial) # Note this model corresponds to the Poisson regression with uniform association

c(deviance(modu),deviance(modbinr))

summary(modbina<-glm(ybin~age,femsmoke[1:14,],family=binomial)) # Note this model corresponds to the Poisson regression with conditional independence (modc)
c(deviance(modc),deviance(modbina))

summary(modbins<-glm(ybin~smoker,femsmoke[1:14,],family=binomial)) # Note this model does not correspond to the Poisson regression with joint independence (modj)
c(deviance(modj),deviance(modbins))


summary(modbinn<-glm(ybin~1,femsmoke[1:14,],family=binomial)) # Note this model does not corresponds to the Poisson regression with mutual independence (modc)
c(deviance(modm),deviance(modbinn))  # not equal


# an example of ordinal variables
data(nes96)
help(nes96)
nes96
xtabs(~PID+educ,nes96)
partyed<-as.data.frame.table(xtabs(~PID+educ,nes96))
partyed

modno<-glm(Freq~PID+educ,partyed,family=poisson) # a nominal model
summary(modno)
pchisq(40.743,36,lower=FALSE)  

(partyed$oPID<-unclass(partyed$PID))
(partyed$oeduc<-unclass(partyed$educ))

modod<-glm(Freq~PID+educ+I(oPID*oeduc),partyed,family=poisson) # linear-by-linear association model
summary(modod)

anova(modno,modod,test="Chi")

modod2<-glm(Freq~PID+educ+oPID:educ,partyed,family=poisson) # an ordinal-by-nominal model
summary(modod2)
anova(modod,modod2,test="Chi") # linear-by-linear association model is preferred