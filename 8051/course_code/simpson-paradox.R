y<-c(360,30,310,390)
n<-c(900,100,450,600)
tr<-c(1,2,1,2) # 1 = 3 mammograms; 2 =  1 mammogram
h<-c(1,1,2,2)  # 1 = With family history; 2 = Without family history
demo<-cbind(Survival=y,Total=n,Treatment=tr,History=h)
xtabs(y~tr,demo)/xtabs(n~tr,demo)
xtabs(y/n~tr+h,demo)
m<-glm(cbind(y,n-y)~tr,family=binomial)
summary(m)


m1<-glm(cbind(y,n-y)~tr+h,family=binomial)
summary(m1) # Missing a lurking variable gives a totally wrong conclusion! The lack of fit is seen from summary(m). But if the history variable is not in the data, no lack of fit can be detected. The new data: y<-c(670,420); n<-c(1350,700); tr<-c(1,2); mm<-glm(cbind(y,n-y)~tr,family=binomial); summary(mm)

mm<-glm(cbind(y,n-y)~tr+h,family=binomial(link=probit))
summary(mm) 


# How did the paradox occur?

a<-xtabs(n~h+tr,demo)

a/(a%*%matrix(1,2,2))     # Reason: Most women with family history went with 3 mammograms; the problme is avoided if women are randomly assigned with 3 mammograms or 1 mammogram. The commond is the same as:  
prop.table(a,1)

m2<-glm(cbind(y,n-y)~tr*h,family=binomial) # no df left.
summary(m2)

# We re-organize the above data and run Poisson regressions

z<-c(360,30,310,390,540,70,140,210)
tr<-c(1,2,1,2,1,2,1,2)
h<-c(1,1,2,2,1,1,2,2)
s<-c(1,1,1,1,2,2,2,2) # 1 = survival; 2 = death

g<-glm(z~s*tr+h*tr,family=poisson)
summary(g) 
summary(m) # Note that the deviances of m and g are the same, though the null deviances are different. Also, the coefficient for s:tr in g is the same as that for tr in m, up to sign change. If we use other link than logit, the deviances are not the same (for cbind(y,n-y)~tr, the deviance does not depend on the choice of the link).

g1<-glm(z~s*tr+s*h+tr*h,family=poisson)
summary(g1)
summary(m1)

g2<-glm(z~s*tr*h,family=poisson)
summary(g2)
summary(m2)


summary(glm(cbind(y,n-y)~1,family=binomial))
summary(glm(z~s+tr*h,family=poisson))



