library(alr4)
data(turk0)
n1 <- nls(Gain ~ th1 + th2*(1 - exp(-th3 * A)), data=turk0,
          start=list(th1=620, th2=200, th3=10))
summary(n1)

plot(Gain ~ A, turk0, xlab="Amount (percent of diet)", ylab="Weight gain, g")
x <- (0:44)/100
lines(x, predict(n1, data.frame(A=x)))

m1 <- lm(Gain ~ as.factor(A), turk0) # the one-way ANOVA model has no bias
anova(n1,m1)    # lack-of-fit test

n2 <- nls(Gain ~ th1 + 150*(1 - exp(-th3 * A)), data=turk0,
          start=list(th1=620, th3=10))

anova(n2,n1)    # compare the two nested models

# May need to create indicator variables to deal with factors

data(turkey)

temp <- model.matrix( ~ as.factor(S) -1, turkey) # "-1" drops intercept
colnames(temp) <- c("S1", "S2", "S3")
turkey <- cbind(turkey, temp)  # the modified data matrix with the indicators

m4 <- nls( Gain ~ (th1 + th2*(1-exp(-th3*A))),
           data=turkey,start=list(th1=620, th2=200, th3=10),weights=m)
# most general
m1 <- nls( Gain ~ S1*(th11 + th21*(1-exp(-th31*A)))+
             S2*(th12 + th22*(1-exp(-th32*A)))+
             S3*(th13 + th23*(1-exp(-th33*A))),data=turkey,
           start= list(th11=620, th12=620, th13=620,
                       th21=200, th22=200, th23=200,th31=10, th32=10, th33=10),weights=m)
# common intercept
m2 <- nls(Gain ~ th1 + S1*(th21*(1-exp(-th31*A)))+
            S2*(th22*(1-exp(-th32*A)))+S3*(th23*(1-exp(-th33*A))),
          data=turkey,start= list(th1=620,th21=200, th22=200, th23=200,
                                  th31=10, th32=10, th33=10),weights=m)
# common intercept and asymptote
m3 <- nls( Gain ~ (th1 + th2 *(S1*(1-exp(-th31*A))+S2*(1-exp(-th32*A))+
                                 S3*(1-exp(-th33*A)))),data=turkey,start= list(th1=620, th2=200,
                                                                               th31=10,th32=10,th33=10),weights=m)

anova(m4, m3, m2 ,m1)
anova(m4,m1)

# boostrap
library(car)
data(segreg)
s1 <- nls(C ~ th0 + th1 * (pmax(0, Temp - gamma)), data=segreg,
          start=list(th0=70, th1=.5, gamma=40))
summary(s1)

s1.boot <- Boot(s1, R=999)
summary(s1.boot)
confint(s1.boot,type="bca")  # The intervals calculated using the adjusted bootstrap percentile (BCa) method. Other choices are "norm", "basic", "stud", and "perc".



