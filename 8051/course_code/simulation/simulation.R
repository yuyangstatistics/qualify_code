setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
seed <- 999
n = 500
# in total 30 variables
p = 30
p0 = 10 # important variables are in the first 10 variables
set.seed(seed)
x = matrix(rnorm(n * p), n, p)
# add collinearity to the features
corr.seq <- c(0.6, 0.7, 0.8, 0.9, 0.96, 0.8, 0.89, 0.9, 0.99, 0.9)
for (i in 1:p0) {
  set.seed(i)
  x[, i + 20] <- x[, i + 10] * corr.seq[i] + rnorm(n) * 0.5
}
# add scale to the features
set.seed(seed+9)
scale.seq <- rnorm(p, mean = 0, sd = 10)
mean.seq <- rnorm(p, mean = 0, sd = 30)
for (i in (p0+1):p) {
  x[, i] <- x[, i] * scale.seq[i] + mean.seq[i]
}
# 6 linear terms
set.seed(seed)
lterms <- sample(1:p0, 6)
set.seed(seed+1)
lbeta = rnorm(6) * 10
# 1 pure quadratic term
set.seed(seed+2)
pqterms <- sample(lterms, 1)
set.seed(seed+3)
qbeta <- rnorm(1) * 5
# 1 interaction terms
set.seed(seed+4)
iterms <- sample(lterms, 2)
set.seed(seed+5)
ibeta <- rnorm(1) * 5

set.seed(seed+6)
fx = x[, lterms] %*% lbeta + x[, pqterms]^2 * qbeta + 
  x[, iterms[1]] * x[, iterms[2]] * ibeta + rnorm(1, mean=10)

px = exp(fx)
px = px/(1 + px)
set.seed(seed+6)
y = rbinom(n = length(px), prob = px, size = 1)
simulated <- data.frame(cbind(x, y))
set.seed(seed + 7)
simulated <- simulated[sample(nrow(simulated)), ]
y.test <- simulated$y[301:500]
simulated$y[301:500] <- NA
head(simulated)
save(simulated, file = 'simulated.rda')
save(y.test, file = "ytest.rda")




###### check my results ########
load('data/simulate_prediction.rda')
load('ytest.rda')
mean(pred == y.test)
table(pred, y.test)

lterms
lbeta
pqterms
qbeta
iterms
ibeta