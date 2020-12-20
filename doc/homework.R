## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo=FALSE---------------------------------------------------------------
GPA <- data.frame(
  Class=c("A+","A","A-","B+","B","B-","C+","C","C-","D+","D","D-","F"), 
  GPA= c(4.3,4.0,3.7,3.3,3.0,2.7,2.3,2.0,1.7,1.5,1.3,1.0,0),
  grade= c("100~95","94~90","89~85","84~82","81~78","77~75","74~72","71~68","67~65","64","63~61","60","<60"))
library(kableExtra)
kable(GPA, "html") %>%
  kable_styling(full_width = F,  position = "left")

## -----------------------------------------------------------------------------
set.seed(12)
n <- 1000
a <- 2
b <- 2
u <- runif(n)
x <- b*(1-u)^(-1/a)
f <- function(x) fx <- 8/(x^3)
hist(x,prob=TRUE,breaks=length(x)/10,main = expression(f(x)==8/(x^3)))
y <- seq(2,max(x),0.1)
lines(y,f(y))

## -----------------------------------------------------------------------------
deliver <- function(x,y,z){
  if(abs(z)>=abs(y)&abs(z)>=abs(x)) return(y)
  else return(z)
} # choose u from u1,u2,u3
f_e <- function(x) 0.75*(1-x^2)
set.seed(123)
n <- 1000
u1 <- runif(n,min = -1,max = 1)
u2 <- runif(n,min = -1,max = 1)
u3 <- runif(n,min = -1,max = 1)
u <- numeric(n)
for (i in 1:n) {
  u[i]<- deliver(u1[i],u2[i],u3[i])
}
x <- seq(-1,1,0.01)
hist(u,prob=TRUE,main = expression(f(x)==3/4 (1-x^2)))
lines(x,f_e(x))

## -----------------------------------------------------------------------------
set.seed(313)
n <- 1000
r <- 4
beta <- 2
lambda <- rgamma(n,r,beta)
x <- rexp(n,lambda)
y <- seq(0,max(x),length.out = n)
hist(x,prob=TRUE,breaks=100,main = expression(f(y)==64/ (2+y)^5))
lines(y,64/(2+y)^5)

## -----------------------------------------------------------------------------
set.seed(51)
m <- 1e4
x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi/3
print(c(theta.hat,cos(0)-cos(pi/3)))

## -----------------------------------------------------------------------------
set.seed(57)
m <- 1e4
# Antithetic variate approach function
AV.appro <- function(R){
  u1 <- runif(R/2,0,1)
  u2 <- 1- u1
  u <- c(u1,u2)
  theta.hat.a <- mean(exp(u))
  return(theta.hat.a)
}

# simple MC
simpleMC <- function(R){
  x <- runif(R,0,1)
  theta.hat.s <- mean(exp(x))
  return(theta.hat.s)
}

# the coresponding estimators
print(c(AV.appro(m),simpleMC(m),exp(1)-1))

# compare sd
n <- 1000
MC1 <- numeric(n)
MC2 <- numeric(n)

for (i in 1:n) {
  MC1[i] <- AV.appro(m)
  MC2[i] <- simpleMC(m)
}
print(c(sd(MC1),sd(MC2),sd(MC1)/sd(MC2),1-sd(MC1)/sd(MC2)))

## -----------------------------------------------------------------------------
set.seed(513)
m <- 10000
theta.hat <- se <- numeric(2)

g <- function(x){x^2*exp(-0.5*x^2)/sqrt(2*pi)}

# f_1
x.1 <- rexp(m,1)
x.1 <- x.1+1
fg <- g(x.1) / exp(-x.1+1)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)

#f_2
u <- runif(m)
x.2 <- tan((1+u)*pi/4)
fg <- g(x.2) / (4 / ((1 + x.2^2) * pi))
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)

#result
rbind(theta.hat, se)
x.1 <- sort(x.1)
x.2 <- sort(x.2)
x <- seq(1,10,length.out = m)
y <- g(x)
plot(x,y,type='l')
lines(x.1,exp(-x.1+1),col=2)
lines(x.2,(4 / ((1 + x.2^2) * pi)),col=3)
legend("topright",legend = c('g','f1','f2'),col=1:3,lty=1)

## -----------------------------------------------------------------------------
set.seed(515)
M <- 10000; 
k <- 5 # what if k is larger?
r <- M/k #replicates per stratum

g <- function(x){exp(-x)/(1+x^2)}

f <- function(x,j){
  C <- exp(-0.2*(j-1))-exp(-0.2*j)
  p <- exp(-x)/C
}

X.trans <- function(u,j){
  C <- exp(-0.2*(j-1))-exp(-0.2*j)
  x <- -log(exp(-0.2*(j-1))-C*u)
  return(x)
}

theta.hat <- sd <- numeric(k)

for (j in 1:5) {
  u <- runif(r)
  x <- X.trans(u,j)
  fg <- g(x)/f(x,j)
  theta.hat[j+1] <- mean(fg)
  sd[j+1] <- sd(fg)
}
theta.total <- sum(theta.hat)
sd.total <- sqrt(sum(sd^2))
#names(theta.total)="theta"
#names(sd.total)="sd"
c(theta.total,sd.total)

## -----------------------------------------------------------------------------
mu <- 0
sigma <- 2
m <- 1000
n <- 10

set.seed(64)
mu.hat <- mu.se <- UCL <- numeric(m)

calcCI <- function(n,mu,sigma){
  x <- rlnorm(n,meanlog = 0,sdlog = 2)
  mu.hat <- sum(log(x))/n
  sigma.hat <- sum((log(x)-mu.hat)^2)/(n-1)
  A <- exp(sqrt(sigma.hat))^{qt(0.975,df=n-1)/sqrt(n)}
  CI.left <- exp(mu.hat)/A
  CI.right <- exp(mu.hat)*A
  names(CI.left) <- "Lowbound"
  names(CI.right) <- "Upbound"
  return(c(CI.left,CI.right))
}

for (i in 1:m) {
  result <- calcCI(n,0,2)
  low <- result[1]
  up <- result[2]
  UCL[i] <- low<1&up>1
}
sum(UCL)/m

## -----------------------------------------------------------------------------
mu <- 0
sigma <- 2
m <- 1000
n <- 20

set.seed(651)
mu.hat <- mu.se <- UCL <- numeric(m)

calcCI <- function(n){
  x <- rchisq(n,df=2)
  mu.hat <- sum(log(x))/n
  sigma.hat <- sum((log(x)-mu.hat)^2)/(n-1)
  A <- exp(sqrt(sigma.hat))^{qt(0.975,df=n-1)/sqrt(n)}
  CI.left <- exp(mu.hat)/A
  CI.right <- exp(mu.hat)*A
  names(CI.left) <- "Lowbound"
  names(CI.right) <- "Upbound"
  return(c(CI.left,CI.right))
}

for (i in 1:m) {
  result <- calcCI(n)
  low <- result[1]
  up <- result[2]
  UCL[i] <- low<2&up>2
}
sum(UCL)/m

## -----------------------------------------------------------------------------
set.seed(652)
m <-1000
n <- 20
alpha <- .05
UCL <- numeric(m)
for (i in 1:m) {
  x <- rnorm(n, mean=0, sd=2)
  Upbound <- (n-1) * var(x) / qchisq(alpha,df=n-1)
  UCL[i] <- Upbound>4
}
sum(UCL/m)

## -----------------------------------------------------------------------------
set.seed(67)

nu <- 2
alpha <- 2

n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv <- qnorm(.975, 0, sqrt(6/n)) #crit. values for each n

sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#n is a vector of sample sizes
#we are doing length(n) different simulations
p.reject2 <- p.reject1 <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
sktests2 <- sktests1 <- numeric(m) #test decisions

# beta(alpha,alpha) distribution
for (j in 1:m) {
x1 <- rbeta(n[i],alpha,alpha)
#test decision is 1 (reject) or 0
sktests1[j] <- as.integer(abs(sk(x1)) >= cv[i] )
}
p.reject1[i] <- mean(sktests1) #proportion rejected
}
# t distribution
for (i in 1:length(n)) {
sktests2 <- numeric(m) #test decisions
for (j in 1:m) {
x2 <- rt(n[i],2)
#test decision is 1 (reject) or 0
sktests2[j] <- as.integer(abs(sk(x2)) >= cv[i] )
}
p.reject2[i] <- mean(sktests2) #proportion rejected
}

p.reject1 # for beta(alpha,alpha)
p.reject2 # for t(nu)

## -----------------------------------------------------------------------------
temp <- seq(0,1,length.out = 1000)
pdf2 <- temp*(1-temp)/beta(2,2)
pdf3 <- temp^2*(1-temp)^2/beta(3,3)
pdf4 <- temp^3*(1-temp)^3/beta(4,4)
plot(temp,pdf4,main = "Beta(alpha,alpha)",type = "l",col=4)
lines(temp,pdf3,col=3)
lines(temp,pdf2,col=2)
legend("topright",legend = c('alpha=2','alpha=3','alpha=4'),col=2:4,lty=1)

## -----------------------------------------------------------------------------
# generate samples under H1 to estimate power
set.seed(68)

sigma1 <- 1
sigma2 <- 1.5
n <- c(10, 20, 30, 50, 100, 500) # sample size
m <- 1e4


# The function count5test returns the value 1
# (reject H0) or 0 (do not reject H0)
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

power <- Ftest <- numeric(m)
por.power <- por.Ftest <-  numeric(length(n))
for (j in 1:length(n)) {
  for (i in 1:m) {
  x <- rnorm(n[j], 0, sigma1)
  y <- rnorm(n[j], 0, sigma2)
  power[i] <- count5test(x, y)
  Ftest[i] <- var.test(x,y,conf.level = 0.945)$p.value < 0.055
  }
  por.power[j] <- mean(power)
  por.Ftest[j] <- mean(Ftest)
}
# count5test
por.power
# F test
por.Ftest

## -----------------------------------------------------------------------------
set.seed(666)
n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv.6C <- qchisq(.95, 1)#crit. values for each n


sk.6C <- function(x){
  #computes the sample skewness coeff by  Mardia's approach.
  len <- length(x)
  x<- as.matrix(x)
  xbar <- mean(x)
  m1 <- sum((x-xbar)^2)
  mat <- (x-xbar)%*%t(x-xbar)
  mat <- mat^3
  m2 <- sum(mat)
  return(len*m2*m1^(-3))
}
# repeat example 6.8
p.reject.6C1 <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
sktests.6C1 <- numeric(m) #test decisions
for (j in 1:m) {
x <- rnorm(n[i])
#test decision is 1 (reject) or 0
sktests.6C1[j] <- as.integer(n[i]*sk.6C(x)/6 >= cv.6C )
}
p.reject.6C1[i] <- mean(sktests.6C1) #proportion rejected
}
p.reject.6C1

#repeat example 6.10


alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-alpha, 1)

for (j in 1:N) {   #for each epsilon
e <- epsilon[j]
sktests.6C2 <- numeric(m)
for (i in 1:m) { #for each replicate
sigma <- sample(c(1, 10), replace = TRUE,
size = n, prob = c(1-e, e))
x <- rnorm(n, 0, sigma)
sktests.6C2[i] <- as.integer(n*sk.6C(x)/6 >= cv)
}
pwr[j] <- mean(sktests.6C2)
}
# the proportion rejected
pwr
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
mat <- matrix(c(6510,3490,10000,6760,3240,10000,13270,6730,20000),3,3,dimnames = list(c("Rejected","Accepted","total"),c("A method","B method","total")))
mat

## -----------------------------------------------------------------------------
library(boot)
data(law, package = "bootstrap")
n <- nrow(law)
y <- law$LSAT
z <- law$GPA
mat <- as.matrix(law)
theta.hat <-corr(mat)
print (theta.hat)

#compute the jackknife replicates, leave-one-out estimates
theta.jack <- numeric(n)
for (i in 1:n){
  mat.i <- mat[-i,]
  theta.jack[i] <- corr(mat.i)
}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))

print(c(bias,se)) # jackknife estimate of the bias and the standard error

## -----------------------------------------------------------------------------
library(boot)
data(aircondit, package = "boot")
air.time <- aircondit$hours
theta.boot <- function(x,i) mean(x[i])
boot.obj <- boot(air.time, statistic = theta.boot, R=2000)
print(boot.obj)
print(boot.ci(boot.obj, type = c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
library(bootstrap)
score <- as.matrix(scor)
n <- nrow(score)
#compute the jackknife replicates, leave-one-out estimates
Sigma <- cov(score)
eigens <- eigen(Sigma)$values
theta.hat <- max(eigens)/sum(eigens)
theta.jack <- numeric(n)
for (i in 1:n){
  Sigma.hat.i <- cov(score[-i,])
  eigens.i <- eigen(Sigma.hat.i)$values
  theta.jack[i] <- max(eigens.i)/sum(eigens.i)
}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))

print(theta.hat) # 
print(c(bias,se)) # jackknife estimate of the bias and the standard error


## -----------------------------------------------------------------------------
library(bootstrap)
library(DAAG); attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
# for n-fold cross validation
# fit models on leave-one-out samples

# functions to be cross-validated.
theta.fit.1 <- function(x,y){
  J1 <- lm(y ~ x)
}
theta.predict.1 <- function(fit,x){
  cbind(1,x)%*%fit$coef
}

theta.fit.2 <- function(x,y){
  J2 <- lm(y ~ x + I(x^2))
}
theta.predict.2 <- function(fit,x){
  cbind(1,x,x^2)%*%fit$coef
}

theta.fit.3 <- function(x,y){
  J3 <- lm(log(y) ~ x)
}
theta.predict.3 <- function(fit,x){
  exp(cbind(1,x)%*%fit$coef)
}

theta.fit.4 <- function(x,y){
  J4 <- lm(log(y) ~ log(x))
}
theta.predict.4 <- function(fit,x){
  exp(cbind(1,log(x))%*%fit$coef)
}
# results of CV.
results.1 <- crossval(chemical,magnetic,theta.fit= theta.fit.1,theta.predict = theta.predict.1,ngroup = n/2)
results.2 <- crossval(chemical,magnetic,theta.fit= theta.fit.2,theta.predict = theta.predict.2,ngroup = n/2)
results.3 <- crossval(chemical,magnetic,theta.fit= theta.fit.3,theta.predict = theta.predict.3,ngroup = n/2)
results.4 <- crossval(chemical,magnetic,theta.fit= theta.fit.4,theta.predict = theta.predict.4,ngroup = n/2)
# errors of each method
e1 <- magnetic - results.1$cv.fit
e2 <- magnetic - results.2$cv.fit
e3 <- magnetic - results.3$cv.fit
e4 <- magnetic - results.4$cv.fit
#compare the four methods
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
# L2 is the best model
lm(magnetic ~ chemical + I(chemical^2))

## -----------------------------------------------------------------------------
library(RANN)
library(boot)
library(energy)
library(Ball)
library(kableExtra)

## -----------------------------------------------------------------------------
# experiments for evaluating the performance of the NN,
# energy, and ball methods in various situations.


# statistic for NN
Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}
eqdist.nn <- function(z,sizes,k){
boot.obj <- boot(data=z,statistic=Tn,R=R,
sim = "permutation", sizes = sizes,k=k)
ts <- c(boot.obj$t0,boot.obj$t)
p.value <- mean(ts>=ts[1])
list(statistic=ts[1],p.value=p.value)
}
# settings
m <- 1e3; k<-3; p<-2; mu <- 0.5;
n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

# power comparison 

# case 1: Unequal variances and equal expectations
set.seed(123)
p.values.1 <- matrix(NA,m,3)
for(i in 1:m){
x1 <- matrix(rnorm(n1*p),ncol=p);
y1 <- cbind(rnorm(n2,sd=2),rnorm(n2,sd=2));
z1 <- rbind(x1,y1)
p.values.1[i,1] <- eqdist.nn(z1,N,k)$p.value
p.values.1[i,2] <- eqdist.etest(z1,sizes=N,R=R)$p.value
p.values.1[i,3] <- bd.test(x=x1,y=y1,R=999,seed=i*12345)$p.va
}

# case 2: Unequal variances and Unequal expectations
p.values.2 <- matrix(NA,m,3)
for(i in 1:m){
x2 <- matrix(rnorm(n1*p),ncol=p);
y2 <- cbind(rnorm(n2,mean=3,sd=2),rnorm(n2,mean=3,sd=2));
z2 <- rbind(x2,y2)
p.values.2[i,1] <- eqdist.nn(z2,N,k)$p.value
p.values.2[i,2] <- eqdist.etest(z2,sizes=N,R=R)$p.value
p.values.2[i,3] <- bd.test(x=x2,y=y2,R=999,seed=i*12345)$p.value
}

# case3: Non-normal distributions
# t 
p.values.3.1 <- matrix(NA,m,3)
for(i in 1:m){
x3 <- matrix(rt(n1*p,df=1),ncol=p);
y3<- cbind(rt(n2,df=1),rt(n2,df=1));
z3 <- rbind(x3,y3)
p.values.3.1[i,1] <- eqdist.nn(z3,N,k)$p.value
p.values.3.1[i,2] <- eqdist.etest(z3,sizes=N,R=R)$p.value
p.values.3.1[i,3] <- bd.test(x=x3,y=y3,R=999,seed=i*12345)$p.value
}
# bimodel  X \sim 0.5*N(0,1) + 0.5*N(0,4) Y\sim 0.5*N(0,1) + 0.5* N(0,16)
p.values.3.2 <- matrix(NA,m,3)
temp <- numeric(n)
for(i in 1:m){
u <- rbinom(n,1,0.5)
for (j in 1:n) {
  if(u[j]==0) temp[j] <- rnorm(1)
  else temp[j] <- rnorm(1,mean = 0,sd=2)
}
x3 <- matrix(temp,ncol=p);
v <- rbinom(n,1,0.5)
for (l in 1:n) {
  if(v[l]==0) temp[l] <- rnorm(1)
  else temp[l] <- rnorm(1,mean = 0,sd=4)
}
y3<- matrix(temp,ncol=p);
z3 <- rbind(x3,y3)
p.values.3.2[i,1] <- eqdist.nn(z3,N,k)$p.value
p.values.3.2[i,2] <- eqdist.etest(z3,sizes=N,R=R)$p.value
p.values.3.2[i,3] <- bd.test(x=x3,y=y3,R=999,seed=i*12345)$p.value
}

# case 4: Unbalanced samples- consider 2-calssification problem
p.values.4 <- matrix(NA,m,3)
for(i in 1:m){
x4 <- matrix(rbinom(n,1,10/11),ncol=p)
y4 <- matrix(rbinom(n,1,10/11),ncol=p)
z4 <- rbind(x4,y4)
p.values.4[i,1] <- eqdist.nn(z4,N,k)$p.value
p.values.4[i,2] <- eqdist.etest(z4,sizes=N,R=R)$p.value
p.values.4[i,3] <- bd.test(x=x4,y=y4,R=999,seed=i*12345)$p.value
}
alpha <- 0.1
col_names <- c("NN","energy","Ball")
row_names <- c("case1","case2","case 3-1","case 3-2","case 4")

pow <- matrix(0,5,3,dimnames=list(row_names,col_names))
pow[1,] <- colMeans(p.values.1<alpha)
pow[2,] <- colMeans(p.values.2<alpha)
pow[3,] <- colMeans(p.values.3.1<alpha)
pow[4,] <- colMeans(p.values.3.2<alpha)
pow[5,] <- colMeans(p.values.4<alpha)


kable(pow, "html") %>%
  kable_styling(full_width = F,  position = "left")

## -----------------------------------------------------------------------------
set.seed(1)
denLaplace <- function(x){
  return(0.5*exp(-abs(x)))
}

RW.MSample <- function(N,sigma,x0) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (denLaplace(y) / denLaplace(x[i-1]))){
      x[i] <- y
    }
    else {
      x[i] <- x[i-1]
      k <- k + 1
      } 
    }
  return(list(x=x, k=k))
}

# start
N <- 2000
sigma <- c(.05, .5, 2, 10)
x0 <- 25
RW1 <- RW.MSample(N,sigma[1], x0)
RW2 <- RW.MSample(N,sigma[2], x0)
RW3 <- RW.MSample(N,sigma[3], x0)
RW4 <- RW.MSample(N,sigma[4], x0)

print(c((2000-RW1$k)/2000, (2000-RW2$k)/2000, (2000-RW4$k)/2000, (2000-RW4$k)/2000))

## -----------------------------------------------------------------------------
# set.seed(11)
# Gelman.Rubin <- function(psi) {
# # psi[i,j] is the statistic psi(X[i,1:j])
# # for chain in i-th row of X
#   psi <- as.matrix(psi)
#   n <- ncol(psi)
#   psi.means <- rowMeans(psi) #row means
#   B <- n * var(psi.means) #between variance est.
#   psi.w <- apply(psi, 1, "var") #within variances
#   W <- mean(psi.w) #within est.
#   v.hat <- W*(n-1)/n + (B/n) #upper variance est.
#   r.hat <- v.hat / W #G-R statistic
#   return(r.hat)
# }
# denLaplace<-function(x){
#   return(0.5*exp(-abs(x)))
# }
# 
# RW.MSample <- function(N,sigma,x0) {
#   x <- numeric(N)
#   x[1] <- x0
#   u <- runif(N)
#   k <- 0
#   for (i in 2:N) {
#     y <- rnorm(1, x[i-1], sigma)
#     if (u[i] <= (denLaplace(y) / denLaplace(x[i-1]))){
#       x[i] <- y
#     }
#     else {
#       x[i] <- x[i-1]
#       k <- k + 1
#       } 
#     }
#   return(list(x=x, k=k))
# }
# 
# 
# sigma <- sqrt(2) #parameter of proposal distribution
# k <- 4 #number of chains to generate
# n <- 10000 #length of chains
# b <- 500 #burn-in length
# 
# x0 <- c(-10, -5, 5, 10)
# 
# #generate the chains
# X <- matrix(0, nrow=k, ncol=n)
# for (i in 1:k)
#   X[i, ] <- RW.MSample(n, sigma, x0[i])$x
# 
# #compute diagnostic statistics
# psi <- t(apply(X, 1, cumsum))
# for (i in 1:nrow(psi))
#   psi[i,] <- psi[i,] / (1:ncol(psi))
# print(Gelman.Rubin(psi))
# # plot
# par(mfrow=c(2,2))
# for (i in 1:k)
#   plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
k<-c(4:25,100,500,1000)
root<-numeric(length(k))
for (i in 1:length(k)) {
  res <- uniroot(function(a){
    pt(sqrt(a^2*(k[i]-1)/(k[i]-a^2)),df=k[i]-1,log.p = T)-pt(sqrt(a^2*(k[i])/(k[i]+1-a^2)),df=k[i],log.p = T)
  },lower = 1e-5,upper = sqrt(k[i]-1e-5))
  root[i]<-unlist(res)[[1]]
}
root


## -----------------------------------------------------------------------------
# complete data likelihood
Clikelihood <- function(pq,nAA,nBB){
  # pq is a vector contains p and q
  # data
  p<- pq[1]
  q<- pq[2]
  r <- 1-p-q
  
  nA. <- 444
  nB. <- 132
  nOO <- 361
  nAB <- 63
  nAO <- nA. - nAA
  nBO <- nB. - nBB
  
  fx1 <- 2*nAA*log(p)+2*nBB*log(q)+2*nOO*log(r)
  fx2 <- nAO*log(2*p*r)+nBO*log(2*q*r)+nAB*log(2*p*q)
  fx <- fx1 + fx2
  return(-fx)
}

# gradient for llikelihood
DLikelihood <- function(pq,nAA,nBB){
  p<- pq[1]
  q<- pq[2]
  r <- 1-p-q
  
  nA. <- 444
  nB. <- 132
  nOO <- 361
  nAB <- 63
  nAO <- nA. - nAA
  nBO <- nB. - nBB
  
  Dfx1 <- (2*nAA+nAO+nAB)/p - (2*nOO+nAO+nBO)/r
  Dfx2 <- (2*nBB+nBO+nAB)/q - (2*nOO+nAO+nBO)/r
  Dfx  <- c(Dfx1, Dfx2)
  return(-Dfx)
}
# calculate M-step
Mstep <- function(pq){
  p <- pq[1]
  q <- pq[2]
  r <- 1-p-q
  nA. <- 444
  nB. <- 132
  p0 <- p^2/(p^2 + 2*p*r)
  q0 <- q^2/(q^2 + 2*q*r)
  x <- nA.*(1-p0)/(1+p0)
  y <- nB.*(1-q0)/(1+q0)
  return(c(x,y))
}


# initial value
pq0 <- c(0.4,0.2)

Mtemp <- Mstep(pq0)
nAA <- Mtemp[1]
nBB <- Mtemp[2]

Ltemp <- optim(pq0,Clikelihood,gr=DLikelihood,nAA=nAA,nBB=nBB)
pq1 <- Ltemp[[1]]

# star
while (sum(abs(pq1-pq0))>1e-5) {
  pq0 <- pq1
  Mtemp <- Mstep(pq0)
  nAA <- Mtemp[1]
  nBB <- Mtemp[2]
  Ltemp <- optim(pq0,Clikelihood,gr=DLikelihood,nAA=nAA,nBB=nBB)
  pq1 <- Ltemp[[1]]
  like <- Clikelihood(pq1,nAA,nBB)
  print(-like)
}
# the  are shown as follow.

## -----------------------------------------------------------------------------
data <- mtcars

formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

looplist <- vector("list",4)
for (i in 1:length(formulas)) {
  looplist[[i]] <- lm(formulas[[i]],data)
}
looplist
lapplylist <- lapply(formulas, lm,data=data)
lapplylist

## -----------------------------------------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)

sapply(trials, function(i) i$p.value)

## -----------------------------------------------------------------------------
testlist <- list(iris, mtcars, cars)
lapply(testlist, function(x) vapply(x, mean, numeric(1)))


## -----------------------------------------------------------------------------
mlapply <- function(X, FUN, FUN.VALUE, simplify = FALSE){
  out <- Map(function(x) vapply(x, FUN, FUN.VALUE), X)
  if(simplify == TRUE){return(simplify2array(out))}
  out
}

mlapply(testlist, mean, numeric(1))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

## -----------------------------------------------------------------------------
# library(StatComp20088)
# #library(Rcpp) 
# #sourceCpp('RW.cpp')
# #test: rw=RWcpp(2,25,2000)
# set.seed(123)
# N <- 2000
# sigma <- c(.05, .5, 5, 50)
# x0 <- 25
# rw1_c <- RWcpp(sigma[1],x0,N)
# rw2_c <- RWcpp(sigma[2],x0,N)
# rw3_c <- RWcpp(sigma[3],x0,N)
# rw4_c <- RWcpp(sigma[4],x0,N)
# #number of candidate points rejected
# reject <- cbind(rw1_c$k, rw2_c$k, rw3_c$k, rw4_c$k)
# accept_rates <- round((N-reject)/N,4)
# sigma <- matrix(sigma,nrow = 1)
# mat <- rbind(sigma,accept_rates)
# rownames(mat) <- c("sigma","Acceptance rates")
# colnames(mat) <- c("1","2","3","4")
# knitr::kable(mat)
# 
# #plot
# par(mfrow = c(2,2))  #display 4 graphs together
# rw <- cbind(rw1_c$x,rw2_c$x, rw3_c$x,rw4_c$x)
# for (j in 1:4) {
#   plot(rw[,j], type="l",
#   xlab=bquote(sigma == .(round(sigma[j],3))),
#   ylab="X", ylim=range(rw[,j]))
#}

## -----------------------------------------------------------------------------
# # R function for Ex. 9.4
# denLaplace <- function(x){
#   return(0.5*exp(-abs(x)))
# }
# 
# RW.MSample <- function(N,sigma,x0) {
#   x <- numeric(N)
#   x[1] <- x0
#   u <- runif(N)
#   k <- 0
#   for (i in 2:N) {
#     y <- rnorm(1, x[i-1], sigma)
#     if (u[i] <= (denLaplace(y) / denLaplace(x[i-1]))){
#       x[i] <- y
#     }
#     else {
#       x[i] <- x[i-1]
#       k <- k + 1
#       } 
#     }
#   return(list(x=x, k=k))
# }
# # by R function
# set.seed(12)
# N <- 2000
# sigma <- 1  # N(0,1)
# x0 <- 0
# rw_R <- RW.MSample(N,sigma, x0)
# rw_R <- rw_R$x
# rw_R <- rw_R[501:N]
# # by Rcpp 
# rw_cpp <- RWcpp(1,x0,N)
# rw_cpp <- rw_cpp$x
# rw_cpp <- rw_cpp[501:N]
# #plot qqplot
# Giao <- ppoints(500)
# Quan_R <- quantile(rw_R,Giao)
# Quan_cpp <- quantile(rw_cpp,Giao)
# qqplot(Quan_R,Quan_cpp,main="",xlab="rw3 quantiles",ylab="rw3_c quantiles")
# qqline(Quan_cpp)

## -----------------------------------------------------------------------------
library(microbenchmark)
# microbenchmark(
#   RW.MSample(N,1, x0),
#   RWcpp(1,x0,N))

