library("MASS")

#' @title inverse of symmetry matrix
#' @description A function to calculating the inverse of a symmetry matrix
#' @param smatrix the symmetry matrix
#' @param sB sB
#' @param det determinant
#' @param log logarithm
#' @param jitter jiiter to aviod singular case
#' @return return the inverse
#' @export
mymatinv=function (smatrix, sB = "sB", det = F, log = T, jitter = 1e-10)
{
  mat = smatrix + diag(jitter, dim(smatrix)[1])
  smatrix = mat
  if (is.character(sB))
    sB = diag(1, dim(mat)[1])
  else sB = as.matrix(sB)
  x = solve(smatrix, sB)
  d = NULL
  if (det == T) {
    L = chol(smatrix)
    if (log == T)
      d = 2 * log(prod(diag(L)))
    if (log == F)
      d = prod(diag(L))
  }
  return(list(res = as.matrix(x), det = d))
}

#' @title Covariance matrix
#' @description A function to calculating the covariance matrix
#' @param hyper the hyper parameters involved in covariance matrix
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @return return covariance matrix
#' @export
cov.func=function(hyper,x1,x2){
  n <- length(x1)
  #x <- c(x[u==0],x[u>0])
  x1.mat <- matrix(rep(x1,n),ncol=n)
  X1 <- (x1.mat - t(x1.mat))^2
  X1 <- -0.5*exp(hyper[2])*X1
  x2.mat <- matrix(rep(x2,n),ncol=n)
  X2 <- (x2.mat - t(x2.mat))^2
  X2 <- -0.5*exp(hyper[3])*X2
  X <- X1 + X2
  Cov.K <- exp(hyper[1])*exp(X)  + diag(1e-4,dim(X)[1])
  Cov.K <- as.matrix(Cov.K)
  return(Cov.K)
}

#' @title Log-likelihood
#' @description A function to calculating the minus log-likelihood
#' @param hyper the hyper parameters involved in covariance matrix
#' @param f.hat the latent process
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param y the response variable
#' @return return minus log-likelihood value
#' @export
y.loglike <- function(hyper,f.hat,x1,x2,y){
  n <- length(x1)
  temp <- lapply(f.hat, exp)
  temp <- as.numeric(temp)
  temp1 <- 1/(1 + temp)
  temp2 <- 1 - temp1
  A <- diag(as.vector(temp1*temp2),n)
  K <- cov.func(hyper,x1,x2)
  invK <- mymatinv(K)$res
  squ <- t(f.hat)%*%invK%*%f.hat
  l <- 0.5*log(det(K)) + 0.5*squ[1,1]+0.5*log(det(invK + A))
  return(l)
}

#' @title The gradient for log-likelihood
#' @description A function to calculating the gradient of log-likelihood respect to the hyper-paramters
#' @param hyper the hyper parameters involved in covariance matrix
#' @param f.hat the latent process
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param y the response variable
#' @return return gradient of log-likelihood respect to the hyper-paramters
#' @export
y.Dloglike <- function(hyper,f.hat,x1,x2,y){
  n <- length(x1)
  temp <- lapply(f.hat, exp)
  temp <- as.numeric(temp)
  temp1 <- 1/(1 + temp)
  temp2 <- 1 - temp1
  A <- diag(as.vector(temp1*temp2),n)
  K <- cov.func(hyper,x1,x2)
  invK <- mymatinv(K)$res
  B <- invK + mymatinv(invK + A)$res
  # partial v
  DK.v <-  K/exp(hyper[1])
  squ.v <- t(f.hat)%*%invK%*%DK.v%*%invK%*%f.hat
  Dlike.v <- 0.5*sum(diag(B%*%DK.v)) - 0.5*squ.v[1,1]
  # partial w1
  x1.mat <- matrix(rep(x1,n),ncol=n)
  X1 <- (x1.mat - t(x1.mat))^2
  DK.w1 <- -0.5*X1*K
  squ.w1 <- t(f.hat)%*%invK%*%DK.w1%*%invK%*%f.hat
  Dlike.w1 <- 0.5*sum(diag(B%*%DK.w1)) - 0.5*squ.w1[1,1]
  # partial w2
  x2.mat <- matrix(rep(x2,n),ncol=n)
  X2 <- (x2.mat - t(x2.mat))^2
  DK.w2 <- -0.5*X2*K
  squ.w2 <- t(f.hat)%*%invK%*%DK.w2%*%invK%*%f.hat
  Dlike.w2 <- 0.5*sum(diag(B%*%DK.w2)) - 0.5*squ.w2[1,1]
  # gradiant
  Dlike <- c(Dlike.v,Dlike.w1,Dlike.w2)
  return(Dlike)
}

#' @title  Laplace approximation
#' @description A function to approximate the marginal log-likelihood by using Laplace method
#' @param f the latent process f
#' @param hyper the hyper parameters involved in covariance matrix
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param y the response variable
#' @return return the gradient
#' @export
gamma.f <- function(f,hyper,x1,x2,y){
  n <- length(y)
  temp1 <- lapply(f*y, exp)
  temp1 <- as.numeric(temp1)
  temp2 <- lapply(f, exp)
  temp2 <- 1 + as.numeric(temp2)
  temp <- lapply(temp1/temp2,log)
  temp <- as.numeric(temp)
  part1 <- sum(temp)
  K <- cov.func(hyper,x1,x2)
  invK <- mymatinv(K)$res
  squ <- t(f)%*%invK%*%f
  part2 <- -0.5*log(det(K))-0.5*squ[1,1]-0.5*n*log(2*pi)
  g <- part1 + part2
  return(-g)
}

#' @title The gradient for Laplace approximation
#' @description A function to calculating the gradient for Laplace approximation
#' @param f the latent process f
#' @param hyper the hyper parameters involved in covariance matrix
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param y the response variable
#' @return return the gradient
#' @export
Dgamma.f <- function(f,hyper,x1,x2,y){
  temp <- lapply(f, exp)
  temp <- 1 + as.numeric(temp)
  temp.pi <- 1 - 1/temp
  K <- cov.func(hyper,x1,x2)
  invK <- mymatinv(K)$res
  # gradient
  D <- -y + temp.pi + invK%*%f
  return(as.vector(D))
}

#' @title The main function for fitting a binary Generalized GPR
#' @description A function for fitting a binary Generalized GPR and return the estimated kernel hyper-parameters
#' @param hyper the hyper parameters involved in covariance matrix
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param y the response variable
#' @return a list contains hyper-parameters and the latent process f
#' @export
fit.GPbinary <- function(hyper,x1,x2,y){
  hyper.0 <- hyper
  n <- length(x1)
  mu.0 <-numeric(n)
  I <- diag(n)
  f.initial <- mvrnorm(1,mu = mu.0,Sigma = I)
  f.0 <- nlminb(f.initial,gamma.f,Dgamma.f,hyper = hyper.0,x1=x1,x2=x2,y=y)
  f.0 <- f.0[[1]]
  hyper.1 <- nlminb(hyper.0,y.loglike,y.Dloglike,f.hat = f.0,x1=x1,x2=x2,y=y)
  hyper.1 <- hyper.1[[1]]
  f.1 <- nlminb(f.0,gamma.f,Dgamma.f,hyper = hyper.1,x1=x1,x2=x2,y=y)
  f.1 <- f.1[[1]]
  while (sum(abs(hyper.1 -hyper.0)) + sum(abs(f.1- f.0)) >1e-6) {
    f.0 <- f.1
    hyper.0 <- hyper.1
    hyper.1 <- nlminb(hyper.0,y.loglike,y.Dloglike,f.hat = f.0,x1=x1,x2=x2,y=y)
    hyper.1 <- hyper.1[[1]]
    f.1 <- nlminb(f.0,gamma.f,Dgamma.f,hyper = hyper.1,x1=x1,x2=x2,y=y)
    f.1 <- f.1[[1]]
  }
  result <- list('hyper'=hyper.1, 'f'=f.1)
  return(result)
}

#' @title The prediction function
#' @description A function for predicting at a new input after we have fitted the binary generalized Gaussian process regression model with \code{fit.GPbinary()}
#' @param hyper the hyper parameters obtained from fit.GPbinary()
#' @param x1 the fitst covariate
#' @param x2 the second covariate
#' @param x1.new the first new covariate
#' @param x2.new the second new covariate
#' @param f.hat the latent process obtained from fit.GPbinary()
#' @return return a prediction
#' @export
GP.predict <- function(hyper,x1,x2,x1.new,x2.new,f.hat){
  n <- length(x1)
  K <- cov.func(hyper,x1,x2)
  invK <- mymatinv(K)$res
  k <- numeric(n)
  for (i in 1:n) {
    temp <- exp(hyper[2])*(x1.new-x1[i])**2 + exp(hyper[3])*(x2.new - x2[i])**2
    k[i] <- exp(hyper[1])*exp(-0.5*temp)
  }
  f.new.hat <- t(k)%*%invK%*%f.hat
  f.new.hat <- f.new.hat[1,1]
  return(f.new.hat)
}
