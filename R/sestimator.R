##Functions for Matlab Sm.m function:

##MADNOC

madnoc <- function(y){
  y <- abs(y)
  if (is.matrix(y)==0) {
    y <- t(as.matrix(y))}
  n <- nrow(y)
  p <- ncol(y)
  if((n==1)&(p>1)){
    y <- t(y)
    n <- nrow(y)
    p <- ncol(y)
  }
  y <- apply(y, 2, sort)                 
  if (floor(n/2)==n/2){
    mad <- (y[n/2,]+y[(n+2)/2,])/2}
  else {
    mad <- y[(n+1)/2,]
  }
  mad <- mad/.6745
  return(mad)
}

##--------------------------------------------------------------------
##MAHALANOBIS
##build in function in S:

##mahalanobis(x, center, cov, inverted=F) 

##--------------------------------------------------------------------
##KSIINT

ksiint <- function(c, s, p){
  (2^s)*gamma(s+p/2)*pgamma(c^2/2, s+p/2)/gamma(p/2)
}
##--------------------------------------------------------------------
## Computes Tukey's biweight psi function with constant c for all 
## values in the vector x.

psibiweight <- function(x, c){
  hulp <- x-2*x^3/(c^2)+x^5/(c^4)
  psi <- hulp*(abs(x)<c)
  return(psi)
}

##--------------------------------------------------------------------
##UNIRAN
## The random generator.

uniran <- function(seed){
  seed <- floor(seed*5761)+999
  quot <- floor(seed/65536)
  seed <- floor(seed)-floor(quot*65536)
  random <- seed/65536     ##.D0???
  return(list(seed, random))
}

##--------------------------------------------------------------------
##RANDOMSET

## This function is called if not all (p+1)-subsets out of n will be 
## considered. It randomly draws a subsample of nel cases out of tot. 

## tot must be greater than nel!!

randomset <- function(tot, nel, seed){
  ranset <- rep(0, nel)
  for (j in 1:nel){
    random <- uniran(seed)[[2]]
    seed <- uniran(seed)[[1]]
    num <- floor(random*tot)+1
    if (j>1){    
      while (sum(ranset==num)>0){
        random <- uniran(seed)[[2]]
        seed <- uniran(seed)[[1]]
        num <- floor(random*tot)+1
      }
    }
    ranset[j:nel] <- num
  }
  return(list(ranset, seed))
}
##--------------------------------------------------------------------
##RHOBIWEIGHT

rhobiweight <- function(x, c){
  hulp <- x^2/2-x^4/(2*c^2)+x^6/(6*c^4)
  rho <- hulp*(abs(x)<c)+c^2/6*(abs(x)>=c)
  return(rho)
}

##--------------------------------------------------------------------
##TBSB

Tbsb <- function(c, p){
  y1 <- ksiint(c, 1, p)*3/c-ksiint(c, 2, p)*3/(c^3)+ksiint(c, 3, p)/(c^5)
  y2 <- c*(1-pchisq(c^2, p))
  res <- y1+y2
  return(res)
}

##--------------------------------------------------------------------
##TBSC

## constant for Tukey Biweight S 

Tbsc <- function(alpha, p){
  talpha <-  sqrt(qchisq(1-alpha, p))
  maxit <- 1000
  eps <- 10^(-8)
  diff <- 10^6
  ctest <- talpha
  iter <- 1
  while ((diff>eps)& (iter<maxit)) {
    cold <- ctest
    ctest <- Tbsb(cold, p)/alpha
    diff <- abs(cold-ctest)
    eter <- iter+1
  }
  res <- ctest
  return(res)
}
##--------------------------------------------------------------------
##SESTCK
## Computes Tukey's biweight objectief function (scale) corresponding 
## with the mahalanobis distances x.  

##x must be entered as a matrix!!

sestck <- function(x, start, c, k, tol){
  if (start>0){
    s <- start
  }
  else { 
    y <- abs(x)
    if (is.matrix(x)==0) {
      x <- t(as.matrix(x))
    }
    n <- nrow(x)
    p <- ncol(x)
    if ((n==1)&(p>1)){
      x <- t(x)
      n <- nrow(x)
      p <- ncol(x)
    }
    x <- apply(x, 2, sort)
    if (floor(n/2)==n/2){
      s <- (x[n/2,]+x[(n+2)/2,])/2
    }
    else {
      s <- x[(n+1)/2,]
    }
    s <- s/.6745
  }
  crit <- 2*tol
  rhoold <- mean(rhobiweight(x/s, c))-k
  while (crit>=tol){
    delta <- rhoold/mean(psibiweight(x/s, c)*(x/(s^2)))
    isqu <- 1
    okay <- 0
    while((isqu<10)&(okay!=1)){
      rhonew <- mean(rhobiweight(x/(s+delta), c))-k
      if (abs(rhonew)<abs(rhoold)){
        s <- s+delta
        okay <- 1
      } 
      else {
        delta <- delta/2 
        isqu <- isqu+1
      }
    }
    if (isqu==10){
      crit <- 0
    }
    else {
      crit <- (abs(rhoold)-abs(rhonew))/max(abs(rhonew), tol)
    }
    rhoold <- rhonew
  }
  scale <- abs(s)
  return(scale)
}
##--------------------------------------------------------------------
## Computes biweight multivariate S-estimator for location/scatter with algorithm of 
## Ruppert
##
## Input:
##
## Required arguments
## x     :a data matrix. Rows of the matrix represent observations, 
##        columns represent variables. 
## nsamp :The number of random p-subsets considered.
## bdp   :Breakdown value of the S-estimator must be 0.15, 0.25 or 0.5.
##
##
## Output (field):
##
## mean        : vector of estimates for the center of the data.
## covariance  : matrix of estimates for the scatter of the data.
## distances   : vector of robust distances versus mean & covariance.
## scale       : distance scale estimate.

##SM
##p+1 must be greater than n!!!

mySm <- function(x, nsamp=2000, bdp=0.25){
  tol <- 10^(-5)
  seed <- 0
  s <- 10^(11)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  c <- Tbsc(bdp, p)
  k <- (c/6)*Tbsb(c, p)
  la <- 1
  for (i in 1:nsamp){
    ranset <- randomset(n, p+1, seed)[[1]]
    seed <- randomset(n, p+1, seed)[[2]]
    xj <- as.matrix(x[ranset,])
    mu <- apply(xj, 2, mean)
    cov <- var(xj)*(nrow(xj)-1)/nrow(xj)
    determ <- det(cov)
    if (determ>10^(-15)) {
      if (determ^(1/p)>10^(-5)) {
        cov <- determ^(-1/p)*cov
        if (i>ceiling(nsamp/5)) {
          if (i==ceiling(nsamp/2)){
            la <- 2
          }
          else {
            if (i==ceiling(nsamp*.8)) {
              la <- 4
            }
          }
          
          random <-  uniran(seed)[[2]]
          seed  <-  uniran(seed)[[1]] 
          random <- random^la
          mu <- random*mu+(1-random)*muopt
          cov <- random*cov+(1-random)*covopt
        }
        determ <- det(cov)
        cov <- determ^(-1/p)*cov
        md <- mahalanobis(x, mu, cov, inverted = FALSE, tol.inv =.Machine$double.eps)
        md <- md^(1/2)
        if (mean(rhobiweight(md/s, c))<k) {
          if (s<5*10^10) {
            s <- sestck(md, s, c, k, tol)
          }
          else {
            s <- sestck(md, 0, c, k, tol)
          }
          muopt <- mu
          covopt <- cov
          mdopt <- md
          psi <- psibiweight(md, s*c)
          u <- psi/md
          ubig <- matrix(t(u), nrow=length(u), ncol=p, byrow=FALSE)
          aux <- (ubig*x)/mean(u)
          mu <- apply(aux, 2, mean)
          xcenter <- t(t(x)-mu)
          cov <- t(ubig*xcenter)%*%xcenter
          cov <- det(cov)^(-1/p)*cov
          okay <- 0
          jj <- 1
          while ((jj<3)&(okay!=1)) {
            jj <- jj+1
            md <- mahalanobis(x, mu, cov, tol.inv =.Machine$double.eps)
            md <- md^(1/2)
            if (mean(rhobiweight(md/s, c))<k){
              muopt <- mu
              covopt  <- cov
              mdopt  <- md
              okay <- 1
              if (s<5*10^10) {
                s <- sestck(md, s, c, k, tol)}
              else {
                s <- sestck(md, 0, c, k, tol)}
            }
            else {
              mu <- (mu+muopt)/2
              cov <- (cov+covopt)/2
              cov <- determ^(-1/p)*cov
            }
          }
        }
      }
    }
  }
  covopt <- s^2*covopt
  mdopt <- mdopt/s
  res.mean <- muopt
  res.covariance <- covopt
  res.distances <- mdopt
  res.scale <- s
  return(list(center=res.mean, cov=res.covariance, res.distances, res.scale, c))
}

##end


## RIV estimator
betaIV <- function(L, V){
  b1 <- V[1,3]/V[2,3]
  b0 <- L[1]-b1*L[2]
  beta <- c(b0,b1)
  return(beta)
}
