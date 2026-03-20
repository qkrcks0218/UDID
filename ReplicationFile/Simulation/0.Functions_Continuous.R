library(MASS)
library(Matrix)
library(pracma)

sigma0.preliminary <- 2
sigma1.preliminary <- 1

expit <- function(x){
  exp(x)/(1+exp(x))
}

MM <- function(X){
  if( is.null(dim(X)) ){
    X <- matrix(X,length(X),1)
  } else {
    X
  }
  return(X)
}

ProbA <- function(X){
  X <- MM(X)
  expit( 0.1*X[,1]+0.1*X[,2] )
}

# Testing This:
alpha.preliminary <- function(X){
  X <- MM(X)
  0.05+0.02*X[,1]+0.02*X[,2]
} 

muY0 <- function(A,X){
  X <- MM(X)
  3 + alpha.preliminary(X)*sigma0.preliminary^2*A + 
    + 0.1*X[,1] + 0.1*X[,2]
}

muY10 <- function(A,X){
  X <- MM(X)
  3.5 +
    alpha.preliminary(X)*sigma1.preliminary^2*A +
    + 0.1*X[,1] + 0.1*X[,2]
}

muY11 <- function(A,X){
  X <- MM(X)
  4 +
    alpha.preliminary(X)*sigma1.preliminary^2*A +
    + 0.1*X[,1] + 0.1*X[,2]
}


odds_A_given_Y10 <- function(y, X) {
  densratio <- dnorm(y, mean = muY10(1, X), sd = sigma1.preliminary) /
    dnorm(y, mean = muY10(0, X), sd = sigma1.preliminary)
  pA <- ProbA(X)
  prior_odds <- pA / (1 - pA)
  densratio*prior_odds
}

density_ratio_Y10_over_Y0 <- function(y, a, X) {
  dnorm(y, mean = muY10(1, X), sd = sigma1.preliminary) /
    dnorm(y, mean = muY0(a, X), sd = sigma0.preliminary)
}


phi0 <- function(Y1,Y0,A,X){
  (1-A)*odds_A_given_Y10(Y1,X)*(Y1-muY10(1,X)) + A*muY10(1,X) + 
    (2*A-1)*density_ratio_Y10_over_Y0(Y0,A,X)*(Y0-muY10(1,X))
}

Oracle.UDID <- function(Y1,Y0,A,X){
  tau <- mean(A*Y1-phi0(Y1,Y0,A,X))/mean(A)
  eif <- (A*Y1-phi0(Y1,Y0,A,X)-A*tau)/mean(A)
  return(list(ATT=tau,SE=sd(eif)/sqrt(length(eif))))
}













