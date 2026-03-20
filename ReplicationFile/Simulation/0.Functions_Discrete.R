library(MASS)
library(Matrix)
library(pracma)

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
  expit(0.1*X[,1]+0.1*X[,2])
}


OR.Original <- function(X){
  1 - 0.2*(X[,1]+X[,2])
}

probY0givenA <- function(A,X){
  X <- MM(X)
  expit( -0.5 + OR.Original(X)*A + 0.1*(X[,1]+X[,2]) )
}


probY10givenA <- function(A,X){
  X <- MM(X)
  expit( +0.5 + OR.Original(X)*A + 0.1*(X[,1]+X[,2]) )
}

probY11givenA <- function(A,X){
  X <- MM(X)
  probY10givenA(A,X)
}

alpha <- function(X){ 
  
  X <- MM(X)
  pA1 <- ProbA(X)
  pY1A0 <- probY0givenA(0,X)
  pY1A1 <- probY0givenA(1,X)
  
  (1-pA1)*(1-pY1A0)*(  pA1)*(  pY1A1)/((1-pA1)*(  pY1A0)*(  pA1)*(1-pY1A1))
}

mu10 <- function(X){
  X <- MM(X)
  probY11givenA(1,X)
}


odds_A_given_Y10 <- function(y, X) {
  X <- MM(X)
  
  # f(Y10=y|A=1,X) / f(Y10=y|A=0,X)
  p1 <- probY10givenA(1, X)
  p0 <- probY10givenA(0, X)
  lik_ratio <- (p1^y * (1 - p1)^(1 - y)) / (p0^y * (1 - p0)^(1 - y))
  
  # prior odds
  pA <- ProbA(X)
  prior_odds <- pA / (1 - pA)
  
  lik_ratio * prior_odds
}

density_ratio_Y10_A1_over_Y0 <- function(y, a, X) {
  X <- MM(X)
  
  # f(Y11=y|A=1,X) * P(A=1|X)
  p_Y10 <- probY10givenA(1, X)
  num <- (p_Y10^y * (1 - p_Y10)^(1 - y)) * ProbA(X)
  
  # f(Y0=y|A=a,X) * P(A=a|X)
  p_Y0 <- probY0givenA(a, X)
  pA <- ProbA(X)
  pAa <- pA * a + (1 - pA) * (1 - a)
  den <- (p_Y0^y * (1 - p_Y0)^(1 - y)) * pAa
  
  num / den
}


phi0 <- function(Y1, Y0, A, X) {
  (1 - A) * odds_A_given_Y10(Y1, X) * (Y1 - mu10(X)) + A * mu10(X) +
    (2*A - 1) * density_ratio_Y10_A1_over_Y0(Y0, A, X) * (Y0 - mu10(X))
}

Oracle.UDID <- function(Y1,Y0,A,X){
  tau <- mean(A*Y1-phi0(Y1,Y0,A,X))/mean(A)
  eif <- (A*Y1-phi0(Y1,Y0,A,X)-A*tau)/mean(A)
  return(list(ATT=tau,SE=sd(eif)/sqrt(length(eif))))
}
