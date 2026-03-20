args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])  #folder number
  set.seed(BATCH)
} else {
  stop()
}

###########################
# Source Files
###########################

library(UDID)
source("0.NPDiD.R")

###########################
# SL Parameters
###########################

SL.list <- 1:9
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 8: gbm
# 9: 1-layer MLP

###########################
# DGP
###########################

PARA.GRID <- expand.grid(1:2,
                         c(1,2,3,4,8)*500,
                         1:10000)
TYPE <- PARA.GRID[BATCH,1]
N <- PARA.GRID[BATCH,2]
SEED <- PARA.GRID[BATCH,3]

set.seed(SEED)
d <- 2

if(TYPE==1){
  
  source("0.Functions_Continuous.R")
  Outcome.Type <- "continuous"
  True.ATT <- 0.5
  
  X <- matrix(rnorm(N*d),N,d)
  A <- rbinom(N,1,ProbA(X))
  Y0 <- rnorm( N, mean=muY0(A,X), sd=2 )
  Y10 <- rnorm( N, mean=muY10(A,X), sd=1 )
  Y11 <- rnorm( N, mean=muY11(A,X) )
  Y1 <- Y10*(1-A)+Y11*A
  
} else {
  
  source("0.Functions_Discrete.R")
  Outcome.Type <- "binary"
  True.ATT <- 0
  
  X <- matrix(rnorm(N*d),N,d)
  A <- rbinom(N,1,ProbA(X))
  Y0 <- rbinom(N,1,probY0givenA(A,X))
  Y10 <- rbinom(N,1,probY10givenA(A,X))
  Y11 <- rbinom(N,1,probY11givenA(A,X))
  Y1 <- Y10*(1-A)+Y11*A
  
}

Data <- as.matrix(cbind(Y0,Y1,Y10,Y11,A,X))
Data <- data.frame(Data)
colnames(Data) <- c("Y0","Y1","Y10","Y11","A",sprintf("X%0.1d",1:d))

X.pos <- which(substring(colnames(Data),1,1)=="X")

LongData <- matrix(0,2*N,4+d)
LongData[,1] <- c(Y0,Y1)
LongData[,2] <- c(rep(0,N),rep(1,N))
LongData[,3] <- c(A,A)
LongData[,3+1:d] <- rbind(X,X)
LongData[,4+d] <- c(1:N,1:N)

LongData <- data.frame(LongData)
colnames(LongData) <- c("Y","Time","Trt",sprintf("X%0.1d",1:d),"ID")

###########################
# Run UDID
###########################

DID.Result <- UDID.Result <- list()

T0 <- Sys.time()

for(ss in 1:3){
  UDID.Result[[ss]] <- UDID_Nonparametric(Y0=Data$Y0,
                                          Y1=Data$Y1,
                                          A=Data$A,
                                          X=Data[,X.pos],
                                          type=Outcome.Type,
                                          log_Gamma_seq=0,
                                          seed=10*SEED+ss,
                                          SL.list = SL.list,
                                          hyperparameter="fast")
}

T1 <- Sys.time()


T2 <- Sys.time()

for(ss in 1:3){
  
  DID.Result[[ss]] <- att_gt(yname="Y", tname="Time", idname="ID", gname="Trt", 
                             xformla=~0+X1+X2, data=LongData,
                             est_method=function(y1, y0, D, covariates, i.weights, inffunc){
                               my_did_MySL(y1, y0, D, covariates, i.weights, inffunc,
                                           type = Outcome.Type,
                                           seed = 10*SEED+ss,
                                           SL.list = SL.list)
                             })
  
}

T3 <- Sys.time()

UDID.Result.Mat <- rbind( UDID.Result[[1]]$Effect,
                          UDID.Result[[2]]$Effect,
                          UDID.Result[[3]]$Effect)

DID.Result.Mat <- rbind( c( DID.Result[[1]]$att, DID.Result[[1]]$se ),
                         c( DID.Result[[2]]$att, DID.Result[[2]]$se ),
                         c( DID.Result[[3]]$att, DID.Result[[3]]$se ) )

UDID.ATT <- median(UDID.Result.Mat[,1], na.rm=TRUE)
UDID.ASE <- median( sqrt( (UDID.Result.Mat[,1] - UDID.ATT)^2 + UDID.Result.Mat[,2]^2 ), na.rm=TRUE)
UDID.BSE <- median( sqrt( (UDID.Result.Mat[,1] - UDID.ATT)^2 + UDID.Result.Mat[,3]^2 ), na.rm=TRUE)

DID.ATT <- median(DID.Result.Mat[,1], na.rm=TRUE)
DID.ASE <- median( sqrt( (DID.Result.Mat[,1] - DID.ATT)^2 + DID.Result.Mat[,2]^2 ), na.rm=TRUE)

UDID.Oracle.Result <- Oracle.UDID(Data$Y1,Data$Y0,Data$A,Data[,X.pos])

Result <- c(
  TYPE,
  N,
  SEED,
  
  UDID.ATT - True.ATT ,
  UDID.ASE, 
  UDID.BSE, 
  as.numeric(abs(UDID.ATT - True.ATT) <= qnorm(0.975)*UDID.ASE) ,
  as.numeric(abs(UDID.ATT - True.ATT) <= qnorm(0.975)*UDID.BSE) ,
  
  DID.ATT - True.ATT ,
  DID.ASE, 
  as.numeric(abs(DID.ATT - True.ATT) <= qnorm(0.975)*DID.ASE) ,
  
  as.numeric(UDID.Oracle.Result$ATT) ,
  as.numeric(UDID.Oracle.Result$SE) ,
  as.numeric(abs(UDID.Oracle.Result$ATT - True.ATT) <= qnorm(0.975)*UDID.Oracle.Result$SE),
  
  as.numeric(difftime(T1, T0, units = "secs")/3),
  as.numeric(difftime(T3, T2, units = "secs")/3)) 

Result <- matrix(Result,1,length(Result))

colnames(Result) <- c("Type","N","Seed",
                      "UDID_ATT","UDID_ASE","UDID_BSE","UDID_Cover_ASE","UDID_Cover_BSE",
                      "DID_ATT","DID_ASE","DID_Cover_ASE",
                      "UDID_Oracle_ATT","UDID_Oracle_ASE","UDID_Oracle_Cover_ASE",
                      "Time_UDID", "Time_DID" )

write.csv(Result,
          sprintf("Result_Simulation_Type%0.1d_N%0.5d_B%0.5d.csv",TYPE,N,SEED),
          row.names=FALSE)


