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

Data <- read.csv("ZikaData.csv")
 
mY.save <- mean(c(Data$Yn1,Data$Y0,Data$Y1))
sY.save <- sd(c(Data$Yn1,Data$Y0,Data$Y1))
Data$Yn1 <- (Data$Yn1 - mY.save)
Data$Y0 <- (Data$Y0 - mY.save)
Data$Y1 <- (Data$Y1 - mY.save)


LongData <- matrix(0,nrow(Data)*3,8)
LongData[,1] <- rep(1:nrow(Data),each=3)
LongData[,2] <- rep(1:3,nrow(Data))
LongData[,3] <- rep(Data$A,each=3)
LongData[,4] <- rep(Data$X1,each=3)
LongData[,5] <- rep(Data$X2,each=3)
LongData[,6] <- rep(Data$X3,each=3)
LongData[,7] <- rep(Data$X4,each=3)
LongData[,8] <- as.numeric(t(cbind(Data$Yn1,Data$Y0,Data$Y1)))

LongData <- data.frame(LongData)
colnames(LongData) <- c("ID",
                        "Time",
                        "Trt",
                        "LogPop",
                        "LogPDen",
                        "PropF",
                        "GDP",
                        "Y")
LongData$Trt <- LongData$Trt*3

################################################################################
# DID
################################################################################

T.DID.0 <- Sys.time()
DID_Nonpara_X <- att_gt(yname="Y", tname="Time", idname="ID", gname="Trt", 
                        xformla=~0+LogPop+LogPDen+PropF+GDP, data=LongData,
                        est_method=function(y1, y0, D, covariates, i.weights, inffunc){
                          my_did_MySL(y1, y0, D, covariates, i.weights, inffunc,
                                      type = "continuous",
                                      seed = BATCH,
                                      SL.list = SL.list)
                        })
T.DID.1 <- Sys.time()

################################################################################
# UDID
################################################################################

###########################
# UDID_Placebo
###########################

T00 <- Sys.time()
UDID_NP_X_Placebo <- UDID_Nonparametric(Y0=Data$Yn1,
                                        Y1=Data$Y0,
                                        A=Data$A,
                                        X=Data[,2:5],
                                        type="continuous",
                                        log_Gamma_seq=0,
                                        seed=BATCH,
                                        SL.list = SL.list,
                                        hyperparameter="fast")
T10 <- T01 <- Sys.time()

UDID_NP_X_Main <- UDID_Nonparametric(Y0=Data$Y0,
                                     Y1=Data$Y1,
                                     A=Data$A,
                                     X=Data[,2:5],
                                     type="continuous",
                                     log_Gamma_seq=0,
                                     seed=BATCH,
                                     SL.list = SL.list,
                                     hyperparameter="fast",
                                     density=TRUE)
T11 <- Sys.time()

UDID_NP_X_Main_Sens <- UDID_Nonparametric(Y0=Data$Y0,
                                          Y1=Data$Y1,
                                          A=Data$A,
                                          X=Data[,2:5],
                                          type="continuous",
                                          log_Gamma_seq=seq(0,1,by=0.01),
                                          seed=BATCH,
                                          SL.list = SL.list,
                                          hyperparameter="fast")

RESULT <- c( BATCH,
             DID_Nonpara_X$att,
             as.numeric(UDID_NP_X_Placebo$Effect[1]),
             as.numeric(UDID_NP_X_Main$Effect[1]),
             
             DID_Nonpara_X$se,
             as.numeric(UDID_NP_X_Placebo$Effect[2]),
             as.numeric(UDID_NP_X_Main$Effect[2]),
             as.numeric(UDID_NP_X_Placebo$Effect[3]),
             as.numeric(UDID_NP_X_Main$Effect[3]),
             
             as.numeric(difftime(T.DID.1, T.DID.0, units = "secs")),
             as.numeric(difftime(T01, T00, units = "secs")),
             as.numeric(difftime(T11, T10, units = "secs")))

RESULT <- matrix(RESULT,1,length(RESULT))
colnames(RESULT) <- c("BATCH",
                      sprintf("Est_DID_NP_X_%0.2d",1:2),
                      sprintf("Est_UDID_NP_X_%0.2d",1:2),
                      
                      sprintf("SE_DID_NP_X_%0.2d",1:2),
                      sprintf("SE_UDID_NP_X_%0.2d",1:2),
                      sprintf("BSE_UDID_NP_X_%0.2d",1:2),
                      
                      "DID_Time",
                      sprintf("UDID_Time_%0.2d",1:2))

write.csv(RESULT,
          sprintf("Result_Zika_B%0.5d.csv",BATCH),
          row.names=FALSE)


write.csv(UDID_NP_X_Main_Sens$Sensitivity,
          sprintf("Result_Zika_Sens_B%0.5d.csv",BATCH),
          row.names=FALSE)


FData <- cbind(Data,matrix(0,nrow(Data),501*4))
FData[, ncol(Data)+0*501+1:501 ] <- UDID_NP_X_Main$Density$fY0A0
FData[, ncol(Data)+1*501+1:501 ] <- UDID_NP_X_Main$Density$fY0A1
FData[, ncol(Data)+2*501+1:501 ] <- UDID_NP_X_Main$Density$fY1A0
FData[, ncol(Data)+3*501+1:501 ] <- UDID_NP_X_Main$Density$fY1A1 

write.csv(FData,
          sprintf("Result_Zika_Density_B%0.5d.csv",BATCH),
          row.names=FALSE)














