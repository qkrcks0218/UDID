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


###########################
# Sensitivity Range
###########################

# Zika_UDID_NP_X_Gamma_Range <-
#   matrix(0,100,6)
# 
# for(ss in 1:100){
#   Zika_UDID_NP_X_Gamma_Range[ss,] <-
#     UDID_Nonparametric_SensPara(Yn1=Data$Yn1,
#                                 Y0=Data$Y0,
#                                 A=Data$A,
#                                 X=Data[,2:5],
#                                 type="continuous",
#                                 hyperparameter="fast",
#                                 quantile.range=c(0,0.01,0.025,0.975,0.99,1),
#                                 SL.list = SL.list,
#                                 seed=ss)
#   
# }
# 
# log(apply(Zika_UDID_NP_X_Gamma_Range,2,median))
# -0.3975529 -0.1446995 -0.1147944  0.1627893  0.2776814  1.0101133





result_sens <- read.csv("RESULT_Zika_Sens_Merge.csv")
Aggregate.ATT <- 
  aggregate(cbind(ATT_LB,ATT_UB)~log_Gamma,
            result_sens,
            "median")

Aggregate.ATT$SE_LB <- Aggregate.ATT$SE_UB <- NULL

LG <- Aggregate.ATT[,1]

for(gi in 1:101){
  Temp.Data <- result_sens[result_sens$log_Gamma==LG[gi],]
  Aggregate.ATT$SE_LB[gi] <- 
    median((Temp.Data$ATT_LB-Aggregate.ATT[gi,2])^2 + Temp.Data$SE_LB^2)^(1/2)
  Aggregate.ATT$SE_UB[gi] <- 
    median((Temp.Data$ATT_UB-Aggregate.ATT[gi,3])^2 + Temp.Data$SE_UB^2)^(1/2)
}

result <- UDID_Sensitivity_Bounds(Aggregate.ATT[1:51,],0.95)$bounds
ylim <- range(c(result$CI_LB, result$CI_UB), na.rm = TRUE)
ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05

png("Sensitivity.png",height=5,width=10,unit="in",res=500)
layout(matrix(1:2,1,2))
par(mar=c(4,4,2,0.5))
UDID_Sensitivity_Plot(Aggregate.ATT[1:51,], main="",legend.pos="topleft")
abline(v=0.1627893, lty=2, col=rgb(0,0,0,0.25))
text(x=0.1627893, y=ylim[1], 
     labels="0.163",
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)

abline(v=0.2776814, lty=2, col=rgb(0,0,0,0.25))
text(x=0.2776814, y=ylim[1], 
     labels="0.278",
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)

Aggregate.ATT.Gamma <- Aggregate.ATT

Aggregate.ATT.Gamma$log_Gamma <- exp(Aggregate.ATT.Gamma$log_Gamma)
UDID_Sensitivity_Plot(Aggregate.ATT.Gamma[1:51,],xlab=expression(Gamma),main="",ylab="",legend.pos="topleft")
abline(v=exp(0.1627893), lty=2, col=rgb(0,0,0,0.25))
text(x=exp(0.1627893), y=ylim[1], 
     labels="1.177",
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)

abline(v=exp(0.2776814), lty=2, col=rgb(0,0,0,0.25))
text(x=exp(0.2776814), y=ylim[1], 
     labels="1.320",
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)

dev.off()


#########################################################


Data_or <- setNames(
  data.frame(cbind(Data$Y0 - Data$Yn1, Data$A, 
                   Data$X1, Data$X2, Data$X3, Data$X4)), 
  c("dY", "A", "X1", "X2", "X3", "X4" ))

Gamma.095.Vec <- Gamma.098.Vec <- rep(0,100)
for(tt in 1:100){
  set.seed(tt)
  SL_or_1 <- UDID::MySL(Data_or[Data_or$A==1,], locY = 1, locX = 2 + 1:4,
                        Ydist = gaussian(),
                        SL.list = SL.list)
  SL_or_0 <- UDID::MySL(Data_or[Data_or$A==0,], locY = 1, locX = 2 + 1:4,
                        Ydist = gaussian(),
                        SL.list = SL.list)
  
  mu_1 <- predict(SL_or_1,
                  newdata = Data_or[,-c(1,2)],
                  onlySL = TRUE)$pred
  
  mu_0 <- predict(SL_or_0,
                  newdata = Data_or[,-c(1,2)],
                  onlySL = TRUE)$pred
  
  Gamma.095.Vec[tt] <- max(abs(quantile(mu_1 - mu_0,c(0.025,0.975))))
  Gamma.098.Vec[tt] <- max(abs(quantile(mu_1 - mu_0,c(0.01,0.99))))
  print(tt)
}

Gamma.095 <- median(Gamma.095.Vec)
# 1.222484

Gamma.098 <- median(Gamma.098.Vec)
# 1.449523

result <- read.csv("RESULT_Zika_Merge.csv")

ATT.DID.1 <- median(result$Est_DID_NP_X_02)
SE.DID.1 <- sqrt(median(result$SE_DID_NP_X_02^2 + (result$Est_DID_NP_X_02-ATT.DID.1)^2))

ATT.DID <- c(ATT.DID.0, ATT.DID.1)
SE.DID <- c(SE.DID.0, SE.DID.1)

Gamma <- seq(0,1.5,by=0.01)
ATT.Sens.DID <- data.frame(log_Gamma=seq(0,1.5,by=0.01),
                           ATT_LB = ATT.DID.1-Gamma,
                           ATT_UB = ATT.DID.1+Gamma,
                           SE_LB = SE.DID.1,
                           SE_UB = SE.DID.1)


ylim <- range(c(ATT.DID.1-1.5-SE.DID.1*qnorm(0.975),
                ATT.DID.1+1.5-SE.DID.1*qnorm(0.975)), na.rm = TRUE)
ylim <- ylim + c(-1, 1) * diff(ylim) * 0.05

png("Sensitivity_DID.png",height=5,width=5,unit="in",res=500)
par(mar=c(4,4,2,0.5))

UDID_Sensitivity_Plot(ATT.Sens.DID,
                      xlab=expression(Gamma),
                      main="",ylab="",legend.pos="topleft")
abline(v=Gamma.095, lty=2, col=rgb(0,0,0,0.25))
text(x=Gamma.095, y=ylim[1], 
     labels=round(Gamma.095,3),
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)

abline(v=Gamma.098, lty=2, col=rgb(0,0,0,0.25))
text(x=Gamma.098, y=ylim[1], 
     labels=sprintf("%0.3f",Gamma.098),
     adj=c(0.5,-0.4),
     col=rgb(0,0,0,0.5),
     cex=0.8,pos=2)
dev.off()


























