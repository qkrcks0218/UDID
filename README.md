# Universal Difference-in-Differences

This Github repository contains UDID R package that implements the methodologies in [Park and Tchetgen Tchetgen (2026+)](https://arxiv.org/abs/2212.13641 "UDID"). 
This package is currently in beta.

## Installation

To install UDID package in R, run the commands below:

```{r}
library(devtools)
install_github("qkrcks0218/UDID")
```

## Example usage

Here are some examples:

```{r}
library(UDID)

###########################################
# SuperLearner Hyperparameters
###########################################

SL.list <- c(1,2)
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

###########################################
# Data
###########################################

Data <- read.csv("ZikaData.csv")

mY.save <- mean(c(Data$Yn1,Data$Y0,Data$Y1))
sY.save <- sd(c(Data$Yn1,Data$Y0,Data$Y1))
Data$Yn1 <- (Data$Yn1 - mY.save)
Data$Y0 <- (Data$Y0 - mY.save)
Data$Y1 <- (Data$Y1 - mY.save)

X.pos <- which(substring(colnames(Data),1,1)=="X")

###########################################
# UDID
###########################################

UDID.Result <- UDID_Nonparametric(Data$Y0,Data$Y1,Data$A,Data[,X.pos],
                                  type="continuous",
                                  SL.list = SL.list,
                                  density.report=TRUE,
                                  hyperparameter.report=TRUE)

UDID.Result$Effect

#          ATT        SE   Boot.SE
# 1 -0.8932693 0.1237851 0.1289639

UDID.Result$EIF[c(1,1000),]

#      Uncentered EIF Centered EIF
# 1       -0.06212923  -0.06212923
# 1000    -2.80384687  -0.63838159

UDID.Result$SL_Library 

#   SuperLearner Library Fold 1: Pr(A|X) Fold 2: Pr(A|X)
# 1         SL.glm_1_All               1       0.0000000
# 2      SL.glmnet_1_All               0       0.7982535
# 3      SL.glmnet_2_All               0       0.0000000
# 4      SL.glmnet_3_All               0       0.2017465
# 5      SL.glmnet_4_All               0       0.0000000
# 6      SL.glmnet_5_All               0       0.0000000
# 7      SL.glmnet_6_All               0       0.0000000

UDID.Result$selected_bws[,1:2]

#    Fold 1: f(Y0|A=0,X) Fold 2: f(Y0|A=0,X)
# Y            1.9451773           1.8820872
# X1           0.5688129           0.5534825
# X2           0.5370118           0.4965944
# X3           0.4508771           0.3800171
# X4           0.2377253           0.2406518

plot(UDID.Result$Density$Y.Grid,
     UDID.Result$Density$fY0A0[1,],
     type='l')
```

![Alt text](./images/Density.png?raw=true "Density.png")

```{r}
UDID.Result.Binary <- 
  UDID_Nonparametric(as.numeric(Data$Y0>0),as.numeric(Data$Y1>0),Data$A,Data[,X.pos],
                     type="binary",
                     SL.list = SL.list,
                     density.report=TRUE,
                     hyperparameter.report=TRUE)

###########################################
# UDID Sensitivity Parameter Range
###########################################

UDID.SensPara <- UDID_Nonparametric_SensPara(Data$Yn1,Data$Y0,Data$A,Data[,X.pos],
                                             type="continuous",
                                             quantile.range=c(0.01,0.025,0.975,0.99))
log(UDID.SensPara)

#           1%       2.5%      97.5%        99% 
#   -0.4193990 -0.4040136  0.5059941  0.6195689 

###########################################
# UDID Sensitivity Analysis
###########################################

UDID.Sensitivity <- UDID_Nonparametric(Data$Y0,Data$Y1,Data$A,Data[,X.pos],
                                       type="continuous",
                                       log_Gamma_seq = seq(0,1,by=0.1))

UDID.Sens.Bound <- UDID_Sensitivity_Bounds(UDID.Sensitivity)

UDID_Sensitivity_Plot(UDID.Sens.Bound)

```

![Alt text](./images/Sensitivity.png?raw=true "Sensitivity.png")


## References

Chan Park & Eric Tchetgen Tchetgen (2026+) **A Universal Nonparametric Framework for Difference-in-Differences Analyses**, _arXiv_ [[link](https://arxiv.org/abs/2212.13641 "UDID")]