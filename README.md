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

SL.hpara <- list()
SL.hpara$SLL <- c(1,2)
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 9: gbm
# 10: 1-layer MLP
SL.hpara$MLPL <- c(2,4)
SL.hpara$MTRY <- c(1,2)
SL.hpara$NMN <- 20
SL.hpara$MLPdecay <- 10^c(-1,-3)

Data <- read.csv("ZikaData.csv")

mY.save <- mean(c(Data$Yn1,Data$Y0,Data$Y1))
sY.save <- sd(c(Data$Yn1,Data$Y0,Data$Y1))
Data$Yn1 <- (Data$Yn1 - mY.save)
Data$Y0 <- (Data$Y0 - mY.save)
Data$Y1 <- (Data$Y1 - mY.save)

X.pos <- which(substring(colnames(Data),1,1)=="X")

UDID.Result <- UDID_Nonparametric(Data$Y0,Data$Y1,Data$A,Data[,X.pos],
                                  type="continuous")

UDID.Result

UDID.SensPara <- UDID_Nonparametric_SensPara(Data$Yn1,Data$Y0,Data$A,Data[,X.pos],
                                             type="continuous",
                                             quantile.range=c(0.01,0.025,0.975,0.99))
log(UDID.SensPara)

UDID.Sensitivity <- UDID_Nonparametric(Data$Y0,Data$Y1,Data$A,Data[,X.pos],
                                       type="continuous",
                                       log_Gamma_seq = seq(0,1,by=0.1))

UDID.Sens.Bound <- UDID_Sensitivity_Bounds(UDID.Sensitivity)

UDID_Sensitivity_Plot(UDID.Sens.Bound)


```

![Alt text](./images/Sensitivity.png?raw=true "Sensitivity.png")


## References

Chan Park & Eric Tchetgen Tchetgen (2026+) **A Universal Nonparametric Framework for Difference-in-Differences Analyses**, _arXiv_ [[link](https://arxiv.org/abs/2212.13641 "UDID")]