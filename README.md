# Universal Difference-in-Differences

This Github repository contains UDID R package that implements the methodologies in [Park and Tchetgen Tchetgen (2026+)](https://arxiv.org/abs/2212.13641 "UDID").
This package is currently in beta.

## Installation

To install UDID package in R, run the commands below:

```{r}
library(devtools)
install_github("qkrcks0218/UDID")
```

## Methodology

Consider a panel study with two time periods $t \in \{0, 1\}$, where $t = 0$ is pre-treatment and $t = 1$ is post-treatment. Let $A \in \{0, 1\}$ be a binary treatment indicator, $Y_t$ the observed outcome at time $t$, $Y_t(a)$ the potential outcome at time $t$ under treatment $A = a$, and $X \in \mathbb{R}^d$ a vector of pre-treatment covariates. The target estimand is the **average treatment effect on the treated (ATT)**:

$$\tau = E \left[Y_1(1) - Y_1(0) \mid A = 1\right].$$

Under consistency ($Y_t = Y_t(A)$ a.s.) and no anticipation ($Y_0(0) = Y_0(1)$ a.s.), the first term $E[Y_1(1) \mid A=1]$ is directly identified from the data. The key challenge is identifying the counterfactual mean $E[Y_1(0) \mid A=1]$.

### Odds Ratio Equi-Confounding (OREC) Assumption

Let $f_t(y \mid a, x)$ denote the conditional density of $Y_t(0)$ given $A = a$ and $X = x$, and fix a reference value $y_R$ satisfying $(y_R, x) \in \mathcal{S}$ (the common support). For each $t \in \{0, 1\}$, define the **generalized odds ratio function**:

$$\alpha_t(y, x) = \frac{f_t(y \mid A=1,   X=x)}{f_t(y_R \mid A=1,   X=x)} \cdot \frac{f_t(y_R \mid A=0,   X=x)}{f_t(y \mid A=0,   X=x)}.$$

By construction, $\alpha_t(y_R, x) = 1$ for all $x$, and $\alpha_t(y, x) = 1$ for all $(y, x)$ whenever $Y_t(0) \perp\perp A \mid X$.

The key assumption of the UDID framework is **odds ratio equi-confounding (OREC)**, which states that

$$\alpha_0(y, x) = \alpha_1(y, x) \quad \text{for all } (y, x) \in \mathcal{S}.$$

In other words, the generalized odds ratio function is stable across time periods. Intuitively, OREC states that the degree of unmeasured confounding between treatment $A$ and the treatment-free potential outcome $Y_t(0)$—measured on the odds ratio scale—does not change from the pre-treatment period to the post-treatment period. This is a distributional analogue of the classical parallel trends assumption, but is scale-invariant and compatible with continuous, binary, or mixed-type outcomes.

### ATT Identification

Under OREC, the counterfactual mean $\mu(x) = E[Y_1(0) \mid A=1, X=x]$ is identified as:

$$\mu(x) = \frac{E \left[Y_1  \alpha_0(Y_1, X) \mid A=0,   X=x\right]}{E \left[\alpha_0(Y_1, X) \mid A=0,   X=x\right]},$$

where $\alpha_1$ is replaced by the identified $\alpha_0$ under OREC. The ATT is then identified via $\tau = E[\mu(X) \mid A=1]$.

### Efficient Influence Function

The **efficient influence function (EIF)** for $\tau$ in the nonparametric model under OREC (Theorem 5.1 in Park and Tchetgen Tchetgen, 2026+) is:

$$\mathrm{IF}(O) = \frac{A Y_1 - \phi_0(O)  - A \tau }{\Pr(A=1)},$$

where

$$\phi_0(O) = \beta_1(X)  \alpha_1(Y_1, X)  (1-A) \left(Y_1 - \mu(X)\right) + A  \mu(X) + (2A-1)  R(Y_0, A, X) \left(Y_0 - \mu(X)\right). $$

Here, $\beta_1(x) = \Pr(A=1 \mid Y_1(0)=y_R,X=x)/ \Pr(A=0 \mid Y_1(0)=y_R,X=x)$ is the baseline odds at time $t=1$ and $R(y, a, x) = f_1(y,1,x)/f_0(y,a,x)$ is the density ratio. 
These functions are identified under OREC; see Park and Tchetgen Tchetgen (2026+) for details.

The ATT is estimated as the sample mean of the estimated uncentered EIF after cross-fitting (Chernozhukov et al., 2017):

$$\hat{\tau} = \frac{\mathbb{P}_n \left[A Y_1 - \hat{\phi}_0(O)\right]}{\mathbb{P}_n(A)},$$

where $\mathbb{P}_n$ denotes the empirical average computed across the cross-fitted folds.

### Nuisance Function Estimation

All nuisance functions are estimated via 2-fold cross-fitting to avoid overfitting bias and ensure valid inference.
In each fold, nuisance functions are estimated on one half of the data and evaluated on the other.

**Continuous outcomes.** Three nuisance components are estimated nonparametrically:

1. **Conditional density** $f_0(Y_0 \mid A=0, X)$: estimated via the nonparametric kernel conditional density estimator of Hall et al. (2004), implemented in the `np` R package (Hayfield and Racine, 2008). The bandwidth is selected either by the rule-of-thumb (Silverman, 1986) or by cross-validation.

2. **Density ratios** $f_0(Y_0 \mid A=1, X) / f_0(Y_0 \mid A=0, X)$ and $f_1(Y_1 \mid A=0, X) / f_0(Y_0 \mid A=0, X)$: estimated directly via the Kullback–Leibler Importance Estimation Procedure (KLIEP; Sugiyama et al., 2007; Nguyen et al., 2007), implemented in the `densratio` R package (Makiyama, 2019). KLIEP models the density ratio as a positive kernel function in a Reproducing Kernel Hilbert Space (RKHS) with a Gaussian kernel, and selects the bandwidth by the median heuristic (Garreau et al., 2017) or cross-validation.

3. **Propensity score** $\Pr(A=1 \mid X)$: estimated via an ensemble of machine learning algorithms using the Super Learner (van der Laan et al., 2007), implemented in the `SuperLearner` R package (Polley et al., 2025). Available base learners include GLM, lasso/ridge, earth, GAM, XGBoost, polynomial splines, random forests, gradient boosting machines, and one-layer neural networks.

**Binary outcomes.** All conditional distributions, $\Pr(A=1 \mid X)$, $\Pr(Y_0=1 \mid A, X)$, and $\Pr(Y_1=1 \mid A=0, X)$, are estimated via Super Learner.

Given the estimated nuisance functions, $\hat{\alpha}_1$, $\hat{\beta}_1$, $\hat{\mu}$, and $\hat{R}$ are constructed from the estimated densities and density ratios, and plugged into the EIF to obtain $\hat{\tau}$.

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

Park, C. & Tchetgen Tchetgen, E. (2026+) **A Universal Nonparametric Framework for Difference-in-Differences Analyses**. _arXiv_ [[link](https://arxiv.org/abs/2212.13641 "UDID")]

Victor Chernozhukov, Denis Chetverikov, Mert Demirer, Esther Duflo, Christian Hansen, Whitney Newey, James Robins (2017) Double/debiased machine learning for treatment and structural parameters. _The Econometric Journal_, 21: C1-C68

Garreau, D., Jitkrittum, W., & Kanagawa, M. (2017). Large sample analysis of the median heuristic. _arXiv_

Hall, P., Racine, J. S., & Li, Q. (2004). Cross-validation and the estimation of conditional probability densities. _Journal of the American Statistical Association_, 99, 1015--1026.

Hayfield, T., & Racine, J. S. (2008). Nonparametric econometrics: The np package. _Journal of Statistical Software_, 27, 1--32.

Makiyama, K. (2019). densratio: Density Ratio Estimation. R package version 0.2.1.

Nguyen, X., Wainwright, M. J., & Jordan, M. (2007). Estimating divergence functionals and the likelihood ratio by penalized convex risk minimization. _Advances in Neural Information Processing Systems_, 20.

Polley, E., LeDell, E., Kennedy, C., Lendle, S., & van der Laan, M. (2025). SuperLearner: Super Learner Prediction. R package version 2.0-40.

Silverman, B. W. (1986). _Density Estimation for Statistics and Data Analysis_. London: Chapman and Hall.

Sugiyama, M., Nakajima, S., Kashima, H., Buenau, P., & Kawanabe, M. (2007). Direct importance estimation with model selection and its application to covariate shift adaptation. _Advances in Neural Information Processing Systems_, 20.

van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. _Statistical Applications in Genetics and Molecular Biology_, 6(1).
