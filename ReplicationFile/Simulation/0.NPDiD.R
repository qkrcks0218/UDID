library(did)

my_did_MySL <- function(y1, y0, D, covariates, i.weights, inffunc,
                        type = "continuous",
                        seed = 42,
                        SL.list = c(1),
                        ...) {
  #-----------------------------------------------------------------
  # Doubly-robust DID (panel) with 2-fold cross-fitting via MySL
  # Mirrors the sample-split structure of UDID_Nonparametric
  #-----------------------------------------------------------------
  
  if(type=="continuous"){
    outcome.type <- gaussian()
  } else if(type=="binary") {
    outcome.type <- binomial()
  }
  
  n      <- length(D)
  deltaY <- y1 - y0
  w      <- if (missing(i.weights) || is.null(i.weights)) rep(1, n) else i.weights
  d      <- ncol(covariates)
  xnms   <- paste0("X", sprintf("%05d", seq_len(d)))
  
  ## ---- Stratified 2-fold split (same as UDID_Nonparametric) ----
  set.seed(seed)
  idx1     <- sample(which(D == 1), floor(sum(D == 1) / 2))
  idx0     <- sample(which(D == 0), ceiling(sum(D == 0) / 2))
  SS.Index <- list(c(idx1, idx0), setdiff(seq_len(n), c(idx1, idx0)))
  
  ## Accumulate per-unit EIF * pr_D across folds
  uncEIFtimesPrD <- numeric(n)
  
  for (ss in 1:2) {
    train <- SS.Index[[ss]]
    eval  <- SS.Index[[3 - ss]]
    n.eval <- length(eval)
    
    D.tr  <- D[train];  D.ev  <- D[eval]
    dY.tr <- deltaY[train]; dY.ev <- deltaY[eval]
    X.tr  <- covariates[train, , drop = FALSE]
    X.ev  <- covariates[eval,  , drop = FALSE]
    w.ev  <- w[eval]
    
    ## --- Propensity score: P(D=1 | X) ---
    Data_ps <- setNames(data.frame(cbind(D.tr, X.tr)), c("D", xnms))
    set.seed(seed * 10L + ss)
    SL_ps <- UDID::MySL(Data_ps, locY = 1, locX = 1 + 1:d,
                  Ydist = binomial(),
                  SL.list = SL.list)
    ps_ev <- pmin(pmax(
      predict(SL_ps, newdata = setNames(data.frame(X.ev), xnms),
              onlySL = TRUE)$pred,
      1e-6), 1 - 1e-6)
    
    ## --- Outcome regression: E[deltaY | X, D=0] ---
    ctrl_tr <- which(D.tr == 0)
    
    if (type == "binary") {
      ## For binary Y: fit E[Y1|X,D=0] and E[Y0|X,D=0] separately
      y1.tr <- (y1[train]); y0.tr <- (y0[train])
      
      Data_or1 <- setNames(
        data.frame(cbind(y1.tr[ctrl_tr], X.tr[ctrl_tr, , drop = FALSE])),
        c("dY", xnms))
      set.seed(seed * 10L + ss + 2L)
      SL_or1 <- UDID::MySL(Data_or1, locY = 1, locX = 1 + 1:d,
                     Ydist = binomial(),
                     SL.list = SL.list)
      
      Data_or0 <- setNames(
        data.frame(cbind(y0.tr[ctrl_tr], X.tr[ctrl_tr, , drop = FALSE])),
        c("dY", xnms))
      set.seed(seed * 10L + ss + 4L)
      SL_or0 <- UDID::MySL(Data_or0, locY = 1, locX = 1 + 1:d,
                     Ydist = binomial(),
                     SL.list = SL.list)
      
      newX_ev <- setNames(data.frame(X.ev), xnms)
      mu1_ev  <- predict(SL_or1, newdata = newX_ev, onlySL = TRUE)$pred
      mu0_0ev <- predict(SL_or0, newdata = newX_ev, onlySL = TRUE)$pred
      mu0_ev  <- mu1_ev - mu0_0ev
      
    } else {
      ## Continuous: fit E[deltaY | X, D=0] directly
      Data_or <- setNames(
        data.frame(cbind(dY.tr[ctrl_tr], X.tr[ctrl_tr, , drop = FALSE])),
        c("dY", xnms))
      set.seed(seed * 10L + ss + 2L)
      SL_or <- UDID::MySL(Data_or, locY = 1, locX = 1 + 1:d,
                    Ydist = gaussian(),
                    SL.list = SL.list)
      mu0_ev <- predict(SL_or,
                        newdata = setNames(data.frame(X.ev), xnms),
                        onlySL = TRUE)$pred
    }
    
    ## --- DR-DID EIF on eval fold ---
    pr_D.ev  <- mean(D.ev)
    ipw_wt   <- ps_ev / (1 - ps_ev)
    
    att_treat <- w.ev * D.ev * (dY.ev - mu0_ev)
    att_ctrl  <- w.ev * (1 - D.ev) * ipw_wt * (dY.ev - mu0_ev)
    
    # uncEIF * pr_D (same convention as .UDID_fold)
    uncEIFtimesPrD[eval] <- att_treat - att_ctrl
  }
  
  ## ---- Aggregate across folds ----
  pr_D <- mean(D)
  uncEIF <- uncEIFtimesPrD / pr_D
  ATT    <- mean(uncEIF)
  
  ## Influence function (centered)
  att.inf.func <- uncEIF - D * ATT / pr_D
  
  return(list(ATT = ATT, att.inf.func = as.numeric(att.inf.func)))
}
