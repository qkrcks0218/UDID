result <- read.csv("RESULT_Zika_Merge.csv")

ATT.UDID.0 <- median(result$Est_UDID_NP_X_01)
SE.UDID.0 <- sqrt(median(result$SE_UDID_NP_X_01^2 + (result$Est_UDID_NP_X_01-ATT.UDID.0)^2))

ATT.UDID.1 <- median(result$Est_UDID_NP_X_02)
SE.UDID.1 <- sqrt(median(result$SE_UDID_NP_X_02^2 + (result$Est_UDID_NP_X_02-ATT.UDID.1)^2))

ATT.UDID <- c(ATT.UDID.0, ATT.UDID.1)
SE.UDID <- c(SE.UDID.0, SE.UDID.1)

ATT.DID.0 <- median(result$Est_DID_NP_X_01)
SE.DID.0 <- sqrt(median(result$SE_DID_NP_X_01^2 + (result$Est_DID_NP_X_01-ATT.DID.0)^2))

ATT.DID.1 <- median(result$Est_DID_NP_X_02)
SE.DID.1 <- sqrt(median(result$SE_DID_NP_X_02^2 + (result$Est_DID_NP_X_02-ATT.DID.1)^2))

ATT.DID <- c(ATT.DID.0, ATT.DID.1)
SE.DID <- c(SE.DID.0, SE.DID.1)




a0 <- ATT.UDID
s0 <- SE.UDID
a1 <- ATT.DID
s1 <- SE.DID

R1 <- sprintf("Estimate & \\multicolumn{1}{c|}{%0.3f} & %0.3f & \\multicolumn{1}{c|}{%0.3f} & %0.3f \\\\ \\hline",a0[1],a0[2],a1[1],a1[2])
R2 <- sprintf("ASE & \\multicolumn{1}{c|}{%0.3f} & %0.3f & \\multicolumn{1}{c|}{%0.3f} & %0.3f \\\\ \\hline",s0[1],s0[2],s1[1],s1[2])
R3 <- sprintf("95\\%% CI & \\multicolumn{1}{c|}{(%0.3f,%0.3f)} & (%0.3f,%0.3f) & \\multicolumn{1}{c|}{(%0.3f,%0.3f)} & (%0.3f,%0.3f) \\\\ \\hline",
        a0[1]-qnorm(0.975)*s0[1],
        a0[1]+qnorm(0.975)*s0[1],
        a0[2]-qnorm(0.975)*s0[2],
        a0[2]+qnorm(0.975)*s0[2],
        a1[1]-qnorm(0.975)*s1[1],
        a1[1]+qnorm(0.975)*s1[1],
        a1[2]-qnorm(0.975)*s1[2],
        a1[2]+qnorm(0.975)*s1[2])

print(data.frame(rbind(R1,R2,R3)),row.names=F)


 