library(png)

addImg <- function( obj,  x = NULL,  y = NULL,  width = NULL,  interpolate = TRUE ){ 
  USR <- par()$usr 
  PIN <- par()$pin 
  DIM <- dim(obj) 
  ARp <- DIM[1]/DIM[2] 
  WIDi <- width/(USR[2]-USR[1])*PIN[1] 
  HEIi <- WIDi * ARp 
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) 
  rasterImage(image = obj, 
              xleft = x, xright = x+(width),
              ybottom = y-0.5*(HEIu), ytop = y+0.5*(HEIu), 
              interpolate = interpolate)
}

addImg2 <- function( obj,  x = NULL,  y = NULL,  height = NULL,  interpolate = TRUE ){ 
  USR <- par()$usr 
  PIN <- par()$pin 
  DIM <- dim(obj) 
  ARp <- DIM[1]/DIM[2] 
  HEIi <- height/(USR[4]-USR[3])*PIN[2]
  WIDi <- HEIi / ARp 
  WIDi <- WIDi/PIN[1]*(USR[2]-USR[1]) 
  rasterImage(image = obj, 
              xleft = x, xright = x+(WIDi),
              ybottom = y-0.5*(height), ytop = y+0.5*(height), 
              interpolate = interpolate)
}

Title1 <- readPNG("Title1.png")
Title2 <- readPNG("Title2.png")

N.Grid <- c(500,1000,2000,4000)

FL   <- "Final"  

Estimand <- c(0,0)

Eff.Obs <- Estimand
Eff.Exp <- Estimand
AXIS2 <- -7.5

WANTTOSEE <- "Eff"

TT <- 1:4

Est <- SE <- SE.B <- 
  Est.Oracle <- SE.Oracle <- SE.B.Oracle <- 
  Est.NoU <- SE.NoU <- SE.B.NoU <- list()

N.Text <- c(sprintf("N=%s",N.Grid[TT]))

XL <- c(-0.5,max(TT)+0.5)
YL <- list()
YL[[1]] <- c(-0.9,0.6)
YL[[2]] <- c(-0.3,0.2)
shift <- 0.125
boxwidth <- 0.175
COL <- list()
COL[[1]] <- rgb(0,0,0,0.5)
COL[[2]] <- rgb(0,0,0,0.25)
COL[[3]] <- rgb(1,0,0,0.25)
XP1 <- 1.2
XP2 <- 3.0
TW  <- 1.2
YP <- seq(0.9,0.1,length=6)

RESULT.Conti <- read.csv("Results/RESULT_Continuous_F.csv")
RESULT.Conti.N <- list()
RESULT.Conti.N[[1]] <- RESULT.Conti[RESULT.Conti$N==N.Grid[1],]
RESULT.Conti.N[[2]] <- RESULT.Conti[RESULT.Conti$N==N.Grid[2],]
RESULT.Conti.N[[3]] <- RESULT.Conti[RESULT.Conti$N==N.Grid[3],]
RESULT.Conti.N[[4]] <- RESULT.Conti[RESULT.Conti$N==N.Grid[4],]

RESULT.Binary <- read.csv("Results/RESULT_Binary_F.csv")
RESULT.Binary.N <- list()
RESULT.Binary.N[[1]] <- RESULT.Binary[RESULT.Binary$N==N.Grid[1],]
RESULT.Binary.N[[2]] <- RESULT.Binary[RESULT.Binary$N==N.Grid[2],]
RESULT.Binary.N[[3]] <- RESULT.Binary[RESULT.Binary$N==N.Grid[3],]
RESULT.Binary.N[[4]] <- RESULT.Binary[RESULT.Binary$N==N.Grid[4],]


Comp.Time <- cbind( RESULT.Conti.N[[1]]$Time_UDID, 
                    RESULT.Conti.N[[1]]$Time_DID,
                    RESULT.Conti.N[[2]]$Time_UDID, 
                    RESULT.Conti.N[[2]]$Time_DID,
                    RESULT.Conti.N[[3]]$Time_UDID, 
                    RESULT.Conti.N[[3]]$Time_DID,
                    RESULT.Conti.N[[4]]$Time_UDID, 
                    RESULT.Conti.N[[4]]$Time_DID,
                    
                    RESULT.Binary.N[[1]]$Time_UDID, 
                    RESULT.Binary.N[[1]]$Time_DID,
                    RESULT.Binary.N[[2]]$Time_UDID, 
                    RESULT.Binary.N[[2]]$Time_DID,
                    RESULT.Binary.N[[3]]$Time_UDID, 
                    RESULT.Binary.N[[3]]$Time_DID,
                    RESULT.Binary.N[[4]]$Time_UDID, 
                    RESULT.Binary.N[[4]]$Time_DID )

round( cbind( matrix(apply(Comp.Time,2,mean)[1:8],4,2,byrow=T) , 
              matrix(apply(Comp.Time,2,mean)[8+1:8],4,2,byrow=T) ), 2)


















png("Simulation.png",height=4.5,width=11,unit="in",res=500)

MAR    <- c(1,0,0.5,0.5)
layout(rbind(c(1,1),
             c(2,3),
             c(4,6),
             c(8,8),
             c(5,7)),
       heights=c(2, 1.5, 10, 2, 5))

par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(Title1, 
        x = 0, 
        y = 0.5, 
        height = 0.6)


par(mar=MAR*c(0,1,0,1))

MAR    <- c(1,0,1,0.5)
par(mar=MAR*c(0,1,0,1))
plot.new()
text(0.5,0.3,"Continuous Y",cex=1.4,font=1)
plot.new()
text(0.5,0.3,"Binary Y",cex=1.4,font=1)

MAR    <- c(1,0,0.5,0.5)

for(TYPE in 1:2){
  
  for(tt in TT){
    
    if(TYPE==1){
      Result.temp <- RESULT.Conti.N[[tt]]
      True_ATT <- 0
      
    } else {
      Result.temp <- RESULT.Binary.N[[tt]]
      True_ATT <- 0
    }
    
    
    print(nrow(Result.temp)) 
    
    Est[[tt]] <- Result.temp$UDID_ATT-True_ATT
    SE[[tt]]  <- Result.temp$UDID_ASE
    SE.B[[tt]]  <- Result.temp$UDID_BSE 
    
    Est.Oracle[[tt]] <- Result.temp$UDID_Oracle_ATT
    SE.Oracle[[tt]]  <- Result.temp$UDID_Oracle_ASE
    SE.B.Oracle[[tt]]  <- Result.temp$UDID_Oracle_BSE
    
    Est.NoU[[tt]] <- Result.temp$DID_ATT-True_ATT
    SE.NoU[[tt]]  <- Result.temp$DID_ASE 
    
  }
  
  par(mar=MAR*c(0,1,1,1))
  for(tt in TT[1]){
    boxplot(Est[[tt]],at=TT[1]-shift,
            xlab="",ylab="",axes=FALSE,
            xlim=XL,
            ylim=YL[[TYPE]],
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[1]],
            pch=19,cex=0.2)
    boxplot(Est.Oracle[[tt]],at=TT[1],
            xlab="",ylab="",axes=FALSE,add=TRUE,
            xlim=XL,
            ylim=YL[[TYPE]],
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[2]],
            pch=19,cex=0.2)
    boxplot(Est.NoU[[tt]],at=TT[1]+shift,
            xlab="",ylab="",axes=FALSE,add=TRUE,
            xlim=XL,
            ylim=YL[[TYPE]],
            medlwd = 1,
            boxwex = boxwidth,
            col=COL[[3]],
            pch=19,cex=0.2)
    
  }
  if(length(TT)>1){
    for(tt in TT[-1]){
      boxplot(Est[[tt]],at=tt-shift,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL[[TYPE]],
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[1]],
              pch=19,cex=0.2)
      boxplot(Est.Oracle[[tt]],at=tt,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL[[TYPE]],
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[2]],
              pch=19,cex=0.2)
      boxplot(Est.NoU[[tt]],at=tt+shift,
              xlab="",ylab="",axes=FALSE,add=TRUE,
              xlim=XL,
              ylim=YL[[TYPE]],
              medlwd = 1,
              boxwex = boxwidth,
              col=COL[[3]],
              pch=19,cex=0.2)
    }
    
  }
  
  if(TYPE==1){
    axis(2,line=AXIS2,at=round(seq(-0.9,0.6,by=0.3),1),cex.axis=0.8)
    
  } else {
    axis(2,line=AXIS2,at=round(seq(-0.3,0.2,by=0.1),1),cex.axis=0.8) 
  }
  title(ylab="Bias",line=-4.5)
  segments(0.5,0,max(TT)+0.5,0,col=c(2),lty=2)
  axis(3,col=1,at=(1:8)[TT],labels=N.Text,line=-2,tick=F)
  
  
  Summary <- function(tt){
    c(mean( Est[[tt]]  )*10,
      sd( Est[[tt]] )*10,
      mean( SE[[tt]] )*10, 
      mean( SE.Oracle[[tt]] )*10, 
      mean( SE[[tt]] )/mean( SE.Oracle[[tt]] ),
      mean( abs( Est[[tt]]  )/SE[[tt]] <= qnorm(0.975) ) 
    )
  }
  
  boxplot(cbind(rep(0,2),rep(0,2),rep(0,2),rep(0,2)),
          ylim=c(10,11),bty="n",xlim=XL,
          axes=F)
  text(0.4,10+YP[1],"Bias (x10)",pos=2,cex=1)
  text(0.4,10+YP[2],"ESE (x10)",pos=2,cex=1)
  text(0.4,10+YP[3],"ASE (x10)",pos=2,cex=1) 
  text(0.4,10+YP[4],"SEB (x10)",pos=2,cex=1)
  text(0.4,10+YP[5],"ASE/SEB",pos=2,cex=1)
  text(0.4,10+YP[6],"Coverage (ASE)",pos=2,cex=1) 
  for(tt in TT){
    vvv <- sapply(tt:tt,Summary)
    text(tt,10+YP[1],sprintf("%0.3f",vvv[1]),cex=1)
    text(tt,10+YP[2],sprintf("%0.3f",vvv[2]),cex=1)
    text(tt,10+YP[3],sprintf("%0.3f",vvv[3]),cex=1)
    text(tt,10+YP[4],sprintf("%0.3f",vvv[4]),cex=1)
    text(tt,10+YP[5],sprintf("%0.3f",vvv[5]),cex=1) 
    text(tt,10+YP[6],sprintf("%0.3f",vvv[6]),cex=1) 
  }
  
}





par(mar=MAR*c(0,1,0,1))
plot.new()
addImg2(Title2, 
        x = 0, 
        y = 0.5, 
        height = 0.6)




dev.off()


