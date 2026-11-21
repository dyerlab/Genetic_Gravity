.packageName <- "GenCondInd"
#
# Does a canonical variate analysis on all baja data
#

library(MASS)


seqDF <- function(G,n){
  for(i in seq(1,n)){
    df <- (n-i+1)*(G-1)
    cat(i,df,"\n")
  }
}

dist2cov <- function( genDistMatrix ){
  N <- dim( genDistMatrix )[1]
  colMean <- matrix(0,N,1)
  totMean <- mean(genDistMatrix)

  for(i in 1:N) 
    colMean[i] <- mean(genDistMatrix[,i])

  totMeanMatrix = matrix(1,N,N)*totMean
  colMeanMatrix = matrix(colMean,N,N,byrow=T)
  rowMeanMatrix = matrix(colMean,N,N,byrow=F)
  C <- -0.5 * (genDistMatrix - colMeanMatrix - rowMeanMatrix + totMeanMatrix) 
  C
}

cov2dist <- function( covMat ) {
  N <- dim( covMat )[1]
  D <- matrix(0,N,N)
  for(i in 1:N)
    for(j in 1:N)
      D[i,j] <- covMat[i,i] + covMat[j,j] - 2*covMat[i,j]
  D
}

cov2cor <- function( covMat ) {
  N <- dim( covMat )[1]
  R <- matrix(1,N,N)
  for( i in 1:N )
    for( j in 1:N )
      if(i!=j)
        R[i,j] = covMat[i,j]/sqrt( covMat[i,i] * covMat[j,j] )
  R
}


runPG <- function(fname){	
  data <- read.table(fname,sep=",",quote="",header=T)
  pops <- as.factor(data$Strat0)
  genos1 <- as.matrix(data[,2:dim(data)[2]])
  pcfit <- prcomp(genos1,retx=T)
  genos <- as.matrix(pcfit$x[,pcfit$sdev > 1.0e-5])
  print(dim(genos))
  fit1 <- lda(x=genos,grouping=pops,method="moment")
  means <- colMeans(fit1$means)
  LDValues <- scale(genos,center=means,scale=FALSE) %*% fit1$scaling
  allLD <- by(LDValues,pops,mean)
  npop <- length(levels(pops))
  nLDAxes <- 20
  PWDist <- matrix(0.0,nrow=npop,ncol=npop)
  for(i in seq(1,npop)){
    for(j in seq(1,npop)){
      p1 <- unlist( allLD[i])
      p2 <- unlist( allLD[j])
      dist <- sqrt( sum( (p1-p2)^2 ))
      PWDist[i,j]=dist
    }
  }
  R <- cov2cor(dist2cov(PWDist))
  N <- length(pops)
  strata <- levels(pops)
  RI <- ginv(R)
  ni <- dim(R)[1]
  SRI <- matrix(1,ni,ni)
  for(i in seq(1,ni)){
    for(j in seq(1,ni)){
      if(i!=j)
        SRI[i,j] <- -1*RI[i,j]/sqrt( RI[i,i]*RI[j,j])
    }
  }
  EED <- matrix(0,ni,ni)
  EdgeStr <- matrix(0,ni,ni)
  for(i in seq(1,ni)){
    for(j in seq(1,ni)){
      if(i!=j) {
        EED[i,j] = -N * log(1-SRI[i,j]^2)
        EdgeStr[i,j] = -0.5 * log(1-SRI[i,j]^2)
      }
    }
  }
  EED[EED<3.84]<-NA
  EdgeStr[EED<3.84] <- NA

  p <- as.character(unlist(unique(pops)))
  cat(length(p),length(EED[!is.na(EED)])/2.0,"\n")
  allSD <-  by(LDValues,pops,sd)
  popSD <- matrix(0.0,nrow=length(p),ncol=1)
  for(i in seq(1,length(popSD)))
    popSD[i] <- sum(unlist(allSD[i]))

  for(i in seq(1,length(p))) 
    cat(p[i],(popSD[i]-min(popSD))/max(popSD)/2+0.05," 0\n")


  for(i in seq(1,ni)){
    for(j in seq(i+1,ni)){
      if(j<=ni){
        if(!is.na(EED[i,j])){
          cat(as.character(strata[i]),as.character(strata[j]),PWDist[i,j],"\n")

        }
      }
    }
  }   
  f <-  list("EED"=EED,"EdgeStr"=EdgeStr, "PWDist"=PWDist)
  f
}


runPGRaw <- function( rawGenos, pops ){	
  pops <- as.factor(pops)
  pcfit <- prcomp(rawGenos,retx=T)
  genos <- as.matrix(pcfit$x[,pcfit$sdev > 1.0e-5])
  print(dim(genos))
  fit1 <- lda(x=genos,grouping=pops,method="moment")
  means <- colMeans(fit1$means)
  LDValues <- scale(genos,center=means,scale=FALSE) %*% fit1$scaling
  allLD <- by(LDValues,pops,mean)
  npop <- length(levels(pops))
  nLDAxes <- 20
  PWDist <- matrix(0.0,nrow=npop,ncol=npop)
  for(i in seq(1,npop)){
    for(j in seq(1,npop)){
      p1 <- unlist( allLD[i])
      p2 <- unlist( allLD[j])
      dist <- sqrt( sum( (p1-p2)^2 ))
      PWDist[i,j]=dist
    }
  }
  R <- cov2cor(dist2cov(PWDist))
  N <- length(pops)
  strata <- levels(pops)
  RI <- ginv(R)
  ni <- dim(R)[1]
  SRI <- matrix(1,ni,ni)
  for(i in seq(1,ni)){
    for(j in seq(1,ni)){
      if(i!=j)
        SRI[i,j] <- -1*RI[i,j]/sqrt( RI[i,i]*RI[j,j])
    }
  }
  EED <- matrix(0,ni,ni)
  EdgeStr <- matrix(0,ni,ni)
  for(i in seq(1,ni)){
    for(j in seq(1,ni)){
      if(i!=j) {
        EED[i,j] = -N * log(1-SRI[i,j]^2)
        EdgeStr[i,j] = -0.5 * log(1-SRI[i,j]^2)
      }
    }
  }
  EED[EED<3.84]<-NA
  p <- as.character(unlist(unique(pops)))
  cat(length(p),length(EED[!is.na(EED)])/2.0,"\n")
  allSD <-  by(LDValues,pops,sd)
  popSD <- matrix(0.0,nrow=length(p),ncol=1)
  for(i in seq(1,length(popSD)))
    popSD[i] <- sum(unlist(allSD[i]))

  for(i in seq(1,length(p))) 
    cat(p[i],(popSD[i]-min(popSD))/max(popSD)/2+0.05," 0\n")


  for(i in seq(1,ni)){
    for(j in seq(i+1,ni)){
      if(j<=ni){
        if(!is.na(EED[i,j])){
          cat(as.character(strata[i]),as.character(strata[j]),PWDist[i,j],"\n")

        }
      }
    }
  }
  
   f <-  list("EED"=EED,"EdgeStr"=EdgeStr, "PWDist"=PWDist)
    f

}


