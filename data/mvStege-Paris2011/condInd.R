#
# Does a canonical variate analysis on all baja data
#

rm(list=ls())

library(mva)
library(MASS)
library(lqs)

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




C <- dist2cov(PWDist)
#cat("COV\n")
#print(C)


R <- cov2cor(C)
#cat("CORR\n")
#print(R)

N <- length(pops)
strata <- names(LD1)

RI <- ginv(R)
#cat("Inverse Correlation Matrix\n")

#print(RI)
ni <- dim(R)[1]

SRI <- matrix(1,ni,ni)
for(i in seq(1,ni)){
  for(j in seq(1,ni)){
    if(i!=j)
      SRI[i,j] <- -1*RI[i,j]/sqrt( RI[i,i]*RI[j,j])
  }
}

#cat("Standardized Inverse Correlation Matrix\n");
#print(SRI)


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
cat("Edge Exclustion Deviance\n");
print(EED)


p <- as.character(unlist(unique(pops)))
cat(length(p),length(EED[!is.na(EED)])/2.0,"\n")
allSD <-  by(LDValues,pops,sd)
popSD <- matrix(0.0,nrow=length(p),ncol=1)
for(i in seq(1,length(popSD)))
  popSD[i] <- sum(unlist(allSD[i]))

for(i in seq(1,length(p))) 
  cat(p[i],(popSD[i]-min(popSD))/max(popSD)/2+0.05,"\n")


for(i in seq(1,ni)){
  for(j in seq(i+1,ni)){
    if(j<=ni){
      if(!is.na(EED[i,j])){
        cat(as.character(strata[i]),as.character(strata[j]),PWDist[i,j],"\n")

      }
    }
  }
}



