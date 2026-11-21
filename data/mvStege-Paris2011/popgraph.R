library(MASS);

testFunc <- function( x, y ) {
    z <- x+2*y
    z
}



popgraphGetEED <- function( data, pops, alpha ) {
    pops <- as.factor( pops )
    N <- dim(data)[1]
    npop <- length(levels(pops))
    R <- matrix(1,npop,npop)
    colMean <- matrix(0,npop,1)
    PWDist <- matrix(0.0,nrow=npop,ncol=npop)
    critVal <- qchisq( 1.0-alpha,1 )
    SRI <- matrix(1,npop,npop)
    EED <- matrix(0,npop,npop)
    EdgeStr <- matrix(0,npop,npop)
    pcfit <- prcomp(data,retx=T)
    genos <- as.matrix(pcfit$x[,pcfit$sdev > 1.0e-3])
    fit1 <- lda(x=genos,grouping=pops,method="moment")
    means <- colMeans(fit1$means)
    LDValues <- scale(genos,center=means,scale=FALSE) %*% fit1$scaling
    allLD <- by(LDValues,pops,mean)
    allSD <-  by(LDValues,pops,sd)
    popSD <- matrix(0.0,nrow=length(names(allSD)),ncol=1)
    for(i in seq(1,length(popSD)))
        popSD[i] <- sum(unlist(allSD[i]))
    popNames <- names( allSD )
    for(i in seq(1,npop)){
        for(j in seq(1,npop)){
            p1 <- unlist( allLD[i])
            p2 <- unlist( allLD[j])
            PWDist[i,j] = sqrt( sum( (p1-p2)^2 ))
        }
    }
    totMean <- mean(PWDist)
    for(i in 1:npop) colMean[i] <- mean(PWDist[,i])
    totMeanMatrix = matrix(1,npop,npop)*totMean
    colMeanMatrix = matrix(colMean,npop,npop,byrow=T)
    rowMeanMatrix = matrix(colMean,npop,npop,byrow=F)
    C <- -0.5 * (PWDist - colMeanMatrix - rowMeanMatrix + totMeanMatrix)
    for( i in 1:npop )
        for( j in 1:npop )
            if(i!=j)
                R[i,j] = C[i,j]/sqrt( C[i,i] * C[j,j] )
    RI <- ginv(R)
    for(i in seq(1,npop))
        for(j in seq(1,npop))
            if(i!=j)  
                SRI[i,j] <- -1*RI[i,j]/sqrt( RI[i,i]*RI[j,j])
    for(i in seq(1,npop))
        for(j in seq(1,npop))
            if(i!=j) {
                EED[i,j] = -N * log(1-SRI[i,j]^2)
                EdgeStr[i,j] = -0.5 * log(1-SRI[i,j]^2)
            }
    EED[EED<critVal]<-NA
    res <- list()
    res$PWDist <- PWDist
    res$EED <- EED
    res$popSD <- popSD
    res$popNames <- popNames
    res$npop <- npop
    res
}


getPopGraphAdjacencyMatrix <- function( data, pops, alpha ) {
    res <- popgraphGetEED( data, pops, alpha )
    A <- res$PWDist
    A[ is.na(res$EED) ] <- 0
    ret <- list()
    ret$popNames <- res$popNames
    ret$adjacency <- A
    ret
}

getFullPopGraph <- function(data, pops, color, alpha) {
        res <- popgraphGetEED( data, pops, alpha )
        numEdges <- length(res$EED[!is.na(res$EED)])/2.0
        cat(res$npop,numEdges,"\n")
        mnSd = min(res$popSD)
        mxSd = max(res$popSD)
        for(i in seq(1,res$npop)) 
            cat(res$popNames[i],(res$popSD[i]-mnSd)/mxSd/2+0.05, color[i],"\n")    
        for(i in seq(1,res$npop))
            for(j in seq(i+1,res$npop))
                if(j<=res$npop && !is.na(res$EED[i,j]))
                    cat(res$popNames[i],res$popNames[j],res$PWDist[i,j],"\n")
}

   
