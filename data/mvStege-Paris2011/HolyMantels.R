
#   matrix2Vector( X )
#
# This function takes an (N,N) matrix and returns a vector of length (N^2)
#
matrix2vector <- function( X ) {
    N = dim(X)[1]
    V <- matrix( 0, nrow=N*N, ncol=1)
    ctr <- 1
    for( i in seq(1,N) ) {
        for( j in seq( 1, N )) {
            V[ctr] <- X[i,j]
            ctr <- ctr+1
        }
    }
    V
}



