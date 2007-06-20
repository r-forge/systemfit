.calcXtOmegaInv <- function( xMat, sigma, nObsEq, solvetol = 1e-5, invertSigma = TRUE ){
   nEq <- length( nObsEq )
   if( invertSigma ) {
      sigmaInv <- solve( sigma, tol = solvetol )
   } else {
      sigmaInv <- sigma
   }
   eqSelect <- rep( 0, nrow( xMat ) )
   for( i in 1:nEq ) {
      eqSelect[ ( sum( nObsEq[ 0:( i - 1 ) ] ) + 1 ):sum( nObsEq[ 1:i ] ) ] <- i
   }
   result <- matrix( 0, nrow = ncol( xMat ), ncol = nrow( xMat ) )
   for( i in 1:nEq ) {
      for( j in 1:nEq ) {
         result[ , eqSelect == i ] <- result[ , eqSelect == i ] +
            t( xMat )[ , eqSelect == j ] * sigmaInv[ i, j ]
      }
   }
   return( result )
}

.calcGLS <- function( xMat, yVec = NULL, xMat2 = xMat, R.restr = NULL, q.restr = NULL,
      sigma, nObsEq, solvetol = 1e-5 ){
   xtOmegaInv <- .calcXtOmegaInv( xMat = xMat, sigma = sigma, nObsEq = nObsEq,
      solvetol = solvetol )
   if( is.null( R.restr ) ) {
      result <- solve( xtOmegaInv %*% xMat2, tol = solvetol )
      if( !is.null( yVec ) ) {
         result <- result %*% xtOmegaInv %*% yVec
      }
   } else {
      W <- rbind( cbind( xtOmegaInv %*% xMat2, t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      Winv <- solve( W, tol=solvetol )
      if( is.null( yVec ) ) {
         result <- Winv[ 1:ncol(xMat), 1:ncol(xMat) ]
      } else{
         V <- rbind( xtOmegaInv %*% yVec , q.restr )
         result <- ( Winv %*% V )[1:ncol( xMat ),]     # restricted coefficients
      }
   }
   return( result )
}

