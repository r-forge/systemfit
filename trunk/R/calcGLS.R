.calcXtOmegaInv <- function( xMat, sigma, nObsEq, invertSigma = TRUE,
      useMatrix = FALSE, solvetol = 1e-5 ){

   nEq <- length( nObsEq )

   if( useMatrix ){
      sigma <- as( sigma, "dspMatrix" )
      xMat <- as( xMat, "dgCMatrix" )
   }

   if( invertSigma ) {
      sigmaInv <- solve( sigma, tol = solvetol )
   } else {
      sigmaInv <- sigma
   }

   if( useMatrix ){
      result <- crossprod( xMat, kronecker( sigmaInv, Diagonal( nObsEq[ 1 ] ) ) )
   } else {
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
   }

   return( result )
}

.calcGLS <- function( xMat, yVec = NULL, xMat2 = xMat, R.restr = NULL,
      q.restr = NULL, sigma, nObsEq, useMatrix = TRUE, solvetol = 1e-5 ){

   if( useMatrix ){
      xMat  <- as( xMat,  "dgCMatrix" )
      xMat2 <- as( xMat2, "dgCMatrix" )
      sigma <- as( sigma, "dspMatrix" )
   }

   xtOmegaInv <- .calcXtOmegaInv( xMat = xMat, sigma = sigma, nObsEq = nObsEq,
      useMatrix = useMatrix, solvetol = solvetol )
   if( is.null( R.restr ) ) {
      result <- solve( xtOmegaInv %*% xMat2, tol = solvetol )
      if( !is.null( yVec ) ) {
         result <- result %*% xtOmegaInv %*% yVec
      }
   } else {
      W <- rbind2( cbind2( xtOmegaInv %*% xMat2, t(R.restr) ),
                  cbind2( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      Winv <- solve( W, tol=solvetol )
      if( is.null( yVec ) ) {
         result <- Winv[ 1:ncol(xMat), 1:ncol(xMat) ]
      } else{
         V <- c( as.numeric( xtOmegaInv %*% yVec ), q.restr )
         result <- ( Winv %*% V )[1:ncol( xMat ),]     # restricted coefficients
      }
   }
   return( result )
}

