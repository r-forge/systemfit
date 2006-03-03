.calcGLS <- function( x, y = NULL, x2 = x, R.restr = NULL, q.restr = NULL, 
      sigma, nObsEq, solvetol = 1e-5 ){
   nEq <- length( nObsEq )
   sigmaInv <- solve( sigma, tol = solvetol )
   # omegaInv <- sigmaInv %x% diag( 1 ,nObsEq[1], nObsEq[1] )
   # result <- solve( t( x ) %*% omegaInv %*% x, tol = solvetol ) %*%
   #    t( x ) %*% omegaInv %*% y
   eqSelect <- rep( 0, nrow( x ) )
   for( i in 1:nEq ) {
      eqSelect[ ( sum( nObsEq[ 0:( i - 1 ) ] ) + 1 ):sum( nObsEq[ 1:i ] ) ] <- i
   }
   xtOmegaInv <- matrix( 0, nrow = ncol( x ), ncol = nrow( x ) )
   for( i in 1:nEq ) {
      for( j in 1:nEq ) {
         xtOmegaInv[ , eqSelect == i ] <-  xtOmegaInv[ , eqSelect == i ] +
            t( x )[ , eqSelect == j ] * sigmaInv[ i, j ]
      }
   }
   if( is.null( R.restr ) ) {
      result <- solve( xtOmegaInv %*% x2, tol = solvetol )
      if( !is.null( y ) ) {
         result <- result %*% xtOmegaInv %*% y
      }
   } else {
      W <- rbind( cbind( xtOmegaInv %*% x2, t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      Winv <- solve( W, tol=solvetol )
      if( is.null( y ) ) {
         result <- Winv[ 1:ncol(x), 1:ncol(x) ]
      } else{
         V <- rbind( xtOmegaInv %*% y , q.restr )
         result <- ( Winv %*% V )[1:ncol( x )]     # restricted coefficients
      }
   }
   return( result )
}