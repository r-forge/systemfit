.calcXtOmegaInv <- function( xMat, sigma, validObsEq, invertSigma = TRUE,
      useMatrix = FALSE, warnMatrix = TRUE, solvetol = 1e-5 ){

   nEq <- ncol( validObsEq )

   if( useMatrix && warnMatrix ){
      if( !inherits( sigma, "dsyMatrix" ) ){
         warning( "class of 'sigma' is '", 
            paste( class( sigma ), collapse = ", " ),
            "', but it should be 'dsyMatrix'" )
      }
      if( !inherits( xMat,"dgCMatrix" ) ){
         warning( "class of 'xMat' is '", 
            paste( class( xMat ), collapse = ", " ),
            "', but it should be 'dgCMatrix'" )
      }
   }

   if( invertSigma ) {
      sigmaInv <- solve( sigma, tol = solvetol )
   } else {
      sigmaInv <- sigma
   }

   if( useMatrix ){
      for( i in 1:nEq ) {
         for( j in 1:nEq ) {
            thisBlock <- sparseMatrix(
               i = which( validObsEq[ validObsEq[ , i ], j ] ),
               j = which( validObsEq[ validObsEq[ , j ], i ] ),
               x = sigmaInv[ i, j ],
               dims = c( sum( validObsEq[ , i ] ), sum( validObsEq[ , j ] ) ) )
            if( j == 1 ) {
               thisRow <- thisBlock
            } else {
               thisRow <- cbind( thisRow, thisBlock )
            }
         }
         if( i == 1 ) {
            omegaInv <- thisRow
         } else {
            omegaInv <- rbind( omegaInv, thisRow )
         }
      }
      result <- crossprod( xMat, omegaInv )
   } else {
      eqSelect <- rep( 0, nrow( xMat ) )
      for( i in 1:nEq ) {
         eqSelect[ ( sum( validObsEq[ , 0:( i - 1 ) ] ) + 1 ):sum( validObsEq[ , 1:i ] ) ] <- i
      }
      result <- matrix( 0, nrow = ncol( xMat ), ncol = nrow( xMat ) )
      for( i in 1:nEq ) {
         for( j in 1:nEq ) {
            colSelectI <- eqSelect == i
            colSelectI[ colSelectI ] <- validObsEq[ validObsEq[ , i ], j ]
            colSelectJ <- eqSelect == j
            colSelectJ[ colSelectJ ] <- validObsEq[ validObsEq[ , j ], i ]
            result[ , colSelectI ] <- result[ , colSelectI ] +
               t( xMat )[ , colSelectJ ] * sigmaInv[ i, j ]
         }
      }
   }

   return( result )
}

.calcGLS <- function( xMat, yVec = NULL, xMat2 = xMat, R.restr = NULL,
      q.restr = NULL, sigma, validObsEq, useMatrix = TRUE, warnMatrix = TRUE,
      solvetol = 1e-5 ){

   if( useMatrix && warnMatrix ){
      if( !inherits( xMat, "dgCMatrix" ) ){
         warning( "class of 'xMat' is '", 
            paste( class( xMat ), collapse = ", " ),
            "', but it should be 'dgCMatrix'" )
      }
      if( !inherits( xMat2, "dgCMatrix" ) ){
         warning( "class of 'xMat2' is '", 
            paste( class( xMat2 ), collapse = ", " ),
            "', but it should be 'dgCMatrix'" )
      }
      if( !inherits( sigma, "dsyMatrix" ) ){
         warning( "class of 'sigma' is '", 
            paste( class( sigma ), collapse = ", " ),
            "', but it should be 'dsyMatrix'" )
      }
   }

   xtOmegaInv <- .calcXtOmegaInv( xMat = xMat, sigma = sigma, validObsEq = validObsEq,
      useMatrix = useMatrix, solvetol = solvetol )
   if( is.null( R.restr ) ) {
      if( is.null( yVec ) ) {
         result <- solve( xtOmegaInv %*% xMat2, tol = solvetol )
      } else {
         result <- as.numeric( solve( xtOmegaInv %*% xMat2, xtOmegaInv %*% yVec,
            tol = solvetol ) )
      }
   } else {
      W <- rbind2( cbind2( xtOmegaInv %*% xMat2, t(R.restr) ),
                  cbind2( R.restr, matrix(0, nrow(R.restr), nrow(R.restr) )))
      if( is.null( yVec ) ) {
         result <- as.matrix(
            solve( W, tol=solvetol )[ 1:ncol(xMat), 1:ncol(xMat) ] )
      } else{
         V <- c( as.numeric( xtOmegaInv %*% yVec ), q.restr )
         result <- solve( W, V, tol=solvetol )[ 1:ncol( xMat ) ]
      }
   }
   return( result )
}

