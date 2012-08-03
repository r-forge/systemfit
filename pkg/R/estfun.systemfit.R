estfun.systemfit <- function ( obj, ... ) {
   if( !is.null( obj$restrict.matrix ) || !is.null( obj$restrict.rhs ) ||
        !is.null( obj$restrict.regMat ) ) {
      stop( "returning the estimation function for models with restrictions",
            " has not yet been implemented.")
   }
   
   # residuals
   res <- unlist(  residuals( obj ) )
   
   # model matrix
   mm <- model.matrix( obj )
   
   if( sum( !is.na( res ) ) != nrow( mm ) ) {
      stop( "internal error: the number of residuals is not equal to the",
         " number of rows of the model matrix. Please contact the maintainer." )
   }

   if( is.null( obj$residCovEst ) ) {
      omegaInvXmat <- mm
   } else {
      omegaInvXmat <- t( .calcXtOmegaInv( xMat = mm, sigma = obj$residCovEst, 
         validObsEq = !is.na( residuals( obj ) ), invertSigma = TRUE ) )
   }
   
   result <- res[ !is.na( res ) ] * omegaInvXmat
   
   dimnames( result ) <- dimnames( mm )
   
   return( result )
}