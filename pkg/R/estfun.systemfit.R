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
   
   mmAll <- matrix( NA, nrow = length( res ), ncol = ncol( mm ) )
   mmAll[ !is.na( res ), ] <- mm
      
   result <- res * mmAll
   
   result <- result[ !is.na( res ), ]
   dimnames( result ) <- dimnames( mm )
   
   return( result )
}