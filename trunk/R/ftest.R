ftest.systemfit <- function( object, R.restr,
   q.restr = rep( 0, nrow( R.restr ) ) ){

   coef <- coef( object )
   vcov <- vcov( object )
   resid <- unlist( residuals( object ) )
   nEq   <- length( object$eq )
   nObsPerEq <- nrow( residuals( object ) )
   if( is.null( object$rcovest ) ) {
      rcov <- diag( nEq )
   } else {
      rcov <- object$rcovest
   }

   result <- list()

   result$nRestr <- nrow( R.restr )
   result$df.residual.sys  <- object$df.residual

   numerator <- t( R.restr %*% coef - q.restr ) %*%
      solve( R.restr %*% vcov %*% t( R.restr ) ) %*%
      ( R.restr %*% coef - q.restr )

   denominator <- t( resid ) %*%
      ( solve( rcov ) %x% diag( nObsPerEq ) ) %*%
      resid
   #print( denominator )

   result$statistic <- ( numerator / result$nRestr ) /
      ( denominator / result$df.residual.sys )

   result$p.value <- 1 - pf( result$statistic, result$nRestr, result$df.residual.sys )

   class( result ) <- "ftest.systemfit"
   return( result )
}

print.ftest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "F-test for linear parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "F-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "degrees of freedom of the numerator:", x$nRestr, "\n" )
   cat( "degrees of freedom of the denominator:", x$df.residual.sys, "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n\n" )
   invisible( x )
}