ftest.systemfit <- function( object, restrictions,
   restrict.rhs = rep( 0, nrow( restrictions ) ) ){

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

   result$nRestr <- nrow( restrictions )
   result$df.residual.sys  <- object$df.residual

   numerator <- t( restrictions %*% coef - restrict.rhs ) %*%
      solve( restrictions %*% vcov %*% t( restrictions ) ) %*%
      ( restrictions %*% coef - restrict.rhs )

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