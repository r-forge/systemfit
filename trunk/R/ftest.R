ftest.systemfit <- function( object, R.restr,
   q.restr = rep( 0, nrow( R.restr ) ) ){

   coef <- coef( object )
   vcov <- vcov( object )
   resid <- unlist( residuals( object ) )
   rcov  <- object$rcov
   nEq   <- object$g
   nObsEq <- object$n / nEq

   result <- list()

   result$nRestr <- nrow( R.restr )
   result$dfSys  <- nEq * nObsEq - length( coef )

   numerator <- t( R.restr %*% coef - q.restr ) %*%
      solve( R.restr %*% vcov %*% t( R.restr ) ) %*%
      ( R.restr %*% coef - q.restr )

   denominator <- t( resid ) %*%
      ( solve( rcov ) %x% diag( nObsEq ) ) %*%
      resid

   result$statistic <- ( numerator / result$nRestr ) /
      ( denominator / result$dfSys )

   result$p.value <- 1 - pf( result$statistic, result$nRestr, result$dfSys )

   class( result ) <- "ftest.systemfit"
   return( result )
}

print.ftest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "F-test for linear parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "F-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n" )
   invisible( x )
}