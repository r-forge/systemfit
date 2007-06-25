.ftest.systemfit <- function( object, restrictions,
   restrict.rhs = NULL, vcov. = NULL ){

   coef <- coef( object )

   # coefficient covariance matrix
   if( is.null( vcov. ) ){
      vcov <- vcov( object )
   } else if( is.function( vcov. ) ){
      vcov <- vcov.( object )
   } else {
      vcov <- vcov.
   }

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

   return( result )
}
