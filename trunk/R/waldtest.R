waldtest.systemfit <- function( object, restrictions,
   restrict.rhs = rep( 0, nrow( restrictions ) ) ){

   coef <- coef( object )
   vcov <- vcov( object )

   result <- list()

   result$nRestr <- nrow( restrictions )

   result$statistic <- t( restrictions %*% coef - restrict.rhs ) %*%
      solve( restrictions %*% vcov %*% t( restrictions ) ) %*%
      ( restrictions %*% coef - restrict.rhs )

   result$p.value <- 1 - pchisq( result$statistic, result$nRestr )

   class( result ) <- "waldtest.systemfit"
   return( result )
}

print.waldtest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "Wald-test for linear parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "Wald-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "degrees of freedom:", x$nRestr, "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n\n" )
   invisible( x )
}