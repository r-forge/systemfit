linear.hypothesis.systemfit <- function( model,
      hypothesis.matrix, rhs = NULL, test = c( "F", "Chisq" ),
      vcov. = NULL, ... ){

   thisCall <- match.call()
   test <- match.arg( test )

   result <- car:::linear.hypothesis.default( model,
      hypothesis.matrix = hypothesis.matrix, rhs = rhs, test = test,
      vcov. = vcov., ... )

   if( "model" %in% names( thisCall ) ) {
      modelPos <- grep( "^Model 1: model\n", attributes( result )$heading )
      if( class( thisCall$model ) == "name" ) {
         modelName <- as.character( thisCall$model )
      } else if( class( thisCall$model ) == "call" ) {
         modelName <- format( thisCall$model )
      } else {
         modelName <- thisCall$model
      }
      attributes( result )$heading[ modelPos[ 1 ] ] <-
         sub( "^Model 1: model\n",
            paste( "Model 1: ", modelName, "\n", sep = "" ),
            attributes( result )$heading[ modelPos[ 1 ] ] )
   }
   if( test == "Chisq" ){
      attributes( result )$heading[ 1 ] <-
         "Linear hypothesis test (Wald-test)\n\nHypothesis:"
   } else if ( test == "F" ) {
      ftest <- .ftest.systemfit( object = model,
         restrictions = hypothesis.matrix, restrict.rhs = rhs,
         vcov. = vcov. )
      attributes( result )$heading[ 1 ] <-
         "Linear hypothesis test (F-test)\n\nHypothesis:"
      colnames( result )[ 3:4 ] <- c( "F", "Pr(>F)" )
      result[ 2, 3 ] <- ftest$statistic
      result[ 2, 4 ] <- ftest$p.value
      if( ftest$nRestr != abs( result[ 2, 2 ] ) ){
         stop( "internal error: wrong degrees of freedom of test statistic" )
      }
      if( ftest$df.residual.sys != result[ 1, 1 ] ){
         stop( "internal error: wrong degrees of freedom of residuals" )
      }
   } else {
      stop( "unknown test statistic '", test, "'. Please use 'F' or 'Chisq'" )
   }

   return( result )
}

