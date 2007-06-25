linear.hypothesis.systemfit <- function( model,
      hypothesis.matrix, rhs = NULL, test = c( "F", "Chisq" ),
      vcov. = NULL, ... ){

   thisCall <- match.call()
   test <- match.arg( test )

   if( "model" %in% names( thisCall ) ) {
      if( class( thisCall$model ) == "name" ) {
         modelName <- as.character( thisCall$model )
      } else if( class( thisCall$model ) == "call" ) {
         modelName <- format( thisCall$model )
      } else {
         modelName <- thisCall$model
      }
   }
   if( test == "Chisq" ){
      result <- car:::linear.hypothesis.default( model,
         hypothesis.matrix = hypothesis.matrix, rhs = rhs, test = test,
         vcov. = vcov., ... )

      attributes( result )$heading[ 1 ] <-
         "Linear hypothesis test (Wald-test)\n\nHypothesis:"

      modelPos <- grep( "^Model 1: model\n", attributes( result )$heading )
      attributes( result )$heading[ modelPos[ 1 ] ] <-
         sub( "^Model 1: model\n",
            paste( "Model 1: ", modelName, "\n", sep = "" ),
            attributes( result )$heading[ modelPos[ 1 ] ] )

   } else if ( test == "F" ) {
      result <- as.data.frame( matrix( NA, nrow = 2, ncol = 4 ) )
      names( result ) <- c( "Res.Df", "Df", "F", "Pr(>F)" )
      rownames( result ) <- c( 1, 2 )

      if( is.null( rhs ) ){
         rhs <- rep( 0, nrow( hypothesis.matrix ) )
      }

      ftest <- .ftest.systemfit( object = model,
         restrictions = hypothesis.matrix, restrict.rhs = rhs,
         vcov. = vcov. )
      result[ 1, 1 ] <- ftest$df.residual.sys
      result[ 2, 1 ] <- ftest$df.residual.sys + ftest$nRestr
      result[ 2, 2 ] <- result[ 1, 1 ] - result[ 2, 1 ]
      result[ 2, 3 ] <- ftest$statistic
      result[ 2, 4 ] <- ftest$p.value

      title <- "Linear hypothesis test (F-test)\n\nHypothesis:"
      topnote <- paste( "Model 1: ", modelName,
         "\nModel 2: restricted model", sep = "" )
      if( is.null( vcov. ) ){
         note <- ""
      } else {
         note <- "\nNote: Coefficient covariance matrix supplied.\n"
      }
      attributes( result )$heading <- c( title,
         car:::printHypothesis( hypothesis.matrix, rhs, names( coef( model ) ) ),
         "", topnote, note )
      class( result ) <- c( "anova", "data.frame" )
   } else {
      stop( "unknown test statistic '", test, "'. Please use 'F' or 'Chisq'" )
   }

   return( result )
}
