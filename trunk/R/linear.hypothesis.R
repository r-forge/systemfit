linear.hypothesis.systemfit <- function( model,
      hypothesis.matrix, rhs = NULL, test = c( "F", "Chisq" ), ... ){

   thisCall <- match.call()

   result <- car:::linear.hypothesis.default( model,
      hypothesis.matrix = hypothesis.matrix, rhs = rhs, test = test, ... )

   if( test == "Chisq" ){
      attributes( result )$heading[ 1 ] <-
         "Linear hypothesis test (Wald-test)\n\nHypothesis:"
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
# attributes( result )$thisCall <- thisCall
      }
   } else if ( test == "F" ) {
      stop( "F-test with linear.hypothesis is not implemented yet" )
   } else {
      stop( "unknown test statistic '", test, "'. Please use 'F' or 'Chisq'" )
   }

   return( result )
}

