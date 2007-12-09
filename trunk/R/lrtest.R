## Likelihood Ratio Test
lrtest.systemfit <- function( object, ... ) {

   thisCall <- match.call()

   object$lrtest.systemfit.name <- deparse( substitute( object ) )
   objectList <- list( ... )
   for( i in length( objectList ) ){
      objectList[[ i ]]$lrtest.systemfit.name <-
         paste( "object", i + 1, sep = "_" )
   }
   extractName <- function( object ){
      return( object$lrtest.systemfit.name )
   }

   result <- do.call( lrtest.default,
      c( list( object = object ), objectList, list( name = extractName ) ) )

   for( i in 2:nrow( result ) ){
      if( ( result[ i, "#Df" ] - result[ i - 1, "#Df" ] ) *
            ( result[ i, "LogLik" ] - result[ i - 1, "LogLik" ] ) < 0 ) {
         if( result[ i, "LogLik" ] > result[ i - 1, "LogLik" ] ) {
            compareLikelihood <- "larger"
            compareDf <- "less"
         } else {
            compareLikelihood <- "smaller"
            compareDf <- "more"
         }
         warning( "model '", i, "' has a ", compareLikelihood,
            " log-likelihood value than the ", compareDf,
            " restricted model '", i - 1, "'" )
      }
   }

   return( result )
}
