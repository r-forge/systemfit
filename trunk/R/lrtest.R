## Likelihood Ratio Test
lrtest.systemfit <- function( object, ... ) {

   thisCall <- match.call()

   object$name <- deparse( substitute( object ) )

   extractName <- function( object ){
      if( !exists( ".lrtestSystemfitNameNumber" ) ) {
         .lrtestSystemfitNameNumber <<- 1
      } else {
         .lrtestSystemfitNameNumber <<- .lrtestSystemfitNameNumber + 1
      }
      objectName <- object$name
      if( is.null( objectName ) ) {
         objectName <- paste( "object", .lrtestSystemfitNameNumber, sep = "_" )
      }
      return( objectName )
   }

   result <- lrtest.default( object = object, ..., name = extractName )

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

   rm( .lrtestSystemfitNameNumber, inherits = TRUE )

   return( result )
}
