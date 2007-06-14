## Likelihood Ratio Test
lrtest.systemfit <- function( object, ... ) {

   thisCall <- match.call()

   if( "object" %in% names( thisCall ) ) {
      if( class( thisCall$object ) == "name" ) {
         object$name <- as.character( thisCall$object )
      } else if( class( thisCall$object ) == "call" ) {
         object$name <- format( thisCall$object )
      } else {
         object$name <- thisCall$object
      }
   }

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

   rm( .lrtestSystemfitNameNumber, inherits = TRUE )

   return( result )
}
