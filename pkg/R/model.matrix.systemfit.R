## return model matrix of the entire system
model.matrix.systemfit <- function( object, ... ){
   result <- matrix( NA, 0, 0 )
   mmRowNames <- NULL
   mmColNames <- NULL
   for( i in 1:length( object$eq ) ) {
      mmi <- model.matrix( object$eq[[ i ]] ) 
      result <- rbind(
         cbind( result, matrix( 0, nrow( result ), ncol( mmi ) ) ),
         cbind( matrix( 0, nrow( mmi ), ncol( result ) ), mmi ) )
      mmRowNames <- c( mmRowNames,
         paste( object$eq[[ i ]]$eqnLabel, "_", rownames( mmi ), sep = "" ) )
      for( j in 1:ncol( mmi ) ){
         cName <- colnames( mmi )[ j ]
         if( object$panelLike && cName != "(Intercept)" ){
            mmColNames <- c( mmColNames, cName )
         } else {
            mmColNames <- c( mmColNames,
               paste( object$eq[[ i ]]$eqnLabel, "_", cName, sep = "" ) )
         }
      }
   }
   rownames( result ) <- mmRowNames
   colnames( result ) <- mmColNames
   return( result )
}

## return model matrix of a single equation
model.matrix.systemfit.equation <- function( object, ... ){
   if( !is.null( object$x ) ) {
      result <- object$x
   } else if( !is.null( model.frame( object ) ) ) {
      result <- model.matrix( object$terms, data = model.frame( object ) )
      attrAssign <- attributes( result )$assign
      result <- result[ !is.na( residuals( object ) ), , drop = FALSE ]
      attributes( result )$assign <- attrAssign
   } else {
      stop( "returning model matrix not possible. Please re-estimate",
         " the system with either control variable",
         "  'x' or 'model' set to TRUE" )
   }
   return( result )
}
