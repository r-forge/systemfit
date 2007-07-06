.stackMatList <- function( matList, way, useMatrix = FALSE ){
   if( way == "diag" ){
      result <- matrix( 0, 0, 0 )
      for( i in 1:length( matList ) ){
         result <- rbind( 
            cbind( result, 
               matrix( 0, nrow( result ), ncol( matList[[ i ]] ) ) ),
            cbind( matrix( 0, nrow( matList[[ i ]] ), ncol( result ) ),
               matList[[ i ]] ) )
      }
   } else if( way == "below" ) {
      result <- NULL
      for( i in 1:length( matList ) ){
         result <- rbind( result, matList[[ i ]] )
      }
   }

   if( useMatrix ){
      result <- as( result, "dgCMatrix" )
   }

   return( result )
}

.prepareWmatrix <- function( upperleft, R.restr ){
   result <- rbind2( cbind2( upperleft, t(R.restr) ),
      cbind2( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
   return( result )
}
