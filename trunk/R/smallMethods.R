## this function returns a vector of the
## cross-equation corrlations between eq i and eq j
## from the results set for equation ij
correlation.systemfit <- function( results, eqni, eqnj ) {
  nCoefEq <- NULL
  for( i in 1:length( results$eq ) ) {
     nCoefEq <- c( nCoefEq, length( coef( results$eq[[ i ]] ) ) )
  }
  cij <- results$bcov[(1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni])),
                      (1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj]))]
  cii <- results$bcov[(1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni])),
                      (1+sum(nCoefEq[1:eqni])-nCoefEq[eqni]):(sum(nCoefEq[1:eqni]))]
  cjj <- results$bcov[(1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj])),
                      (1+sum(nCoefEq[1:eqnj])-nCoefEq[eqnj]):(sum(nCoefEq[1:eqnj]))]
  rij <- NULL

  for( i in 1:nrow( residuals( results ) ) ) {
    xik    <- results$eq[[eqni]]$modelMatrix[i,]
    xjk    <- results$eq[[eqnj]]$modelMatrix[i,]
    top    <- xik %*% cij %*% xjk
    bottom <- sqrt( ( xik %*% cii %*% xik ) * ( xjk %*% cjj %*% xjk ) )
    rijk   <- top / bottom
    rij    <- rbind( rij, rijk )
  }
  rij
}

## determines the improvement of resultsj (3sls) over
## resultsi (2sls) for equation i and returns a matrix
## of the values, so you can examine the range, mean, etc
se.ratio.systemfit <- function( resultsi, resultsj, eqni ) {
  ratio <- NULL
  for( i in 1:nrow( residuals( resultsi ) ) ) {
    xik    <- resultsi$eq[[eqni]]$modelMatrix[i,]
    top    <- sqrt( xik %*% resultsi$eq[[eqni]]$bcov %*% xik )
    bottom <- sqrt( xik %*% resultsj$eq[[eqni]]$bcov %*% xik )
    rk     <- top / bottom
    ratio  <- rbind( ratio, rk )
  }
  ratio
}


## return all coefficients
coef.systemfit <- function( object, transformed = FALSE, ... ) {
   if( transformed ){
      if( is.null( object$restrict.reg ) ){
         stop( "transformed coefficients are not available,",
            " because argument 'restrict.reg' has not been used in this estimation." )
      } else {
         return( drop( solve( crossprod( object$restrict.reg ),
            t( object$restrict.reg ) %*% coef( object ) ) ) )
      }
   } else {
      return( object$coefficients )
   }
}

## return all coefficients, std.errors, t-values and p-values
coef.summary.systemfit <- function( object, transformed = FALSE, ... ) {
   if( transformed ){
      if( is.null( object$coefTrans ) ){
         stop( "transformed coefficients are not available,",
            " because argument 'restrict.reg' has not been used in this estimation." )
      } else {
         return( object$coefTrans )
      }
   } else {
      return( object$coefficients )
   }
}

## return the coefficients of a single equation
coef.systemfit.equation <- function( object, ... ) {
   object$coefficients
}

## return coefficients, std.errors, t-values and p-values of a single equation
coef.summary.systemfit.equation <- function( object, ... ) {
   object$coefficients
}

## return all residuals
residuals.systemfit <- function( object, ... ) {
   result <- data.frame( obsNo = c( 1:length( residuals( object$eq[[1]] ) ) ) )
   for( i in 1:length( object$eq ) ) {
      result[[ object$eq[[i]]$eqnLabel ]] <- residuals( object$eq[[i]] )
   }
   result$obsNo <- NULL
   return( result )
}

## return residuals of a single equation
residuals.systemfit.equation <- function( object, ... ) {
   object$residuals
}

## return the variance covariance matrix of the coefficients
vcov.systemfit <- function( object, transformed = FALSE, ... ) {
   if( transformed ){
      if( is.null( object$restrict.reg ) ){
         stop( "transformed coefficients and their covariance matrix",
            " are not available,",
            " because argument 'restrict.reg' has not been used in this estimation." )
      } else {
         txtxInv <- solve( crossprod( object$restrict.reg ) )
         result <- txtxInv %*% t( object$restrict.reg ) %*% vcov( object ) %*%
            object$restrict.reg %*% txtxInv
         return( result )
      }
   } else {
      return( object$bcov )
   }
}

## return the variance covariance matrix of the coefficients of a single equation
vcov.systemfit.equation <- function( object, ... ) {
   object$bcov
}


## return the fitted values
fitted.systemfit <- function( object, ... ) {
   nEq <- length( object$eq )
   fitted.values <- matrix( NA, length( object$eq[[1]]$fitted.values ), nEq )
   colnames( fitted.values ) <- as.character( 1:ncol( fitted.values ) )
   for(i in 1:nEq )  {
      fitted.values[ , i ]           <- object$eq[[ i ]]$fitted.values
      colnames( fitted.values )[ i ] <- object$eq[[ i ]]$eqnLabel
   }
   fitted.values
}

## return the fitted values of e single euation
fitted.systemfit.equation <- function( object, ... ) {
   object$fitted.values
}

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
   if( !is.null( object$modelMatrix ) ) {
      result <- object$modelMatrix
   } else if( !is.null( model.frame( object ) ) ) {
      result <- model.matrix( object$terms, data = model.frame( object ) )
   } else {
      stop( "returning model matrix not possible. Please re-estimate",
         " the system with either control variable",
         "  'returnModelMatrix' or 'returnModelFrame' set to TRUE" )
   }
   return( result )
}

## return model frame of the entire system
model.frame.systemfit <- function( formula, ... ){
   mfColNames <- NULL
   for( i in 1:length( formula$eq ) ) {
      mfi <- model.frame( formula$eq[[ i ]] )
      if( i == 1 ) {
         result <- mfi
      } else {
         result <- cbind( result, mfi[ , ! names( mfi ) %in% names( result ) ] )
      }
   }
   return( result )
}

## return model frame of a single equation
model.frame.systemfit.equation <- function( formula, ... ){
   if( !is.null( formula$modelFrame ) ) {
      result <- formula$modelFrame
   } else {
      stop( "returning model frame not possible. Please re-estimate",
         " the system with control variable 'returnModelFrame'",
         " set to TRUE" )
   }
   return( result )
}
