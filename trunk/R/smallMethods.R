## this function returns a vector of the
## cross-equation corrlations between eq i and eq j
## from the results set for equation ij
correlation.systemfit <- function( results, eqni, eqnj ) {
  nExogEq <- results$nExogEq
  cij <- results$bcov[(1+sum(nExogEq[1:eqni])-nExogEq[eqni]):(sum(nExogEq[1:eqni])),
                      (1+sum(nExogEq[1:eqnj])-nExogEq[eqnj]):(sum(nExogEq[1:eqnj]))]
  cii <- results$bcov[(1+sum(nExogEq[1:eqni])-nExogEq[eqni]):(sum(nExogEq[1:eqni])),
                      (1+sum(nExogEq[1:eqni])-nExogEq[eqni]):(sum(nExogEq[1:eqni]))]
  cjj <- results$bcov[(1+sum(nExogEq[1:eqnj])-nExogEq[eqnj]):(sum(nExogEq[1:eqnj])),
                      (1+sum(nExogEq[1:eqnj])-nExogEq[eqnj]):(sum(nExogEq[1:eqnj]))]
  rij <- NULL

  for(i in 1:results$eq[[1]]$nObs ) {
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
  for(i in 1:resultsi$eq[[eqni]]$nObs ) {
    xik    <- resultsi$eq[[eqni]]$modelMatrix[i,]
    top    <- sqrt( xik %*% resultsi$eq[[eqni]]$bcov %*% xik )
    bottom <- sqrt( xik %*% resultsj$eq[[eqni]]$bcov %*% xik )
    rk     <- top / bottom
    ratio  <- rbind( ratio, rk )
  }
  ratio
}


## return all coefficients
coef.systemfit <- function( object, ... ) {
   object$coef
}

## return all coefficients, std.errors, t-values and p-values
coef.summary.systemfit <- function( object, ... ) {
   object$coefTab
}

## return the coefficients of a single equation
coef.systemfit.equation <- function( object, ... ) {
   object$coef
}

## return coefficients, std.errors, t-values and p-values of a single equation
coef.summary.systemfit.equation <- function( object, ... ) {
   object$coef
}

## return all residuals
residuals.systemfit <- function( object, ... ) {
   if( is.null( colnames( object$rcov ) ) ) {
      eqNames <- paste( "eq", c( 1:object$nEq ) )
   } else {
      eqNames <- colnames( object$rcov )
   }
   residuals <- data.frame( obsNo = c( 1:length( object$eq[[1]]$residuals ) ) )
   for( i in 1:object$nEq ) {
      residuals[[ eqNames[ i ] ]] <- object$eq[[i]]$residuals
   }
   residuals$obsNo <- NULL
   return( residuals )
}

## return residuals of a single equation
residuals.systemfit.equation <- function( object, ... ) {
   object$residuals
}

## return the variance covariance matrix of the coefficients
vcov.systemfit <- function( object, ... ) {
   object$bcov
}

## return the variance covariance matrix of the coefficients of a single equation
vcov.systemfit.equation <- function( object, ... ) {
   object$bcov
}


## return the fitted values
fitted.systemfit <- function( object, ... ) {
   fitted.values <- matrix( NA, length( object$eq[[1]]$fitted.values ), object$nEq )
   colnames( fitted.values ) <- as.character( 1:ncol( fitted.values ) )
   for(i in 1:object$nEq )  {
      fitted.values[ , i ]           <- object$eq[[ i ]]$fitted.values
      colnames( fitted.values )[ i ] <- paste( "eq", as.character(i), sep="" )
   }
   fitted.values
}

## return the fitted values of e single euation
fitted.systemfit.equation <- function( object, ... ) {
   object$fitted.values
}




