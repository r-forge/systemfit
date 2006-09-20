## Likelihood Ratio Test
lrtest.systemfit <- function( resultc, resultu ) {
  lrtest <- list()
  if( resultc$method %in% c( "SUR", "WSUR" ) &
      resultu$method %in% c( "SUR", "WSUR" ) ) {
    nObs <- resultu$nObsAll / resultu$nEq
    lrtest$nRestr  <- resultu$nExogLiAll - resultc$nExogLiAll
    if( resultc$control$methodRCov != resultu$control$methodRCov ) {
      stop( paste( "both estimations must use the same formula to calculate",
                   "the residual covariance matrix!" ) )
    }
    if( resultc$control$methodRCov == 0 ) {
      lrtest$statistic  <- nObs * ( log( resultc$drcov ) - log( resultu$drcov ) )
    } else {
      residc <- as.matrix( residuals( resultc ) )
      residu <- as.matrix( residuals( resultu ) )
      lrtest$statistic <- nObs * ( log( det( (t(residc) %*% residc)) ) -
                         log( det( (t(residu) %*% residu))))
    }
    lrtest$p.value <- 1 - pchisq( lrtest$statistic, lrtest$nRestr )
  }
  class( lrtest ) <- "lrtest.systemfit"
  return( lrtest )
}

print.lrtest.systemfit <- function( x, digits = 4, ... ){
   cat( "\n", "Likelihood-Ratio-test for parameter restrictions",
      " in equation systems\n", sep = "" )
   cat( "LR-statistic:", formatC( x$statistic, digits = digits ), "\n" )
   cat( "degrees of freedom:", x$nRestr, "\n" )
   cat( "p-value:", formatC( x$p.value, digits = digits ), "\n\n" )
   invisible( x )
}

