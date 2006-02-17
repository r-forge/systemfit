## Likelihood Ratio Test
lrtest.systemfit <- function( resultc, resultu ) {
  lrtest <- list()
  if( resultc$method %in% c( "SUR", "WSUR" ) &
      resultu$method %in% c( "SUR", "WSUR" ) ) {
    n   <- resultu$eq[[1]]$n
    lrtest$df  <- resultu$ki - resultc$ki
    if(resultc$rcovformula != resultu$rcovformula) {
      stop( paste( "both estimations must use the same formula to calculate",
                   "the residual covariance matrix!" ) )
    }
    if(resultc$rcovformula == 0) {
      lrtest$lr  <- n * ( log( resultc$drcov ) - log( resultu$drcov ) )
    } else {
      residc <- array(resultc$resids,c(n,resultc$g))
      residu <- array(resultu$resids,c(n,resultu$g))
      lrtest$lr <- n * ( log( det( (t(residc) %*% residc)) ) -
                         log( det( (t(residu) %*% residu))))
    }
    lrtest$p <- 1-pchisq( lrtest$lr, lrtest$df )
  }
  lrtest
}
