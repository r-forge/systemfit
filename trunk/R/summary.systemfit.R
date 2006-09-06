## prepare summary results that belong to the whole system
summary.systemfit <- function( object, probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   # preparing objects that will be returned
   result <- list()
   result$call <- object$call
   result$nEq <- object$nEq
   result$method <- object$method
   result$iter <- object$iter
   result$maxiter <- object$maxiter
   result$residuals <- residuals( object )
   result$residCovEst <- object$rcovest
   result$residCov <- object$rcov
   rownames( result$residCov ) <- colnames( result$residCov ) <- object$eqnLabels
   if( !is.null( result$residCovEst ) ) {
      dimnames( result$residCovEst ) <- dimnames( result$residCov )
   }
   result$residCor <- cor( residuals( object ) )
   dimnames( result$residCor ) <- dimnames( result$residCov )
   result$detResidCov <- object$drcov

   # coefficients, standard errors, ... 
   result$coefCov <- object$bcov
   coef <- object$coef
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( probDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df ) )
   } else {
      pVal <- rep( NA, length( coef ) )
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( object$nExogAll, object$nObsAll - object$nExogAll )
   result$ols.r.squared <- object$olsr2
   result$mcelroy.r.squared <- object$mcelr2

   # now prepare summury results for the individual equations
   result$eq <- list()
   for( i in 1:object$nEq ) {
       result$eq[[ i ]] <- summary( object$eq[[i]], probDfSys = probDfSys )
   }

   class( result ) <- "summary.systemfit"
   return( result )
}

## print summary results of the whole system
print.summary.systemfit <- function( x, digits=6,... ) {

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("systemfit results \n")
  cat("method: ")
  if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
  cat( paste( x$method, "\n\n"))
  if(!is.null(x$iter)) {
    if(x$iter>1) {
      if(x$iter<x$maxiter) {
        cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
      } else {
        cat( paste( "warning: convergence not achieved after", x$iter,
                    "iterations\n\n" ) )
      }
    }
  }
  for(i in 1:x$nEq) {
    row <- NULL
    row <- cbind( round( x$eq[[i]]$nObs,  digits ),
                  round( x$eq[[i]]$df[2], digits ),
                  round( x$eq[[i]]$ssr,   digits ),
                  round( x$eq[[i]]$sigma^2, digits ),
                  round( x$eq[[i]]$sigma,   digits ),
                  round( x$eq[[i]]$r.squared,     digits ),
                  round( x$eq[[i]]$adj.r.squared, digits ))
    table  <- rbind( table, row )
    if( is.null( x$eq[[i]]$eqnLabel ) ) {
      labels <- rbind( labels, paste( "equation", i ) )
    } else {
      labels <- rbind( labels, x$eq[[i]]$eqnLabel )
    }
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )

  cat("\n")

  if(!is.null(x$residCovEst)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    print( x$residCovEst )
    cat("\n")
    if( min(eigen( x$residCov, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  print( x$residCov )
  cat("\n")

  cat("The correlations of the residuals\n")
  print( x$residCor )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat( x$detResidCov )
  cat("\n")

  cat("OLS R-squared value of the system: ")
  cat( x$ols.r.squared )
  cat("\n")

  if(!is.null(x$mcelroy.r.squared)) {
    cat("McElroy's R-squared value for the system: ")
    cat(x$mcelroy.r.squared)
    cat("\n")
  }
  ## now print the individual equations
  for(i in 1:x$nEq) {
      print( x$eq[[i]], digits )
  }
  invisible( x )
}


## prepare summary results for a single equation
summary.systemfit.equation <- function( object, probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   # preparing objects that will be returned
   result <- list()
   result$eqnLabel <- object$eqnlabel
   result$eqnNo <- object$i
   result$formula <- object$formula
   result$instruments <- object$inst
   result$method <- object$method
   result$residuals <- object$residuals

   # coefficients, standard errors, ... 
   result$coefCov <- object$bcov
   coef <- object$coef
   stdEr <- diag( result$coefCov )^0.5  # standard errors
   tStat <- coef / stdEr                # t-statistic
   if( probDfSys ) {             # p-values
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$dfSys ) )
   } else {
      pVal <- 2 * ( 1 - pt( abs( tStat ), object$df ) )
   }
   result$coefficients <- cbind( coef, stdEr, tStat, pVal )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( object$nExog, object$nObs - object$nExog )
   result$nObs <- object$nObs
   result$sigma <- object$sigma
   result$ssr <- object$ssr
   result$r.squared <- object$r2
   result$adj.r.squared <- object$adjr2
   class( result ) <- "summary.systemfit.equation"
   return( result )
}


## print summary results for a single equation
print.summary.systemfit.equation <- function( x, digits=6, ... ) {

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  cat("\n")
  if( is.null( x$eqnLabel ) ) {
    cat( x$method, " estimates for equation ", x$eqnNo, "\n", sep = "" )
  } else {
    cat( x$method, " estimates for '", x$eqnLabel,
         "' (equation ", x$eqnNo, ")\n", sep = "" )
  }

  cat("Model Formula: ")
  print(x$formula)
  if(!is.null(x$inst)) {
    cat("Instruments: ")
    print(x$inst)
  }
  cat("\n")

  Signif <- symnum(x$coefficients[,4], corr = FALSE, na = FALSE,
                   cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
                   symbols   = c( "***", "**", "*", "." ," " ))

  table <- cbind(round( x$coefficients,  digits ),
                 Signif)

  rownames(table) <- rownames(x$coefficients)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round( x$sigma, digits ),
            "on", x$df[ 2 ], "degrees of freedom\n" ))

  cat( paste( "Number of observations:", round( x$nObs, digits ),
              "Degrees of Freedom:", round( x$df[ 2 ], digits ),"\n" ) )

  cat( paste( "SSR:", round( x$ssr, digits ),
              "MSE:", round( x$sigma^2, digits ),
              "Root MSE:", round(x$sigma, digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round( x$r.squared, digits ),
              "Adjusted R-Squared:", round( x$adj.r.squared, digits ),
              "\n" ) )
  cat("\n")
  invisible( x )
}
