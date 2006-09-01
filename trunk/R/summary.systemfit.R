## prepare summary results that belong to the whole system
summary.systemfit <- function(object,...) {
   object$coefTab <- cbind( object$coef, object$se, object$t, object$p )
   colnames( object$coefTab ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   object$rcor <- cor( residuals( object ) )
   dimnames( object$rcor ) <- dimnames( object$rcov )
   class( object ) <- "summary.systemfit"
   return( object )
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
                  round( x$eq[[i]]$df,    digits ),
                  round( x$eq[[i]]$ssr,   digits ),
                  round( x$eq[[i]]$mse,   digits ),
                  round( x$eq[[i]]$rmse,  digits ),
                  round( x$eq[[i]]$r2,    digits ),
                  round( x$eq[[i]]$adjr2, digits ))
    table  <- rbind( table, row )
    if( is.null( x$eq[[i]]$eqnlabel ) ) {
      labels <- rbind( labels, paste( "equation", i ) )
    } else {
      labels <- rbind( labels, x$eq[[i]]$eqnlabel )
    }
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )

  cat("\n")

  if(!is.null(x$rcovest)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    rcov <- x$rcovest
    rownames(rcov) <- labels
    colnames(rcov) <- labels
    print( rcov )
    cat("\n")
    if( min(eigen( x$rcov, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  rcov <- x$rcov
  rownames(rcov) <- labels
  colnames(rcov) <- labels
  print( rcov )
  cat("\n")

  cat("The correlations of the residuals\n")
  rcor <- x$rcor
  rownames(rcor) <- labels
  colnames(rcor) <- labels
  print( rcor )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat(x$drcov)
  cat("\n")

  cat("OLS R-squared value of the system: ")
  cat(x$olsr2)
  cat("\n")

  if(!is.null(x$mcelr2)) {
    cat("McElroy's R-squared value for the system: ")
    cat(x$mcelr2)
    cat("\n")
  }
  ## now print the individual equations
  for(i in 1:x$nEq) {
      print( summary( x$eq[[i]] ), digits )
  }
  invisible( x )
}


## prepare summary results for a single equation
summary.systemfit.equation <- function(object,...) {
   result <- list()
   result$eqnLabel <- object$eqnlabel
   result$eqnNo <- object$i
   result$formula <- object$formula
   result$instruments <- object$inst
   result$method <- object$method
   result$residuals <- object$residuals
   result$coefficients <- cbind( object$b, object$se, object$t, object$p )
   colnames( result$coefficients ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   result$df <- c( object$nExog, object$nObs - object$nExog )
   result$nObs <- object$nObs
   result$sigma <- object$s
   result$ssr <- object$ssr
   result$mse <- object$mse
   result$rmse <- object$rmse
   result$r.squared <- object$r2
   result$adj.r.squared <- object$adjr2
   result$coefCov <- object$bcov
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
            # s ist the variance, isn't it???

  cat( paste( "Number of observations:", round( x$nObs, digits ),
              "Degrees of Freedom:", round( x$df[ 2 ], digits ),"\n" ) )

  cat( paste( "SSR:", round( x$ssr, digits ),
              "MSE:", round( x$mse, digits ),
              "Root MSE:", round(x$rmse, digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round( x$r.squared, digits ),
              "Adjusted R-Squared:", round( x$adj.r.squared, digits ),
              "\n" ) )
  cat("\n")
  invisible( x )
}
