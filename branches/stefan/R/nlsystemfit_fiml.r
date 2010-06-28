#man müsste da jetzt noch rausfinden, ob das mit den 
#tests passt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



knls1 <- function( theta, eqns, data, fitmethod="OLS", parmnames,gra=NULL ) {

 r  <- matrix()   
 r <- NULL
  
 residi  <- list()          # residuals equation wise
 lhs <- list()
 rhs <- list()
 neqs <- length( eqns )
 nobs <- dim( data )[[1]]    
 npa<-length(theta)        

 lhs <- NULL
 rhs <- NULL
 residi <- NULL
 # partial derivatives of the residuals with respect to the parameters
 dResidTheta <- NULL
 dResidThetai <- list()
   
 ## get the values of the parameters
 for( i in 1:length( parmnames ) ) {
   name <- names( parmnames )[i]
   val <- theta[i]
   storage.mode( val ) <-  "double"
   assign( name, val )
  }

 ## build the residual vector...
 for( i in 1:length( eqns ) ) {
    lhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[2]],data ) )
    rhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[3]] ,data ) )
    residi[[i]] <- lhs[[i]] - rhs[[i]]
    r <- cbind( r, as.matrix( residi[[i]] ) )

 
  ## implement analytic gradient
  #without analytic gradient
  if(is.null(gra)){
    dResidThetai[[ i ]] <- - attributes( with( data, with( as.list( theta ),
    eval( deriv( eqns[[i]], names( parmnames ))))))$gradient}else{

  ##with analytic gradient
  #analytic gradient:
  #gra must be a list of a lists 
  #1st dimension parameters
  #2nd dimension number of equations
  ##################
    dResidThetai<-list()
    dResidThetama<-matrix(NA,dim(data)[1],length(parmnames))
    for(i2 in 1: length( parmnames )){
      dResidThetama[,i2]<- -  with( data, with( as.list( parmnames ),
      eval( gra[[i2]][[i]][[2]])))}
  
    dResidThetai[[ i ]]<-dResidThetama}

    
 #erste ableitung ist dann ne matrix mit anzahl der parameter spalten
 dResidTheta <- rbind( dResidTheta, dResidThetai[[ i ]] )
 }
 
 
  ## these are the objective functions for the various fitting methods
  if( fitmethod == "FIML" ) {
    ma<-matrix(NA,neqs,neqs)
       for (k in 1:neqs){
       for (l in 1:neqs){
       ma[k,l]<-1/nobs*(r[,k]%*%r[,l])
       }}

       cma <- chol(ma )
       ma_m1<-chol2inv(cma)
   obj <- nobs/2*log(det(ma))
  
  ## calculate the gradients
 
   gradient<-theta
   gradient[1:npa]<-NA
    for (m in 1:npa){
     gra_su<-rep(NA,neqs*neqs)
     kl1<-0
     for(k in 1:neqs){
      for(l in 1:neqs){
       kl1<-kl1+1
       ein<-((l-1)*nobs+1):(l*nobs)
       zwe<-((k-1)*nobs+1):(k*nobs)
       gra_su[kl1]<-ma_m1[l,k]*(sum(dResidTheta[ein,m]*r[,k]+dResidTheta[zwe,m]*r[,l]))  
       }}
  gradient[m]<-1/2*sum(gra_su)}


  attributes( obj ) <- list( gradient = gradient )
  return(obj)
  }
 }



nlsystemfit_fiml <- function( method="FIML",
                        eqns,
                        startvals,
                        eqnlabels=c(as.character(1:length(eqns))),
                        inst=NULL,
                        data=list(),gra=NULL,
                        solvtol=.Machine$double.eps,
                        maxiter=1000, par_tol=1e-6, ... ) {# par_tol is the 
#stopping tolerance for ITSUR


  attach( data )

  ## some tests
  if(!(method=="OLS" | method=="SUR" | method=="2SLS" | method=="3SLS" | method=="GMM" |
method=="FIML" | method=="ITSUR")){
    stop("The method must be 'OLS', 'SUR', '2SLS', '3SLS', 'FIML' or 'ITSUR'")}
  if((method=="2SLS" | method=="3SLS" | method=="GMM") & is.null(inst)) {
    stop("The methods '2SLS', '3SLS' and 'GMM' need instruments!")}

  lhs <- list()
  rhs <- list()
  derivs <- list()

  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  resulti <- list()               # results of the ith equation
  residi  <- list()               # residuals equation wise
  iter    <- NULL                 # number of iterations
  G       <- length( eqns )       # number of equations
  n       <- array( 0, c(G))      # number of observations in each equation
  k       <- array( 0, c(G) )     # number of (unrestricted) coefficients/regressors in each equation
  df       <- array( 0, c(G) )     # degrees of freedom in each equation
  instl   <- list()               # list of the instruments for each equation
  ssr     <- array( 0, c(G))      # sum of squared residuals of each equation
  mse     <- array( 0, c(G))      # mean square error (residuals) of each equation
  rmse    <- array( 0, c(G))      # root of mse
  r2      <- array( 0, c(G))      # R-squared value
  adjr2   <- array( 0, c(G))      # adjusted R-squared value
  nobs <- dim( data )[[1]]
  S <- matrix( 0, G, G )               # covariance matrix of the residuals
  X <- array()
  x <- list()

  resids <- array()
  resids <- NULL

  if( method == "FIML" ) {
    est <- nlm( knls1, startvals,
       typsize=abs(startvals),iterlim=1000,eqns,
       data,gradtol = 1e-6,
       stepmax=20,
       method,parmnames=startvals )
     }
  
  ## done with the fitting...
  ## now, part out the results from the nlm function
  ## to rebuild the equations and return object
  ## get the parameters for each of the equations and


  ## evaluate the residuals for eqn
  ## get the values of the final parameters
  if( TRUE ) {
    estimate <- est$estimate
  } else {
    estimate <- est$par
  }
  names( estimate ) <- names( startvals )
  for( i in 1:length( estimate ) ) {
    name <- names( estimate )[i]
    ### I wonder if I need to clear out the variables before assigning them for good measure...
    assign( name, NULL )
    val <- estimate[i]
    storage.mode( val ) <-  "double"
    assign( name, val )
  }

  ## get the rank for the eqns, compute the first-stage
  ## cov matrix to finish the SUR and 3SLS methods
  X <- NULL
  results$resids <- array()
  results$resids <- NULL

  ## you're working on parsing out the parameters and the estiamtes for the return structure...
  for(i in 1:G) {
    lhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[2]] ) )
    rhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[3]] ) )
    residi[[i]] <- lhs[[i]] - rhs[[i]]
    derivs[[i]] <- deriv( as.formula( eqns[[i]] ), names( startvals ) )

    ## computing the jacobian to get the rank to get the number of variables...
    jacobian <- attr( eval( derivs[[i]] ), "gradient" )
    n[i]   <-  length( lhs[[i]] )
    k[i] <- qr( jacobian )$rank
    df[i] <- n[i] - k[i]
    ssr[i] <- crossprod( residi[[i]] )
    mse[i] <- ssr[i] / ( n[i] - k[i] )
    rmse[i] <- sqrt( mse[i] )

    X <- rbind( X, jacobian )
    results$resids <- cbind( results$resids, as.matrix( residi[[i]] ) )
  }

  ## compute the final covariance matrix
  ## you really should use the code below to handle weights...
  rcovformula <- 1
  for(i in 1:G) {
    for(j in 1:G) {
      S[i,j] <- sum(residi[[i]]*residi[[j]])/(
                                              sqrt((n[i]-rcovformula*k[i])*(n[j]-rcovformula*k[j])))
    }
  }

### for when you get the weights working...
#     vardef <- 1
#     if( vardef == 1 ) {
#       D <- diag( G ) * 1 / sqrt( nrow( data ) )
#     }
#     if( vardef == 2 ) {
#       D <- diag( G ) * 1 / sqrt( sum( weights ) )
#     }
#     if( vardef == 3 ) {
#       D <- diag( G ) * 1 / sqrt( sum( weights ) - ( sum( n ) - sum( k ) ) )
#     }
#     if( vardef == 4 ) {
#       for(i in 1:G) {
#         D <- diag( G )
#         D[i,i] <- D[i,i] * 1 / sqrt( nrow( data ) - n[i] )
#       }
#     }
#     ## the docs have this, but the table contains the above equations
#     R <- crossprod( results$resids )
#     S <- D %*% R %*% D
#     SI <- qr.solve( S, tol=solvtol ) %x% diag( nrow( data ) )
#     covb <- qr.solve(t(X) %*% SI %*% X, tol=solvtol )



  ## get the variance-covariance matrix
  if( method == "FIML" ) {
    SI <- diag( diag( qr.solve( S, tol=solvtol ) ) ) %x% diag( nrow( data ) )
    covb <- qr.solve(t(X) %*% SI %*% X, tol=solvtol )
  }

  colnames( covb ) <- rownames( covb )


  ## bind the standard errors to the parameter estimate matrix
  se2 <- sqrt( diag( covb ) )
  t.val <- estimate / se2
  prob  <- 2*( 1 - pt( abs( t.val ), sum( n ) - sum( k ) ) ) ### you better check this...

  results$method       <- method
  results$n <- sum( n )
  results$k <- sum( k )
  results$b <- estimate
  results$se <- se2
  results$t <- t.val
  results$p <- prob

  ## build the results structure...
  for(i in 1:G) {
    ## you may be able to shrink this up a little and write the values directly to the return structure...
    eqn.terms <- vector()
    eqn.est <- vector()
    eqn.se <- vector()
    jacob <- attr( eval( deriv( as.formula( eqns[[i]] ), names( startvals ) ) ), "gradient" )
    for( v in 1:length( estimate ) ) {
      if( qr( jacob[,v] )$rank > 0 ) {
        eqn.terms <- rbind( eqn.terms, name <- names( estimate )[v] )
        eqn.est <- rbind( eqn.est, estimate[v] )
        eqn.se <- rbind( eqn.se, se2[v] )
      }
    }


    ## build the "return" structure for the equations
    resulti$method       <- method
    resulti$i            <- i               # equation number
    resulti$eqnlabel     <- eqnlabels[[i]]
    resulti$formula      <- eqns[[i]]
    resulti$b <- as.vector( eqn.est )
    names( resulti$b )   <- eqn.terms
    resulti$se           <- eqn.se
    resulti$t            <- resulti$b / resulti$se
    resulti$p            <- 2*( 1-pt(abs(resulti$t), n[i] - k[i] ))
    resulti$n            <- n[i]            # number of observations
    resulti$k            <- k[i]            # number of coefficients/regressors
    resulti$df           <- df[i]           # degrees of freedom of residuals
    resulti$predicted    <- rhs[[i]]           # predicted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr[i]             # sum of squared errors/residuals
    resulti$mse          <- mse[i]             # estimated variance of the residuals (mean squared error)
    resulti$s2           <- mse[i]             #        the same (sigma hat squared)
    resulti$rmse         <- rmse[i]            # estimated standard error of the residuals
    resulti$s            <- rmse[i]            #        the same (sigma hat)

#     ## you'll need these to compute the correlations...
#     print( paste( "eqn ", i ) )
    coefNames <- rownames( covb )[ rownames( covb ) %in%
      strsplit( as.character( eqns[[ i ]] )[ 3 ], "[^a-zA-Z0-9.]" )[[ 1 ]] ]
    resulti$covb <- covb[ coefNames, coefNames ]

#     resulti$x <- model.frame( as.formula( eqns[[i]] )[[3]] )
#     print( resulti$x )
#    print( model.frame( eval( eqns[[i]] ) ) )



    ## fix this to allow for multiple instruments?
    if(method=="2SLS" | method=="3SLS" | method=="GMM") {
      resulti$inst         <- inst
      ##resulti$inst         <- inst[[i]]
      ##resulti$inst         <- instl[[i]]
      ## resulti$h            <- h[[i]]          # matrix of instrumental variables
    }

    resulti$r2     <- 1 - ssr[i] / ( ( crossprod( lhs[[i]]) ) - mean( lhs[[i]] )^2 * nobs )
    resulti$adjr2  <- 1 - ((n[i]-1)/df[i])*(1-resulti$r2)

    class(resulti)        <- "nlsystemfit.equation"
    results$eq[[i]]      <- resulti
  }

  results$solvtol <- solvtol
  results$covb <- covb
  results$rcov <- S
  results$rcor <- cor( results$resids )
  results$drcov <- det(results$rcov)          # det(rcov, tol=solvetol)

  if(method=="2SLS" || method=="3SLS") {
    ##      results$h       <- H            # matrix of all (diagonally stacked) instrumental variables
  }
  if(method=="SUR" || method=="3SLS" || method=="GMM" || method=="ITSUR") {
    results$rcovest <- Sols      # residual covarance matrix used for estimation
    ##results$mcelr2  <- mcelr2       # McElroy's R-squared value for the equation system
  }

  ## build the "return" structure for the whole system
  results$method  <- method
  results$g       <- G              # number of equations
  results$nlmest  <- est

  class(results)  <- "nlsystemfit.system"

  detach(data)

  if( results$nlmest$code >= 4 ) {
    warning( "Estimation did not converge!" )
  }

  results
}


## print the (summary) results that belong to the whole system
summary.nlsystemfit.system <- function(object,...) {
  summary.nlsystemfit.system <- object
  summary.nlsystemfit.system
}


## print the results that belong to the whole system
print.nlsystemfit.system <- function( x, digits=6,... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("nlsystemfit results \n")
  cat("method: ")

#  if(!is.null(object$iter)) if(object$iter>1) cat("iterated ")
  cat( paste( object$method, "\n\n"))
#   if(!is.null(object$iter)) {
#     if(object$iter>1) {
#       if(object$iter<object$maxiter) {
#         cat( paste( "convergence achieved after",object$iter,"iterations\n\n" ) )
#       } else {
#         cat( paste( "warning: convergence not achieved after",object$iter,"iterations\n\n" ) )
#       }
#     }
#   }

  cat( paste( "convergence achieved after",object$nlmest$iterations,"iterations\n" ) )
  cat( paste( "nlsystemfit objective function value:",object$nlmest$minimum,"\n\n" ) )


  for(i in 1:object$g) {
    row <- NULL
    row <- cbind( round( object$eq[[i]]$n,     digits ),
                  round( object$eq[[i]]$df,    digits ),
                  round( object$eq[[i]]$ssr,   digits ),
                  round( object$eq[[i]]$mse,   digits ),
                  round( object$eq[[i]]$rmse,  digits ),
                  round( object$eq[[i]]$r2,    digits ),
                  round( object$eq[[i]]$adjr2, digits ))
    table  <- rbind( table, row )
    labels <- rbind( labels, object$eq[[i]]$eqnlabel )
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("\n")

  ## check this code before release...
  if(!is.null(object$rcovest)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    rcov <- object$rcovest
    rownames(rcov) <- labels
    colnames(rcov) <- labels
    print( rcov )
    cat("\n")
    if( min(eigen( object$rcovest, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  rcov <- object$rcov
  rownames(rcov) <- labels
  colnames(rcov) <- labels
  print( rcov )
  cat("\n")

  cat("The correlations of the residuals\n")
  rcor <- object$rcor
  rownames(rcor) <- labels
  colnames(rcor) <- labels
  print( rcor )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat(object$drcov)
  cat("\n")

### check this code before release
#   cat("OLS R-squared value of the system: ")
#   cat(object$olsr2)
#   cat("\n")

#   if(!is.null(object$mcelr2)) {
#     cat("McElroy's R-squared value for the system: ")
#     cat(object$mcelr2)
#     cat("\n")
#   }

  ## now print the individual equations
  for(i in 1:object$g) {
    print( object$eq[[i]], digits )
  }

}


## print the (summary) results for a single equation
summary.nlsystemfit.equation <- function(object,...) {
  summary.nlsystemfit.equation <- object
  summary.nlsystemfit.equation
}


## print the results for a single equation
print.nlsystemfit.equation <- function( x, digits=6, ... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  cat("\n")
  cat( paste( object$method, " estimates for ", object$eqnlabel, " (equation ", object$i, ")\n", sep="" ) )

  cat("Model Formula: ")
  print(object$formula)
  if(!is.null(object$inst)) {
    cat("Instruments: ")
    print(object$inst)
  }
  cat("\n")

  Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols   = c("***","**","*","."," "))

  table <- cbind(round( object$b,  digits ),
                 round( object$se, digits ),
                 round( object$t,  digits ),
                 round( object$p,  digits ),
                 Signif)

  rownames(table) <- names(object$b)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round(object$s, digits),  ## s ist the variance, isn't it???
            "on", object$df, "degrees of freedom\n"))

  cat( paste( "Number of observations:", round(object$n, digits),
              "Degrees of Freedom:", round(object$df, digits),"\n" ) )

  cat( paste( "SSR:",     round(object$ssr,    digits),
              "MSE:", round(object$mse, digits),
              "Root MSE:",   round(object$rmse,  digits), "\n" ) )

   cat( paste( "Multiple R-Squared:", round(object$r2,    digits),
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
  cat("\n")
}

# from Model Selection and Inference: A Practical Information-Theoretic Approach
# Kenneth P. Burnham and David R. Anderson, 1998. Springer-Verlag, New York, New York.

## Akaike's Information Criterion
## AIC = n * log( sigmahat^2 ) + 2K
## n = number of obs
## sigmahat^2 = sum( error^2 ) / n == residual sums of squares
## K is the total number if estimated parameters, including the intercept and sigma^2 (nparams + 1)
## second order AIC
## AICc = AIC + (2K*(K+1))/(n-K-1)




