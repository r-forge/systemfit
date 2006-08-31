###   $Id$
###
###            Simultaneous Equation Estimation for R
###
### Copyright 2002-2004 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###                     Arne Henningsen <http://www.arne-henningsen.de>
###
### This file is part of the nlsystemfit library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


systemfit <- function(  eqns,
                        method = "OLS",
                        inst=NULL,
                        data=list(),
                        R.restr=NULL,
                        q.restr=matrix(0,max(nrow(R.restr),0),1),
                        TX=NULL,
                        maxiter=1,
                        tol=1e-5,
                        methodRCov="geomean",
                        centerResiduals = FALSE,
                        method3sls="GLS",
                        probdfsys=!(is.null(R.restr) & is.null(TX)),
                        single.eq.sigma=(is.null(R.restr) & is.null(TX)),
                        solvetol=.Machine$double.eps,
                        saveMemory = ( nrow( data ) * length( eqns ) > 1000 &&
                           length( data ) > 0 ) )
{

   ## some tests
   if(!( method=="OLS" | method=="WLS" | method=="SUR" | method=="WSUR" |
         method=="2SLS" | method=="W2SLS" | method=="3SLS" | method=="W3SLS" |
         method=="LIML" | method=="FIML")){
      stop( "The method must be 'OLS', 'WLS', 'SUR', 'WSUR',",
         " '2SLS', 'W2SLS', '3SLS', or 'W3SLS'" )
   }
   if( ( method=="2SLS" | method=="W2SLS" | method=="3SLS" | method=="W3SLS" ) &
         is.null(inst) ) {
      stop( "The methods '2SLS', 'W2SLS', '3SLS', and 'W3SLS' need instruments!" )
   }

  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  resulti <- list()               # results of the ith equation
  residi  <- list()               # residuals equation wise
  iter    <- NULL                 # number of iterations
  nEq     <- length( eqns )       # number of equations
  yVecEq  <- list()               # list for vectors of endogenous variables in each equation
  yVecAll <- matrix( 0, 0, 1 )    # stacked endogenous variables of all equations
  xMatEq  <- list()               # list for matrices of regressors in each equation
  xMatAll <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
  nObsEq  <- array( 0, c(nEq))    # number of observations in each equation
  nExogEq <- array( 0, c(nEq) )   # number of exogenous variables /(unrestricted) coefficients
                                     # in each equation
  instl   <- list()               # list of the instruments for each equation
  ssr     <- array( 0, c(nEq))    # sum of squared residuals of each equation
  mse     <- array( 0, c(nEq))    # mean square error (residuals) of each equation
  rmse    <- array( 0, c(nEq))    # root of mse
  r2      <- array( 0, c(nEq))    # R-squared value
  adjr2   <- array( 0, c(nEq))    # adjusted R-squared value
  xnames  <- NULL                 # names of regressors

   if( is.null( names( eqns ) ) ) {
      eqnlabels <- paste( "eq", c( 1:nEq ), sep = "" )
   } else {
      eqnlabels <- names(eqns)
   }

#   for(i in 1:nEq )  {
#     yVecEq[[i]] <-  eval( attr( terms( eqns[[i]] ), "variables" )[[2]] )
#     yVecAll      <-  c( yVecAll, yVecEq[[i]] )
#     xMatEq[[i]] <-  model.matrix( eqns[[i]] )
#     xMatAll      <-  rbind( cbind( xMatAll, matrix( 0, nrow( xMatAll ), ncol( xMatEq[[i]] ))),
#                        cbind( matrix( 0, nrow( xMatEq[[i]] ), ncol( xMatAll )), xMatEq[[i]]))
#     nObsEq[i]   <-  length( yVecEq[[i]] )
#     nExogEq[i]   <-  ncol(xMatEq[[i]])
#     for(j in 1:nExogEq[i]) {
#       xnames <- c( xnames, paste("eq",as.character(i),colnames( xMatEq[[i]] )[j] ))
#     }
#   }

   # the previous lines are subtituted by the following,
   # because Ott Toomet reported that they might lead to
   # problems with special data sets. He suggested the
   # following lines, which are copied from the survreg
   # package
   # how were we called?
   call <- match.call() # get the original call
   m0 <- match.call( expand.dots = FALSE ) #-"- without ...-expansion
   pos <- which( names( m0 ) == "data" )
   if( length( pos ) == 1 ) {
      results$data.name <- as.character( m0[ pos ] )
   } else {
      results$data.name <- "unknown"
   }
   rm( pos )
   temp <- c("", "data", "weights", "subset", "na.action")
                  # arguments for model matrices
   m0 <- m0[match(temp, names(m0), nomatch = 0)]
            # positions of temp-arguments
   m0[[1]] <- as.name("model.frame")
            # find matrices for individual models
   for(i in 1:nEq ) {
      m <- m0
      Terms <- terms(eqns[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      weights <- model.extract(m, "weights")
      yVecEq[[i]] <- model.extract(m, "response")
      xMatEq[[i]] <- model.matrix(Terms, m)
      yVecAll <- c(yVecAll,yVecEq[[i]])
      xMatAll <- rbind( cbind( xMatAll, matrix( 0, nrow( xMatAll ), ncol( xMatEq[[i]] ))),
                  cbind( matrix( 0, nrow( xMatEq[[i]] ), ncol( xMatAll )), xMatEq[[i]]))
      nObsEq[i] <- length( yVecEq[[i]] )
      nExogEq[i] <- ncol(xMatEq[[i]])
      for(j in 1:nExogEq[i]) {
         xnames <- c( xnames, paste("eq",as.character(i),colnames( xMatEq[[i]] )[j] ))
      }
   }
   if( nEq > 1 ) {
      if( var ( nObsEq ) != 0 ) {
         stop( "Systems with unequal numbers of observations are not supported yet." )
      }
   }
   if( sum( nObsEq ) > 1000 && !saveMemory ) {
      warning( paste( "You have more than 1000 observations.",
         "Setting argument 'saveMemory' to TRUE speeds up",
         "the estimation. Estimation of larger data sets might even",
         "require this setting.\n" ) )
   }

   nObsAll  <- sum( nObsEq )  # total number of observations of all equations
   nExogAll <- sum( nExogEq ) # total number of exogenous variables/(unrestricted) coefficients in all equations
   nExogLiAll <- nExogAll     # total number of linear independent coefficients in all equations
   nExogLiEq  <- nExogEq      # total number of linear independent coefficients in each equation
   if(!is.null(TX)) {
      XU <- xMatAll
      xMatAll  <- XU %*% TX
      nExogLiAll <- nExogLiAll - ( nrow( TX ) - ncol( TX ) )
      for(i in 1:nEq) {
         nExogLiEq[i] <- ncol(xMatAll)
         for(j in 1: ncol(xMatAll) ) {
            if(sum(xMatAll[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])),j]^2)==0) {
               nExogLiEq[i] <- nExogLiEq[i]-1
            }
         }
      }
   }
   if(!is.null(R.restr)) {
      nExogLiAll  <- nExogLiAll - nrow(R.restr)
      if(is.null(TX)) {
         for(j in 1:nrow(R.restr)) {
            for(i in 1:nEq) {  # search for restrictions that are NOT cross-equation
               if( sum( R.restr[ j, (1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))]^2) ==
                   sum(R.restr[j,]^2)) {
                  nExogLiEq[i] <- nExogLiEq[i]-1
               }
            }
         }
      }
    }
    df <- nObsEq - nExogLiEq    # degress of freedom of each equation

  ## only for OLS, WLS and SUR estimation
  if(method=="OLS" | method=="WLS" | method=="SUR" | method=="WSUR") {
    if(is.null(R.restr)) {
      b <- solve( crossprod( xMatAll ), crossprod( xMatAll, yVecAll ), tol=solvetol )
               # estimated coefficients
    } else {
      W <- rbind( cbind( t(xMatAll) %*% xMatAll, t(R.restr) ),
                  cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
      V <- rbind( t(xMatAll) %*% yVecAll , q.restr )
      b <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(xMatAll)]
    }
  }

  ## only for OLS estimation
  if(method=="OLS") {
    resids <- yVecAll - xMatAll %*% b                                        # residuals
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )               # residual covariance matrix
      bcov <- .calcGLS( x = xMatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
                    # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nExogLiAll, methodRCov = methodRCov )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( xMatAll ), tol=solvetol )
                          # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coefficient covariance matrix
      }
    }
  }

  ## only for WLS estimation
  if( method == "WLS" | method == "WSUR" ) {
    bl    <- b   # coefficients of previous step
    bdif  <- b   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "WLS" ) ) {
      iter  <- iter+1
      bl    <- b                # coefficients of previous step
      resids <- yVecAll - xMatAll %*% b     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      b  <- .calcGLS( x = xMatAll, y = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol ) # coefficients
      bdif <- b-bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
       # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% b                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for SUR estimation
  if( method == "SUR" | method == "WSUR" ) {
    bl    <- b    # coefficients of previous step
    bdif  <- b    # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%b                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = centerResiduals,
         solvetol = solvetol )
      b <- .calcGLS( x = xMatAll, y = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )     # coefficients
      bdif <- b-bl  # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
            # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% b                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
      method == "W3SLS" ) {
    for(i in 1:nEq) {
      if(is.list(inst)) {
         instl[[i]] <- inst[[i]]
      } else {
         instl[[i]] <- inst
      }
    }
    Xf <- array(0,c(0,ncol(xMatAll)))       # fitted X values
    H  <- matrix( 0, 0, 0 )           # stacked matrices of all instruments
    h  <- list()
    for(i in 1:nEq) {
      Xi <- xMatAll[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])),]
            # regressors of the ith equation (including zeros)
      #h[[i]] <- model.matrix( instl[[i]] )
      # the following lines have been substituted for the previous
      # line due to changes in the data handling.
      # code provided by Ott Toomet
      m <- m0
      Terms <- terms(instl[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      h[[i]] <- model.matrix(Terms, m)
      if( nrow( h[[ i ]] ) != nrow( Xi ) ) {
         stop( paste( "The instruments and the regressors of equation",
            as.character( i ), "have different numbers of observations." ) )
      }
      # extract instrument matrix
      Xf <- rbind(Xf, h[[i]] %*% solve( crossprod( h[[i]]) , tol=solvetol )
              %*% crossprod( h[[i]], Xi ))       # 'fitted' X-values
      H  <-  rbind( cbind( H, matrix( 0, nrow( H ), ncol( h[[i]] ))),
                         cbind( matrix( 0, nrow( h[[i]] ), ncol( H )), h[[i]]))

    }
    if(is.null(R.restr)) {
      b <- solve( crossprod( Xf ), crossprod( Xf, yVecAll ), tol=solvetol )
         # 2nd stage coefficients
    } else {
      W <- rbind( cbind( crossprod(Xf), t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
      V <- rbind( t(Xf) %*% yVecAll , q.restr )
      b <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(xMatAll)] # restricted coefficients
    }
    b2 <- b
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- yVecAll - xMatAll %*% b                        # residuals
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nExogLiAll, methodRCov = methodRCov )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( Xf ), tol=solvetol )
                  # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if( method == "W2SLS" | method == "W3SLS" ) {
    bl     <- b   # coefficients of previous step
    bdif   <- b   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "W2LS" ) ) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%b                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      b <- .calcGLS( x = Xf, y = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )          # (unrestr.) coeffic.
      bdif <- b - bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% b                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 3SLS estimation
  if( method == "3SLS" | method == "W3SLS" ) {
    bl     <- b  # coefficients of previous step
    bdif   <- b  # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- b                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%b                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = centerResiduals, solvetol = solvetol )
      if(method3sls=="GLS") {
         b <- .calcGLS( x = Xf, y = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      if(method3sls=="IV") {
         b <- .calcGLS( x = Xf, x2 = xMatAll, y = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )   # (unrestr.) coeffic.
      }
      if(method3sls=="GMM") {
        HtOmega <- .calcXtOmegaInv( x = H, sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
        if(is.null(R.restr)) {
          b <- solve(t(xMatAll) %*% H %*% solve( HtOmega %*%
                 H, tol=solvetol) %*% t(H) %*% xMatAll, tol=solvetol) %*% t(xMatAll) %*% H %*%
                 solve( HtOmega %*%
                 H, tol=solvetol) %*% t(H) %*% yVecAll  #(unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(xMatAll) %*% H %*% solve( HtOmega
                              %*% H, tol=solvetol) %*% t(H) %*% xMatAll, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(xMatAll) %*% H %*% solve( HtOmega
                      %*% H, tol=solvetol) %*% t(H) %*% yVecAll , q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(method3sls=="Schmidt") {
        XftOmegaInv <- .calcXtOmegaInv( x = Xf, sigma = rcov, nObsEq = nObsEq,
           solvetol = solvetol )
        if(is.null(R.restr)) {
          b <- solve( t(Xf) %*% t( XftOmegaInv ), tol=solvetol) %*% ( XftOmegaInv
                      %*% H %*% solve( crossprod( H ), tol=solvetol ) %*% crossprod(H, yVecAll) )
                           # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(Xf) %*% t( XftOmegaInv ), t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( XftOmegaInv %*% H %*% solve( crossprod( H ), tol=solvetol ) %*%
                      crossprod( H, yVecAll ), q.restr )
          Winv <- solve( W, tol=solvetol )
          b <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(method3sls=="EViews") {
         b <- b2 + .calcGLS( x = Xf, y = ( yVecAll -  xMatAll %*% b2 ),
            R.restr = R.restr, q.restr = q.restr, sigma = rcov,
            nObsEq = nObsEq, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      bdif <- b - bl # difference of coefficients between this and previous step
    }
    if(method3sls=="GLS") {
       bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    }
    if(method3sls=="IV") {
       bcov <- .calcGLS( x = Xf, x2 = xMatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
    }
    if(method3sls=="GMM") {
      if(is.null(R.restr)) {
        bcov <- solve( t(xMatAll) %*% H %*% solve( HtOmega %*% H, tol=solvetol ) %*%
           t(H) %*% xMatAll, tol=solvetol )
                # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(xMatAll),1:ncol(xMatAll)] # coefficient covariance matrix
      }
    }
    if(method3sls=="Schmidt") {
      XftOmegaInv <- .calcXtOmegaInv( x = Xf, sigma = rcov, nObsEq = nObsEq,
         solvetol = solvetol )
      PH <- H %*%  solve( t(H) %*% H, tol=solvetol ) %*% t(H)
      PHOmega <- .calcXtOmegaInv( x = t( PH ), sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
      if(is.null(R.restr)) {
         bcov <- solve( XftOmegaInv %*% Xf, tol=solvetol ) %*%
            XftOmegaInv %*% PHOmega %*%
            PH %*% t( XftOmegaInv ) %*% solve( XftOmegaInv %*% Xf, tol=solvetol )
                  # final step coefficient covariance matrix
      } else {
         VV <- XftOmegaInv %*% PHOmega %*%
            PH %*% t( XftOmegaInv )
         VV <- rbind( cbind( VV, matrix( 0, nrow( VV ), nrow( R.restr ) ) ),
            matrix( 0, nrow( R.restr ), nrow( VV ) + nrow( R.restr ) ) )
         bcov <- ( Winv %*% VV %*% Winv )[ 1:ncol(xMatAll), 1:ncol(xMatAll) ]
                  # coefficient covariance matrix
      }
    }
    if(method3sls=="EViews") {
       bcov <- .calcGLS( x = Xf, R.restr = R.restr, q.restr = q.restr, 
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # final step coefficient covariance matrix
    }
    resids <- yVecAll - xMatAll %*% b                        # residuals
  }

  ## FIML estimation
  if( method == "FIML" ) {
    fimlResult <- .systemfitFiml( systemfitCall = call, nObsEq = nObsEq,
      nCoefEq = nExogLiEq, yVec = yVecAll, xMat = xMatAll, xEq = xMatEq, methodRCov = methodRCov,
      centerResiduals = centerResiduals, solvetol = solvetol )
    #print( fimlResult )
    b <- fimlResult$coef
    bcov <- fimlResult$coefCov
    resids <- fimlResult$resids
  }

  ## for all estimation methods
  fitted <- xMatAll %*% b                              # fitted endogenous values
  bt     <- NULL
  btcov  <- NULL
  if(!is.null(TX)) {
    bt <- b
    b  <- TX %*% bt
    btcov <- bcov
    bcov  <- TX %*% btcov %*% t(TX)
  }
  se     <- diag(bcov)^0.5                       # standard errors of all estimated coefficients
  t      <- b/se                                 # t-values of all estimated coefficients
  if(probdfsys) {
    prob <- 2*( 1-pt(abs(t), nObsAll - nExogLiAll))
       # p-values of all estimated coefficients
  } else {
    prob <- matrix( 0, 0, 1 )                    # p-values of all estimated coefficients
  }



  ## equation wise results
  for(i in 1:nEq) {
    residi[[i]] <- resids[ ( 1 + sum(nObsEq[1:i]) -nObsEq[i] ):( sum(nObsEq[1:i]) ) ]
    bi     <- b[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))]
              # estimated coefficients of equation i
    sei    <- c(se[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))])
              # std. errors of est. param. of equation i
    ti     <- c(t[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))])
              # t-values of estim. param. of equation i
    bcovi  <- bcov[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i])),(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))]
              # covariance matrix of estimated coefficients of equation i
    if(probdfsys) {
      probi <- c(prob[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))])
               # p-values of estim. param. of equation i
    } else {
      probi <- 2*( 1 - pt(abs(ti), df[i] ))
               # p-values of estim. param. of equation i
      prob <- c(prob,probi) # p-values of all estimated coefficients
    }

    # set names
    names( bi ) <- colnames( xMatEq[[i]] )
    names( sei ) <- colnames( xMatEq[[i]] )
    names( ti ) <- colnames( xMatEq[[i]] )
    names( probi ) <- colnames( xMatEq[[i]] )
    colnames( bcovi ) <- colnames( xMatEq[[i]] )
    rownames( bcovi ) <- colnames( xMatEq[[i]] )

    ssr    <- sum(residi[[i]]^2)                         # sum of squared residuals
    mse    <- ssr/df[i]                                  # estimated variance of residuals
    rmse   <- sqrt( mse )                                # estimated standard error of residuals
    r2     <- 1 - ssr/(t(yVecEq[[i]])%*%yVecEq[[i]]-nObsEq[i]*mean(yVecEq[[i]])^2)
    adjr2  <- 1 - ((nObsEq[i]-1)/df[i])*(1-r2)
    fittedi <- fitted[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
    #datai  <- model.frame( eqns[[i]] )
    # the following lines have to be substituted for the previous
    # line due to changes in the data handling.
    # code provided by Ott Toomet
    m <- m0
    Terms <- terms( eqns[[i]], data = data)
    m$formula <- Terms
    m <- eval(m, parent.frame())
    # datai <- model.frame(Terms, m)
    # the following lines are substituted for the previous line to
    # allow transformed variables in the formulas (e.g. "log(x1)")
    # code provided by Mikko Pakkanen
    resp <- model.extract( m, "response" )
    # using model.matrix instead of model.frame, need to get the output
    # variable separately
    datai <- data.frame( cbind( resp, ( model.matrix( Terms, m ) )[ , -1 ] ) )
    # I guess there's a better way to extract the name of the output variable?
    names( datai )[1] <- as.character( terms( eqns[[ i ]] ) )[2]
    rm( resp )
    if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
        method == "W3SLS" ) {
      #datai <- cbind( datai, model.frame( instl[[i]] ))
      # the following lines have to be substituted for the previous
      # line due to changes in the data handling.
      # code provided by Ott Toomet
      m <- m0
      Terms <- terms(instl[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      # datai <- cbind( datai, model.frame(Terms, m))
      # the following lines are substituted for the previous line to
      # allow transformed variables in the formulas (e.g. "log(x1)")
      # code provided by Mikko Pakkanen
      datai <- cbind( datai,
         as.data.frame( ( model.matrix( Terms, m ) )[ , -1 ] ) )
    }

    if(i==1) {
      alldata <- datai                    # dataframe for all data used for estimation
    } else {
      alldata <- cbind( alldata, datai )  # dataframe for all data used for estimation
    }

    ## build the "return" structure for the equations
    resulti$method       <- method
    resulti$i            <- i               # equation number
    resulti$eqnlabel     <- eqnlabels[[i]]
    resulti$formula      <- eqns[[i]]
    resulti$nObs         <- nObsEq[i]       # number of observations
    resulti$nExog        <- nExogEq[i]      # number of exogenous variables/coefficients
    resulti$nExogLi      <- nExogLiEq[i]    # number of linear independent coefficients
    resulti$df           <- df[i]           # degrees of freedom of residuals
    resulti$dfSys        <- nObsAll- nExogLiAll
       # degrees of freedom of residuals of the whole system
    resulti$probdfsys    <- probdfsys       #
    resulti$b            <- c( bi )         # estimated coefficients
    resulti$se           <- c( sei )        # standard errors of estimated coefficients
    resulti$t            <- c( ti )         # t-values of estimated coefficients
    resulti$p            <- c( probi )      # p-values of estimated coefficients
    resulti$bcov         <- bcovi           # covariance matrix of estimated coefficients
    resulti$yVec         <- yVecEq[[i]]     # vector of endogenous variables
    resulti$xMat         <- xMatEq[[i]]     # matrix of regressors
    resulti$data         <- datai           # data frame of this equation (incl. instruments)
    resulti$fitted       <- fittedi         # fitted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$mse          <- mse             # estimated variance of the residuals (mean squared error)
    resulti$s2           <- mse             #        the same (sigma hat squared)
    resulti$rmse         <- rmse            # estimated standard error of the residuals
    resulti$s            <- rmse            #        the same (sigma hat)
    resulti$r2           <- r2              # R-sqared value
    resulti$adjr2        <- adjr2           # adjusted R-squared value
    if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
        method == "W3SLS" ) {
      resulti$inst         <- instl[[i]]
      resulti$h            <- h[[i]]          # matrix of instrumental variables
    }
    class(resulti)        <- "systemfit.equation"
    results$eq[[i]]      <- resulti
  }

  ## results of the total system
  #olsr2 <- 1 - t(resids) %*% resids / ( t(yVecAll) %*% ( diag(1,nEq,nEq)     # OLS system R2
  # %x% ( diag( 1, nObsEq[1], nObsEq[1]) - rep(1, nObsEq[1]) %*% t(rep(1,nObsEq[1])) / nObsEq[1])) %*% yVecAll)
  # the following lines are substituted for the previous 2 lines to increase
  # speed ( idea suggested by Ott Toomet )
   meanY <- numeric(length(yVecAll)) # compute mean of Y by equations
   for(i in 1:nEq) {
      meanY[ (1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])) ] <-
         mean( yVecAll[ (1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])) ])
   }
   olsr2 <- 1 - t(resids) %*% resids / sum( ( yVecAll - meanY )^2 )
                        # OLS system R2
  if( method == "SUR" | method == "WSUR" | method == "3SLS" |
      method == "W3SLS" ) {
    rcovest <- rcov                   # residual covariance matrix used for estimation
  }
  rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
      nCoefEq = nExogLiEq, xEq = xMatEq, centered = centerResiduals, solvetol = solvetol )
  drcov <- det(rcov, tol=solvetol)
  if( !saveMemory ) {
#       # original formula from McElroy (1977)
#       mcelr2 <- 1 - ( t(resids) %*% ( solve(rcov, tol=solvetol) %x%
#                 diag(1, nObsEq[1],nObsEq[1])) %*% resids ) /
#                 ( t(yVecAll) %*% ( solve(rcov, tol=solvetol ) %x%
#                 ( diag(1,nObsEq[1],nObsEq[1] ) - rep(1,nObsEq[1]) %*%
#                 t(rep(1,nObsEq[1])) / nObsEq[1] )) %*% yVecAll )   # McElroy's (1977a) R2
      # first formula from Greene (2003, p. 345) (numerator modified to save memory)
      rtOmega <- .calcXtOmegaInv( x = resids, sigma = rcov, nObsEq = nObsEq,
         solvetol = solvetol )
      yCov <- .calcRCov( yVecAll, methodRCov = "noDfCor", nObsEq = nObsEq, centered = TRUE,
         solvetol = solvetol )
      residCovInv <- solve( rcov, tol = solvetol )
      denominator <- 0
      for( i in 1:nEq ) {
         for( j in 1:nEq ) {
            denominator <- denominator + residCovInv[ i, j ] * yCov[ i, j ] * nObsEq[1]
         }
      }
      mcelr2 <- 1 - ( rtOmega %*% resids ) / denominator
#       # second formula from Greene (2003, p. 345)
#        yCov <- sum(diag(.calcRCov( yVecAll, methodRCov = "noDfCor", nObsEq = nObsEq, centered = TRUE,
#           solvetol = solvetol )))
#        yCov <- drop( t(yVecAll-mean(yVecAll)) %*% (yVecAll-mean(yVecAll)) / sum(nObsEq) )
#       yCov <- .calcRCov( yVecAll, methodRCov = "geomean", nObsEq = nObsEq, nCoefEq=rep(1,nEq),
#          centered = TRUE, solvetol = solvetol )
#        mcelr2 <- 1 - nEq / sum( diag( solve( rcov, tol = solvetol ) * yCov ) )
  } else {
     mcelr2 <- NA
  }

  b              <- c(b)
  names(b)       <- xnames
  se             <- c(se)
  names(se)      <- xnames
  t              <- c(t)
  names(t)       <- xnames
  prob           <- c(prob)
  names(prob)    <- xnames
  colnames( bcov ) <- xnames
  rownames( bcov ) <- xnames
  colnames( rcov ) <- eqnlabels
  rownames( rcov ) <- eqnlabels

  ## build the "return" structure for the whole system
  results$method  <- method
  results$nEq     <- nEq            # number of equations
  results$nObsAll <- nObsAll        # total number of observations of all equations
  results$nObsEq  <- nObsEq         # number of observations in each equation
  results$nExogAll <- nExogAll      # total number of exogenous variables/coefficients in all equations
  results$nExogEq <- nExogEq        # number of exogenous variables/coefficients
                                       # in each equation
  results$nExogLiAll <- nExogLiAll
     # total number of linear independent coefficients of all equations
  results$nExogLiEq  <- nExogLiEq
     # number of linear independent coefficients in each equation
  results$df      <- nObsAll - nExogLiAll
     # degrees of freedom of the whole system
  results$b       <- b              # all estimated coefficients
  results$bt      <- bt             # transformed vector of estimated coefficients
  results$se      <- se             # standard errors of estimated coefficients
  results$t       <- t              # t-values of estimated coefficients
  results$p       <- prob           # p-values of estimated coefficients
  results$bcov    <- bcov           # coefficients covariance matrix
  results$btcov   <- btcov          # covariance matrix for transformed coeff. vector
  results$rcov    <- rcov           # residual covarance matrix
  results$drcov   <- drcov          # determinant of residual covarance matrix
  results$olsr2   <- olsr2          # R-squared value of the equation system
  results$iter    <- iter           # residual correlation matrix
  results$yVec    <- yVecAll        # vector of all (stacked) endogenous variables
  results$xMat    <- xMatAll        # matrix of all (diagonally stacked) regressors
  results$resids  <- resids         # vector of all (stacked) residuals
  results$data    <- alldata        # data frame for all data used in the system
  if( method == "2SLS" | method == "W2SLS" | method == "3SLS" |
      method == "W3SLS" ) {
    results$h       <- H            # matrix of all (diagonally stacked) instr. variables
    results$xHat    <- Xf           # matrix of "fitted" regressors
  }
  if( method == "SUR" | method == "WSUR" | method == "3SLS" |
      method == "W3SLS" ) {
    results$rcovest <- rcovest      # residual covarance matrix used for estimation
    results$mcelr2  <- mcelr2       # McElroy's R-squared value for the equation system
  }
  results$R.restr <- R.restr
  results$q.restr <- q.restr
  results$TX      <- TX
  results$maxiter <- maxiter
  results$tol     <- tol
  results$methodRCov     <- methodRCov
  results$method3sls     <- method3sls
  results$probdfsys       <- probdfsys
  results$single.eq.sigma <- single.eq.sigma
  results$solvetol        <- solvetol
  class(results)  <- "systemfit"

  results
}


## print the (summary) results that belong to the whole system
summary.systemfit <- function(object,...) {
   object$coef <- cbind( object$b, object$se, object$t, object$p )
   colnames( object$coef ) <- c( "Estimate", "Std. Error",
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

## print a few results of the whole system
print.systemfit <- function( x, digits=6,... ) {

   save.digits <- unlist(options(digits=digits))
   on.exit(options(digits=save.digits))

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
   cat( "Coefficients:\n" )
   print( x$b )
   invisible( x )
}

## print the (summary) results for a single equation
summary.systemfit.equation <- function(object,...) {
   object$coef <- cbind( object$b, object$se, object$t, object$p )
   colnames( object$coef ) <- c( "Estimate", "Std. Error",
      "t value", "Pr(>|t|)" )
   class( object ) <- "summary.systemfit.equation"
   return( object )
}


## print summary results for a single equation
print.summary.systemfit.equation <- function( x, digits=6, ... ) {

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  cat("\n")
  if( is.null( x$eqnlabel ) ) {
    cat( x$method, " estimates for equation ", x$i, "\n", sep = "" )
  } else {
    cat( x$method, " estimates for '", x$eqnlabel,
         "' (equation ", x$i, ")\n", sep = "" )
  }

  cat("Model Formula: ")
  print(x$formula)
  if(!is.null(x$inst)) {
    cat("Instruments: ")
    print(x$inst)
  }
  cat("\n")

  Signif <- symnum(x$p, corr = FALSE, na = FALSE,
                   cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
                   symbols   = c( "***", "**", "*", "." ," " ))

  table <- cbind(round( x$b,  digits ),
                 round( x$se, digits ),
                 round( x$t,  digits ),
                 round( x$p,  digits ),
                 Signif)

  rownames(table) <- names(x$b)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round(x$s, digits),
            "on", x$df, "degrees of freedom\n"))
            # s ist the variance, isn't it???

  cat( paste( "Number of observations:", round(x$nObs, digits),
              "Degrees of Freedom:", round(x$df, digits),"\n" ) )

  cat( paste( "SSR:",     round(x$ssr,    digits),
              "MSE:", round(x$mse, digits),
              "Root MSE:",   round(x$rmse,  digits), "\n" ) )

  cat( paste( "Multiple R-Squared:", round(x$r2,    digits),
              "Adjusted R-Squared:", round(x$adjr2, digits),
              "\n" ) )
  cat("\n")
  invisible( x )
}

## print a few results for a single equation
print.systemfit.equation <- function( x, digits=6, ... ) {

   save.digits <- unlist(options(digits=digits))
   on.exit(options(digits=save.digits))

   cat("\n")
   if( is.null( x$eqnlabel ) ) {
      cat( x$method, " estimates for equation ", x$i, "\n", sep = "" )
   } else {
      cat( x$method, " estimates for '", x$eqnlabel,
            "' (equation ", x$i, ")\n", sep = "" )
   }

   cat("Model Formula: ")
   print(x$formula)
   if(!is.null(x$inst)) {
      cat("Instruments: ")
      print(x$inst)
   }

   cat("\nCoefficients:")
   print( x$b )
   invisible( x )
}


## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, data=object$data,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95, ... ) {
   attach(data); on.exit( detach( data ) )

   predicted <- data.frame( obs=seq( nrow( data ) ) )
   colnames( predicted ) <- as.character( 1:ncol( predicted ) )
   g       <- object$nEq
   nObsEq       <- array(NA,c(g))
   eqns    <- list()
   xMatEq  <- list()               # regressors equation-wise
   xMatAll <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
   for(i in 1:g )  {
      eqns[[i]] <- object$eq[[i]]$formula
      xMatEq[[i]] <-  model.matrix( eqns[[i]] )
      xMatAll      <-  rbind( cbind( xMatAll, matrix( 0, nrow( xMatAll ), ncol( xMatEq[[i]] ))),
                       cbind( matrix( 0, nrow( xMatEq[[i]] ), ncol( xMatAll )), xMatEq[[i]]))
      nObsEq[i]   <-  nrow( xMatEq[[i]] )
   }
   yVecAll <- xMatAll %*% object$b
   if( object$method == "SUR" | object$method == "WSUR" |
       object$method == "3SLS" | object$method == "W3SLS") {
      if( se.fit | interval == "confidence" ) {
         ycovc <- xMatAll %*% object$bcov %*% t(xMatAll)
      }
      if( se.pred | interval == "prediction" ) {
         ycovp <- xMatAll %*% object$bcov %*% t(xMatAll) + object$rcov %x% diag(1,nObsEq[1],nObsEq[1])
      }
   }
   for(i in 1:g) {
      # fitted values
      Yi <- yVecAll[(1+sum(nObsEq[1:i])-nObsEq[i]):sum(nObsEq[1:i]),]
      predicted <- cbind( predicted, Yi )
      names( predicted )[ length( predicted ) ] <- paste( "eq", as.character(i),
                                                          ".pred", sep="" )
      # calculate variance covariance matrices
      if( se.fit | interval == "confidence" ) {
         if( object$method == "SUR" | object$method == "WSUR" |
             object$method == "3SLS" | object$method== "W3SLS" ) {
            ycovci <- ycovc[ ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ),
                             ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ) ]
         } else {
            ycovci <- xMatEq[[i]] %*% object$eq[[i]]$bcov %*% t(xMatEq[[i]])
         }
      }
      if( se.pred | interval == "prediction" ) {
         if( object$method == "SUR" | object$method == "WSUR" |
             object$method == "3SLS" | object$method == "W3SLS" ) {
            ycovpi <- ycovp[ ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ),
                            ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ) ]
         } else {
            ycovpi <- xMatEq[[i]] %*% object$eq[[i]]$bcov %*% t(xMatEq[[i]]) +
                                 object$eq[[i]]$s2
         }
      }
      # standard errors of fitted values
      if( se.fit ) {
         if(nrow(data)==1) {
            predicted <- cbind( predicted, sqrt( ycovci ) )
         } else {
            predicted <- cbind( predicted, sqrt( diag( ycovci ) ) )
         }
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".se.fit", sep="" )
      }
      # standard errors of prediction
      if( se.pred ) {
         if(nrow(data)==1) {
            predicted <- cbind( predicted, sqrt( ycovpi ) )
         } else {
            predicted <- cbind( predicted, sqrt( diag( ycovpi ) ) )
         }
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".se.pred", sep="" )
      }

      # confidence intervals
      if( interval == "confidence" ) {
         if( object$probdfsys ) {
            tval   <- qt( 1 - ( 1- level )/2, object$df )
         } else {
            tval   <- qt( 1 - ( 1- level )/2, object$eq[[i]]$df )
         }
         if( nrow(data)==1 ) {
            seci    <- sqrt( ycovci )
         } else {
            seci    <- sqrt( diag( ycovci ) )
         }
         predicted <- cbind( predicted, Yi - ( tval * seci ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".lwr", sep="" )
         predicted <- cbind( predicted, Yi + ( tval * seci ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".upr", sep="" )
      }
      # prediction intervals
      if( interval == "prediction" ) {
         if( object$probdfsys ) {
            tval   <- qt( 1 - ( 1- level )/2, object$df )
         } else {
            tval   <- qt( 1 - ( 1- level )/2, object$eq[[i]]$df )
         }
         if(nrow(data)==1) {
            sepi <- sqrt( ycovpi )
         } else {
            sepi <- sqrt( diag( ycovpi ) )
         }
         predicted <- cbind( predicted, Yi - ( tval * sepi ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".lwr", sep="" )
         predicted <- cbind( predicted, Yi + ( tval * sepi ) )
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".upr", sep="" )
      }
   }
   predicted[ 2: length( predicted ) ]
}

## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit.equation <- function( object, data=object$data, ... ) {
   attach( data ); on.exit( detach( data ) )
   x <-  model.matrix( object$formula )
   predicted <- drop( x %*% object$b )
   predicted
}

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
    xik    <- results$eq[[eqni]]$xMat[i,]
    xjk    <- results$eq[[eqnj]]$xMat[i,]
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
    xik    <- resultsi$eq[[eqni]]$xMat[i,]
    top    <- sqrt( xik %*% resultsi$eq[[eqni]]$bcov %*% xik )
    bottom <- sqrt( xik %*% resultsj$eq[[eqni]]$bcov %*% xik )
    rk     <- top / bottom
    ratio  <- rbind( ratio, rk )
  }
  ratio
}


## return all coefficients
coef.systemfit <- function( object, ... ) {
   object$b
}

## return all coefficients, std.errors, t-values and p-values
coef.summary.systemfit <- function( object, ... ) {
   object$coef
}

## return the coefficients of a single equation
coef.systemfit.equation <- function( object, ... ) {
   object$b
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

## return the variance covariance matrix of the coefficients
confint.systemfit <- function( object, parm = NULL, level = 0.95, ... ) {
   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- array( NA, dim = c( length( object$b ), 2),
            dimnames = list( names( object$b ), pct ) )
   j <- 1
   for( i in 1:object$nEq ) {
      object$eq[[i]]$dfSys <- object$df
      object$eq[[i]]$probdfsys <- object$probdfsys
      ci[ j:(j+object$eq[[ i ]]$nExog-1), ] <- confint( object$eq[[ i ]] )
      j <- j + object$eq[[ i ]]$nExog
   }
   class( ci ) <- "confint.systemfit"
   ci
}

## return the variance covariance matrix of the coefficients of a single equation
confint.systemfit.equation <- function( object, parm = NULL, level = 0.95, ... ) {
   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- array( NA, dim = c( length( object$b ), 2),
            dimnames = list( names( object$b ), pct ) )
   if( object$ probdfsys ) {
      fac <- qt( a, object$dfSys )
   } else {
      fac <- qt( a, object$df )
   }
   ci[] <- object$b + object$se %o% fac
   class( ci ) <- "confint.systemfit"
   ci
}

## print the confidence intervals of the coefficients
print.confint.systemfit <- function( x, digits = 3, ... ) {
   print( unclass( round( x, digits = digits, ...) ) )
   invisible(x)
}

## return the fitted values
fitted.systemfit <- function( object, ... ) {
   fitted <- array( NA, c( length( object$eq[[1]]$fitted ), object$nEq ) )
   colnames( fitted ) <- as.character( 1:ncol( fitted ) )
   for(i in 1:object$nEq )  {
      fitted[ , i ]           <- object$eq[[ i ]]$fitted
      colnames( fitted )[ i ] <- paste( "eq", as.character(i), sep="" )
   }
   fitted
}

## return the fitted values of e single euation
fitted.systemfit.equation <- function( object, ... ) {
   object$fitted
}




