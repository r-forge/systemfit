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
                        single.eq.sigma=(is.null(R.restr) & is.null(TX)),
                        solvetol=.Machine$double.eps,
                        saveMemory = ( nrow( data ) * length( eqns ) > 1000 &&
                           length( data ) > 0 ) )
{

   ## some tests
   if(!( method %in% c( "OLS", "WLS", "SUR", "WSUR", "2SLS", "W2SLS", "3SLS",
         "W3SLS", "LIML", "FIML" ) ) ){
      stop( "The method must be 'OLS', 'WLS', 'SUR', 'WSUR',",
         " '2SLS', 'W2SLS', '3SLS', or 'W3SLS'" )
   }
   if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) &
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
  nObsEq  <- numeric( nEq ) # number of observations in each equation
  nExogEq <- numeric( nEq ) # number of exogenous variables /(unrestricted) coefficients
                                     # in each equation
  instl   <- list()         # list of the instruments for each equation
  ssr     <- numeric( nEq ) # sum of squared residuals of each equation
  sigma   <- numeric( nEq ) # estimated sigma (std. dev. of residuals) of each equation
  r2      <- numeric( nEq ) # R-squared value
  adjr2   <- numeric( nEq ) # adjusted R-squared value
  xnames  <- NULL           # names of regressors

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
  if( method %in% c( "OLS", "WLS", "SUR", "WSUR" ) ) {
    if(is.null(R.restr)) {
      coef <- solve( crossprod( xMatAll ), crossprod( xMatAll, yVecAll ), tol=solvetol )
               # estimated coefficients
    } else {
      W <- rbind( cbind( t(xMatAll) %*% xMatAll, t(R.restr) ),
                  cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
      V <- rbind( t(xMatAll) %*% yVecAll , q.restr )
      coef <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(xMatAll)]
    }
  }

  ## only for OLS estimation
  if(method=="OLS") {
    resids <- yVecAll - xMatAll %*% coef                                        # residuals
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )               # residual covariance matrix
      bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
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
  if( method %in% c( "WLS", "WSUR" ) ) {
    bl    <- coef   # coefficients of previous step
    bdif  <- coef   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "WLS" ) ) {
      iter  <- iter+1
      bl    <- coef                # coefficients of previous step
      resids <- yVecAll - xMatAll %*% coef     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      coef  <- .calcGLS( xMat = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol ) # coefficients
      bdif <- coef-bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
       # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for SUR estimation
  if( method %in% c( "SUR", "WSUR" ) ) {
    bl    <- coef    # coefficients of previous step
    bdif  <- coef    # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = centerResiduals,
         solvetol = solvetol )
      coef <- .calcGLS( xMat = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )     # coefficients
      bdif <- coef-bl  # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
            # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
    for(i in 1:nEq) {
      if(is.list(inst)) {
         instl[[i]] <- inst[[i]]
      } else {
         instl[[i]] <- inst
      }
    }
    xMatHatAll <- matrix( 0, 0, ncol( xMatAll ) ) # fitted X values
    hMatAll  <- matrix( 0, 0, 0 )           # stacked matrices of all instruments
    hMatEq  <- list()
    for(i in 1:nEq) {
      Xi <- xMatAll[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])),]
            # regressors of the ith equation (including zeros)
      #hMatEq[[i]] <- model.matrix( instl[[i]] )
      # the following lines have been substituted for the previous
      # line due to changes in the data handling.
      # code provided by Ott Toomet
      m <- m0
      Terms <- terms(instl[[i]], data = data)
      m$formula <- Terms
      m <- eval(m, parent.frame())
      hMatEq[[i]] <- model.matrix(Terms, m)
      if( nrow( hMatEq[[ i ]] ) != nrow( Xi ) ) {
         stop( paste( "The instruments and the regressors of equation",
            as.character( i ), "have different numbers of observations." ) )
      }
      # extract instrument matrix
      xMatHatAll <- rbind(xMatHatAll, hMatEq[[i]] %*% solve( crossprod( hMatEq[[i]]) , tol=solvetol )
              %*% crossprod( hMatEq[[i]], Xi ))       # 'fitted' X-values
      hMatAll  <-  rbind( cbind( hMatAll, matrix( 0, nrow( hMatAll ), ncol( hMatEq[[i]] ))),
                         cbind( matrix( 0, nrow( hMatEq[[i]] ), ncol( hMatAll )), hMatEq[[i]]))

    }
    if(is.null(R.restr)) {
      coef <- solve( crossprod( xMatHatAll ), crossprod( xMatHatAll, yVecAll ), tol=solvetol )
         # 2nd stage coefficients
    } else {
      W <- rbind( cbind( crossprod(xMatHatAll), t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
      V <- rbind( t(xMatHatAll) %*% yVecAll , q.restr )
      coef <- ( solve( W, tol=solvetol ) %*% V )[1:ncol(xMatAll)] # restricted coefficients
    }
    b2 <- coef
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    if(single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nExogLiAll, methodRCov = methodRCov )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( xMatHatAll ), tol=solvetol )
                  # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if( method %in% c( "W2SLS", "W3SLS" ) ) {
    bl     <- coef   # coefficients of previous step
    bdif   <- coef   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter^( method == "W2LS" ) ) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = centerResiduals,
         solvetol = solvetol )
      coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )          # (unrestr.) coeffic.
      bdif <- coef - bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 3SLS estimation
  if( method %in% c( "3SLS", "W3SLS" ) ) {
    bl     <- coef  # coefficients of previous step
    bdif   <- coef  # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>tol & iter < maxiter) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = centerResiduals, solvetol = solvetol )
      if(method3sls=="GLS") {
         coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      if(method3sls=="IV") {
         coef <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )   # (unrestr.) coeffic.
      }
      if(method3sls=="GMM") {
        HtOmega <- .calcXtOmegaInv( xMat = hMatAll, sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
        if(is.null(R.restr)) {
          coef <- solve(t(xMatAll) %*% hMatAll %*% solve( HtOmega %*%
                 hMatAll, tol=solvetol) %*% t(hMatAll) %*% xMatAll, tol=solvetol) %*% t(xMatAll) %*% hMatAll %*%
                 solve( HtOmega %*%
                 hMatAll, tol=solvetol) %*% t(hMatAll) %*% yVecAll  #(unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(xMatAll) %*% hMatAll %*% solve( HtOmega
                              %*% hMatAll, tol=solvetol) %*% t(hMatAll) %*% xMatAll, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(xMatAll) %*% hMatAll %*% solve( HtOmega
                      %*% hMatAll, tol=solvetol) %*% t(hMatAll) %*% yVecAll , q.restr )
          Winv <- solve( W, tol=solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(method3sls=="Schmidt") {
        xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov, nObsEq = nObsEq,
           solvetol = solvetol )
        if(is.null(R.restr)) {
          coef <- solve( t(xMatHatAll) %*% t( xMatHatOmegaInv ), tol=solvetol) %*% ( xMatHatOmegaInv
                      %*% hMatAll %*% solve( crossprod( hMatAll ), tol=solvetol ) %*% crossprod(hMatAll, yVecAll) )
                           # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(xMatHatAll) %*% t( xMatHatOmegaInv ), t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( xMatHatOmegaInv %*% hMatAll %*% solve( crossprod( hMatAll ), tol=solvetol ) %*%
                      crossprod( hMatAll, yVecAll ), q.restr )
          Winv <- solve( W, tol=solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(method3sls=="EViews") {
         coef <- b2 + .calcGLS( xMat = xMatHatAll, yVec = ( yVecAll -  xMatAll %*% b2 ),
            R.restr = R.restr, q.restr = q.restr, sigma = rcov,
            nObsEq = nObsEq, solvetol = solvetol )  # (unrestr.) coeffic.
      }
      bdif <- coef - bl # difference of coefficients between this and previous step
    }
    if(method3sls=="GLS") {
       bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # coefficient covariance matrix
    }
    if(method3sls=="IV") {
       bcov <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )
    }
    if(method3sls=="GMM") {
      if(is.null(R.restr)) {
        bcov <- solve( t(xMatAll) %*% hMatAll %*% solve( HtOmega %*% hMatAll, tol=solvetol ) %*%
           t(hMatAll) %*% xMatAll, tol=solvetol )
                # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(xMatAll),1:ncol(xMatAll)] # coefficient covariance matrix
      }
    }
    if(method3sls=="Schmidt") {
      xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov, nObsEq = nObsEq,
         solvetol = solvetol )
      PH <- hMatAll %*%  solve( t(hMatAll) %*% hMatAll, tol=solvetol ) %*% t(hMatAll)
      PHOmega <- .calcXtOmegaInv( xMat = t( PH ), sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
      if(is.null(R.restr)) {
         bcov <- solve( xMatHatOmegaInv %*% xMatHatAll, tol=solvetol ) %*%
            xMatHatOmegaInv %*% PHOmega %*%
            PH %*% t( xMatHatOmegaInv ) %*% solve( xMatHatOmegaInv %*% xMatHatAll, tol=solvetol )
                  # final step coefficient covariance matrix
      } else {
         VV <- xMatHatOmegaInv %*% PHOmega %*%
            PH %*% t( xMatHatOmegaInv )
         VV <- rbind( cbind( VV, matrix( 0, nrow( VV ), nrow( R.restr ) ) ),
            matrix( 0, nrow( R.restr ), nrow( VV ) + nrow( R.restr ) ) )
         bcov <- ( Winv %*% VV %*% Winv )[ 1:ncol(xMatAll), 1:ncol(xMatAll) ]
                  # coefficient covariance matrix
      }
    }
    if(method3sls=="EViews") {
       bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = solvetol )  # final step coefficient covariance matrix
    }
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## FIML estimation
  if( method == "FIML" ) {
    fimlResult <- .systemfitFiml( systemfitCall = call, nObsEq = nObsEq,
      nCoefEq = nExogLiEq, yVec = yVecAll, xMat = xMatAll, xEq = xMatEq, methodRCov = methodRCov,
      centerResiduals = centerResiduals, solvetol = solvetol )
    #print( fimlResult )
    coef <- fimlResult$coef
    bcov <- fimlResult$coefCov
    resids <- fimlResult$resids
  }

  ## for all estimation methods
  fitted <- xMatAll %*% coef                              # fitted endogenous values
  bt     <- NULL
  btcov  <- NULL
  if(!is.null(TX)) {
    bt <- coef
    coef  <- TX %*% bt
    btcov <- bcov
    bcov  <- TX %*% btcov %*% t(TX)
  }


  ## equation wise results
  for(i in 1:nEq) {
    residi[[i]] <- resids[ ( 1 + sum(nObsEq[1:i]) -nObsEq[i] ):( sum(nObsEq[1:i]) ) ]
    coefEqI <- drop( coef[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))] )
              # estimated coefficients of equation i
    bcovi  <- bcov[(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i])),(1+sum(nExogEq[1:i])-nExogEq[i]):(sum(nExogEq[1:i]))]
              # covariance matrix of estimated coefficients of equation i

    # set names
    names( coefEqI ) <- colnames( xMatEq[[i]] )
    colnames( bcovi ) <- colnames( xMatEq[[i]] )
    rownames( bcovi ) <- colnames( xMatEq[[i]] )

    ssr    <- sum(residi[[i]]^2)                         # sum of squared residuals
    sigma  <- sqrt( ssr / df[i] ) # estimated standand deviation of residuals
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
    if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
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
    resulti$nExogAll     <- nExogAll        # number of exogenous variables/coefficients
    resulti$nExogLiAll   <- nExogLiAll      # number of linear independent coefficients
    resulti$df           <- df[i]           # degrees of freedom of residuals
    resulti$dfSys        <- nObsAll- nExogLiAll
       # degrees of freedom of residuals of the whole system
    resulti$coef         <- coefEqI         # estimated coefficients
    resulti$bcov         <- bcovi           # covariance matrix of estimated coefficients
    resulti$yVec         <- yVecEq[[i]]     # vector of endogenous variables
    resulti$xMat         <- xMatEq[[i]]     # matrix of regressors
    resulti$data         <- datai           # data frame of this equation (incl. instruments)
    resulti$fitted       <- fittedi         # fitted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$sigma        <- sigma           # estimated standard error of the residuals
    resulti$r2           <- r2              # R-sqared value
    resulti$adjr2        <- adjr2           # adjusted R-squared value
    if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
      resulti$inst         <- instl[[i]]
      resulti$hMat         <- hMatEq[[i]]          # matrix of instrumental variables
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
  if( method %in% c(  "SUR", "WSUR", "3SLS", "W3SLS" ) ) {
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
      rtOmega <- .calcXtOmegaInv( xMat = resids, sigma = rcov, nObsEq = nObsEq,
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

  coef           <- drop(coef)
  names(coef)    <- xnames
  colnames( bcov ) <- xnames
  rownames( bcov ) <- xnames
  colnames( rcov ) <- eqnlabels
  rownames( rcov ) <- eqnlabels

  ## build the "return" structure for the whole system
  results$method  <- method
  results$eqnLabels <- eqnlabels
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
  results$coef    <- coef           # all estimated coefficients
  results$bt      <- bt             # transformed vector of estimated coefficients
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
  if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ){
    results$hMat    <- hMatAll            # matrix of all (diagonally stacked) instr. variables
    results$xHat    <- xMatHatAll           # matrix of "fitted" regressors
  }
  if( method %in% c( "SUR", "WSUR", "3SLS", "W3SLS" ) ){
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
  results$single.eq.sigma <- single.eq.sigma
  results$solvetol        <- solvetol
  class(results)  <- "systemfit"

  results
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
   print( x$coef )
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
   print( x$coef )
   invisible( x )
}


## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, data=object$data,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   attach(data); on.exit( detach( data ) )

   predicted <- data.frame( obs=seq( nrow( data ) ) )
   colnames( predicted ) <- as.character( 1:ncol( predicted ) )
   nObsEq  <- numeric( object$nEq )
   eqns    <- list()
   xMatEq  <- list()               # regressors equation-wise
   xMatAll <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
   for(i in 1:object$nEq )  {
      eqns[[i]] <- object$eq[[i]]$formula
      xMatEq[[i]] <-  model.matrix( eqns[[i]] )
      xMatAll      <-  rbind( cbind( xMatAll, matrix( 0, nrow( xMatAll ), ncol( xMatEq[[i]] ))),
                       cbind( matrix( 0, nrow( xMatEq[[i]] ), ncol( xMatAll )), xMatEq[[i]]))
      nObsEq[i]   <-  nrow( xMatEq[[i]] )
   }
   yVecAll <- xMatAll %*% object$coef
   if( object$method %in% c( "SUR", "WSUR", "3SLS", "W3SLS" ) ){
      if( se.fit | interval == "confidence" ) {
         ycovc <- xMatAll %*% object$bcov %*% t(xMatAll)
      }
      if( se.pred | interval == "prediction" ) {
         ycovp <- xMatAll %*% object$bcov %*% t(xMatAll) + object$rcov %x% diag(1,nObsEq[1],nObsEq[1])
      }
   }
   for(i in 1:object$nEq) {
      # fitted values
      Yi <- yVecAll[(1+sum(nObsEq[1:i])-nObsEq[i]):sum(nObsEq[1:i]),]
      predicted <- cbind( predicted, Yi )
      names( predicted )[ length( predicted ) ] <- paste( "eq", as.character(i),
                                                          ".pred", sep="" )
      # calculate variance covariance matrices
      if( se.fit | interval == "confidence" ) {
         if( object$method %in% c( "SUR", "WSUR", "3SLS", "W3SLS" ) ){
            ycovci <- ycovc[ ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ),
                             ( 1 + sum( nObsEq[1:i] ) - nObsEq[i] ) : sum( nObsEq[1:i] ) ]
         } else {
            ycovci <- xMatEq[[i]] %*% object$eq[[i]]$bcov %*% t(xMatEq[[i]])
         }
      }
      if( se.pred | interval == "prediction" ) {
         if( object$method %in% c( "SUR", "WSUR", "3SLS", "W3SLS" ) ){
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
         if( probDfSys ) {
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
         if( probDfSys ) {
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
   xMat <-  model.matrix( object$formula )
   predicted <- drop( xMat %*% object$coef )
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

## return the variance covariance matrix of the coefficients
confint.systemfit <- function( object, parm = NULL, level = 0.95,
      probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- matrix( NA, length( object$coef ), 2,
            dimnames = list( names( object$coef ), pct ) )
   j <- 1
   for( i in 1:object$nEq ) {
      object$eq[[i]]$dfSys <- object$df
      ci[ j:(j+object$eq[[ i ]]$nExog-1), ] <- confint( object$eq[[ i ]],
         probDfSys = probDfSys )
      j <- j + object$eq[[ i ]]$nExog
   }
   class( ci ) <- "confint.systemfit"
   ci
}

## return the variance covariance matrix of the coefficients of a single equation
confint.systemfit.equation <- function( object, parm = NULL, level = 0.95,
   probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   a <- ( 1 - level ) / 2
   a <- c( a, 1 - a )
   pct <- paste( round( 100 * a, 1 ), "%" )
   ci <- matrix( NA, length( object$coef ), 2,
            dimnames = list( names( object$coef ), pct ) )
   if( probDfSys ) {
      fac <- qt( a, object$dfSys )
   } else {
      fac <- qt( a, object$df )
   }
   coef <- summary( object )$coefficients
   ci[] <- coef[ , 1 ] + coef[ , 2 ] %o% fac
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
   fitted <- matrix( NA, length( object$eq[[1]]$fitted ), object$nEq )
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




