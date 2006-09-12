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
                        control = systemfit.control( ... ),
                        ... )
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

   if( is.null( control$single.eq.sigma ) ) {
      control$single.eq.sigma <- ( is.null( R.restr ) & is.null( TX ) )
   }

  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  resulti <- list()               # results of the ith equation
  residi  <- list()               # residuals equation wise
  iter    <- NULL                 # number of iterations
  nEq     <- length( eqns )       # number of equations
  ssr     <- numeric( nEq ) # sum of squared residuals of each equation
  sigma   <- numeric( nEq ) # estimated sigma (std. dev. of residuals) of each equation
  r2      <- numeric( nEq ) # R-squared value
  adjr2   <- numeric( nEq ) # adjusted R-squared value

   if( is.null( names( eqns ) ) ) {
      eqnlabels <- paste( "eq", c( 1:nEq ), sep = "" )
   } else {
      eqnlabels <- names(eqns)
   }

   results$call <- match.call() # get the original call
   callNoDots <- match.call( expand.dots = FALSE ) #-"- without ...-expansion
   if( "data" %in% names( callNoDots ) ) {
      results$data.name <- callNoDots$data
   } else {
      results$data.name <- "unknown"
   }

   # prepare data
   preparedData <- .prepareData.systemfit( data = data, eqns = eqns,
      inst = inst, TX = TX, control = control )
   # list of terms objects of each equation
   termsEq <- preparedData$termsEq
   # list of evaluated model frames of each equation 
   evalModelFrameEq <- preparedData$evalModelFrameEq
   # list for vectors of endogenous variables in each equation
   yVecEq <- preparedData$yVecEq
   # stacked endogenous variables of all equations
   yVecAll <- preparedData$yVecAll
   # list for matrices of regressors in each equation
   xMatEq <- preparedData$xMatEq
   # stacked matrices of all regressors (multiplied with TX)
   xMatAll <- preparedData$xMatAll
   # number of observations in each equation
   nObsEq <- preparedData$nObsEq
   # number of exogenous variables /(unrestricted) coefficients in each equation
   nExogEq <- preparedData$nExogEq
   # names of regressors
   xnames <- preparedData$xnames
   if( !is.null( inst ) ) {
      # list of formulas for instruments of each equation
      instEq <- preparedData$instEq
      # list of terms objects of instruments of each equation
      termsInst  <- preparedData$termsInst
      # model frame of instruments
      evalModelFrameInst  <- preparedData$evalModelFrameInst
      # stacked matrices of all instruments
      hMatAll  <- preparedData$hMatAll      # stacked matrices of all instruments
      # list for matrices of instruments in each equation
      hMatEq  <- preparedData$hMatEq
      # fitted values of all regressors
      xMatHatAll <- preparedData$xMatHatAll
   }

   if( is.null( control$saveMemory ) ) {
      control$saveMemory <- sum( nObsEq ) > 1000
   } else if( sum( nObsEq ) > 1000 && !control$saveMemory ) {
      warning( paste( "You have more than 1000 observations.",
         "Setting control variable 'saveMemory' to TRUE speeds up",
         "the estimation. Estimation of larger data sets might even",
         "require this setting." ) )
   }

   nObsAll  <- sum( nObsEq )  # total number of observations of all equations
   nExogAll <- sum( nExogEq ) # total number of exogenous variables/(unrestricted) coefficients in all equations
   nExogLiAll <- nExogAll     # total number of linear independent coefficients in all equations
   nExogLiEq  <- nExogEq      # total number of linear independent coefficients in each equation
   if(!is.null(TX)) {
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
      coef <- solve( crossprod( xMatAll ), crossprod( xMatAll, yVecAll ), tol=control$solvetol )
               # estimated coefficients
    } else {
      W <- rbind( cbind( t(xMatAll) %*% xMatAll, t(R.restr) ),
                  cbind( R.restr, matrix( 0, nrow(R.restr), nrow(R.restr) )))
      V <- rbind( t(xMatAll) %*% yVecAll , q.restr )
      coef <- ( solve( W, tol=control$solvetol ) %*% V )[1:ncol(xMatAll)]
    }
  }

  ## only for OLS estimation
  if(method=="OLS") {
    resids <- yVecAll - xMatAll %*% coef                                        # residuals
    if(control$single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )               # residual covariance matrix
      bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )
                    # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nExogLiAll, methodRCov = control$methodRCov )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( xMatAll ), tol=control$solvetol )
                          # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=control$solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coefficient covariance matrix
      }
    }
  }

  ## only for WLS estimation
  if( method %in% c( "WLS", "WSUR" ) ) {
    bl    <- coef   # coefficients of previous step
    bdif  <- coef   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>control$tol & iter < control$maxiter^( method == "WLS" ) ) {
      iter  <- iter+1
      bl    <- coef                # coefficients of previous step
      resids <- yVecAll - xMatAll %*% coef     # residuals
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )
      coef  <- .calcGLS( xMat = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol ) # coefficients
      bdif <- coef-bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )
       # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for SUR estimation
  if( method %in% c( "SUR", "WSUR" ) ) {
    bl    <- coef    # coefficients of previous step
    bdif  <- coef    # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>control$tol & iter < control$maxiter) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = control$centerResiduals,
         solvetol = control$solvetol )
      coef <- .calcGLS( xMat = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )     # coefficients
      bdif <- coef-bl  # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )
            # final step coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 2SLS, W2SLS and 3SLS estimation
  if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
    if(is.null(R.restr)) {
      coef <- solve( crossprod( xMatHatAll ), crossprod( xMatHatAll, yVecAll ), tol=control$solvetol )
         # 2nd stage coefficients
    } else {
      W <- rbind( cbind( crossprod(xMatHatAll), t(R.restr) ),
                  cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
      V <- rbind( t(xMatHatAll) %*% yVecAll , q.restr )
      coef <- ( solve( W, tol=control$solvetol ) %*% V )[1:ncol(xMatAll)] # restricted coefficients
    }
    b2 <- coef
  }

  ## only for 2SLS estimation
  if(method=="2SLS") {
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    if(control$single.eq.sigma) {
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )
      bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nExogLiAll, methodRCov = control$methodRCov )
                           # sigma squared
      if(is.null(R.restr)) {
        bcov   <- s2 * solve( crossprod( xMatHatAll ), tol=control$solvetol )
                  # coefficient covariance matrix
      } else {
        bcov   <- s2 * solve( W, tol=control$solvetol )[1:ncol(xMatAll),1:ncol(xMatAll)]
                    # coeff. covariance matrix
      }
    }
  }

  ## only for W2SLS estimation
  if( method %in% c( "W2SLS", "W3SLS" ) ) {
    bl     <- coef   # coefficients of previous step
    bdif   <- coef   # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>control$tol & iter < control$maxiter^( method == "W2LS" ) ) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )
      coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )          # (unrestr.) coeffic.
      bdif <- coef - bl # difference of coefficients between this and previous step
    }
    bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
       sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # coefficient covariance matrix
    resids <- yVecAll - xMatAll %*% coef                        # residuals
    for(i in 1:nEq) residi[[i]] <- resids[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
  }

  ## only for 3SLS estimation
  if( method %in% c( "3SLS", "W3SLS" ) ) {
    bl     <- coef  # coefficients of previous step
    bdif   <- coef  # difference of coefficients between this and previous step
    iter  <- 0
    while((sum(bdif^2)/sum(bl^2))^0.5>control$tol & iter < control$maxiter) {
      iter  <- iter+1
      bl    <- coef                           # coefficients of previous step
      resids <- yVecAll-xMatAll%*%coef                     # residuals
      rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
         nCoefEq = nExogLiEq, xEq = xMatEq, centered = control$centerResiduals, solvetol = control$solvetol )
      if(control$method3sls=="GLS") {
         coef <- .calcGLS( xMat = xMatHatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # (unrestr.) coeffic.
      }
      if(control$method3sls=="IV") {
         coef <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll, yVec = yVecAll, R.restr = R.restr, q.restr = q.restr,
            sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )   # (unrestr.) coeffic.
      }
      if(control$method3sls=="GMM") {
        HtOmega <- .calcXtOmegaInv( xMat = hMatAll, sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
        if(is.null(R.restr)) {
          coef <- solve(t(xMatAll) %*% hMatAll %*% solve( HtOmega %*%
                 hMatAll, tol=control$solvetol) %*% t(hMatAll) %*% xMatAll, tol=control$solvetol) %*% t(xMatAll) %*% hMatAll %*%
                 solve( HtOmega %*%
                 hMatAll, tol=control$solvetol) %*% t(hMatAll) %*% yVecAll  #(unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(xMatAll) %*% hMatAll %*% solve( HtOmega
                              %*% hMatAll, tol=control$solvetol) %*% t(hMatAll) %*% xMatAll, t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( t(xMatAll) %*% hMatAll %*% solve( HtOmega
                      %*% hMatAll, tol=control$solvetol) %*% t(hMatAll) %*% yVecAll , q.restr )
          Winv <- solve( W, tol=control$solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(control$method3sls=="Schmidt") {
        xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov, nObsEq = nObsEq,
           solvetol = control$solvetol )
        if(is.null(R.restr)) {
          coef <- solve( t(xMatHatAll) %*% t( xMatHatOmegaInv ), tol=control$solvetol) %*% ( xMatHatOmegaInv
                      %*% hMatAll %*% solve( crossprod( hMatAll ), tol=control$solvetol ) %*% crossprod(hMatAll, yVecAll) )
                           # (unrestr.) coeffic.
        } else {
          W <- rbind( cbind( t(xMatHatAll) %*% t( xMatHatOmegaInv ), t(R.restr) ),
                      cbind( R.restr, matrix(0, nrow(R.restr), nrow(R.restr))))
          V <- rbind( xMatHatOmegaInv %*% hMatAll %*% solve( crossprod( hMatAll ), tol=control$solvetol ) %*%
                      crossprod( hMatAll, yVecAll ), q.restr )
          Winv <- solve( W, tol=control$solvetol )
          coef <- ( Winv %*% V )[1:ncol(xMatAll)]     # restricted coefficients
        }
      }
      if(control$method3sls=="EViews") {
         coef <- b2 + .calcGLS( xMat = xMatHatAll, yVec = ( yVecAll -  xMatAll %*% b2 ),
            R.restr = R.restr, q.restr = q.restr, sigma = rcov,
            nObsEq = nObsEq, solvetol = control$solvetol )  # (unrestr.) coeffic.
      }
      bdif <- coef - bl # difference of coefficients between this and previous step
    }
    if(control$method3sls=="GLS") {
       bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # coefficient covariance matrix
    }
    if(control$method3sls=="IV") {
       bcov <- .calcGLS( xMat = xMatHatAll, xMat2 = xMatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )
    }
    if(control$method3sls=="GMM") {
      if(is.null(R.restr)) {
        bcov <- solve( t(xMatAll) %*% hMatAll %*% solve( HtOmega %*% hMatAll, tol=control$solvetol ) %*%
           t(hMatAll) %*% xMatAll, tol=control$solvetol )
                # final step coefficient covariance matrix
      } else {
        bcov   <- Winv[1:ncol(xMatAll),1:ncol(xMatAll)] # coefficient covariance matrix
      }
    }
    if(control$method3sls=="Schmidt") {
      xMatHatOmegaInv <- .calcXtOmegaInv( xMat = xMatHatAll, sigma = rcov, nObsEq = nObsEq,
         solvetol = control$solvetol )
      PH <- hMatAll %*%  solve( t(hMatAll) %*% hMatAll, tol=control$solvetol ) %*% t(hMatAll)
      PHOmega <- .calcXtOmegaInv( xMat = t( PH ), sigma = rcov, nObsEq = nObsEq,
           invertSigma = FALSE )
      if(is.null(R.restr)) {
         bcov <- solve( xMatHatOmegaInv %*% xMatHatAll, tol=control$solvetol ) %*%
            xMatHatOmegaInv %*% PHOmega %*%
            PH %*% t( xMatHatOmegaInv ) %*% solve( xMatHatOmegaInv %*% xMatHatAll, tol=control$solvetol )
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
    if(control$method3sls=="EViews") {
       bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
          sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # final step coefficient covariance matrix
    }
    resids <- yVecAll - xMatAll %*% coef                        # residuals
  }

  ## FIML estimation
  if( method == "FIML" ) {
    fimlResult <- .systemfitFiml( systemfitCall = results$call, nObsEq = nObsEq,
      nCoefEq = nExogLiEq, yVec = yVecAll, xMat = xMatAll, xEq = xMatEq, methodRCov = control$methodRCov,
      centerResiduals = control$centerResiduals, solvetol = control$solvetol )
    #print( fimlResult )
    coef <- fimlResult$coef
    bcov <- fimlResult$coefCov
    resids <- fimlResult$resids
  }

  ## for all estimation methods
  fitted.values <- xMatAll %*% coef   # fitted endogenous values
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
    fitted.values.i <- fitted.values[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]
    resp <- model.extract( evalModelFrameEq[[ i ]], "response" )
    datai <- data.frame( cbind( resp, ( model.matrix( termsEq[[ i ]],
      evalModelFrameEq[[ i ]] ) )[ , -1 ] ) )
    names( datai )[1] <- as.character( terms( eqns[[ i ]] ) )[2]
    rm( resp )
    if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
      datai <- cbind( datai, as.data.frame( model.matrix( termsInst[[ i ]],
         evalModelFrameInst[[ i ]] )[ , -1 ] ) )
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
    resulti$terms        <- termsEq[[ i ]]
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
    resulti$fitted.values <- fitted.values.i # fitted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$sigma        <- sigma           # estimated standard error of the residuals
    resulti$r2           <- r2              # R-sqared value
    resulti$adjr2        <- adjr2           # adjusted R-squared value
    if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
      resulti$inst         <- instEq[[i]]
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
  rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
      nCoefEq = nExogLiEq, xEq = xMatEq, centered = control$centerResiduals, solvetol = control$solvetol )
  drcov <- det(rcov, tol=control$solvetol)
  if( !control$saveMemory ) {
#       # original formula from McElroy (1977)
#       mcelr2 <- 1 - ( t(resids) %*% ( solve(rcov, tol=control$solvetol) %x%
#                 diag(1, nObsEq[1],nObsEq[1])) %*% resids ) /
#                 ( t(yVecAll) %*% ( solve(rcov, tol=control$solvetol ) %x%
#                 ( diag(1,nObsEq[1],nObsEq[1] ) - rep(1,nObsEq[1]) %*%
#                 t(rep(1,nObsEq[1])) / nObsEq[1] )) %*% yVecAll )   # McElroy's (1977a) R2
      # first formula from Greene (2003, p. 345) (numerator modified to save memory)
      rtOmega <- .calcXtOmegaInv( xMat = resids, sigma = rcov, nObsEq = nObsEq,
         solvetol = control$solvetol )
      yCov <- .calcRCov( yVecAll, methodRCov = "noDfCor", nObsEq = nObsEq, centered = TRUE,
         solvetol = control$solvetol )
      residCovInv <- solve( rcov, tol = control$solvetol )
      denominator <- 0
      for( i in 1:nEq ) {
         for( j in 1:nEq ) {
            denominator <- denominator + residCovInv[ i, j ] * yCov[ i, j ] * nObsEq[1]
         }
      }
      mcelr2 <- 1 - ( rtOmega %*% resids ) / denominator
#       # second formula from Greene (2003, p. 345)
#        yCov <- sum(diag(.calcRCov( yVecAll, methodRCov = "noDfCor", nObsEq = nObsEq, centered = TRUE,
#           solvetol = control$solvetol )))
#        yCov <- drop( t(yVecAll-mean(yVecAll)) %*% (yVecAll-mean(yVecAll)) / sum(nObsEq) )
#       yCov <- .calcRCov( yVecAll, methodRCov = "geomean", nObsEq = nObsEq, nCoefEq=rep(1,nEq),
#          centered = TRUE, solvetol = control$solvetol )
#        mcelr2 <- 1 - nEq / sum( diag( solve( rcov, tol = control$solvetol ) * yCov ) )
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
  results$control <- control
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
         if(x$iter<x$control$maxiter) {
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
   print( formula( x$terms ) )
   if(!is.null(x$inst)) {
      cat("Instruments: ")
      print(x$inst)
   }

   cat("\nCoefficients:")
   print( x$coef )
   invisible( x )
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




