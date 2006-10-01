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

   if( is.null( names( eqns ) ) ) {
      eqnLabels <- paste( "eq", c( 1:nEq ), sep = "" )
   } else {
      eqnLabels <- names(eqns)
      if( sum( regexpr( " |_", eqnLabels ) != -1 ) > 0 ) {
         stop( "equation labels may not contain blanks (' ') or underscores ('_')" )
      }
   }

   results$call <- match.call() # get the original call

   # prepare data
   preparedData <- .prepareData.systemfit( data = data, eqns = eqns,
      inst = inst, TX = TX, control = control, eqnLabels = eqnLabels )
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
   nCoefEq <- preparedData$nCoefEq
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

   nObsAll  <- sum( nObsEq )  # total number of observations of all equations
   nCoefAll <- sum( nCoefEq ) # total number of exogenous variables/(unrestricted) coefficients in all equations
   nCoefLiAll <- nCoefAll     # total number of linear independent coefficients in all equations
   nCoefLiEq  <- nCoefEq      # total number of linear independent coefficients in each equation
   if(!is.null(TX)) {
      nCoefLiAll <- nCoefLiAll - ( nrow( TX ) - ncol( TX ) )
      for(i in 1:nEq) {
         nCoefLiEq[i] <- ncol(xMatAll)
         for(j in 1: ncol(xMatAll) ) {
            if(sum(xMatAll[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])),j]^2)==0) {
               nCoefLiEq[i] <- nCoefLiEq[i]-1
            }
         }
      }
   }
   if(!is.null(R.restr)) {
      nCoefLiAll  <- nCoefLiAll - nrow(R.restr)
      if(is.null(TX)) {
         for(j in 1:nrow(R.restr)) {
            for(i in 1:nEq) {  # search for restrictions that are NOT cross-equation
               if( sum( R.restr[ j, (1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))]^2) ==
                   sum(R.restr[j,]^2)) {
                  nCoefLiEq[i] <- nCoefLiEq[i]-1
               }
            }
         }
      }
    }
    df <- nObsEq - nCoefLiEq    # degress of freedom of each equation

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
         nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )               # residual covariance matrix
      bcov <- .calcGLS( xMat = xMatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )
                    # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nCoefLiAll, methodRCov = control$methodRCov )
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
         nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
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
         nCoefEq = nCoefLiEq, xEq = xMatEq, centered = control$centerResiduals,
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
         nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
         solvetol = control$solvetol )
      bcov <- .calcGLS( xMat = xMatHatAll, R.restr = R.restr, q.restr = q.restr,
         sigma = rcov, nObsEq = nObsEq, solvetol = control$solvetol )  # coefficient covariance matrix
    } else {
      s2 <- .calcSigma2( resids, nObs = nObsAll, nCoef = nCoefLiAll, methodRCov = control$methodRCov )
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
         nCoefEq = nCoefLiEq, xEq = xMatEq, diag = TRUE, centered = control$centerResiduals,
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
         nCoefEq = nCoefLiEq, xEq = xMatEq, centered = control$centerResiduals, solvetol = control$solvetol )
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
      nCoefEq = nCoefLiEq, yVec = yVecAll, xMat = xMatAll, xEq = xMatEq, methodRCov = control$methodRCov,
      centerResiduals = control$centerResiduals, solvetol = control$solvetol )
    #print( fimlResult )
    coef <- fimlResult$coefficients
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
    coefEqI <- drop( coef[(1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))] )
              # estimated coefficients of equation i
    bcovi  <- bcov[(1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i])),(1+sum(nCoefEq[1:i])-nCoefEq[i]):(sum(nCoefEq[1:i]))]
              # covariance matrix of estimated coefficients of equation i

    # set names
    names( coefEqI ) <- colnames( xMatEq[[i]] )
    colnames( bcovi ) <- colnames( xMatEq[[i]] )
    rownames( bcovi ) <- colnames( xMatEq[[i]] )

    ssr    <- sum(residi[[i]]^2)                         # sum of squared residuals
    sigma  <- sqrt( ssr / df[i] ) # estimated standand deviation of residuals
    fitted.values.i <- fitted.values[(1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i]))]

    ## build the "return" structure for the equations
    resulti$method       <- method
    resulti$i            <- i               # equation number
    resulti$eqnLabel     <- eqnLabels[[i]]
    resulti$terms        <- termsEq[[ i ]]
    resulti$nObs         <- nObsEq[i]       # number of observations
    resulti$nCoef        <- nCoefEq[i]      # number of exogenous variables/coefficients
    resulti$nCoefLi      <- nCoefLiEq[i]    # number of linear independent coefficients
    resulti$nCoef.sys    <- nCoefAll        # number of exogenous variables/coefficients
    resulti$nCoefLi.sys  <- nCoefLiAll      # number of linear independent coefficients
    resulti$df.residual  <- df[i]           # degrees of freedom of residuals
    resulti$df.residual.sys  <- nObsAll- nCoefLiAll
       # degrees of freedom of residuals of the whole system
    resulti$coefficients <- coefEqI         # estimated coefficients
    resulti$bcov         <- bcovi           # covariance matrix of estimated coefficients
    if( control$returnResponse ){
      resulti$response     <- yVecEq[[i]]     # vector of endogenous variables
    }
    if( control$returnModelMatrix ){
      resulti$modelMatrix  <- xMatEq[[i]]     # matrix of regressors
    }
    if( control$returnModelFrame ){
      resulti$modelFrame   <- evalModelFrameEq[[ i ]] # model frame of this equation
    }
    resulti$fitted.values <- fitted.values.i # fitted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr             # sum of squared errors/residuals
    resulti$sigma        <- sigma           # estimated standard error of the residuals
    if( method %in% c( "2SLS", "W2SLS", "3SLS", "W3SLS" ) ) {
      resulti$inst         <- instEq[[i]]
      if(  control$returnInstMatrix ) {
         resulti$instMatrix   <- hMatEq[[i]]  # matrix of instrumental variables
      }
    }
    class(resulti)        <- "systemfit.equation"
    results$eq[[i]]      <- resulti
  }

  ## results of the total system
  if( method %in% c(  "SUR", "WSUR", "3SLS", "W3SLS" ) ) {
    rcovest <- rcov                   # residual covariance matrix used for estimation
  }
  rcov <- .calcRCov( resids, methodRCov = control$methodRCov, nObsEq = nObsEq,
      nCoefEq = nCoefLiEq, xEq = xMatEq, centered = control$centerResiduals, solvetol = control$solvetol )

  coef           <- drop(coef)
  names(coef)    <- xnames
  colnames( bcov ) <- xnames
  rownames( bcov ) <- xnames
  colnames( rcov ) <- eqnLabels
  rownames( rcov ) <- eqnLabels

  ## build the "return" structure for the whole system
  results$method  <- method
  results$nObs    <- nObsAll        # total number of observations of all equations
  results$nCoef   <- nCoefAll      # total number of exogenous variables/coefficients in all equations
  results$nCoefLi <- nCoefLiAll
     # total number of linear independent coefficients of all equations
  results$df.residual <- nObsAll - nCoefLiAll
     # degrees of freedom of the whole system
  results$coefficients <- coef           # all estimated coefficients
  results$bt      <- bt             # transformed vector of estimated coefficients
  results$bcov    <- bcov           # coefficients covariance matrix
  results$btcov   <- btcov          # covariance matrix for transformed coeff. vector
  results$rcov    <- rcov           # residual covarance matrix
  results$iter    <- iter           # residual correlation matrix
  if( method %in% c( "SUR", "WSUR", "3SLS", "W3SLS" ) ){
    results$rcovest <- rcovest      # residual covarance matrix used for estimation
  }
  results$R.restr <- R.restr
  results$q.restr <- q.restr
  results$TX      <- TX
  results$control <- control
  class(results)  <- "systemfit"

  results
}
