## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, newdata=object$data,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   predicted <- data.frame( obs=seq( nrow( newdata ) ) )
   colnames( predicted ) <- as.character( 1:ncol( predicted ) )

   # get equations
   eqns <- list()
   for(i in 1:object$nEq )  {
      eqns[[i]] <- formula( object$eq[[i]]$terms )
   }

   # prepare data
   preparedData <- .prepareData.systemfit( data = newdata, eqns = eqns )
   # number of observations in each equation
   nObsEq  <- preparedData$nObsEq
   # list for matrices of regressors in each equation
   xMatEq <- preparedData$xMatEq
   # stacked matrices of all regressors
   xMatAll <- preparedData$xMatAll

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
                                 object$eq[[i]]$sigma^2
         }
      }
      # standard errors of fitted values
      if( se.fit ) {
         if(nrow(newdata)==1) {
            predicted <- cbind( predicted, sqrt( ycovci ) )
         } else {
            predicted <- cbind( predicted, sqrt( diag( ycovci ) ) )
         }
         names( predicted )[ length( predicted ) ] <-
            paste( "eq", as.character(i), ".se.fit", sep="" )
      }
      # standard errors of prediction
      if( se.pred ) {
         if(nrow(newdata)==1) {
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
         if( nrow(newdata)==1 ) {
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
         if(nrow(newdata)==1) {
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
predict.systemfit.equation <- function( object, newdata=object$data,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   xMat <-  model.matrix( formula( object$terms ), data = newdata )
   # fitted values
   predicted <- data.frame( fit = drop( xMat %*% object$coef ) )

   # calculate variance covariance matrices
   if( se.fit | interval == "confidence" ) {
      yCovConf <- xMat %*% object$bcov %*% t( xMat )
   }
   if( se.pred | interval == "prediction" ) {
      yCovPred <- xMat %*% object$bcov %*% t( xMat ) + object$sigma^2
   }
   # standard errors of fitted values
   if( se.fit ) {
      if(nrow(newdata)==1) {
         predicted[[ "se.fit" ]] <- sqrt( yCovConf )
      } else {
         predicted[[ "se.fit" ]] <- sqrt( diag( yCovConf ) )
      }
   }
   # standard errors of prediction
   if( se.pred ) {
      if(nrow(newdata)==1) {
         predicted[[ "se.pred" ]] <- sqrt( yCovPred )
      } else {
         predicted[[ "se.pred" ]] <- sqrt( diag( yCovPred ) )
      }
   }

   # confidence intervals
   if( interval == "confidence" ) {
      if( probDfSys ) {
         tval   <- qt( 1 - ( 1- level )/2, object$dfSys )
      } else {
         tval   <- qt( 1 - ( 1- level )/2, object$df )
      }
      if( nrow(newdata)==1 ) {
         stdErConf <- sqrt( yCovConf )
      } else {
         stdErConf <- sqrt( diag( yCovConf ) )
      }
      predicted[[ "lwr" ]] <- predicted$fit - ( tval * stdErConf )
      predicted[[ "upr" ]] <- predicted$fit + ( tval * stdErConf )
   }
   # prediction intervals
   if( interval == "prediction" ) {
      if( probDfSys ) {
         tval   <- qt( 1 - ( 1- level )/2, object$dfSys )
      } else {
         tval   <- qt( 1 - ( 1- level )/2, object$df )
      }
      if(nrow(newdata)==1) {
         stdErPred <- sqrt( yCovPred )
      } else {
         stdErPred <- sqrt( diag( yCovPred ) )
      }
      predicted[[ "lwr" ]] <- predicted$fit - ( tval * stdErPred )
      predicted[[ "upr" ]] <- predicted$fit + ( tval * stdErPred )
   }

   return( predicted )
}

