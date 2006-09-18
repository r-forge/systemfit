## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit <- function( object, newdata = NULL,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   for(i in 1:object$nEq )  {
      predicted.i <- predict( object$eq[[ i ]], newdata = newdata,
         se.fit = se.fit, se.pred = se.pred, interval = interval,
         level = level, probDfSys = probDfSys )
      names( predicted.i ) <- paste( "eq", i, ".",
         names( predicted.i ), sep = "" )
      if( i == 1 ) {
         predicted <- predicted.i
      } else {
         predicted <- cbind( predicted, predicted.i )
      }
   }
   names( predicted ) <- sub( "([0-9])\.fit", "\\1.pred",
         names( predicted ) )

   return( predicted )
}

## calculate predicted values, its standard errors and the prediction intervals
predict.systemfit.equation <- function( object, newdata = NULL,
                               se.fit=FALSE, se.pred=FALSE,
                               interval="none", level=0.95,
                               probDfSys = NULL, ... ) {

   if( is.null( probDfSys ) ) {
      probDfSys <- object$nExogAll != object$nExogLiAll
         # TRUE if there are restrictions imposed
   }

   if( is.null( newdata ) ) {
      xMat <-  model.matrix( object )
   } else {
      xMat <-  model.matrix( formula( object$terms ), data = newdata )
   }

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
      if( length( yCovConf ) == 1 ) {
         predicted[[ "se.fit" ]] <- sqrt( yCovConf )
      } else {
         predicted[[ "se.fit" ]] <- sqrt( diag( yCovConf ) )
      }
   }
   # standard errors of prediction
   if( se.pred ) {
      if( length( yCovPred ) == 1 ) {
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
      if(  length( yCovConf ) == 1 ) {
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
      if( length( yCovPred ) == 1 ) {
         stdErPred <- sqrt( yCovPred )
      } else {
         stdErPred <- sqrt( diag( yCovPred ) )
      }
      predicted[[ "lwr" ]] <- predicted$fit - ( tval * stdErPred )
      predicted[[ "upr" ]] <- predicted$fit + ( tval * stdErPred )
   }

   return( predicted )
}

