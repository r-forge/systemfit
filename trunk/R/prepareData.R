.prepareData.systemfit <- function( data, eqns, inst = NULL, TX = NULL, control = NULL,
      eqnLabels )
{
   # list for results
   result <- list()

   nEq     <- length( eqns )       # number of equations
   yVecEq  <- list()               # list for vectors of endogenous variables in each equation
   yVecAll <- matrix( 0, 0, 1 )    # stacked endogenous variables of all equations
   xMatEq  <- list()               # list for matrices of regressors in each equation
   xMatAll <- matrix( 0, 0, 0 )    # stacked matrices of all regressors (unrestricted)
   nObsEq  <- numeric( nEq ) # number of observations in each equation
   nExogEq <- numeric( nEq ) # number of exogenous variables /(unrestricted) coefficients
                                     # in each equation
   instEq  <- list()         # list of the instruments for each equation
   xnames  <- NULL           # names of regressors

   callNoDots <- match.call( expand.dots = FALSE ) #-"- without ...-expansion

   # model frame (without formula)
   modelFrame <- callNoDots[ c( 1, match( "data", names( callNoDots ), 0 ) ) ]
   modelFrame[[1]] <- as.name( "model.frame" )
   # terms and model frames for the individual equations
   termsEq <- list()
   modelFrameEq <- list()
   evalModelFrameEq <- list()
   # prepare data for individual equations
   for(i in 1:nEq ) {
      modelFrameEq[[ i ]] <- modelFrame
      modelFrameEq[[ i ]]$formula <- eqns[[ i ]]
      evalModelFrameEq[[ i ]] <- eval( modelFrameEq[[ i ]], parent.frame() )
      termsEq[[ i ]] <- attr( evalModelFrameEq[[ i ]], "terms" )
      weights <- model.extract( evalModelFrameEq[[ i ]], "weights" )
      yVecEq[[i]] <- model.extract( evalModelFrameEq[[ i ]], "response" )
      xMatEq[[i]] <- model.matrix( termsEq[[ i ]], evalModelFrameEq[[ i ]] )
      yVecAll <- c(yVecAll,yVecEq[[i]])
      xMatAll <- rbind( cbind( xMatAll, matrix( 0, nrow( xMatAll ), ncol( xMatEq[[i]] ))),
                  cbind( matrix( 0, nrow( xMatEq[[i]] ), ncol( xMatAll )), xMatEq[[i]]))
      nObsEq[i] <- length( yVecEq[[i]] )
      nExogEq[i] <- ncol(xMatEq[[i]])
      for(j in 1:nExogEq[i]) {
         xnames <- c( xnames, paste( eqnLabels[ i ],colnames( xMatEq[[i]] )[j],
            sep = "_" ))
      }
   }
   if( nEq > 1 ) {
      if( var ( nObsEq ) != 0 ) {
         stop( "Systems with unequal numbers of observations are not supported yet." )
      }
   }
   if( !is.null( TX ) ) {
      XU <- xMatAll
      xMatAll  <- XU %*% TX
   }
   result$termsEq    <- termsEq
   result$evalModelFrameEq <- evalModelFrameEq
   result$yVecEq     <- yVecEq
   result$xMatEq     <- xMatEq
   result$yVecAll    <- yVecAll
   result$xMatAll    <- xMatAll
   result$nObsEq     <- nObsEq
   result$nExogEq    <- nExogEq
   result$xnames     <- xnames

   ## preparing instruments
   if( !is.null( inst ) ) {
      for(i in 1:nEq) {
         if(is.list(inst)) {
            instEq[[i]] <- inst[[i]]
         } else {
            instEq[[i]] <- inst
         }
      }
      xMatHatAll <- matrix( 0, 0, ncol( xMatAll ) ) # fitted X values
      hMatAll  <- matrix( 0, 0, 0 )           # stacked matrices of all instruments
      hMatEq  <- list()
      # terms and model frames for the instruments of the individual equations
      termsInst <- list()
      modelFrameInst <- list()
      evalModelFrameInst <- list()
      # prepare data for individual equations
      for(i in 1:nEq) {
         rowsEq <- c( (1+sum(nObsEq[1:i])-nObsEq[i]):(sum(nObsEq[1:i])) )
            # rows that belong to the ith equation
         modelFrameInst[[ i ]] <- modelFrame
         modelFrameInst[[ i ]]$formula <- instEq[[ i ]]
         evalModelFrameInst[[ i ]] <- eval( modelFrameInst[[ i ]], parent.frame() )
         termsInst[[ i ]] <- attr( evalModelFrameInst[[ i ]], "terms" )
         hMatEq[[i]] <- model.matrix( termsInst[[ i ]], evalModelFrameInst[[ i ]] )
         if( nrow( hMatEq[[ i ]] ) != nrow( xMatAll[ rowsEq, ] ) ) {
            stop( paste( "The instruments and the regressors of equation",
               as.character( i ), "have different numbers of observations." ) )
         }
         # extract instrument matrix
         xMatHatAll <- rbind(xMatHatAll, hMatEq[[i]] %*% solve( crossprod( hMatEq[[i]]) , tol=control$solvetol )
              %*% crossprod( hMatEq[[i]], xMatAll[ rowsEq, ] ))       # 'fitted' X-values
         hMatAll  <-  rbind( cbind( hMatAll, matrix( 0, nrow( hMatAll ), ncol( hMatEq[[i]] ))),
                         cbind( matrix( 0, nrow( hMatEq[[i]] ), ncol( hMatAll )), hMatEq[[i]]))
      }
      result$instEq     <- instEq
      result$termsInst  <- termsInst
      result$evalModelFrameInst <- evalModelFrameInst
      result$xMatHatAll <- xMatHatAll
      result$hMatAll    <- hMatAll
      result$hMatEq     <- hMatEq
   }

   return( result )
}




