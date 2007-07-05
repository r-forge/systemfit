nEq <- 3
nExog <- 4
nObs <- 1000

# myData <- data.frame( obsNo = c( 1:nObs ) )
# yNames <- paste( "y", 1:nEq, sep = "." )
# xNames <- matrix(
#    paste( "x", rep( 1:nEq, nExog ), rep( 1:nExog, each = nEq ), sep = "." ),
#    nrow = nEq, ncol = nExog )
# eqns <- list()

xMat   <- NULL
xMatEq <- list()
for( eqNo in 1:nEq ){
   yVecEq[[ eqNo ]] <- rep( eqNo, nObs )
   xMatEq[[ eqNo ]] <- matrix( 1, nrow = nObs, ncol = 1 )
   for( exogNo in 1:nExog ){
      xMatEq[[ eqNo ]] <- cbind( xMatEq[[ eqNo ]], rnorm( nObs ) )
   }
   if( eqNo == 1 ){
      xMat <- xMatEq[[ eqNo ]]
   } else {
      xMat <- rbind(
         cbind( xMat, matrix( 0, nrow = nrow( xMat ), ncol = nExog + 1 ) ),
         cbind( matrix( 0, nrow = nObs, ncol = ncol( xMat ) ), xMatEq[[ eqNo ]] ) )
   }
#    eqns[[ eqNo ]] <- as.formula( paste( "y.", eqNo, " ~ ",
#       paste( paste( "x.", eqNo, ".", sep = "" ),
#          c( 1:nExog ), sep = "", collapse = " + " ), sep = "" ) )
}
myCoef <- c( 1:( nEq * ( nExog + 1 ) ) )

yVec <- xMat %*% myCoef

myCoefEst <- solve( crossprod( xMat ), crossprod( xMat, yVec ) )


