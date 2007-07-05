library( MASS )

nEq <- 6
nExog <- 10
nObs <- 2000

xMat   <- NULL
xMatEq <- list()
for( eqNo in 1:nEq ){
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
}
myCoef <- c( 10:( 9 + nEq * ( nExog + 1 ) ) ) / 10

sigma <- diag( c( 1:nEq ) ) * 20
sigma[ upper.tri( sigma ) ] <- sigma[ lower.tri( sigma ) ] <- 12

disturbances <- mvrnorm( nObs, rep( 0, nEq ), sigma )

yVec <- xMat %*% myCoef + c( disturbances )

## OLS
print( system.time(
   olsNaive <- solve( t( xMat ) %*% xMat ) %*% t( xMat ) %*%  yVec
) )

print( system.time(
   olsSolve <- solve( t( xMat ) %*% xMat, t( xMat ) %*%  yVec )
) )

print( system.time(
   olsCross <- solve( crossprod( xMat ) ) %*% crossprod( xMat, yVec )
) )

print( system.time(
   olsCrossSolve <- solve( crossprod( xMat ), crossprod( xMat, yVec ) )
) )


## residuals
residVec <- yVec - xMat %*% olsNaive
residMat <- matrix( residVec, nrow = nObs, ncol = nEq )
residCov <- crossprod( residMat ) / nObs


## SUR ( GLS )
print( system.time(
   omegaInv <- solve( residCov ) %x% diag( nObs )
))
print( system.time(
   surNaive <- solve( t( xMat ) %*% omegaInv %*% xMat ) %*%
      t( xMat ) %*% omegaInv %*%  yVec
) )


