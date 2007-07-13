library( MASS )
set.seed( 20070705 )

nEq <- 6
nExog <- 10
nObs <- 2500

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
# Naive
system.time(
   olsNaive <- solve( t( xMat ) %*% xMat ) %*% t( xMat ) %*%  yVec
)

# with solve( , )
system.time(
   olsSolve <- solve( t( xMat ) %*% xMat, t( xMat ) %*%  yVec )
)
all.equal( olsNaive, olsSolve )

# with crossprod
system.time(
   olsCross <- solve( crossprod( xMat ) ) %*% crossprod( xMat, yVec )
)
all.equal( olsNaive, olsCross )

#with
system.time(
   olsCrossSolve <- solve( crossprod( xMat ), crossprod( xMat, yVec ) )
)
all.equal( olsNaive, olsCrossSolve )


## residuals
residVec <- yVec - xMat %*% olsNaive
residMat <- matrix( residVec, nrow = nObs, ncol = nEq )
residCov <- crossprod( residMat ) / nObs


## SUR ( GLS )
# Naive
system.time( {
   surNaive <- solve( t( xMat ) %*% ( solve( residCov ) %x% diag( nObs ) ) %*%
      xMat ) %*% t( xMat ) %*% ( solve( residCov ) %x% diag( nObs ) ) %*%  yVec
} )

# omegaInv
system.time( {
   omegaInv <- solve( residCov ) %x% diag( nObs )
   surOmegaInv <- solve( t( xMat ) %*% omegaInv %*% xMat,
      t( xMat ) %*% omegaInv %*%  yVec )
} )
all.equal( surNaive, surOmegaInv )

# Loop
system.time( {
   sigmaInv <- solve( residCov )
   xtOmegaInv <- matrix( 0, nrow = nEq * ( nExog + 1 ), ncol = nEq * nObs )
   eqSelect <- rep( NA, nEq * nObs )
   for( i in 1:nEq ){
      eqSelect[ ( ( i - 1 ) * nObs + 1 ):( i * nObs ) ] <- i
   }
   for( i in 1:nEq ){
      for( j in 1:nEq ){
         xtOmegaInv[ , eqSelect == i ] <- xtOmegaInv[ , eqSelect == i ] +
            t( xMat )[ , eqSelect == j ] * sigmaInv[ i, j ]
      }
   }
   surLoop <- solve( xtOmegaInv %*% xMat, xtOmegaInv %*%  yVec )
} )
all.equal( surNaive, surLoop )

# Matrix package
library( "Matrix" )
xMatM <- as( xMat, "dgCMatrix" )
residCovM <- as( residCov, "dspMatrix" )
system.time( {
   sigmaInvM <- solve( residCovM )
   omegaInvM <- kronecker( sigmaInvM, Diagonal( nObs ) )
   surMatrix <- solve( t( xMatM ) %*% omegaInvM %*% xMatM,
      t( xMatM ) %*% omegaInvM %*%  yVec )
} )
all.equal( surNaive, surMatrix )

xMatM <- as( xMat, "dgCMatrix" )
residCovM <- as( residCov, "dspMatrix" )
system.time( {
   sigmaInvM <- solve( residCovM )
   xtOmegaInvM <- crossprod( xMatM, kronecker( sigmaInvM, Diagonal( nObs ) ) )
   surMatrix2 <- solve( xtOmegaInvM %*% xMatM, xtOmegaInvM %*%  yVec )
} )
all.equal( surNaive, surMatrix )

system.time( sigmaInvM <- solve( residCovM ) )
system.time( omegaInvM <- kronecker( sigmaInvM, Diagonal( nObs ) ) )
system.time( surMatrix <- solve( t( xMatM ) %*% omegaInvM %*% xMatM,
      t( xMatM ) %*% omegaInvM %*%  yVec ) )

system.time( xtOmegaInvM <- t( xMatM ) %*% omegaInvM )

system.time( omegaInvM <- as( kronecker( sigmaInvM, Diagonal( nObs ) ), "dsCMatrix" ) )
system.time( xtOmegaInvM <- t( xMatM ) %*% omegaInvM )

system.time( xtOmegaInvM <- t( xMatM ) %*% as( kronecker( sigmaInvM, Diagonal( nObs ) ), "dsCMatrix" ) )

system.time( sigmaInvM <- solve( residCovM ) )
system.time( xtOmegaInvM <- t( xMatM ) %*% kronecker( sigmaInvM, Diagonal( nObs ) ) )
system.time( surMatrix2 <- solve( xtOmegaInvM %*% xMatM, xtOmegaInvM %*%  yVec ) )

system.time( xtOmegaInvM <- crossprod( xMatM, kronecker( sigmaInvM, Diagonal( nObs ) ) ) )
