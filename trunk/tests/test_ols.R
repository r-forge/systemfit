
library( systemfit )
data( kmenta )

demand <- q ~ p + d
supply <- q ~ p + f + a
inst   <- ~ d + f + a
inst1  <- ~ d + f
instlist <- list( inst1, inst )
labels <- list( "demand", "supply" )
system <- list( demand, supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2q <- matrix(0,2,1)  # restriction vector "q" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q[2,1] <-  0.5
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3q <- matrix(0,1,1)  # restriction vector "q" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q[1,1] <-  0.5


## *************** OLS estimation ************************
fitols1 <- systemfit( "OLS", system, labels, data = kmenta )
print( fitols1 )
print( fitols1$bcov )

## ********** OLS estimation (no single.eq sigma) ******************
fitols1n <- systemfit( "OLS", system, labels, data = kmenta,
   single.eq.sigma = FALSE )
print( fitols1n )
print( fitols1n$bcov )

## ****************  OLS (probdfsys = TRUE) ***********************
fitols1e <- systemfit( "OLS", system, labels, data = kmenta,
   probdfsys = TRUE )
print( fitols1e )
print( fitols1e$bcov )

## ****************  OLS (rcovformula = 0) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = kmenta,
   rcovformula = 0 )
print( fitols1r )
print( fitols1r$bcov )
# It is not possible to estimate OLS without any restrictions
# with systemfit exactly as EViews does, because in the absence
# of cross-equation restrictions EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## ********* OLS with cross-equation restriction ************
fitols2 <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm )
print( fitols2 )
print( fitols2$bcov )

## ****** OLS with cross-equation restriction (EViews-like) *******
fitols2e <- systemfit( "OLS", system, labels, data = kmenta,
   R.restr = restrm, rcovformula = 0 )
print( fitols2e )
print( fitols2e$bcov )

## *** OLS with cross-equation restriction via TX ***
fitols3 <- systemfit( "OLS", system, labels, data = kmenta, TX = tc )
print( fitols3 )
print( fitols3$bcov )

## *** OLS with cross-equation restriction via TX (EViews-like) ***
fitols3e <- systemfit( "OLS", system, labels, data = kmenta,
   TX = tc, rcovformula = 0 )
print( fitols3e )
print( fitols3e$bcov )

## ********* OLS with 2 cross-equation restrictions ***********
fitols4 <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( fitols4 )
print( fitols4$bcov )

## ****** OLS with 2 cross-equation restrictions (EViews-like) *******
fitols4e <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0 )
print( fitols4e )
print( fitols4e$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX ****
fitols5 <- systemfit( "OLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0)
print( fitols5 )
print( fitols5$bcov )

## ***** OLS with 2 cross-equation restrictions via R and TX (EViews-like) ****
fitols5e <- systemfit( "OLS", system, labels, data = kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0 )
print( fitols5e )
print( fitols5e$bcov )

