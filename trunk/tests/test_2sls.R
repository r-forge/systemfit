
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

## *************** 2SLS estimation ************************
fit2sls1 <- systemfit( "2SLS", system, labels, data = kmenta, inst = inst )
print( fit2sls1 )
print( fit2sls1$bcov )

## ********************* 2SLS (probdfsys = TRUE) *****************
fit2sls1e <- systemfit( "2SLS", system, labels, data = kmenta, inst = inst,
   probdfsys = TRUE )
print( fit2sls1e )
print( fit2sls1e$bcov )

## ********************* 2SLS (rcovformula = 0) *****************
fit2sls1e <- systemfit( "2SLS", system, labels, data = kmenta, inst = inst,
   rcovformula = 0 )
print( fit2sls1e )
print( fit2sls1e$bcov )
# It is not possible to estimate 2SLS without any restrictions
# with systemfit exactly as EViews does, because in the absence
# of cross-equation restrictions EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## ********************* 2SLS with restriction ********************
fit2sls2 <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = inst )
print( fit2sls2 )
print( fit2sls2$bcov )

## ********************* 2SLS with restriction (EViews-like) **************
fit2sls2e <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fit2sls2e )
print( fit2sls2e$bcov )

## ********************* 2SLS with restriction via TX ******************
fit2sls3 <- systemfit( "2SLS", system, labels, data = kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fit2sls3 )
print( fit2sls3$bcov )

## ********************* 2SLS with restriction via TX (EViews-like) *******
fit2sls3e <- systemfit( "2SLS", system, labels, data = kmenta, TX = tc,
   inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fit2sls3e )
print( fit2sls3e$bcov )

## ***************** 2SLS with 2 restrictions *******************
fit2sls4 <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( fit2sls4 )
print( fit2sls4$bcov )

## ***************** 2SLS with 2 restrictions (EViews-like) **************
fit2sls4e <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fit2sls4e )
print( fit2sls4e$bcov )

## ************* 2SLS with 2 restrictions via R and TX ******************
fit2sls5 <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( fit2sls5 )
print( fit2sls5$bcov )

## ************* 2SLS with 2 restrictions via R and TX (EViews-like) *********
fit2sls5e <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, rcovformula = 0, probdfsys = TRUE )
print( fit2sls5e )
print( fit2sls5e$bcov )

## *********** 2SLS estimation with different instruments **************
fit2slsd1 <- systemfit( "2SLS", system, labels, data = kmenta, inst = instlist )
print( fit2slsd1 )
print( fit2slsd1$bcov )

## *********** 2SLS estimation with different instruments (probdfsys = TRUE)********
fit2slsd1e <- systemfit( "2SLS", system, labels, data = kmenta, inst = instlist,
   probdfsys = TRUE )
print( fit2slsd1e )
print( fit2slsd1e$bcov )

## *********** 2SLS estimation with different instruments (rcovformula = 0)********
fit2slsd1e <- systemfit( "2SLS", system, labels, data = kmenta, inst = instlist,
   rcovformula = 0 )
print( fit2slsd1e )
print( fit2slsd1e$bcov )
# It is not possible to estimate 2SLS without any restrictions
# with systemfit exactly as EViews does, because in the absence
# of cross-equation restrictions EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## **** 2SLS estimation with different instruments and restriction *******
fit2slsd2 <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = instlist )
print( fit2slsd2 )
print( fit2slsd2$bcov )

## **** 2SLS estimation with different instruments and restriction (EViews-like)*
fit2slsd2e <- systemfit( "2SLS", system, labels, data = kmenta, R.restr = restrm,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( fit2slsd2e )
print( fit2slsd2e$bcov )

## **** 2SLS estimation with different instruments and restriction via TX *
fit2slsd3 <- systemfit( "2SLS", system, labels, data = kmenta, TX = tc,
   inst = instlist )
print( fit2slsd3 )
print( fit2slsd3$bcov )

## **** 2SLS estimation with different instruments and restriction via TX (EViews-like)*
fit2slsd3e <- systemfit( "2SLS", system, labels, data = kmenta, TX = tc,
   inst = instlist, rcovformula = 0, probdfsys = TRUE )
print( fit2slsd3e )
print( fit2slsd3e$bcov )
