
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
inst1  <- ~ income + farmPrice
instlist <- list( inst1, inst )
system <- list( demand = demand, supply = supply )
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

# It is not possible to estimate 2SLS with systemfit exactly
# as EViews does, because EViews uses
# methodRCov == "geomean" for the coefficient covariance matrix and
# methodRCov == "noDfCor" for the residual covariance matrix.
# systemfit uses always the same formulas for both calculations.

## *************** 2SLS estimation ************************
## ************ 2SLS estimation (default)*********************
fit2sls1 <- systemfit( system, "2SLS", data = Kmenta, inst = inst )
print( summary( fit2sls1 ) )
print( round( fit2sls1$bcov, digits = 6 ) )

## *************** 2SLS estimation (single.eq.sigma=F)*******************
fit2sls1s <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   single.eq.sigma = FALSE )
print( summary( fit2sls1s ) )
print( round( fit2sls1s$bcov, digits = 6 ) )

## ********************* 2SLS (probDfSys = TRUE) *****************
fit2sls1p <- systemfit( system, "2SLS", data = Kmenta, inst = inst )
print( summary( fit2sls1p, probDfSys = TRUE ) )
print( round( fit2sls1p$bcov, digits = 6 ) )

## ********************* 2SLS (methodRCov = "noDfCor" ) *****************
fit2sls1r <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   methodRCov = "noDfCor" )
print( summary( fit2sls1r ) )
print( round( fit2sls1r$bcov, digits = 6 ) )

## *************** 2SLS (methodRCov="noDfCor", single.eq.sigma=F) *************
fit2sls1rs <- systemfit( system, "2SLS", data = Kmenta, inst = inst,
   methodRCov = "noDfCor", single.eq.sigma = FALSE )
print( summary( fit2sls1rs ) )
print( round( fit2sls1rs$bcov, digits = 6 ) )

## ********************* 2SLS with restriction ********************
## **************** 2SLS with restriction (default)********************
fit2sls2 <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fit2sls2 ) )
print( round( fit2sls2$bcov, digits = 6 ) )

## ************* 2SLS with restriction (single.eq.sigma=T) *****************
fit2sls2s <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls2s ) )
print( round( fit2sls2s$bcov, digits = 6 ) )

## ********************* 2SLS with restriction (probDfSys=T) **************
fit2sls2p <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fit2sls2p, probDfSys = TRUE ) )
print( round( fit2sls2p$bcov, digits = 6 ) )

## ********************* 2SLS with restriction (methodRCov = "noDfCor") **************
fit2sls2r <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = inst, methodRCov = "noDfCor" )
print( summary( fit2sls2r ) )
print( round( fit2sls2r$bcov, digits = 6 ) )

## ******** 2SLS with restriction (methodRCov="noDfCor", single.eq.sigma=TRUE) *********
fit2sls2rs <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = inst, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fit2sls2rs ) )
print( round( fit2sls2rs$bcov, digits = 6 ) )

## ********************* 2SLS with restriction via TX ******************
## *************** 2SLS with restriction via TX (default )***************
fit2sls3 <- systemfit( system, "2SLS", data = Kmenta, TX = tc,
   inst = inst, methodRCov = "noDfCor" )
print( summary( fit2sls3, probDfSys = TRUE ) )
print( round( fit2sls3$bcov, digits = 6 ) )

## ********************* 2SLS with restriction via TX (EViews-like) *******
fit2sls3e <- systemfit( system, "2SLS", data = Kmenta, TX = tc,
   inst = inst, methodRCov = "noDfCor" )
print( summary( fit2sls3e, probDfSys = TRUE ) )
print( round( fit2sls3e$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions *******************
## ************** 2SLS with 2 restrictions (default) **************
fit2sls4 <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fit2sls4 ) )
print( round( fit2sls4$bcov, digits = 6 ) )

## ************ 2SLS with 2 restrictions (single.eq.sigma=T) **************
fit2sls4s <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls4s ) )
print( round( fit2sls4s$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions (probDfSys=T) **************
fit2sls4p <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fit2sls4p, probDfSys = TRUE ) )
print( round( fit2sls4p$bcov, digits = 6 ) )

## ***************** 2SLS with 2 restrictions (methodRCov="noDfCor") **************
fit2sls4r <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, methodRCov = "noDfCor" )
print( summary( fit2sls4r ) )
print( round( fit2sls4r$bcov, digits = 6 ) )

## ***** 2SLS with 2 restrictions (methodRCov="noDfCor", single.eq.sigma=T) *******
fit2sls4rs <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fit2sls4rs ) )
print( round( fit2sls4rs$bcov, digits = 6 ) )

## ************* 2SLS with 2 restrictions via R and TX ******************
## ******** 2SLS with 2 restrictions via R and TX (default) *************
fit2sls5 <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fit2sls5 ) )
print( round( fit2sls5$bcov, digits = 6 ) )

## ******* 2SLS with 2 restrictions via R and TX (single.eq.sigma=T) ******
fit2sls5s <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, single.eq.sigma = TRUE )
print( summary( fit2sls5s ) )
print( round( fit2sls5s$bcov, digits = 6 ) )

## ********** 2SLS with 2 restrictions via R and TX (probDfSys=T) *******
fit2sls5p <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fit2sls5p, probDfSys = TRUE ) )
print( round( fit2sls5p$bcov, digits = 6 ) )

## ************* 2SLS with 2 restrictions via R and TX (methodRCov="noDfCor") *********
fit2sls5r <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, methodRCov = "noDfCor" )
print( summary( fit2sls5r ) )
print( round( fit2sls5r$bcov, digits = 6 ) )

## ** 2SLS with 2 restrictions via R and TX (methodRCov="noDfCor", single.eq.sigma=T) **
fit2sls5rs <- systemfit( system, "2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fit2sls5rs ) )
print( round( fit2sls5rs$bcov, digits = 6 ) )

## *********** 2SLS estimation with different instruments **************
## ******* 2SLS estimation with different instruments (default) *********
fit2slsd1 <- systemfit( system, "2SLS", data = Kmenta, inst = instlist )
print( summary( fit2slsd1 ) )
print( round( fit2slsd1$bcov, digits = 6 ) )

## *********** 2SLS estimation with different instruments (single.eq.sigma=F)*****
fit2slsd1s <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   single.eq.sigma = FALSE )
print( summary( fit2slsd1s ) )
print( round( fit2slsd1s$bcov, digits = 6 ) )

## ********* 2SLS estimation with different instruments (probDfSys=T) *******
fit2slsd1p <- systemfit( system, "2SLS", data = Kmenta, inst = instlist )
print( summary( fit2slsd1p, probDfSys = TRUE ) )
print( round( fit2slsd1p$bcov, digits = 6 ) )

## ********* 2SLS estimation with different instruments (methodRCov="noDfCor") ******
fit2slsd1r <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   methodRCov = "noDfCor" )
print( summary( fit2slsd1r ) )
print( round( fit2slsd1r$bcov, digits = 6 ) )

## 2SLS estimation with different instruments (methodRCov="noDfCor",single.eq.sigma=F)
fit2slsd1r <- systemfit( system, "2SLS", data = Kmenta, inst = instlist,
   methodRCov = "noDfCor", single.eq.sigma = FALSE )
print( summary( fit2slsd1r ) )
print( round( fit2slsd1r$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction *******
## ** 2SLS estimation with different instruments and restriction (default) ****
fit2slsd2 <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fit2slsd2 ) )
print( round( fit2slsd2$bcov, digits = 6 ) )

## 2SLS estimation with different instruments and restriction (single.eq.sigma=T)
fit2slsd2s <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist, single.eq.sigma = TRUE )
print( summary( fit2slsd2s ) )
print( round( fit2slsd2s$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction (probDfSys=F)
fit2slsd2p <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fit2slsd2p, probDfSys = FALSE ) )
print( round( fit2slsd2p$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction (methodRCov="noDfCor")
fit2slsd2r <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fit2slsd2r ) )
print( round( fit2slsd2r$bcov, digits = 6 ) )

## 2SLS estimation with different instr. and restr. (methodRCov="noDfCor", single.eq.sigma=T)
fit2slsd2rs <- systemfit( system, "2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fit2slsd2rs ) )
print( round( fit2slsd2rs$bcov, digits = 6 ) )

## **** 2SLS estimation with different instruments and restriction via TX *
## 2SLS estimation with different instruments and restriction via TX (default)
fit2slsd3 <- systemfit( system, "2SLS", data = Kmenta, TX = tc,
   inst = instlist )
print( summary( fit2slsd3 ) )
print( round( fit2slsd3$bcov, digits = 6 ) )

## **** 2SLS estimation with different instr. and restr. via TX (methodRCov="noDfCor")
fit2slsd3r <- systemfit( system, "2SLS", data = Kmenta, TX = tc,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fit2slsd3r ) )
print( round( fit2slsd3r$bcov, digits = 6 ) )


## ****************** residuals **************************
print( residuals( fit2sls1 ) )
print( residuals( fit2sls1$eq[[ 1 ]] ) )

print( residuals( fit2sls2s ) )
print( residuals( fit2sls2s$eq[[ 2 ]] ) )

print( residuals( fit2sls3e ) )
print( residuals( fit2sls3e$eq[[ 1 ]] ) )

print( residuals( fit2sls4r ) )
print( residuals( fit2sls4r$eq[[ 2 ]] ) )

print( residuals( fit2sls5rs ) )
print( residuals( fit2sls5rs$eq[[ 1 ]] ) )

print( residuals( fit2slsd1p ) )
print( residuals( fit2slsd1p$eq[[ 2 ]] ) )

print( residuals( fit2slsd2r ) )
print( residuals( fit2slsd2r$eq[[ 1 ]] ) )


## *********** confidence intervals of coefficients *************
print( confint( fit2sls1 ) )
print( confint( fit2sls1$eq[[ 1 ]], level = 0.9 ) )

print( confint( fit2sls2s, level = 0.9 ) )
print( confint( fit2sls2s$eq[[ 2 ]], level = 0.99 ) )

print( confint( fit2sls3e, level = 0.99, probDfSys = TRUE ) )
print( confint( fit2sls3e$eq[[ 1 ]], level = 0.5, probDfSys = TRUE ) )

print( confint( fit2sls4r, level = 0.5 ) )
print( confint( fit2sls4r$eq[[ 2 ]], level = 0.25 ) )

print( confint( fit2sls5rs, level = 0.25 ) )
print( confint( fit2sls5rs$eq[[ 1 ]], level = 0.975 ) )

print( confint( fit2slsd1p, level = 0.975, probDfSys = TRUE ) )
print( confint( fit2slsd1p$eq[[ 2 ]], level = 0.999, probDfSys = TRUE ) )

print( confint( fit2slsd2r, level = 0.999 ) )
print( confint( fit2slsd2r$eq[[ 1 ]] ) )

