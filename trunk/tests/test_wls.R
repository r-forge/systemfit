
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
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


## *************** WLS estimation ************************
fitwls1 <- systemfit( system, "WLS", data = Kmenta )
print( summary( fitwls1 ) )
print( round( fitwls1$bcov, digits = 6 ) )

## *************** WLS estimation (EViews-like) ************************
fitwls1e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor" )
print( summary( fitwls1e, probDfSys = TRUE ) )
print( round( fitwls1e$bcov, digits = 6 ) )

## ************** WLS with cross-equation restriction ***************
fitwls2 <- systemfit( system, "WLS", data = Kmenta, R.restr = restrm )
print( summary( fitwls2 ) )
print( round( fitwls2$bcov, digits = 6 ) )

## ************** WLS with cross-equation restriction (EViews-like) *******
fitwls2e <- systemfit( system, "WLS", data = Kmenta, R.restr = restrm,
   methodRCov = "noDfCor" )
print( summary( fitwls2e ) )
print( round( fitwls2e$bcov, digits = 6 ) )

## ******* WLS with cross-equation restriction via TX **********
fitwls3 <- systemfit( system,"WLS", data = Kmenta, TX = tc,)
print( summary( fitwls3 ) )
print( round( fitwls3$bcov, digits = 6 ) )

## ******* WLS with cross-equation restriction via TX (EViews-like) *****
fitwls3e <- systemfit( system,"WLS", data = Kmenta, TX = tc,
   methodRCov = "noDfCor" )
print( summary( fitwls3e ) )
print( round( fitwls3e$bcov, digits = 6 ) )

## ***** WLS with 2 cross-equation restrictions ***************
fitwls4 <- systemfit( system,"WLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitwls4 ) )
print( round( fitwls4$bcov, digits = 6 ) )

## ***** WLS with 2 cross-equation restrictions (EViews-like) **********
fitwls4e <- systemfit( system,"WLS", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitwls4e ) )
print( round( fitwls4e$bcov, digits = 6 ) )

## *********** WLS with 2 cross-equation restrictions via R and TX ******
fitwls5 <- systemfit( system, "WLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( summary( fitwls5 ) )
print( round( fitwls5$bcov, digits = 6 ) )

## *********** WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwls5e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr3m, q.restr = restr3q, TX = tc )
print( summary( fitwls5e ) )
print( round( fitwls5e$bcov, digits = 6 ) )

## *************** iterated WLS estimation *********************
fitwlsi1 <- systemfit( system, "WLS", data = Kmenta,
   maxit = 100 )
print( summary( fitwlsi1, probDfSys = TRUE ) )
print( round( fitwlsi1$bcov, digits = 6 ) )

## *************** iterated WLS estimation (EViews-like) ************
fitwlsi1e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   maxit = 100 )
print( summary( fitwlsi1e, probDfSys = TRUE ) )
print( round( fitwlsi1e$bcov, digits = 6 ) )

## ****** iterated WLS with cross-equation restriction ***************
fitwlsi2 <- systemfit( system, "WLS", data = Kmenta, R.restr = restrm,
   maxit = 100 )
print( summary( fitwlsi2 ) )
print( round( fitwlsi2$bcov, digits = 6 ) )

## ****** iterated WLS with cross-equation restriction (EViews-like) ********
fitwlsi2e <- systemfit( system, "WLS", data = Kmenta, R.restr = restrm,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitwlsi2e ) )
print( round( fitwlsi2e$bcov, digits = 6 ) )

## ******* iterated WLS with cross-equation restriction via TX **********
fitwlsi3 <- systemfit( system, "WLS", data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitwlsi3 ) )
print( round( fitwlsi3$bcov, digits = 6 ) )

## ******* iterated WLS with cross-equation restriction via TX (EViews-like) ***
fitwlsi3e <- systemfit( system, "WLS", data = Kmenta, TX = tc,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitwlsi3e ) )
print( round( fitwlsi3e$bcov, digits = 6 ) )

## ******* iterated WLS with 2 cross-equation restrictions ***********
fitwlsi4 <- systemfit( system, "WLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, maxit = 100 )
print( summary( fitwlsi4 ) )
print( round( fitwlsi4$bcov, digits = 6 ) )

## ******* iterated WLS with 2 cross-equation restrictions (EViews-like) *****
fitwlsi4e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr2m, q.restr = restr2q, maxit = 100 )
print( summary( fitwlsi4e ) )
print( round( fitwlsi4e$bcov, digits = 6 ) )

## ***** iterated WLS with 2 cross-equation restrictions via R and TX ******
fitwlsi5 <- systemfit( system, "WLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5 ) )
print( round( fitwlsi5$bcov, digits = 6 ) )

## *** iterated WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwlsi5e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5e ) )
print( round( fitwlsi5e$bcov, digits = 6 ) )


## ****************** residuals **************************
print( residuals( fitwls1 ) )
print( residuals( fitwls1$eq[[ 2 ]] ) )

print( residuals( fitwls2e ) )
print( residuals( fitwls2e$eq[[ 1 ]] ) )

print( residuals( fitwls3 ) )
print( residuals( fitwls3$eq[[ 2 ]] ) )

print( residuals( fitwls4e ) )
print( residuals( fitwls4e$eq[[ 1 ]] ) )

print( residuals( fitwls5 ) )
print( residuals( fitwls5$eq[[ 2 ]] ) )

print( residuals( fitwlsi1e ) )
print( residuals( fitwlsi1e$eq[[ 1 ]] ) )

print( residuals( fitwlsi2 ) )
print( residuals( fitwlsi2$eq[[ 2 ]] ) )

print( residuals( fitwlsi3e ) )
print( residuals( fitwlsi3e$eq[[ 1 ]] ) )

print( residuals( fitwlsi4 ) )
print( residuals( fitwlsi4$eq[[ 2 ]] ) )

print( residuals( fitwlsi5e ) )
print( residuals( fitwlsi5e$eq[[ 1 ]] ) )


## *********** confidence intervals of coefficients *************
print( confint( fitwls1 ) )
print( confint( fitwls1$eq[[ 2 ]], level = 0.9 ) )

print( confint( fitwls2e, level = 0.9 ) )
print( confint( fitwls2e$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitwls3, level = 0.99 ) )
print( confint( fitwls3$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitwls4e, level = 0.5 ) )
print( confint( fitwls4e$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitwls5, level = 0.25 ) )
print( confint( fitwls5$eq[[ 2 ]], level = 0.975 ) )

print( confint( fitwlsi1e, level = 0.975, probDfSys = TRUE ) )
print( confint( fitwlsi1e$eq[[ 1 ]], level = 0.999, probDfSys = TRUE ) )

print( confint( fitwlsi2, level = 0.999 ) )
print( confint( fitwlsi2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitwlsi3e, level = 0.1 ) )
print( confint( fitwlsi3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitwlsi4, level = 0.01 ) )
print( confint( fitwlsi4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitwlsi5e, level = 0.33 ) )
print( confint( fitwlsi5e$eq[[ 1 ]] ) )

