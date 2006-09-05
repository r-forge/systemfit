
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

# the standard equations do not converge and lead to a singular weighting matrix
# both in R and in EViews, since both equations have the same endogenous variable
supply2 <- price ~ income + farmPrice + trend
system2 <- list( demand = demand, supply = supply2 )


## *************** SUR estimation ************************
fitsur1 <- systemfit( system, "SUR", data = Kmenta )
print( summary( fitsur1 ) )
print( round( fitsur1$bcov, digits = 6 ) )

## ********************* SUR (EViews-like) *****************
fitsur1e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor" )
print( summary( fitsur1e, probDfSys = TRUE ) )
print( round( fitsur1e$bcov, digits = 6 ) )

## ********************* SUR (methodRCov="Theil") *****************
fitsur1c <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1c ) )
print( round( fitsur1c$bcov, digits = 6 ) )

## *************** SUR (methodRCov="Theil", probDfSys = TRUE ) ***************
fitsur1cp <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1cp, probDfSys = TRUE ) )
print( round( fitsur1cp$bcov, digits = 6 ) )

## ********************* SUR (methodRCov="max") *****************
fitsur1c <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max" )
print( summary( fitsur1c ) )
print( round( fitsur1c$bcov, digits = 6 ) )

## *************** SUR with cross-equation restriction **************
fitsur2 <- systemfit( system, "SUR", data = Kmenta, R.restr = restrm )
print( summary( fitsur2 ) )
print( round( fitsur2$bcov, digits = 6 ) )

## *************** SUR with cross-equation restriction (EViews-like) **
fitsur2e <- systemfit( system, "SUR", data = Kmenta, R.restr = restrm,
   methodRCov = "noDfCor" )
print( summary( fitsur2e ) )
print( round( fitsur2e$bcov, digits = 6 ) )

## *************** SUR with restriction via TX *******************
fitsur3 <- systemfit( system, "SUR", data = Kmenta, TX = tc )
print( summary( fitsur3 ) )
print( round( fitsur3$bcov, digits = 6 ) )

## *************** SUR with restriction via TX (EViews-like) **************
fitsur3e <- systemfit( system, "SUR", data = Kmenta, TX = tc,
   methodRCov = "noDfCor" )
print( summary( fitsur3e ) )
print( round( fitsur3e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions ***************************
fitsur4 <- systemfit( system, "SUR", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitsur4 ) )
print( round( fitsur4$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (EViews-like) **************
fitsur4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (methodRCov = "Theil") **************
fitsur4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (methodRCov = "max") **************
fitsur4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4e ) )
print( round( fitsur4e$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions via R and TX ****************
fitsur5 <- systemfit( system, "SUR", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( summary( fitsur5 ) )
print( round( fitsur5$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions via R and TX (EViews-like) **************
fitsur5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr3m, q.restr = restr3q, TX = tc )
print( summary( fitsur5e ) )
print( round( fitsur5e$bcov, digits = 6 ) )

## ************** iterated SUR ****************************
fitsuri1 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100 )
print( summary( fitsuri1 ) )
print( round( fitsuri1$bcov, digits = 6 ) )

## ************** iterated SUR (EViews-like) *****************
fitsuri1e <- systemfit( system2, "SUR", data = Kmenta, methodRCov = "noDfCor",
   maxit = 100 )
print( summary( fitsuri1e, probDfSys = TRUE ) )
print( round( fitsuri1e$bcov, digits = 6 ) )

## ************** iterated SUR (methodRCov = "Theil") ****************************
fitsuri1c <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "Theil" )
print( summary( fitsuri1c ) )
print( round( fitsuri1c$bcov, digits = 6 ) )

## ************** iterated SUR (methodRCov="Theil", probDfSys=TRUE) *****************
fitsuri1cp <- systemfit( system2, "SUR", data = Kmenta, methodRCov = "Theil",
   maxit = 100 )
print( summary( fitsuri1cp, probDfSys = TRUE ) )
print( round( fitsuri1cp$bcov, digits = 6 ) )

## ************** iterated SUR (methodRCov = "max") ****************************
fitsuri1c <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "max" )
print( summary( fitsuri1c ) )
print( round( fitsuri1c$bcov, digits = 6 ) )

## *********** iterated SUR with restriction *******************
fitsuri2 <- systemfit( system2, "SUR", data = Kmenta, R.restr = restrm,
   maxit = 100 )
print( summary( fitsuri2 ) )
print( round( fitsuri2$bcov, digits = 6 ) )

## *********** iterated SUR with restriction (EViews-like) ***************
fitsuri2e <- systemfit( system2, "SUR", data = Kmenta, R.restr = restrm,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitsuri2e ) )
print( round( fitsuri2e$bcov, digits = 6 ) )

## *********** iterated SUR with restriction via TX ********************
fitsuri3 <- systemfit( system2, "SUR", data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitsuri3 ) )
print( round( fitsuri3$bcov, digits = 6 ) )

## *********** iterated SUR with restriction via TX (EViews-like) ***************
fitsuri3e <- systemfit( system2, "SUR", data = Kmenta, TX = tc,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitsuri3e ) )
print( round( fitsuri3e$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions ***************************
fitsuri4 <- systemfit( system, "SUR", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, maxit = 100 )
print( summary( fitsuri4 ) )
print( round( fitsuri4$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions (EViews-like) **************
fitsuri4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr2m, q.restr = restr2q, maxit = 100 )
print( summary( fitsuri4e ) )
print( round( fitsuri4e$bcov, digits = 6 ) )

## *************** iterated SUR with 2 restrictions via R and TX ****************
fitsuri5 <- systemfit( system, "SUR", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5 ) )
print( round( fitsuri5$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (EViews-like) **********
fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5e ) )
print( round( fitsuri5e$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (methodRCov="Theil") **********
fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5e ) )
print( round( fitsuri5e$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (methodRCov="max") **********
# fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
#    R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
# print( summary( fitsuri5e ) )
# print( round( fitsuri5e$bcov, digits = 6 ) )
# disabled, because the estimation does not converge


## ****************** residuals **************************
print( residuals( fitsur1cp ) )
print( residuals( fitsur1cp$eq[[ 2 ]] ) )

print( residuals( fitsur2e ) )
print( residuals( fitsur2e$eq[[ 1 ]] ) )

print( residuals( fitsur3 ) )
print( residuals( fitsur3$eq[[ 2 ]] ) )

print( residuals( fitsur4e ) )
print( residuals( fitsur4e$eq[[ 1 ]] ) )

print( residuals( fitsur5 ) )
print( residuals( fitsur5$eq[[ 2 ]] ) )

print( residuals( fitsuri1c ) )
print( residuals( fitsuri1c$eq[[ 1 ]] ) )

print( residuals( fitsuri2 ) )
print( residuals( fitsuri2$eq[[ 2 ]] ) )

print( residuals( fitsuri3e ) )
print( residuals( fitsuri3e$eq[[ 1 ]] ) )

print( residuals( fitsuri4 ) )
print( residuals( fitsuri4$eq[[ 2 ]] ) )

print( residuals( fitsuri5e ) )
print( residuals( fitsuri5e$eq[[ 1 ]] ) )


## *********** confidence intervals of coefficients *************
print( confint( fitsur1cp, probDfSys = TRUE ) )
print( confint( fitsur1cp$eq[[ 2 ]], level = 0.9, probDfSys = TRUE ) )

print( confint( fitsur2e, level = 0.9 ) )
print( confint( fitsur2e$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitsur3, level = 0.99 ) )
print( confint( fitsur3$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitsur4e, level = 0.5 ) )
print( confint( fitsur4e$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitsur5, level = 0.25 ) )
print( confint( fitsur5$eq[[ 2 ]], level = 0.975 ) )

print( confint( fitsuri1c, level = 0.975 ) )
print( confint( fitsuri1c$eq[[ 1 ]], level = 0.999 ) )

print( confint( fitsuri2, level = 0.999 ) )
print( confint( fitsuri2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitsuri3e, level = 0.1 ) )
print( confint( fitsuri3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitsuri4, level = 0.01 ) )
print( confint( fitsuri4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitsuri5e, level = 0.33 ) )
print( confint( fitsuri5e$eq[[ 1 ]] ) )
