
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
fitsur1r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1r2 ) )
print( round( fitsur1r2$bcov, digits = 6 ) )

## *************** SUR (methodRCov="Theil", probDfSys = TRUE ) ***************
fitsur1e2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1e2, probDfSys = TRUE ) )
print( round( fitsur1e2$bcov, digits = 6 ) )

## ********************* SUR (methodRCov="max") *****************
fitsur1r3 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max" )
print( summary( fitsur1r3 ) )
print( round( fitsur1r3$bcov, digits = 6 ) )

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
fitsur4r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4r2 ) )
print( round( fitsur4r2$bcov, digits = 6 ) )

## *************** SUR with 2 restrictions (methodRCov = "max") **************
fitsur4r3 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
   R.restr = restr2m, q.restr = restr2q )
print( summary( fitsur4r3 ) )
print( round( fitsur4r3$bcov, digits = 6 ) )

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
fitsuri1r2 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "Theil" )
print( summary( fitsuri1r2 ) )
print( round( fitsuri1r2$bcov, digits = 6 ) )

## ************** iterated SUR (methodRCov="Theil", probDfSys=TRUE) *****************
fitsuri1e2 <- systemfit( system2, "SUR", data = Kmenta, methodRCov = "Theil",
   maxit = 100 )
print( summary( fitsuri1e2, probDfSys = TRUE ) )
print( round( fitsuri1e2$bcov, digits = 6 ) )

## ************** iterated SUR (methodRCov = "max") ****************************
fitsuri1r3 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "max" )
print( summary( fitsuri1r3 ) )
print( round( fitsuri1r3$bcov, digits = 6 ) )

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
fitsuri5r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5r2 ) )
print( round( fitsuri5r2$bcov, digits = 6 ) )

## ********* iterated SUR with 2 restrictions via R and TX (methodRCov="max") **********
# fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
#    R.restr = restr3m, q.restr = restr3q, TX = tc, maxit = 100 )
# print( summary( fitsuri5e ) )
# print( round( fitsuri5e$bcov, digits = 6 ) )
# disabled, because the estimation does not converge


## ****************** residuals **************************
print( residuals( fitsur1e2 ) )
print( residuals( fitsur1e2$eq[[ 2 ]] ) )

print( residuals( fitsur2e ) )
print( residuals( fitsur2e$eq[[ 1 ]] ) )

print( residuals( fitsur3 ) )
print( residuals( fitsur3$eq[[ 2 ]] ) )

print( residuals( fitsur4r3 ) )
print( residuals( fitsur4r3$eq[[ 1 ]] ) )

print( residuals( fitsur5 ) )
print( residuals( fitsur5$eq[[ 2 ]] ) )

print( residuals( fitsuri1r3 ) )
print( residuals( fitsuri1r3$eq[[ 1 ]] ) )

print( residuals( fitsuri2 ) )
print( residuals( fitsuri2$eq[[ 2 ]] ) )

print( residuals( fitsuri3e ) )
print( residuals( fitsuri3e$eq[[ 1 ]] ) )

print( residuals( fitsuri4 ) )
print( residuals( fitsuri4$eq[[ 2 ]] ) )

print( residuals( fitsuri5r2 ) )
print( residuals( fitsuri5r2$eq[[ 1 ]] ) )


## *********** confidence intervals of coefficients *************
print( confint( fitsur1e2, probDfSys = TRUE ) )
print( confint( fitsur1e2$eq[[ 2 ]], level = 0.9, probDfSys = TRUE ) )

print( confint( fitsur2e, level = 0.9 ) )
print( confint( fitsur2e$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitsur3, level = 0.99 ) )
print( confint( fitsur3$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitsur4r3, level = 0.5 ) )
print( confint( fitsur4r3$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitsur5, level = 0.25 ) )
print( confint( fitsur5$eq[[ 2 ]], level = 0.975 ) )

print( confint( fitsuri1r3, level = 0.975 ) )
print( confint( fitsuri1r3$eq[[ 1 ]], level = 0.999 ) )

print( confint( fitsuri2, level = 0.999 ) )
print( confint( fitsuri2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitsuri3e, level = 0.1 ) )
print( confint( fitsuri3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitsuri4, level = 0.01 ) )
print( confint( fitsuri4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitsuri5r2, level = 0.33 ) )
print( confint( fitsuri5r2$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitsur1e2, se.fit = TRUE, interval = "prediction",
   probDfSys = TRUE ) )
print( predict( fitsur1e2$eq[[ 2 ]] ) )

print( predict( fitsur2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fitsur2e$eq[[ 1 ]] ) )

print( predict( fitsur3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitsur3$eq[[ 2 ]] ) )

print( predict( fitsur4r3, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitsur4r3$eq[[ 1 ]] ) )

print( predict( fitsur5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fitsur5$eq[[ 2 ]] ) )

print( predict( fitsuri1r3, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitsuri1r3$eq[[ 1 ]] ) )

print( predict( fitsuri2, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fitsuri2$eq[[ 2 ]] ) )

print( predict( fitsuri3e, interval = "prediction", level = 0.925 ) )
print( predict( fitsuri3e$eq[[ 1 ]] ) )

print( predict( fitsuri4, interval = "confidence", data = predictData ) )
print( predict( fitsuri4$eq[[ 2 ]] ) )

print( predict( fitsuri5r2 ) )
print( predict( fitsuri5r2$eq[[ 1 ]] ) )


## *********** likelihood ratio tests *************
# testing first restriction
# non-iterating, methodRCov = 1
print( lrtest.systemfit( fitsur2, fitsur1 ) )
print( lrtest.systemfit( fitsur3, fitsur1 ) )
# non-iterating, methodRCov = 0
print( lrtest.systemfit( fitsur2e, fitsur1e ) )
print( lrtest.systemfit( fitsur3e, fitsur1e ) )
# iterating, methodRCov = 1
print( lrtest.systemfit( fitsuri2, fitsuri1 ) )
print( lrtest.systemfit( fitsuri3, fitsuri1 ) )
# iterating, methodRCov = 0
print( lrtest.systemfit( fitsuri2e, fitsuri1e ) )
print( lrtest.systemfit( fitsuri3e, fitsuri1e ) )

# testing second restriction
# non-iterating, methodRCov = 1
print( lrtest.systemfit( fitsur4, fitsur2 ) )
print( lrtest.systemfit( fitsur4, fitsur3 ) )
print( lrtest.systemfit( fitsur5, fitsur2 ) )
print( lrtest.systemfit( fitsur5, fitsur3 ) )
# non-iterating, methodRCov = 0
print( lrtest.systemfit( fitsur4e, fitsur2e ) )
print( lrtest.systemfit( fitsur4e, fitsur3e ) )
print( lrtest.systemfit( fitsur5e, fitsur2e ) )
print( lrtest.systemfit( fitsur5e, fitsur3e ) )
# iterating, methodRCov = 1
print( lrtest.systemfit( fitsuri4, fitsuri2 ) )
print( lrtest.systemfit( fitsuri4, fitsuri3 ) )
print( lrtest.systemfit( fitsuri5, fitsuri2 ) )
print( lrtest.systemfit( fitsuri5, fitsuri3 ) )
# iterating, methodRCov = 0
print( lrtest.systemfit( fitsuri4e, fitsuri2e ) )
print( lrtest.systemfit( fitsuri4e, fitsuri3e ) )
print( lrtest.systemfit( fitsuri5e, fitsuri2e ) )
print( lrtest.systemfit( fitsuri5e, fitsuri3e ) )

# testing both of the restrictions
# non-iterating, methodRCov = 1
print( lrtest.systemfit( fitsur4, fitsur1 ) )
print( lrtest.systemfit( fitsur5, fitsur1 ) )
# non-iterating, methodRCov = 0
print( lrtest.systemfit( fitsur4e, fitsur1e ) )
print( lrtest.systemfit( fitsur5e, fitsur1e ) )
# iterating, methodRCov = 1
print( lrtest.systemfit( fitsuri4, fitsuri1 ) )
print( lrtest.systemfit( fitsuri5, fitsuri1 ) )
# iterating, methodRCov = 0
print( lrtest.systemfit( fitsuri4e, fitsuri1e ) )
print( lrtest.systemfit( fitsuri5e, fitsuri1e ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fitsur1, restrm ) )
print( ftest.systemfit( fitsur1r2, restrm ) )
print( ftest.systemfit( fitsuri1e2, restrm ) )
print( ftest.systemfit( fitsuri1r3, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fitsur1e2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitsuri1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fitsur2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitsur3, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitsuri2e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitsuri3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fitsur1r3, restr2m, restr2q ) )
print( ftest.systemfit( fitsuri1e2, restr2m, restr2q ) )
