
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

## ********************* SUR (EViews-like) *****************
fitsur1e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor" )
print( summary( fitsur1e, useDfSys = TRUE ) )

## ********************* SUR (methodRCov="Theil") *****************
fitsur1r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1r2 ) )

## *************** SUR (methodRCov="Theil", useDfSys = TRUE ) ***************
fitsur1e2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil" )
print( summary( fitsur1e2, useDfSys = TRUE ) )

## ********************* SUR (methodRCov="max") *****************
fitsur1r3 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max" )
print( summary( fitsur1r3 ) )

## *************** SUR with cross-equation restriction **************
fitsur2 <- systemfit( system, "SUR", data = Kmenta, restrictions = restrm )
print( summary( fitsur2 ) )

## *************** SUR with cross-equation restriction (EViews-like) **
fitsur2e <- systemfit( system, "SUR", data = Kmenta, restrictions = restrm,
   methodRCov = "noDfCor" )
print( summary( fitsur2e ) )

## *************** SUR with restriction via TX *******************
fitsur3 <- systemfit( system, "SUR", data = Kmenta, TX = tc )
print( summary( fitsur3 ) )

## *************** SUR with restriction via TX (EViews-like) **************
fitsur3e <- systemfit( system, "SUR", data = Kmenta, TX = tc,
   methodRCov = "noDfCor" )
print( summary( fitsur3e ) )

## *************** SUR with 2 restrictions ***************************
fitsur4 <- systemfit( system, "SUR", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q )
print( summary( fitsur4 ) )

## *************** SUR with 2 restrictions (EViews-like) **************
fitsur4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr2m, restrict.rhs = restr2q )
print( summary( fitsur4e ) )

## *************** SUR with 2 restrictions (methodRCov = "Theil") **************
fitsur4r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   restrictions = restr2m, restrict.rhs = restr2q )
print( summary( fitsur4r2 ) )

## *************** SUR with 2 restrictions (methodRCov = "max") **************
fitsur4r3 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
   restrictions = restr2m, restrict.rhs = restr2q )
print( summary( fitsur4r3 ) )

## *************** SUR with 2 restrictions via R and TX ****************
fitsur5 <- systemfit( system, "SUR", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc )
print( summary( fitsur5 ) )

## *************** SUR with 2 restrictions via R and TX (EViews-like) **************
fitsur5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr3m, restrict.rhs = restr3q, TX = tc )
print( summary( fitsur5e ) )

## ************** iterated SUR ****************************
fitsuri1 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100 )
print( summary( fitsuri1 ) )

## ************** iterated SUR (EViews-like) *****************
fitsuri1e <- systemfit( system2, "SUR", data = Kmenta, methodRCov = "noDfCor",
   maxit = 100 )
print( summary( fitsuri1e, useDfSys = TRUE ) )

## ************** iterated SUR (methodRCov = "Theil") ****************************
fitsuri1r2 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "Theil" )
print( summary( fitsuri1r2 ) )

## ************** iterated SUR (methodRCov="Theil", useDfSys=TRUE) *****************
fitsuri1e2 <- systemfit( system2, "SUR", data = Kmenta, methodRCov = "Theil",
   maxit = 100 )
print( summary( fitsuri1e2, useDfSys = TRUE ) )

## ************** iterated SUR (methodRCov = "max") ****************************
fitsuri1r3 <- systemfit( system2, "SUR", data = Kmenta, maxit = 100,
   methodRCov = "max" )
print( summary( fitsuri1r3 ) )

## *********** iterated SUR with restriction *******************
fitsuri2 <- systemfit( system2, "SUR", data = Kmenta, restrictions = restrm,
   maxit = 100 )
print( summary( fitsuri2 ) )

## *********** iterated SUR with restriction (EViews-like) ***************
fitsuri2e <- systemfit( system2, "SUR", data = Kmenta, restrictions = restrm,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitsuri2e ) )

## *********** iterated SUR with restriction via TX ********************
fitsuri3 <- systemfit( system2, "SUR", data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitsuri3 ) )

## *********** iterated SUR with restriction via TX (EViews-like) ***************
fitsuri3e <- systemfit( system2, "SUR", data = Kmenta, TX = tc,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitsuri3e ) )

## *************** iterated SUR with 2 restrictions ***************************
fitsuri4 <- systemfit( system, "SUR", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q, maxit = 100 )
print( summary( fitsuri4 ) )

## *************** iterated SUR with 2 restrictions (EViews-like) **************
fitsuri4e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr2m, restrict.rhs = restr2q, maxit = 100 )
print( summary( fitsuri4e ) )

## *************** iterated SUR with 2 restrictions via R and TX ****************
fitsuri5 <- systemfit( system, "SUR", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5 ) )

## ********* iterated SUR with 2 restrictions via R and TX (EViews-like) **********
fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr3m, restrict.rhs = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5e ) )

## ********* iterated SUR with 2 restrictions via R and TX (methodRCov="Theil") **********
fitsuri5r2 <- systemfit( system, "SUR", data = Kmenta, methodRCov = "Theil",
   restrictions = restr3m, restrict.rhs = restr3q, TX = tc, maxit = 100 )
print( summary( fitsuri5r2 ) )

## ********* iterated SUR with 2 restrictions via R and TX (methodRCov="max") **********
# fitsuri5e <- systemfit( system, "SUR", data = Kmenta, methodRCov = "max",
#    restrictions = restr3m, restrict.rhs = restr3q, TX = tc, maxit = 100 )
# print( summary( fitsuri5e ) )
# print( round( vcov( fitsuri5e ), digits = 6 ) )
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


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitsur1e2 ), digits = 6 ) )
print( round( vcov( fitsur1e2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur1r3 ), digits = 6 ) )
print( round( vcov( fitsur1r3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsur2e ), digits = 6 ) )
print( round( vcov( fitsur2e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur3 ), digits = 6 ) )
print( round( vcov( fitsur3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsur4r2 ), digits = 6 ) )
print( round( vcov( fitsur4r2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsur5e ), digits = 6 ) )
print( round( vcov( fitsur5e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsuri1r3 ), digits = 6 ) )
print( round( vcov( fitsuri1r3$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsuri2 ), digits = 6 ) )
print( round( vcov( fitsuri2$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsuri3e ), digits = 6 ) )
print( round( vcov( fitsuri3e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitsuri4e ), digits = 6 ) )
print( round( vcov( fitsuri4e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitsuri5r2 ), digits = 6 ) )
print( round( vcov( fitsuri5r2$eq[[ 1 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitsur1e2, useDfSys = TRUE ) )
print( confint( fitsur1e2$eq[[ 2 ]], level = 0.9, useDfSys = TRUE ) )

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


## *********** fitted values *************
print( fitted( fitsur1e2 ) )
print( fitted( fitsur1e2$eq[[ 2 ]] ) )

print( fitted( fitsur2e ) )
print( fitted( fitsur2e$eq[[ 1 ]] ) )

print( fitted( fitsur3 ) )
print( fitted( fitsur3$eq[[ 2 ]] ) )

print( fitted( fitsur4r3 ) )
print( fitted( fitsur4r3$eq[[ 1 ]] ) )

print( fitted( fitsur5 ) )
print( fitted( fitsur5$eq[[ 2 ]] ) )

print( fitted( fitsuri1r3 ) )
print( fitted( fitsuri1r3$eq[[ 1 ]] ) )

print( fitted( fitsuri2 ) )
print( fitted( fitsuri2$eq[[ 2 ]] ) )

print( fitted( fitsuri3e ) )
print( fitted( fitsuri3e$eq[[ 1 ]] ) )

print( fitted( fitsuri4 ) )
print( fitted( fitsuri4$eq[[ 2 ]] ) )

print( fitted( fitsuri5r2 ) )
print( fitted( fitsuri5r2$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitsur1e2, se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )
print( predict( fitsur1e2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )

print( predict( fitsur2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitsur2e$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitsur3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitsur3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fitsur4r3, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitsur4r3$eq[[ 1 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fitsur5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitsur5$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitsuri1r3, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitsuri1r3$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )

print( predict( fitsuri2, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitsuri2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fitsuri3e, interval = "prediction", level = 0.925 ) )
print( predict( fitsuri3e$eq[[ 1 ]], interval = "prediction", level = 0.925 ) )

print( predict( fitsuri4, interval = "confidence", newdata = predictData ) )
print( predict( fitsuri4$eq[[ 2 ]], interval = "confidence",
   newdata = predictData ) )

print( predict( fitsuri5r2 ) )
print( predict( fitsuri5r2$eq[[ 1 ]] ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitsur1e2, newdata = smallData ) )
print( predict( fitsur1e2$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitsur2e, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitsur2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitsur3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitsur3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitsur4r3, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitsur4r3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitsur5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitsur5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitsuri5r2, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitsuri5r2$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitsur1e2, 2, 1 ) )

print( correlation.systemfit( fitsur2e, 1, 2 ) )

print( correlation.systemfit( fitsur3, 2, 1 ) )

print( correlation.systemfit( fitsur4r3, 1, 2 ) )

print( correlation.systemfit( fitsur5, 2, 1 ) )

print( correlation.systemfit( fitsuri1r3, 1, 2 ) )

print( correlation.systemfit( fitsuri2, 2, 1 ) )

print( correlation.systemfit( fitsuri3e, 1, 2 ) )

print( correlation.systemfit( fitsuri4, 2, 1 ) )

print( correlation.systemfit( fitsuri5r2, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fitsur1e2 ) )

print( logLik( fitsur2e ) )

print( logLik( fitsur3 ) )

print( logLik( fitsur4r3 ) )

print( logLik( fitsur5 ) )

print( logLik( fitsuri1r3 ) )

print( logLik( fitsuri2 ) )

print( logLik( fitsuri3e ) )

print( logLik( fitsuri4 ) )

print( logLik( fitsuri5r2 ) )


## *********** likelihood ratio tests *************
# testing first restriction
# non-iterating, methodRCov = 1
print( lrtest( fitsur2, fitsur1 ) )
print( lrtest( fitsur3, fitsur1 ) )
# non-iterating, methodRCov = 0
print( lrtest( fitsur2e, fitsur1e ) )
print( lrtest( fitsur3e, fitsur1e ) )
# iterating, methodRCov = 1
print( lrtest( fitsuri2, fitsuri1 ) )
print( lrtest( fitsuri3, fitsuri1 ) )
# iterating, methodRCov = 0
print( lrtest( fitsuri2e, fitsuri1e ) )
print( lrtest( fitsuri3e, fitsuri1e ) )

# testing second restriction
# non-iterating, methodRCov = 1
print( lrtest( fitsur4, fitsur2 ) )
print( lrtest( fitsur4, fitsur3 ) )
print( lrtest( fitsur5, fitsur2 ) )
print( lrtest( fitsur5, fitsur3 ) )
# non-iterating, methodRCov = 0
print( lrtest( fitsur4e, fitsur2e ) )
print( lrtest( fitsur4e, fitsur3e ) )
print( lrtest( fitsur5e, fitsur2e ) )
print( lrtest( fitsur5e, fitsur3e ) )
# iterating, methodRCov = 1
print( lrtest( fitsuri4, fitsuri2 ) )
print( lrtest( fitsuri4, fitsuri3 ) )
print( lrtest( fitsuri5, fitsuri2 ) )
print( lrtest( fitsuri5, fitsuri3 ) )
# iterating, methodRCov = 0
print( lrtest( fitsuri4e, fitsuri2e ) )
print( lrtest( fitsuri4e, fitsuri3e ) )
print( lrtest( fitsuri5e, fitsuri2e ) )
print( lrtest( fitsuri5e, fitsuri3e ) )

# testing both of the restrictions
# non-iterating, methodRCov = 1
print( lrtest( fitsur4, fitsur1 ) )
print( lrtest( fitsur5, fitsur1 ) )
# non-iterating, methodRCov = 0
print( lrtest( fitsur4e, fitsur1e ) )
print( lrtest( fitsur5e, fitsur1e ) )
# iterating, methodRCov = 1
print( lrtest( fitsuri4, fitsuri1 ) )
print( lrtest( fitsuri5, fitsuri1 ) )
# iterating, methodRCov = 0
print( lrtest( fitsuri4e, fitsuri1e ) )
print( lrtest( fitsuri5e, fitsuri1e ) )


## ************** F tests ****************
# testing first restriction
print( linear.hypothesis( fitsur1, restrm ) )
print( linear.hypothesis( fitsur1r2, restrm ) )
print( linear.hypothesis( fitsuri1e2, restrm ) )
print( linear.hypothesis( fitsuri1r3, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( linear.hypothesis( fitsur1e2, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitsuri1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( linear.hypothesis( fitsur2, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitsur3, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitsuri2e, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitsuri3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( linear.hypothesis( fitsur1r3, restr2m, restr2q ) )
print( linear.hypothesis( fitsuri1e2, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( linear.hypothesis( fitsur1, restrm, test = "Chisq" ) )
print( linear.hypothesis( fitsur1r2, restrm, test = "Chisq" ) )
print( linear.hypothesis( fitsuri1e2, restrm, test = "Chisq" ) )
print( linear.hypothesis( fitsuri1r3, restrm, test = "Chisq" ) )

# testing second restriction
# first restriction not imposed
print( linear.hypothesis( fitsur1e2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitsuri1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
# first restriction imposed
print( linear.hypothesis( fitsur2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitsur3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitsuri2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitsuri3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )

# testing both of the restrictions
print( linear.hypothesis( fitsur1r3, restr2m, restr2q, test = "Chisq" ) )
print( linear.hypothesis( fitsuri1e2, restr2m, restr2q, test = "Chisq" ) )


## ****************** model frame **************************
print( mf <- model.frame( fitsur1e2 ) )
print( mf1 <- model.frame( fitsur1e2$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitsur1e2$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitsur2e ) ) )
print( all.equal( mf1, model.frame( fitsur2e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsur3 ) ) )
print( all.equal( mf2, model.frame( fitsur3$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsur4r3 ) ) )
print( all.equal( mf1, model.frame( fitsur4r3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsur5 ) ) )
print( all.equal( mf2, model.frame( fitsur5$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri1r3 ) ) )
print( all.equal( mf1, model.frame( fitsuri1r3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri2 ) ) )
print( all.equal( mf1, model.frame( fitsuri2$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri3e ) ) )
print( all.equal( mf1, model.frame( fitsuri3e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri4 ) ) )
print( all.equal( mf2, model.frame( fitsuri4$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitsuri5r2 ) ) )
print( all.equal( mf1, model.frame( fitsuri5r2$eq[[ 1 ]] ) ) )


## **************** model matrix ************************
print( mm <- model.matrix( fitsur1e2 ) )
print( mm1 <- model.matrix( fitsur1e2$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitsur1e2$eq[[ 2 ]] ) )
fitsur1e2$eq[[ 1 ]]$modelMatrix <- NULL
fitsur1e2$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitsur1e2 ) ) )
print( all.equal( mm1, model.matrix( fitsur1e2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur1e2$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitsur2 ) ) )
print( all.equal( mm1, model.matrix( fitsur2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2$eq[[ 2 ]] ) ) )
fitsur2$eq[[ 1 ]]$modelMatrix <- NULL
fitsur2$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitsur2 ) ) )
print( all.equal( mm1, model.matrix( fitsur2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur2$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitsur3e ) ) )
print( all.equal( mm1, model.matrix( fitsur3e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3e$eq[[ 2 ]] ) ) )
fitsur3e$eq[[ 1 ]]$modelMatrix <- NULL
fitsur3e$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitsur3e ) ) )
print( all.equal( mm1, model.matrix( fitsur3e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur3e$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitsur4r3 ) ) )
print( all.equal( mm1, model.matrix( fitsur4r3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur4r3$eq[[ 2 ]] ) ) )
fitsur4r3$eq[[ 1 ]]$modelMatrix <- NULL
fitsur4r3$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitsur4r3 ) ) )
print( all.equal( mm1, model.matrix( fitsur4r3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitsur4r3$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitsuri5r2 ) ) )
print( all.equal( mm1, model.matrix( fitsuri5r2$eq[[ 1 ]] ) ) )
#print( all.equal( mm2, model.matrix( fitsuri5r2$eq[[ 2 ]] ) ) )
fitsuri5r2$eq[[ 1 ]]$modelMatrix <- NULL
fitsuri5r2$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitsuri5r2 ) ) )
print( all.equal( mm1, model.matrix( fitsuri5r2$eq[[ 1 ]] ) ) )
#print( all.equal( mm2, model.matrix( fitsuri5r2$eq[[ 1 ]] ) ) )
