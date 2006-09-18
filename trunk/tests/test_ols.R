
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

# It is not possible to estimate OLS with systemfit
# exactly as EViews does, because EViews uses
# methodRCov == "geomean" for the coefficient covariance matrix and
# methodRCov == "noDfCor" for the residual covariance matrix, while
# systemfit uses always the same formulas for both calculations.

## *************** OLS estimation ************************
## ********** OLS estimation (default) ********************
fitols1 <- systemfit( system, "OLS", data = Kmenta )
print( summary( fitols1 ) )

## ********** OLS estimation (no single.eq sigma=F) ******************
fitols1s <- systemfit( system, "OLS", data = Kmenta,
   single.eq.sigma = FALSE )
print( summary( fitols1s ) )

## ****************  OLS (probDfSys=T) ***********************
fitols1p <- systemfit( system, "OLS", data = Kmenta )
print( summary( fitols1p, probDfSys = TRUE ) )

## ****************  OLS (methodRCov="noDfCor") ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "noDfCor" )
print( summary( fitols1r ) )

## ********  OLS (methodRCov="noDfCor", single.eq.sigma=F) ***********
fitols1rs <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "noDfCor", single.eq.sigma = FALSE )
print( summary( fitols1rs ) )

## ****************  OLS (methodRCov="Theil" ) ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "Theil" )
print( summary( fitols1r ) )

## ****************  OLS (methodRCov="max") ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "max" )
print( summary( fitols1r ) )

## ********  OLS (methodRCov="max", single.eq.sigma=F) ***********
fitols1rs <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "max", single.eq.sigma = FALSE )
print( summary( fitols1rs ) )


## ********* OLS with cross-equation restriction ************
## ****** OLS with cross-equation restriction (default) *********
fitols2 <- systemfit( system, "OLS", data = Kmenta,
   R.restr = restrm )
print( summary( fitols2 ) )

## ****** OLS with cross-equation restriction (single.eq.sigma=T) *******
fitols2s <- systemfit( system, "OLS", data = Kmenta,
   R.restr = restrm, single.eq.sigma = TRUE )
print( summary( fitols2s ) )

## ****** OLS with cross-equation restriction (probDfSys=F) *******
fitols2p <- systemfit( system, "OLS", data = Kmenta,
   R.restr = restrm )
print( summary( fitols2p, probDfSys = FALSE ) )

## ****** OLS with cross-equation restriction (methodRCov="noDfCor") *******
fitols2r <- systemfit( system, "OLS", data = Kmenta,
   R.restr = restrm, methodRCov = "noDfCor" )
print( summary( fitols2r ) )

## ** OLS with cross-equation restriction (methodRCov="noDfCor",single.eq.sigma=T) ***
fitols2rs <- systemfit( system, "OLS", data = Kmenta,
   R.restr = restrm, methodRCov = "noDfCor" )
print( summary( fitols2rs ) )

## *** OLS with cross-equation restriction via TX ***
## *** OLS with cross-equation restriction via TX (default) ***
fitols3 <- systemfit( system, "OLS", data = Kmenta, TX = tc )
print( summary( fitols3 ) )

## *** OLS with cross-equation restriction via TX (single.eq.sigma=T) ***
fitols3s <- systemfit( system, "OLS", data = Kmenta,
   TX = tc, single.eq.sigma = TRUE )
print( summary( fitols3s ) )

## *** OLS with cross-equation restriction via TX (probDfSys=F) ***
fitols3p <- systemfit( system, "OLS", data = Kmenta,
   TX = tc )
print( summary( fitols3p, probDfSys = FALSE ) )

## *** OLS with cross-equation restriction via TX (methodRCov="noDfCor") ***
fitols3r <- systemfit( system, "OLS", data = Kmenta,
   TX = tc, methodRCov = "noDfCor" )
print( summary( fitols3r ) )

## OLS with cross-equation restriction via TX (methodRCov="noDfCor",single.eq.sigma=T)
fitols3rs <- systemfit( system, "OLS", data = Kmenta,
   TX = tc, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fitols3rs ) )

## ********* OLS with 2 cross-equation restrictions ***********
## ********* OLS with 2 cross-equation restrictions (default) ***********
fitols4 <- systemfit( system, "OLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitols4 ) )

## ****** OLS with 2 cross-equation restrictions (single.eq.sigma=T) *******
fitols4s <- systemfit( system, "OLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, single.eq.sigma = T )
print( summary( fitols4s ) )

## ****** OLS with 2 cross-equation restrictions (probDfSys=F) *******
fitols4p <- systemfit( system, "OLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitols4p, probDfSys = FALSE ) )

## ****** OLS with 2 cross-equation restrictions (methodRCov="noDfCor") *******
fitols4r <- systemfit( system, "OLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, methodRCov = "noDfCor" )
print( summary( fitols4r ) )

## OLS with 2 cross-equation restrictions (methodRCov="noDfCor", single.eq.sigma=T) *
fitols4rs <- systemfit( system, "OLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, methodRCov = "noDfCor", single.eq.sigma = T )
print( summary( fitols4rs ) )

## ***** OLS with 2 cross-equation restrictions via R and TX ****
## ***** OLS with 2 cross-equation restrictions via R and TX (default) ****
fitols5 <- systemfit( system, "OLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, methodRCov = "noDfCor")
print( summary( fitols5 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (single.eq.sigma=T) ****
fitols5s <- systemfit( system, "OLS", data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, single.eq.sigma = T )
print( summary( fitols5s ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (probDfSys=F) ****
fitols5p <- systemfit( system, "OLS", data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc )
print( summary( fitols5p, probDfSys = FALSE ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (methodRCov="noDfCor") ****
fitols5r <- systemfit( system, "OLS", data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, methodRCov = "noDfCor" )
print( summary( fitols5r ) )

## OLS with 2 cross-equation restr. via R and TX (methodRCov="noDfCor",single.eq.sigma=T)
fitols5rs <- systemfit( system, "OLS", data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, methodRCov = "noDfCor", single.eq.sigma = T )
print( summary( fitols5rs ) )


## ****************** residuals **************************
print( residuals( fitols1p ) )
print( residuals( fitols1p$eq[[ 2 ]] ) )

print( residuals( fitols2r ) )
print( residuals( fitols2r$eq[[ 1 ]] ) )

print( residuals( fitols3s ) )
print( residuals( fitols3s$eq[[ 2 ]] ) )

print( residuals( fitols4rs ) )
print( residuals( fitols4rs$eq[[ 1 ]] ) )

print( residuals( fitols5 ) )
print( residuals( fitols5$eq[[ 2 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitols1rs ), digits = 6 ) )
print( round( vcov( fitols1rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols2s ), digits = 6 ) )
print( round( vcov( fitols2s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols3p ), digits = 6 ) )
print( round( vcov( fitols3p$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols4r ), digits = 6 ) )
print( round( vcov( fitols4r$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols5 ), digits = 6 ) )
print( round( vcov( fitols5$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitols1p, probDfSys = TRUE ) )
print( confint( fitols1p$eq[[ 2 ]], level = 0.9, probDfSys = TRUE ) )

print( confint( fitols2r, level = 0.9 ) )
print( confint( fitols2r$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitols3s, level = 0.99 ) )
print( confint( fitols3s$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitols4rs, level = 0.5 ) )
print( confint( fitols4rs$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitols5, level = 0.25 ) )
print( confint( fitols5$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitols3p, level = 0.999, probDfSys = FALSE ) )
print( confint( fitols3p$eq[[ 1 ]], probDfSys = FALSE ) )


## *********** fitted values *************
print( fitted( fitols1p ) )
print( fitted( fitols1p$eq[[ 2 ]] ) )

print( fitted( fitols2r ) )
print( fitted( fitols2r$eq[[ 1 ]] ) )

print( fitted( fitols3s ) )
print( fitted( fitols3s$eq[[ 2 ]] ) )

print( fitted( fitols4rs ) )
print( fitted( fitols4rs$eq[[ 1 ]] ) )

print( fitted( fitols5 ) )
print( fitted( fitols5$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitols1p, se.fit = TRUE, interval = "prediction",
   probDfSys = TRUE ) )
print( predict( fitols1p$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   probDfSys = TRUE ) )

print( predict( fitols2r, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitols2r$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitols3s, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitols3s$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitols4rs, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitols4rs$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )

print( predict( fitols5, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitols5$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitols1p, newdata = smallData ) )
print( predict( fitols1p$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitols2r, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitols2r$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitols3s, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitols3s$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitols4rs, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitols4rs$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitols5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitols5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitols5rs, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitols5rs$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fitols1, restrm ) )
print( ftest.systemfit( fitols1s, restrm ) )
print( ftest.systemfit( fitols1p, restrm ) )
print( ftest.systemfit( fitols1r, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fitols1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fitols2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitols3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fitols1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fitols1, restrm ) )
print( waldtest.systemfit( fitols1s, restrm ) )
print( waldtest.systemfit( fitols1p, restrm ) )
print( waldtest.systemfit( fitols1r, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fitols1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fitols2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitols3, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fitols1, restr2m, restr2q ) )


## ****************** model frame **************************
print( mf <- model.frame( fitols1p ) )
print( mf1 <- model.frame( fitols1p$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitols1p$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitols2r ) ) )
print( all.equal( mf1, model.frame( fitols2r$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitols3s ) ) )
print( all.equal( mf2, model.frame( fitols3s$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitols4rs ) ) )
print( all.equal( mf1, model.frame( fitols4rs$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitols5 ) ) )
print( all.equal( mf2, model.frame( fitols5$eq[[ 2 ]] ) ) )


## **************** model matrix ************************
print( mm <- model.matrix( fitols1r ) )
print( mm1 <- model.matrix( fitols1r$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitols1r$eq[[ 2 ]] ) )
fitols1r$eq[[ 1 ]]$modelMatrix <- NULL
fitols1r$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitols1r ) ) )
print( all.equal( mm1, model.matrix( fitols1r$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols1r$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitols2rs ) ) )
print( all.equal( mm1, model.matrix( fitols2rs$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols2rs$eq[[ 2 ]] ) ) )
fitols2rs$eq[[ 1 ]]$modelMatrix <- NULL
fitols2rs$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitols2rs ) ) )
print( all.equal( mm1, model.matrix( fitols2rs$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols2rs$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitols3 ) ) )
print( all.equal( mm1, model.matrix( fitols3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols3$eq[[ 2 ]] ) ) )
fitols3$eq[[ 1 ]]$modelMatrix <- NULL
fitols3$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitols3 ) ) )
print( all.equal( mm1, model.matrix( fitols3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols3$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitols4p ) ) )
print( all.equal( mm1, model.matrix( fitols4p$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols4p$eq[[ 2 ]] ) ) )
fitols4p$eq[[ 1 ]]$modelMatrix <- NULL
fitols4p$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitols4p ) ) )
print( all.equal( mm1, model.matrix( fitols4p$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols4p$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitols5s ) ) )
print( all.equal( mm1, model.matrix( fitols5s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols5s$eq[[ 2 ]] ) ) )
fitols5s$eq[[ 1 ]]$modelMatrix <- NULL
fitols5s$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitols5s ) ) )
print( all.equal( mm1, model.matrix( fitols5s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols5s$eq[[ 2 ]] ) ) )
