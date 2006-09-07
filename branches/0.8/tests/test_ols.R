
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
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

# It is not possible to estimate OLS with systemfit
# exactly as EViews does, because EViews uses
# rcovformula == 1 for the coefficient covariance matrix and
# rcovformula == 0 for the residual covariance matrix, while
# systemfit uses always the same formulas for both calculations.

## *************** OLS estimation ************************
## ********** OLS estimation (default) ********************
fitols1 <- systemfit( "OLS", system, labels, data = Kmenta )
print( summary( fitols1 ) )
print( round( fitols1$bcov, digits = 6 ) )

## ********** OLS estimation (no single.eq sigma=F) ******************
fitols1s <- systemfit( "OLS", system, labels, data = Kmenta,
   single.eq.sigma = FALSE )
print( summary( fitols1s ) )
print( round( fitols1s$bcov, digits = 6 ) )

## ****************  OLS (probdfsys=T) ***********************
fitols1p <- systemfit( "OLS", system, labels, data = Kmenta,
   probdfsys = TRUE )
print( summary( fitols1p ) )
print( round( fitols1p$bcov, digits = 6 ) )

## ****************  OLS (rcovformula=0) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 0 )
print( summary( fitols1r ) )
print( round( fitols1r$bcov, digits = 6 ) )

## ********  OLS (rcovformula=0, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 0, single.eq.sigma = FALSE )
print( summary( fitols1rs ) )
print( round( fitols1rs$bcov, digits = 6 ) )

## ****************  OLS (rcovformula=2 ) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 2 )
print( summary( fitols1r ) )
print( round( fitols1r$bcov, digits = 6 ) )

## ****************  OLS (rcovformula=3) ***********************
fitols1r <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 3 )
print( summary( fitols1r ) )
print( round( fitols1r$bcov, digits = 6 ) )

## ********  OLS (rcovformula=3, single.eq.sigma=F) ***********
fitols1rs <- systemfit( "OLS", system, labels, data = Kmenta,
   rcovformula = 3, single.eq.sigma = FALSE )
print( summary( fitols1rs ) )
print( round( fitols1rs$bcov, digits = 6 ) )


## ********* OLS with cross-equation restriction ************
## ****** OLS with cross-equation restriction (default) *********
fitols2 <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm )
print( summary( fitols2 ) )
print( round( fitols2$bcov, digits = 6 ) )

## ****** OLS with cross-equation restriction (single.eq.sigma=T) *******
fitols2s <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, single.eq.sigma = TRUE )
print( summary( fitols2s ) )
print( round( fitols2s$bcov, digits = 6 ) )

## ****** OLS with cross-equation restriction (probdfsys=F) *******
fitols2p <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, probdfsys = FALSE )
print( summary( fitols2p ) )
print( round( fitols2p$bcov, digits = 6 ) )

## ****** OLS with cross-equation restriction (rcovformula=0) *******
fitols2r <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, rcovformula = 0 )
print( summary( fitols2r ) )
print( round( fitols2r$bcov, digits = 6 ) )

## ** OLS with cross-equation restriction (rcovformula=0,single.eq.sigma=T) ***
fitols2rs <- systemfit( "OLS", system, labels, data = Kmenta,
   R.restr = restrm, rcovformula = 0 )
print( summary( fitols2rs ) )
print( round( fitols2rs$bcov, digits = 6 ) )

## *** OLS with cross-equation restriction via TX ***
## *** OLS with cross-equation restriction via TX (default) ***
fitols3 <- systemfit( "OLS", system, labels, data = Kmenta, TX = tc )
print( summary( fitols3 ) )
print( round( fitols3$bcov, digits = 6 ) )

## *** OLS with cross-equation restriction via TX (single.eq.sigma=T) ***
fitols3s <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, single.eq.sigma = TRUE )
print( summary( fitols3s ) )
print( round( fitols3s$bcov, digits = 6 ) )

## *** OLS with cross-equation restriction via TX (probdfsys=F) ***
fitols3p <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, probdfsys = FALSE )
print( summary( fitols3p ) )
print( round( fitols3p$bcov, digits = 6 ) )

## *** OLS with cross-equation restriction via TX (rcovformula=0) ***
fitols3r <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, rcovformula = 0 )
print( summary( fitols3r ) )
print( round( fitols3r$bcov, digits = 6 ) )

## OLS with cross-equation restriction via TX (rcovformula=0,single.eq.sigma=T)
fitols3rs <- systemfit( "OLS", system, labels, data = Kmenta,
   TX = tc, rcovformula = 0, single.eq.sigma = TRUE )
print( summary( fitols3rs ) )
print( round( fitols3rs$bcov, digits = 6 ) )

## ********* OLS with 2 cross-equation restrictions ***********
## ********* OLS with 2 cross-equation restrictions (default) ***********
fitols4 <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q )
print( summary( fitols4 ) )
print( round( fitols4$bcov, digits = 6 ) )

## ****** OLS with 2 cross-equation restrictions (single.eq.sigma=T) *******
fitols4s <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, single.eq.sigma = T )
print( summary( fitols4s ) )
print( round( fitols4s$bcov, digits = 6 ) )

## ****** OLS with 2 cross-equation restrictions (probdfsys=F) *******
fitols4p <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, probdfsys = F )
print( summary( fitols4p ) )
print( round( fitols4p$bcov, digits = 6 ) )

## ****** OLS with 2 cross-equation restrictions (rcovformula=0) *******
fitols4r <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0 )
print( summary( fitols4r ) )
print( round( fitols4r$bcov, digits = 6 ) )

## OLS with 2 cross-equation restrictions (rcovformula=0, single.eq.sigma=T) *
fitols4rs <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, rcovformula = 0, single.eq.sigma = T )
print( summary( fitols4rs ) )
print( round( fitols4rs$bcov, digits = 6 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX ****
## ***** OLS with 2 cross-equation restrictions via R and TX (default) ****
fitols5 <- systemfit( "OLS", system, labels, data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0)
print( summary( fitols5 ) )
print( round( fitols5$bcov, digits = 6 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (single.eq.sigma=T) ****
fitols5s <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, single.eq.sigma = T )
print( summary( fitols5s ) )
print( round( fitols5s$bcov, digits = 6 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (probdfsys=F) ****
fitols5p <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, probdfsys = F )
print( summary( fitols5p ) )
print( round( fitols5p$bcov, digits = 6 ) )

## ***** OLS with 2 cross-equation restrictions via R and TX (rcovformula=0) ****
fitols5r <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0 )
print( summary( fitols5r ) )
print( round( fitols5r$bcov, digits = 6 ) )

## OLS with 2 cross-equation restr. via R and TX (rcovformula=0,single.eq.sigma=T)
fitols5rs <- systemfit( "OLS", system, labels, data = Kmenta,R.restr = restr3m,
   q.restr = restr3q, TX = tc, rcovformula = 0, single.eq.sigma = T )
print( summary( fitols5rs ) )
print( round( fitols5rs$bcov, digits = 6 ) )


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


## *********** confidence intervals of coefficients *************
print( confint( fitols1p ) )
print( confint( fitols1p$eq[[ 2 ]], level = 0.9 ) )

print( confint( fitols2r, level = 0.9 ) )
print( confint( fitols2r$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitols3s, level = 0.99 ) )
print( confint( fitols3s$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitols4rs, level = 0.5 ) )
print( confint( fitols4rs$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitols5, level = 0.25 ) )
print( confint( fitols5$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitols3p, level = 0.999 ) )
print( confint( fitols3p$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitols1p, se.fit = TRUE, interval = "prediction" ) )
print( predict( fitols1p$eq[[ 2 ]] ) )

print( predict( fitols2r, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fitols2r$eq[[ 1 ]] ) )

print( predict( fitols3s, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fitols3s$eq[[ 2 ]] ) )

print( predict( fitols4rs, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitols4rs$eq[[ 1 ]] ) )

print( predict( fitols5, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData ) )
print( predict( fitols5$eq[[ 2 ]] ) )


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
