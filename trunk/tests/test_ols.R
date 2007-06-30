
library( systemfit )
data( "Kmenta" )

demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
system <- list( demand = demand, supply = supply )
restrm <- matrix(0,1,7)  # restriction matrix "R"
restrm[1,3] <-  1
restrm[1,7] <- -1
restrict <- "demand_income - supply_trend = 0"
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )    # restriction vector "q" 2
restrict2 <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
tc <- matrix(0,7,6)
tc[1,1] <- 1
tc[2,2] <- 1
tc[3,3] <- 1
tc[4,4] <- 1
tc[5,5] <- 1
tc[6,6] <- 1
tc[7,3] <- 1
restr3m <- matrix(0,1,6)  # restriction matrix "R" 2
restr3m[1,2] <- -1
restr3m[1,5] <-  1
restr3q <- c( 0.5 )       # restriction vector "q" 2
restrict3 <- "- C2 + C5 = 0.5"

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

## ****************  OLS (useDfSys=T) ***********************
fitols1p <- systemfit( system, "OLS", data = Kmenta )
print( summary( fitols1p, useDfSys = TRUE ) )

## ****************  OLS (methodRCov="noDfCor") ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "noDfCor", returnModelMatrix = TRUE )
print( summary( fitols1r ) )

## ********  OLS (methodRCov="noDfCor", single.eq.sigma=F) ***********
fitols1rs <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "noDfCor", single.eq.sigma = FALSE )
print( summary( fitols1rs ) )

## ****************  OLS (methodRCov="Theil" ) ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "Theil", returnModelMatrix = TRUE )
print( summary( fitols1r ) )

## ****************  OLS (methodRCov="max") ***********************
fitols1r <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "max", returnModelMatrix = TRUE )
print( summary( fitols1r ) )

## ********  OLS (methodRCov="max", single.eq.sigma=F) ***********
fitols1rs <- systemfit( system, "OLS", data = Kmenta,
   methodRCov = "max", single.eq.sigma = FALSE )
print( summary( fitols1rs ) )


## ********* OLS with cross-equation restriction ************
## ****** OLS with cross-equation restriction (default) *********
fitols2 <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrm )
print( summary( fitols2 ) )
# the same with symbolically specified restrictions
fitols2Sym <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrict )
all.equal( fitols2, fitols2Sym )

## ****** OLS with cross-equation restriction (single.eq.sigma=T) *******
fitols2s <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrm, single.eq.sigma = TRUE )
print( summary( fitols2s ) )

## ****** OLS with cross-equation restriction (useDfSys=F) *******
fitols2p <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrm )
print( summary( fitols2p, useDfSys = FALSE ) )

## ****** OLS with cross-equation restriction (methodRCov="noDfCor") *******
fitols2r <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrm, methodRCov = "noDfCor" )
print( summary( fitols2r ) )

## ** OLS with cross-equation restriction (methodRCov="noDfCor",single.eq.sigma=T) ***
fitols2rs <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrm, methodRCov = "noDfCor", returnModelMatrix = TRUE )
print( summary( fitols2rs ) )

## *** OLS with cross-equation restriction via restrict.regMat ***
## *** OLS with cross-equation restriction via restrict.regMat (default) ***
fitols3 <- systemfit( system, "OLS", data = Kmenta, restrict.regMat = tc,
   returnModelMatrix = TRUE )
print( summary( fitols3 ) )

## *** OLS with cross-equation restriction via restrict.regMat (single.eq.sigma=T) ***
fitols3s <- systemfit( system, "OLS", data = Kmenta,
   restrict.regMat = tc, single.eq.sigma = TRUE )
print( summary( fitols3s ) )

## *** OLS with cross-equation restriction via restrict.regMat (useDfSys=F) ***
fitols3p <- systemfit( system, "OLS", data = Kmenta,
   restrict.regMat = tc )
print( summary( fitols3p, useDfSys = FALSE ) )

## *** OLS with cross-equation restriction via restrict.regMat (methodRCov="noDfCor") ***
fitols3r <- systemfit( system, "OLS", data = Kmenta,
   restrict.regMat = tc, methodRCov = "noDfCor" )
print( summary( fitols3r ) )

## OLS with cross-equation restriction via restrict.regMat (methodRCov="noDfCor",single.eq.sigma=T)
fitols3rs <- systemfit( system, "OLS", data = Kmenta,
   restrict.regMat = tc, methodRCov = "noDfCor", single.eq.sigma = TRUE )
print( summary( fitols3rs ) )

## ********* OLS with 2 cross-equation restrictions ***********
## ********* OLS with 2 cross-equation restrictions (default) ***********
fitols4 <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q )
print( summary( fitols4 ) )
# the same with symbolically specified restrictions
fitols4Sym <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrict2 )
all.equal( fitols4, fitols4Sym )

## ****** OLS with 2 cross-equation restrictions (single.eq.sigma=T) *******
fitols4s <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, single.eq.sigma = T )
print( summary( fitols4s ) )

## ****** OLS with 2 cross-equation restrictions (useDfSys=F) *******
fitols4p <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, returnModelMatrix = TRUE )
print( summary( fitols4p, useDfSys = FALSE ) )

## ****** OLS with 2 cross-equation restrictions (methodRCov="noDfCor") *******
fitols4r <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, methodRCov = "noDfCor" )
print( summary( fitols4r ) )

## OLS with 2 cross-equation restrictions (methodRCov="noDfCor", single.eq.sigma=T) *
fitols4rs <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr2m,
   restrict.rhs = restr2q, methodRCov = "noDfCor", single.eq.sigma = T )
print( summary( fitols4rs ) )

## ***** OLS with 2 cross-equation restrictions via R and restrict.regMat ****
## ***** OLS with 2 cross-equation restrictions via R and restrict.regMat (default) ****
fitols5 <- systemfit( system, "OLS", data = Kmenta, restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, methodRCov = "noDfCor")
print( summary( fitols5 ) )
# the same with symbolically specified restrictions
fitols5Sym <- systemfit( system, "OLS", data = Kmenta,
   restrict.matrix = restrict3, restrict.regMat = tc, methodRCov = "noDfCor")
all.equal( fitols5, fitols5Sym )

## ***** OLS with 2 cross-equation restrictions via R and restrict.regMat (single.eq.sigma=T) ****
fitols5s <- systemfit( system, "OLS", data = Kmenta,restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, single.eq.sigma = T,
   returnModelMatrix = TRUE )
print( summary( fitols5s ) )

## ***** OLS with 2 cross-equation restrictions via R and restrict.regMat (useDfSys=F) ****
fitols5p <- systemfit( system, "OLS", data = Kmenta,restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc )
print( summary( fitols5p, useDfSys = FALSE ) )

## ***** OLS with 2 cross-equation restrictions via R and restrict.regMat (methodRCov="noDfCor") ****
fitols5r <- systemfit( system, "OLS", data = Kmenta,restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, methodRCov = "noDfCor" )
print( summary( fitols5r ) )

## OLS with 2 cross-equation restr. via R and restrict.regMat (methodRCov="noDfCor",single.eq.sigma=T)
fitols5rs <- systemfit( system, "OLS", data = Kmenta,restrict.matrix = restr3m,
   restrict.rhs = restr3q, restrict.regMat = tc, methodRCov = "noDfCor", single.eq.sigma = T )
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


## *************** coefficients *********************
print( round( coef( fitols1rs ), digits = 6 ) )
print( round( coef( fitols1rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitols2s ), digits = 6 ) )
print( round( coef( fitols2s$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitols3p ), digits = 6 ) )
print( round( coef( fitols3p, modified.reg = TRUE ), digits = 6 ) )
print( round( coef( fitols3p$eq[[ 2 ]] ), digits = 6 ) )

print( round( coef( fitols4r ), digits = 6 ) )
print( round( coef( fitols4r$eq[[ 1 ]] ), digits = 6 ) )

print( round( coef( fitols5 ), digits = 6 ) )
print( round( coef( fitols5, modified.reg = TRUE ), digits = 6 ) )
print( round( coef( fitols5$eq[[ 2 ]] ), digits = 6 ) )


## *************** coefficients with stats *********************
print( round( coef( summary( fitols1rs, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitols1rs$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fitols2s ) ), digits = 6 ) )
print( round( coef( summary( fitols2s$eq[[ 1 ]] ) ), digits = 6 ) )

print( round( coef( summary( fitols3p, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitols3p, useDfSys = FALSE ), modified.reg = TRUE ),
   digits = 6 ) )
print( round( coef( summary( fitols3p$eq[[ 2 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fitols4r, useDfSys = FALSE ) ), digits = 6 ) )
print( round( coef( summary( fitols4r$eq[[ 1 ]], useDfSys = FALSE ) ),
   digits = 6 ) )

print( round( coef( summary( fitols5 ) ), digits = 6 ) )
print( round( coef( summary( fitols5 ), modified.reg = TRUE ), digits = 6 ) )
print( round( coef( summary( fitols5$eq[[ 2 ]] ) ), digits = 6 ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitols1rs ), digits = 6 ) )
print( round( vcov( fitols1rs$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols2s ), digits = 6 ) )
print( round( vcov( fitols2s$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols3p ), digits = 6 ) )
print( round( vcov( fitols3p, modified.reg = TRUE ), digits = 6 ) )
print( round( vcov( fitols3p$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitols4r ), digits = 6 ) )
print( round( vcov( fitols4r$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitols5 ), digits = 6 ) )
print( round( vcov( fitols5, modified.reg = TRUE ), digits = 6 ) )
print( round( vcov( fitols5$eq[[ 2 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitols1p, useDfSys = TRUE ) )
print( confint( fitols1p$eq[[ 2 ]], level = 0.9, useDfSys = TRUE ) )

print( confint( fitols2r, level = 0.9 ) )
print( confint( fitols2r$eq[[ 1 ]], level = 0.99 ) )

print( confint( fitols3s, level = 0.99 ) )
print( confint( fitols3s$eq[[ 2 ]], level = 0.5 ) )

print( confint( fitols4rs, level = 0.5 ) )
print( confint( fitols4rs$eq[[ 1 ]], level = 0.25 ) )

print( confint( fitols5, level = 0.25 ) )
print( confint( fitols5$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitols3p, level = 0.999, useDfSys = FALSE ) )
print( confint( fitols3p$eq[[ 1 ]], useDfSys = FALSE ) )


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
   useDfSys = TRUE ) )
print( predict( fitols1p$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )

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


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitols1p, 1, 2 ) )

print( correlation.systemfit( fitols2r, 2, 1 ) )

print( correlation.systemfit( fitols3s, 1, 2 ) )

print( correlation.systemfit( fitols4rs, 2, 1 ) )

print( correlation.systemfit( fitols5, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fitols1p ) )

print( logLik( fitols2r ) )

print( logLik( fitols3s ) )

print( logLik( fitols4rs ) )

print( logLik( fitols5 ) )


## ************** F tests ****************
# testing first restriction
print( linear.hypothesis( fitols1, restrm ) )
linear.hypothesis( fitols1, restrict )

print( linear.hypothesis( fitols1s, restrm ) )
linear.hypothesis( fitols1s, restrict )

print( linear.hypothesis( fitols1p, restrm ) )
linear.hypothesis( fitols1p, restrict )

print( linear.hypothesis( fitols1r, restrm ) )
linear.hypothesis( fitols1r, restrict )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
restrictOnly2 <- "- demand_price + supply_price = 0.5"
# first restriction not imposed 
print( linear.hypothesis( fitols1, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fitols1, restrictOnly2 )

# first restriction imposed
print( linear.hypothesis( fitols2, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fitols2, restrictOnly2 )

print( linear.hypothesis( fitols3, restrOnly2m, restrOnly2q ) )
linear.hypothesis( fitols3, restrictOnly2 )

# testing both of the restrictions
print( linear.hypothesis( fitols1, restr2m, restr2q ) )
linear.hypothesis( fitols1, restrict2 )


## ************** Wald tests ****************
# testing first restriction
print( linear.hypothesis( fitols1, restrm, test = "Chisq" ) )
linear.hypothesis( fitols1, restrict, test = "Chisq" )

print( linear.hypothesis( fitols1s, restrm, test = "Chisq" ) )
linear.hypothesis( fitols1s, restrict, test = "Chisq" )

print( linear.hypothesis( fitols1p, restrm, test = "Chisq" ) )
linear.hypothesis( fitols1p, restrict, test = "Chisq" )

print( linear.hypothesis( fitols1r, restrm, test = "Chisq" ) )
linear.hypothesis( fitols1r, restrict, test = "Chisq" )

# testing second restriction
# first restriction not imposed
print( linear.hypothesis( fitols1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fitols1, restrictOnly2, test = "Chisq" )
# first restriction imposed
print( linear.hypothesis( fitols2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fitols2, restrictOnly2, test = "Chisq" )

print( linear.hypothesis( fitols3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
linear.hypothesis( fitols3, restrictOnly2, test = "Chisq" )

# testing both of the restrictions
print( linear.hypothesis( fitols1, restr2m, restr2q, test = "Chisq" ) )
linear.hypothesis( fitols1, restrict2, test = "Chisq" )


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
# with returnModelMatrix = TRUE
print( !is.null( fitols1r$eq[[ 1 ]]$modelMatrix ) )
print( mm <- model.matrix( fitols1r ) )
print( mm1 <- model.matrix( fitols1r$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitols1r$eq[[ 2 ]] ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fitols1rs ) ) )
print( all.equal( mm1, model.matrix( fitols1rs$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols1rs$eq[[ 2 ]] ) ) )
print( !is.null( fitols1rs$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fitols2rs$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fitols2rs ) ) )
print( all.equal( mm1, model.matrix( fitols2rs$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols2rs$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fitols2p ) ) )
print( all.equal( mm1, model.matrix( fitols2p$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols2p$eq[[ 2 ]] ) ) )
print( !is.null( fitols2p$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fitols3$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fitols3 ) ) )
print( all.equal( mm1, model.matrix( fitols3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols3$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fitols3r ) ) )
print( all.equal( mm1, model.matrix( fitols3r$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols3r$eq[[ 2 ]] ) ) )
print( !is.null( fitols3r$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fitols4p$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fitols4p ) ) )
print( all.equal( mm1, model.matrix( fitols4p$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols4p$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fitols4Sym ) ) )
print( all.equal( mm1, model.matrix( fitols4Sym$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols4Sym$eq[[ 2 ]] ) ) )
print( !is.null( fitols4Sym$eq[[ 1 ]]$modelMatrix ) )

# with returnModelMatrix = TRUE
print( !is.null( fitols5s$eq[[ 1 ]]$modelMatrix ) )
print( all.equal( mm, model.matrix( fitols5s ) ) )
print( all.equal( mm1, model.matrix( fitols5s$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols5s$eq[[ 2 ]] ) ) )

# with returnModelMatrix = FALSE
print( all.equal( mm, model.matrix( fitols5 ) ) )
print( all.equal( mm1, model.matrix( fitols5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitols5$eq[[ 2 ]] ) ) )
print( !is.null( fitols5$eq[[ 1 ]]$modelMatrix ) )


## **************** formulas ************************
formula( fitols1p )
formula( fitols1p$eq[[ 2 ]] )

formula( fitols2r )
formula( fitols2r$eq[[ 1 ]] )

formula( fitols3s )
formula( fitols3s$eq[[ 2 ]] )

formula( fitols4rs )
formula( fitols4rs$eq[[ 1 ]] )

formula( fitols5 )
formula( fitols5$eq[[ 2 ]] )


## **************** model terms *******************
terms( fitols1p )
terms( fitols1p$eq[[ 2 ]] )

terms( fitols2r )
terms( fitols2r$eq[[ 1 ]] )

terms( fitols3s )
terms( fitols3s$eq[[ 2 ]] )

terms( fitols4rs )
terms( fitols4rs$eq[[ 1 ]] )

terms( fitols5 )
terms( fitols5$eq[[ 2 ]] )
