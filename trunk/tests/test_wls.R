
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

## *************** WLS estimation (EViews-like) ************************
fitwls1e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor" )
print( summary( fitwls1e, useDfSys = TRUE ) )

## ************** WLS with cross-equation restriction ***************
fitwls2 <- systemfit( system, "WLS", data = Kmenta, restrictions = restrm )
print( summary( fitwls2 ) )

## ************** WLS with cross-equation restriction (EViews-like) *******
fitwls2e <- systemfit( system, "WLS", data = Kmenta, restrictions = restrm,
   methodRCov = "noDfCor" )
print( summary( fitwls2e ) )

## ******* WLS with cross-equation restriction via TX **********
fitwls3 <- systemfit( system,"WLS", data = Kmenta, TX = tc,)
print( summary( fitwls3 ) )

## ******* WLS with cross-equation restriction via TX (EViews-like) *****
fitwls3e <- systemfit( system,"WLS", data = Kmenta, TX = tc,
   methodRCov = "noDfCor" )
print( summary( fitwls3e ) )

## ***** WLS with 2 cross-equation restrictions ***************
fitwls4 <- systemfit( system,"WLS", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q )
print( summary( fitwls4 ) )

## ***** WLS with 2 cross-equation restrictions (EViews-like) **********
fitwls4e <- systemfit( system,"WLS", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr2m, restrict.rhs = restr2q )
print( summary( fitwls4e ) )

## *********** WLS with 2 cross-equation restrictions via R and TX ******
fitwls5 <- systemfit( system, "WLS", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc )
print( summary( fitwls5 ) )

## *********** WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwls5e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr3m, restrict.rhs = restr3q, TX = tc )
print( summary( fitwls5e ) )

## *************** iterated WLS estimation *********************
fitwlsi1 <- systemfit( system, "WLS", data = Kmenta,
   maxit = 100 )
print( summary( fitwlsi1, useDfSys = TRUE ) )

## *************** iterated WLS estimation (EViews-like) ************
fitwlsi1e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   maxit = 100 )
print( summary( fitwlsi1e, useDfSys = TRUE ) )

## ****** iterated WLS with cross-equation restriction ***************
fitwlsi2 <- systemfit( system, "WLS", data = Kmenta, restrictions = restrm,
   maxit = 100 )
print( summary( fitwlsi2 ) )

## ****** iterated WLS with cross-equation restriction (EViews-like) ********
fitwlsi2e <- systemfit( system, "WLS", data = Kmenta, restrictions = restrm,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitwlsi2e ) )

## ******* iterated WLS with cross-equation restriction via TX **********
fitwlsi3 <- systemfit( system, "WLS", data = Kmenta, TX = tc,
   maxit = 100 )
print( summary( fitwlsi3 ) )

## ******* iterated WLS with cross-equation restriction via TX (EViews-like) ***
fitwlsi3e <- systemfit( system, "WLS", data = Kmenta, TX = tc,
   methodRCov = "noDfCor", maxit = 100 )
print( summary( fitwlsi3e ) )

## ******* iterated WLS with 2 cross-equation restrictions ***********
fitwlsi4 <- systemfit( system, "WLS", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q, maxit = 100 )
print( summary( fitwlsi4 ) )

## ******* iterated WLS with 2 cross-equation restrictions (EViews-like) *****
fitwlsi4e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr2m, restrict.rhs = restr2q, maxit = 100 )
print( summary( fitwlsi4e ) )

## ***** iterated WLS with 2 cross-equation restrictions via R and TX ******
fitwlsi5 <- systemfit( system, "WLS", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5 ) )

## *** iterated WLS with 2 cross-equation restrictions via R and TX (EViews-like)
fitwlsi5e <- systemfit( system, "WLS", data = Kmenta, methodRCov = "noDfCor",
   restrictions = restr3m, restrict.rhs = restr3q, TX = tc, maxit = 100 )
print( summary( fitwlsi5e ) )


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


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitwls1e ), digits = 6 ) )
print( round( vcov( fitwls1e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwls2 ), digits = 6 ) )
print( round( vcov( fitwls2$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwls3e ), digits = 6 ) )
print( round( vcov( fitwls3e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwls4 ), digits = 6 ) )
print( round( vcov( fitwls4$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwls5 ), digits = 6 ) )
print( round( vcov( fitwls5$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi1 ), digits = 6 ) )
print( round( vcov( fitwlsi1$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi2e ), digits = 6 ) )
print( round( vcov( fitwlsi2e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi3 ), digits = 6 ) )
print( round( vcov( fitwlsi3$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi4e ), digits = 6 ) )
print( round( vcov( fitwlsi4e$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitwlsi5e ), digits = 6 ) )
print( round( vcov( fitwlsi5e$eq[[ 2 ]] ), digits = 6 ) )


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

print( confint( fitwlsi1e, level = 0.975, useDfSys = TRUE ) )
print( confint( fitwlsi1e$eq[[ 1 ]], level = 0.999, useDfSys = TRUE ) )

print( confint( fitwlsi2, level = 0.999 ) )
print( confint( fitwlsi2$eq[[ 2 ]], level = 0.1 ) )

print( confint( fitwlsi3e, level = 0.1 ) )
print( confint( fitwlsi3e$eq[[ 1 ]], level = 0.01 ) )

print( confint( fitwlsi4, level = 0.01 ) )
print( confint( fitwlsi4$eq[[ 2 ]], level = 0.33 ) )

print( confint( fitwlsi5e, level = 0.33 ) )
print( confint( fitwlsi5e$eq[[ 1 ]] ) )


## *********** fitted values *************
print( fitted( fitwls1 ) )
print( fitted( fitwls1$eq[[ 2 ]] ) )

print( fitted( fitwls2e ) )
print( fitted( fitwls2e$eq[[ 1 ]] ) )

print( fitted( fitwls3 ) )
print( fitted( fitwls3$eq[[ 2 ]] ) )

print( fitted( fitwls4e ) )
print( fitted( fitwls4e$eq[[ 1 ]] ) )

print( fitted( fitwls5 ) )
print( fitted( fitwls5$eq[[ 2 ]] ) )

print( fitted( fitwlsi1e ) )
print( fitted( fitwlsi1e$eq[[ 1 ]] ) )

print( fitted( fitwlsi2 ) )
print( fitted( fitwlsi2$eq[[ 2 ]] ) )

print( fitted( fitwlsi3e ) )
print( fitted( fitwlsi3e$eq[[ 1 ]] ) )

print( fitted( fitwlsi4 ) )
print( fitted( fitwlsi4$eq[[ 2 ]] ) )

print( fitted( fitwlsi5e ) )
print( fitted( fitwlsi5e$eq[[ 1 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitwls1, se.fit = TRUE, interval = "prediction" ) )
print( predict( fitwls1$eq[[ 2 ]], se.fit = TRUE, interval = "prediction" ) )

print( predict( fitwls2e, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitwls2e$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitwls3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitwls3$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fitwls4e, se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )
print( predict( fitwls4e$eq[[ 1 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25 ) )

print( predict( fitwls5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitwls5$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitwlsi1e, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )
print( predict( fitwlsi1e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99, useDfSys = TRUE ) )

print( predict( fitwlsi2, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )
print( predict( fitwlsi2$eq[[ 2 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData ) )

print( predict( fitwlsi3e, interval = "prediction", level = 0.925 ) )
print( predict( fitwlsi3e$eq[[ 1 ]], interval = "prediction", level = 0.925 ) )

print( predict( fitwlsi4, interval = "confidence", newdata = predictData ) )
print( predict( fitwlsi4$eq[[ 2 ]], interval = "confidence",
   newdata = predictData ) )

print( predict( fitwlsi5e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01 ) )
print( predict( fitwlsi5e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01 ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitwls1, newdata = smallData ) )
print( predict( fitwls1$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitwls2e, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitwls2e$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitwls3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitwls3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitwls4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitwls4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitwls5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitwls5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitwlsi3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitwlsi3e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitwls1, 2, 1 ) )

print( correlation.systemfit( fitwls2e, 1, 2 ) )

print( correlation.systemfit( fitwls3, 2, 1 ) )

print( correlation.systemfit( fitwls4e, 1, 2 ) )

print( correlation.systemfit( fitwls5, 2, 1 ) )

print( correlation.systemfit( fitwlsi1e, 1, 2 ) )

print( correlation.systemfit( fitwlsi2, 2, 1 ) )

print( correlation.systemfit( fitwlsi3e, 1, 2 ) )

print( correlation.systemfit( fitwlsi4, 2, 1 ) )

print( correlation.systemfit( fitwlsi5e, 1, 2 ) )


## ************ Log-Likelihood values ***************
print( logLik( fitwls1 ) )

print( logLik( fitwls2e ) )

print( logLik( fitwls3 ) )

print( logLik( fitwls4e ) )

print( logLik( fitwls5 ) )

print( logLik( fitwlsi1e ) )

print( logLik( fitwlsi2 ) )

print( logLik( fitwlsi3e ) )

print( logLik( fitwlsi4 ) )

print( logLik( fitwlsi5e ) )


## ************** F tests ****************
# testing first restriction
print( linear.hypothesis( fitwls1, restrm ) )
print( linear.hypothesis( fitwlsi1e, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( linear.hypothesis( fitwls1e, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitwlsi1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( linear.hypothesis( fitwls2, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitwls3, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitwlsi2e, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitwlsi3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( linear.hypothesis( fitwls1e, restr2m, restr2q ) )
print( linear.hypothesis( fitwlsi1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( linear.hypothesis( fitwls1, restrm, test = "Chisq" ) )
print( linear.hypothesis( fitwlsi1e, restrm, test = "Chisq" ) )

# testing second restriction
# first restriction not imposed
print( linear.hypothesis( fitwls1e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitwlsi1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
# first restriction imposed
print( linear.hypothesis( fitwls2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitwls3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitwlsi2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitwlsi3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )

# testing both of the restrictions
print( linear.hypothesis( fitwls1e, restr2m, restr2q, test = "Chisq" ) )
print( linear.hypothesis( fitwlsi1, restr2m, restr2q, test = "Chisq" ) )


## ****************** model frame **************************
print( mf <- model.frame( fitwls1 ) )
print( mf1 <- model.frame( fitwls1$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitwls1$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitwls2e ) ) )
print( all.equal( mf1, model.frame( fitwls2e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwls3 ) ) )
print( all.equal( mf2, model.frame( fitwls3$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwls4e ) ) )
print( all.equal( mf1, model.frame( fitwls4e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwls5 ) ) )
print( all.equal( mf2, model.frame( fitwls5$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi1e ) ) )
print( all.equal( mf1, model.frame( fitwlsi1e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi2 ) ) )
print( all.equal( mf2, model.frame( fitwlsi2$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi3e ) ) )
print( all.equal( mf1, model.frame( fitwlsi3e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi4 ) ) )
print( all.equal( mf2, model.frame( fitwlsi4$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitwlsi5e ) ) )
print( all.equal( mf1, model.frame( fitwlsi5e$eq[[ 1 ]] ) ) )


## **************** model matrix ************************
print( mm <- model.matrix( fitwlsi1e ) )
print( mm1 <- model.matrix( fitwlsi1e$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitwlsi1e$eq[[ 2 ]] ) )
fitwlsi1e$eq[[ 1 ]]$modelMatrix <- NULL
fitwlsi1e$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitwlsi1e ) ) )
print( all.equal( mm1, model.matrix( fitwlsi1e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi1e$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitwls2 ) ) )
print( all.equal( mm1, model.matrix( fitwls2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls2$eq[[ 2 ]] ) ) )
fitwls2$eq[[ 1 ]]$modelMatrix <- NULL
fitwls2$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitwls2 ) ) )
print( all.equal( mm1, model.matrix( fitwls2$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls2$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitwlsi3 ) ) )
print( all.equal( mm1, model.matrix( fitwlsi3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi3$eq[[ 2 ]] ) ) )
fitwlsi3$eq[[ 1 ]]$modelMatrix <- NULL
fitwlsi3$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitwlsi3 ) ) )
print( all.equal( mm1, model.matrix( fitwlsi3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwlsi3$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitwls4e ) ) )
print( all.equal( mm1, model.matrix( fitwls4e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls4e$eq[[ 2 ]] ) ) )
fitwls4e$eq[[ 1 ]]$modelMatrix <- NULL
fitwls4e$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitwls4e ) ) )
print( all.equal( mm1, model.matrix( fitwls4e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls4e$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitwls5 ) ) )
print( all.equal( mm1, model.matrix( fitwls5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls5$eq[[ 2 ]] ) ) )
fitwls5$eq[[ 1 ]]$modelMatrix <- NULL
fitwls5$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitwls5 ) ) )
print( all.equal( mm1, model.matrix( fitwls5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitwls5$eq[[ 2 ]] ) ) )
