
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
restrict <- "demand_income - supply_trend = 0"
restr2m <- matrix(0,2,7)  # restriction matrix "R" 2
restr2m[1,3] <-  1
restr2m[1,7] <- -1
restr2m[2,2] <- -1
restr2m[2,5] <-  1
restr2q <- c( 0, 0.5 )  # restriction vector "q" 2
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
restr3q <- c( 0.5 )  # restriction vector "q" 2
restrict3 <- "demand_income - supply_price = 0"


## ********************* W2SLS *****************
fitw2sls1 <- systemfit( system, "W2SLS", data = Kmenta, inst = inst )
print( summary( fitw2sls1 ) )

## ********************* W2SLS (EViews-like) *****************
fitw2sls1e <- systemfit( system, "W2SLS", data = Kmenta, inst = inst,
   methodRCov = "noDfCor" )
print( summary( fitw2sls1e, useDfSys = TRUE ) )

## ********************* W2SLS with restriction *******************
fitw2sls2 <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restrm,
   inst = inst )
print( summary( fitw2sls2 ) )
# the same with symbolically specified restrictions
fitw2sls2Sym <- systemfit( system, "W2SLS", data = Kmenta,
   restrictions = restrict, inst = inst )
all.equal( fitw2sls2, fitw2sls2Sym )

## ********************* W2SLS with restriction (EViews-like) **************
fitw2sls2e <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restrm,
   inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls2e, useDfSys = TRUE ) )

## ********************* W2SLS with restriction via TX *******************
fitw2sls3 <- systemfit( system, "W2SLS", data = Kmenta, TX = tc, inst = inst )
print( summary( fitw2sls3 ) )

## ********************* W2SLS with restriction via TX (EViews-like) **************
fitw2sls3e <- systemfit( system, "W2SLS", data = Kmenta, TX = tc, inst = inst,
   methodRCov = "noDfCor" )
print( summary( fitw2sls3e, useDfSys = TRUE ) )

## ***************** W2SLS with 2 restrictions ********************
fitw2sls4 <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q, inst = inst )
print( summary( fitw2sls4 ) )
# the same with symbolically specified restrictions
fitw2sls4Sym <- systemfit( system, "W2SLS", data = Kmenta,
   restrictions = restrict2, inst = inst )
all.equal( fitw2sls4, fitw2sls4Sym )

## ***************** W2SLS with 2 restrictions (EViews-like) **************
fitw2sls4e <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restr2m,
   restrict.rhs = restr2q, inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls4e, useDfSys = TRUE ) )

## ***************** W2SLS with 2 restrictions via R and TX ******************
fitw2sls5 <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc, inst = inst )
print( summary( fitw2sls5 ) )

## ***************** W2SLS with 2 restrictions via R and TX (EViews-like) **************
fitw2sls5e <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restr3m,
   restrict.rhs = restr3q, TX = tc, inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls5e, useDfSys = TRUE ) )

## ****** 2SLS estimation with different instruments **********************
fitw2slsd1 <- systemfit( system, "W2SLS", data = Kmenta, inst = instlist )
print( summary( fitw2slsd1 ) )

## ****** 2SLS estimation with different instruments (EViews-like)******************
fitw2slsd1e <- systemfit( system, "W2SLS", data = Kmenta, inst = instlist,
   methodRCov = "noDfCor" )
print( summary( fitw2slsd1e, useDfSys = TRUE ) )

## **** W2SLS estimation with different instruments and restriction ********
fitw2slsd2 <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restrm,
   inst = instlist )
print( summary( fitw2slsd2 ) )

## **** W2SLS estimation with different instruments and restriction (EViews-like)*
fitw2slsd2e <- systemfit( system, "W2SLS", data = Kmenta, restrictions = restrm,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fitw2slsd2e, useDfSys = TRUE ) )

## ** W2SLS estimation with different instruments and restriction via TX ****
fitw2slsd3 <- systemfit( system, "W2SLS", data = Kmenta, TX = tc,
   inst = instlist)
print( summary( fitw2slsd3 ) )

## W2SLS estimation with different instruments and restriction via TX (EViews-like)
fitw2slsd3e <- systemfit( system, "W2SLS", data = Kmenta, TX = tc,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fitw2slsd3e, useDfSys = TRUE ) )


## ****************** residuals **************************
print( residuals( fitw2sls1e ) )
print( residuals( fitw2sls1e$eq[[ 1 ]] ) )

print( residuals( fitw2sls2 ) )
print( residuals( fitw2sls2$eq[[ 2 ]] ) )

print( residuals( fitw2sls3 ) )
print( residuals( fitw2sls3$eq[[ 1 ]] ) )

print( residuals( fitw2sls4e ) )
print( residuals( fitw2sls4e$eq[[ 2 ]] ) )

print( residuals( fitw2sls5 ) )
print( residuals( fitw2sls5$eq[[ 1 ]] ) )

print( residuals( fitw2slsd1 ) )
print( residuals( fitw2slsd1$eq[[ 2 ]] ) )

print( residuals( fitw2slsd2e ) )
print( residuals( fitw2slsd2e$eq[[ 1 ]] ) )

print( residuals( fitw2slsd3e ) )
print( residuals( fitw2slsd3e$eq[[ 2 ]] ) )


## *********** variance covariance matrix of the coefficients *******
print( round( vcov( fitw2sls1e ), digits = 6 ) )
print( round( vcov( fitw2sls1e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls2 ), digits = 6 ) )
print( round( vcov( fitw2sls2$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls3e ), digits = 6 ) )
print( round( vcov( fitw2sls3e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls4 ), digits = 6 ) )
print( round( vcov( fitw2sls4$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2sls5 ), digits = 6 ) )
print( round( vcov( fitw2sls5$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd1 ), digits = 6 ) )
print( round( vcov( fitw2slsd1$eq[[ 1 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd2e ), digits = 6 ) )
print( round( vcov( fitw2slsd2e$eq[[ 2 ]] ), digits = 6 ) )

print( round( vcov( fitw2slsd3 ), digits = 6 ) )
print( round( vcov( fitw2slsd3$eq[[ 1 ]] ), digits = 6 ) )


## *********** confidence intervals of coefficients *************
print( confint( fitw2sls1e, useDfSys = TRUE ) )
print( confint( fitw2sls1e$eq[[ 1 ]], level = 0.9, useDfSys = TRUE ) )

print( confint( fitw2sls2, level = 0.9 ) )
print( confint( fitw2sls2$eq[[ 2 ]], level = 0.99 ) )

print( confint( fitw2sls3, level = 0.99 ) )
print( confint( fitw2sls3$eq[[ 1 ]], level = 0.5 ) )

print( confint( fitw2sls4e, level = 0.5, useDfSys = TRUE ) )
print( confint( fitw2sls4e$eq[[ 2 ]], level = 0.25, useDfSys = TRUE ) )

print( confint( fitw2sls5, level = 0.25 ) )
print( confint( fitw2sls5$eq[[ 1 ]], level = 0.975 ) )

print( confint( fitw2slsd1, level = 0.975 ) )
print( confint( fitw2slsd1$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitw2slsd2e, level = 0.999, useDfSys = TRUE ) )
print( confint( fitw2slsd2e$eq[[ 1 ]], level = 0.01, useDfSys = TRUE ) )

print( confint( fitw2slsd3e, level = 0.01, useDfSys = TRUE ) )
print( confint( fitw2slsd3e$eq[[ 2 ]], useDfSys = TRUE ) )


## *********** fitted values *************
print( fitted( fitw2sls1e ) )
print( fitted( fitw2sls1e$eq[[ 1 ]] ) )

print( fitted( fitw2sls2 ) )
print( fitted( fitw2sls2$eq[[ 2 ]] ) )

print( fitted( fitw2sls3 ) )
print( fitted( fitw2sls3$eq[[ 1 ]] ) )

print( fitted( fitw2sls4e ) )
print( fitted( fitw2sls4e$eq[[ 2 ]] ) )

print( fitted( fitw2sls5 ) )
print( fitted( fitw2sls5$eq[[ 1 ]] ) )

print( fitted( fitw2slsd1 ) )
print( fitted( fitw2slsd1$eq[[ 2 ]] ) )

print( fitted( fitw2slsd2e ) )
print( fitted( fitw2slsd2e$eq[[ 1 ]] ) )

print( fitted( fitw2slsd3e ) )
print( fitted( fitw2slsd3e$eq[[ 2 ]] ) )


## *********** predicted values *************
predictData <- Kmenta
predictData$consump <- NULL
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitw2sls1e, se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )
print( predict( fitw2sls1e$eq[[ 1 ]], se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE ) )

print( predict( fitw2sls2, se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )
print( predict( fitw2sls2$eq[[ 2 ]], se.pred = TRUE, interval = "confidence",
   level = 0.999, newdata = predictData ) )

print( predict( fitw2sls3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitw2sls3$eq[[ 1 ]], se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )

print( predict( fitw2sls4e, se.fit = TRUE, interval = "confidence",
   level = 0.25, useDfSys = TRUE ) )
print( predict( fitw2sls4e$eq[[ 2 ]], se.fit = TRUE, interval = "confidence",
   level = 0.25, useDfSys = TRUE ) )

print( predict( fitw2sls5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )
print( predict( fitw2sls5$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = predictData ) )

print( predict( fitw2slsd1, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitw2slsd1$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )

print( predict( fitw2slsd2e, se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData, useDfSys = TRUE ) )
print( predict( fitw2slsd2e$eq[[ 1 ]], se.fit = TRUE, interval = "prediction",
   level = 0.9, newdata = predictData, useDfSys = TRUE ) )

print( predict( fitw2slsd3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01, useDfSys = TRUE ) )
print( predict( fitw2slsd3e$eq[[ 2 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01, useDfSys = TRUE ) )

# predict just one observation
smallData <- data.frame( price = 130, income = 150, farmPrice = 120,
   trend = 25 )

print( predict( fitw2sls1e, newdata = smallData ) )
print( predict( fitw2sls1e$eq[[ 1 ]], newdata = smallData ) )

print( predict( fitw2sls2, se.fit = TRUE, level = 0.9,
   newdata = smallData ) )
print( predict( fitw2sls2$eq[[ 1 ]], se.pred = TRUE, level = 0.99,
   newdata = smallData ) )

print( predict( fitw2sls3, interval = "prediction", level = 0.975,
   newdata = smallData ) )
print( predict( fitw2sls3$eq[[ 1 ]], interval = "confidence", level = 0.8,
   newdata = smallData ) )

print( predict( fitw2sls4e, se.fit = TRUE, interval = "confidence",
   level = 0.999, newdata = smallData ) )
print( predict( fitw2sls4e$eq[[ 2 ]], se.pred = TRUE, interval = "prediction",
   level = 0.75, newdata = smallData ) )

print( predict( fitw2sls5, se.fit = TRUE, interval = "prediction",
   newdata = smallData ) )
print( predict( fitw2sls5$eq[[ 1 ]], se.pred = TRUE, interval = "confidence",
   newdata = smallData ) )

print( predict( fitw2slsd2e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, newdata = smallData ) )
print( predict( fitw2slsd2e$eq[[ 1 ]], se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.25, newdata = smallData ) )


## ************ correlation of predicted values ***************
print( correlation.systemfit( fitw2sls1e, 1, 2 ) )

print( correlation.systemfit( fitw2sls2, 2, 1 ) )

print( correlation.systemfit( fitw2sls3, 1, 2 ) )

print( correlation.systemfit( fitw2sls4e, 2, 1 ) )

print( correlation.systemfit( fitw2sls5, 1, 2 ) )

print( correlation.systemfit( fitw2slsd1, 2, 1 ) )

print( correlation.systemfit( fitw2slsd2e, 1, 2 ) )

print( correlation.systemfit( fitw2slsd3e, 2, 1 ) )


## ************ LOG-Likelihood values ***************
print( logLik( fitw2sls1e ) )

print( logLik( fitw2sls2 ) )

print( logLik( fitw2sls3 ) )

print( logLik( fitw2sls4e ) )

print( logLik( fitw2sls5 ) )

print( logLik( fitw2slsd1 ) )

print( logLik( fitw2slsd2e ) )

print( logLik( fitw2slsd3e ) )


## ************** F tests ****************
# testing first restriction
print( linear.hypothesis( fitw2sls1, restrm ) )
print( linear.hypothesis( fitw2slsd1e, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( linear.hypothesis( fitw2sls1e, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitw2slsd1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( linear.hypothesis( fitw2sls2, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitw2sls3, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitw2slsd2e, restrOnly2m, restrOnly2q ) )
print( linear.hypothesis( fitw2slsd3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( linear.hypothesis( fitw2sls1e, restr2m, restr2q ) )
print( linear.hypothesis( fitw2slsd1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( linear.hypothesis( fitw2sls1, restrm, test = "Chisq" ) )
print( linear.hypothesis( fitw2slsd1e, restrm, test = "Chisq" ) )

# testing second restriction
# first restriction not imposed
print( linear.hypothesis( fitw2sls1e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitw2slsd1, restrOnly2m, restrOnly2q, test = "Chisq" ) )
# first restriction imposed
print( linear.hypothesis( fitw2sls2, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitw2sls3, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitw2slsd2e, restrOnly2m, restrOnly2q, test = "Chisq" ) )
print( linear.hypothesis( fitw2slsd3e, restrOnly2m, restrOnly2q, test = "Chisq" ) )

# testing both of the restrictions
print( linear.hypothesis( fitw2sls1e, restr2m, restr2q, test = "Chisq" ) )
print( linear.hypothesis( fitw2slsd1, restr2m, restr2q, test = "Chisq" ) )


## ****************** model frame **************************
print( mf <- model.frame( fitw2sls1e ) )
print( mf1 <- model.frame( fitw2sls1e$eq[[ 1 ]] ) )
print( attributes( mf1 )$terms )
print( mf2 <- model.frame( fitw2sls1e$eq[[ 2 ]] ) )
print( attributes( mf2 )$terms )

print( all.equal( mf, model.frame( fitw2sls2 ) ) )
print( all.equal( mf2, model.frame( fitw2sls2$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitw2sls3 ) ) )
print( all.equal( mf1, model.frame( fitw2sls3$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitw2sls4e ) ) )
print( all.equal( mf2, model.frame( fitw2sls4e$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitw2sls5 ) ) )
print( all.equal( mf1, model.frame( fitw2sls5$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitw2slsd1 ) ) )
print( all.equal( mf2, model.frame( fitw2slsd1$eq[[ 2 ]] ) ) )

print( all.equal( mf, model.frame( fitw2slsd2e ) ) )
print( all.equal( mf1, model.frame( fitw2slsd2e$eq[[ 1 ]] ) ) )

print( all.equal( mf, model.frame( fitw2slsd3e ) ) )
print( all.equal( mf2, model.frame( fitw2slsd3e$eq[[ 2 ]] ) ) )


## **************** model matrix ************************
print( mm <- model.matrix( fitw2sls1e ) )
print( mm1 <- model.matrix( fitw2sls1e$eq[[ 1 ]] ) )
print( mm2 <- model.matrix( fitw2sls1e$eq[[ 2 ]] ) )
fitw2sls1e$eq[[ 1 ]]$modelMatrix <- NULL
fitw2sls1e$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitw2sls1e ) ) )
print( all.equal( mm1, model.matrix( fitw2sls1e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls1e$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitw2sls2e ) ) )
print( all.equal( mm1, model.matrix( fitw2sls2e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls2e$eq[[ 2 ]] ) ) )
fitw2sls2e$eq[[ 1 ]]$modelMatrix <- NULL
fitw2sls2e$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitw2sls2e ) ) )
print( all.equal( mm1, model.matrix( fitw2sls2e$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls2e$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitw2slsd3 ) ) )
print( all.equal( mm1, model.matrix( fitw2slsd3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2slsd3$eq[[ 2 ]] ) ) )
fitw2slsd3$eq[[ 1 ]]$modelMatrix <- NULL
fitw2slsd3$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitw2slsd3 ) ) )
print( all.equal( mm1, model.matrix( fitw2slsd3$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2slsd3$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitw2sls4 ) ) )
print( all.equal( mm1, model.matrix( fitw2sls4$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls4$eq[[ 2 ]] ) ) )
fitw2sls4$eq[[ 1 ]]$modelMatrix <- NULL
fitw2sls4$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitw2sls4 ) ) )
print( all.equal( mm1, model.matrix( fitw2sls4$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls4$eq[[ 2 ]] ) ) )

print( all.equal( mm, model.matrix( fitw2sls5 ) ) )
print( all.equal( mm1, model.matrix( fitw2sls5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls5$eq[[ 2 ]] ) ) )
fitw2sls5$eq[[ 1 ]]$modelMatrix <- NULL
fitw2sls5$eq[[ 2 ]]$modelMatrix <- NULL
print( all.equal( mm, model.matrix( fitw2sls5 ) ) )
print( all.equal( mm1, model.matrix( fitw2sls5$eq[[ 1 ]] ) ) )
print( all.equal( mm2, model.matrix( fitw2sls5$eq[[ 2 ]] ) ) )
