
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


## ********************* W2SLS *****************
fitw2sls1 <- systemfit( system, "W2SLS", data = Kmenta, inst = inst )
print( summary( fitw2sls1 ) )
print( round( fitw2sls1$bcov, digits = 6 ) )

## ********************* W2SLS (EViews-like) *****************
fitw2sls1e <- systemfit( system, "W2SLS", data = Kmenta, inst = inst,
   methodRCov = "noDfCor" )
print( summary( fitw2sls1e, probDfSys = TRUE ) )
print( round( fitw2sls1e$bcov, digits = 6 ) )

## ********************* W2SLS with restriction *******************
fitw2sls2 <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restrm,
   inst = inst )
print( summary( fitw2sls2 ) )
print( round( fitw2sls2$bcov, digits = 6 ) )

## ********************* W2SLS with restriction (EViews-like) **************
fitw2sls2e <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restrm,
   inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls2e, probDfSys = TRUE ) )
print( round( fitw2sls2e$bcov, digits = 6 ) )

## ********************* W2SLS with restriction via TX *******************
fitw2sls3 <- systemfit( system, "W2SLS", data = Kmenta, TX = tc, inst = inst )
print( summary( fitw2sls3 ) )
print( round( fitw2sls3$bcov, digits = 6 ) )

## ********************* W2SLS with restriction via TX (EViews-like) **************
fitw2sls3e <- systemfit( system, "W2SLS", data = Kmenta, TX = tc, inst = inst,
   methodRCov = "noDfCor" )
print( summary( fitw2sls3e, probDfSys = TRUE ) )
print( round( fitw2sls3e$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions ********************
fitw2sls4 <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst )
print( summary( fitw2sls4 ) )
print( round( fitw2sls4$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions (EViews-like) **************
fitw2sls4e <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restr2m,
   q.restr = restr2q, inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls4e, probDfSys = TRUE ) )
print( round( fitw2sls4e$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions via R and TX ******************
fitw2sls5 <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst )
print( summary( fitw2sls5 ) )
print( round( fitw2sls5$bcov, digits = 6 ) )

## ***************** W2SLS with 2 restrictions via R and TX (EViews-like) **************
fitw2sls5e <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restr3m,
   q.restr = restr3q, TX = tc, inst = inst, methodRCov = "noDfCor" )
print( summary( fitw2sls5e, probDfSys = TRUE ) )
print( round( fitw2sls5e$bcov, digits = 6 ) )

## ****** 2SLS estimation with different instruments **********************
fitw2slsd1 <- systemfit( system, "W2SLS", data = Kmenta, inst = instlist )
print( summary( fitw2slsd1 ) )
print( round( fitw2slsd1$bcov, digits = 6 ) )

## ****** 2SLS estimation with different instruments (EViews-like)******************
fitw2slsd1e <- systemfit( system, "W2SLS", data = Kmenta, inst = instlist,
   methodRCov = "noDfCor" )
print( summary( fitw2slsd1e, probDfSys = TRUE ) )
print( round( fitw2slsd1e$bcov, digits = 6 ) )

## **** W2SLS estimation with different instruments and restriction ********
fitw2slsd2 <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist )
print( summary( fitw2slsd2 ) )
print( round( fitw2slsd2$bcov, digits = 6 ) )

## **** W2SLS estimation with different instruments and restriction (EViews-like)*
fitw2slsd2e <- systemfit( system, "W2SLS", data = Kmenta, R.restr = restrm,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fitw2slsd2e, probDfSys = TRUE ) )
print( round( fitw2slsd2e$bcov, digits = 6 ) )

## ** W2SLS estimation with different instruments and restriction via TX ****
fitw2slsd3 <- systemfit( system, "W2SLS", data = Kmenta, TX = tc,
   inst = instlist)
print( summary( fitw2slsd3 ) )
print( round( fitw2slsd3$bcov, digits = 6 ) )

## W2SLS estimation with different instruments and restriction via TX (EViews-like)
fitw2slsd3e <- systemfit( system, "W2SLS", data = Kmenta, TX = tc,
   inst = instlist, methodRCov = "noDfCor" )
print( summary( fitw2slsd3e, probDfSys = TRUE ) )
print( round( fitw2slsd3e$bcov, digits = 6 ) )


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


## *********** confidence intervals of coefficients *************
print( confint( fitw2sls1e, probDfSys = TRUE ) )
print( confint( fitw2sls1e$eq[[ 1 ]], level = 0.9, probDfSys = TRUE ) )

print( confint( fitw2sls2, level = 0.9 ) )
print( confint( fitw2sls2$eq[[ 2 ]], level = 0.99 ) )

print( confint( fitw2sls3, level = 0.99 ) )
print( confint( fitw2sls3$eq[[ 1 ]], level = 0.5 ) )

print( confint( fitw2sls4e, level = 0.5, probDfSys = TRUE ) )
print( confint( fitw2sls4e$eq[[ 2 ]], level = 0.25, probDfSys = TRUE ) )

print( confint( fitw2sls5, level = 0.25 ) )
print( confint( fitw2sls5$eq[[ 1 ]], level = 0.975 ) )

print( confint( fitw2slsd1, level = 0.975 ) )
print( confint( fitw2slsd1$eq[[ 2 ]], level = 0.999 ) )

print( confint( fitw2slsd2e, level = 0.999, probDfSys = TRUE ) )
print( confint( fitw2slsd2e$eq[[ 1 ]], level = 0.01, probDfSys = TRUE ) )

print( confint( fitw2slsd3e, level = 0.01, probDfSys = TRUE ) )
print( confint( fitw2slsd3e$eq[[ 2 ]], probDfSys = TRUE ) )


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
predictData$price <- Kmenta$price * 0.9
predictData$income <- Kmenta$income * 1.1

print( predict( fitw2sls1e, se.fit = TRUE, interval = "prediction",
   probDfSys = TRUE ) )
print( predict( fitw2sls1e$eq[[ 1 ]] ) )

print( predict( fitw2sls2, se.pred = TRUE, interval = "confidence",
   level = 0.999, data = predictData ) )
print( predict( fitw2sls2$eq[[ 2 ]] ) )

print( predict( fitw2sls3, se.pred = TRUE, interval = "prediction",
   level = 0.975 ) )
print( predict( fitw2sls3$eq[[ 1 ]] ) )

print( predict( fitw2sls4e, se.fit = TRUE, interval = "confidence",
   level = 0.25, probDfSys = TRUE ) )
print( predict( fitw2sls4e$eq[[ 2 ]] ) )

print( predict( fitw2sls5, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.5, data = predictData ) )
print( predict( fitw2sls5$eq[[ 1 ]] ) )

print( predict( fitw2slsd1, se.fit = TRUE, se.pred = TRUE,
   interval = "confidence", level = 0.99 ) )
print( predict( fitw2slsd1$eq[[ 2 ]] ) )

print( predict( fitw2slsd2e, se.fit = TRUE, interval = "prediction",
   level = 0.9, data = predictData, probDfSys = TRUE ) )
print( predict( fitw2slsd2e$eq[[ 1 ]] ) )

print( predict( fitw2slsd3e, se.fit = TRUE, se.pred = TRUE,
   interval = "prediction", level = 0.01, probDfSys = TRUE ) )
print( predict( fitw2slsd3e$eq[[ 2 ]] ) )


## ************** F tests ****************
# testing first restriction
print( ftest.systemfit( fitw2sls1, restrm ) )
print( ftest.systemfit( fitw2slsd1e, restrm ) )

# testing second restriction
restrOnly2m <- matrix(0,1,7)
restrOnly2q <- 0.5
restrOnly2m[1,2] <- -1
restrOnly2m[1,5] <-  1
# first restriction not imposed 
print( ftest.systemfit( fitw2sls1e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( ftest.systemfit( fitw2sls2, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2sls3, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd2e, restrOnly2m, restrOnly2q ) )
print( ftest.systemfit( fitw2slsd3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( ftest.systemfit( fitw2sls1e, restr2m, restr2q ) )
print( ftest.systemfit( fitw2slsd1, restr2m, restr2q ) )


## ************** Wald tests ****************
# testing first restriction
print( waldtest.systemfit( fitw2sls1, restrm ) )
print( waldtest.systemfit( fitw2slsd1e, restrm ) )

# testing second restriction
# first restriction not imposed
print( waldtest.systemfit( fitw2sls1e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd1, restrOnly2m, restrOnly2q ) )
# first restriction imposed
print( waldtest.systemfit( fitw2sls2, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2sls3, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd2e, restrOnly2m, restrOnly2q ) )
print( waldtest.systemfit( fitw2slsd3e, restrOnly2m, restrOnly2q ) )

# testing both of the restrictions
print( waldtest.systemfit( fitw2sls1e, restr2m, restr2q ) )
print( waldtest.systemfit( fitw2slsd1, restr2m, restr2q ) )
