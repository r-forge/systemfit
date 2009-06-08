library( "systemfit" )

data( "KleinI" )
eqConsump  <- consump ~ corpProf + corpProfLag + wages
eqInvest   <- invest ~ corpProf + corpProfLag + capitalLag
eqPrivWage <- privWage ~ gnp + gnpLag + trend
inst <- ~ govExp + taxes + govWage + trend + capitalLag + corpProfLag + gnpLag
system <- list( Consumption = eqConsump, Investment = eqInvest,
   PrivateWages = eqPrivWage )

# OLS
kleinOls <- systemfit( system, data = KleinI )
summary( kleinOls )
residuals( kleinOls )
fitted( kleinOls )
predict( kleinOls, se.fit = TRUE, interval = "prediction",
   useDfSys = TRUE )
model.frame( kleinOls )
model.matrix( kleinOls )

# 2SLS
klein2sls <- systemfit( system, "2SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
summary( klein2sls )
residuals( klein2sls )
fitted( klein2sls )
predict( klein2sls, se.fit = TRUE, interval = "confidence",
   useDfSys = FALSE )
model.frame( klein2sls )
all.equal( model.matrix( kleinOls ), model.matrix( klein2sls ) )

# SUR
kleinSur <- systemfit( system, "SUR", data = KleinI,
   methodResidCov = "noDfCor" )
summary( kleinSur )
residuals( kleinSur )
fitted( kleinSur )
predict( kleinSur, se.fit = TRUE, interval = "confidence",
   useDfSys = TRUE )
all.equal( model.frame( kleinOls ), model.frame( kleinSur ) )
all.equal( model.matrix( kleinOls ), model.matrix( kleinSur ) )

# 3SLS
klein3sls <- systemfit( system, "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor" )
summary( klein3sls )
residuals( klein3sls )
fitted( klein3sls )
predict( klein3sls, se.fit = TRUE, interval = "prediction",
   useDfSys = FALSE )
all.equal( model.frame( klein2sls ), model.frame( klein3sls ) )
all.equal( model.matrix( kleinOls ), model.matrix( klein3sls ) )

# I3SLS
kleinI3sls <- systemfit( system, "3SLS", inst = inst, data = KleinI,
   methodResidCov = "noDfCor", maxit = 500 )
summary( kleinI3sls )
residuals( kleinI3sls )
fitted( kleinI3sls )
predict( kleinI3sls, se.fit = TRUE, interval = "confidence",
   useDfSys = TRUE )
all.equal( model.frame( klein2sls ), model.frame( kleinI3sls ) )
all.equal( model.matrix( kleinOls ), model.matrix( kleinI3sls ) )
