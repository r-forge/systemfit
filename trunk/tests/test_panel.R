library( systemfit )
library( plm )

## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
pdata.frame( GrunfeldTheil, "firm", "year" )
formulaGrunfeld <- invest ~ value + capital

# OLS
theilOls <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldTheil )
print( theilOls )
summary( theilOls )
coef( theilOls )
coef( summary(theilOls ) )
vcov( theilOls )
residuals( theilOls )
confint( theilOls )
fitted(theilOls  )
logLik( theilOls )
model.frame( theilOls )
model.matrix( theilOls )

# SUR
theilSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldTheil, methodRCov = "noDfCor" )
print( theilSur )
summary( theilSur )
coef( theilSur )
coef( summary( theilSur ) )
vcov( theilSur )
residuals( theilSur )
confint( theilSur )
fitted( theilSur )
logLik( theilSur )
model.frame( theilSur )
model.matrix( theilSur )


## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
data( "GrunfeldGreene" )
pdata.frame( GrunfeldGreene, "firm", "year" )
formulaGrunfeld <- invest ~ value + capital

# OLS
greeneOls <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldGreene )
print( greeneOls )
summary( greeneOls )
sapply( greeneOls$eq, function(x){return(x$ssr/20)} ) # sigma^2
coef( greeneOls )
coef( summary( greeneOls ) )
vcov( greeneOls )
residuals( greeneOls )
confint(greeneOls  )
fitted( greeneOls )
logLik( greeneOls )
model.frame( greeneOls )
model.matrix( greeneOls )

# OLS Pooled
greeneOlsPooled <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldGreene, pooled = TRUE )
print( greeneOlsPooled )
summary( greeneOlsPooled )
sum( sapply( greeneOlsPooled$eq, function(x){return(x$ssr)}) )/97 # sigma^2
coef( greeneOlsPooled )
coef( greeneOlsPooled, modified.reg = TRUE )
coef( summary( greeneOlsPooled ) )
coef( summary( greeneOlsPooled ), modified.reg = TRUE )
vcov( greeneOlsPooled )
vcov( greeneOlsPooled, modified.reg = TRUE )
residuals( greeneOlsPooled )
confint( greeneOlsPooled )
fitted( greeneOlsPooled )
logLik( greeneOlsPooled )
model.frame( greeneOlsPooled )
model.matrix( greeneOlsPooled )

# SUR
greeneSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldGreene, methodRCov = "noDfCor" )
print( greeneSur )
summary( greeneSur )
coef( greeneSur )
coef( summary( greeneSur ) )
vcov( greeneSur )
residuals( greeneSur )
confint( greeneSur )
fitted( greeneSur )
logLik( greeneSur )
model.frame( greeneSur )
model.matrix( greeneSur )

# SUR Pooled
greeneSurPooled <- systemfit( formulaGrunfeld, "WSUR",
   data = GrunfeldGreene, pooled = TRUE, methodRCov = "noDfCor" )
print( greeneSurPooled )
summary( greeneSurPooled )
coef( greeneSurPooled )
coef( greeneSurPooled, modified.reg = TRUE )
coef( summary( greeneSurPooled ) )
coef( summary( greeneSurPooled ), modified.reg = TRUE )
vcov( greeneSurPooled )
vcov( greeneSurPooled, modified.reg = TRUE )
residuals( greeneSurPooled )
confint( greeneSurPooled )
fitted( greeneSurPooled )
logLik( greeneSurPooled )
model.frame( greeneSurPooled )
model.matrix( greeneSurPooled )
