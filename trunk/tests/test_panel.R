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

# SUR
theilSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldTheil, methodRCov = "noDfCor" )
print( theilSur )
summary( theilSur )


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

# OLS Pooled
greeneOlsPooled <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldGreene, pooled = TRUE )
print( greeneOlsPooled )
summary( greeneOlsPooled )
sum( sapply( greeneOlsPooled$eq, function(x){return(x$ssr)}) )/97 # sigma^2

# SUR
greeneSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldGreene, methodRCov = "noDfCor" )
print( greeneSur )
summary( greeneSur )

# SUR Pooled
greeneSurPooled <- systemfit( formulaGrunfeld, "WSUR",
   data = GrunfeldGreene, pooled = TRUE, methodRCov = "noDfCor" )
print( greeneSurPooled )
summary( greeneSurPooled )
