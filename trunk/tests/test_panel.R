library( systemfit )

## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
formulaGrunfeld <- invest ~ value + capital

# OLS
theilOls <- systemfitClassic( formulaGrunfeld, "OLS", "firm", "year",
   data = GrunfeldTheil )
print( theilOls )
summary( theilOls )

# SUR
theilSur <- systemfitClassic( formulaGrunfeld, "SUR", "firm", "year",
   data = GrunfeldTheil, methodRCov = "noDfCor" )
print( theilSur )
summary( theilSur )


## Repeating the OLS and SUR estimations in Greene (2003, pp. 351)
data( "GrunfeldGreene" )
formulaGrunfeld <- invest ~ value + capital

# OLS
greeneOls <- systemfitClassic( formulaGrunfeld, "OLS", "firm", "year",
   data = GrunfeldGreene )
print( greeneOls )
summary( greeneOls )
sapply( greeneOls$eq, function(x){return(x$ssr/20)} ) # sigma^2

# OLS Pooled
greeneOlsPooled <- systemfitClassic( formulaGrunfeld, "OLS", "firm", "year",
   data = GrunfeldGreene, pooled = TRUE )
print( greeneOlsPooled )
summary( greeneOlsPooled )
sum( sapply( greeneOlsPooled$eq, function(x){return(x$ssr)}) )/97 # sigma^2

# SUR
greeneSur <- systemfitClassic( formulaGrunfeld, "SUR", "firm", "year",
   data = GrunfeldGreene, methodRCov = "noDfCor" )
print( greeneSur )
summary( greeneSur )

# SUR Pooled
greeneSurPooled <- systemfitClassic( formulaGrunfeld, "WSUR", "firm", "year",
   data = GrunfeldGreene, pooled = TRUE, methodRCov = "noDfCor" )
print( greeneSurPooled )
summary( greeneSurPooled )
