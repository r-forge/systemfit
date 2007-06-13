systemfit.control <- function(
      maxiter = 1,
      tol = 1e-5,
      methodRCov = "geomean",
      centerResiduals = FALSE,
      method3sls = "GLS",
      single.eq.sigma = NULL,
      solvetol = .Machine$double.eps,
      residCovRestricted = TRUE,
      returnModelFrame = TRUE,
      returnModelMatrix = TRUE,
      returnInstMatrix = TRUE,
      returnResponse = TRUE )
{
   result <- list()

   ## maxiter
   if( maxiter <= 0 || round( maxiter ) != maxiter ) {
      stop( "control parameter 'maxiter' must be a positive integer" )
   }
   result$maxiter <- maxiter

   ## tol
   if( tol <= 0 || !is.numeric( tol ) || length( tol ) != 1 ) {
      stop( "control parameter 'tol' must be a positive scalar" )
   }
   result$tol <- tol

   ## methodRCov
   if( !( methodRCov %in% c( "noDfCor", "geomean", "max", "Theil" ) ) ) {
      stop( "control parameter 'methodRCov' must be either",
         " 'noDfCor', 'geomean', 'max', or 'Theil'" )
   }
   result$methodRCov <- methodRCov

   ## centerResiduals
   if( !is.logical( centerResiduals ) || length( centerResiduals ) != 1 ) {
      stop( "control parameter 'centerResiduals' must be logical" )
   }
   result$centerResiduals <- centerResiduals

   ## method3sls
   if( !( method3sls %in% c( "GLS", "IV", "GMM", "Schmidt", "EViews" ) ) ) {
      stop( "control parameter 'method3sls' must be either",
         " 'GLS', 'IV', 'GMM', 'Schmidt', or 'EViews'" )
   }
   result$method3sls <- method3sls

   ## single.eq.sigma
   if( ( !is.logical( single.eq.sigma ) || length( single.eq.sigma ) != 1 )
         && !is.null( single.eq.sigma ) ) {
      stop( "control parameter 'single.eq.sigma' must be logical or NULL" )
   }
   result$single.eq.sigma <- single.eq.sigma

   ## solvetol
   if( solvetol <= 0 || !is.numeric( solvetol ) || length( solvetol ) != 1 ) {
      stop( "control parameter 'solvetol' must be a positive scalar" )
   }
   result$solvetol <- solvetol

   ## residCovRestricted
   if( !is.logical( residCovRestricted ) || length( residCovRestricted ) != 1 ) {
      stop( "control parameter 'residCovRestricted' must be logical" )
   }
   result$residCovRestricted <- residCovRestricted

   ## returnModelFrame
   if( !is.logical( returnModelFrame ) || length( returnModelFrame ) != 1 ) {
      stop( "control parameter 'returnModelFrame' must be logical" )
   }
   result$returnModelFrame <- returnModelFrame

   ## returnModelMatrix
   if( !is.logical( returnModelMatrix ) || length( returnModelMatrix ) != 1 ) {
      stop( "control parameter 'returnModelMatrix' must be logical" )
   }
   result$returnModelMatrix <- returnModelMatrix

   ## returnInstMatrix
   if( !is.logical( returnInstMatrix ) || length( returnInstMatrix ) != 1 ) {
      stop( "control parameter 'returnInstMatrix' must be logical" )
   }
   result$returnInstMatrix <- returnInstMatrix

   ## returnResponse
   if( !is.logical( returnResponse ) || length( returnResponse ) != 1 ) {
      stop( "control parameter 'returnResponse' must be logical" )
   }
   result$returnResponse <- returnResponse

   return( result )
}