systemfit.control <- function(
      maxiter = 1,
      tol = 1e-5,
      methodResidCov = "geomean",
      centerResiduals = FALSE,
      residCovRestricted = TRUE,
      residCovWeighted = FALSE,
      method3sls = "GLS",
      single.eq.sigma = NULL,
      useMatrix = TRUE,
      solvetol = .Machine$double.eps,
      returnModelFrame = TRUE,
      returnModelMatrix = FALSE,
      returnInstMatrix = FALSE,
      returnResponse = FALSE )
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

   ## methodResidCov
   if( !( methodResidCov %in% c( "noDfCor", "geomean", "max", "Theil" ) ) ) {
      stop( "control parameter 'methodResidCov' must be either",
         " 'noDfCor', 'geomean', 'max', or 'Theil'" )
   }
   result$methodResidCov <- methodResidCov

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

   ## useMatrix
   if( !is.logical( useMatrix ) || length( useMatrix ) != 1 ) {
      stop( "control parameter 'useMatrix' must be logical" )
   }
   result$useMatrix <- useMatrix

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

   ## residCovWeighted
   if( !is.logical( residCovWeighted ) || length( residCovWeighted ) != 1 ) {
      stop( "control parameter 'residCovWeighted' must be logical" )
   }
   result$residCovWeighted <- residCovWeighted

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