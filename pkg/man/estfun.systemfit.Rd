\name{estfun.systemfit}
\alias{estfun.systemfit}

\title{Extract Gradients of the Objective Function at each Observation}

\description{
   Extract the gradients of the objective function 
   with respect to the coefficients 
   evaluated at each observation 
   (\sQuote{Empirical Estimating Function}, 
   see \code{\link[sandwich]{estfun}}).
}

\usage{
\method{estfun}{systemfit}( obj, ... )
}

\arguments{
   \item{obj}{an object of class \code{systemfit}.}
   \item{\dots}{further arguments (currently ignored).}
}

\value{
   Matrix of gradients of the objective function
   with respect to the coefficients
   evaluated at each observation.
}

\section{Warnings}{
   The \pkg{sandwich} package must be loaded before this method
   can be used.
   
   In some cases with unbalanced systems estimated by 3SLS 
   (and perhaps also by W2SLS),
   the columns of the matrix returned by this \code{estfun} method
   do not sum up to zero.
   This indicates a bug.
   Unfortunately, the reason for this has not yet been found.
}

\author{
   Arne Henningsen
}

\seealso{\code{\link[sandwich]{estfun}}, \code{\link{systemfit}}.}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )
inst <- ~ income + farmPrice + trend

## OLS estimation
fitols <- systemfit( system, "OLS", data = Kmenta )

## obtain the estimation function
library( "sandwich" )
estfun( fitols )

## this is only true for OLS models
all.equal( estfun( fitols ),
   unlist( residuals( fitols ) ) * model.matrix( fitols ) )

# each column should sum up to (approximately) zero
colSums( estfun( fitols ) )


## 2SLS estimation
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )

## obtain the estimation function
estfun( fit2sls )

## this is only true for 2SLS models
all.equal( estfun( fit2sls ),
   unlist( residuals( fit2sls ) ) * model.matrix( fit2sls, which = "xHat" ) )

# each column should sum up to (approximately) zero
colSums( estfun( fit2sls ) )


## iterated SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta, maxit = 100 )

## obtain the estimation function
estfun( fitsur )

## this should be true for SUR and WLS models
all.equal( estfun( fitsur ),
   unlist( residuals( fitsur ) ) * 
   ( ( solve( fitsur$residCovEst ) \%x\% diag( nrow( Kmenta ) ) ) \%*\% 
      model.matrix( fitsur ) ), check.attributes = FALSE )

# each column should sum up to (approximately) zero
colSums( estfun( fitsur ) )


## 3SLS estimation
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta )

## obtain the estimation function
estfun( fit3sls )

## this should be true for 3SLS and W2SLS models
all.equal( estfun( fit3sls ),
   unlist( residuals( fit3sls ) ) * 
   ( ( solve( fit3sls$residCovEst ) \%x\% diag( nrow( Kmenta ) ) ) \%*\% 
      model.matrix( fit3sls, which = "xHat" ) ), check.attributes = FALSE )

# each column should sum up to (approximately) zero
colSums( estfun( fit3sls ) )
}

\keyword{methods}
