\name{coef.systemfit}
\alias{coef.systemfit}
\alias{coef.systemfit.equation}
\alias{coef.summary.systemfit}
\alias{coef.summary.systemfit.equation}
\title{Coefficients of systemfit object}

\description{
   These functions extract the coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{coef}{systemfit}( object, transformed = FALSE, \dots )

   \method{coef}{systemfit.equation}( object, \dots )

   \method{coef}{summary.systemfit}( object, transformed = FALSE, \dots )

   \method{coef}{summary.systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit},
      \code{systemfit.equation}, \code{summary.systemfit}, or
      \code{summary.systemfit.equation}.}
   \item{transformed}{logical. If \code{TRUE},
      the transformed coefficients
      (due to transforming the regressor matrix by argument \code{restrict.reg})
      rather than the re-transformed original coefficients
      are returned.}
   \item{\dots}{other arguments.}
}

\value{
   \code{coef.systemfit} returns a vector of all estimated coefficients.

   \code{coef.systemfit.equation} returns a vector of the estimated coefficients
      of a single equation.

   \code{coef.summary.systemfit} returns a matrix of all estimated coefficients,
      their standard errors, t-values, and p-values.

   \code{coef.summary.systemfit.equation} returns a matrix of the estimated
      coefficients of a single equation, their standard errors, t-values,
      and p-values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{coef}}
}


\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## all coefficients
coef( fitols )
coef( summary( fitols ) )

## coefficients of the first equation
coef( fitols$eq[[1]] )
coef( summary( fitols$eq[[1]] ) )

## coefficients of the second equation
coef( fitols$eq[[2]] )
coef( summary( fitols$eq[[2]] ) )

## transformed coefficients
coefMap <- matrix( 0, 7, 6 )
colnames( coefMap ) <- c( "demIntercept", "demPrice", "demIncome",
   "supIntercept", "supPrice2", "supTrend" )
coefMap[ 1, "demIntercept" ] <- 1
coefMap[ 2, "demPrice" ]     <- 1
coefMap[ 3, "demIncome" ]    <- 1
coefMap[ 4, "supIntercept" ] <- 1
coefMap[ 5, "supPrice2" ]    <- 1
coefMap[ 6, "supPrice2" ]    <- 1
coefMap[ 7, "supTrend" ]     <- 1
fitols3 <- systemfit( system, data = Kmenta, restrict.reg = coefMap )
coef( fitols3, transformed  = TRUE )
coef( fitols3 )
}

\keyword{models}



