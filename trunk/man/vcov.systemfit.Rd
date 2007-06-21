\name{vcov.systemfit}
\alias{vcov.systemfit}
\alias{vcov.systemfit.equation}
\title{Variance covariance matrix of coefficients}

\description{
   These functions extract the variance covariance matrix of the
   coefficients from an object returned by
   \code{\link{systemfit}}.
}

\usage{
   \method{vcov}{systemfit}( object, transformed = FALSE, \dots )

   \method{vcov}{systemfit.equation}( object, \dots )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{transformed}{logical. If \code{TRUE},
      the covariance matrix of the transformed coefficients
      (due to transforming the regressor matrix by argument \code{TX})
      is returned.}
   \item{\dots}{other arguments.}
}

\value{
   \code{vcov.systemfit} returns the variance covariance matrix
      of all estimated coefficients.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{
   \code{\link{systemfit}}, \code{\link{vcov}}
}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## perform OLS on each of the equations in the system
fitols <- systemfit( system, data = Kmenta )

## variance covariance matrix of all coefficients
vcov( fitols )

## variance covariance matrix of the coefficients in the first equation
vcov( fitols$eq[[1]] )

## variance covariance matrix of the coefficients in the second equation
vcov( fitols$eq[[2]] )

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
fitsur <- systemfit( system, "SUR", data = Kmenta, TX = coefMap )
vcov( fitsur, transformed  = TRUE )
vcov( fitsur )
}

\keyword{models}



