\name{summary.systemfit}
\alias{summary.systemfit}
\title{Print output of systemfit estimation}

\description{
  This function returns a summary of the system of equations.
}

\usage{
   \method{summary}{systemfit}( object, ... )
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.systemfit} returns an object of type \code{systemfit}.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit}} and
  \code{\link{print.systemfit}}}


\examples{
data( "kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = kmenta )

## print the results
summary.systemfit( fitols )
}

\keyword{models}



