\name{summary.systemfit}
\alias{summary.systemfit}
\alias{summary.systemfit.equation}
\title{Summary of systemfit estimation}

\description{
   These functions print a summary of the estimated equation system.
}

\usage{
   \method{summary}{systemfit}( object, ... )

   \method{summary}{systemfit.equation}( object, ... )
}

\arguments{
   \item{object}{an object of class \code{systemfit} or
      \code{systemfit.equation}.}
   \item{...}{not used by user.}
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de} }

\seealso{\code{\link{systemfit}}, \code{\link{print.systemfit}}}


\examples{
data( "kmenta" )
demand <- consump ~ price + income
supply <- consump ~ price + farmPrice + trend
inst <- ~ income + farmPrice + trend
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## perform OLS on each of the equations in the system
fitols <- systemfit( "OLS", system, labels, data = kmenta )

## results of the whole system
summary( fitols )

## results of the first equation
summary( fitols$eq[[1]] )

## results of the second equation
summary( fitols$eq[[2]] )
}

\keyword{models}



