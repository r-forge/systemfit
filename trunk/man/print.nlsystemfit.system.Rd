
% 	$Id$	

\name{print.nlsystemfit.system}
\alias{print.nlsystemfit.system}
\title{print.nlsystemfit.system}

\description{
  This function prints a summary of the estimated equation system.
}

\usage{
print.nlsystemfit.system(x, digits=6,...)
}

\arguments{
  \item{x}{an object of type \code{nlsystemfit.system}.}
  \item{digits}{number of digits to print.}
  \item{...}{not used by user.}
}

\value{
  \code{print.nlsystemfit.system} returns nothing.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{nlsystemfit}},
  \code{\link{print.nlsystemfit.equation}},
  \code{\link{summary.nlsystemfit.system}} and
  }


\examples{
library( systemfit )
data( ritchie )

hg.formula <- hg ~ exp( h0 + h1 * log( ht ) + h2 * ht^2 )
dg.formula <- dg ~ exp( d0 + d1 * log( bd ) + d2 * log( hg ) )
labels <- list( "height.growth", "diameter.growth" )
inst <- ~ ht + bd
start.values <- c(h0=-1.1, h1=0.55, h2=-0.0008,
                  d0=-5.51072, d1=0.832112, d2=0.165224 )
model <- list( hg.formula, dg.formula )

## perform Three-Stage Least Squares fit for the system
model.3sls <- nlsystemfit( "3SLS", model, start.values, data=ritchie,
eqnlabels=labels, inst=inst, pl=2, solvtol=1e-16 )

# print the results
print( model.3sls )
}


\keyword{models}
\keyword{regression}
\keyword{nonlinear}

