\name{predict.systemfit}
\alias{predict.systemfit}

\title{Predictions from Equation System Estimation}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction.
}

\usage{
\method{predict}{systemfit}( object, data=object$data, level=0.05 )
}

\arguments{
  \item{object}{an object of type \code{systemfit}.}
  \item{data}{data frame in which to predict.}
  \item{level}{level of significance.}
}

\value{
   \code{predict.systemfit} returns a matrix. The columns of this matrix
      contain for each equation the predicted values, the standard errors of
      the predicted values and the lower and upper confidence limits of prediction.
}

\references{

  Greene, W. H. (1993)
  \emph{Econometric Analysis, Second Edition}, Macmillan.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{systemfit}}, \code{\link{predict.systemfit.equation}}
   and \code{\link{predict}}
}

\examples{
\dontrun{library( systemfit )}
data( kmenta )
demand <- q ~ p + d
supply <- q ~ p + f + a
labels <- list( "demand", "supply" )
system <- list( demand, supply )

## OLS estimation
fitols <- systemfit("OLS", system, labels, data=kmenta )

## calculate predicted values and limits
predict( fitols )
}

\keyword{models}



