
% 	$Id$	

\keyword{datasets}
\name{ritchie}
\alias{ritchie}
\non_function{}
\title{Young Tree Growth Data for Southwestern Oregon and Northern California}
\usage{data(ritchie)}
\description{
This data is a subset from a much larger dataset of young tree growth
observations (Ritchie, put the correct year in here...).

The \code{ritchie} data frame has 428 rows and 7 columns.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{bd}{
      basal diameter (the diameter of the tree at 6 inches above the ground.
    }
    \item{ht}{
      total stem height for the tree.
    }
    \item{hg}{
      two-year height increment.
    }
    \item{dg}{
      two-year basal diameter increment.
    }
    \item{whc}{
      water holding capacity of the soil on the plot.
    }
    \item{cr}{
      the ratio of the stem in live crown.
    }
    \item{bainlh}{
      hardwood basal area in larger (fatter) stems than the subject tree.
    }
  }
}

\details{
 The exogenous variables are \code{bd}, \code{ht}, \code{whc},
 \code{bainlh}; the endogenous variables \code{dg} and \code{hg}. There
 are no lagged variables in the dataset and the observations are for a
 single remeasurement. 
}
\source{
 Ritchie, M. W. (19XX)
 \emph{Get the pub these things come from, Second Edition}, Macmillan.
}
\examples{
data(ritchie)
}

\keyword{datasets}


