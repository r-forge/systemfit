
% 	$Id$	

\name{summary.nlsystemfit.equation}
\alias{summary.nlsystemfit.equation}
\title{summary.nlsystemfit.equation}

\description{
  This function returns a summary of one equation of the estimated equation system.
}

\usage{
summary.nlsystemfit.equation(object,...)
}

\arguments{
  \item{object}{an object of type \code{nlsystemfit.equation}.}
  \item{...}{not used by user.}
}

\value{
  \code{summary.nlsystemfit.equation} returns an object of type \code{nlsystemfit.equation}.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com}}

\seealso{\code{\link{nlsystemfit}},
  \code{\link{summary.nlsystemfit.system}}, and 
  \code{\link{print.nlsystemfit.equation}} }


\examples{
library( systemfit )
data( ppine )

hg.formula <- hg ~ exp( h0 + h1*log(tht) + h2*tht^2 + h3*cr + h4*elev )
dg.formula <- dg ~ exp( d0 + d1*log(dbh) + d2*hg + d3*cr + d4*ba  )
hd.formula <-  tht ~ 4.5 + exp( b0 + b1 * dbh^b2 )
labels <- list( "height.growth", "diameter.growth", "tht-dbh" )
inst <- ~ dbh + ba + elev + cr 
start.values <- c(h0=-0.5, h1=0.5, h2=-0.0001, h3=0.1, h4=0.0001,
                  d0=-0.1, d1=0.2, d2=0.05, d3=0.07, d4=-0.025,
                  b0=7.0, b1=-6.0, b2=-0.3 )
model <- list( hg.formula, dg.formula, hd.formula )

model.ols <- nlsystemfit( "OLS", model, start.values, data=ppine, eqnlabels=labels )
model.3sls <- nlsystemfit( "3SLS", model, start.values, data=ppine,
                                   eqnlabels=labels, inst=inst )

#print the results
print( model.ols )
print( model.3sls )
}



\keyword{models}
\keyword{regression}
\keyword{nonlinear}
