
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

