
% 	$Id$	

\name{systemfit-internal}
\alias{knls}

\title{Internal systemfit functions}
\description{
  Internal systemfit functions
}
\usage{
knls( theta, eqns, data, fitmethod="OLS", parmnames, instr=NULL, S=NULL )
}

\details{

  These functions are no the called by the user.
  
  The knls is the driver function that is passed to the nlm call in
  \link{\code{systemfit}}. 
}

\keyword{internal}
