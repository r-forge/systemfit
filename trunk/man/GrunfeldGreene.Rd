\name{GrunfeldGreene}
\alias{GrunfeldGreene}
\title{Grunfeld Data as published by Greene (2003)}
\description{
   Panel data on 5 US firms for the years 1935-1954.
}

\usage{
   data("GrunfeldGreene")
}
\format{
A data frame containing 20 annual observations on 3 variables for 5 firms.
\describe{
   \item{invest}{gross investment.}
   \item{value}{market value of the firm (at the end of the previous year).}
   \item{capital}{capital stock of the firm (at the end of the previous year).}
   \item{firm}{name of the firm ("General Motors", "Chrysler",
      "General Electric", "Westinghouse" or "US Steel").}
   \item{year}{year.}
}
}

\source{
   \url{http://www.stern.nyu.edu/~wgreene/Text/tables/TableF13-1.txt}
}

\references{
Greene, W.H. (2003). \emph{Econometric Analysis}, 5th edition. Prentice Hall,
Upper Saddle River (NJ).

Grunfeld, Y. (1958). \emph{The Determinants of Corporate Investment},
Unpublished Ph.D. Dissertation, University of Chicago.
}

\keyword{datasets}
