\name{linear.hypothesis.systemfit}
\alias{linear.hypothesis.systemfit}

\title{Test Linear Hypothesis}

\description{
   Testing linear hypothesis on parameters of a system of equations
   by an F-test or Wald-test.
}

\usage{
   \method{linear.hypothesis}{systemfit}( model,
      hypothesis.matrix, rhs = NULL, test = c( "F", "Chisq" ), ... )
}

\arguments{
   \item{model}{a fitted object of type \code{systemfit}.}
   \item{hypothesis.matrix}{matrix (or vector) giving linear combinations
      of coefficients by rows,
      or a character vector giving the hypothesis in symbolic form
      (see documentation of \code{\link[car]{linear.hypothesis}} for details).}
   \item{rhs}{optional right-hand-side vector for hypothesis,
      with as many entries as rows in the hypothesis matrix;
      if omitted, it defaults to a vector of zeroes.}
   \item{test}{character string, "\code{F}" or "\code{Chisq}",
      specifying whether to compute the
      finite-sample F statistic (with approximate F distribution) or
      the large-sample Wald statistic
      (with asymptotic Chi-squared distribution).}
   \item{\dots}{further arguments passed to
      \code{\link[car]{linear.hypothesis.default}}.}
}

\details{
   The Wald-statistic for sytems of equations is
   \deqn{ W =
      ( R \hat{b} - q )'
      ( R \widehat{Cov} [ \hat{b} ] R' )^{-1}
      ( R \hat{b} - q )
   }
   Asymptotically, \eqn{W} has a \eqn{\chi^2}
   distribution with \eqn{j} degrees of freedom
   under the null hypothesis
   (Greene, 2003, p. 347).
}

\value{
   see documentation of \code{\link[car]{linear.hypothesis}}
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{systemfit}}, \code{\link[car]{linear.hypothesis}}
   \code{\link{ftest.systemfit}},
   \code{\link{lrtest.systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## unconstrained SUR estimation
fitsur <- systemfit( system, method = "SUR", data=Kmenta )

# create restriction matrix to test whether \eqn{beta_2 = \beta_6}
R1 <- matrix( 0, nrow = 1, ncol = 7 )
R1[ 1, 2 ] <- 1
R1[ 1, 6 ] <- -1

## perform Wald-test
waldTest1 <- linear.hypothesis( fitsur, R1, test = "Chisq" )
print( waldTest1 )   # rejected

# create restriction matrix to test whether \eqn{beta_2 = - \beta_6}
R2 <- matrix( 0, nrow = 1, ncol = 7 )
R2[ 1, 2 ] <- 1
R2[ 1, 6 ] <- 1

## perform Wald-test
waldTest2 <- linear.hypothesis( fitsur, R2, test = "Chisq" )
print( waldTest2 )   # accepted
}

\keyword{models}



