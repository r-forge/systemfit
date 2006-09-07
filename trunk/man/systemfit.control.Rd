\name{systemfit.control}
\alias{systemfit.control}

\title{Create list of control parameters for systemfit}

\description{
   Create a list of control pararameters for function \code{\link{systemfit}}.
   All control parameters that are not passed to this function
   are set to default values.
}

\usage{
systemfit.control(
      maxiter = 1,
      tol = 1e-5,
      methodRCov = "geomean",
      centerResiduals = FALSE,
      method3sls = "GLS",
      single.eq.sigma = NULL,
      solvetol = .Machine$double.eps,
      saveMemory = NULL )
}

\arguments{
   \item{maxiter}{maximum number of iterations for WLS, SUR, W2SLS and
      3SLS estimations.}
   \item{tol}{tolerance level indicating when to stop the iteration (only
      WLS, SUR, W2SLS and 3SLS estimations).}
   \item{methodRCov}{method for calculating the estimated residual covariance
      matrix, one of "noDfCor", "geomean", "max", or "Theil" (see details).}
   \item{centerResiduals}{logical. Subtract the means from the residuals
      of each equation before calculating the estimated residual covariance matrix.}
   \item{method3sls}{method for calculating the 3SLS estimator,
      one of "GLS", "IV", "GMM", "Schmidt", or "EViews" (see details).}
   \item{single.eq.sigma}{logical. use different \eqn{\sigma^2}s for each
      single equation to calculate the covariance matrix and the
      standard errors of the coefficients (only OLS and 2SLS)?
      If \code{single.eq.sigma} is \code{NULL}, it is automatically determined:
      It is set to \code{TRUE}, if there are any parameter restrictions imposed,
      and it is set to \code{FALSE} otherwise.}
   \item{solvetol}{tolerance level for detecting linear dependencies
      when inverting a matrix or calculating a determinant (see
      \code{\link{solve}} and \code{\link{det}}).}
   \item{saveMemory}{logical. Save memory by omitting some calculation that
      are not crucial for the basic estimation (e.g McElroy's
      \eqn{R^2})?
      If \code{saveMemory} is \code{NULL}, it is automatically determined:
      It is set to \code{TRUE}, if there are more than 1000 observation
      for all equations together,
      and it is set to \code{FALSE} otherwise.}
}

\details{
   If the estimation is iterated
      (WLS, SUR, W2SLS or 3SLS estimation with \code{maxiter}>1),
      the convergence criterion is
   \deqn{\sqrt{ \frac{ \sum_i (b_{i,g} - b_{i,g-1})^2 }{ \sum_i b_{i,g-1}^2 }}
      < \code{tol}}
   (\eqn{b_{i,g}} is the ith coefficient of the gth iteration step).

   The method for calculating the estimated covariance matrix of the residuals
   (\eqn{\hat{\Sigma}}) can be one of the following
   (see Judge et al., 1985, p. 469): \cr
   if methodRCov='noDfCor':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T}}
   if methodRCov='geomean':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {\sqrt{(T - k_i)*(T - k_j)}}}
   if methodRCov='Theil':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}{T - k_i - k_j
      + tr[X_i(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j']}}
   if methodRCov='max':
      \deqn{\hat{\sigma}_{ij} = \frac{\hat{e}_i' \hat{e}_j}
         {T - \max( k_i, k_j)}}
   If \eqn{ i = j}, the formulas 'geomean', 'Theil', and 'max' are equal.
   All these three formulas yield unbiased estimators
   for the diagonal elements of the residual covariance matrix.
   If \eqn{i \neq j}, only formula 'Theil' yields an unbiased estimator for the residual
   covariance matrix, but it is not neccessarily positive semidefinit. Thus, it is
   doubtful whether formula 'Theil' is really superior to formula 'noDfCor'
   (Theil, 1971, p. 322).

   The methods for calculating the 3SLS estimator lead to identical results
   if the same instruments are used in all equations. If different instruments
   are used in the different equations, only the GMM-3SLS estimator ("GMM")
   and the 3SLS estimator proposed by Schmidt (1990) ("Schmidt") are consistent,
   whereas "GMM" is efficient relative to "Schmidt" (see Schmidt, 1990).
}

\value{
   A list of the above components.
}

\references{
  Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Lütkepohl and Tsoung-Chao Lee (1985)
  \emph{The Theory and Practice of Econometrics, Second Edition}, Wiley.

  Schmidt, P. (1990)
  \emph{Three-Stage Least Squares with different Instruments for different equations},
  Journal of Econometrics 43, p. 389-394.

  Theil, H. (1971)
  \emph{Principles of Econometrics}, Wiley, New York.
}

\author{
  Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}
}

\seealso{\code{\link{systemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
eqSystem <- list( demand = eqDemand, supply = eqSupply )

## SUR estimation: calculation of residual covariance
## matrix without correction for degrees of freedom
fitsur <- systemfit( eqSystem, "SUR", data = Kmenta,
    control = systemfit.control( methodRCov = "noDfCor" ) )
print( fitsur )
}

\keyword{models}
\keyword{regression}



