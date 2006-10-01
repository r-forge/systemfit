
%    $Id$


\name{systemfit}
\alias{systemfit}
\alias{systemfitClassic}

\title{Linear Equation System Estimation}

\description{
  Fits a set of linear structural equations using Ordinary Least
  Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Weighted Two-Stage Least Squares (W2SLS)
  or Three-Stage Least Squares (3SLS).
}

\usage{
systemfit( eqns, method = "OLS",
           inst=NULL, data=list(), R.restr=NULL,
           q.restr=matrix(0,max(nrow(R.restr),0),1),
           TX=NULL, control = systemfit.control( ... ), ... )

systemfitClassic( formula, method = "OLS", eqnVar, timeVar, data,
                  pooled = FALSE, ... )
}

\arguments{
   \item{eqns}{a list of structural equations to be estimated; a regression
      constant is implied if not explicitly omitted.}
   \item{method}{the estimation method, one of "OLS", "WLS", "SUR", "WSUR",
      "2SLS", "W2SLS", "3SLS", or "W3SLS" (see details);
      iterated estimation methods can be specified by setting control parameter
      \code{maxiter} larger than 1 (e.g. 500).}
   \item{inst}{one-sided model formula specifying instrumental variables
      or a list of one-sided model formulas if different instruments should
      be used for the different equations (only needed for 2SLS, W2SLS and
      3SLS estimations).}
   \item{data}{an optional data frame containing the variables in the model.
      By default the variables are taken from the environment from which
      systemfit is called.}
   \item{R.restr}{an optional j x k matrix to impose linear
      restrictions on the parameters by \code{R.restr} * \eqn{b} = \code{q.restr}
      (j = number of restrictions, k = number of all parameters,
      \eqn{b} = vector of all parameters).}
   \item{q.restr}{an optional j x 1 matrix to impose linear
      restrictions (see \code{R.restr}); default is a j x 1 matrix
      that contains only zeros.}
   \item{TX}{an optional matrix to transform the regressor matrix and,
      hence, also the coefficient vector (see details).}
   \item{control}{list of control parameters.
      The default is constructed by the function \code{\link{systemfit.control}}.
      See the documentation of \code{\link{systemfit.control}} for details.}
   \item{formula}{formula to be estimated (for each equation).}
   \item{eqnVar}{variable name indicating the equation to which the
      observation belongs.}
   \item{timeVar}{variable name indicating the time.}
   \item{pooled}{logical, restrict coefficients to be equal in all equations.}
   \item{...}{arguments passed to \code{\link{systemfit.control}}
      (from  \code{\link{systemfit}})
      or passed to \code{\link{systemfit}}
      (from  \code{\link{systemfitClassic}}).}
}

\details{
   \code{systemfitClassic} is a wrapper function for \code{systemfit}
   that can be applied to panel-like data in long format
   if the regressors are the same for all equations.

   If argument \code{method} is "WSUR" or "W3SLS",
   the "SUR" or "3SLS" estimation uses a residual variance covariance matrix
   that is calculated from a "WLS" or "W2SLS" estimation, respectively
   (and not from an "OLS" or "2SLS" estimation as for a standard "SUR" or "3SLS"
   estimation).
   The "WSUR" method is the default method of command "TSCS"
   in the software LIMDEP that carries out "SUR" estimations
   in which all coefficient vectors are constrained to be equal
   (personal information from W.H. Greene, 2006/02/16).
   If no cross-equation restrictions are imposed, "WSUR" and "W3SLS"
   generate identical results compared to "SUR" and "3SLS", respectively.

   The matrix \code{TX} transforms the regressor matrix (\eqn{X}) by
   \eqn{X^{*} = X *} \code{TX}. Thus, the vector of coefficients is now
   \eqn{b =} \code{TX} \eqn{\cdot b^{*}} , where \eqn{b} is the original (stacked) vector
   of all coefficients and \eqn{b^{*}} is the new coefficient vector that is
   estimated instead. Thus, the elements of vector \eqn{b} are
   \eqn{b_i = \sum_j TX_{ij} \cdot b^{*}_j} \cr
   The \code{TX} matrix can be used to change the order of the
   coefficients and also to restrict coefficients (if \code{TX} has less
   columns than it has rows). However restricting coefficients
   by the \code{TX} matrix is less powerfull and flexible than the
   restriction by providing the \code{R.restr} matrix and the
   \code{q.restr} vector. The advantage of restricting the coefficients
   by the \code{TX} matrix is that the matrix that is inverted for
   estimation gets smaller by this procedure, while it gets larger
   if the restrictions are imposed by \code{R.restr} and \code{q.restr}.
}

\value{
  \code{systemfit} returns a list of the class \code{systemfit} and
  contains all results that belong to the whole system.
  This list contains one special object: "eq". It is a list and contains
  one object for each estimated equation. These objects are of the class
  \code{systemfit.equation} and contain the results that belong only to the
  regarding equation.

  The objects of the class \code{systemfit} and
  \code{systemfit.equation} have the following components (the elements of
  the latter are marked with an asterisk (\eqn{*})):

  \item{call}{the matched call.}
  \item{method}{estimation method.}
  \item{nObs}{total number of observations of all equations.}
  \item{nCoefAll}{total number of exogenous variables/coefficients
      in all equations (including the constants).}
  \item{nCoefEq}{number of exogenous variables/coefficients in each equation
      (including the constants).}
  \item{nCoefLiAll}{total number of linear independent coefficients
      in all equations.}
  \item{nCoefLiEq}{number of linear independent coefficients in each equation.}
  \item{df.residual}{degrees of freedom of the whole system.}
  \item{iter}{number of iteration steps.}
  \item{coefficients}{vector of all estimated coefficients.}
  \item{bt}{coefficient vector transformed by \code{TX}.}
  \item{bcov}{estimated covariance matrix of \code{coefficients}.}
  \item{btcov}{covariance matrix of \code{bt}.}
  \item{rcov}{estimated residual covariance matrix.}
  \item{rcovest}{residual covariance matrix used for estimation (only SUR and 3SLS).}
  \item{R.restr}{the restriction matrix.}
  \item{q.restr}{the restriction vector.}
  \item{TX}{matrix used to transform the regressor matrix.}
  \item{control}{list of control parameters used for the estimation.}

  ## elements of the class systemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnLabel*}{the equation label of the ith equation (from the labels list).}
  \item{terms*}{the 'terms' object used for the ith equation.}
  \item{inst*}{instruments of the ith equation (only 2SLS and 3SLS).}
  \item{nObs*}{number of observations of the ith equation.}
  \item{nCoef*}{number of exogenous variables/coefficients in the ith equation
     (including the constant).}
  \item{nCoefLi*}{number of linear independent coefficients in the ith equation
     (including the constant differs from \code{nCoef*} only if there are
     restrictions that are not cross-equation).}
  \item{nCoefAll*}{total number of exogenous variables/coefficients
      in all equation (including constants).}
  \item{nCoefLiAll*}{total number of linear independent coefficients
      in all equations (including constants).}
  \item{df.residual*}{degrees of freedom of the ith equation.}
  \item{df.residual.sys*}{degrees of freedom of the whole system.}
  \item{coefficients*}{estimated coefficients of the ith equation.}
  \item{covb*}{estimated covariance matrix of \code{coefficients}.}

  \item{response*}{if requested, the response of the ith equation.}
  \item{modelMatrix*}{if requested, the model matrix of the ith equation.}
  \item{instMatrix*}{if requested, the matrix of instrumental variables
      of the ith equation (only 2SLS and 3SLS).}
  \item{modelFrame*}{if requested (the default), the model frame of the ith equation.}
  \item{fitted.values*}{vector of fitted values of the ith equation.}
  \item{residuals*}{vector of residuals of the ith equation.}

  \item{ssr*}{sum of squared residuals of the ith equation.}
  \item{sigma*}{estimated standard error of the residuals (\eqn{\hat{\sigma}}) of the ith equation.}
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Judge, George G.; W. E. Griffiths; R. Carter Hill; Helmut Lütkepohl and Tsoung-Chao Lee (1985)
  \emph{The Theory and Practice of Econometrics, Second Edition}, Wiley.

  Kmenta, J. (1997)
  \emph{Elements of Econometrics, Second Edition}, University of
  Michigan Publishing.

  Schmidt, P. (1990)
  \emph{Three-Stage Least Squares with different Instruments for different equations},
  Journal of Econometrics 43, p. 389-394.

  Theil, H. (1971)
  \emph{Principles of Econometrics}, Wiley, New York.
}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
  Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}
}

\seealso{\code{\link{lm}} and \code{\link{nlsystemfit}}}

\examples{
data( "Kmenta" )
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list( demand = eqDemand, supply = eqSupply )

## OLS estimation
fitols <- systemfit( system, data=Kmenta )
print( fitols )

## OLS estimation with 2 restrictions
Rrestr <- matrix(0,2,7)
qrestr <- matrix(0,2,1)
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr[2,1] <-  0.5
fitols2 <- systemfit( system, data = Kmenta,
                      R.restr = Rrestr, q.restr = qrestr )
print( fitols2 )

## iterated SUR estimation
fitsur <- systemfit( system, "SUR", data = Kmenta, maxit = 100 )
print( fitsur )

## 2SLS estimation
inst <- ~ income + farmPrice + trend
fit2sls <- systemfit( system, "2SLS", inst = inst, data = Kmenta )
print( fit2sls )

## 2SLS estimation with different instruments in each equation
inst1 <- ~ income + farmPrice
inst2 <- ~ income + farmPrice + trend
instlist <- list( inst1, inst2 )
fit2sls2 <- systemfit( system, "2SLS", inst = instlist, data = Kmenta )
print( fit2sls2 )

## 3SLS estimation with GMM-3SLS formula
inst <- ~ income + farmPrice + trend
fit3sls <- systemfit( system, "3SLS", inst = inst, data = Kmenta,
   method3sls = "GMM" )
print( fit3sls )


## Examples how to use systemfitClassic()
## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
formulaGrunfeld <- invest ~ value + capital
# OLS
theilOls <- systemfitClassic( formulaGrunfeld, "OLS", "firm", "year",
   data = GrunfeldTheil )
summary( theilOls )
# SUR
theilSur <- systemfitClassic( formulaGrunfeld, "SUR", "firm", "year",
   data = GrunfeldTheil, methodRCov = "noDfCor" )
summary( theilSur )


## Further examples are in the documentation to the data sets
## 'KleinI' and 'GrunfeldGreene'.
}

\keyword{models}
\keyword{regression}



