
%    $Id$


\name{systemfit}
\alias{systemfit}

\title{Linear Equation System Estimation}

\description{
  Fits a set of linear structural equations using Ordinary Least
  Squares (OLS), Weighted Least Squares (WLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Weighted Two-Stage Least Squares (W2SLS)
  or Three-Stage Least Squares (3SLS).
}

\usage{
systemfit( eqns, method = "OLS",
           inst=NULL, data=list(),
           restrictions = NULL, restrict.rhs = NULL,
           TX=NULL, pooled = FALSE, control = systemfit.control( ... ), ... )
}

\arguments{
   \item{eqns}{a list of structural equations
      (objects of class \code{formula}) to be estimated or
      (if argument \code{data} is of class \code{pdata.frame})
      a single equation (object of class \code{formula})
      that is estimated for all individuals.}
   \item{method}{the estimation method, one of "OLS", "WLS", "SUR", "WSUR",
      "2SLS", "W2SLS", "3SLS", or "W3SLS" (see details);
      iterated estimation methods can be specified by setting control parameter
      \code{maxiter} larger than 1 (e.g. 500).}
   \item{inst}{one-sided model formula specifying instrumental variables
      or a list of one-sided model formulas if different instruments should
      be used for the different equations (only needed for 2SLS, W2SLS and
      3SLS estimations).}
   \item{data}{an optional data frame (or pdata.frame for panel-like data)
      containing the variables in the model.
      By default the variables are taken from the environment from which
      systemfit is called.}
   \item{restrictions}{an optional j x k matrix to impose linear
      restrictions on the parameters by \code{restrictions} * \eqn{b} = \code{restrict.rhs}
      (j = number of restrictions, k = number of all parameters,
      \eqn{b} = vector of all parameters).}
   \item{restrict.rhs}{an optional vector with j elements to impose linear
      restrictions (see \code{restrictions}); default is a vector
      that contains j zeros.}
   \item{TX}{an optional matrix to transform the regressor matrix and,
      hence, also the coefficient vector (see details).}
   \item{control}{list of control parameters.
      The default is constructed by the function \code{\link{systemfit.control}}.
      See the documentation of \code{\link{systemfit.control}} for details.}
   \item{pooled}{logical, restrict coefficients to be equal in all equations
      (only for panel-like data).}
   \item{...}{arguments passed to \code{\link{systemfit.control}}.}
}

\details{
   If argument \code{data} is of class \code{pdata.frame}
   (i.e. contains panel-like data in long format),
   argument \code{eqns} must be a single equation
   that is applied to all individuals.
   In this case, argument \code{pooled} specifies
   whether the coefficients are restricted to be equal for all
   individuals.

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
   restriction by providing the \code{restrictions} matrix and the
   \code{restrict.rhs} vector. The advantage of restricting the coefficients
   by the \code{TX} matrix is that the matrix that is inverted for
   estimation gets smaller by this procedure, while it gets larger
   if the restrictions are imposed by \code{restrictions} and \code{restrict.rhs}.
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
  \item{rank}{total number of linear independent coefficients
      = number of coefficients minus number of linear restrictions.}
  \item{df.residual}{degrees of freedom of the whole system.}
  \item{iter}{number of iteration steps.}
  \item{coefficients}{vector of all estimated coefficients.}
  \item{bcov}{estimated covariance matrix of \code{coefficients}.}
  \item{rcov}{estimated residual covariance matrix.}
  \item{rcovest}{residual covariance matrix used for estimation (only SUR and 3SLS).}
  \item{restrictions}{the restriction matrix,
      or a character vector giving the hypothesis in symbolic form
      (see \code{\link[car]{linear.hypothesis}} in package "car").}
  \item{restrict.rhs}{the restriction vector.}
  \item{TX}{matrix used to transform the regressor matrix.}
  \item{control}{list of control parameters used for the estimation.}

  ## elements of the class systemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnLabel*}{the equation label of the ith equation (from the labels list).}
  \item{terms*}{the 'terms' object used for the ith equation.}
  \item{inst*}{instruments of the ith equation (only 2SLS and 3SLS).}
  \item{rank*}{number of linear independent coefficients in the ith equation
     (differs from the number of coefficients only if there are
     restrictions that are not cross-equation).}
  \item{nCoef.sys*}{total number of coefficients in all equations.}
  \item{rank.sys*}{total number of linear independent coefficients
      in all equations.}
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
Rrestr[1,3] <-  1
Rrestr[1,7] <- -1
Rrestr[2,2] <- -1
Rrestr[2,5] <-  1
qrestr <- c( 0, 0.5 )
fitols2 <- systemfit( system, data = Kmenta,
                      restrictions = Rrestr, restrict.rhs = qrestr )
print( fitols2 )

## OLS estimation with the same 2 restrictions in symbolic form
restrict <- c( "demand_income - supply_trend = 0",
   "- demand_price + supply_price = 0.5" )
fitols2b <- systemfit( system, data = Kmenta, restrictions = restrict )
print( fitols2b )

# test whether both restricted estimators are identical
all.equal( coef( fitols2 ), coef( fitols2b ) )

## OLS with parameter restriction by argument TX
coefMap <- matrix( 0, 7, 6 )
colnames( coefMap ) <- c( "demIntercept", "demPrice", "demIncome",
   "supIntercept", "supPrice2", "supTrend" )
coefMap[ 1, "demIntercept" ] <- 1
coefMap[ 2, "demPrice" ]     <- 1
coefMap[ 3, "demIncome" ]    <- 1
coefMap[ 4, "supIntercept" ] <- 1
coefMap[ 5, "supPrice2" ]    <- 1
coefMap[ 6, "supPrice2" ]    <- 1
coefMap[ 7, "supTrend" ]     <- 1
fitols3 <- systemfit( system, data = Kmenta, TX = coefMap )
print( fitols3 )


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


## Examples how to use systemfit() with panel-like data
## Repeating the OLS and SUR estimations in Theil (1971, pp. 295, 300)
data( "GrunfeldTheil" )
library( plm )
pdata.frame( GrunfeldTheil, "firm", "year" )
formulaGrunfeld <- invest ~ value + capital
# OLS
theilOls <- systemfit( formulaGrunfeld, "OLS",
   data = GrunfeldTheil )
summary( theilOls )
# SUR
theilSur <- systemfit( formulaGrunfeld, "SUR",
   data = GrunfeldTheil, methodRCov = "noDfCor" )
summary( theilSur )


## Further examples are in the documentation to the data sets
## 'KleinI' and 'GrunfeldGreene'.
}

\keyword{models}
\keyword{regression}



