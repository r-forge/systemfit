
% 	$Id$	


\name{nlsystemfit}
\alias{nlsystemfit}

\title{Nonlinear Equation System Estimation}

\description{
  Fits a set of structural nonlinear equations using Ordinary Least
  Squares (OLS), Seemingly Unrelated Regression (SUR),
  Two-Stage Least Squares (2SLS), Three-Stage Least Squares (3SLS).
}


\usage{ nlsystemfit( method="OLS", eqns, startvals, eqnlabels, inst,
data, solvtol=.Machine$double.eps, pl=0, maxiter=1000 ) }


\arguments{
  \item{method}{the estimation method, one of "OLS", "SUR", "2SLS",
    "3SLS".}
  \item{eqns}{a list of structural equations to be estimated.}
  \item{startvals}{a list of starting parameter values for the minimization.}
  \item{eqnlabels}{an optional list of character vectors of names
    for the equation labels.}
  \item{inst}{one-sided model formula specifying instrumental variables
    or a list of one-sided model formulas if different instruments should
    be used for the different equations (only needed for 2SLS, 3SLS and
    GMM estimations).}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment from which
    nlsystemfit is called.}
  \item{pl}{print.level argument that is passed to the \code{\link{nlm}} function. }
  \item{solvtol}{tolerance level used in the step tolerance and gradient
    tolerances arguments for \code{\link{nlm}} and the tolerance for
    detecting linear dependencies in the columns of X in the
    \code{\link{qr}} function calls.}
  \item{maxiter}{the maximum number of iterations for the \code{\link{nlm}} function.}
%   \item{rcovformula}{formula to calculate the estimated residual covariance
%     matrix (see details).}
%   \item{probdfsys}{use the degrees of freedom of the whole system
%     (in place of the degrees of freedom of the single equation)
%     to calculate prob values for the t-test of individual parameters.}
%   \item{single.eq.sigma}{use different \eqn{\sigma^2}s for each
%     single equation to calculate the covariance matrix and the
%     standard errors of the coefficients (only OLS and 2SLS).}
%   \item{solvetol}{tolerance level for detecting linear dependencies
%     when inverting a matrix or calculating a determinant (see
%     \code{\link{solve}} and \code{\link{det}}).}
}

\details{

  The nlsystemfit function relies on the \code{\link{nlm}} perform the minimization of
  the objective functions and the \code{\link{qr}} set of functions.

  A system of nonlinear equations can be written as:

  \deqn{\epsilon_{t} = q( y_t, x_t, \theta )}
  \deqn{z_{t} = Z( x_t )}

  where \eqn{\epsilon_{t}} are the residuals from the y observations and
  the function evaluated at the parameter estimates. 

  The objective functions for the methods are:
  
%  \tabular{|l|c|c|c|}{
  \tabular{lccc}{
%    \hline
    Method \tab Instruments \tab Objective Function \tab Covariance of
    \eqn{\theta}\cr %\hline  
    OLS \tab no \tab \eqn{r'r} \tab \eqn{(X(diag(S)^{-1}\bigotimes
      I)X)^{-1}}\cr %\hline
    SUR \tab no \tab \eqn{r'(diag(S)_{OLS}^{-1}\bigotimes I)r} \tab
    \eqn{(X(S^{-1}\bigotimes I)X)^{-1}}\cr %\hline
    2SLS \tab yes \tab \eqn{r'(I \bigotimes W)r} \tab
    \eqn{(X(diag(S)^{-1}\bigotimes I)X)^{-1}}\cr %\hline 
    3SLS \tab yes \tab \eqn{r'(S_{2SLS}^{-1} \bigotimes W)r} \tab
    \eqn{(X(diag(S)^{-1}\bigotimes W)X)^{-1}} %\hline
  }

  where, r is a column vector for the residuals for each equation, S is
  variance-covariance matrix between the equations
  (\eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
    \sqrt{(T - k_i)*(T - k_j)}}), X is matrix of the
  partial derivates with respect to the parameters, W is a matrix of the
  instrument variables \eqn{Z(Z'Z)^{-1}Z}, Z is a matrix of the
  instrument variables, and I is an nxn identity matrix.

  The SUR and 3SLS methods requires two solutions. The first solution
  for the SUR is an OLS solution to obtain the variance-covariance
  matrix. The 3SLS uses the variance-covatiance from a 2SLS solution,
  then fits all the equations simultaneously.  

  The user should be aware that the function is \bold{VERY} sensative to
  the starting values and the nlm function may not converge. The nlm
  function will be called with the \code{typsize} argument set the
  absolute values of the starting values for the OLS and 2SLS
  methods. For the SUR and 3SLS methods, the \code{typsize} argument is
  set to the absolute values of the resulting OLS and 2SLS parameter
  estimates from the nlm result structre. In addition, the starting
  values for the SUR and 3SLS methods are obtained from the OLS and 2SLS
  parameter estimates to shorten the number of iterations. The number of
  iterations reported in the summary are only those used in the last
  call to nlm, thus the number of iterations in the OLS portion of the
  SUR fit and the 2SLS portion of the 3SLS fit are not included. 
  
%   \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
%     \sqrt{(T - k_i)*(T - k_j)}}
  
%   The formula to calculate the estimated covariance matrix of the
%   residuals, S,
%   (\eqn{\hat{\Sigma}}) can be one of the following (see Judge et al., 1985, p. 469): \cr
%   if rcovformula=0: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / T}; \cr
%   if rcovformula=1: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
%     \sqrt{(T - k_i)*(T - k_j)}}; \cr
%   if rcovformula=2: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / (T - k_i - k_j
%     + tr[(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j'X_i]}. \cr
%   If \eqn{k_i = k_j}, formula 1 and 2 are equal and yield an unbiased estimator for the
%   residual covariance matrix.
%   If \eqn{k_i \neq k_j}, only formula 2 yields an unbiased estimator for the residual
%   covariance matrix, but it is not neccessarily positive semidefinit and
%   its inverse is \bold{not} an unbiased estimator for the inverse of the residual
%   covariance matrix. Thus, it is doubtful whether formula 2 is really superior
%   to formula 1 (see Theil, 1971, p. 322).

  
  

  
%   The matrix \code{TX} transforms the regressor matrix (\eqn{X}) by
%    \eqn{X^{*} = X *} \code{TX}. Thus, the vector of coefficients is now
%    \eqn{b =} \code{TX} \eqn{\cdot b^{*}} , where \eqn{b} is the original (stacked) vector
%    of all coefficients and \eqn{b^{*}} is the new coefficient vector that is
%    estimated instead. Thus, the elements of vector \eqn{b} are
%    \eqn{b_i = \sum_j TX_{ij} \cdot b^{*}_j} \cr
%    The \code{TX} matrix can be used to change the order of the
%    coefficients and also to restrict coefficients (if \code{TX} has less
%    columns than it has rows). However restricting coefficients
%    by the \code{TX} matrix is less powerfull and flexible than the
%    restriction by providing the \code{R.restr} matrix and the
%    \code{q.restr} vector. The advantage of restricting the coefficients
%    by the \code{TX} matrix is that the matrix that is inverted for
%    estimation gets smaller by this procedure, while it gets larger
%    if the restrictions are imposed by \code{R.restr} and \code{q.restr}.

%    If iterated (WLS, SUR, W2SLS or 3SLS estimation with \code{maxit}>1),
%      the convergence criterion is \cr
%    \eqn{\sqrt{\sum_i (b_{i,g} - b_{i,g-1})^2 \left/ \sum_i b_{i,g-1}^2 \right.}}
%    < \code{tol} (\eqn{b_{i,g}} is the ith coefficient of the gth iteration step).

%    The formula to calculate the estimated covariance matrix of the residuals
%    (\eqn{\hat{\Sigma}}) can be one of the following (see Judge et al., 1985, p. 469): \cr
%    if rcovformula=0: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / T}; \cr
%    if rcovformula=1: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) /
%        \sqrt{(T - k_i)*(T - k_j)}}; \cr
%    if rcovformula=2: \eqn{\hat{\sigma}_{ij} = (\hat{e}_i' \hat{e}_j) / (T - k_i - k_j
%    + tr[(X_i'X_i)^{-1}X_i'X_j(X_j'X_j)^{-1}X_j'X_i]}. \cr
%    If \eqn{k_i = k_j}, formula 1 and 2 are equal and yield an unbiased estimator for the
%    residual covariance matrix.
%    If \eqn{k_i \neq k_j}, only formula 2 yields an unbiased estimator for the residual
%    covariance matrix, but it is not neccessarily positive semidefinit and
%    its inverse is \bold{not} an unbiased estimator for the inverse of the residual
%    covariance matrix. Thus, it is doubtful whether formula 2 is really superior
%    to formula 1 (see Theil, 1971, p. 322).

%    The formulas to calculate the 3SLS estimator lead to identical results
%    if the same instruments are used in all equations. If different instruments
%    are used in the different equations, only the GMM-3SLS estimator ("GMM")
%    and the 3SLS estimator proposed by Schmidt (1990) ("Schmidt") are consistent,
%    whereas "GMM" is efficient relative to "Schmidt" (see Schmidt, 1990).
}

\value{
  \code{nlsystemfit} returns a list of the class \code{nlsystemfit.system} and
  contains all results that belong to the whole system.
  This list contains one special object: "eq". It is a list and contains
  one object for each estimated equation. These objects are of the class
  \code{nlsystemfit.equation} and contain the results that belong only to the
  regarding equation.

  The objects of the class \code{nlsystemfit.system} and
  \code{nlsystemfit.equation} have the following components (the elements of
  the latter are marked with an asterisk (\eqn{*})):

  \item{eq}{a list object that contains a list object for each equation.}
  \item{method}{estimation method.}
  \item{resids}{an \eqn{n \times g} matrix of the residuals.}
  \item{g}{number of equations.}
  \item{n}{total number of observations.}
  \item{k}{total number of coefficients.}
  \item{b}{vector of all estimated coefficients.}
  \item{se}{estimated standard errors of \code{b}.}
  \item{t}{t values for \code{b}.}
  \item{p}{p values for \code{b}.}
  \item{bcov}{estimated covariance matrix of \code{b}.}
  \item{rcov}{estimated residual covariance matrix.}
  \item{drcov}{determinant of \code{rcov}.}
  \item{rcovest}{residual covariance matrix used for estimation (only SUR and 3SLS).}
  \item{rcor}{estimated residual correlation matrix.}
  \item{nlmest}{results from the nlm function call}
  
  %  \item{olsr2}{System OLS R-squared value.}
%  \item{mcelr2}{McElroys R-squared value for the system (only SUR and 3SLS).}
%  \item{y}{vector of all (stacked) endogenous variables}
%  \item{x}{matrix of all (diagonally stacked) regressors}
%  \item{h}{matrix of all (diagonally stacked) instrumental variables (only 2SLS and 3SLS)}
%  \item{data}{data frame of the whole system (including instruments)}
%  \item{R.restr}{the restriction matrix.}
%  \item{q.restr}{the restriction vector.}
%  \item{TX}{matrix used to transform the regressor matrix.}
%  \item{maxiter}{maximum number of iterations.}
%  \item{tol}{tolerance level indicating when to stop the iteration}
%  \item{rcovformula}{formula to calculate the estimated residual covariance
%    matrix}
%  \item{formula3sls}{formula for calculating the 3SLS estimator.}
%  \item{probdfsys}{system degrees of freedom to calculate prob values?.}
%  \item{single.eq.sigma}{different \eqn{\sigma^2}s for each single equation?.}
  \item{solvetol}{tolerance level when inverting a matrix or calculating a determinant.}

  ## elements of the class nlsystemfit.eq
  \item{eq}{a list that contains the results that belong to the individual equations.}
  \item{eqnlabel*}{the equation label of the ith equation (from the labels list).}
  \item{formula*}{model formula of the ith equation.}
%  \item{inst*}{instruments of the ith equation (only 2SLS and 3SLS).}
  \item{n*}{number of observations of the ith equation.}
  \item{k*}{number of coefficients/regressors in the ith equation (including the constant).}
%   \item{ki*}{number of linear independent coefficients in the ith equation (including the constant
%               differs from \code{k} only if there are restrictions that are not cross-equation).}
  \item{df*}{degrees of freedom of the ith equation.}
  \item{b*}{estimated coefficients of the ith equation.}
  \item{se*}{estimated standard errors of \code{b}.}
  \item{t*}{t values for \code{b}.}
  \item{p*}{p values for \code{b}.}
  \item{covb*}{estimated covariance matrix of \code{b}.}
%  \item{y*}{vector of endogenous variable (response values) of the ith equation.}
%  \item{x*}{matrix of regressors (model matrix) of the ith equation.}
%  \item{h*}{matrix of instrumental variables of the ith equation (only 2SLS and 3SLS).}
%  \item{data*}{data frame (including instruments) of the ith equation.}
  \item{predicted*}{vector of predicted values of the ith equation.}
  \item{residuals*}{vector of residuals of the ith equation.}
  \item{ssr*}{sum of squared residuals of the ith equation.}
  \item{mse*}{estimated variance of the residuals (mean of squared errors) of the ith equation.}
  \item{s2*}{estimated variance of the residuals (\eqn{\hat{\sigma}^2}) of the ith equation.}
  \item{rmse*}{estimated standard error of the residulas (square root of mse) of the ith equation.}
  \item{s*}{estimated standard error of the residuals (\eqn{\hat{\sigma}}) of the ith equation.}
  \item{r2*}{R-squared (coefficient of determination).}
  \item{adjr2*}{adjusted R-squared value.}
}

\references{

  Gallant, R. H. (1987)
  \emph{Nonlinear Equation Estimation}, John Wiley and Sons, 610 pp.

  SAS Institute (1999)
  \emph{SAS/ETS User's Guide, Version 8}, Cary NC: SAS Institute 1546 pp.

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com} }
}

\seealso{\code{\link{nlm}},\code{\link{qr}}, and \code{\link{systemfit}} }

\examples{
library( systemfit )
data( ritchie )

hg.formula <- hg ~ exp( h0 + h1*log(ht) + h2*ht^2 + h3*cr + h4*whc )
dg.formula <- dg ~ exp( d0 + d1*log(bd) + d2*log(hg) + d3*bainlh )
cr.formula <- cr ~ 1.0 / ( 1.0 + exp( c0 + c1*log(bd) + c2*ht^2 + c3*bainlh ) )
labels <- list( "height.growth", "diameter.growth", "crown.ratio" )
inst <- ~ ht + bd + whc + bainlh
start.values <- c(h0=-1.0, h1=0.5, h2=-0.0005, h3=0.05, h4=-0.2,
                  d0=-5.0, d1=1.2, d2=0.5, d3=-0.002,
                  c0=-1.5, c1=-0.5, c2=0.0005, c3=0.002 )
model <- list( hg.formula, dg.formula, cr.formula )
model.ols <- nlsystemfit( "OLS", model, start.values, data=ritchie, eqnlabels=labels )

# print the results
print( model.ols )

}

\keyword{models}
\keyword{regression}
\keyword{nonlinear}



