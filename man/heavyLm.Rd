\name{heavyLm}
\alias{heavyLm}

\title{(Multivariate) linear regression under heavy-tailed distributions}
\description{
  This function is used to fit linear models considering heavy-tailed errors.
  It can be used to carry out univariate or multivariate regression.
}
\usage{
heavyLm(formula, data, family = Student(df = 4), subset, na.action,
  control, model = TRUE, x = FALSE, y = FALSE, contrasts = NULL)
}

\arguments{
  \item{formula}{ an object of class \code{"formula"}: a symbolic description of
  the model to be fitted.}
  \item{data}{ an optional data frame containing the variables in the model. If
  not found in \code{data}, the variables are taken from \code{environment(formula)},
  typically the environment from which \code{heavyLm} is called.}
  \item{family}{ a description of the error distribution to be used in the model.
  By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{subset}{ an optional expression indicating the subset of the rows of
  data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace
  the default values returned by the function \code{\link{heavy.control}}.}
  \item{model, x, y}{ logicals. If \code{TRUE} the corresponding components of
  the fit (the model frame, the model matrix, the response) are returned.}
  \item{contrasts}{ an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}
}

\details{
  Models for \code{heavyLm} are specified symbolically (for additional information see the "Details"
  section from \code{\link{lm}} function). If \code{response} is a matrix, then a multivariate linear
  model is fitted.
}

\value{
  An object of class \code{"heavyLm"} or \code{"heavyMLm"} for multiple responses
  which represents the fitted model. Generic functions \code{print} and \code{summary},
  show the results of the fit.

  The following components must be included in a legitimate \code{"heavyLm"} object.
  \item{call}{  a list containing an image of the \code{heavyLm} call that produced the object.}
  \item{family}{  the \code{\link{heavy.family}} object used, with the estimated shape parameters (if requested).}
  \item{coefficients}{  final estimate of the coefficients vector.}
  \item{sigma2}{  final scale estimate of the random error (only available for univariate regression models).}
  \item{Sigma}{  estimate of scatter matrix for each row of the response matrix (only available for objects of class \code{"heavyMLm"}).}
  \item{fitted.values}{  the fitted mean values.}
  \item{residuals}{  the residuals, that is response minus fitted values.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  squared of scaled residuals or Mahalanobis distances.}
  \item{acov}{  asymptotic covariance matrix of the coefficients estimates.}
}

\references{
  Dempster, A.P., Laird, N.M., and Rubin, D.B. (1980).
  Iteratively reweighted least squares for linear regression when errors are Normal/Independent distributed.
  In P.R. Krishnaiah (Ed.), \emph{Multivariate Analysis V}, p. 35-57. North-Holland.

  Lange, K., and Sinsheimer, J.S. (1993).
  Normal/Independent distributions and their applications in robust regression.
  \emph{Journal of Computational and Graphical Statistics} \bold{2}, 175-198.
}

\examples{
# univariate linear regression
data(ereturns)
fit <- heavyLm(m.marietta ~ CRSP, data = ereturns, family = Student(df = 5))
summary(fit)

# multivariate linear regression
data(dialyzer)
fit <- heavyLm(cbind(y1,y2,y3,y4) ~ -1 + centre, data = dialyzer, family = slash(df = 4))
fit

# fixing the degrees of freedom at df = 5
fit <- heavyLm(m.marietta ~ CRSP, data = ereturns, family = Student(df = 5),
  control = heavy.control(fix.shape = TRUE))
summary(fit)
}
\keyword{regression}
\keyword{models}
