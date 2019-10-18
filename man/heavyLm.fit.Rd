\name{heavyLm.fit}
\alias{heavyLm.fit}
\alias{heavyMLm.fit}

\title{Fitter functions for linear models under heavy-tailed distributions}

\description{
  These are basic computing functions called by \code{\link{heavyLm}} used to
  fit linear models considering heavy-tailed errors. These should usually
  \emph{not} be used directly unless by experienced users.
}

\usage{
heavyLm.fit(x, y, family = Student(df = 4), control = heavy.control())

heavyMLm.fit(x, y, family = Student(df = 4), control = heavy.control())
}

\arguments{
  \item{x}{ design matrix of dimension \code{n}-by-\code{p}.}
  \item{y}{ vector of observations of length \code{n}, or a matrix with
    \code{n} rows.}
  \item{family}{ a description of the error distribution to be used in the model.
    By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{control}{ a list of control values for the estimation algorithm to replace
    the default values returned by the function \code{\link{heavy.control}}.}
}
\value{
  a list with components
  \item{family}{ the \code{\link{heavy.family}} object used, with the estimated
    shape parameters (if requested).}
  \item{coefficients}{ \code{p} vector or a matrix with \code{p} rows.}
  \item{sigma2}{ scale estimate of the random error (only available for univariate
    regression models).}
  \item{Sigma}{ estimate of scatter matrix for each row of the response matrix
    (only available for objects of class \code{"heavyMLm"}).}
  \item{residuals}{ \code{n} vector or matrix}
  \item{fitted.values}{ \code{n} vector or matrix}
  \item{weights}{ estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{ squared of scaled residuals or Mahalanobis distances.}
  \item{acov}{ asymptotic covariance matrix of the coefficients estimates.}
  \item{logLik}{  the log-likelihood at convergence.}
}
\seealso{
  \code{\link{heavyLm}} which you should use for multivariate or univariate linear
  regression under heavy-tailed distributions unless you know better.
}
\examples{
# univariate linear regression
data(ereturns)
x <- cbind(1, ereturns$CRSP)
colnames(x) <- c("Intercept", "CRSP")
y <- ereturns$m.marietta
z <- heavyLm.fit(x = x, y = y)

# multivariate linear regression
data(dialyzer)
y <- as.matrix(dialyzer[,1:4])
n <- nrow(y)
x <- matrix(1, nrow = n, ncol = 1) # a vector of ones
z <- heavyMLm.fit(x = x, y = y)
}
\keyword{regression}
\keyword{array}
