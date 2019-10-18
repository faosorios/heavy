\name{heavyFit}
\alias{heavyFit}

\title{Estimation of Location and Scatter using multivariate heavy-tailed distributions}
\description{
  Estimates the center vector and scatter matrix assuming that the data came from
  a multivariate heavy-tailed distribution. This provides some degree of robustness
  to outliers without giving a high breakdown point.
}
\usage{
heavyFit(x, data, family = Student(df = 4), subset, na.action, control)
}
\arguments{
  \item{x}{  a formula or a numeric matrix or an object that can be coerced to a numeric matrix.}
  \item{data}{  an optional data frame (or similar: see \code{\link{model.frame}}), used only if
  \code{x} is a formula. By default the variables are taken from \code{environment(formula)}.}
  \item{family}{ a description of the error distribution to be used in the model.
  By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{subset}{ an optional expression indicating the subset of the rows of
  data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace
  the default values returned by the function \code{\link{heavy.control}}.}
}

\value{
  A list with class \code{"heavyFit"} containing the following components:
  \item{call}{  a list containing an image of the \code{heavyFit} call that produced the object.}
  \item{family}{  the \code{\link{heavy.family}} object used, with the estimated shape parameters (if requested).}
  \item{center}{  final estimate of the location vector.}
  \item{Scatter}{  final estimate of the scale matrix.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  \item{acov}{  asymptotic covariance matrix of the center estimates.}
}

\references{
  Kent, J.T., Tyler, D.E., and Vardi, Y. (1994).
  A curious likelihood identity for the multivariate t-distribution.
  \emph{Communications in Statistics - Simulation and Computation} \bold{23}, 441-453.

  Lange, K., Little, R.J.A., and Taylor, J.M.G. (1989).
  Robust statistical modeling using the t distribution.
  \emph{Journal of the American Statistical Association} \bold{84}, 881-896.

  Little, R.J.A. (1988).
  Robust estimation of the mean and covariance matrix from data with missing values.
  \emph{Applied Statistics} \bold{37}, 23-38.
}
\seealso{
\code{\link{cov}}
}
\examples{
data(creatinine)
fit <- heavyFit(~ log(WT) + log(SC) + log(140 - Age) + log(CR), data = creatinine,
  family = Student(df = 4))
summary(fit)
}
\keyword{multivariate}
