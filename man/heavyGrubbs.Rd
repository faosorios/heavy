\name{heavyGrubbs}
\alias{heavyGrubbs}

\title{Estimation for the Grubbs' model under heavy-tailed distributions}
\description{
  This function fits a Grubbs' model (using the parameterization of Theobald and
  Mallison, 1978) considering heavy-tailed distributions. This provides some degree
  of robustness to outliers.
}
\usage{
heavyGrubbs(y, data, family = Student(df = 4), subset, na.action, control)
}
\arguments{
  \item{y}{  a formula or a numeric matrix or an object that can be coerced to a numeric matrix.}
  \item{data}{  an optional data frame (or similar: see \code{\link{model.frame}}), used only if
  \code{y} is a formula. By default the variables are taken from \code{environment(formula)}.}
  \item{family}{ a description of the error distribution to be used in the model.
  By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{subset}{ an optional expression indicating the subset of the rows of
  data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace
  the default values returned by the function \code{\link{heavy.control}}.}
}

\value{
  A list with class \code{"heavyGrubbs"} containing the following components:
  \item{call}{  a list containing an image of the \code{heavyFit} call that produced the object.}
  \item{family}{  the \code{\link{heavy.family}} object used, with the estimated shape parameters (if requested).}
  \item{center}{  final estimate of the center parameters (related with the additive biases).}
  \item{phi}{  final estimate of the dispersion parameters.}
  \item{z}{  estimated latent variables.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  \item{acov}{  asymptotic covariance matrix of the center estimate.}
}

\references{
  Osorio, F., Paula, G.A., Galea, M. (2009).
  On estimation and influence diagnostics for the Grubbs' model under heavy-tailed distributions.
  \emph{Computational Statistics and Data Analysis} \bold{53}, 1249-1263.
}
\examples{
data(thermocouples)
fit <- heavyGrubbs(100 * thermocouples, family = Student(df = 4))
fit
}
\keyword{regression}
\keyword{models}
