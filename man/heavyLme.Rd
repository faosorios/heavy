\name{heavyLme}
\alias{heavyLme}
\alias{heavyLme.formula}

\title{Linear mixed-effects models under heavy-tailed distributions}
\description{
  This function fits a linear mixed-effects model under heavy-tailed errors
  using the formulation described in Pinheiro et al. (2001).
}
\usage{
heavyLme(fixed, random, groups, data, family = Student(df = 4),
  subset, na.action, control)
}

\arguments{
  \item{fixed}{ a two-sided linear formula object describing the fixed-effects part of the model,
  with the response on the left of a ~ operator and the terms, separated by + operators.}
  \item{random}{  a one-sided formula of the form ~x1+...+xn specifying the model for the random effects.}
  \item{groups}{ a one-sided formula for specifying the grouping variable.}
  \item{data}{ an optional data frame containing the variables named in \code{fixed}, \code{random} and \code{group}.
  By default the variables are taken from the environment from which heavy.lme is called.}
  \item{family}{ a description of the error distribution to be used in the model.
  By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{subset}{ an optional expression indicating the subset of the rows of data that should be used in the fitting process.}
  \item{na.action}{ a function that indicates what should happen when the data contain NAs.}
  \item{control}{ a list of control values for the estimation algorithm to replace the default
  values returned by the function \code{\link{heavy.control}}.}
}

\value{
  An object of class \code{heavyLme} representing the linear mixed-effects model fit. Generic function \code{print}
  and \code{summary}, show the results of the fit.

  The following components must be included in a legitimate \code{heavyLme} object.
  \item{lmeData}{  an object representing a list of mixed-effects model components.}
  \item{call}{  a list containing an image of the \code{heavyLme} call that produced the object.}
  \item{family}{  the \code{\link{heavy.family}} object used, with the estimated shape parameters (if requested).}
  \item{coefficients}{  final estimate of the fixed effects.}
  \item{theta}{  final estimate of the scale parameters associated to the random effects.}
  \item{scale}{  final scale estimate of the random error.}
  \item{logLik}{  the log-likelihood at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{ranef}{  a matrix with the estimated random effects.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  estimated squared Mahalanobis distances.}
  \item{Fitted}{  a data frame with the \code{"marginal"} and \code{"conditional"} fitted values as columns.}
  \item{Resid}{  a data frame with the \code{"marginal"} and \code{"conditional"} residuals as columns.}
}

\references{
  Pinheiro, J.C., Liu, C., and Wu, Y.N. (2001).
  Efficient algorithms for robust estimation in linear mixed-effects models using the multivariate \emph{t} distribution.
  \emph{Journal of Computational and Graphical Statistics} \bold{10}, 249--276.
}

\examples{
data(dental)
fm0 <- heavyLme(distance ~ age * Sex, random = ~ age, groups = ~ Subject,
  data = dental, family = Student(df = 4))
summary(fm0)

# fitting model with fixed degrees of freedom
fm1 <- heavyLme(distance ~ age * Sex, random = ~ age, groups = ~ Subject,
  data = dental, family = Student(df = 4), control = heavy.control(fix.shape = TRUE))
summary(fm1) # fixed at df = 4
}
\keyword{models}
