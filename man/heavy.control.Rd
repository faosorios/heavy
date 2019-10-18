\name{heavy.control}
\alias{heavy.control}

\title{Set control parameters}
\description{
  Allows users to set control parameters for routines available in \code{heavy}.
}
\usage{
heavy.control(maxIter = 2000, tolerance = 1e-6, fix.shape = FALSE,
  ndraws = 500, algorithm = c("EM", "NEM"), ncycles = 5)
}

\arguments{
  \item{maxIter}{  maximum number of iterations. The default is 2000.}
  \item{tolerance}{  the relative tolerance in the iterative algorithm.}
  \item{fix.shape}{  whether the shape parameters should be kept fixed in the fitting
  processes. The default is \code{fix.shape = FALSE}. Currently, only is available for
  the Student-t and slash distributions.}
  \item{ndraws}{  independent draws for Monte Carlo integration. The default is 500.}
  \item{algorithm}{  character string specifying the algorithm to use. The default
  algorithm is an EM algorithm. Other possible value is "NEM" for the nested EM algorithm.}
  \item{ncycles}{  number of cycles used in the nested EM algorithm, not used if \code{algorithm="EM"}.
  The default is 5.}
}

\value{
  A list of control arguments to be used in a call to the fitter (i.e., \code{\link{heavyFit}},
  \code{\link{heavyGrubbs}}, \code{\link{heavyLm}}, \code{\link{heavyLme}} or \code{\link{heavyPS}}
  functions).

  A call to \code{heavy.control} can be used directly in the control argument of the call to the fitter.
}

\seealso{\code{\link{heavyFit}}, \code{\link{heavyGrubbs}}, \code{\link{heavyLm}}, \code{\link{heavyLme}}, \code{\link{heavyPS}}}

\examples{
ctrl <- heavy.control(maxIter = 200, tol = 1e-04, ndraws = 9000)
data(dental)
heavyLme(distance ~ age * Sex, random = ~ age, groups = ~ Subject,
  data = dental, control = ctrl)
}
\keyword{models}
