\name{pgamma.deriv}
\alias{pgamma.deriv}

\title{Derivatives of the regularized incomplete gamma function}

\description{
  Computes the incomplete gamma function and its first and second derivatives with
  respect to 'shape' parameter.
}
\usage{
  pgamma.deriv(x, shape, scale, deriv = 0:2)
}
\arguments{
  \item{x}{a numeric argument, \code{x > 0}.}
  \item{shape, scale}{shape and scale parameters, must be positive.}
  \item{deriv}{integer vector, order of the required derivatives.}
}
\value{
  Function \code{pgamma.deriv} returns the value, first and second derivatives of
  the regularized incomplete gamma gamma function with respect to 'shape' parameter.
}
\details{
  The regularized (standard) incomplete gamma function is given by
  \deqn{P(a, x) = \frac{1}{\Gamma(a)} \int_0^x t^{a - 1}\,e^{-t}\,dt,}
  which is related with the CDF of Gamma distribution (see \code{\link{pgamma}}).

  \code{pgamma.deriv} is a C translation of the original Fortran AS 187 subroutine
  by R.J. Moore.
}
\references{
  Abramowitz, M., and Stegun, I.A. (1970).
  \emph{Handbook of Mathematical Functions}.
  Dover, New York.

  Moore, R.J. (1982).
  Algorithm AS 187: Derivatives of the incomplete gamma integral.
  \emph{Applied Statistics} \bold{31}, 330-335.
}
\seealso{
  \code{\link{pgamma}}, \code{\link{psigamma}}
}
\examples{
  x <- 3
  pgamma.deriv(x, 1, 1, deriv = 0:2)

  # only to check 1st element
  pgamma(x, 1, 1)

}
\keyword{math}
