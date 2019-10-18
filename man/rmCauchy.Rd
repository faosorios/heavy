\name{rmCauchy}
\alias{rmCauchy}

\title{Multivariate Cauchy Random Deviates}
\description{
  Random number generation from the multivariate Cauchy distribution.
}
\usage{
rmCauchy(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)))
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{center}{ a vector giving the positions of each variable}
  \item{Scatter}{ a positive-definite dispersion matrix}
}

\value{
  If \code{n = 1} a vector of the same length as \code{center}, otherwise a
  matrix of \code{n} rows of random vectors.
}

\details{
  The function \code{rmCauchy} is an interface to C routines, which make calls to
  subroutines from LAPACK. The matrix decomposition is internally done using
  the Cholesky decomposition. If \code{Scatter} is not non-negative definite then
  there will be a warning message.
}

\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  Springer-Verlag, New York.
}

\seealso{
  \code{\link{rcauchy}}
}

\examples{
# dispersion parameters
Scatter <- matrix(c(10,3,3,2), ncol = 2)
Scatter

# generate the sample
y <- rmCauchy(n = 1000, Scatter = Scatter)

# scatterplot of a random bivariate Cauchy sample with center
# vector zero and scale matrix 'Scatter'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate Cauchy sample", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
