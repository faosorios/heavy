\name{rmt}
\alias{rmt}

\title{Multivariate Student-t Random Deviates}
\description{
  Random number generation from the multivariate Student-t distribution.
}
\usage{
rmt(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)),
  df = 4)
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{center}{ a vector giving the locations of each variable}
  \item{Scatter}{ a positive-definite dispersion matrix}
  \item{df}{ degrees of freedom (positive, maybe non-integer)}
}

\value{
  If \code{n = 1} a vector of the same length as \code{center}, otherwise a
  matrix of \code{n} rows of random vectors.
}

\details{
  The function \code{rmt} is an interface to C routines, which make calls to
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
  \code{\link{rt}}
}

\examples{
# dispersion parameters
Scatter <- matrix(c(10,3,3,2), ncol = 2)
Scatter

# generate the sample
y <- rmt(n = 1000, Scatter = Scatter)

# scatterplot of a random bivariate Student-t sample with center
# vector zero and scale matrix 'Scatter'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate Student-t sample (df = 4)", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
