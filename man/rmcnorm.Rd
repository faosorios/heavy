\name{rmcnorm}
\alias{rmcnorm}

\title{Multivariate Contaminated Normal Random Deviates}
\description{
  Random number generation from the multivariate contaminated normal distribution.
}
\usage{
rmcnorm(n = 1, center = rep(0, nrow(Scatter)), Scatter = diag(length(center)),
  epsilon = 0.05, vif = 0.25)
}

\arguments{
  \item{n}{ the number of samples requested}
  \item{center}{ a vector giving the locations of each variable}
  \item{Scatter}{ a positive-definite dispersion matrix}
  \item{epsilon, vif}{ percentage of outliers contamination and variance inflation factor}
}

\value{
  If \code{n = 1} a vector of the same length as \code{center}, otherwise a
  matrix of \code{n} rows of random vectors.
}

\details{
  The function \code{rmcnorm} is an interface to C routines, which make calls to
  subroutines from LAPACK. The matrix decomposition is internally done using
  the Cholesky decomposition. If \code{Scatter} is not non-negative definite then
  there will be a warning message.
}

\references{
  Devroye, L. (1986).
  \emph{Non-Uniform Random Variate Generation}.
  Springer-Verlag, New York.
}

\examples{
# dispersion parameters
Scatter <- matrix(c(10,3,3,2), ncol = 2)
Scatter

# generate the sample
y <- rmcnorm(n = 1000, Scatter = Scatter, epsilon = .05, vif = .25)

# scatterplot of a random bivariate contaminated normal sample
# with center vector zero and scale matrix 'Scatter'
par(pty = "s")
plot(y, xlab = "", ylab = "")
title("bivariate contaminated normal sample", font.main = 1)
}
\keyword{distribution}
\keyword{multivariate}
