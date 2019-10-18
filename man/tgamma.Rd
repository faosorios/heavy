\name{tgamma}
\alias{tgamma}
\alias{dtgamma}
\alias{ptgamma}
\alias{qtgamma}
\alias{rtgamma}
\title{The right truncated gamma distribution}
\description{
  Density, distribution function, quantile function and random generation for the
  right truncated gamma distribution with shape (\code{shape}), scale (\code{scale})
  parameters and right truncation point (\code{truncation}).
}
\usage{
  dtgamma(x, shape, scale = 1, truncation = 1, log = FALSE)
  ptgamma(q, shape, scale = 1, truncation = 1, lower.tail = TRUE)
  qtgamma(p, shape, scale = 1, truncation = 1, lower.tail = TRUE)
  rtgamma(n, shape, scale = 1, truncation = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{shape, scale}{shape and scale parameters, must be positive.}
  \item{truncation}{right truncation point, must be positive.}
  \item{log}{logical; if TRUE, the log-density is returned.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \eqn{P[X \le x]},
    otherwise, \eqn{P[X > x]}.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is taken
    to be the number of required deviates.}
}
\value{
  \code{dtgamma}, \code{ptgamma}, and \code{qtgamma} are respectively the density,
  distribution function and quantile function of the right truncated gamma distribution.
  \code{rtgamma} generates random deviates from the right truncated gamma distribution.

  The length of the result is determined by \code{n} for \code{rtgamma}, and is
  the maximum of the lengths of the numerical parameters for the other functions.
}
\details{
  If \code{scale} or \code{truncation} are not specified, they assume the default
  values.

  The right truncated gamma distribution with shape \eqn{a}, scale \eqn{b} and
  right truncation point \eqn{t > 0} has density
  \deqn{
    f(x) =
    \frac{b^a}{\gamma(a,bt)} \exp(-bx)x^{a-1}}{
    f(x) = b^a/gamma(a,bt) exp(-bx)x^(a-1)
  }
  con \eqn{x < t} and \eqn{\gamma(a,b)} denotes the incomplete gamma function
  (see Abramowitz and Stegun, 1970, pp. 260).
}
\references{
  Abramowitz, M., and Stegun, I.A. (1970).
  \emph{Handbook of Mathematical Functions}.
  Dover, New York.

  Phillippe, A. (1997).
  Simulation of right and left truncated gamma distribution by mixtures.
  \emph{Statistics and Computing} \bold{7}, 173-181.
}
\seealso{
  \link{Distributions} for other standard distributions.
}
\examples{
x <- seq(0, 2, by = 0.1)
y <- dtgamma(x, shape = 1, truncation = 1)
z <- dgamma(x, shape = 1) # standard gamma pdf
plot(x, z, type = "l", xlab = "x", ylab = "density", ylim = range(y, z), lty = 2)
lines(x, y)

x <- rtgamma(1000, shape = 1)
## Q-Q plot for the right truncated gamma data against true theoretical distribution:
qqplot(qtgamma(ppoints(1000), shape = 1), x, main = "Truncated Gamma Q-Q plot",
  xlab = "Theoretical quantiles", ylab = "Sample quantiles", font.main = 1)
abline(c(0,1), col = "red", lwd = 2)
}
\keyword{distribution}
