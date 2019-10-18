\name{heavy.family}
\alias{heavy.family}
\alias{normal}
\alias{Cauchy}
\alias{Student}
\alias{slash}
\alias{contaminated}

\title{Family objects for heavy-tailed distributions}
\description{
  Provide a convenient way to specify the details of the models used by the fitter functions.
  (i.e., \code{\link{heavyFit}}, \code{\link{heavyGrubbs}}, \code{\link{heavyLm}},
  \code{\link{heavyLme}} and \code{\link{heavyPS}} functions).
}
\usage{
normal()
Cauchy()
Student(df = 4)
slash(df = 2)
contaminated(epsilon = 0.05, vif = 0.25)
}

\arguments{
  \item{df}{  degrees of freedom for Student-t and slash distributions.}
  \item{epsilon, vif}{  percentage of outliers contamination and variance inflation factor.}
}
\details{
  \code{heavy.family} is a generic function to create information about the selected heavy-tailed distribution
  which is passed to the estimation algorithm.
}

\seealso{\code{\link{heavyFit}}, \code{\link{heavyGrubbs}}, \code{\link{heavyLm}}, \code{\link{heavyLme}}, \code{\link{heavyPS}}}

\examples{
CauchyFmly <- Student(df = 1)
CauchyFmly
}
\keyword{models}
