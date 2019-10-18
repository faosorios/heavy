\name{heavyPS}
\alias{heavyPS}

\title{Fit a penalized spline under heavy-tailed distributions}
\description{
  Fits a penalized spline to the supplied data.
}
\usage{
heavyPS(x, y, family = Student(df = 4), nseg = 20, deg = 3, ord = 2,
  lambda = 1, method = c("GCV", "none"), ngrid = 200, control)
}

\arguments{
  \item{x, y}{ vectors giving the coordinates of the points in the scatter plot.
  Missing values are deleted.}
  \item{family}{ a description of the error distribution to be used in the model.
  By default the Student-t distribution with 4 degrees of freedom is considered.}
  \item{nseg}{ number of segments used to divide the domain, this information is
  required to construct the sequence of knots. Default value is 20.}
  \item{deg}{ the degree of the spline transformation. Must be a nonnegative integer.
  The default value is 3. The polynomial degree should be a small integer, usually 0,
  1, 2, or 3. Larger values are rarely useful.}
  \item{ord}{ the order of the roughness penalty. Default value is 2.}
  \item{lambda}{ specifies the smoothing parameter for the fit. It is fixed if \code{method="none"}.
  If \code{method="GCV"} then weighted generalized cross validation is used to choose an "optimal"
  smoothing parameter. The default value is 1.}
  \item{method}{ the method for choosing the smoothing parameter \code{lambda}. If \code{method="none"},
  then \code{lambda} is 'fixed'. If \code{method="GCV"} (the default) then the smoothing parameter is chosen
  automatically using the weighted generalized cross validation criterion.}
  \item{ngrid}{ number of elements in the grid used to compute the smoother.
  Only required to plot the fitted P-spline.}
  \item{control}{ a list of control values for the estimation algorithm to replace
  the default values returned by the function \code{\link{heavy.control}}.}
}

\value{
  an object of class \code{heavyPS} representing the fitted model. Generic
  functions \code{print} and \code{summary}, show the results of the fit.

  The following components must be included in a legitimate \code{heavyPS} object.
  \item{call}{  a list containing an image of the \code{heavyPS} call that produced the object.}
  \item{design}{  a list containing the B-spline basis matrix, the triangular factor of the penalty matrix and a numeric vector of knot positions with non-decreasing values.}
  \item{method}{  one of "GCV" or "none", depending on the fitting criterion used.}
  \item{family}{  the \code{\link{heavy.family}} object used in the fitting process.}
  \item{coefficients}{  final estimate of the coefficients vector.}
  \item{scale}{  final scale estimate of the random error.}
  \item{lambda}{  estimated smoothing parameter for the model (if requested).}
  \item{fitted.values}{  fitted model predictions of expected value for each datum.}
  \item{residuals}{  the residuals for the fitted model.}
  \item{plogLik}{  the penalized log-likelihood at convergence.}
  \item{edf}{  the effective number of parameters.}
  \item{gcv}{  the minimized smoothing parameter selection score (weighted GCV).}
  \item{pen}{  the penalty term at convergence.}
  \item{numIter}{  the number of iterations used in the iterative algorithm.}
  \item{weights}{  estimated weights corresponding to the assumed heavy-tailed distribution.}
  \item{distances}{  squared of scaled residuals.}
  \item{xgrid}{  grid of x-values used to fit the P-spline.}
  \item{ygrid}{  estimated curve on the x-grid, required to plot the fitted P-spline.}
  \item{shape}{  estimated shape parameters, only available if requested.}
}

\references{
  Eilers, P.H.C., and Marx, B.D. (1996).
  Flexible smoothing using B-splines and penalties (with discussion).
  \emph{Statistical Science} \bold{11}, 89-121.

  Osorio, F. (2016).
  Influence diagnostics for robust P-splines using scale mixture of normal distributions.
  \emph{Annals of the Institute of Statistical Mathematics} \bold{68}, 589-619.
}

\examples{
data(life)
x <- life$income
y <- life$life
fit <- heavyPS(x, y, family = Student(df = 5), method = "GCV")
summary(fit)
par(pty = "s")
plot(x, y, xlab = "Per Capita Income", ylab = "Life Expectancy")
lines(fit$xgrid, fit$ygrid)
}
\keyword{smooth}
