pgamma.deriv <- function(x, shape, scale, deriv = 0:2)
{
  if (missing(shape))
    stop("argument 'shape' is missing.")

  if (length(shape) > 1)
    stop("'shape' is not allowed to be a vector.")

  if (length(scale) > 1)
    stop("'scale' is not allowed to be a vector.")

  if (shape <= 0.0)
    stop("'shape' must be non-negative.")

  if (scale <= 0.0)
    stop("'scale' must be non-negative.")

  which <- deriv
  ok <- x > 0
  x <- x[ok]

  z <- .C("cdf_gamma_derivatives",
          x = as.double(x),
          shape = as.double(shape),
          scale = as.double(scale),
          deriv = double(3))$deriv

  idx <- pmatch(which, 0:2)
  ok  <- complete.cases(idx)
  idx <- idx[ok]

  z[idx]
}
