# ID: distn.R, last updated 2019/09/06, F.Osorio */

dtgamma <- function(x, shape, scale = 1, truncation = 1, log = FALSE)
{
  if (missing(shape))
    stop("argument 'shape' is missing.")

  if (any(shape <= 0.0))
    stop("'shape' must be non-negative.")

  if (any(scale <= 0.0))
    stop("'scale' must be non-negative.")

  if (any(truncation <= 0.0))
    stop("'truncation' must be non-negative.")

  x0 <- y0 <- x
  ok <- x >= 0 & x <= truncation
  x <- x[ok]

  n <- length(x)
  nshape <- length(shape)
  nscale <- length(scale)
  ntrunc <- length(truncation)

  y <- .C("pdf_tgamma",
          n = as.integer(n),
          y = double(n),
          x = as.double(x),
          shape  = as.double(shape),
          nshape = as.integer(nshape),
          scale  = as.double(scale),
          nscale = as.integer(nscale),
          truncation = as.double(truncation),
          ntrunc = as.integer(ntrunc),
          give.log = as.integer(log))$y
  y0[ok] <- y
  y0[!ok] <- 0.0
  y0
}

ptgamma <- function(q, shape, scale = 1, truncation = 1, lower.tail = TRUE)
{
  if (missing(shape))
    stop("argument 'shape' is missing.")

  if (any(shape <= 0.0))
    stop("'shape' must be non-negative.")

  if (any(scale <= 0.0))
    stop("'scale' must be non-negative.")

  if (any(truncation <= 0.0))
    stop("'truncation' must be non-negative.")

  n <- length(q)
  nshape <- length(shape)
  nscale <- length(scale)
  ntrunc <- length(truncation)

  y <- .C("cdf_tgamma",
          n = as.integer(n),
          y = double(n),
          q = as.double(q),
          shape  = as.double(shape),
          nshape = as.integer(nshape),
          scale  = as.double(scale),
          nscale = as.integer(nscale),
          truncation = as.double(truncation),
          ntrunc = as.integer(ntrunc),
          lower.tail = as.integer(lower.tail))$y
  y
}

qtgamma <- function(p, shape, scale = 1, truncation = 1, lower.tail = TRUE)
{
  if (missing(shape))
    stop("argument 'shape' is missing.")

  if (any(shape <= 0.0))
    stop("'shape' must be non-negative.")

  if (any(scale <= 0.0))
    stop("'scale' must be non-negative.")

  if (any(truncation <= 0.0))
    stop("'truncation' must be non-negative.")

  n <- length(p)
  nshape <- length(shape)
  nscale <- length(scale)
  ntrunc <- length(truncation)

  y <- .C("quantile_tgamma",
          n = as.integer(n),
          y = double(n),
          p = as.double(p),
          shape  = as.double(shape),
          nshape = as.integer(nshape),
          scale  = as.double(scale),
          nscale = as.integer(nscale),
          truncation = as.double(truncation),
          ntrunc = as.integer(ntrunc),
          lower.tail = as.integer(lower.tail))$y
  y
}

rtgamma <- function(n, shape, scale = 1, truncation = 1)
{
  if (missing(shape))
    stop("argument 'shape' is missing.")

  if (any(shape <= 0.0))
    stop("'shape' must be non-negative.")

  if (any(scale <= 0.0))
    stop("'scale' must be non-negative.")

  if (any(truncation <= 0.0))
    stop("'truncation' must be non-negative.")

  nshape <- length(shape)
  nscale <- length(scale)
  ntrunc <- length(truncation)

  x <- .C("rand_tgamma",
          n = as.integer(n),
          x = double(n),
          shape  = as.double(shape),
          nshape = as.integer(nshape),
          scale  = as.double(scale),
          nscale = as.integer(nscale),
          truncation = as.double(truncation),
          ntrunc = as.integer(ntrunc))$x
  x
}
