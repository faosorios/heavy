heavyGrubbs <-
function(y, data, family = Student(df = 4), subset, na.action, control = heavy.control())
{
  ## initial values using ML estimates under normality (Gleser, 1981)
  grubbs.start <- function(z) {
    o  <- cov.wt(z, method = "unbiased")
    rs <- eigen(o$cov)
    b  <- rs$vectors[,1]
    d1 <- rs$values
    d0 <- d1[1]
    d1 <- d1[-1]
    mu <- o$center
    p  <- length(mu) - 1
    b  <- b[-1] / b[1]
    c0 <- 1.0 + sum(b^2)
    phi <- mean(d1)
    scale <- (p * d0 - sum(d1)) / (p * c0)
    start <- list(center = mu, scale = scale, phi = rep(phi, p +1))
    start
  }

  Call <- match.call()
  if (missing(y))
    stop("'y' is not supplied")
  if (inherits(y, "formula")) {
    mt <- terms(y, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "y"] <- "formula"
    mf$family <- mf$control <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  }
  else {
    z <- as.matrix(y)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("heavyGrubbs applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  if (p < 2)
    stop("Grubbs model requires at least two variables")

  ## extract family info
  if (!inherits(family, "heavy.family"))
    stop("Use only with 'heavy.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$which
  if ((kind < 0) || (kind > 4))
    stop("not valid 'family' object")
  settings <- c(kind, family$npars, unlist(family$pars))

  ## set control values
  if (missing(control))
    control <- heavy.control()
  if (!control$algorithm)
    control$ncycles <- 1
  ctrl <- unlist(control)[1:4]
  ctrl <- c(ctrl, 0)

  ## initial estimates
  start <- grubbs.start(z)

  ## Call fitter
  now <- proc.time()
  fit <- .C("grubbs_fit",
            y = as.double(t(z)),
            dims = as.integer(dz),
            settings = as.double(settings),
            center = as.double(start$center),
            phi = as.double(start$phi),
            scale = as.double(start$scale),
            z = double(n),
            distances = double(n),
            weights = as.double(rep(1, n)),
            residuals = double(n * p),
            logLik = double(1),
            acov = double(p^2),
            control = as.double(ctrl))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              y = z,
              dims = dz,
              family = family,
              settings = fit$settings,
              center = fit$center,
              phi = fit$phi,
              scale = fit$scale,
              residuals = fit$residuals,
              z = fit$z,
              logLik = fit$logLik,
              numIter = fit$control[5],
              control = control,
              weights = fit$weights,
              distances = fit$distances,
              acov = matrix(fit$acov, ncol = p),
              speed = speed,
              start = start,
              converged = FALSE)
  names(out$center) <- znames
  names(out$phi) <- znames
  dimnames(out$acov) <- list(znames, znames)
  if (!control$fix.shape) {
    if ((kind > 1) && (kind < 4)) {
      df <- signif(out$settings[3], 6)
      out$family$call <- call(out$family$family, df = df)
    }
  }
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  class(out) <- "heavyGrubbs"
  out
}

print.heavyGrubbs <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded\n")
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nDispersion:\n ")
  print(format(round(x$phi, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}
