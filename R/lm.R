heavyLm <-
function(formula, data, family = Student(df = 4), subset, na.action,
  control = heavy.control(), model = TRUE, x = FALSE, y = FALSE,
  contrasts = NULL)
{
  ret.x <- x
  ret.y <- y
  Call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$family <- mf$control <- mf$model <- mf$x <- mf$y <- mf$contrasts <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Terms <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(Terms, mf, contrasts)
  ny <- NCOL(y)
  is.multivariate <- is.matrix(y) && (ny > 1)

  ## call fitter
  z <- if (is.multivariate) heavyMLm.fit(x, y, family, control)
       else heavyLm.fit(x, y, family, control)

  ## output
  z$call <- Call
  z$na.action <- attr(mf, "na.action")
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(Terms, mf)
  z$terms <- Terms
  if (model)
    z$model <- mf
  if (ret.y)
    z$y <- y
  if (ret.x)
    z$x <- x
  class(z) <- if (is.multivariate) "heavyMLm" else "heavyLm"
  z
}

heavyLm.fit <-
function(x, y, family = Student(df = 4), control = heavy.control())
{
  ## validating arguments
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L)
    stop("trying to fit a null model")
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (NROW(y) != n)
    stop("incompatible dimensions")
  dx <- c(n,p)

  ## extract family info
  if (!inherits(family, "heavy.family"))
    stop("Use only with 'heavy.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$which
  if ((kind < 0) || (kind > 4))
    stop("not valid 'family' object")
  settings <- c(kind, family$npars, unlist(family$pars))

  ## initial estimates
  fit <- lsfit(x, y, intercept = FALSE)[1:2]
  resid <- fit$residuals
  coef <- fit$coefficients
  sigma2 <- sum(resid^2) / n

  ## set control values
  if (missing(control))
    control <- heavy.control()
  if (!control$algorithm)
    control$ncycles <- 1
  ctrl <- unlist(control)[1:4]
  ctrl <- c(ctrl, 0)

  ## call fitter
  now <- proc.time()
  fit <- .C("lm_fit",
            y = as.double(y),
            x = as.double(x),
            dims = as.integer(dx),
            settings = as.double(settings),
            coefficients = as.double(coef),
            sigma2 = as.double(sigma2),
            fitted = double(n),
            residuals = as.double(resid),
            distances = double(n),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            acov = double(p^2),
            control = as.double(ctrl))
  speed <- proc.time() - now

  coef <- fit$coefficients
  acov <- matrix(fit$acov, ncol = p)
  xnames <- colnames(x)
  if (is.null(xnames))
    xnames <- paste0("x", 1L:p)
  names(coef) <- xnames
  dimnames(acov) <- list(xnames, xnames)

  ## creating the output object
  out <- list(dims = dx, family = family, settings = fit$settings,
              coefficients = coef, sigma2 = fit$sigma2, fitted.values = fit$fitted,
              residuals = fit$residuals, weights = fit$weights, distances = fit$distances,
              acov = acov, logLik = fit$logLik, numIter = fit$control[5],
              control = control, speed = speed)
  if (!control$fix.shape) {
    if ((kind > 1) && (kind < 4)) {
      df <- signif(out$settings[3], 6)
      out$family$call <- call(out$family$family, df = df)
    }
  }
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  else
    out$converged <- FALSE
  out
}

heavyMLm.fit <-
function(x, y, family = Student(df = 4), control = heavy.control())
{
  ## validating arguments
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L)
    stop("trying to fit a null model")
  ny <- NCOL(y)
  if (NROW(y) != n)
    stop("incompatible dimensions")
  dx <- c(n,p,ny)

  ## extract family info
  if (!inherits(family, "heavy.family"))
    stop("Use only with 'heavy.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$which
  if ((kind < 0) || (kind > 4))
    stop("not valid 'family' object")
  settings <- c(kind, family$npars, unlist(family$pars))

  ## initial estimates
  fit <- lsfit(x, y, intercept = FALSE)[1:2]
  resid <- fit$residuals
  coef <- fit$coefficients
  Sigma <- crossprod(resid) / n

  ## set control values
  if (missing(control))
    control <- heavy.control()
  if (!control$algorithm)
    control$ncycles <- 1
  ctrl <- unlist(control)[1:4]
  ctrl <- c(ctrl, 0)

  ## call fitter
  now <- proc.time()
  fit <- .C("mlm_fit",
            y = as.double(y),
            x = as.double(x),
            dims = as.integer(dx),
            settings = as.double(settings),
            coefficients = as.double(coef),
            Sigma = as.double(Sigma),
            fitted = double(n * ny),
            residuals = as.double(resid),
            distances = double(n),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            acov = double(p^2),
            control = as.double(ctrl))
  speed <- proc.time() - now

  coef <- matrix(fit$coefficients, ncol = ny)
  Sigma <- matrix(fit$Sigma, ncol = ny)
  fitted <- matrix(fit$fitted, ncol = ny)
  resid <- matrix(fit$resid, ncol = ny)
  acov <- matrix(fit$acov, ncol = p)
  xnames <- colnames(x)
  if (is.null(xnames))
    xnames <- paste0("x", 1L:p)
  ynames <- colnames(y)
  dimnames(coef) <- list(xnames, ynames)
  dimnames(Sigma) <- list(ynames, ynames)
  dimnames(acov) <- list(xnames, xnames)
  colnames(fitted) <- ynames
  colnames(resid) <- ynames

  ## creating the output object
  out <- list(dims = dx, family = family, settings = fit$settings,
              coefficients = coef, Sigma = Sigma, fitted.values = fitted,
              residuals = resid, weights = fit$weights, distances = fit$distances,
              acov = kronecker(acov, Sigma), logLik = fit$logLik,
              numIter = fit$control[5], control = control, speed = speed)
  if (!control$fix.shape) {
    if ((kind > 1) && (kind < 4)) {
      df <- signif(out$settings[3], 6)
      out$family$call <- call(out$family$family, df = df)
    }
  }
  if (out$numIter < control$maxIter)
    out$converged <- TRUE
  else
    out$converged <- FALSE
  out
}

print.heavyLm <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded\n")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$sigma2), "\n")
  invisible(x)
}

print.heavyMLm <-
function(x, digits = 4, ...)
{
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded\n")
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  nobs <- x$dims[1]
  rdf <- nobs - x$dims[2]
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\n")
  invisible(x)
}

summary.heavyLm <-
function (object, ...)
{
  z <- object
  se <- sqrt(diag(z$acov))
  est <- z$coefficients
  zval <- est / se
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$family <- z$family
  ans$logLik <- z$logLik
  ans$sigma2 <- z$sigma2
  ans$residuals <- z$residuals
  ans$coefficients <- cbind(est, se, zval, 2 * pnorm(abs(zval), lower.tail = FALSE))
  dimnames(ans$coefficients) <- list(names(z$coefficients),
        c("Estimate", "Std.Error", "Z value", "p-value"))
  ans$correlation <- z$acov / outer(se, se)
  dimnames(ans$correlation) <- dimnames(ans$coefficients)[c(1,1)]
  class(ans) <- "summary.heavyLm"
  ans
}

summary.heavyMLm <-
function (object, ...)
{
  z <- object
  ans <- z[c("call", "terms")]
  ans$dims <- z$dims
  ans$family <- z$family
  ans$logLik <- z$logLik
  ans$coefficients <- z$coefficients
  ans$Sigma <- z$Sigma
  ans$residuals <- z$residuals
  ans$acov <- z$acov
  ans$correlation <- cov2cor(z$acov)
  class(ans) <- "summary.heavyMLm"
  ans
}

print.summary.heavyLm <-
function(x, digits = 4, ...)
{
  cat("Linear model under heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  resid <- x$residuals
  nobs <- x$dims[1]
  p <- x$dims[2]
  rdf <- nobs - p
  if (rdf > 5) {
    cat("\nResiduals:\n")
		rq <- quantile(resid)
		names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
		print(rq, digits = digits, ...)
	}
	else if(rdf > 0) {
	 cat("\nResiduals:\n")
	 print(resid, digits = digits, ...)
  }
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nScale estimate:", format(x$sigma2))
  cat("\nLog-likelihood:", format(x$logLik), "on", p + 1, "degrees of freedom\n")
  invisible(x)
}

print.summary.heavyMLm <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Multivariate regression under heavy-tailed distributions\n")
  cat(" Data:", paste(as.name(x$call$data), ";", sep = ""))
  print(x$family)
  nobs <- x$dims[1]
  p <- x$dims[2]
  rdf <- nobs - p
  ny <- x$dims[3]
  npars <- ny * p + ny * (ny + 1) / 2
  cat("\nCoefficients:\n ")
  print(format(round(x$coef, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  print.symmetric(x$Sigma, digits = digits)
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual")
  cat("\nLog-likelihood:", format(x$logLik), "on", npars, "degrees of freedom\n")
  invisible(x)
}
