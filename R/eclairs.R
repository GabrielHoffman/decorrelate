# Gabriel Hoffman March 12, 2021 Eclairs: Estimate covariance/correlation with
# low rank and shrinkage


#' Class eclairs
#'
#' Class \code{eclairs}
#'
#' @details Object storing:
#' \describe{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{sigma: }{standard deviations of input columns}
#'  \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#'  \item{rownames: }{sample names from the original matrix}
#'  \item{colnames: }{features names from the original matrix}
#'  \item{method: }{method used for decomposition}
#'  \item{call: }{the function call}
#' }
#' @name eclairs-class
#' @rdname eclairs-class
#' @exportClass eclairs
# setClass('eclairs', representation('list'))
setClass("eclairs", contains = "list")


#' @importFrom methods show
setMethod("print", "eclairs", function(x) {
  show(x)
})

setMethod("show", "eclairs", function(object) {
  isCorrMatrix <- all(object$sigma == 1)

  if (isCorrMatrix) {
    cat("       Estimate correlation with low rank and shrinkage\n\n")
  } else {
    cat("       Estimate covariance with low rank and shrinkage\n\n")
  }

  cat("  Original data dims:", object$n, "x", object$p, "\n")
  cat("  Low rank component:", length(object$dSq), "\n")
  cat("  lambda:            ", format(object$lambda, digits = 3), "\n")
  cat("  nu:                ", format(object$nu, digits = 3), "\n")

  if (!isCorrMatrix) {
    # print standard deviations of features
    if (length(object$sigma) > 2) {
      values <- object$sigma[c(1, 2, length(object$sigma))]
      values <- format(values, digits = 3)
      values <- paste(values[1], values[2], "...", values[3])
    } else {
      values <- format(object$sigma, digits = 3)
      values <- paste(values, collapse = " ")
    }
    startText <- paste("  sd(", length(object$sigma), "): ", sep = "")
    n <- 21 - nchar(startText) + 1
    startText <- paste0(startText, paste(rep(" ", n), collapse = ""))
    cat(startText, values, "\n", sep = "")
  }

  if (isCorrMatrix) {
    rho_eb <- averageCorr(object, "EB")
    cat("  Avg corr (EB):     ", format(rho_eb, digits = 3), "\n")

    rhoSq_eb <- averageCorrSq(object, "EB")
    cat("  Avg corrSq (EB):   ", format(rhoSq_eb, digits = 3), "\n")

    # peff_eb <- sumInverseCorr(object, "EB")
    # digits <- ceiling(log10(object$p) + 2)
    # cat("  Eff # feat. (EB):  ", format(peff_eb, digits = digits), "\n")
  }

  cat(
    "  logLik:            ", format(object$logLik, digits = 1, scientific = FALSE),
    "\n"
  )
})


#' Get full covariance/correlation matrix from \link{eclairs}
#'
#' Get full covariance/correlation matrix from \link{eclairs} decomposition
#'
#' @param ecl eclairs decomposition
#' @param lambda shrinkage parameter for the convex combination.
#'
#' @param ... other arguments
#'
#' @return p x p covariance/correlation matrix
#'
#' @details The full matrix is computationally expensive to compute and uses a lot of memory for large p.  So it is better to use \link{decorrelate} or \link{mult_eclairs} to perform projections in \eqn{O(np)} time.
#'
#' @examples
#' library(Rfast)
#'
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y)
#'
#' # extract covariance implied by eclairs decomposition
#' getCov(ecl)[1:3, 1:3]
#'
#' @rdname getCov
#' @export
setGeneric("getCov", function(ecl, lambda, ...) standardGeneric("getCov"))

#' @rdname getCov
#' @export
setGeneric("getCor", function(ecl, lambda, ...) standardGeneric("getCor"))



#' @rdname getCov
#' @export
setMethod("getCov", c(ecl = "eclairs"), function(
    ecl, lambda,
    ...) {
  # if disableWarn is specified and is TRUE, set disableWarn = TRUE else
  # FALSE
  lst <- list(...)
  disableWarn <- FALSE
  if (!is.null(lst[["disableWarn"]])) {
    disableWarn <- lst[["disableWarn"]]
  }

  if (!disableWarn & all(ecl$sigma == 1)) {
    warning("eclairs estimated correlation, so the correlation matrix is returned")
  }

  .getCovCor(ecl,
    lambda = lambda,
    method = "covariance"
  )
})


#' @importFrom stats cov2cor
#' @rdname getCov
#' @export
setMethod("getCor", c(ecl = "eclairs"), function(
    ecl, lambda,
    ...) {
  .getCovCor(ecl,
    lambda = lambda,
    method = "correlation"
  )
})


.getCovCor <- function(ecl, lambda, method = c("covariance", "correlation")) {
  method <- match.arg(method)

  if (missing(lambda)) {
    lambda <- ecl$lambda
  }

  # reconstruct correlation matrix
  # C <- ecl$U %*% ((ecl$dSq * (1 - lambda)) * t(ecl$U)) +
    # diag(ecl$nu * lambda, ecl$p)

  A <- with(ecl, dmult(U[,seq(k),drop=FALSE], sqrt(dSq[seq(k)]), "right"))
  C <- (1 - lambda)*tcrossprod(A)
  diag(C) <- diag(C) + ecl$nu * lambda

  if (method == "covariance" & !all(ecl$sigma == 1)) {
    # scale by standard deviations
    D <- diag(ecl$sigma)
    C <- D %*% C %*% D
  }

  rownames(C) <- ecl$colnames
  colnames(C) <- ecl$colnames

  C
}



#' Estimate covariance/correlation with low rank and shrinkage
#'
#' Estimate the covariance/correlation between columns as the weighted sum of a low rank matrix and a scaled identity matrix.  The weight acts to shrink the sample correlation matrix towards the identity matrix or the sample covariance matrix towards a scaled identity matrix with constant variance.  An estimate of this form is useful because it is fast, and enables fast operations downstream.  The method is based on the Gaussian Inverse Wishart Empirical Bayes (GIW-EB) model.
#'
#' @param X data matrix with n samples as rows and p features as columns
#' @param k the rank of the low rank component
#' @param lambda shrinkage parameter. If not specified, it is estimated from the data.
#' @param compute evaluate either the \code{"covariance"} or \code{"correlation"} of \code{X}
#' @param n.samples number of samples data is from.  Usually \code{nrow(X)}, but can be other values in special cases.
#'
#' @return \link{eclairs} object storing:
#' \describe{
#'  \item{U: }{orthonormal matrix with k columns representing the low rank component}
#'  \item{dSq: }{eigen-values so that \eqn{U diag(d^2) U^T} is the low rank component}
#'  \item{lambda: }{shrinkage parameter \eqn{\lambda} for the scaled diagonal component}
#'  \item{sigma: }{standard deviations of input columns}
#'  \item{nu: }{diagonal value, \eqn{\nu}, of target matrix in shrinkage}
#'  \item{n: }{number of samples (i.e. rows) in the original data}
#'  \item{p: }{number of features (i.e. columns) in the original data}
#'  \item{k: }{rank of low rank component}
#'  \item{rownames: }{sample names from the original matrix}
#'  \item{colnames: }{features names from the original matrix}
#'  \item{method: }{method used for decomposition}
#'  \item{call: }{the function call}
#' }
#'
#' @details
#' Compute \eqn{U}, \eqn{d^2} to approximate the correlation matrix between columns of data matrix X by \eqn{U diag(d^2 (1-\lambda)) U^T + I\nu\lambda}.  When computing the covariance matrix, scale by the standard deviation of each feature.
#'
#' @examples
#' library(Rfast)
#'
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1, seed = 1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition: covariance
#' ecl <- eclairs(Y, compute = "covariance")
#'
#' ecl
#'
#' # eclairs decomposition: correlation
#' ecl <- eclairs(Y, compute = "correlation")
#'
#' ecl
#'
#' @importFrom Rfast standardise colVars eachrow
#' @importFrom irlba irlba
#' @importFrom methods new
#' @importFrom stats cor
#'
#' @export
eclairs <- function(X, k, lambda = NULL, compute = c("covariance", "correlation"), n.samples = nrow(X)) {
  stopifnot(is.matrix(X))
  compute <- match.arg(compute)

  # check value of lambda
  if (!missing(lambda) & !is.null(lambda)) {
    if (lambda < 0 || lambda > 1) {
      stop("lambda must be in (0,1)")
    }
  }

  n <- nrow(X)
  p <- ncol(X)

  # save row and columns names since X is overwritten
  rn <- rownames(X)
  cn <- colnames(X)

  # scale so that cross-product gives correlation features are *columns*
  X <- .standardise(X) / sqrt(n - 1)

  # save standard deviation of input features
  sigma <- attr(X, "sd")

  # if computing the correlation matrix,
  # set the scale of each feature to 1
  if (compute == "correlation") {
    sigma[] <- 1
  }

  if (missing(k)) {
    k <- min(n, p)
  } else {
    # k cannot exceed n or p
    k <- min(c(k, p, n))
  }

  # SVD of X to get low rank estimate of Sigma
  if (k < min(p, n) / 3) {
    dcmp <- irlba(X, k)
  } else {
    dcmp <- tryCatch({
      svd(X)
      },
      error = function(e){
      # very rarely, svd() above can fail
      # fall back on interative 
      k <<- min(c(k, dim(X)-1))
      suppressWarnings(irlba(X, k ))
      })

    # if k < min(n,p) truncate spectrum
    if (k < length(dcmp$d)) {
      dcmp$u <- dcmp$u[, seq_len(k), drop = FALSE]
      dcmp$v <- dcmp$v[, seq_len(k), drop = FALSE]
      dcmp$d <- dcmp$d[seq_len(k)]
    }
  }

  # Modify sign of dcmp$v and dcmp$u so principal components are consistant
  # This is motivated by whitening:::makePositivDiagonal() but here adjust
  # both U and V so reconstructed data is correct
  values <- sign0(diag(dcmp$v))

  # faster version
  dcmp$v <- eachrow(dcmp$v, values, "*")
  dcmp$u <- eachrow(dcmp$u, values, "*")

  ecl <- list(
    U = dcmp$v,
    dSq = dcmp$d^2,
    V = dcmp$u,
    lambda = NA,
    logLik = NA,
    sigma = sigma,
    nu = NA,
    n = n.samples,
    p = p,
    k = k,
    rownames = rn,
    colnames = cn,
    # compute is each column in X is increasing or decreasing
    direction = sign(c(cor(X, seq(nrow(X))))),
    method = "svd",
    call = match.call()
  )
  ecl <- new("eclairs", ecl)

  # estimate lambda and nu values
  res <- getShrinkageParams( ecl, k, lambda = lambda)
  ecl$lambda <- res$lambda
  ecl$nu <- res$nu
  ecl$logLik <- res$logLik

  ecl
}



# Like standard sign function, except sign(x) giving 0 is reset to give 1
sign0 <- function(x) {
  # use standard sign function
  res <- sign(x)

  # get entries that equal 0 and set them to 1
  i <- which(res == 0)
  if (length(i) > 0) {
    res[i] <- 1
  }

  res
}

# like Rfast::standardise()
# but returns SD of each column
.standardise <- function(x) {
  y <- t(x) - Rfast::colmeans(x)
  s <- sqrt(Rfast::rowsums(y^2) / (nrow(x) - 1))
  y <- y / s
  y <- t(y)

  attr(y, "sd") <- s
  y
}


#' Plot eclairs object
#'
#' Plot eclairs object
#' @param x eclairs object
#' @param y extra argument, not used
#' @param ... additional arguments
#'
#' @importFrom graphics points legend arrows par
#' @export
setMethod("plot", "eclairs", function(x, y, ...) {
  args <- match.call(expand.dots = TRUE)

  col <- args$col
  if (!is(col, "character")) {
    col <- "deepskyblue"
  }

  main <- args$main
  if (!is(main, "character")) {
    main <- "Eigen-values"
  }

  # plot observed eigen-values
  plot(x$dSq, ylab = "Eigen-values", main = main, ylim = c(0, max(x$dSq)))

  # plot shrunk eigen-values
  ev_shrunk <- with(x, (1 - lambda) * dSq + nu * lambda)
  points(ev_shrunk, col = col, pch = "+")

  # add legend
  legend("topright", legend = c("Observed", "Shrinkage estimate"), fill = c(
    "black",
    col
  ), border = "white", bty = "n")

  # get max x and y values
  maxy <- max(x$dSq)
  maxx <- length(ev_shrunk)

  # add info about dataset
  legend(0.7 * maxx, 0.8 * maxy, ncol = 2, legend = c(
    expression(lambda ~ ":"),
    expression(nu ~ ":"), expression(n ~ ":"), expression(p ~ ":"), format(x$lambda,
      digits = 3
    ), format(x$nu, digits = 3), format(x$n, big.mark = ","), format(x$p,
      big.mark = ","
    )
  ), bty = "n", text.width = 0)

  # if SVD is low rank,
  isLowRank <- length(x$dSq) < x$p

  if (isLowRank) {
    # start at last eigen-value, and end at right side of plot
    xvals <- c(maxx + 0.2, par("usr")[2])
    yvals <- rep(x$lambda * x$nu, 2)

    arrows(xvals[1], yvals[1], xvals[2], yvals[2], length = 0.1, col = col)
  }
})


