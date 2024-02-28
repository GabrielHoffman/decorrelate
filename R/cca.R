# Sept 9, 2021

#' Canonical correlation analysis
#'
#' Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
#'
#' @param X first matrix (n x p1)
#' @param Y first matrix (n x p2)
#' @param k number of canonical components to return
#' @param lambda.x optional shrinkage parameter for estimating covariance of X. If NULL, estimate from data.
#' @param lambda.y optional shrinkage parameter for estimating covariance of Y. If NULL, estimate from data.
#'
#' @return statistics summarizing CCA
#' @details
#' Results from standard CCA are based on the SVD of \eqn{\Sigma_{xx}^{-\frac{1}{2}} \Sigma_{xy} \Sigma_{yy}^{-\frac{1}{2}}}.
#'
#' Avoids computation of \eqn{\Sigma_{xx}^{-\frac{1}{2}}} by using eclairs.  Avoids cov(X,Y) by framing this as a matrix product that can be distributed. Uses low rank SVD.
#' Other regularized CCA adds lambda to covariance like Ridge. Here it is a mixture
#'
#' @examples
#' pop <- LifeCycleSavings[, 2:3]
#' oec <- LifeCycleSavings[, -(2:3)]
#' fastcca(pop, oec)
#'
#' @importFrom irlba irlba
#' @export
cca <- function(X, Y, k = min(dim(X), dim(Y)), lambda.x = NULL, lambda.y = NULL) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  n1 <- nrow(X)
  n2 <- nrow(Y)
  p1 <- ncol(X)
  p2 <- ncol(Y)

  if (n1 != n2) {
    stop("X and Y must have same number of rows")
  }

  k <- min(dim(X), dim(Y), k)

  if (k < 0) {
    stop("k must be > 1")
  }

  n.comp <- min(ncol(X), ncol(Y), nrow(X), k)

  # mean center columns
  X.center <- scale(X, scale = FALSE)
  Y.center <- scale(Y, scale = FALSE)
  n <- nrow(X)

  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }

  # Evaluate: Sig_{xx}^{-.5} Sig_{xy} Sig_{yy}^{-.5} Sig = minvsqrt(cov(X))
  # %*% cov(X,Y) %*% minvsqrt(cov(Y)) But using faster algorithm

  # eclairs decomposition
  ecl.x <- eclairs(X, lambda = lambda.x, compute = "cov")
  ecl.y <- eclairs(Y, lambda = lambda.y, compute = "cov")

  # creates a matrix that is ncol(X) * ncol(Y) tmp1 = decorrelate(cov(X,Y),
  # ecl.x, transpose=TRUE)

  # evaluates the same quantity in linear time X and Y are mean centered by
  # cov()
  tmp1 <- crossprod(Y.center, decorrelate(X.center / (n - 1), ecl.x))

  Sig <- decorrelate(t(tmp1), ecl.y)

  # scale by shrinkage parameters
  Sig <- Sig * sqrt(1 - ecl.x$lambda) * sqrt(1 - ecl.y$lambda)

  # SVD on
  if (k > min(dim(Sig)) / 3) {
    dcmp <- svd(Sig)

    if (k < min(dim(Sig))) {
      dcmp$d <- dcmp$d[seq_len(n.comp)]
      dcmp$u <- dcmp$u[, seq_len(n.comp), drop = FALSE]
      dcmp$v <- dcmp$v[, seq_len(n.comp), drop = FALSE]
    }
  } else {
    dcmp <- irlba(Sig, nv = n.comp)
  }

  # set diagonals to be positive
  dcmp$v <- eachrow(dcmp$v, sign(diag(dcmp$v)), "*")
  dcmp$u <- eachrow(dcmp$u, sign(diag(dcmp$u)), "*")

  x.coefs <- decorrelate(dcmp$u, ecl.x, transpose = TRUE)
  x.vars <- X %*% x.coefs
  rownames(x.vars) <- rownames(X)

  y.coefs <- decorrelate(dcmp$v, ecl.y, transpose = TRUE)
  y.vars <- Y %*% y.coefs
  rownames(y.vars) <- rownames(X)

  rho <- diag(cor(x.vars, y.vars))[seq(n.comp)]
  names(rho) <- paste("can.comp", seq(n.comp), sep = "")

  rho.mod <- dcmp$d
  names(rho.mod) <- paste("can.comp", seq_len(n.comp), sep = "")

  if (k < min(dim(X), dim(Y))) {
    idx <- seq(1, k)
  } else {
    idx <- seq(1, k - 1)
  }

  # Cramer's V-statistic for CCA
  cramer.V <- sqrt(mean(rho.mod[idx]^2))

  # based on CRAN::yacca compute loadings
  load.x <- cor(X, x.vars)
  load.y <- cor(Y, y.vars)

  # compute redundancy index
  ri.x <- rho.mod^2 * apply(load.x, 2, function(x) mean(x^2))
  ri.y <- rho.mod^2 * apply(load.y, 2, function(x) mean(x^2))

  names(ri.x) <- paste("can.comp", seq_len(n.comp), sep = "")
  names(ri.y) <- paste("can.comp", seq_len(n.comp), sep = "")

  res <- list(
    n.comp = n.comp, rho.mod = rho.mod, cor = rho, cramer.V = cramer.V,
    x.coefs = x.coefs, x.vars = x.vars, x.ri = ri.x, y.coefs = y.coefs, y.vars = y.vars,
    y.ri = ri.y, lambdas = c(x = ecl.x$lambda, y = ecl.y$lambda), dims = c(
      n1 = n1,
      n2 = n2, p1 = p1, p2 = p2
    )
  )

  new("fastcca", res)
}
