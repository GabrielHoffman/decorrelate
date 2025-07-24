

#' Fast canonical correlation analysis
#'
#' Fast Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
#'
#' @param X first matrix (n x p1)
#' @param Y first matrix (n x p2)
#' @param k number of canonical components to return
#' @param lambda.x optional shrinkage parameter for estimating covariance of X. If NULL, estimate from data.
#' @param lambda.y optional shrinkage parameter for estimating covariance of Y. If NULL, estimate from data.
#'
#' @examples
#' pop <- LifeCycleSavings[, 2:3]
#' oec <- LifeCycleSavings[, -(2:3)]
#'
#' fastcca2(pop, oec)
#'
#' @export
fastcca2 = function(X, Y, k = NULL, k.x=min(dim(X)), k.y=min(dim(Y)), lambda.x=NULL, lambda.y=NULL, svd.method = c("svd", "irlba", "pcaone") ){

  # checks
  stopifnot("k must be positive" = k > 0)
  stopifnot("k.x must be positive" = k.x > 0)
  stopifnot("k.y must be positive" = k.y > 0)
  stopifnot("k.y must be positive" = k.y > 0)
  stopifnot("lambda.x must be in (0,1)" = (lambda.x >= 0) & (lambda.x <=1))
  stopifnot("lambda.y must be in (0,1)" = (lambda.y >= 0) & (lambda.y <=1))
  svd.method <- match.arg(svd.method)

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }

  # scale columns
  X_scaled <- .standardise(X)
  Y_scaled <- .standardise(Y)

  # SVD and shrinkage
  ecl.x <- eclairs(X_scaled, 
    k = k.x, 
    lambda = lambda.x, 
    svd.method = svd.method)

  ecl.y <- eclairs(Y_scaled, 
    k = k.y, 
    lambda = lambda.y, 
    svd.method = svd.method)

  # shrink singular values
  gamma_x <- with(ecl.x, sqrt(dSq/((1-lambda)*dSq + lambda*nu)))
  gamma_y <- with(ecl.y, sqrt(dSq/((1-lambda)*dSq + lambda*nu)))

  # cross-product
  Sig <- crossprod(dmult(ecl.x$V, gamma_x, "right"), 
                  dmult(ecl.y$V, gamma_y, "right"))

  if( is.null(k) ){
    k <- min(dim(Sig))
  }

  # if "svd" selected, but k indicates partial SVD
  # then use irlba
  if(k < min(dim(Sig)) / 3 & svd.method == "svd"){
    svd.method = "irlba"
  }

  # SVD of cross-product
  dcmp = run_svd(Sig, k, svd.method) 

  # Compute coefs
  # Note that coefs are not unique, 
  # but latent variables are???
  ax <- with(ecl.x, 1/sqrt(dSq * (1 - lambda) + lambda * nu))
  x.coefs <- dmult(ecl.x$U, ax, "right") %*% dcmp$u

  ay <- with(ecl.y, 1/sqrt(dSq * (1 - lambda) + lambda * nu))
  y.coefs <- dmult(ecl.y$U, ay, "right") %*% dcmp$v

  # Compute latent variables
  x.vars <- X_scaled %*% x.coefs
  rownames(x.vars) <- rownames(X)

  y.vars <- Y_scaled %*% y.coefs
  rownames(y.vars) <- rownames(Y)

  rho.mod = dcmp$d

  res = list( 
        dims = c(n=nrow(X), p1 = ncol(X), p2 = ncol(Y)),
        n.comp = k,
        # Cramer's V-statistic for CCA
        cramer.V = sqrt(mean(rho.mod^2)),
        lambdas = c(x = ecl.x$lambda, y = ecl.y$lambda),
        x.coefs = x.coefs, 
        y.coefs = y.coefs, 
        rho.mod = rho.mod,
        x.vars = x.vars, 
        y.vars = y.vars)

  new("fastcca", res)

}