# July 24, 2025

#' Class fastcca
#'
#' Class \code{fastcca}
#'
#' @details Object storing:
#' \describe{
#'  \item{n.comp: }{number of canonical components}
#'  \item{cors: }{canonical correlations}
#'  \item{x.coefs: }{canonical coefficients for X}
#'  \item{x.vars: }{canonical variates for X}
#'  \item{y.coefs: }{canonical coefficients for Y}
#'  \item{y.vars: }{canonical variates for Y}
#'  \item{lambdas: }{shrinkage parameters from \code{eclairs}}
#' }
#' @name fastcca-class
#' @rdname fastcca-class
#' @exportClass fastcca
setClass("fastcca", contains = "list")


setMethod("print", "fastcca", function(x) {
  show(x)
})


setMethod("show", "fastcca", function(object) {
  x <- object

  cat("       Fast regularized canonical correlation analysis\n\n")

  k <- min(3, x$n.comp)

  cat("  Original data rows:", x$dims["n"], "\n")
  cat("  Original data cols: ", x$dims["p1"], ", ", x$dims["p2"], "\n", sep = "")
  cat("  Num components:    ", x$n.comp, "\n")
  cat("  Cor:               ", round(x$cor[seq_len(k)], digits = 3), "...\n")
  cat("  rho.mod:           ", round(x$rho.mod[seq_len(k)], digits = 3), "...\n")
  cat("  Cramer's V:        ", round(x$cramer.V, digits = 3), "\n")
  cat("  lambda:            ", format(x$lambdas, digits = 3), "\n")
})


#' Fast canonical correlation analysis
#'
#' Fast Canonical correlation analysis that is scalable to high dimensional data.  Uses covariance shrinkage and algorithmic speed ups to be linear time in p when p > n.
#'
#' @param X first matrix (n x p1)
#' @param Y first matrix (n x p2)
#' @param k number of canonical components to return
#' @param k.x number of singular vectors of X to use
#' @param k.y number of singular vectors of Y to use
#' @param lambda.x shrinkage parameter for X. If \code{NULL}, estimate from data
#' @param lambda.y shrinkage parameter for Y. If \code{NULL}, estimate from data
#' @param svd.method SVD algorithm string "svd", "irlba", or "pcaone" 
#'
#' @examples
#' pop <- LifeCycleSavings[, 2:3]
#' oec <- LifeCycleSavings[, -(2:3)]
#'
#' # fit CCA
#' fit <- fastcca(pop, oec)
#' fit
#' 
#' # predict second dataset given first
#' y.pred <- predict(fit, X = pop)
#
#' @importFrom Rfast cora colsums
#' @export
fastcca <- function(X, Y, k = NULL, k.x=min(dim(X)), k.y=min(dim(Y)), lambda.x=NULL, lambda.y=NULL, svd.method = c("svd", "irlba", "pcaone") ){

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
  Sig <- crossprod( dmult(ecl.x$V, gamma_x, "right"), 
                    dmult(ecl.y$V, gamma_y, "right"))

  if( is.null(k) ){
    k <- min(dim(Sig))
  }

  # if "svd" selected, but k indicates partial SVD
  # then use irlba
  if(k < min(dim(Sig)) / 3 & svd.method == "svd"){
    svd.method <- "irlba"
  }

  # SVD of cross-product
  dcmp <- run_svd(Sig, k, svd.method) 

  # Compute coefs
  # Note that coefs are not unique, 
  # but latent variables are???
  ax <- with(ecl.x, 1/sqrt(dSq * (1 - lambda) + lambda * nu))
  x.coefs <- dmult(ecl.x$U, ax, "right") %*% dcmp$u
  rownames(x.coefs) <- colnames(X)
  colnames(x.coefs) <- paste0("comp_", seq(k))

  ay <- with(ecl.y, 1/sqrt(dSq * (1 - lambda) + lambda * nu))
  y.coefs <- dmult(ecl.y$U, ay, "right") %*% dcmp$v
  rownames(y.coefs) <- colnames(Y)
  colnames(y.coefs) <- paste0("comp_", seq(k))

  # Compute latent variables
  x.vars <- X_scaled %*% x.coefs
  y.vars <- Y_scaled %*% y.coefs

  # rho <- diag(cor(x.vars, y.vars))[seq(k)]
  # Faster way to eval diag of correlation
  rho <- colsums(.standardise(x.vars) * .standardise(y.vars)) / (nrow(X)-1)
  names(rho) <- paste("comp", seq(k), sep = "")

  rho.mod <- dcmp$d

  res <- list( 
        dims = c(n=nrow(X), p1 = ncol(X), p2 = ncol(Y)),
        n.comp = k,
        rho.mod = rho.mod,
        cor = rho, 
        # Cramer's V-statistic for CCA
        cramer.V = sqrt(mean(rho.mod^2)),
        lambdas = c(x = ecl.x$lambda, y = ecl.y$lambda),
        x.coefs = x.coefs, 
        y.coefs = y.coefs, 
        x.vars = x.vars, 
        y.vars = y.vars,
        x.mu = attr(X_scaled, "mu"),
        y.mu = attr(Y_scaled, "mu"),
        x.sd = attr(X_scaled, "sd"),
        y.sd = attr(Y_scaled, "sd"))

  new("fastcca", res)
}

#' Predict one dataset from another
#' 
#' Using CCA fit, predict one dataset from another
#' 
#' @param object model fit from \code{fastcca()}
#' @param X first dataset 
#' @param Y second dataset 
#' @param ... other arguments, not used
#'
#' @details Specify which dataset to map from using \code{X} and \code{Y}.  
#'
#' @examples
#' pop <- LifeCycleSavings[, 2:3]
#' oec <- LifeCycleSavings[, -(2:3)]
#'
#' # fit CCA
#' fit <- fastcca(pop, oec)
#' 
#' # predict Y given X
#' y.pred <- predict(fit, X = pop)
#' head(y.pred)
#' 
#' # predict X given Y
#' x.pred <- predict(fit, Y = oec)
#' head(x.pred)
#
#' @export
setMethod("predict", "fastcca", function(object, X, Y,...) {

  if( !missing(X) & !missing(Y) ){
    stop("Only one of X and Y can be specified")
  }

  if( !missing(X) ){
    res <- predict_from_X(object, X)
  }

  if( !missing(Y) ){
    res <- predict_from_Y(object, Y)
  }

  res
})



#' @importFrom MASS ginv
#' @importFrom Rfast eachrow
predict_from_X <- function( fit, X ){

  # check that column names of X are the same as in 
  # the training data
  if( ! identical(rownames(fit$x.coefs), colnames(X)) ){
    stop("column names of X doesn't the training data ")
  }

  X <- as.matrix(X)

  # scale input data
  X_scaled <- eachrow(X, fit$x.mu, "-")
  X_scaled <- eachrow(X_scaled, fit$x.sd, "/")

  # project X onto latent space
  # project from latent space to Y
  y.pred <- X_scaled %*% (fit$x.coefs %*% ginv(fit$y.coefs))
  
  # rescale to Y
  y.pred <- eachrow(y.pred, fit$y.sd, "*")
  y.pred <- eachrow(y.pred, fit$y.mu, "+")

  colnames(y.pred) <- rownames(fit$y.coefs)
  rownames(y.pred) <- rownames(X)

  y.pred
}


#' @importFrom MASS ginv
#' @importFrom Rfast eachrow
predict_from_Y <- function( fit, Y ){

  # check that column names of Y are the same as in 
  # the training data
  if( ! identical(rownames(fit$y.coefs), colnames(Y)) ){
    stop("column names of Y doesn't the training data ")
  }

  Y <- as.matrix(Y)

  # scale input data
  Y_scaled <- eachrow(Y, fit$y.mu, "-")
  Y_scaled <- eachrow(Y_scaled, fit$y.sd, "/")

  # project Y onto latent space
  # project from latent space to X
  x.pred <- Y_scaled %*% (fit$y.coefs %*% ginv(fit$x.coefs))
  
  # rescale to Y
  x.pred <- eachrow(x.pred, fit$x.sd, "*")
  x.pred <- eachrow(x.pred, fit$x.mu, "+")

  colnames(x.pred) <- rownames(fit$x.coefs)
  rownames(x.pred) <- rownames(Y)

  x.pred
}






