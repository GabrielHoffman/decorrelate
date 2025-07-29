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
#' @param X \code{matrix}, \code{data.frame}, \code{eclairs} decomposition or \code{list} from \code{svd()}
#' @param Y \code{matrix}, \code{data.frame}, \code{eclairs} decomposition or \code{list} from \code{svd()}
#' @param k number of canonical components to return
#' @param k.x number of singular vectors of X to use
#' @param k.y number of singular vectors of Y to use
#' @param lambda.x shrinkage parameter for X. If \code{NULL}, estimate from data
#' @param lambda.y shrinkage parameter for Y. If \code{NULL}, estimate from data
#' @param svd.method SVD algorithm string "svd", "irlba", or "pcaone" 
#'
#' @details
#' Objects in \code{X} and \code{Y} are converted to \code{eclairs} decomposition from \code{matrix}, \code{data.frame} or \code{list} from \code{svd()}.  The \code{predict()} function maps from one input to another, except when input object is \code{list} from \code{svd()}.
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
#'
#' # Advanced
#' # fastcca() also works with eclairs() objects
#' ecl.x <- eclairs(as.matrix(pop))
#' ecl.y <- eclairs(as.matrix(oec))
#' fit2 = fastcca(ecl.x, ecl.y)
#' 
#' # fastcca() also works mixed types
#' fit3 = fastcca(pop, ecl.y)
#' 
#' # it also works with results of svd()
#
#' @importFrom Rfast cora colsums
#' @export
fastcca <- function( X, Y, 
  k = NULL, 
  k.x = NULL, 
  k.y = NULL, 
  lambda.x = NULL, 
  lambda.y = NULL, 
  svd.method = c("svd", "irlba", "pcaone") ){

  # checks
  stopifnot("k must be positive" = k > 0)
  stopifnot("k.x must be positive" = k.x > 0)
  stopifnot("k.y must be positive" = k.y > 0)
  stopifnot("k.y must be positive" = k.y > 0)
  stopifnot("lambda.x must be in (0,1)" = (lambda.x >= 0) & (lambda.x <=1))
  stopifnot("lambda.y must be in (0,1)" = (lambda.y >= 0) & (lambda.y <=1))
  svd.method <- match.arg(svd.method)

  # if (nrow(X) != nrow(Y)) {
  #   stop("X and Y must have the same number of rows")
  # }

  # Convert X to eclairs decomposition
  if( is.matrix(X) | is.data.frame(X) ){

    if( is.null(k.x) ) k.x <- min(dim(X))

    # SVD and shrinkage
    X_scaled <- .standardise(X)
    ecl.x <- eclairs(X_scaled, 
              k = k.x, 
              lambda = lambda.x, 
              svd.method = svd.method)

  }else if( is(X, "eclairs")){
    ecl.x <- X
  }else if( is(X, "list") ){
    dcmp <- X;
    ecl.x <- as.eclairs(dcmp, 
                k = ncol(dcmp$u), 
                n = nrow(dcmp$u), 
                p = ncol(dcmp$v),
                svd.method = "svd",
                lambda = lambda.x)
  }
 
  # Convert Y to eclairs decomposition
  if( is.matrix(Y) | is.data.frame(Y) ){

    if( is.null(k.y) ) k.y <- min(dim(Y))

    # SVD and shrinkage
    Y_scaled <- .standardise(Y)
    ecl.y <- eclairs(Y_scaled, 
              k = k.y, 
              lambda = lambda.y, 
              svd.method = svd.method)

  }else if( is(Y, "eclairs")){
    ecl.y <- Y
  }else if( is(Y, "list") ){
    dcmp <- Y;
    ecl.y <- as.eclairs(dcmp, 
                k = ncol(dcmp$u), 
                n = nrow(dcmp$u), 
                p = ncol(dcmp$v),
                svd.method = "svd",
                lambda = lambda.y)
  }

  # Perform CCA
  res = .fastcca( ecl.x, ecl.y, k = k, svd.method = svd.method)

  # Assign variable names
  rownames(res$y.coefs) <- colnames(Y)

  # Compute latent variables
  if( exists("X_scaled") ){
    rownames(res$x.coefs) <- colnames(X_scaled)
    res$x.vars <- X_scaled %*% res$x.coefs
    res$x.mu = attr(X_scaled, "mu")
    res$x.sd = attr(X_scaled, "sd")
  }else{
    rownames(res$x.coefs) <- ecl.x$colnames

    # X_recon <- with(ecl.x, V %*% diag(sqrt(dSq)) %*% t(U) )
    # res$x.vars <- X_recon %*% res$x.coef
    res$x.vars <- dmult(ecl.x$V, sqrt(ecl.x$dSq), "right") %*% 
      crossprod(ecl.x$U, res$x.coef)

    res$x.mu = ecl.x$mu
    res$x.sd = ecl.x$sigma
  }
  if( exists("Y_scaled") ){
    rownames(res$y.coefs) <- colnames(Y_scaled)
    res$y.vars <- Y_scaled %*% res$y.coefs
    res$y.mu = attr(Y_scaled, "mu")
    res$y.sd = attr(Y_scaled, "sd")
  }else{
    rownames(res$y.coefs) <- ecl.y$colnames

    # Y_recon <- with(ecl.y, V %*% diag(sqrt(dSq)) %*% t(U) )
    # res$y.vars <- Y_recon %*% res$y.coef
    res$y.vars <- dmult(ecl.y$V, sqrt(ecl.y$dSq), "right") %*% 
      crossprod(ecl.y$U, res$y.coef)

    res$y.mu = ecl.y$mu
    res$y.sd = ecl.y$sigma
  }

  # Faster way to eval diag of correlation
  # diag(cor(x.vars, y.vars))[seq(k)]
  res$cor <- colsums(.standardise(res$x.vars) * .standardise(res$y.vars)) / (nrow(res$y.vars)-1)
  names(res$cor) <- paste("comp", seq(k), sep = "")

  res
}





# CCA on eclairs objects
.fastcca <- function(ecl.x, ecl.y, k = NULL, svd.method = c("svd", "irlba", "pcaone") ){

  # checks
  stopifnot("k must be positive" = k > 0)
  svd.method <- match.arg(svd.method)

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
  colnames(x.coefs) <- paste0("comp_", seq(k))

  ay <- with(ecl.y, 1/sqrt(dSq * (1 - lambda) + lambda * nu))
  y.coefs <- dmult(ecl.y$U, ay, "right") %*% dcmp$v
  colnames(y.coefs) <- paste0("comp_", seq(k))

  rho.mod <- dcmp$d

  res <- list( 
        dims = c(n=ecl.x$n, p1 = ecl.x$p, p2 = ecl.y$p),
        n.comp = k,
        rho.mod = rho.mod,
        # cor = rho, 
        # Cramer's V-statistic for CCA
        cramer.V = sqrt(mean(rho.mod^2)),
        lambdas = c(x = ecl.x$lambda, y = ecl.y$lambda),
        x.coefs = x.coefs, 
        y.coefs = y.coefs)

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
    stop("column names of X doesn't match the training data ")
  }

  if( ! is.matrix(X) ){
    X <- as.matrix(X)
  }

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
    stop("column names of Y doesn't match the training data ")
  }  

  if( ! is.matrix(Y) ){
    Y <- as.matrix(Y)
  }

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






