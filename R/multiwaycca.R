


#' Class multiwaycca
#'
#' Class \code{multiwaycca}
#'
#' @name multiwaycca-class
#' @rdname multiwaycca-class
#' @exportClass multiwaycca
setClass("multiwaycca", contains = "list")


setMethod("print", "multiwaycca", function(x) {
  show(x)
})


setMethod("show", "multiwaycca", function(object) {
  x <- object

  cat("       Fast regularized canonical correlation analysis\n\n")

  k <- min(3, x$n.comp)

  cat("  Original data rows:", x$n, "\n")
  cat("  Original data cols: ", paste0(x$p, collpse=" "), "\n", sep = "")
  cat("  Num components:    ", x$n.comp, "\n")
  # cat("  Cor:               ", round(x$cor[seq_len(k)], digits = 3), "...\n")
  # cat("  rho.mod:           ", round(x$rho.mod[seq_len(k)], digits = 3), "...\n")
  # cat("  Cramer's V:        ", round(x$cramer.V, digits = 3), "\n")
  cat("  lambda:            ", format(x$lambdas, digits = 3), "\n")
})



multiwaycca = function(datasets, 
  k = NULL, 
  ks = rep(NULL, length(datasets)),  
  lambdas = rep(NULL, length(datasets)), 
  svd.method = c("svd", "irlba", "pcaone")){

  stopifnot("k must be positive" = k > 0)
  stopifnot("ks must be positive" = all(ks > 0))
  stopifnot("lambdas must be in (0,1)" = all(lambdas >= 0) & (lambdas <=1))
  svd.method <- match.arg(svd.method)

  # put k and lambda as attributes for each dataset
  for(i in seq(length(datasets))) {
    attr(datasets[[i]], "k") <- ks[i]
    attr(datasets[[i]], "lambda") <- lambdas[i]
    attr(datasets[[i]], "svd.method") <- svd.method
  }

  if( is.null(names(datasets)) ){
    names(datasets) <- paste0("dataset_", seq(length(datasets)))
  }

  # convert to eclairs decomp
  datasets <- lapply(datasets, to_eclairs)

  # shrink singular values
  # whiten data
  data_whitened <- lapply(datasets, function(x){
    gamma <- with(x$ecl, sqrt(dSq/((1-lambda)*dSq + lambda*nu)))
    dmult(x$ecl$V, gamma, "right")
    })

  # Concatenate whitened data
  dataCat <- do.call(cbind, data_whitened)

  # SVD
  dcmp <- run_svd(dataCat, k, svd.method) 

  n = nrow(dataCat)

  # get coefficients
  betas <- lapply( datasets, function(x){ 

    a <- with(x$ecl, 1/sqrt(dSq * (1 - lambda) + lambda * nu)/ sqrt(n-1))

    B <- dmult(x$ecl$U, a, "right") %*% crossprod(x$ecl$V, dcmp$u)
    colnames(B) <- paste0("comp_", seq(ncol(B)))
    rownames(B) <- x$ecl$colnames
    B
  })

  # construct latent variables
  vars = lapply(seq(length(datasets)), function(i){
    if( !is.null(datasets[[i]]$X_scaled) ){      
      v <- datasets[[i]]$X_scaled %*% betas[[i]]
    }else{
      ecl = datasets[[i]]$ecl
      v <- dmult(ecl$V, sqrt(ecl$dSq), "right") %*% 
      crossprod(ecl$U, betas[[i]])
    }
    v
    })
  names(vars) = names(betas)

  # return results
  as.multiwaycca( datasets, betas, vars, dcmp$d, k, lambdas)
}

as.multiwaycca = function(datasets, betas, vars, d, k, lambdas){

  mu <- lapply(datasets, function(x){
    if( ! is.null(x$X_scaled) ){
      mu <- attr(x$X_scaled, "mu")
    }else{
      mu <- x$ecl$mu
      names(mu) <- x$ecl$colnames
    }
    mu
  })

  sds <- lapply(datasets, function(x){
    if( ! is.null(x$X_scaled) ){
      s <- attr(x$X_scaled, "sd")
    }else{
      s <- x$ecl$sigma
      names(s) <- x$ecl$colnames
    }
    s
  }) 

  obj <- list(datasets = datasets, 
            coef = betas, 
            vars = vars,
            d = d[seq(k)],
            n = nrow(datasets[[1]]$X_scaled),
            p = sapply(datasets, function(x) ncol(x$X_scaled)),
            n.comp = k,
            lambdas = lambdas,
            mu = mu,
            sd = sds)
  new("multiwaycca", obj)
}




#' @importFrom MASS ginv
#' @importFrom Rfast eachrow
#' @export
setMethod("predict", "multiwaycca", function(object, X, given, target,...) {

  if( !given %in% names(object$datasets) ){
    stop(cat("Not a valid datset name:", given), "\n")
  }
  if( !target %in% names(object$datasets) ){
    stop(cat("Not a valid datset name:", target), "\n")
  }

  if( ! is.matrix(X) ){
    X <- as.matrix(X)
  }

  if( ! identical(colnames(X), names(object$mu[[given]])) ){
    stop("Column names in X and training data from 'given' must match ")
  }

  # scale input data
  X_scaled <- eachrow(X, object$mu[[given]], "-")
  X_scaled <- eachrow(X_scaled, object$sd[[given]], "/")

  # project inout onto latent space
  # project from latent space to target space
  y.pred <- X_scaled %*% (object$coef[[given]] %*% ginv(object$coef[[target]] ))
  
  # rescale to Y
  y.pred <- eachrow(y.pred, object$sd[[target]], "*")
  y.pred <- eachrow(y.pred, object$mu[[target]], "+")

  colnames(y.pred) <- rownames(object$coef[[target]])
  rownames(y.pred) <- rownames(X)

  y.pred
})





to_eclairs = function(X){

  k.x <- attr(X, "k")
  lambda.x <- attr(X, "lambda")
  svd.method <- attr(X, "svd.method")

  X_scaled <- NULL

  # Convert X to eclairs decomposition
  if( is.matrix(X) | is.data.frame(X) ){

    if( is.null(k.x) ) k.x <- min(dim(X))

    # SVD and shrinkage
    X_scaled <- decorrelate:::.standardise(X)
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
  }else{
    stop("invalid datatype")
  }

  list(ecl = ecl.x, X_scaled = X_scaled)
}


