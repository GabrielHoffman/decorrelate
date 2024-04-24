# # Gabriel Hoffman April 12, 2021 Recompute eclairs after reforming original
# # data, and droppping features

#' Recompute eclairs after dropping features
#'
#' Recompute eclairs after dropping features
#'
#' @param ecl covariance/correlation matrix as an \link{eclairs} object
#' @param k the rank of the low rank component
#' @param drop array of variable names to drop.
#'
#' @details Reform the dataset from the eclairs decomposition, drop features, then recompute the eclairs decomposition.  If the original SVD/eigen was truncated, then the reconstruction of the original data will be approximate.  Note that the target shrinkage matrix is the same as in \code{ecl}, so \eqn{\nu} is not recomputed from the retained features.
#'
#' @return \link{eclairs} decomposition for a subset of features
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
#' # Correlation
#' #------------
#'
#' # eclairs decomposition
#' Sigma.eclairs <- eclairs(Y, compute = "correlation")
#'
#' # features to drop
#' drop <- paste0("gene_", 1:100)
#'
#' # Compute SVD on subset of eclairs decomposition
#' ecl1 <- reform_decomp(Sigma.eclairs, drop = drop)
#'
#' ecl1
#'
#' @importFrom utils head
#' @export
reform_decomp <- function(ecl, k = ecl$k, drop = NULL) {
  stopifnot(is(ecl, "eclairs"))

  if (is.null(drop) | length(drop) == 0) {
    warning("No variables are dropped.  Returning original eclairs")
    return(ecl)
  }

  if (!(ecl$method %in% c("svd", "eigen"))) {
    stop("method must be 'svd' or 'eigen'")
  }

  # check if there are entries in drop that are not in ecl
  problematic <- which(!(drop %in% ecl$colnames))

  # if yes, throw warning
  if (length(problematic) > 0) {
    warning(length(problematic), " entries in drop are not in ecl\nThe first few entries not found are: ",
      paste0("'", paste(head(drop[problematic], 3), collapse = "', '"), "'"),
      immediate. = TRUE
    )
  }

  # keep features that are not in drop
  keep <- which(!ecl$colnames %in% drop)

  V <- dSq <- U <- NULL # pass R CMD BiocCheck

  n <- ecl$n
  p <- length(keep)

  # k cannot exceed n or p
  k <- min(c(k, p, n))

  # reconstruct original dataset from SVD
  if (ecl$method == "svd") {
    X_reconstruct <- with(ecl, V %*% (sqrt(dSq) * t(U)))
    X_reconstruct <- .standardise(X_reconstruct) %*% diag(ecl$sigma)

    compute <- ifelse(all(ecl$sigma == 1), "correlation", "covariance")
    ecl <- eclairs(X_reconstruct[, keep, drop = FALSE],
      k = k,
      # lambda = ecl$lambda,
      compute = compute
    )
  } else {
    # if( ecl$method == 'eigen' ) Reconstruct original dataset
    # from eigen decomposition

    # C_reconstruct <- with(ecl, U %*% (dSq * t(U)))
    # C_reconstruct <- diag(ecl$sigma) %*% C_reconstruct %*% diag(ecl$sigma)
    dM <- with(ecl, dmult(U, dSq, "right"))
    Cr <- with(ecl, tcrossprod(U, dM))
    dM2 <- dmult(Cr, ecl$sigma, "right")
    C_reconstruct <- dmult(dM2, ecl$sigma, "left")

    ecl <- eclairs_corMat(C_reconstruct[keep, keep, drop = FALSE],
      n = n,
      k = k
    )
    # lambda = ecl$lambda
  }

  ecl$nu <- ecl$nu
  ecl$call <- match.call()

  ecl
}
