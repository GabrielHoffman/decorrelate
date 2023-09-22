#' Compute condition number
#'
#' Compute condition number of matrix from \code{eclairs} decomposition
#'
#' @param z \code{eclairs()} decomposition
#' @param lambda specify lambda to override value from \code{z}
#'
#' @examples
#' library(Rfast)
#' set.seed(1)
#' n <- 800 # number of samples
#' p <- 200 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # eclairs decomposition
#' ecl <- eclairs(Y, compute = "correlation")
#'
#' # compute condition number
#' kappa(ecl)
#'
#' @rdname kappa
#' @export
setMethod("kappa", c(z = "eclairs"), function(z, lambda = NULL) {
  if (!missing(lambda) & !is.null(lambda)) {
    if (lambda < 0 || lambda > 1) {
      stop("lambda must be in (0,1)")
    }
  } else {
    lambda <- z$lambda
  }

  # compute max eigen-values
  l_max <- (1 - lambda) * z$dSq[1] + lambda * z$nu

  # compute min eigen-values
  #--------------------------

  # if k >= 0 use last computed eigen-value else sample cov is low rank so
  # last sample eigen-value is zero
  e_min <- ifelse(z$k >= z$p, z$dSq[length(z$dSq)], 0)

  l_min <- (1 - lambda) * e_min + lambda * z$nu

  # condition number is ratio of largest and smallest eigen-value
  l_max / l_min
})
