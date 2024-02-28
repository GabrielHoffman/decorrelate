# Gabriel Hoffman April 25, 2022 eclairs_sq: exclairs to estimate 2*cor(Y^2)





#' Compute eclairs decomp of squared correlation matrix
#'
#' Given the eclairs decomp of C, compute the eclairs decomp of C^2
#'
#' @param ecl estimate of correlation matrix from \code{eclairs()} storing \eqn{U}, \eqn{d_1^2}, \eqn{\lambda} and \eqn{\nu}

#' @param rank1 use the first 'rank' singular vectors from the SVD.  Using increasing rank1 will increase the accuracy of the estimation.  But note that the computationaly complexity is O(P choose(rank, 2)), where P is the number of features in the dataset
#' @param rank2 rank of \code{eclairs()} decomposition returned
#' @param varianceFraction fraction of variance to retain after truncating trailing eigen values
#'
#' @details Consider a data matrix X_{N x P} of P features and N samples where N << P. Let the columns of X be scaled so that C_{P x P} = XX^T.  C is often too big to compute directly since it is O(P^2) and O(P^3) to invert.  But we can compute the SVD of X in O(PN^2).
#' The goal is to compute the SVD of the matrix C^2, given only the SVD of C in less than O(P^2 time).  Here we compute this SVD of C^2 in O(PN^4) time, which is tractible for small N.
#' Moreover, if we use an SVD X = UDV^T with of rank R, we can approximate the SVD of C^2 in O(PR^4) using only D and V
#' In practice, this can be reduced to O(P (choose(R,2) + R)^2)
#'
#' @examples
#' # Compute correlations directly and using eclairs decomp
#'
#' n <- 600 # number of samples
#' p <- 100 # number of features
#'
#' # create correlation matrix
#' Sigma <- autocorr.mat(p, .9)
#'
#' # draw data from correlation matrix Sigma
#' Y <- Rfast::rmvnorm(n, rep(0, p), sigma = Sigma, seed = 1)
#' rownames(Y) <- paste0("sample_", seq(n))
#' colnames(Y) <- paste0("gene_", seq(p))
#'
#' # correlation computed directly
#' C <- cor(Y)
#'
#' # correlation from eclairs decomposition
#' ecl <- eclairs(Y, compute = "cor")
#' C.eclairs <- getCor(ecl, lambda = 0)
#'
#' all.equal(C, C.eclairs)
#'
#' # Correlation of Y^2
#' #-------------------
#'
#' # exact quadratic way
#' C <- 2 * cor(Y)^2
#'
#' # faster low rank
#' ecl2 <- eclairs_sq(ecl)
#' C.eclairs <- 2 * getCor(ecl2, lambda = 0)
#'
#' all.equal(C.eclairs, C)
#'
#' @return compute the eclairs of C^2
#' @importFrom Rfast eachrow colVars
#' @importFrom irlba irlba
#' @export
eclairs_sq <- function(
    ecl, rank1 = ecl$k, rank2 = Inf, varianceFraction = 1) {
  if (!is(ecl, "eclairs")) {
    stop("ecl must be of class eclairs")
  }

  d <- sqrt(ecl$dSq)
  V <- ecl$U

  rank1 <- min(rank1, ncol(V))

  if (rank1 < ncol(V)) {
    V <- V[, seq_len(rank1), drop = FALSE]
    d <- d[seq_len(rank1), drop = FALSE]
  }

  # compute scaled eigen vectors mat = sweep(V, 2, d, FUN='*')
  mat <- eachrow(V, d, "*")
  rm(V)

  # compute matrix of element-wise products for all pairs of columns
  n <- ncol(mat)
  j1 <- rep.int(seq(1, n), seq(n, 1))
  j2 <- sequence(seq(n, 1)) - 1L + j1
  G <- mat[, j1, drop = FALSE] * mat[, j2, drop = FALSE]
  rm(mat)

  # scale if j1 != j2\t\t
  s <- (j1 == j2)
  s[s] <- 1
  s[!s] <- sqrt(2)
  # G = sweep(G, 2, c(s)*sqrt(2), FUN='*')\t G = eachrow(G, c(s)*sqrt(2),
  # '*')\t
  G <- eachrow(G, c(s), "*")

  # retain columns with highest variance

  # due to precision isues, make sure variance >= 0
  df <- data.frame(index = seq_len(ncol(G)), variance = pmax(colVars(G), 0))

  df <- df[order(df$variance, decreasing = TRUE), ]
  df$cumsum <- cumsum(df$variance / sum(df$variance))

  # set cutoff to the specified varianceFraction if no columns pass this
  # cutoff, use new cutoff
  varianceFraction <- max(varianceFraction, df$cumsum[min(nrow(df), 10)])

  idx <- which(df$cumsum <= varianceFraction)
  G <- G[, df$index[idx], drop = FALSE]

  # SVD of G

  n <- ecl$n
  p <- ecl$p
  nu <- 1

  # SVD of X to get low rank estimate of Sigma
  if (rank2 < min(p, n) / 3) {
    # if( is.null(warmStart) ){ dcmp = svds(G, k, isreal=TRUE)
    dcmp <- irlba(G, rank2) # should be faster thatn PRIMME::svds
    # }else{\t\t\t \tdcmp = svds(G, rank2,
    # u0=ecl$U[,seq_len(rank2)], isreal=TRUE) }
  } else {
    dcmp <- svd(G)

    # if rank2 < min(n,p) truncate spectrum
    if (rank2 < length(dcmp$d)) {
      dcmp$u <- dcmp$u[, seq_len(rank2), drop = FALSE]
      dcmp$v <- dcmp$v[, seq_len(rank2), drop = FALSE]
      dcmp$d <- dcmp$d[seq_len(rank2)]
    }
  }

  # eclairs

  # Estimate lambda by empirical Bayes, using nu as scale of target Since
  # data is scaled to have var 1 (instead of n), multiply by n
  res <- estimate_lambda_eb(n * dcmp$d^2, n, p, nu)

  # Modify sign of dcmp$v and dcmp$u so principal components are consistant
  # This is motivated by whitening:::makePositivDiagonal() but here adjust
  # both U and V so reconstructed data is correct
  values <- sign(diag(dcmp$u))

  # faster version
  dcmp$v <- eachrow(dcmp$v, values, "*")
  dcmp$u <- eachrow(dcmp$u, values, "*")

  ecl <- list(
    U = dcmp$u, dSq = dcmp$d^2, V = dcmp$v, lambda = res$lambda, logML = res$logML,
    nu = nu, n = n, p = p, k = length(dcmp$d), rownames = ecl$rownames,
    colnames = ecl$colnames, method = "svd", call = match.call()
  )

  new("eclairs", ecl)
}
