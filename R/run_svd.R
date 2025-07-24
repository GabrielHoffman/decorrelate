
#' Compute SVD with selected algorithm
#' 
#' Compute SVD using exact, approximate or randomized algorithms
#' 
#' @param X data matrix
#' @param k number of left and right singular values to compute
#' @param method SVD algorithm string "svd", "irlba", or "pcaone" 
#' @param ... other arguments passed to SVD functions
#' 
#' @details Compute SVD using select algorithm.  Let data matrix \code{X} have \code{n} rows and \code{p} columns.  At most \code{min(n,p)} singular vectors can be computed.
#' \itemize{
#'  \item "svd": exact algorithm using standard \code{svd()} function.  Fastest choice if computing greater than about \code{min(n,p)/3} singular vectors.
#' 
#'  \item "irlba": approximate algorithm using implicitly restarted Lanczos bidiagonalization in \code{irlba::irlba()}.  Much faster than "svd" for computing less than about \code{min(n,p)/3}  singular vectors, while matching \code{svd()} to moderate precision
#' 
#'  \item "pcaone": randomized SVD algorithm implemented in \code{pcaone::pcaone()}.  Fastest method for large datasets, when computing about \code{min(n,p)/3}  singular vectors.  Singular vectors and values are approximate and can differ from \code{svd()} by a few percent.  
#' }
#' @examples
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
#' X <- hilbert(9)[, 1:6]
#' dcmp1 <- run_svd(X, method = "svd")
#' dcmp2 <- run_svd(X, method = "irlba")
#' dcmp3 <- run_svd(X, method = "pcaone")
#' 
#' @importFrom irlba irlba 
#' @importFrom pcaone pcaone
#' @importFrom Rfast eachrow
#' @export
run_svd = function(X, k = min(dim(X)), method = c("svd", "irlba", "pcaone"),...){

  method <- match.arg(method)
  stopifnot("k must be numeric" = is.numeric(k))
  stopifnot("k must be positive" = k > 0)

  if( k > min(dim(X))-1){
    method = "svd"
  }

  if( method == "svd" ){
    dcmp <- tryCatch({svd(X, k, k,...)},
      error = function(e) {
        # very rarely, svd() above can fail
        # fall back on interative
        k <<- min(c(k, dim(X) - 1))
        suppressWarnings(irlba(X, k))
      }
    ) 
  }else if( method == "irlba"){    
    dcmp <- irlba(X, k, k,...) 
  }else if( method == "pcaone"){    
    dcmp <- pcaone(X, k,...) 
  }

  # Ensure diagonals of v are positive
  # so PCs are always pointed in the same direction
  values <- sign0(diag(dcmp$v))
  dcmp$v <- eachrow(dcmp$v, values, "*")
  dcmp$u <- eachrow(dcmp$u, values, "*")

  dcmp 
}
