# Gabriel Hoffman
# April 16, 2024

#' Multiply by diagonal matrix
#' 
#' Multiply by diagonal matrix using efficient algorithm
#' 
#' @param M matrix 
#' @param v vector with entries forming a diagonal matrix matching the dimensions of \code{M} depending on the value of \code{side}
#' @param side is the matrix \code{M} \code{"left"} or \code{"right"} multiplied by the diagonal matrix 
#' 
#' @details 
#' Naively multiplying by a diagonal matrix with \code{p} entries takes \eqn{\mathcal{O}(p^2)}, even though must values in the diagonal matrix are zero.  R has built in syntax to scale \emph{columns}, so \code{diag(v) \%*\% M} can be evaluated with \code{v * M}.  This is often difficult to remember, hard to look up, and scaling \emph{rows} requires \code{t(t(M) * v)}.  This is hard to read and write, and requires two transpose operations.  
#'
#' Here, \code{dmult()} uses \code{Rcpp} to evaluate the right multiplication efficiently, and uses \code{v * M} for left multiplication.  This gives good performance and readability.
#' 
#' In principle, the \code{Rcpp} code can be modified to use BLAS for better performance by treating a \code{NumericMatrix} as a \code{C} array.  But this is not currently a bottleneck
#' 
#' @return matrix product
#' @examples
#' # right multiply
#' # mat %*% diag(v)
#' n = 30
#' p = 20
#' mat = matrix(n*p, n,p)
#' v = rnorm(p)
#' 
#' A = dmult(mat, v, side="right")
#' B = mat %*% diag(v)
#' range(A - B)
#' 
#' # left multiply
#' # diag(v) %*% mat
#' n = 30
#' p = 20
#' mat = matrix(n*p, p,n)
#' v = rnorm(p)
#' 
#' A = dmult(mat, v, side="left")
#' B = diag(v) %*% mat
#' range(A - B)
#
#' @import Rcpp
#' @export
dmult = function(M, v, side = c( "left", "right")){

	side = match.arg(side)

	stopifnot(is.matrix(M))
	stopifnot(is.vector(v))

	if( side == "right"){

		stopifnot(length(v) == ncol(M))

		# mat %*% diag(v)
		# t( t(M) * v )
		res = dmult_(M, v, dleft=FALSE)
	}else{

		stopifnot(length(v) == nrow(M))

		# diag(v) %*% mat
		res = v * M
	}

	res
}