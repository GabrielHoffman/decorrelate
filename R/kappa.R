

#' Compute condition number
#' 
#' Compute condition number of matrix from \code{eclairs} decomposition
#' 
#' @param z \code{eclairs} decomposition
#' @param lambda specify lambda to override value from \code{z}
#' 
#' @rdname kappa
#' @export
setMethod('kappa', c(z = "eclairs"),
	function(z, lambda=NULL){

	if (!missing(lambda) & !is.null(lambda)) {
        if (lambda < 0 || lambda > 1) {
            stop("lambda must be in (0,1)")
        }
    }else{
    	lambda = z$lambda
    }

	# compute max eigen-values
	l_max = (1-lambda) * z$dSq[1] + lambda * z$nu

	# compute min eigen-values
	#--------------------------

	# if k >= 0 use last computed eigen-value
	# else sample cov is low rank so last sample eigen-value is zero
    e_min = ifelse( z$k >= z$p, z$dSq[length(z$dSq)], 0)

	l_min = (1-lambda) * e_min + lambda * z$nu

	# condition number is ratio of largest and smallest eigen-value
	l_max / l_min
})



# ecl = eclairs(Y)
# C = getCov(ecl)
# # eigen(C)$values
# kappa(C, exact=TRUE)

# kappa(ecl)

