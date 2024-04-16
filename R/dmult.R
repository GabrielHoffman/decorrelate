
#' @export
dmult = function(M, v, side = c( "left", "right")){

	side = match.arg(side)

	stopifnot(is.matrix(M))
	stopifnot(is.vector(v))

	if( side == "right"){

		stopifnot(length(v) == ncol(M))

		# t( t(M) * v )
		res = dmult_(M, v, dleft=FALSE)
	}else{

		stopifnot(length(v) == nrow(M))

		res = v * M
	}

	res
}