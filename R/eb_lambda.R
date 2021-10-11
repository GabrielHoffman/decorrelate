# Gabriel Hoffman
# May 24, 2021
#
# Estimate shrinkage parameter by emprical Bayes
# Simple R adaptation of the Rcpp code of CRAN's beam package
# Current code only works when target matrix is identity







# cholD = chol(D);
# logdetD = 2*sum(log(cholD.diag()));
# logdetD = 2*sum(log(diag(chol(D))))

# # Optimal shrinkage
# deltaOpt = getDeltaOpt(n, p, eigs, logdetD);
# alphaOpt = deltaToAlpha(deltaOpt, n, p);   


# double alphaToDelta(double alpha, int n, int p){
#   return (alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha);
# }
alphaToDelta = function(alpha, n, p){
	(alpha*n+(1-alpha)*p+(1-alpha))/(1-alpha)
}

# double deltaToAlpha(double delta, int n, int p){
#   return (delta-p-1)/(n+delta-p-1);
# }
deltaToAlpha = function(delta, n, p){
	(delta-p-1)/(n+delta-p-1)
}


# double lpvarGamma(const double x, const int p) {
#   double ans = (p * (p - 1) * 0.25) * log(datum::pi);
#   for(int j = 1; j < (p + 1); ++j){
#     ans += std::lgamma(x - ((j - 1.0) * 0.5));
#   }
#   return ans;
# }
# CholWishart::lmvgamma

# lpvarGamma = function(x,p){
#   ans = (p * (p - 1) * 0.25) * log(pi);
#   for(j in 1:(p) ){
#     ans = ans + lgamma(x - ((j - 1.0) * 0.5));
#   }
#   ans
# }

# lpvarGamma(4, 2)
# CholWishart::lmvgamma(4, 2)



# double logML(const double delta, const int p, const int n, colvec eigs, const double logdetD){
#   double out = -0.5*n*p*std::log(datum::pi);
#   out += lpvarGamma((delta+n)*0.5, p);
#   out -= lpvarGamma(delta*0.5, p);
#   out += 0.5*delta*p*std::log(delta-p-1);
#   out -= 0.5*(delta+n)*sum(arma::log((delta-p-1)+eigs));
#   if(logdetD!=0){
#     out -= 0.5*n*logdetD;
#   }
#   return(out);
# }


# Integrated log-likelihood for empirical Bayes
# eigs stores non-zero eigen-values of p total eigen values
#' @importFrom CholWishart lmvgamma
logML = function(delta, p, n, eigs, logdetD){

	out = -0.5*n*p*log(pi);
	out = out + lmvgamma((delta+n)*0.5, p);
	out = out - lmvgamma(delta*0.5, p);
	out = out + 0.5*delta*p*log(delta-p-1);

	# in bream Rcpp code assums eigs are zero after rank n
	# eigs[(n+1):p] =0
	if( n > p){
		out = out - 0.5*(delta+n)*sum(log((delta-p-1)+eigs))
	}else{	
		# if full svd, length(eigs) is 0
		# but for low rank svd, use lenght(eigs)
		# out = out - 0.5*(delta+n)*(sum(log((delta-p-1)+eigs)))+ (p-n)*sum(log(delta-p-1))) 
		out = out - 0.5*(delta+n)*(sum(log((delta-p-1)+eigs)) + (p-length(eigs))*sum(log(delta-p-1))) 
	}

	out = out - 0.5*n*logdetD;
	
	out
}


# logML(6.010101, p, n, eigs, logdetD)



# double getDeltaOpt(const int n, const int p, colvec eigs, const double logdetD){
#   const double lowerVal = alphaToDelta(0.001, n, p);
#   const double upperVal = alphaToDelta(0.999, n, p);
#   const auto obj = [p, n, eigs, logdetD](double x) { return -logML(x, p, n, eigs, logdetD); };
#   boost::uintmax_t it = 1000;
#   const auto result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
#   //std::pair<double, double> result = brent_find_minima(obj, lowerVal, upperVal, 1000, it);
#   auto deltaOpt = 0.0, valOpt = 0.0;
#   std::tie(deltaOpt, valOpt) = result;
#   return(deltaOpt);
# }

getShrinkageValue = function(n, p, eigs, logdetD, minimum = 1e-4, lambda=NULL){

	# if lambda is NULL, estimate it
	if( is.null(lambda) ){
		# get upper and lower values of range
		lowerVal = alphaToDelta(minimum, n, p);
		upperVal = alphaToDelta(1-minimum, n, p);

		# get optimal estiamte of delta using log-likelihood
		res = optimize( logML, lower = lowerVal, upper=upperVal,
			p=p, n=n, eigs=eigs, logdetD=logdetD, 
			maximum=TRUE)

		# convert delta to alpha between 0 and 1
		lambda = deltaToAlpha( res$maximum, n, p)
		value = res$objective
	}else{
		# compute logML
		delta = alphaToDelta(lambda,n,p)
		value = logML(delta, p=p, n=n, eigs=eigs, logdetD=logdetD)
	}

	list(lambda = lambda, logML = value)
}




#' Estimate shrinkage parameter by empirical Bayes
#'
#' Estimate shrinkage parameter by empirical Bayes
#'
#' @param ev array of eigen values 
#' @param n number of samples
#' @param p number of features
#' @param nu scale of prior covariance matrix
#' @param lambda (default: NULL) If NULL, estimate lambda from data. Else evaluate logML using specified lambda value.
#'
# @seealso \link{getShrinkageValue}
#'
#' @details Estimate shrinkage parameter for covariance matrix estimation using empirical Bayes method (Hannart and Naveau, 2014; Leday and Richardson, 2019).  The shrinage estimate of the covariance matrix is \eqn{(1-\lambda)\hat\Sigma + \lambda \nu I}, where \eqn{\hat\Sigma} is the sample covariance matrix, given a value of \eqn{lambda}.  A large value of \eqn{\lambda} indicates more weight on the prior.
#'
#' @return value \eqn{\lambda} indicating the shrinkage between sample and prior covariance matrices.
#'
#' @examples
#' ev = c(10, 2, 1) # eigen values
#' n = 12 # samples
#' p = 3 # features
#' nu = 2 # scale of target covariance
#' 
#' estimate_lambda_eb(ev, n, p, nu)
#' 
#' @references{
#'   \insertRef{leday2019fast}{decorrelate}
#'
#'   \insertRef{hannart2014estimating}{decorrelate}
#' }
#'
#' @import Rdpack
#' @export
estimate_lambda_eb = function(ev, n, p, nu, lambda=NULL){

	# when D = diag(nu,p),
	# logdet(D) is 2*p*log(sqrt(nu))
	# and the eigen values of X %*% Dinv %*% t(X) equal
	# eigen(tcrossprod(X)) / nu
	# 	so divide the eigen-values by nu

	# D = diag(nu, p)
	# logdetD = 2*sum(log(diag(chol(D))))	

	# if is low rank
	if( length(ev) < min(n,p)){

		# if eigen-values are truncated, include additial eigen-values
		# so that sum(ev) equals the total variance, regardless of rank
		# Note that eclairs calls 
		# 	estimate_lambda_eb( n*dcmp$d^2, n, p, nu)
		# so eigen-values are already scaled by n
		totalVar = p * n * nu

		# if the emprical sum of variance is larger then the theoretical
		# then the specified n is likely mis-specified
		if( sum(ev) > totalVar ){
			warning("The sample size n is likely mis-specified")
		}

		# min(n,p) so works regardless of n > p or m < p
		# can give over-estimate of lambda for small k
		# 	but the deviation from full k can be huge
		# idx = seq(length(ev)+1, min(n,p))
		# ev[idx] = (totalVar-sum(ev)) / length(idx)

		# can give under-estimate of lambda for small k
		# but much smaller difference than above
		len = min(n,p) - length(ev) 
		remainingVar = totalVar-sum(ev)
		ev = append(ev, series_start_total(min(ev), remainingVar, len))
	}

	# estimate optimal lambda (i.e. alpha) value
	getShrinkageValue(n, p, ev / nu, logdetD=2*p*log(sqrt(nu)), minimum = 1e-4, lambda=lambda)
}




#' Create decreasing series of values
#' 
#' Create series of \eqn{n} values decreasing from \eqn{start} with a constant ratio between terms where the values sum to \eqn{totalSum}
#' 
#' @param start maximum value
#' @param totalSum specify target sum of values produces
#' @param n number of terms
#' 
# @examples
# start = 10
# totalSum = 40
# n = 4
# 
# values = series_start_total(start, totalSum, n)
# 
series_start_total = function(start, totalSum, n){

	# get alpha so that:
	# sum(start*alpha^seq(1,n)) == totalSum
	
	f = function(alpha, start, n, totalSum){

		(sum(start*alpha^seq(1,n)) - totalSum)^2
	}

	fit = optimize(f, c(1e-12, 1-1e-12), start=start, n=n, totalSum=totalSum, tol=1e-10)
	alpha = fit$minimum

	start*alpha^seq(1,n)
}




# #' Create decreasing series of values
# #' 
# #' Create series of \eqn{n} values decreasing from \eqn{start} to \eqn{end} with a constant ratio between terms.
# #' 
# #' @param start maximum value
# #' @param end minimum value
# #' @param n number of terms
# #' 
# #' @examples
# #' start = 10
# #' end = 1
# #' n = 43
# #' 
# #' values = series(start, end, n)
# #' 
# series_start_end = function(start, end, n){

# 	# get alpha so that start*alpha^n = end
# 	alpha = (end/ start)^(1/n)

# 	# create series
# 	sapply( seq(0, n), function(i){
# 		start*alpha^i
# 		})
# }







