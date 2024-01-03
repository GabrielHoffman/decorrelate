

test_mahalanobisDistance = function(){
	library(RUnit)
	library(Rfast)
	set.seed(1)
	n <- 800 # number of samples
	p <- 200 # number of features

	# create correlation matrix
	Sigma <- autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y <- rmvnorm(n, rep(0, p), sigma = Sigma * 5.1)

	# eclairs decomposition
	ecl <- eclairs(Y)

	# no centering
	a = mahalanobis(Y, FALSE, cov = cov(Y))
	b = mahalanobisDistance(ecl, Y, lambda=0)
	checkEqualsNumeric(a,b)

	# centering
	Y_center = scale(Y, scale=FALSE)
	mu = colMeans(Y)
	a = mahalanobis(Y, mu, cov = cov(Y))
	b = mahalanobisDistance(ecl, Y_center, lambda=0)
	checkEqualsNumeric(a,b)

	# compare to quadForm that comptues full square
	a = diag(quadForm(ecl, Y_center))
	b = mahalanobisDistance(ecl, Y_center)
	checkEqualsNumeric(a,b)


	a = mahalanobis(Y_center[1:2,], FALSE, cov=cov(Y))
	b = mahalanobis(Y_center[1:3,], FALSE, cov=cov(Y))
	c = mahalanobisDistance(ecl, Y_center[1:2,], lambda=0)
	d = mahalanobisDistance(ecl, Y_center[1:3,], lambda=0)
	checkEqualsNumeric(a,b[1:2])
	checkEqualsNumeric(a,c[1:2])
	checkEqualsNumeric(a,d[1:2])
}