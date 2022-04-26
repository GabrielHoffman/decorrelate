
library(RUnit)



test_eclairs_sq = function(){

	# Compute correlations directly and using eclairs decomp
	set.seed(1)
	n = 600 # number of samples
	p = 100 # number of features

	# create correlation matrix
	Sigma = autocorr.mat(p, .9)

	# draw data from correlation matrix Sigma
	Y = Rfast::rmvnorm(n, rep(0, p), sigma=Sigma)
	rownames(Y) = paste0("sample_", 1:n)
	colnames(Y) = paste0("gene_", 1:p)

	# correlation computed directly
	C = cor(Y)

	# correlation from eclairs decomposition
	ecl = eclairs(Y, compute="cor")
	C.eclairs = getCor(ecl, lambda=0)

	checkEqualsNumeric(C, C.eclairs)

	# Correlation of Y^2
	#-------------------

	# exact quadratic way
	C = 2*cor(Y)^2 

	# faster low rank
	ecl2 = eclairs_sq(ecl)
	C.eclairs = 2*getCor(ecl2, lambda=0)

	checkEqualsNumeric(C.eclairs, C)
}


# source("/Users/gabrielhoffman/workspace/repos/decorrelate/R/eclairs_sq.R") 

# # exact quadratic way
# C = 2*cor(Y)^2 

# # faster low rank
# ecl2 = eclairs_sq(ecl, rank2=Inf)
# C.eclairs = 2*getCor(ecl2, lambda=0)

# plot(C, C.eclairs)
# abline(0,1, col='red')

