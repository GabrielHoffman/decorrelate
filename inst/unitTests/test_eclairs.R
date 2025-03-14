
test_eclairs_eclairs_corMat = function(){

	# library(decorrelate)
	# library(RUnit)
	library(Matrix)
	library(Rfast)
	set.seed(1)
	n = 800 # number of samples
	p = 8*200 # number of features

	# Create correlation matrix with autocorrelation
	autocorr.mat <- function(p = 100, rho = 0.9) {
	mat <- diag(p)
	return(rho^abs(row(mat)-col(mat)))
	}

	# create correlation matrix
	Sigma = autocorr.mat(p/8, .9)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)
	Sigma = bdiag(Sigma, Sigma)

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
	rownames(Y) = paste0("sample_", 1:n)
	colnames(Y) = paste0("gene_", 1:p)
	Y = scale(Y)

	# eclairs decomposition
	ecl = eclairs(Y)
	Y.decorr = decorrelate(Y, ecl)

	# source("eclair_corMat.R")
	# from correlation matrix
	ecl2 = eclairs_corMat(cora(Y), n=n)
	Y.decorr2 = decorrelate(Y, ecl2)

	checkEqualsNumeric( Y.decorr, Y.decorr2, tol=1e-7)

	# source("eclairs_reform.R")
	# source("eclair_corMat.R")

	ecl2 = eclairs_corMat(cora(Y), n=n)
	# Compute SVD on subset of eclairs decomposition
	drop = paste0("gene_",1:1)

	# reform
	eclA = reform_decomp( ecl2, drop=drop)

	# eclairs using original data
	eclB = eclairs_corMat(cora(Y[,!colnames(Y) %in% drop]), n=n, lambda=ecl2$lambda)
	

	# since the last singular value is very small (below machine precision)
	# the large singular vectors of U and V are unstable 
	# and should be ommitted from the check
	idx = seq(1:50)
	res1 = checkEqualsNumeric( eclA$V[,idx], eclB$V[,idx])
	res2 = checkEqualsNumeric( eclA$U[,idx], eclB$U[,idx], tol=1e-5)
	res3 = checkEqualsNumeric( eclA$dSq, eclB$dSq)

	# Comparison of irlba::partial_eigen and PRIMME::eigs_sym
	# the leading eigen vectors have little relative error, 
	# 	and it increase further down the spectrum
	d = sapply(1:ncol(eclA$U), function(i){
		mean(abs(eclA$U[,i] -  eclB$U[,i]))
	})
	# plot(log10(d))


}


# system.time({
# a <- eigs_sym(Sigma, 1200, isreal=TRUE)
# })

# system.time({
# b <- eigs_sym(Sigma, 1200, isreal=TRUE, tol=1e-8)
# })

# a$stats$numMatvecs
# b$stats$numMatvecs

# source("eclair_corMat.R")
# C = cor(Y)

# system.time({
# eclA <- eclairs_corMat(C, n=n)
# })


# system.time({
# eclB <- eclairs_corMat(C, n=n, warmStart=eclA)
# })




# system.time({
# a <- reform_decomp(eclA)
# })


# system.time({
# b <- reform_decomp(eclA, drop=colnames(C)[1:10])
# })
















