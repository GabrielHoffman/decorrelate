
test_summary_statistics = function(){

	set.seed(1)
	n = 200 # number of samples
	p = 300 # number of features

	# create correlation matrix
	Sigma = matrix(.2, p, p)
	diag(Sigma) = 1

	# draw data from correlation matrix Sigma
	Y = rmvnorm(n, rep(0, p), sigma=Sigma)
	rownames(Y) = paste0("sample_", 1:n)
	colnames(Y) = paste0("gene_", 1:p)

	# eclairs decomposition
	Sigma.eclairs = eclairs(Y, compute="cor")

	# Average correlation value
	#--------------------------
	a = averageCorr( Sigma.eclairs, method="EB" )

	# check
	S = getCor(Sigma.eclairs)
	b = mean(S[lower.tri(S)])

	checkEqualsNumeric(a,b)


	# Average squared correlation value
	#--------------------------
	a = averageCorrSq( Sigma.eclairs, method="EB" )

	# check
	S = getCor(Sigma.eclairs)
	b = mean(S[lower.tri(S)]^2)

	checkEqualsNumeric(a,b)

	# Sum elements in inverse correlation matrix
	# Gives the effective number of independent features
	#--------------------------
	a = sumInverseCorr( Sigma.eclairs, method="EB" )

	b = sum(solve(getCor(Sigma.eclairs)))

	checkEqualsNumeric(a,b)


	# Effective variance
	#-------------------

	a = effVariance( Sigma.eclairs, method="EB" )

	b = det(getCor(Sigma.eclairs))^(1/p) 

	checkEqualsNumeric(a,b)



}