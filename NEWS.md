# decorrelate 0.0.8

	July 3-12, 2021
	-  getCov() and getCor() now have lambda argument
	-  plot() for eclairs shows arrow on right for zero eigen-values
	-  estimate_lambda_eb() now returns logML for estimated or specified lambda
		- this is stored by eclairs()

# decorrelate 0.0.7
	June 1, 2021
	-  add whiten() that combines eclairs() and decorrelate() into one function call
	-  add eclairs_corMat() to perform decomposition on correlation matrix
	-  extend reform_decomp2() to work with result of eclairs_corMat() 

# decorrelate 0.0.6
	May 25, 2021
	-  add estimate_lambda_eb() to perform empirical Bayes estimation of lambda
		- works for truncated SVD by setting missing ev's to the mean residual variance
		- this makes estimation of lambda stable for varying ranks
	-  add plot() for eclairs

# decorrelate 0.0.5
	-  add reform_decomp()
	-  improve documentation

# decorrelate 0.0.4
	-  add lm_eclairs() and lm_each_eclairs()
	-  add R/lm_projection.R

# decorrelate 0.0.1
	-  initial version
	
