# decorrelate 0.1.5
* April 2, 2025
* update docs
* update `sumInverseCorr()` 

# decorrelate 0.1.4
* March 20, 2025
* update docs
* update `sumInverseCorr()` based on `direction`

# decorrelate 0.1.3
* March 11, 2025
* in `eclairs()`, if `svd()` fails fall back on `irlba()`
* `sumInverseCorr()` has upper bound of `p`
* fix docs

# decorrelate 0.1.2
* April 23, 2024
* use `dmult()` instead of transposing
* estimate `nu` to give correlation matrix close to having diagonals 1

# decorrelate 0.1.1
* Feb 27, 2024
* update all functions to process covariance by shrinking correlation and retaining variance information

# decorrelate 0.1.0
* Feb 27, 2024
* update all functions to process covariance by shrinking correlation and retaining variance information

# decorrelate 0.0.19
* Jan 02, 2024
* add `mahalanobisDistance()`

# decorrelate 0.0.18
* Sept 18, 2023
* add `averageCorr()`, `averageCorrSq()` and `sumInverseCorr()`
 * add tests for these functions
* add `alpha` parameter to `quadForm()`

# decorrelate 0.0.17
* add `n.samples` argument to \code{eclairs()}


# decorrelate 0.0.16
* add `eclairs_sq`

# decorrelate 0.0.15
* use `irlba` for SVD instead of `PRIMME`
* improve documentation

# decorrelate 0.0.14
* add `series_start_total()` and use it in `estimate_lambda_eb()` for partial SVD
 * this substantially improves estimation of lambda when partial SVD is used
* add `averageCorr()`

# decorrelate 0.0.13
* add redundancy index as `x.ri`, `y.ri`
* add effective variance and effective dependency metrics

# decorrelate 0.0.12
* `fastcca()` and `cca()` give equivalent results
 * but canonical variates are rotated with respect to each other

# decorrelate 0.0.11
* Scale Cxy by `sqrt(1-lambda.x)*sqrt(1-lambda.y)`
* Add Cramer's V statistic

# decorrelate 0.0.10
* Sept 13, 2021
* add ``fastcca()` and `cca()`

# decorrelate 0.0.9
* Sept 7, 2021
* add `kappa()` to compute condition number
* add `logDet()` to compute log determinant
* add `cca()` for canonical correlation analysis
* simplify examples

# decorrelate 0.0.8
* July 3-12, 2021
* `getCov()` and `getCor()` now have lambda argument
* `plot()` for eclairs shows arrow on right for zero eigen-values
* `estimate_lambda_eb()` now returns logML for estimated or specified lambda
 - this is stored by `eclairs()`

# decorrelate 0.0.7
* June 1, 2021
* add `whiten()` that combines `eclairs()` and `decorrelate()` into one function call
* add `eclairs_corMat()` to perform decomposition on correlation matrix
* extend `reform_decomp2()` to work with result of `eclairs_corMat()`

# decorrelate 0.0.6
* May 25, 2021
* add `estimate_lambda_eb()` to perform empirical Bayes estimation of lambda
 - works for truncated SVD by setting missing ev's to the mean residual variance
- this makes estimation of lambda stable for varying ranks
* add `plot()` for eclairs

# decorrelate 0.0.5
-  add `reform_decomp()`
-  improve documentation

# decorrelate 0.0.4
-  add `lm_eclairs()` and `lm_each_eclairs()`
-  add R/lm_projection.R

# decorrelate 0.0.1
-  initial version
	
