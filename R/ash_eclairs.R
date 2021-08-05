





# library(ashr)
# library(emulator)


# 1) Joint fine-mapping
# 	run ash on decorrelated data

# 2) Association



# #' Adaptive shrinkage in decorrelated space
# #'
# #' Adaptive shrinkage in decorrelated space
# #'
# #' @param betahat estimated regression coefficients
# #' @param sebetahat estimated standard error of regression coefficients
# #' @param Sigma.eclairs covariance/correlation matrix as an \link{eclairs} object
# #'
# #' @importFrom Matrix Diagonal
# #' @importFrom ashr ash
# #' @export

# ash_eclairs = function(betahat, sebetahat, Sigma.eclairs, ...){

# 	# decorrelate betahat
# 	betahat.decor = decorrelate(betahat, Sigma.eclairs)

# 	# decorrelate sebetahat
# 	# Currently quadratic time, but maybe possible in linear time?
# 	# A = getWhiteningMatrix( Sigma.eclairs )
# 	# step1 = decorrelate( Diagonal(length(sebetahat), sebetahat^2), Sigma.eclairs)
# 	# sebetahat.decor = sqrt(colSums(step1 * A))

# 	# same as:
# 	A = getWhiteningMatrix( Sigma.eclairs )
# 	sebetahat.decor = sqrt(quad.diag(diag(sebetahat^2), A))

# 	# k = 30 # number of components
# 	# a = seq(0, max(betahat)*1.1, length.out=k)
# 	# g.init = unimix(rep(1,k)/n, a, -a)
# 	# learn prior from decorrelated data
# 	ash.fit = ash(betahat.decor, sebetahat.decor, alpha=1)#), mixcompdist="normal")

# 	# pi1
# 	pi1 = 1 - ash.fit$fitted_g$pi[1]

# 	# sample beat values from posterior
# 	beta_post = get_post_sample(ash.fit, 1000)

# 	# back-transform into original space
# 	beta_post.back = mult_eclairs(beta_post, 
# 		U1 		= Sigma.eclairs$U, 
# 		dSq1 	= Sigma.eclairs$dSq, 
# 		lambda 	= Sigma.eclairs$lambda, 
# 		nu 		= Sigma.eclairs$nu, 
# 		alpha 	= 1/2)

# 	# local false sign rate from 	
# 	NegativeProb = apply(beta_post.back, 2, function(x){
# 		sum(x < 0) / length(x)		
# 		})

# 	PositiveProb = apply(beta_post.back, 2, function(x){
# 		sum(x > 0) / length(x)
# 		})

# 	ZeroProb = apply(beta_post.back, 2, function(x){
# 		sum(x == 0) / length(x)
# 		})

# 	df = data.frame( betahat 	= betahat, 
# 				    sebetahat 	= sebetahat,
# 				    NegativeProb= NegativeProb,
# 				    PositiveProb= PositiveProb, 
# 				    lfsr 		= compute_lfsr(NegativeProb, ZeroProb),
# 					PosteriorMean = colMeans(beta_post.back),
# 					PosteriorSD = apply(beta_post.back, 2, sd))

# 	list( joint = ash.fit, assoc = df)
# }







# 	# # For association, need back transform!!!!!!
# 	# # so use normal model
# 	# # need to back-transform full posterior distribution
# 	# # back transform
# 	# betahat.back = mult_eclairs(matrix(ash.fit$result$PosteriorMean, ncol=1), 
# 	# 	U1 		= Sigma.eclairs$U, 
# 	# 	dSq1 	= Sigma.eclairs$dSq, 
# 	# 	lambda 	= Sigma.eclairs$lambda, 
# 	# 	nu 		= Sigma.eclairs$nu, 
# 	# 	alpha 	= 1/2, 
# 	# 	transpose = TRUE)

# 	# beta_post = get_post_sample(ash.fit, 1000)


# 	# plot(colMeans(beta_post))
# 	# tstats = apply(betahat.back, 2, function(x) mean(x) / sd(x))
# 	# plot(colMeans(betahat.back))
# 	# plot(tstats)


# 	# Joint fine-mapping
# 	# Learn posterior distribution of original data 
# 	# using priors used from decorrelated data
# 	ash.fit2 = ash( betahat, sebetahat, g = ash.fit$fitted_g, fixg=TRUE, alpha=1)
# 	ash.fit2





# library(ashr)
# library(emulator)
# library(decorrelate)
# library(Matrix)
# library(Rfast)
# library(limma)

# set.seed(1)
# n = 800 # number of samples
# p = 8*200 # number of features

# # Create correlation matrix with autocorrelation
# autocorr.mat <- function(p = 100, rho = 0.9) {
# mat <- diag(p)
# return(rho^abs(row(mat)-col(mat)))
# }

# # create correlation matrix
# Sigma = autocorr.mat(p/8, .9)
# Sigma = bdiag(Sigma, Sigma)
# Sigma = bdiag(Sigma, Sigma)
# Sigma = bdiag(Sigma, Sigma)

# # draw data from correlation matrix Sigma
# Y = Rfast::rmvnorm(n, rep(0, p), sigma=Sigma*5.1)
# rownames(Y) = paste0("sample_", 1:n)
# colnames(Y) = paste0("gene_", 1:p)

# # eclairs decomposition


# x = Y[,1] + rnorm(n)*6


# fit = lmFit(t(Y), model.matrix(~x))
# fit = eBayes(fit)
# tab = topTable(fit, coef='x', sort.by='none', number=Inf)


# betahat = tab$t
# sebetahat = rep(1,p)



# par(mfrow=c(3,1))
# fit1 = ash(tab$t, rep(1,p))
# plot(get_lfsr(fit1), main="ash")

# Sigma.eclairs = eclairs(Y, compute="corr")
# fit2 = ash_eclairs(tab$t, rep(1,p), Sigma.eclairs)
# plot(get_lfsr(fit2$joint), main="ash_eclairs")

# plot(fit2$assoc$lfsr, main="ash_eclairs")






# head(fit2$result)




# W = getWhiteningMatrix(Sigma.eclairs)
# betahat.decor = c(tcrossprod(betahat, W))

# # same as:
# sebetahat.decor = sqrt(quad.diag(diag(sebetahat^2), W))

# fit1 = ash(betahat, sebetahat, pointmass=FALSE, alpha=1)

# fit2 = ash( betahat, sebetahat, g = fit1$fitted_g, fixg=TRUE, alpha=1)

# head(fit2$result)


# head(ash(betahat, sebetahat, pointmass=FALSE, alpha=1)$result)


# par(mfrow=c(2,1))
# fit1 = ash(tab$t, rep(1,p))
# plot(get_lfsr(fit1), main="ash")


# Sigma.eclairs = eclairs(Y, compute="corr", fastLambda=FALSE)
# fit2 = ash_eclairs(tab$t, rep(1,p), Sigma.eclairs)
# plot(get_lfsr(fit2), main="ash_eclairs")









# beta = c(rep(0,100),rnorm(100))
# sebetahat = abs(rnorm(200,0,1))
# betahat = rnorm(200,beta,sebetahat)

# beta.ash = ash(betahat, sebetahat)



# fit = ash_eclairs(betahat, sebetahat, Sigma.eclairs)



# dcmp = eigen(getCov(Sigma.eclairs))
# A = sweep(dcmp$vectors, 2, dcmp$values^-0.5, FUN='*')

# se.transform = as.numeric(sqrt(colSums((A * sebetahat^2) * A)))
# se.transform[1:5]

# a1 = sqrt(diag(A %*% diag(sebetahat^2) %*% t(A)))
# head(a1)


# colSums((A * sebetahat^2))[1:3]



# b = sqrt(quad.diag(diag(sebetahat^2), A))
# head(b)



# Sigma.eclairs = eclairs(Y, compute="corr", lambda=0)
# A = getWhiteningMatrix( Sigma.eclairs )

# a = sqrt(colSums(crossprod(diag(sebetahat^2), A) * A))
# head(a)

# step1 = decorrelate( Diagonal(length(sebetahat), sebetahat^2), Sigma.eclairs)
# a = sqrt(colSums(step1 * A))
# head(a)








# crossprod(diag(sebetahat^2), A)[1:3, 1:3]
# step1[1:3, 1:3]









