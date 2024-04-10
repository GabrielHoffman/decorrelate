

test_lm_proj = function(){

	# standard lm()
	fit = lm(Petal.Width ~ Sepal.Length + Sepal.Width, iris)
	res1 = coef(summary(fit))['Sepal.Width',]

	# Hypothesis test of Sepal.Width using pre-fit model 
	obj = decorrelate:::lm.projection(iris$Petal.Width, model.matrix(~Sepal.Length, iris))
	res2 = decorrelate:::lm.test(obj, iris$Sepal.Width)

	checkEqualsNumeric(res1, res2)
}