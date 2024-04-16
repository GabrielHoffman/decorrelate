




test_dmult = function(){

	n = 30
	p = 20
	mat = Rfast::matrnorm(n,p)
	v = rnorm(p)

	# right multiply
	A = dmult(mat, v, side="right")
	B = mat %*% diag(v)
	checkEqualsNumeric(A,B)

	# left multiply
	n = 30
	p = 20
	mat = Rfast::matrnorm(p,n)
	v = rnorm(p)

	A = dmult(mat, v, side="left")
	B = diag(v) %*% mat
	checkEqualsNumeric(A,B)
}

# Timings

# library(Matrix)
# library(RhpcBLASctl)
# library(decorrelate)
# library(RUnit)

# RhpcBLASctl::omp_set_num_threads(1)

# n = 3000
# p = 6000
# mat = Rfast::matrnorm(n,p)
# v = rnorm(p)

# library(microbenchmark)
# m <- microbenchmark( 
# 	# mat %*% diag(v) , 
# 	as.matrix(mat %*% Diagonal(length(v), v)) , 
# 	dmult( mat , v, side="right" ) ,
# 	decorrelate:::dmult_arma( mat , v, FALSE ) ,
# 	sweep(mat, 2, v, FUN = "*") ,   
# 	t( t(mat) * v ) , 
# 	times = 20 )
# print( m , "relative" , order = "median" , digits = 3 )


# n = 3000
# p = 6000
# mat = Rfast::matrnorm(p,n)
# v = rnorm(p)

# library(microbenchmark)
# m <- microbenchmark( 
# 	diag(v) %*% mat , 
# 	as.matrix(Diagonal(length(v), v) %*% mat) , 
# 	dmult( mat , v, side="left" ) ,
# 	decorrelate:::dmult_arma( mat , v, TRUE ) ,
# 	sweep(mat, 1, v, FUN = "*") ,   
# 	v * mat , 
# 	times = 20 )
# print( m , "relative" , order = "median" , digits = 3 )







# mmult( mat , v, FALSE  ) 
# fun(mat, v, FALSE )




# q()
# R
# library(decorrelate)
# library(microbenchmark)


# func <- '
# NumericMatrix mmult( const NumericMatrix m, const NumericVector v, bool dleft = true ){

#   if( dleft && m.nrow() != v.size() )
#     stop("Non-conformable arrays") ;

#   if( ! dleft && m.ncol() != v.size() )
#     stop("Non-conformable arrays") ;

#   NumericMatrix out(m.nrow(), m.ncol()) ;

#   if( dleft ){
#     for (int i = 0; i < m.nrow(); i++) {
#       for (int j = 0; j < m.ncol(); j++) {
#         out(i,j) = m(i,j) * v[i];
#       }
#     }
#   }
#   if( ! dleft ){
#     for (int j = 0; j < m.ncol(); j++) {
#       for (int i = 0; i < m.nrow(); i++) {
#         out(i,j) = m(i,j) * v[j];
#       }
#     }
#   }
#   return out ;
# }
# '

# #  Make it available
# Rcpp::cppFunction( func, verbose=TRUE, rebuild=TRUE )

# mmult2 = function(M, v, dleft){
# 	mmult(M, v, dleft)
# }


# Rcpp:::cppFunction(
#     "arma::mat fun(arma::mat M, arma::rowvec v, bool dleft=true) 
#     { 
#   if( dleft )
#     M.each_col() %= v.t();
#   else
#     M.each_row() %= v;

#   return M;
#     }", depends = "RcppArmadillo"
# )



# n = 3
# p = 2
# mat = Rfast::matrnorm(n,p)
# v = rnorm(p)
# mmult( mat , v, FALSE  ) 
# fun(mat, v, FALSE )


# decorrelate:::dmult_arma(mat, v, FALSE )





# # A = dmult( mat , v, FALSE )
# A = mat %*% diag(v)
# B = mmult( mat , v, FALSE  ) 

# A[1:2, 1:2]
# B[1:2, 1:2]

# range(A - B)



# n = 300
# p = 20000
# mat = Rfast::matrnorm(n,p)
# v = rnorm(p)

# m <- microbenchmark( 
# 	# mat %*% diag(v) , 
# 	decorrelate:::dmult_arma(mat, v, FALSE),
# 	dmult( mat , v, FALSE ) ,
# 	fun( mat , v, FALSE ) ,
# 	decorrelate:::dmult_( mat , v, FALSE ) ,
# 	decorrelate:::dmult2_( mat , v, FALSE ) ,
# 	 mmult( mat , v,FALSE ), 
# 	 mmult2( mat , v,FALSE ), 
# 	sweep(mat, 2, v, FUN = "*") ,   
# 	t( t(mat) * v ) , 
# 	times = 20 )
# print( m , "relative" , order = "median" , digits = 3 )



# n = 3000
# p = 2000
# mat = Rfast::matrnorm(p,n)
# v = rnorm(p)

# m <- microbenchmark( 
# 	# diag(v) %*% mat , 
# 	decorrelate:::dmult_( mat , v, TRUE ) ,
# 	decorrelate:::dmult2_( mat , v, TRUE ) ,
# 	# decorrelate:::dmult_arma(mat, v, TRUE),
# 	fun( mat , v, TRUE ) ,
# 	 mmult( mat , v,TRUE ), 
# 	sweep(mat, 1, v, FUN = "*") ,   
# 	v * mat , 
# 	times = 10L )
# print( m , "relative" , order = "median" , digits = 3 )


# library(RUnit)

# # Right multiply
# n = 3000
# p = 20
# mat = Rfast::matrnorm(n,p)
# v = rnorm(p)

# mat.in = mat

# A = decorrelate:::dmult( mat , v, FALSE )
# B = mat %*% diag(v) 

# checkEqualsNumeric(mat.in, mat)
# checkEqualsNumeric(A,B)


# # left multiply
# n = 3000
# p = 20
# mat = Rfast::matrnorm(p,n)
# v = rnorm(p)

# mat.in = mat

# A = decorrelate:::dmult( mat , v, TRUE )
# B = diag(v) %*% mat

# checkEqualsNumeric(mat.in, mat)
# checkEqualsNumeric(A,B)


















