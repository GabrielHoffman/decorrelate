// #include  "Rcpp.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dmult_( const NumericMatrix m, const NumericVector v, bool dleft){

  if( dleft && m.nrow() != v.size() )
    stop("Non-conformable arrays") ;

  if( ! dleft && m.ncol() != v.size() )
    stop("Non-conformable arrays") ;

  // initialize new matrix
  NumericMatrix out(m.nrow(), m.ncol()) ;

  // diag(v) %*% mat
  if( dleft ){
    for (int i = 0; i < m.nrow(); i++) {
      for (int j = 0; j < m.ncol(); j++) {
        out(i,j) = m(i,j) * v[i];
      }
    }
  }

  // mat %*% diag(v)
  if( ! dleft ){
    for (int j = 0; j < m.ncol(); j++) {
      for (int i = 0; i < m.nrow(); i++) {
        out(i,j) = m(i,j) * v[j];
      }
    }
  }

  return out ;
}
  


// Use RcppArmadaillo, but not faster
// [[Rcpp::export]]
arma::mat dmult_arma(arma::mat M, arma::rowvec v, bool dleft=true ){ 
  if( dleft )
    M.each_col() %= v.t();
  else
    M.each_row() %= v;

  return M;
}









// // Use Rcpp sugar to vectorize
// // [[Rcpp::export]]
// NumericMatrix dmult2_( const NumericMatrix m, const NumericVector v, bool dleft=true ){

//   if( dleft && m.nrow() != v.size() )
//     stop("Non-conformable arrays") ;

//   if( ! dleft && m.ncol() != v.size() )
//     stop("Non-conformable arrays") ;

//   // initialize new matrix
//   NumericMatrix out(m.nrow(), m.ncol()) ;

//   if( dleft ){
//     for (int i = 0; i < m.nrow(); i++) {
//         out(i,_) = m(i,_) * v[i];
//     }
//   }
//   if( ! dleft ){
//     for (int j = 0; j < m.ncol(); j++) {
//       out(_,j) = m(_,j) * v[j];
//     }
//   }
//   return out ;
// }


  

  