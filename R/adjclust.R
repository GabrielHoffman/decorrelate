
library(inline)
library(sparseMatrixStats)




matL <- rcpp(signature(Csq_ = "dsCMatrix", h_ = "integer"),
               " 
               arma::SpMat<double> Csq = as< arma::SpMat<double> >(Csq_);
              
              int p = Csq.n_rows;
               int h = Rcpp::as<int>(h_);
               arma::SpMat<double> out(p, h+1);

                int k;
               double value;

               omp_set_num_threads(4) ;

               #pragma omp parallel for               
               for( int i=0; i<p; i++){
                k = 0;
                for( int j=i; j<std::min(i+h+1,p); j++){
                  if(k == 0){
                    out(i,k) = 1;
                  }else{
                    value = Csq(i,j);
                    if( value!= 0.0){
                        out(i,k) = 2.0*value;
                    }
                  }
                  k++;
                }
               }

               return wrap(out);
               ", plugin = "RcppArmadillo", settings=settings)


# Rcpp::Rcout << i-h-1 << ' ' << i << ' ' << p-i-1 << ' ' << j << ' '  << k << ' ' <<std::endl;
matR <- rcpp(signature(Csq_ = "dsCMatrix", h_ = "integer"),
               " 
               arma::SpMat<double> Csq = as< arma::SpMat<double> >(Csq_);
              
              int p = Csq.n_rows;
              int h = Rcpp::as<int>(h_);
              arma::SpMat<double> out(p, h+1);

              int k;
               double value;
              omp_set_num_threads(4) ;

               #pragma omp parallel for   
               for( int i=0; i<p; i++){
                k = 0;
                for( int j=i; j>=std::max(i-h, (int) 0); j--){
                  if(k == 0){
                    out(p-i-1,k) = 1;
                  }else{
                    value = Csq(i,j);
                    if( value != 0.0){
                        out(p-i-1,k) = 2.0*value;
                    }
                  }
                  k++;
                }
               }

               return wrap(out);
               ", plugin = "RcppArmadillo")

run.adjclust <- function(mat, type = c("similarity", "dissimilarity"), h) {
  # sanity checks
  type <- match.arg(type)
  if (!(nrow(mat) == ncol(mat)))
    stop("Input matrix is not a square matrix")
  if (any(is.na(mat)))
    stop("Missing values in the input are not allowed")
  
  p <- nrow(mat)
  
  # 'h'
  if (!is.numeric(h))
    stop("Input band width 'h' must be numeric")
  if (h != as.integer(h))
    stop("Input band width 'h' must be an integer")
  if (h < 0)
    stop("Input band width 'h' must be non negative")
  if (h >= p) 
    stop("Input band width 'h' must be strictly less than dimensions of matrix")
  
  # data preprocessing
  if (type == "dissimilarity") {
    mat <- 1 - 0.5*(mat^2)
  }
  
  # res_cc <- checkCondition(mat)
  # if (is.numeric(res_cc)) {
  #   message(paste("Note: modifying similarity to ensure positive heights...
  #     added", res_cc, "to diagonal (merges will not be affected)"))
  #   mat <- mat + diag(rep(res_cc, ncol(mat)))
  # }
  
  out_matL <- matL(mat, h)
  out_matR <- matR(mat, h)
  
  ## computing pencils
  rCumL <- rowCumsums(out_matL) # p x (h+1) matrix
  rcCumL <- colCumsums(rCumL) # p x (h+1) matrix
  rm(rCumL)
    
  rCumR <- rowCumsums(out_matR) # p x (h+1) matrix
  rcCumR <- colCumsums(rCumR) # p x (h+1) matrix
  rm(rCumR)
  
  ## Initialization:
  ##
  maxSize <- 3*p
  ## NB: The size of heap, D and chainedL are set to the maximum size required,
  ## ie 3*p. This comes from: size p at initialization; at each of the p
  ## iterations, at most 2 new merges are added (only one if the last merge
  ## involves the first or last cluster).

  gains <- rep(0, p-1)
  merge <- matrix(0, nrow = p-1, ncol = 2)  # matrix of the merges
  traceW <- matrix(0, nrow = p-1, ncol = 2) # matrix of traceW
  sd1 <- out_matL[1:(p-1),2]/2 # similarity of objects with their right neighbors
  sii <- out_matL[1:p,1] # auto-similarity of objects
  rm(out_matL)
    
  ## initialization of the heap
  heap <- as.integer(rep(-1, maxSize))
  lHeap <- length(heap)
  v <- 1:(p-1)
  heap[v] <- v
  D <- rep(-1, maxSize)
  
  ## linkage value of objects with their right neighbors
  D[v] <- (sii[v] + sii[v+1]) / 2 - sd1 
  ## initialization of the length of the Heap
  lHeap <- p-1
  ## each element contains a vector: c(cl1, cl2, label1, label2, posL, posR, valid)
  chainedL <- matrix(-1, nrow = 12, ncol = maxSize)
  rownames(chainedL) <- c("minCl1", "maxCl1", "minCl2", "maxCl2", "lab1", 
                          "lab2", "posL", "posR", "Sii", "Sjj", "Sij", "valid")
  w <- as.integer(v + 1)
  chainedL[1,v] <- v
  chainedL[2,v] <- v
  chainedL[3,v] <- w
  chainedL[4,v] <- w
  chainedL[5,v] <- -v
  chainedL[6,v] <- -w
  chainedL[7,v] <- v - 1
  chainedL[8,v] <- w
  chainedL[9,v] <- sii[v]
  chainedL[10,v] <- sii[v+1]
  chainedL[11,v] <- sd1
  chainedL[12,v] <- 1
  chainedL[7,1] <- -1
  chainedL[8,p-1] <- -1

  heap <- adjclust:::buildHeap(heap, D, lHeap)
  # performing clustering  
  res <- .Call("cWardHeaps", rcCumR, rcCumL, as.integer(h), as.integer(p), 
               chainedL, heap, D, as.integer(lHeap), merge, gains, traceW, 
               PACKAGE = "adjclust")
    
  # formatting outputs (as 'hclust') and checking if decreasing heights are present
  height <- gains
  if (is.null(rownames(mat))) {
    labels <- as.character(1:p)
  } else {
    labels <- rownames(mat)
  }
  tree <- list(merge = res,
               height = height,
               order = 1:p,
               labels = labels,
               call = match.call(),
               method = "adjClust",
               dist.method = attr(D, "method"),
               data = mat)
  class(tree) <- c("chac")
    
  if (any(diff(tree$height) < 0)) 
    message(paste("Note:", sum(diff(tree$height) < 0),
                  "merges with non increasing heights."))
    
  return(tree)
}