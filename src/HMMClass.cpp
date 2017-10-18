#include <Rcpp.h>
using namespace Rcpp;

// Define numericList type
typedef Rcpp::ListOf<Rcpp::NumericVector> NumericVectorList;

double logSumExpCpp(const NumericVector& x){
  double a=max(x);
  return a+log(sum(exp(x-a)));
}

// [[Rcpp::export]] 
NumericMatrix forwardAlgorithmCpp(NumericMatrix& e, NumericMatrix& t, NumericVector& prior) {
  int nobs = e.nrow();
  int nsta = e.ncol();
  // Initialize alpha table
  NumericMatrix aMat(nobs,nsta); 
  //Fill first row of alpha table
  for(int j=0;j < nsta; j++){
    aMat(0,j)=prior(j)+e(0,j);
  }
  // Loop over alpha table to fill in whole thing
  if( nobs > 1 ){
    // iterate over loci
    for(int i=1;i < nobs; i++){
      // iterate over states
      for(int j=0;j < nsta; j++){
	NumericVector foo(nsta);
	// iterate over prior states
	for(int k=0; k<nsta;k++){
	  foo(k)=aMat(i-1,k)+t(k,j);
	}
	aMat(i,j) = logSumExpCpp(foo) + e(i,j);
      }
    }
  }
  return(aMat);
}

// [[Rcpp::export]] 
NumericMatrix forwardAlgorithmSparseCpp(NumericMatrix& e, NumericMatrix& t, NumericVector& prior,  IntegerVector permTrans, IntegerVector tLen) {
  int nobs = e.nrow();
  int nsta = e.ncol();
  // Initialize alpha table
  NumericMatrix aMat(nobs,nsta);
  //Fill first row of alpha table
  for(int j=0;j < nsta; j++){
    aMat(0,j)=prior(j)+e(0,j);
  }
  // Loop over alpha table to fill in whole thing
  if( nobs > 1 ){
    // iterate over loci
    for(int i=1;i < nobs; i++){
      int position=0; //track position in acceptable transition vector
      // iterate over states
      for(int j=0;j < nsta; j++){
	NumericVector foo(tLen(j));
	// iterate over permissible prior states
	for(int k=position; k<position+tLen(j);k++){
	  foo(k-position)=aMat(i-1,permTrans(k))+t(permTrans(k),j);
	}
	aMat(i,j) = logSumExpCpp(foo) + e(i,j);
	position=position+tLen(j);
      }
    }
  }
  return(aMat);
}


// [[Rcpp::export]] 
NumericMatrix backwardAlgorithmCpp(NumericMatrix& e, NumericMatrix& t) {
  int nobs = e.nrow();
  int nsta = e.ncol();
  // Initialize alpha table
  NumericMatrix bMat(nobs,nsta); 
  //Fill first row of beta table
  for(int j=0;j < nsta; j++){
    bMat(nobs-1,j)=0;
  }
  // Loop over alpha table to fill in whole thing
  if( nobs > 1 ){
    // iterate over loci
    for(int i=nobs-2;i >= 0; i--){
      // iterate over states
      for(int j=0;j < nsta; j++){
	NumericVector foo(nsta);
	// iterate over prior states
	for(int k=0; k < nsta;k++){
	  foo(k)=bMat(i+1,k) + e(i+1,k) + t(j,k) ;
	}
	bMat(i,j) = logSumExpCpp(foo);
      }
    }
  }
  return(bMat);
}

// [[Rcpp::export]] 
NumericMatrix backwardAlgorithmSparseCpp(NumericMatrix& e, NumericMatrix& t,  IntegerVector permTrans, IntegerVector tLen) {
  int nobs = e.nrow();
  int nsta = e.ncol();
  // Initialize beta table
  NumericMatrix bMat(nobs,nsta); 
  //Fill first row of beta table
  for(int j=0;j < nsta; j++){
    bMat(nobs-1,j)=0;
  }
  // Loop over beta table to fill in whole thing
  if( nobs > 1 ){
    // iterate over loci
    for(int i=nobs-2;i >= 0; i--){
      int position=0; //track position in acceptable transition vector
      // iterate over states
      for(int j=0;j < nsta; j++){
	NumericVector foo(tLen(j));
	// iterate over prior states
	for(int k=position; k<position+tLen(j);k++){
	  foo(k-position)=bMat(i+1,permTrans(k)) + e(i+1,permTrans(k)) + t(j,permTrans(k)) ;
	}
	bMat(i,j) = logSumExpCpp(foo);
	position=position+tLen(j);
      }
    }
  }
  return(bMat);
}



