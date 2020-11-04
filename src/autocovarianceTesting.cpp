// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function that, given two multivariate time series, computes the DEMEANED autocovariance function
// of the concatenated series

//' Compute computes the DEMEANED autocovariance function of two concatenated multivariate time series
//' 
//' @param X a time series
//' @param Y a time series
//' @export
// [[Rcpp::export]]
arma::cube calculateAutocovariance(const arma::mat & X, const arma::mat & Y){
    // Get some details from the inputs
    int n = X.n_rows; // time series length
    int k = X.n_cols; // time series dimension
    
    // Concatenate the series
    arma::mat Z(n, 2 * k);
    Z(span(0, n - 1), span(0, k - 1)) = X; 
    Z(span(0, n - 1), span(k, 2*k - 1)) = Y; 
    
    // Center the time series
    arma::mat ZCenter = Z.each_row() - arma::mean(Z, 0);
    
    // Preallocate acvf matrices
    arma::cube ZZt(2 * k, 2 * k, n, fill::zeros);
    arma::cube gammaZ(2 * k, 2 * k, n, fill::zeros);
    
    for (int h = 0; h < n; h++){
        for(int t = 0; t < n - h; t++){
            // Compute summand
            ZZt.slice(t) = ZCenter.row(t + h).t() * ZCenter.row(t);
        }
        // Compute acf
        gammaZ.slice(h) = arma::mean(ZZt, 2);
        
        // Reset summand to zero
        ZZt.zeros();
    }
    
    // Bind the covariance matrices into one large list of covariances
    arma::cube acf(2 * k, 2 * k, 2 * (n - 1) + 1, fill::zeros);
    // positive lagged acf
    acf(span(0, 2 * k - 1), span(0, 2 * k - 1), span(n - 1, 2*n - 2)) = gammaZ;
    
    // loop to set up negative lag acf
    for (int h = 0; h < n - 1; h++){
        acf.slice(h) = acf.slice(2 * n - 2 - h).t();
    }
    
    return acf;
}

// Function that computes the asymptotic covariance matrix of the autocovariance differences 
// when the two time series are dependent
Rcpp::List dependentCovariance(const arma::mat & X, const arma::mat & Y, const double & L) {
    // Get some details from the inputs
    int n = X.n_rows; // time series length
    int k = X.n_cols; // time series dimension
}

