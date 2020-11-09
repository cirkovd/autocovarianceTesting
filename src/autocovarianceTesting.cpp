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
arma::cube calculateAutocovariance(const arma::mat & X, const arma::mat & Y, const int & maxLag){
    // Get some details from the inputs
    int n = X.n_rows; // time series length
    int k = X.n_cols; // time series dimension
    
    // Concatenate the series
    arma::mat Z(n, 2 * k);
    Z(span(0, n - 1), span(0, k - 1)) = X; 
    Z(span(0, n - 1), span(k, 2 * k - 1)) = Y; 
    
    // Center the time series
    arma::mat ZCenter = Z.each_row() - arma::mean(Z, 0);
    
    // Preallocate acvf matrices
    arma::cube ZZt(2 * k, 2 * k, n, fill::zeros);
    arma::cube gammaZ(2 * k, 2 * k, maxLag, fill::zeros);
    
    for (int h = 0; h < maxLag ; h++){
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
    arma::cube acf(2 * k, 2 * k, 2 * (maxLag - 1) + 1, fill::zeros);
    // positive lagged acf
    acf(span(0, 2 * k - 1), span(0, 2 * k - 1), span(maxLag - 1, 2*maxLag - 2)) = gammaZ;
    
    // loop to set up negative lag acf
    for (int h = 0; h < maxLag - 1; h++){
        acf.slice(h) = acf.slice(2 * maxLag - 2 - h).t();
    }
    
    return acf;
}

double bartlettForumla(const arma::cube & acf_H0, const int trunc, const int a, const int b, const int c, const int d, const int p, const int q, const int n){
    arma::colvec cov = acf_H0(span(a, a), span(c, c), span(n - 1 - trunc,  n - 1 + trunc)) % acf_H0(span(b, b), span(d, d), span(n - 1 - trunc - p + q,  n - 1 + trunc - p + q)) + acf_H0(span(a, a), span(d, d), span(n - 1 - trunc + q,  n - 1 + trunc + q)) % acf_H0(span(b, b), span(c, c), span(n - 1 - trunc - p,  n - 1 + trunc - p));
    double covar = sum(cov);
    return covar;
}

// Function that computes the asymptotic covariance matrix of the autocovariance differences 
// when the two time series are dependent
// [[Rcpp::export]]
Rcpp::List dependentCovariance(const arma::mat & X, const arma::mat & Y, const double & L) {
    // Get some details from the inputs
    int n = X.n_rows; // time series length
    int k = X.n_cols; // time series dimension
    
    // truncation for asymptotic covariance
    int trunc = floor(cbrt(n));
    
    int maxLag = trunc + L + 2;
    // Get autocovariance function
    arma::cube acf = calculateAutocovariance(X, Y, maxLag);
    
    // Get our autocovariances up to lag L
    arma::colvec eta = vectorise(acf.slices(maxLag - 1, maxLag + L - 1));
    
    // Create the contrast matrix used in the test statistic
    arma::mat cont_mat(2 * k, 2 * k, fill::zeros);
    arma::mat A_sub(pow(k, 2), pow(2 * k, 2), fill::zeros);
    int counter = 0;
    
    // Create contrast matrix submatrix
    for (int j = 0; j < k; j++){
        for (int i = 0; i < k; i++){
            // Create mini matrices that when vectorized give rows of A
            cont_mat(i, j) = 1;
            cont_mat(i + k, j + k) = -1;
            A_sub.row(counter) = trans(vectorise(cont_mat));
            
            // reset matrix
            cont_mat.zeros();
            counter = counter + 1;
        }
    }
    // Create large contrast matrix 
    arma::mat A = kron(arma::eye(L + 1, L + 1), A_sub);
    
    // Get autocovariances under the null hypothesis (average the block diagonals)
    arma::cube acf_H0 = acf;
    // Sum and replace the blog diagonals
    arma::cube acf_sum = acf_H0(span(0, k - 1), span(0, k - 1), span(0, 2 * maxLag - 2)) + acf_H0(span(k, 2 * k - 1), span(k, 2 * k - 1), span(0, 2 * maxLag - 2)); 
    acf_H0(span(0, k - 1), span(0, k - 1), span(0, 2 * maxLag - 2)) = acf_sum;
    acf_H0(span(k, 2 * k - 1), span(k, 2 * k - 1), span(0, 2 * maxLag - 2)) = acf_sum;
    
    // Create a matrix with 0.5 to average the autocovariances
    arma::mat half(k, k, fill::zeros);
    half.fill(0.5);
    arma::mat one(k, k, fill::ones);
    half = kron(arma::eye(2, 2), half) + kron(arma::fliplr(arma::eye(2, 2)), one);
    
    // Divide block diagonals by 2
    for (int h = 0; h < 2 * maxLag - 2; h++){
        acf_H0.slice(h) = acf_H0.slice(h) % half;
    }
    
    // Compute asymptotic covariance matrix
    int Wdim = pow(2 * k, 2) * (L + 1);
    arma::mat W(Wdim, Wdim, fill::zeros);
    int row = 0;
    int col = 0;
    
    // Fill only upper diagonal
    for (int lag1 = 0; lag1 < L + 1; lag1++){
        for (int b = 0; b < 2 * k; b++){
            for (int a = 0; a < 2 * k; a++){
                for (int lag2 = 0; lag2 < L + 1; lag2++){
                    for (int d = 0; d < 2 * k; d++){
                        for (int c = 0; c < 2 * k; c++){
                            // Apply Bartlett's Formula
                            W(row, col) = bartlettForumla(acf_H0, trunc, a, b, c, d, lag1, lag2, maxLag);
                            col = (col + 1) % (Wdim);
                        }
                    }
                }
                row = (row + 1) % (Wdim);
            }
        }
    }
    
    // Get our delta vector (difference in autocovariance up to lag L)
    arma::colvec delta = A * eta;
    // delta covariance matrix
    arma::mat dep_cov = A * W * trans(A);
    
    if ( k > 1 ){
        // Find indices with duplicate lag 0 autocovariances
        arma::mat dup(k, k, fill::ones);
        dup = arma::trimatu(dup, 1);
        arma::uvec dup_ind = find(vectorise(dup));
        
        // Remove duplicates
        delta.shed_rows(dup_ind);
        dep_cov.shed_rows(dup_ind);
        dep_cov.shed_cols(dup_ind);
    }
    
    return Rcpp::List::create(delta, dep_cov);
    
}

