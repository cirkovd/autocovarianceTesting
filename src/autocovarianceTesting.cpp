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
Rcpp::List calculateCovariance(const arma::mat & X, const arma::mat & Y, const double & L) {
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
    
    // Compute independent covariance matrix
    arma::mat onemat(k, k, fill::ones);
    arma::colvec ind_vec = vectorise(kron(arma::eye(2, 2), onemat));
    arma::mat ind_mat = ind_vec * trans(ind_vec);
    arma::mat ind_mat2(L + 1, L + 1, fill::ones);
    arma::mat large_ind_mat = kron(ind_mat2, ind_mat);
    arma::mat onemat2(2 * pow(k, 2), 2 * pow(k, 2), fill::ones);
    large_ind_mat = kron(ind_mat2, kron(arma::eye(2, 2), onemat2)) % large_ind_mat;
        
    // Get our delta vector (difference in autocovariance up to lag L)
    arma::colvec delta = A * eta;
    // delta covariance matrix
    arma::mat dep_cov = A * W * trans(A);
    // independent covariance matrix
    arma::mat ind_cov = A * (large_ind_mat % W) * trans(A);
    
    if ( k > 1 ){
        // Find indices with duplicate lag 0 autocovariances
        arma::mat dup(k, k, fill::ones);
        dup = arma::trimatu(dup, 1);
        arma::uvec dup_ind = find(vectorise(dup));
        
        // Remove duplicates
        delta.shed_rows(dup_ind);
        dep_cov.shed_rows(dup_ind);
        dep_cov.shed_cols(dup_ind);
        ind_cov.shed_rows(dup_ind);
        ind_cov.shed_cols(dup_ind);
    }
    
    return Rcpp::List::create(Rcpp::Named("delta") = delta, 
                              Rcpp::Named("dep_cov") = dep_cov,
                              Rcpp::Named("ind_cov") = ind_cov);
    
}

// Function that computes test statistics for fixed lag tests given autocovariance differences
// and asymptotic covariance matrix
// [[Rcpp::export]]
Rcpp::List calculateTestStat(const arma::colvec & delta, const arma::mat & covar, const int & n, const int & L, const int & k) {
    // Compute unweighted stat
    double stat = n * as_scalar(trans(delta) * solve(covar, delta));  
    
    // Compute weighted stat
    arma::colvec weights = arma::linspace<arma::colvec>(1, 1/(L + 1), L + 2);
    weights.shed_row(L + 1);
    weights = repelem(weights, k * k, 1);
    // Get rid of duplicates
    int p1 = delta.n_rows;
    int p2 = weights.n_rows;
    weights = weights.rows(p2 - p1, p2 - 1);
    
    // Weighted statistic
    double weight_stat = n * as_scalar(trans(weights) * (delta % delta));
    double K1 = trace(covar * diagmat(weights));
    double K2 = 2 * trace(covar * diagmat(weights) * covar * diagmat(weights));
    double alpha = pow(K1, 2) / K2;
    double beta = K2 / K1;
    
    return Rcpp::List::create(Rcpp::Named("stat") = stat,
                              Rcpp::Named("df") = p1,
                              Rcpp::Named("weight_stat") = weight_stat,
                              Rcpp::Named("alpha") = alpha,
                              Rcpp::Named("beta") = beta);
}


// Function that prewhitens time series data via an AR(p) process for the bootstrap algorithm
arma::mat prewhitenData(const arma::mat & Z){
    // Get some details from the inputs
    int n = Z.n_rows; // time series length
    int k = Z.n_cols / 2; // time series dimension
    
    // AR order
    int p = ceil(std::min(10 * log(n), 0.8 * sqrt(n)));
    // Create AR(p) matrix for linear regression
    arma::mat V(n, 2 * k * p, fill::zeros);
    for (int i = 1; i < p + 1; i++){
        V.cols(2 * k * (i - 1), 2 * k * i - 1) = shift(Z, +(i - 1), 0);
    }
    V.shed_row(n - 1);
    V.shed_rows(0, p - 2);
    arma::colvec intercept(n - p, fill::ones);
    V = join_rows(intercept, V);
    
    // Compute AR coefficients
    arma::mat betas = solve(trans(V) * V, trans(V) * Z.rows(p, n - 1));  
    
    // Incorporate H0 for beta coefficients
    // Average Mean
    arma::rowvec mu = betas(span(0), span(0, k - 1)) + betas(span(0), span(k, 2 * k - 1));
    arma::rowvec mu_avg = join_rows(mu, mu);
    betas.shed_row(0);
    
    // Average Phis
    arma::mat betas_sub(2 * k, 2 * k, fill::zeros);
    arma::mat beta_col(2 * k, k, fill::zeros);
    int lower = 0;
    int upper = 0;
    // Loop to sum beta matrices
    for (int i = 0; i < p; i++){
        upper = lower + 2 * k - 1;
        betas_sub = betas.rows(lower, upper);
        betas_sub.cols(k, 2 * k - 1) = arma::shift(betas_sub.cols(k, 2 * k - 1), +k, 0);
        beta_col = betas_sub.cols(0, k - 1) + betas_sub.cols(k, 2 * k - 1);
        betas.rows(lower, upper) = join_rows(beta_col, arma::shift(beta_col, +k, 0));
        lower = upper + 1;
    }
    // Divide whole matrix by two to get average
    betas = join_cols(mu_avg, betas);
    arma::colvec half_vec(betas.n_rows, fill::zeros);
    half_vec.fill(0.5);
    betas = betas.each_col() % half_vec;
    // Finally prewhiten data
    arma::mat Zwhitened = Z.rows(p, n - 1) - V * betas;
    
    return Zwhitened;
}

// Function that computes test statistics for order lag tests 
// [[Rcpp::export]]
Rcpp::List calculateBootTestStat(const arma::mat & X, const arma::mat & Y, const double & L, double const & alpha, int const & B, bool const & prewhiten) {
    // Get some details from the inputs
    int n = X.n_rows; // time series length
    int k = X.n_cols; // time series dimension
    
    // Concatenate the series
    arma::mat Z(n, 2 * k);
    Z(span(0, n - 1), span(0, k - 1)) = X; 
    Z(span(0, n - 1), span(k, 2 * k - 1)) = Y; 
    
    // Prewhiten data if asked for
    if (prewhiten == true){
        Z = prewhitenData(Z);
        n = Z.n_rows;
    }
    
    // Set some parameters for block of blocks bootstrap
    int b = floor(std::max(0.5 * cbrt(n), 2.0));
    // no blocks
    int K = ceil(n/b);
    int L_max = 0;
    if ( L == NULL ){
        L_max = ceil(pow(log2(n), 0.999));
    } else {
        L_max = L;
    }
    
    // Do our method but bootstrapping the covariance
    Rcpp::List covar_stats = calculateCovariance(Z.cols(0, k - 1), Z.cols(k, 2 * k - 1), L);
    arma::colvec stats(L_max + 1, fill::zeros);
    arma::colvec delta = covar_stats["delta"];
    arma::mat covar = covar_stats["dep_cov"];
    for (int i = 0; i < L_max + 1; i++){
        stats(i) = n * as_scalar( trans(delta.rows(0, i)) * solve(covar(span(0, i), span(0, i)), delta.rows(0, i)) );
    }
    
    // Do bootstrap algorithm
    arma::mat X_Z = trans(Z.cols(0, k - 1));
    arma::mat Y_Z = trans(Z.cols(k, 2 * k - 1));
    arma::mat shiftedX(k, n);
    arma::mat shiftedY(k, n);
    arma::mat D_X(k * k * (L_max + 1), n);
    arma::mat D_Y(k * k * (L_max + 1), n);
    
    
    arma::colvec one_vec(k, fill::ones);
    arma::colvec zero_vec(k, fill::zeros);
    D_X.rows(0, k * k - 1) = (kron(one_vec, X_Z) % kron(X_Z, one_vec)).eval();
    D_Y.rows(0, k * k - 1) = (kron(one_vec, Y_Z) % kron(Y_Z, one_vec)).eval();
    
    int lower = k * k;
    int upper = 0;
    double epsX = max(vectorise(X_Z));
    double epsY = max(vectorise(Y_Z));
    // Compute D matrices
    for (int l = 1; l < L_max + 1; l++){
        upper = lower + k * k - 1;
        
        // Compute D_X for every l
        shiftedX = shift(X_Z, -l, 1);
        shiftedX.cols(n - l, n - 1).clean(epsX * epsX + 10000);
        D_X.rows(lower, upper) = (kron(one_vec, shiftedX) % kron(X_Z, one_vec)).eval();
            
        // Compute D_Y for every l
        shiftedY = shift(Y_Z, -l, 1);
        shiftedY.cols(n - l, n - 1).clean(epsY * epsY + 10000);
        D_Y.rows(lower, upper) = (kron(one_vec, shiftedY) % kron(Y_Z, one_vec)).eval();
        
        lower = upper + 1; 
    }
    
    arma::mat D = D_X - D_Y;
    
    
    
    return Rcpp::List::create(Rcpp::Named("D") = D,
                              Rcpp::Named("eps") = epsX);
}



