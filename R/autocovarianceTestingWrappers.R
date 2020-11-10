## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#' @useDynLib autocovarianceTesting, .registration=TRUE
#' 

# Function to run compatibility checks
compatibilityChecks <- function(X, Y, L, test, alpha, B, prewhiten){
  
  # X and Y must be matrices
  if ( (!(is.matrix(X))) | (!(is.matrix(Y))) ){
    stop(paste("X and Y must be matrices"))
  }
  
  # X and Y have the same number of rows
  if ( nrow(X) != nrow(Y) ){
    stop(paste("X and Y must be matrices with the same number of rows"))
  }
  
  # X and Y have the same number of columns
  if ( ncol(X) != ncol(Y) ){
    stop(paste("X and Y must be matrices with the same number of columns"))
  }
  
  # L must be numeric
  if ( !is.numeric(L) ){
    stop(paste("L must be numeric"))
  }
  
  # L must be less than n
  if ( L >= nrow(X) ){
    stop(paste("L must be strictly less than n the length of the series"))
  }
  
  # L must be larger than n
  if ( ncol(X) >= nrow(X) ){
    stop(paste("The dimension of the series must be smaller than its length"))
  }
  
  # test must be one of c(Independent, Dependent, bootDependent, bootBartlett)
  if (!(all(test %in% c("Independent", "Dependent", "bootDependent", "bootBartlett")))){
    stop(paste("test must be one of Independent, Dependent, bootDependent, or bootBartlett"))
  }
  
  # alpha must be numeric
  if ( !is.numeric(alpha) ){
    stop(paste("alpha must be numeric"))
  }
  
  # alpha between 0 and 1
  if ( alpha >= 1 | alpha <= 0 ){
    stop(paste("alpha must be between 0 and 1"))
  }
  
  # B must be numeric
  if ( !is.numeric(B) ){
    stop(paste("B must be numeric"))
  }
  
  # prewhiten must be boolean
  if ( !is.logical(prewhiten) ){
    stop(paste("prewhiten must be logical"))
  }
  
}

#' Title
#'
#' @param X 
#' @param Y 
#' @param L 
#' @param test 
#' @param alpha 
#' @param B 
#' @param prewhiten 
#'
#' @return
#' @export
#'
#' @examples
autocovarianceTest <- function(X, Y, L, test = "Dependent", alpha = 0.05, B = 500, prewhiten = FALSE){
  
  # Compatibility Tests
  compatibilityChecks(X, Y, L, test, alpha, B, prewhiten)
  
  # Initialize output list
  out <- list()
  
  # Compute independent and dependent covariance if need be
  if ("Independent" %in% test | "Dependent" %in% test){
    # Compute dependent covariance
    testStatistics <- calculateCovariance(X, Y, L)
  }
  
  
}