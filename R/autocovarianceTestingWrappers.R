## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#' @useDynLib autocovarianceTesting, .registration=TRUE
#' 

# Function to run compatibility checks
compatibilityChecks <- function(X, Y, L, test, trunc, B, b, prewhiten, dependent_series, bootstrap_fixed_stats, weights, average_bartlett_cov){
  
  # X and Y must have the same class
  if ( class(X)[[1]] != class(Y)[[1]] ){
    stop(paste("X and Y must have the same class"))
  }
  
  # X and Y must be matrices
  if ( (!(is.matrix(X))) ){
    if ( (!(stats::is.ts(X))) ){
      if ( (!(is.vector(X))) ){
        stop(paste("X and Y must be matrices, ts objects, or vectors"))
      }
    }
  }
  
  # If X and Y are vectors, put them in the correct format
  X <- t(t(X))
  Y <- t(t(Y))
  
  # Get ts length
  n <- nrow(X)
  
  # X and Y have the same number of rows
  if ( nrow(X) != nrow(Y) ){
    stop(paste("X and Y must be matrices with the same number of rows"))
  }
  
  # X and Y have the same number of rows
  if ( nrow(X) < 3 ){
    stop(paste("X and Y must be must be of length greater than 2"))
  }
  
  # X and Y have the same number of columns
  if ( ncol(X) != ncol(Y) ){
    stop(paste("X and Y must be matrices with the same number of columns"))
  }
  
  # X and Y have no missing values
  if ( sum(is.na(X)) > 0 | sum(is.na(Y)) > 0 ){
    stop(paste("X and Y can not have missing values"))
  }
  
  # AR(p) order for prewhitening
  p <- ifelse(prewhiten == TRUE & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test), ceiling(min(10 * log(n), 0.8 * sqrt(n))), 0)
  
  # Set L if not supplied
  if ( is.null(L) ){
    L <- ceiling((log2(n - p))^0.9999)
  }
  
  # L must be numeric
  if ( !is.numeric(L) ){
    stop(paste("L must be numeric and a positive integer"))
  }
  
  # L must be nonnegative
  if ( !((L == round(L)) & (L >= 0)) ){
    stop(paste("L must be numeric and a positive integer"))
  }
  
  # L must be less than n
  if ( L >= (n - p - 1) ){
    stop(paste("L must be strictly less than n the length of the series (minus AR order + 1 if prewhitening)"))
  }
  
  # Set b if not supplied
  if ( is.null(b) ){
    b <- max(floor(0.5 * (n - p)^(1/3)), 2)
  }
  
  # b must be numeric
  if ( !is.numeric(b) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    stop(paste("b must be an integer and larger than 2"))
  }
  
  # b must be larger than 2
  if ( !((b == round(b)) & (b >= 2)) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    stop(paste("b must be an integer and larger than 2"))
  }
  
  # b must be less than n
  if ( b >= (n - p - 1) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    stop(paste("b must be strictly less than n the length of the series (minus AR order + 1 if prewhitening)"))
  }
  
  # Dimension of time series must be smaller than n
  if ( ncol(X) >= (n - p - 1)){
    stop(paste("The dimension of the series must be smaller than its length (minus AR order + 1 if prewhitening)"))
  }
  
  # test must be one of c(Fixed_Lag, Weighted_Fixed_Lag, Auto_Lag_Jin, Auto_Lag_Bartlett)
  if (!(all(test %in% c("Fixed_Lag", "Weighted_Fixed_Lag", "Auto_Lag_Jin", "Auto_Lag_Bartlett")))){
    stop(paste("test must be one of Fixed_Lag, Weighted_Fixed_Lag, Auto_Lag_Jin, or Auto_Lag_Bartlett"))
  }
  
  # B must be numeric
  if ( !is.numeric(B) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    stop(paste("B must be an integer and larger than 1"))
  }
  
  # B must be larger than 1
  if (!((B == round(B)) & (B >= 1)) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    stop(paste("B must be an integer and larger than 1"))
  }
  
  # Set trunc if not supplied
  if ( is.null(trunc) ){
    trunc <- floor((n - p)^(1/3))
  }
  
  # trunc must be numeric
  if ( !is.numeric(trunc) ){
    stop(paste("trunc must be numeric and a positive integer"))
  }
  
  # trunc must be greater than 0
  if ( !((trunc == round(trunc)) & (trunc >= 0)) ){
    stop(paste("trunc must be numeric and a positive integer"))
  }
  
  # trunc must be less than n
  if ( trunc >= (n - p - L - 1) ){
    stop(paste("trunc must be strictly less than n - L (minus AR order + 1 if prewhitening)"))
  }
  
  # prewhiten must be boolean
  if ((prewhiten == 0 | prewhiten == 1) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    prewhiten = ifelse(prewhiten == 1, TRUE, FALSE)
  } else {
    if ( !is.logical(prewhiten) & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
      stop(paste("prewhiten must be logical"))
    }
  }
  
  # b Warning
  if ( b < 3 & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    warning("b is small (b < 3), bootstrap resamples may not be representative")
  }
  
  # B warning
  if ( B < 100 & ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    warning("B is small (B < 100), results may not be reliable")
  }
  
  # dependent_series must be boolean
  if ((dependent_series == 0 | dependent_series == 1) & ("Fixed_Lag" %in% test | "Weighted_Fixed_Lag" %in% test | "Auto_Lag_Bartlett" %in% test) ){
    dependent_series = ifelse(dependent_series == 1, TRUE, FALSE)
  } else {
    if ( !is.logical(dependent_series) & ("Fixed_Lag" %in% test | "Weighted_Fixed_Lag" %in% test | "Auto_Lag_Bartlett" %in% test) ){
      stop(paste("dependent_series must be logical"))
    }
  }
  
  # bootstrap_fixed_stats must be boolean
  if ((bootstrap_fixed_stats == 0 | bootstrap_fixed_stats == 1) & ("Fixed_Lag" %in% test | "Weighted_Fixed_Lag" %in% test) ){
    bootstrap_fixed_stats = ifelse(bootstrap_fixed_stats == 1, TRUE, FALSE)
  } else {
    if ( !is.logical(bootstrap_fixed_stats) & ("Fixed_Lag" %in% test | "Weighted_Fixed_Lag" %in% test) ){
      stop(paste("bootstrap_fixed_stats must be logical"))
    }
  }
  
  # Weight compatibility checks
  if ( is.null(weights) ){
    weights <- seq(1, 1/(L + 1), by = -1/(L + 1))
  }
  
  if ( !is.vector(weights) | !is.numeric(weights) ){
    stop(paste("weights must be a numeric vector"))
  }
  
  if ( !(length(weights) == L + 1) ){
    stop(paste("weights must be of length max_lag + 1"))
  }
  
  if ( !all(weights <= 1) | !all(weights > 0)){
    stop(paste("weights must be of length max_lag + 1"))
  }
  
  # average_bartlett_cov must be boolean
  if ((average_bartlett_cov == 0 | average_bartlett_cov == 1) & ("Auto_Lag_Bartlett" %in% test) ){
    average_bartlett_cov = ifelse(average_bartlett_cov == 1, TRUE, FALSE)
  } else {
    if ( !is.logical(average_bartlett_cov) & ("Auto_Lag_Bartlett" %in% test) ){
      stop(paste("average_bartlett_cov must be logical"))
    }
  }
  
  return(list(L, b, trunc, X, Y, prewhiten, B, dependent_series, bootstrap_fixed_stats, weights, average_bartlett_cov))
}


#' Summarizing autocovariance_test function output
#'
#' @param x a \code{"acvfTest"} object given by \code{"autocovariance_test()"}.
#' @param ... additional arguments to \link[base]{print}.
#'
#' @export
#'
#' @examples
#' #' ## Simulate two marignal AR(1) processes that are correlated at lag zero and run tests
#' set.seed(12345)
#' 
#' # Simulate correlated errors
#' Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2) # covariance of errors
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X <- matrix(arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1])) 
#' Y <- matrix(arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 2]))
#' 
#' # Run tests for equality of autocovariance functions up to lag 5
#' output <- autocovariance_test(X, Y, max_lag = 5, 
#' test = c("Fixed_Lag","Weighted_Fixed_Lag","Auto_Lag_Jin","Auto_Lag_Bartlett"))
#' # All tests fail to reject the null
#' print(output)
print.acvfTest <- function(x, ...){
  # Create test vector
  test <- c()
  listnames <- names(x)
  test <- if("fixed_lag_stat" %in% listnames | "boot_fixed_lag_stat" %in% listnames){ c(test, "Fixed_Lag") }else{ test }
  test <- if("fixed_lag_weight_stat" %in% listnames | "boot_fixed_lag_weight_stat" %in% listnames){ c(test, "Weighted_Fixed_Lag") }else{ test }
  test <- if("auto_lag_jin_stat" %in% listnames){ c(test, "Auto_Lag_Jin") }else{ test }
  test <- if("auto_lag_bart_stat" %in% listnames){ c(test, "Auto_Lag_Bartlett") }else{ test }
  
  if( "Fixed_Lag" %in% test & x$bootstrap_fixed_stats == FALSE){
    # Initialize Fixed Lag Table
    out_table <- data.frame(0, 1, 2, 3)
    out_table <- out_table[-1, ]
  }
  
  if( "Weighted_Fixed_Lag" %in% test & x$bootstrap_fixed_stats == FALSE){
    # Initialize Weighted Table
    out_table_weight <- data.frame(0, 1, 2, 3, 4)
    out_table_weight <- out_table_weight[-1, ]
  }
  
  if( "Fixed_Lag" %in% test & x$bootstrap_fixed_stats == TRUE){
    # Initialize Fixed Lag Table
    out_table <- data.frame(0, 1, 2)
    out_table <- out_table[-1, ]
  }
  
  if( "Weighted_Fixed_Lag" %in% test & x$bootstrap_fixed_stats == TRUE){
    # Initialize Weighted Table
    out_table_weight <- data.frame(0, 1, 2)
    out_table_weight <- out_table_weight[-1, ]
  }
  
  if ( "Fixed_Lag" %in% test & x$bootstrap_fixed_stats == FALSE){
    # Add Fixed_Lag tests
    out_table <- rbind(out_table, c(ifelse(x$dependent_series == FALSE, "Independent", "Dependent"), round(x$fixed_lag_stat, 3), x$fixed_lag_df, round(x$fixed_lag_pval, 3)))
    colnames(out_table) <- c("Test", "Statistic", "df", "p-value")
  }
  
  if ( "Weighted_Fixed_Lag" %in% test & x$bootstrap_fixed_stats == FALSE){
    # Add Weighted_Fixed_Lag tests
    out_table_weight <- rbind(out_table_weight, c(ifelse(x$dependent_series == FALSE, "Independent", "Dependent"), round(x$fixed_lag_weight_stat, 3), round(x$fixed_lag_weight_alpha, 3), round(x$fixed_lag_weight_beta, 3), round(x$fixed_lag_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Statistic", "alpha", "beta", "p-value")
  }
  
  # Bootstrapped Fixed Lag Tables
  
  if ( "Fixed_Lag" %in% test & x$bootstrap_fixed_stats == TRUE){
    # Add Fixed_Lag tests
    out_table <- rbind(out_table, c(ifelse(x$dependent_series == FALSE, "Independent", "Dependent"), round(x$boot_fixed_lag_stat, 3), round(x$boot_fixed_lag_pval, 3)))
    colnames(out_table) <- c("Test", "Statistic", "p-value")
  }
  
  if ( "Weighted_Fixed_Lag" %in% test & x$bootstrap_fixed_stats == TRUE){
    # Add Weighted_Fixed_Lag tests
    out_table_weight <- rbind(out_table_weight, c(ifelse(x$dependent_series == FALSE, "Independent", "Dependent"), round(x$boot_fixed_lag_weight_stat, 3), round(x$boot_fixed_lag_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Statistic", "p-value")
  }
  
  
  # Initialize bootstrapped tables
  if ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test){
    out_table_boot <- data.frame(0, 1, 2, 3)
    out_table_boot <- out_table_boot[-1, ]
  }
  
  if ("Auto_Lag_Jin" %in% test){
    # Add jin tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Jin", round(x$auto_lag_jin_stat, 3), round(x$auto_lag_jin_L, 3), round(x$auto_lag_jin_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statistic", "L hat", "p-value")
  }
  
  if ("Auto_Lag_Bartlett" %in% test){
    # Add boot Bartlett tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Bartlett", round(x$auto_lag_bart_stat, 3), round(x$auto_lag_bart_L, 3), round(x$auto_lag_bart_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statistic", "L hat", "p-value")
  }
  
  test_string <- ifelse(x$bootstrap_fixed_stats == TRUE, "\n Boot Fixed Lag Tests:\n \n", "\n Fixed Lag Tests:\n \n")
  weighted_test_string <- ifelse(x$bootstrap_fixed_stats == TRUE, "\n Boot Weighted Fixed Lag Tests:\n \n", "\n Weighted Fixed Lag Tests:\n \n")
  
  if ("Fixed_Lag" %in% test ){
    cat(test_string)
    print(out_table, row.names = FALSE, ...)
  }
  if ("Weighted_Fixed_Lag" %in% test ){
    cat(weighted_test_string)
    print(out_table_weight, row.names = FALSE, ...)
  }
  if ("Auto_Lag_Jin" %in% test | "Auto_Lag_Bartlett" %in% test ){
    cat("\n Automatic Lag Selection Tests:\n \n")
    print(out_table_boot, row.names = FALSE, ...)
  }
  
}

#' Test for equality of autocovariance functions for two (linearly Weighted_Fixed_Lag) stationary time series
#'
#' @description Perform a hypothesis test for equality of autocovariance functions for two time series with one or more of the following methods: (Weighted) Fixed_Lag, (Weighted) Weighted_Fixed_Lag, Bootstrapped Weighted_Fixed_Lag and Bootstrapped Bartlett. The former two tests assess equality up to a fixed lag, while the latter two select the "optimal lag" for testing the hypothesis using an AIC-like penalty at each lag. The tests can handle multivariate time series, but the computations become considerably more intense with an increase in dimension. 
#'
#' @param X a \eqn{n x m} mean zero (column-wise) stationary time series with \eqn{m < n}. Must be a matrix.
#' @param Y a \eqn{n x m} mean zero (column-wise) stationary time series with \eqn{m < n}. Must be a matrix.
#' @param max_lag the maximum lag to be considered. Must be a positive integer less than \eqn{n}. Note that the automatic lag selection tests may choose \code{max_lag} less than the one supplied. If not supplied, \code{max_lag = ceiling((log2(n))^0.9999)}.
#' @param test the tests to be performed. Must be a vector containing a subset of \code{"Fixed_Lag"}, \code{"Weighted_Fixed_Lag"}, \code{"Auto_Lag_Jin"} and \code{"Auto_Lag_Bartlett"}.
#' @param trunc for the \code{"Fixed_Lag"}, \code{"Weighted_Fixed_Lag"} and \code{"Auto_Lag_Bartlett"} tests, the truncation rule used in Bartlett's formula. If not supplied, \code{trunc = floor(n^(1/3))}.
#' @param num_bootstrap for the \code{"Auto_Lag_Jin"} and \code{"Auto_Lag_Bartlett"} tests or the \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests with \code{bootstrap_fixed_stats = TRUE}, the number of bootstrap resamples to be used in the Block of Blocks algorithm. Must be a positive integer.
#' @param block_size for the \code{"Auto_Lag_Jin"} and \code{"Auto_Lag_Bartlett"} tests or the \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests with \code{bootstrap_fixed_stats = TRUE}, the block length to be used in the Block of Blocks algorithm. Must be a positive integer less than \eqn{n}. If not supplied, \code{block_size = max(floor(0.5 * n^(1/3)), 2)}.
#' @param prewhiten for the \code{"Auto_Lag_Jin"} and \code{"Auto_Lag_Bartlett"} tests or the \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests with \code{bootstrap_fixed_stats = TRUE}, should the supplied time series be prewhitened? \code{prewhiten = TRUE} is strongly recommended.
#' @param dependent_series for the \code{"Fixed_Lag"}, \code{"Weighted_Fixed_Lag"} and \code{"Auto_Lag_Bartlett"} tests, should the series be assumed dependent?
#' @param bootstrap_fixed_stats for the \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests, should the test null distribution of the statistics be bootstrapped?
#' @param weights for the \code{"Weighted_Fixed_Lag"} test, the weight placed on each lag. Should be a vector of length \code{max_lag + 1} with values between \eqn{0} and \eqn{1}. By default, weights are \code{1, max_lag/(max_lag + 1), \dots, 1/(max_lag + 1)}.
#' @param average_bartlett_cov for the \code{"Auto_Lag_Bartlett"} test, should the covariance be estimated on the original sample or be averaged over each bootstrap resample?
#'
#' @return A named list containing many relevant statistics pertaining to each test
#' \itemize{
#'  \item \code{delta} - the sample difference in autocovariance functions.
#'  \item \code{dep_cov} - the asymptotic covariance matrix of the autocovariance differences computed assuming dependence of the series.
#'  \item \code{ind_cov} - - the asymptotic covariance matrix of the autocovariance differences computed assuming independence of the series.
#'  \item \code{fixed_lag_stat} - for the \code{"Fixed_Lag"} test, the test statistic.
#'  \item \code{fixed_lag_df} - for the \code{"Fixed_Lag"} test, the degrees of freedom associated with the test statistic.
#'  \item \code{fixed_lag_weight_stat} - for the \code{"Weighted_Fixed_Lag"} test, the weighted test statistic.
#'  \item \code{fixed_lag_weight_alpha} - for the \code{"Weighted_Fixed_Lag"} test, the scale parameter associated with the weighted test statistic.
#'  \item \code{fixed_lag_weight_beta} - for the \code{"Weighted_Fixed_Lag"} test, the shape parameter associated with the weighted test statistic.
#'  \item \code{fixed_lag_pval} - for the \code{"Fixed_Lag"} test, the associated p-value.
#'  \item \code{fixed_lag_weight_pval} - for the \code{"Weighted_Fixed_Lag"} test, the associated p-value of the weighted test.
#'  \item \code{auto_lag_jin_stat} - for the \code{"Auto_Lag_Jin"} test, the test statistic.
#'  \item \code{auto_lag_bart_stat} - for the \code{"Auto_Lag_Bartlett"} test, the test statistic.
#'  \item \code{auto_lag_jin_L} - for the \code{"Auto_Lag_Jin"} test, the automatically selected lag.
#'  \item \code{auto_lag_bart_L} - for the \code{"Auto_Lag_Bartlett"} test, the automatically selected lag.
#'  \item \code{auto_lag_jin_cov} - for the \code{"Auto_Lag_Jin"} test, the bootstrapped asymptotic covariance matrix.
#'  \item \code{auto_lag_bart_cov} - for the \code{"Auto_Lag_Bartlett"} test, the bootstrapped asymptotic covariance matrix.
#'  \item \code{auto_lag_jin_pval} - for the \code{"Auto_Lag_Jin"} test, the pseudo p-value.
#'  \item \code{auto_lag_bart_pval} - for the \code{"Auto_Lag_Bartlett"} test, the pseudo p-value.
#'  \item \code{dependent_series} - the provided \code{dependent_series} argument.
#'  \item \code{bootstrap_fixed_stats} - the provided \code{bootstrap_fixed_stats} argument.
#'  \item \code{max_lag} - the provided \code{max_lag} argument.
#'  \item \code{k} - the time series dimension.
#'  \item \code{n} - the time series length.
#'  \item \code{var_names} - the shared column names of \code{X} and \code{Y}.
#' }
#' @export 
#'
#' @details Consider two \eqn{m}-dimensional, stationary time series with mean zero. \code{autocovariance_test} tests for equality of autocovariance functions of the respective series. The independent \code{"Fixed_Lag"} test is given by Lund et. al (2009). Their test assumes
#' independence of the two series and quantifies the distribution of the autocovariance differences through Bartlett's formula (see Brockwell and Davids (1991)). \code{autocovariance_test} provides the corresponding asymptotic covariance matrix under the null, \code{ind_cov}, 
#' the chi-square test statistic, \code{fixed_lag_stat} and its associated degrees of freedom and p-value, \code{fixed_lag_df} and \code{fixed_lag_pval}. An extension of the \code{"Fixed_Lag"} test to linearly dependent series can be run by setting \code{dependent_series = TRUE}. 
#' If the independence of the series is not obvious, the \code{dependent_series = TRUE} test can be used without penalty. Both tests assess equality of autocovariances up to a fixed lag \code{max_lag}. The \code{trunc} parameter controls the truncation of the Bartlett formula infinite 
#' sum (see the vignette for more details). Any truncation rule must follow the assumptions of Theorem A.1 in Berkes et al. (2006).
#' 
#' \code{autocovariance_test} also computes weighted variants of the \code{"Fixed_Lag"} tests using \code{test = "Weighted_Fixed_Lag"}. The weighting idea is borrowed from goodness-of-fit testing (see Fisher and Gallagher (2012)). By default, we assign to lags \eqn{0, 1, \dots , L} weights
#' \eqn{1, max_lag/(max_lag + 1), \dots , 1/(max_lag + 1)}. This can be controlled via the \code{weights} argument. The corresponding test statistics have distributions that can be approximated by a Gamma distribution. For the \code{"Weighted_Fixed_Lag"} test, \code{autocovariance_test} outputs the test statistic, \code{fixed_lag_weight_stat}, its gamma distribtion
#' parameters, \code{weighted_fixed_lag_alpha} and \code{weighted_fixed_lag_beta}, and the associated p-value \code{weighted_fixed_lag_pval}. All previously mentioned tests require Gaussianity of the series. The following two tests can be applied to series with non-normal distributions. 
#' 
#' The \code{"Auto_Lag_Jin"} (see Jin et. al (2019)) and \code{"Auto_Lag_Bartlett"} tests compute a similar test statistic to the \code{"Fixed_Lag"} test but through a block of blocks bootstrap. The original test statistic is monotone increasing with lag, so the tests apply
#' an AIC-like penalty to choose the "best" lag (up to \code{max_lag}) to perform the hypothesis test with. This provides some desirable properties under the alternative hypothesis. \code{prewhiten = TRUE} fits a large VAR(p) process to \code{X} and \code{Y} and then applies 
#' the tests to the residuals of said process. It is highly recommended as prewhitening improves empirical rejection rates. \code{num_bootstrap} and \code{block_size} control the number of bootstrap resamples and block length used in the \code{"Auto_Lag_Jin"} and \code{"Auto_Lag_Bartlett"} tests. \code{"Auto_Lag_Bartlett"} replaces the Jin 
#' et. al (2019) covariance estimate with the \code{"Fixed_Lag"} estimate averaged over each bootstrap resample. For the \code{"Auto_Lag_Jin"} test, \code{autocovariance_test} outputs the test statistic, \code{auto_lag_jin_stat}, the optimally selected lag, \code{auto_lag_jin_L}, the 
#' estimated covariance under the null, \code{auto_lag_jin_cov}, and psuedo p-value \code{auto_lag_jin_pval}. The same statistics for \code{"Auto_Lag_Bartlett"} are given by \code{auto_lag_bart_stat}, \code{auto_lag_bart_L}, \code{auto_lag_bart_cov} and \code{auto_lag_bart_pval}. Alternatively, one may choose to bootstrap the
#' \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests by setting \code{bootstrap_fixed_stats = TRUE}. This uses the same covariance estimate as the \code{"Fixed_Lag"} and \code{"Weighted_Fixed_Lag"} tests, but approximates the null distribution via the Block of Blocks bootstrap.
#' 
#' For more on key assumptions behind each of the tests and a thorough example, please see the vignette. 
#' 
#' @references Berkes, I., Horvath, L., Kokoszka, P. and Shao, Q. M. (2006). \emph{On discriminating between longrange dependence and changes in the mean. Annals of Statistics, 34, 1140–65}\cr
#' \cr Brockwell, P. J., and Davis, R. A. (1991). \emph{Time series: theory and methods (Second ed.). New York: Springer-Verlag. }\cr 
#' \cr Fisher, T. J., and Gallagher, C. M. (2012).  \emph{New weighted portmanteau statistics for time series goodness  of  fit  testing, J.  Amer.  Statist.  Assoc., 107(498),  777–787}\cr 
#' \cr Jin, L., Cai, L., and Wang, S. (2019).  \emph{A computational bootstrap procedure to compare two Weighted_Fixed_Lag  time  series, Journal  of  Statistical  Computation  and  Simulation, 89(15), 2831–2847}\cr 
#' \cr Lund, R., Bassily, H., and Vidakovic, B. (2009). \emph{Testing equality of stationary autocovariances, Journal of Time Series Analysis, 30(3), 332–348}\cr
#'  
#'
#' @examples 
#' ## Simulate two marignal AR(1) processes that are correlated at lag zero and run tests
#' set.seed(12345)
#' 
#' # Simulate correlated errors
#' Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2) # covariance of errors
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X <- matrix(arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1])) 
#' Y <- matrix(arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 2]))
#' 
#' # Run tests for equality of autocovariance functions up to lag 5
#' output <- autocovariance_test(X, Y, max_lag = 5, 
#' test = c("Fixed_Lag","Weighted_Fixed_Lag","Auto_Lag_Jin","Auto_Lag_Bartlett"))
#' # All tests fail to reject the null
#' output
#' 
#' ## In males and females, do monthly lung deaths in the UK from 1974-1979
#' ## have the same autocovariance structure?
#' ## Here we carefully go through the process of making sure the series 
#' ## adhere to the test assumptions
#' ?UKLungDeaths 
#' 
#' n <- length(mdeaths) # series length
#' # plot series
#' plot(mdeaths)
#' plot(fdeaths)
#' # Obvious non-zero means and not stationary 
#' 
#' # Standardize by monthly means and standard deviations to ensure stationarity
#' male_matrix <- matrix(mdeaths, 6, 12, byrow = TRUE)
#' female_matrix <- matrix(fdeaths, 6, 12, byrow = TRUE)
#' # standardize series
#' male_stnd <- c(t(apply(male_matrix, 2, scale)))
#' female_stnd <- c(t(apply(female_matrix, 2, scale)))
#' 
#' # Plot again to check stationarity
#' t <- 1:n
#' plot(t, male_stnd, type = "l") # definitely has a trend
#' plot(t, female_stnd, type = "l") # maybe a slight trend
#' 
#' # Detrend both series to adhere to mean zero assumption + stationarity
#' male_stnd <- lm(male_stnd ~ t)$residuals
#' female_stnd <- lm(female_stnd ~ t)$residuals
#' plot(t, male_stnd, type = "l") # reasonably stationary, 
#' plot(t, female_stnd, type = "l") # reasonably stationary
#' 
#' # Both series are reasonably Gaussian
#' qqnorm(male_stnd)
#' qqline(male_stnd)
#' qqnorm(female_stnd)
#' qqline(female_stnd)
#'
#' # Check for correlation
#' ccf(female_stnd, male_stnd) # Clear lag zero correlation
#' 
#' # The series are likely gaussian and Weighted_Fixed_Lag, hence it 
#' # makes most sense to use "Weighted_Fixed_Lag" test
#' output2 <- autocovariance_test(matrix(male_stnd), matrix(female_stnd), 5, 
#' test = "Weighted_Fixed_Lag")
#' output2
#' 
#' # Both the weighted and unweighted tests fail to reject the null hypothesis
#' # The weighted test comes close to rejecting, which makes sense given the large difference in lag 0
#' acf(cbind(male_stnd, female_stnd), type = "covariance", lag.max = 5)
#' 
autocovariance_test <- function(X, Y, max_lag = NULL, test = "Auto_Lag_Bartlett", trunc = NULL, num_bootstrap = 500, block_size = NULL, prewhiten = TRUE, dependent_series = TRUE, bootstrap_fixed_stats = FALSE, weights = NULL, average_bartlett_cov = TRUE){
  
  # Compatibility Tests and reassign L if need be
  Lb <- compatibilityChecks(X, Y, max_lag, test, trunc, num_bootstrap, block_size, prewhiten, dependent_series, bootstrap_fixed_stats, weights, average_bartlett_cov)
  L <- Lb[[1]]
  b <- Lb[[2]]
  trunc <- Lb[[3]]
  X <- Lb[[4]]
  Y <- Lb[[5]]
  prewhiten <- Lb[[6]]
  B <- Lb[[7]]
  dependent_series <- Lb[[8]]
  bootstrap_fixed_stats <- Lb[[9]]
  weights <- Lb[[10]]
  average_bartlett_cov <- Lb[[11]]
  
  # Get some info
  n <- nrow(X)
  k <- ncol(Y)
  
  # Initialize output list
  out <- list()
  
  # Compute Bartlett covariance
  deltaCovar <- calculateCovariance(X, Y, L, trunc)
  out <- c(out, deltaCovar)
  
  # Compute Independent Fixed Lag test
  if ("Fixed_Lag" %in% test & dependent_series == FALSE & bootstrap_fixed_stats == FALSE){
    indTest <- calculateTestStat(deltaCovar$delta, deltaCovar$ind_cov, n, L, k, weights)
    # Compute p-values
    indTest$pval <- stats::pchisq(indTest$stat, df = indTest$df, lower.tail = FALSE)
    names(indTest) <- paste0("fixed_lag_",names(indTest))
    out <- c(out, indTest[c(1, 2, 6)])
  }
  # Compute Independent Weighted Fix Lag test
  if ("Weighted_Fixed_Lag" %in% test & dependent_series == FALSE & bootstrap_fixed_stats == FALSE){
    indWeightTest <- calculateTestStat(deltaCovar$delta, deltaCovar$ind_cov, n, L, k, weights)
    indWeightTest$weight_pval <- stats::pgamma(indWeightTest$weight_stat, shape = indWeightTest$weight_alpha, scale = indWeightTest$weight_beta, lower.tail = FALSE)
    names(indWeightTest) <- paste0("fixed_lag_",names(indWeightTest))
    out <- c(out, indWeightTest[3:6])
  }
  # Compute Dependent Fixed Lag test
  if ("Fixed_Lag" %in% test & dependent_series == TRUE & bootstrap_fixed_stats == FALSE){
    depTest <- calculateTestStat(deltaCovar$delta, deltaCovar$dep_cov, n, L, k, weights)
    # Compute p-values
    depTest$pval <- stats::pchisq(depTest$stat, df = depTest$df, lower.tail = FALSE)
    names(depTest) <- paste0("fixed_lag_",names(depTest))
    out <- c(out, depTest[c(1, 2, 6)])
  }
  # Compute Dependent Weighted Fix Lag test
  if ("Weighted_Fixed_Lag" %in% test & dependent_series == TRUE & bootstrap_fixed_stats == FALSE){
    depWeightTest <- calculateTestStat(deltaCovar$delta, deltaCovar$dep_cov, n, L, k, weights)
    depWeightTest$weight_pval <- stats::pgamma(depWeightTest$weight_stat, shape = depWeightTest$weight_alpha, scale = depWeightTest$weight_beta, lower.tail = FALSE)
    names(depWeightTest) <- paste0("fixed_lag_",names(depWeightTest))
    out <- c(out, depWeightTest[3:6])
  }
  
  
  # Fixed Lag Bootstrapped tests 
  if ("Fixed_Lag" %in% test & "Weighted_Fixed_Lag" %in% test & bootstrap_fixed_stats == TRUE){
    fixedBoot <- calculateBootTestStatFixed(X, Y, L, B, b, prewhiten, trunc, dependent_series, weights)
    out <- c(out, fixedBoot)
  } else {
    
    if ("Fixed_Lag" %in% test & bootstrap_fixed_stats == TRUE){
      fixedUnweightedBoot <- calculateBootTestStatFixedUnweighted(X, Y, L, B, b, prewhiten, trunc, dependent_series, weights)
      out <- c(out, fixedUnweightedBoot)
    }
    
    if ("Weighted_Fixed_Lag" %in% test & bootstrap_fixed_stats == TRUE){
      fixedWeightedBoot <- calculateBootTestStatFixedWeighted(X, Y, L, B, b, prewhiten, trunc, dependent_series, weights)
      out <- c(out, fixedWeightedBoot)
    }
  }
  
  # Automatic Lag Selection tests
  if ("Auto_Lag_Jin" %in% test & "Auto_Lag_Bartlett" %in% test){
    bootTest <- calculateBootTestStat(X, Y, L, B, b, prewhiten, trunc, dependent_series, average_bartlett_cov)
    out <- c(out, bootTest)
  } else {
    # Compute bootstrapped tests if need be
    if ("Auto_Lag_Jin" %in% test){
      bootTest <- calculateBootTestStatJin(X, Y, L, B, b, prewhiten, trunc)
      out <- c(out, bootTest)
    }
    
    if ("Auto_Lag_Bartlett" %in% test){
      bootTest <- calculateBootTestStatBartlett(X, Y, L, B, b, prewhiten, trunc, dependent_series, average_bartlett_cov)
      out <- c(out, bootTest)
    }
  }
  
  out$dependent_series <- dependent_series
  out$bootstrap_fixed_stats <- bootstrap_fixed_stats
  out$max_lag <- L
  out$k <- k
  out$n <- n
  
  if (is.null(colnames(X))){
    out$var_names <- colnames(Y)
  } else {
    out$var_names <- colnames(X)
  }
  
  class(out) <- "acvfTest"

  return(out)
  
}

#' Monthly High Temperatures in New York and Cleveland from 1960 to 2019
#' 
#' @description Monthly high temperatures in New York City and Cleveland from 1960 to 2019. The stations that record the data are located at Cleveland Hopkins International Airport (ICAO: KCLE; GHCND: USW00014820)
#' and John F. Kennedy International Airport (ICAO: KJFK; GHCND: USW00094789).
#'
#' @docType data
#'
#' @usage data(cityTemps)
#'
#' @format An object of class \code{"data.frame"} with 1440 rows and 4 columns.
#' \describe{
#'   \item{City}{the city, Cleveland or New York.}
#'   \item{Year}{the year of the high temperature.}
#'   \item{Month }{the month of the high temperature.}
#'   \item{TMAX}{the high temperature (degrees Farenheit).}
#' }
#'
#' @details Aggregated Monthly high temperatures at the KCLE and KJFK airports. Daily high temperatures were collected from the National Oceanic and Atmospheric Administration's National Centers for Environmental Information database, the Global Historical Climatological Network (GHCND IDs: USW00014820 and USW00094789).
#'
#' @references Menne,  M.,  Durre,  I.,  Korzeniewski,  B.,  McNeal,  S.,  Thomas,  K.,  Yin,  X., Anthony, S. and  Ray,  R. (2012). \emph{Global  historical  climatology  network  -  daily  (ghcn-daily),  version  3.27}\cr
#' \cr Menne, M., Durre, I., Vose, R., Gleason, B., and Houston, T. (2012). \emph{An overview of the global historical  climatology  network-daily  database. Journal  of  Atmospheric  and  Oceanic Technology, 29, 897-910.}\cr
#'
#' @source https://www.ncdc.noaa.gov/cdo-web/datatools/findstation
#'
"cityTemps"

#' Concentration of Hourly Air Pollutants in two London Locations from 2015 to 2018
#' 
#' @description Concentration of daily average air pollutants in the North Kensington and Marylebone Road stations from 2015 to 2018. Data is taken from the Air Information Resource (UK-AIR) of the Department for Environment, Food and Rural Affairs in the United Kingdom (see https://uk-air.defra.gov.uk/).
#'
#' @docType data
#'
#' @usage data(londonAirQuality)
#'
#' @format An object of class \code{"data.frame"} with 140256 rows and 4 columns.
#' \describe{
#'   \item{Site}{the location of measurement, North Kensington or Marylebone Road.}
#'   \item{Parameter}{the pollutant, Nitric Oxide or Nitrogen Dioxide}
#'   \item{Datetime}{the date and time of the measurement}
#'   \item{Value}{concentration of pollutant (in ug/m3)}
#' }
#' 
#' @details The Air Information Resource (UK-AIR) data was accessed with the \code{rdefra} package in R.
#' 
#' @references Vitolo, Claudia, Russell, Andrew, Tucker, Allan (2016). “rdefra: Interact with the UK AIR Pollution Database from DEFRA.” \emph{The Journal of Open Source Software}, 1(4). doi: 10.21105/joss.00051, http://dx.doi.org/10.21105/joss.00051.
#'
#' @source https://uk-air.defra.gov.uk/
#' 
"londonAirQuality"