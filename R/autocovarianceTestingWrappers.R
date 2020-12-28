## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#' @useDynLib autocovarianceTesting, .registration=TRUE
#' 

# Function to run compatibility checks
compatibilityChecks <- function(X, Y, L, test, trunc, B, b, prewhiten, plot){
  
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
  p <- ifelse(prewhiten == TRUE & ("bootDependent" %in% test | "bootBartlett" %in% test), ceiling(min(10 * log(n), 0.8 * sqrt(n))), 0)
  
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
  if ( !is.numeric(b) ){
    stop(paste("b must be an integer and larger than 2"))
  }
  
  # b must be larger than 2
  if ( !((b == round(b)) & (b >= 2)) ){
    stop(paste("b must be an integer and larger than 2"))
  }
  
  # b must be less than n
  if ( b >= (n - p - 1) ){
    stop(paste("b must be strictly less than n the length of the series (minus AR order + 1 if prewhitening)"))
  }
  
  # Dimension of time series must be smaller than n
  if ( ncol(X) >= (n - p - 1)){
    stop(paste("The dimension of the series must be smaller than its length (minus AR order + 1 if prewhitening)"))
  }
  
  # test must be one of c(Independent, Dependent, bootDependent, bootBartlett)
  if (!(all(test %in% c("Independent", "Dependent", "bootDependent", "bootBartlett")))){
    stop(paste("test must be one of Independent, Dependent, bootDependent, or bootBartlett"))
  }
  
  # B must be numeric
  if ( !is.numeric(B)){
    stop(paste("B must be an integer and larger than 1"))
  }
  
  # B must be larger than 1
  if (!((B == round(B)) & (B >= 1)) ){
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
  if ( !is.logical(prewhiten) ){
    stop(paste("prewhiten must be logical"))
  }
  
  # plot must be boolean
  if ( !is.logical(plot) ){
    stop(paste("plot must be logical"))
  }
  
  # b Warning
  if ( b < 3 ){
    warning("b is small (b < 3), bootstrap resamples may not be representative")
  }
  
  # B warning
  if ( B < 100 ){
    warning("B is small (B < 100), results may not be reliable")
  }
  
  return(list(L, b, trunc, X, Y))
}


#' Summarizing autocovarianceTest function output
#'
#' @param x a \code{"acvfTest"} object given by \code{"autocovarianceTest()"}.
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
#' output <- autocovarianceTest(X, Y, L = 5, 
#' test = c("Independent","Dependent","bootDependent","bootBartlett"))
#' # All tests fail to reject the null
#' print(output)
print.acvfTest <- function(x, ...){
  # Create test vector
  test <- c()
  listnames <- names(x)
  test <- if("ind_stat" %in% listnames){ c(test, "Independent") }else{ test }
  test <- if("dep_stat" %in% listnames){ c(test, "Dependent") }else{ test }
  test <- if("jin_stat" %in% listnames){ c(test, "bootDependent") }else{ test }
  test <- if("bart_stat" %in% listnames){ c(test, "bootBartlett") }else{ test }
  
  if( "Independent" %in% test | "Dependent" %in% test ){
    # Initialize Fixed Lag Table
    out_table <- data.frame(0, 1, 2, 3)
    out_table <- out_table[-1, ]
    
    # Initialize Weighted Table
    out_table_weight <- data.frame(0, 1, 2, 3, 4)
    out_table_weight <- out_table_weight[-1, ]
  }
  
  if ("Independent" %in% test){
    # Add independent tests
    out_table <- rbind(out_table, c("Independent", round(x$ind_stat, 3), x$ind_df, round(x$ind_pval, 3)))
    colnames(out_table) <- c("Test", "Chi-Sq", "df", "p-value")
    out_table_weight <- rbind(out_table_weight, c("Weighted Independent", round(x$ind_weight_stat, 3), round(x$ind_alpha, 3), round(x$ind_beta, 3), round(x$ind_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Gamma", "alpha", "beta", "p-value")
  }
  
  if ("Dependent" %in% test){
    # Add dependent tests
    out_table <- rbind(out_table, c("Dependent", round(x$dep_stat, 3), x$dep_df, round(x$dep_pval, 3)))
    colnames(out_table) <- c("Test", "Chi-Sq", "df", "p-value")
    out_table_weight <- rbind(out_table_weight, c("Weighted Dependent", round(x$dep_weight_stat, 3), round(x$dep_alpha, 3), round(x$dep_beta, 3), round(x$dep_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Gamma", "alpha", "beta", "p-value")
  }
  
  # Initialize bootstrapped tables
  if ("bootDependent" %in% test | "bootBartlett" %in% test){
    out_table_boot <- data.frame(0, 1, 2, 3)
    out_table_boot <- out_table_boot[-1, ]
  }
  
  if ("bootDependent" %in% test){
    # Add jin tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Jin", round(x$jin_stat, 3), round(x$jin_L, 3), round(x$jin_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statitic", "L hat", "p-value")
  }
  
  if ("bootBartlett" %in% test){
    # Add boot Bartlett tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Bartlett", round(x$bart_stat, 3), round(x$bart_L, 3), round(x$bart_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statitic", "L hat", "p-value")
  }
  
  if ("Independent" %in% test | "Dependent" %in% test ){
    cat("\n Fixed Lag Tests:\n \n")
    print(out_table, row.names = FALSE, ...)
    cat("\n Weighted Fixed Lag Tests:\n \n")
    print(out_table_weight, row.names = FALSE, ...)
  }
  if ("bootDependent" %in% test | "bootBartlett" %in% test ){
    cat("\n Automatic Lag Selection Tests:\n \n")
    print(out_table_boot, row.names = FALSE, ...)
  }
  
}

# Plot helper function
barplotACVF <- function(dataset, ymin, ymax, k, time, L){
  if (k == 1){
    graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = c("#000033","#800000"), beside = TRUE, ylim = c(ymin, ymax), data = dataset)
    graphics::lines(x = c(0, 3.5 * L), y = c(0, 0))
  } else {
    if (!(time %in% 1:k) & !(time %% k == 0)){
      graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = c("#000033","#800000"), beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", yaxt = "n", data = dataset)
      graphics::lines(x = c(0, 3.5 * L), y = c(0, 0))
    } else {
      if (time %in% 1:k & !(time %% k == 0)){
        graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = c("#000033","#800000"), beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", data = dataset)
        graphics::lines(x = c(0, 3.5 * L), y = c(0, 0))
      } else {
        if (!(time %in% 1:k) & time %% k == 0){
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = c("#000033","#800000"), beside = TRUE, ylim = c(ymin, ymax), yaxt = "n", data = dataset)
          graphics::lines(x = c(0, 3.5 * L), y = c(0, 0))
        } else {
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = c("#000033","#800000"), beside = TRUE, ylim = c(ymin, ymax), data = dataset)
          graphics::lines(x = c(0, 3.5 * L), y = c(0, 0))
        }
      }
    }
  }
}

# Function to create a autocovariance difference plot
acvfPlot <- function(X, Y, L){
  # Get time series dimension
  k <- ncol(X)
  
  # Get autocovariance up to lag 6
  acvf <- stats::acf(cbind(X, Y), lag.max = L, type = "covariance", plot = FALSE)$acf
  
  # Shape data for plotting
  acvfX <- data.frame(acvf = c(acvf[ , 1:k, 1:k]), lags = 0:L, dim1 = rep(paste0("Dim", 1:k), each = (L + 1) * k ), dim2 = rep(rep(paste0("Dim", 1:k), each = L + 1), k), timeseries = "X")
  acvfY <- data.frame(acvf = c(acvf[ , (k + 1):(2 * k), (k + 1):(2 * k)]), lags = 0:L, dim1 = rep(paste0("Dim", 1:k), each = (L + 1) * k ), dim2 = rep(rep(paste0("Dim", 1:k), each = L + 1), k), timeseries = "Y")
  plot_data <- rbind(acvfX, acvfY)
  plot_data$lags <- factor(plot_data$lags)
  
  # Y axis limits
  ymax <- max(plot_data$acvf)
  ymin <- min(plot_data$acvf)
  
  # Split data by dimensions
  data_list <- split(plot_data, paste(plot_data$dim1, plot_data$dim2))
  
  # Adjust margins and such
  graphics::par(mfcol = c(k, k), oma = c(3, 3, 2, 2), mar = c(1, 1, 0, 0), mgp = c(1, 1, 0), xpd = NA)
  # Plot 
  for (i in 1:(k^2)){
    barplotACVF(data_list[[i]], ymin, ymax, k, i, L)
  }
  graphics::par(mfrow = c(1, 1))
  
  # print the overall labels
  graphics::mtext('Lag', side = 1, outer = TRUE, line = 2)
  graphics::mtext('ACVF', side = 2, outer = TRUE, line = 2)
  at <- seq(1 / (2 * k), 1 - 1 / (2 * k), by = 1 / k)
  graphics::mtext(paste0(1:k), side = 3, outer = TRUE, line = 0, at = at)
  graphics::mtext(paste0(1:k), side = 4, outer = TRUE, line = 0, at = at)
  
  # add legend
  graphics::legend("topright", legend = c("X", "Y"), fill = c("#000033","#800000"), cex = 1/k)
  # Reset graphics
  graphics::par(oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), xpd = FALSE)
  
}

#' Test for equality of autocovariance functions for two (linearly dependent) stationary time series
#'
#' @description Perform a hypothesis test for equality of autocovariance functions for two time series with one or more of the following methods: (Weighted) Independent, (Weighted) Dependent, Bootstrapped Dependent and Bootstrapped Bartlett. The former two tests assess equality up to a fixed lag, while the latter two select the "optimal lag" for testing the hypothesis using an AIC-like penalty at each lag. The tests can handle multivariate time series, but the computations become considerably more intense with an increase in dimension. 
#'
#' @param X a \eqn{n x m} mean zero (column-wise) stationary time series with \eqn{m < n}. Must be a matrix.
#' @param Y a \eqn{n x m} mean zero (column-wise) stationary time series with \eqn{m < n}. Must be a matrix.
#' @param L the maximum lag to be considered. Must be a positive integer less than \eqn{n}. Note that the automatic lag selection tests may choose \eqn{L} less than the one supplied. If not supplied, \code{L = ceiling((log2(n))^0.9999)}.
#' @param test the tests to be performed. Must be a vector containing a subset of \code{"Independent"}, \code{"Dependent"}, \code{"bootDependent"} and \code{"bootBartlett"}.
#' @param trunc for the \code{"Independent"}, \code{"Dependent"} and \code{"bootBartlett"} tests, the truncation rule used in Bartlett's formula. If not supplied, \code{trunc = floor(n^(1/3))}.
#' @param B for the \code{"bootDependent"} and \code{"bootBartlett"} tests, the number of bootstrap resamples to be used in the Block of Blocks algorithm. Must be a positive integer.
#' @param b for the \code{"bootDependent"} and \code{"bootBartlett"} tests, the block length to be used in the Block of Blocks algorithm. Must be a positive integer less than \eqn{n}. If not supplied, \code{b = max(floor(0.5 * n^(1/3)), 2)}.
#' @param prewhiten for the \code{"bootDependent"} and \code{"bootBartlett"} tests, should the supplied time series be prewhitened? \code{prewhiten = TRUE} is strongly recommended.
#' @param plot should a plot of the tested sample autocovariances be made?
#'
#' @return A named list containing many relevant statistics pertaining to each test
#' \itemize{
#'  \item \code{delta} - for the \code{"Independent"} and \code{"Dependent"} tests, the sample difference in autocovariance functions.
#'  \item \code{dep_cov} - for the \code{"Dependent"} test, the asymptotic covariance matrix of the autocovariance differences.
#'  \item \code{ind_cov} - for the \code{"Independent"} test, the asymptotic covariance matrix of the autocovariance differences.
#'  \item \code{ind_stat} - for the \code{"Independent"} test, the test statistic.
#'  \item \code{ind_df} - for the \code{"Independent"} test, the degrees of freedom associated with the test statistic.
#'  \item \code{ind_weight_stat} - for the \code{"Independent"} test, the weighted test statistic.
#'  \item \code{ind_alpha} - for the \code{"Independent"} test, the scale parameter associated with the weighted test statistic.
#'  \item \code{ind_beta} - for the \code{"Independent"} test, the shape parameter associated with the weighted test statistic.
#'  \item \code{ind_pval} - for the \code{"Independent"} test, the associated p-value.
#'  \item \code{ind_weight_pval} - for the \code{"Independent"} test, the associated p-value of the weighted test.
#'  \item \code{dep_stat} - for the \code{"Dependent"} test, the test statistic.
#'  \item \code{dep_df} - for the \code{"Dependent"} test, the degrees of freedom associated with the test statistic.
#'  \item \code{dep_weight_stat} - for the \code{"Dependent"} test, the weighted test statistic.
#'  \item \code{dep_alpha} - for the \code{"Dependent"} test, the scale parameter associated with the weighted test statistic.
#'  \item \code{dep_beta} - for the \code{"Dependent"} test, the shape parameter associated with the weighted test statistic.
#'  \item \code{dep_pval} - for the \code{"Dependent"} test, the associated p-value.
#'  \item \code{dep_weight_pval} - for the \code{"Dependent"} test, the associated p-value of the weighted test.
#'  \item \code{jin_stat} - for the \code{"bootDependent"} test, the test statistic.
#'  \item \code{bart_stat} - for the \code{"bootBartlett"} test, the test statistic.
#'  \item \code{jin_L} - for the \code{"bootDependent"} test, the automatically selected lag.
#'  \item \code{bart_L} - for the \code{"bootBartlett"} test, the automatically selected lag.
#'  \item \code{jin_cov} - for the \code{"bootDependent"} test, the bootstrapped asymptotic covariance matrix.
#'  \item \code{bart_cov} - for the \code{"bootBartlett"} test, the bootstrapped asymptotic covariance matrix.
#'  \item \code{jin_pval} - for the \code{"bootDependent"} test, the pseudo p-value.
#'  \item \code{bart_pval} - for the \code{"bootBartlett"} test, the pseudo p-value. 
#' }
#' @export 
#'
#' @details Consider two \eqn{m}-dimensional, stationary time series with mean zero. \code{autocovarianceTest} tests for equality of autocovariance functions of the respective series. The \code{"Independent"} test is given by Lund et. al (2009). Their test assumes
#' independence of the two series and quantifies the distribution of the autocovariance differences through Bartlett's formula (see Brockwell and Davids (1991)). \code{autocovarianceTest} provides the corresponding asymptotic covariance matrix under the null, \code{dep_cov}, the
#' chi-square test statistic, \code{ind_stat} and its associated degrees of freedom and p-value, \code{ind_df} and \code{ind_pval}. The \code{"Dependent"} test is an extension of the \code{"Independent"} test to linearly dependent series. If the independence of the series is 
#' not obvious, the \code{"Dependent"} test can be used without penalty. Both tests assess equality of autocovariances up to a fixed lag \code{L}. The \code{trunc} parameter controls the truncation of the Bartlett formula infinite sum (see the vignette for more details). Any truncation rule 
#' must follow the assumptions of Theorem A.1 in Berkes et al. (2006).
#' 
#' When running either the \code{"Independent"} or \code{"Dependent"} test, \code{autocovarianceTest} also computes weighted variants. The weighting idea is borrowed from goodness-of-fit testing (see Fisher and Gallagher (2012)). We assign to lags \eqn{0, 1, \dots , L} weights
#' \eqn{1, L/(L + 1), \dots , 1/(L + 1)}. The corresponding test statistics have distributions that can be approximated by a Gamma distribution. For the \code{"Independent"} test, \code{autocovarianceTest} outputs the test statistic, \code{ind_weight_stat}, its gamma distribtion
#' parameters, \code{ind_alpha} and \code{ind_beta}, and the associated p-value \code{ind_weight_pval}. The same statistics for \code{"Dependent"} are given by \code{dep_weight_stat}, \code{dep_alpha}, \code{dep_beta} and \code{dep_weight_pval}. All previously mentioned tests
#' require Gaussianity of the series. The following two tests can be applied to series with non-normal distributions. 
#' 
#' The \code{"bootDependent"} (see Jin et. al (2019)) and \code{"bootBartlett"} tests compute a similar test statistic to the \code{"Dependent"} test but through a block of blocks bootstrap. The original test statistic is monotone increasing with lag, so the tests apply
#' an AIC-like penalty to choose the "best" lag (up to \code{L}) to perform the hypothesis test with. This provides some desirable properties under the alternative hypothesis. \code{prewhiten = TRUE} fits a large VAR(p) process to \code{X} and \code{Y} and then applies 
#' the tests to the residuals of said process. It is highly recommended as prewhitening improves empirical rejection rates. \code{B} and \code{b} control the number of bootstrap resamples and block length used in the \code{"bootDependent"} and \code{"bootBartlett"} tests. \code{"bootBartlett"} replaces the Jin 
#' et. al (2019) covariance estimate with the \code{"Dependent"} estimate averaged over each bootstrap resample. For the \code{"bootDependent"} test, \code{autocovarianceTest} outputs the test statistic, \code{jin_stat}, the optimally selected lag, \code{jin_L}, the 
#' estimated covariance under the null, \code{jin_cov}, and psuedo p-value \code{jin_pval}. The same statistics for \code{"bootBartlett"} are given by \code{bart_stat}, \code{bart_L}, \code{bart_cov} and \code{bart_pval}.
#' 
#' For more on key assumptions behind each of the tests and a thorough example, please see the vignette. 
#' 
#' @references Berkes, I., Horvath, L., Kokoszka, P. and Shao, Q. M. (2006). \emph{On discriminating between longrange dependence and changes in the mean. Annals of Statistics, 34, 1140–65}\cr
#' \cr Brockwell, P. J., and Davis, R. A. (1991). \emph{Time series: theory and methods (Second ed.). New York: Springer-Verlag. }\cr 
#' \cr Fisher, T. J., and Gallagher, C. M. (2012).  \emph{New weighted portmanteau statistics for time series goodness  of  fit  testing, J.  Amer.  Statist.  Assoc., 107(498),  777–787}\cr 
#' \cr Jin, L., Cai, L., and Wang, S. (2019).  \emph{A computational bootstrap procedure to compare two dependent  time  series, Journal  of  Statistical  Computation  and  Simulation, 89(15), 2831–2847}\cr 
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
#' output <- autocovarianceTest(X, Y, L = 5, 
#' test = c("Independent","Dependent","bootDependent","bootBartlett"))
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
#' # The series are likely gaussian and dependent, hence it makes most sense to use "Dependent" test
#' output2 <- autocovarianceTest(matrix(male_stnd), matrix(female_stnd), 5, test = "Dependent")
#' output2
#' 
#' # Both the weighted and unweighted tests fail to reject the null hypothesis
#' # The weighted test comes close to rejecting, which makes sense given the large difference in lag 0
#' acf(cbind(male_stnd, female_stnd), type = "covariance", lag.max = 5)
#' 
autocovarianceTest <- function(X, Y, L = NULL, test = "bootBartlett", trunc = NULL, B = 500, b = NULL, prewhiten = TRUE, plot = FALSE){
  
  # Compatibility Tests and reassign L if need be
  Lb <- compatibilityChecks(X, Y, L, test, trunc, B, b, prewhiten, plot)
  L <- Lb[[1]]
  b <- Lb[[2]]
  trunc <- Lb[[3]]
  X <- Lb[[4]]
  Y <- Lb[[5]]
  
  # Get some info
  n <- nrow(X)
  k <- ncol(Y)
  
  # Initialize output list
  out <- list()
  
  # Compute independent and dependent tests if need be
  if ("Independent" %in% test | "Dependent" %in% test){
    # Compute dependent covariance
    deltaCovar <- calculateCovariance(X, Y, L, trunc)
    out <- c(out, deltaCovar)
    # Compute test statistics
    if ("Independent" %in% test){
      indTest <- calculateTestStat(deltaCovar$delta, deltaCovar$ind_cov, n, L, k)
      # Compute p-values
      indTest$pval <- stats::pchisq(indTest$stat, df = indTest$df, lower.tail = FALSE)
      indTest$weight_pval <- stats::pgamma(indTest$weight_stat, shape = indTest$alpha, scale = indTest$beta, lower.tail = FALSE)
      names(indTest) <- paste0("ind_",names(indTest))
      out <- c(out, indTest)
    }
    if ("Dependent" %in% test){
      depTest <- calculateTestStat(deltaCovar$delta, deltaCovar$dep_cov, n, L, k)
      # Compute p-values
      depTest$pval <- stats::pchisq(depTest$stat, df = depTest$df, lower.tail = FALSE)
      depTest$weight_pval <- stats::pgamma(depTest$weight_stat, shape = depTest$alpha, scale = depTest$beta, lower.tail = FALSE)
      names(depTest) <- paste0("dep_",names(depTest))
      out <- c(out, depTest)
    }
  }
  
  if ("bootDependent" %in% test & "bootBartlett" %in% test){
    bootTest <- calculateBootTestStat(X, Y, L, B, b, prewhiten, trunc)
    out <- c(out, bootTest)
  } else {
    # Compute bootstrapped tests if need be
    if ("bootDependent" %in% test){
      bootTest <- calculateBootTestStatJin(X, Y, L, B, b, prewhiten, trunc)
      out <- c(out, bootTest)
    }
    
    if ("bootBartlett" %in% test){
      bootTest <- calculateBootTestStatBartlett(X, Y, L, B, b, prewhiten, trunc)
      out <- c(out, bootTest)
    }
  }
  
  # Plot 
  if( plot == TRUE ){
    acvfPlot(X, Y, L)
  }
  
  class(out) <- "acvfTest"

  return(out)
  
}

#' Monthly High Temperatures in New York and Cleveland from 1960 to 2019
#' 
#' @description Monthly high temperatures in New York City and Cleveland from 1960 to 2019. The stations that record the data are located at Cleveland Hopkins International Airport and John F. Kennedy International Airport.
#'
#' @docType data
#'
#' @usage data(cityTemps)
#'
#' @format An object of class \code{"data.frame"} with 1440 rows and 4 columns.
#' \describe{
#'   \item{City}{the city, CLE or NYC.}
#'   \item{Year}{the year of the high temperature.}
#'   \item{Month }{the month of the high temperature.}
#'   \item{TMAX}{the high temperature (degrees Farenheit).}
#' }
#'
#' @source https://www.ncdc.noaa.gov/cdo-web/datatools/findstation
#' 
"cityTemps"

#' Concentration of Daily Average Air Pollutants in two London Locations from 2015 to 2018
#' 
#' @description Concentration of daily average air pollutants in the London N. Kensington and London Marylebone Road stations from 2015 to 2018. Data is taken from the Air Information Resource (UK-AIR) of the Department for Environment, Food and Rural Affairs in the United Kingdom (see https://uk-air.defra.gov.uk/).
#'
#' @docType data
#'
#' @usage data(ukAir)
#'
#' @format An object of class \code{"data.frame"} with 5844 rows and 5 columns.
#' \describe{
#'   \item{Site}{the location of measurement, London N. Kensington or London Marylebone Road.}
#'   \item{Parameter}{the pollutant, Nitric Oxide or Nitrogen Dioxide}
#'   \item{Date}{the date of the measurement}
#'   \item{Value}{average daily concentration of pollutant (in ug/m3)}
#'   \item{logValue}{log of average daily concentration of pollutant (in log(ug/m3))}
#' }
#'
#' @source https://uk-air.defra.gov.uk/
#' 
"ukAir"