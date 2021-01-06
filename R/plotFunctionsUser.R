#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

################
## These are the functions the user will have access to
##


# This file has the background functions
#source("plotFunctionsBackground.R")



#' Construct a plot of the Autocovariance Function for two time series
#'
#' @description A convenience function that will construct a plot of the sample Autocovariance function for two (potentially multivariate) time series. The function will be a $k x k$ grid of side-by-side bar graphs showing the autocovariance value at each lag for each series.
#'
#' @param X a \eqn{n x m} time series with \eqn{m < n}. Must be a matrix, vector or ts object.
#' @param Y a \eqn{n x m} time series with \eqn{m < n}. Must be a matrix, vector or ts object.
#' @param max_lag the maximum lag to be considered. Must be a positive integer less than \eqn{n}.  If not supplied, \code{max_lag = ceiling((log2(n))^0.9999)}.
#' @param ggplot2 a flag on whether to construct a ggplot2 version of the plot.
#' @param plot_options a list containing elements to customize the plot, see details.
#'
#' @return If \code{ggplot2=FALSE}, an unnamed list containing base R graphic plot elements is supplied. If \code{ggplot2=TRUE}, a \code{gg} or \code{ggplot} object is returned.
#' 
#' @details The parameter \code{plot_options} is a list that contains three elements. Defaults are supplied. If the user attempts to change them error checking is performed and the user may receive a warning.
#' \itemize{
#'  \item \code{series_name} - a vector of two strings specifying names for the two series, defaults as "Series 1" and "Series 2"
#'  \item \code{bar_colors} - a vector of two strings specifying the colors for the bars of the autocovariance values, defaults as "gray10" and "gray50". The code also checks to make sure legitimate color choices are made.
#'  \item \code{var_names} - an optional vector of strings with length matching the dimensionality of the series specifying the names of the time series within the multivariate series. If not supplied, or incorrectly supplied, defaults to the convention Var1, Var2, ...
#'  }
#'  
#' @export 
#' 
#' @examples
#' library(ggplot2)
#' 
#' ## Simulate two marignal AR(1) processes that are correlated at lag zero and run tests
#' set.seed(12345)
#' 
#' # Simulate correlated errors
#' Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2) # covariance of errors
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1]) 
#' Y <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 2])
#' 
#' ## Base R Graphics
#' ## Default options
#' plot_autocovariances(X, Y)
#' 
#' ## Plot options
#' plot_autocovariances(X, Y, max_lag=7,
#'                      plot_options=list(series_names = c("My AR1", "Your AR1"),
#'                                        bar_colors = c("forestgreen", "gray20")))
#'                                        
#' ## ggplot2 Graphics
#' ## Default options
#' plot_autocovariances(X, Y, ggplot2=TRUE)
#' 
#' ## Plot options
#' plot_autocovariances(X, Y, max_lag=7, ggplot2=TRUE,
#'                      plot_options=list(series_names = c("My AR1", "Your AR1"),
#'                                        bar_colors = c("forestgreen", "gray20")))
#'                                      
#' #########
#' ## A mulitviariate example
#' set.seed(12345)
#' Sigma11 <- matrix(c(1, 0.6, 0.6, 1), 2, 2) # covariance of errors
#' Sigma12 <- matrix(c(0.1, 0.2, 0.2, 0.1), 2, 2) # Cross-covariance
#' Sigma <- rbind(cbind(Sigma11,Sigma12), cbind(Sigma12,Sigma11))
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0, 0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X1 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1]) 
#' X2 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 2]) 
#' Y1 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 3])
#' Y2 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 4])
#' X <- cbind(X1, X2)
#' Y <- cbind(Y1, Y2)
#' 
#' plot_autocovariances(X, Y, max_lag=7,
#'                      plot_options=list(series_names = c("My VAR1", "Your VAR1"),
#'                                        bar_colors = c("forestgreen", "gray20"),
#'                                        var_names = c("Dim One", "Dim 2")))
#' 
#' plot_autocovariances(X, Y, max_lag=7, ggplot2=TRUE,
#'                      plot_options=list(series_names = c("My VAR1", "Your VAR1"),
#'                                        bar_colors = c("forestgreen", "gray20"),
#'                                        var_names = c("Dim One", "Dim 2")))

plot_autocovariances <- function(X, Y, max_lag=NULL, 
                                 ggplot2=FALSE,
                                 plot_options=list(series_names=c("Series 1", "Series 2"),
                                                   bar_colors = c("gray10", "gray50"),
                                                   var_names = NULL)) {
  #####################################################
  ## Start with some error checking
  #####################################################
  ## Dimensionality of supplied data
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if(ncol(X) != ncol(Y) ) 
    stop("Series X and Y must be of the same dimension")
  k <- ncol(X)
  if(nrow(X) != nrow(Y) ) 
    stop("Series X and Y must be of the same length")
  N <- nrow(X)
  
  
  #####################
  ## Check lag
  if(is.null(max_lag) ) 
    max_lag <- ceiling((log2(N))^0.9999)
  if(max_lag<0) 
    stop("A non-negative maximum lag must be supplied")

  #################
  ## For ggplot2 into a logical (if the user does something silly)
  ggplot2 = as.logical(ggplot2)
  
  ######################################
  ## Check the supplied series names
  ##  Must be of dimension 2
  if(!is.null(plot_options$series_names) && length(plot_options$series_names)!=2)
    warning("Length of series names must match dimension of two, using defaults")

  ## If no series_names is supplied or not length 2,
  ## user the defaults: Series 1, Series 2
  if(is.null(plot_options$series_names) || length(plot_options$series_names)!=2)  
    plot_options$series_names = c("Series 1", "Series 2")
  
  
  ##########################
  ## Check the supplied variable names
  ## Need to the be correct dimensions
  if(!is.null(plot_options$var_names) && length(plot_options$var_names)!=k)
    warning("Length of variable names must match dimension of series, using defaults")

  ## If no variable names are supplied, or not correct dimension
  ## first try and get the names from the X or Y
  ## If that fails, then do a default naming structure
  ##      Var1, Var2, etc...
  if(is.null(plot_options$var_names) || length(plot_options$var_names)!=k) {
    if (is.null(colnames(X))){
      data_var_names <- colnames(Y)
    } else {
      data_var_names <- colnames(X)
    }
    
    if(is.null(data_var_names) ) {
      plot_options$var_names <- paste0("Var", 1:k) 
    } else {
      plot_options$var_names <- data_var_names
    }
  }


  #############
  ## Error check supplied colors
  ##
  ## If no colors are supplied or not enough, reset the default
  if(is.null(plot_options$bar_colors) || length(plot_options$bar_colors)<2 ) {
    plot_options$bar_colors = c("gray10", "gray50")
    warning("No colors, or incorrect number of colors supplied, using defaults")
  }
  ## Check for valid colors
  if(!tryCatch(is.matrix(grDevices::col2rgb(plot_options$bar_colors[1])), 
           error = function(e) FALSE) || !tryCatch(is.matrix(grDevices::col2rgb(plot_options$bar_colors[2])), 
                                                error = function(e) FALSE) ) {
    warning("Not valid color choices, using defaults")
    plot_options$bar_colors = c("gray10", "gray50")
  }
  
  ########################
  ## compute the ACVFs and call
  ## the internal functions that 
  ## actually does the plotting
  acvf <- stats::acf(cbind(X, Y), lag.max = max_lag, type = "covariance", plot = FALSE)$acf
  
  if(ggplot2) {
    pack_attached <- search()
    ggp_attached <- sum(ifelse(pack_attached == "package:ggplot2", 1, 0) )
    if(!ggp_attached)
      stop("The package ggplot2 must be installed and attached to the current R session")
    return(acvfPlot_ggplot(acvf=acvf, N=N, k=k, max_lag=max_lag, plot_options=plot_options))
  } else {
    acvfPlot_base(acvf=acvf, N=N, k=k, max_lag=max_lag, plot_options=plot_options)
  }
}




#' Plot Difference of Autocovariance Functions
#' 
#' @description Will built a grid of \eqn{m x m} plotting the difference in the autocovariance function of two \eqn{m}-dimensional time series based on the output of autocovariance_test()
#' 
#' @param x acvfTest object, result of autocovariance_test
#' @param bar_color The color for the bars in the plot of differences, defaults at \code{"gray10"}.
#' @param var_names An optional vector of length matching the dimension of the time series. Only used when greater than one to label the grid plot.
#' @param include_ci A flag to indicate whether to include 95\% Simultaneous Confidence Intervals on the differences of autocovariances. See details.
#' 
#' @return An unnamed list containing base R graphic plot elements is supplied.
#' 
#' @details Confidence bands are plotted when \code{include_ci=TRUE} based on 95\% simultaneous confidence intervals. 
#' 
#' Specifically, a vector \eqn{\boldsymbol{\Delta}} is calculated as the difference of the autocovariance function between two series. Based on Bartlett's formula we calculate a covariance matrix for \eqn{\boldsymbol{\Delta}}, called \eqn{\mathbf{W}}.
#' The quadratic form \eqn{n\boldsymbol{\Delta}'\mathbf{W}^{-1}\boldsymbol{\Delta}} is asymptotically a Chi-square random variable. 95\% critical values are calculated from the corresponding chi square distribution along with marginal variance terms from \code{diag}\eqn{(\mathbf{W})} to construct a 95\% simultanous confidence interval for each autocovarince difference in \eqn{\boldsymbol{\Delta}}. 
#' 
#' It is worth noting the multivariate nature of the proposed method. It is possible to have a significant result in \code{autocovariance_test} but for all differences to be within their respective simultaneous intervals.
#' 
#' @export 
#' 
#' @example 
#' 
#' set.seed(12345)
#' Sigma11 <- matrix(c(1, 0.6, 0.6, 1), 2, 2) # covariance of errors
#' Sigma12 <- matrix(c(0.1, 0.2, 0.2, 0.1), 2, 2) # Cross-covariance
#' Sigma <- rbind(cbind(Sigma11,Sigma12), cbind(Sigma12,Sigma11))
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0, 0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X1 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1]) 
#' X2 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 2]) 
#' Y1 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 3])
#' Y2 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 4])
#' X <- cbind(X1, X2)
#' Y <- cbind(Y1, Y2)
#' 
#' out <- autocovariance_test(X, Y, max_lag=5, test="Fixed_Lag")
#' plot(out)
#' plot(out, bar_color="blue")
#' plot(out, bar_color="blue", var_names=c("Dim 1", "Dim two"))
#' plot(out, bar_color="blue", var_names=c("Dim 1", "Dim two"), include_ci=TRUE)

plot.acvfTest <- function(x, 
                          bar_color="gray20", 
                          var_names=NULL,
                          include_ci=FALSE) {
  
  # Shape data for plotting
  
  k_sum <- x$k*(x$k+1)/2
  
  lag0_mat <- matrix(nrow=x$k, ncol=x$k, 0)
  cov0_mat <- matrix(nrow=x$k, ncol=x$k, 0)
  cnt <- 0
  for(i in 1:x$k) {
    for(j in i:x$k) {
      cnt <- cnt+1
      lag0_mat[i,j] <- x$delta[ cnt,1 ]
      cov0_mat[i,j] <- x$dep_cov[cnt,cnt]
    }
  }
  
  if(x$k>1) {
    lag0_mat <- lag0_mat + t(lag0_mat) - diag(diag(lag0_mat))
    cov0_mat <- cov0_mat + t(cov0_mat) - diag(diag(cov0_mat))
  }
  
  acvf_diff <- c(lag0_mat, x$delta[(k_sum+1):length(x$delta)] )
  
  ### Quick error check on var_names
  if(is.null(var_names) ) {
    if(is.null(x$var_names)) {
      var_names <- paste0("Var", 1:x$k)
    } else {
      var_names <- x$var_names
    }
  }
  ## Error check on supplied color
  if(!tryCatch(is.matrix(grDevices::col2rgb(bar_color)), 
               error = function(e) FALSE) ) {
    warning("Not a valid color choice, using default")
    bar_color = c("gray10")
  }
  
  data_acvf <- data.frame(
    lags= factor(rep(0:x$max_lag, each=x$k*x$k) ),
    parameter1 = rep(var_names, x$k*(x$max_lag+1)),
    parameter2 = rep(rep(var_names, each= x$k),(x$max_lag+1)),
    acvf = acvf_diff,
    acvf_cov <- c(cov0_mat, diag(x$dep_cov)[(k_sum+1):length(diag(x$dep_cov))] ),
    timeseries="One")
  
  # Y axis limits
  ymax <- max(data_acvf$acvf)
  ymin <- min(data_acvf$acvf)
  
  # Split data by dimensions
  data_list <- split(data_acvf, paste(data_acvf$parameter2, data_acvf$parameter1))
  
  # Adjust margins and such
  graphics::par(mfcol = c(x$k, x$k), oma = c(3, 3, 2, 2), mar = c(1, 1, 0, 0), mgp = c(1, 1, 0), xpd = NA)
  # Plot 
  if(include_ci) {
    for (i in 1:(x$k^2)){
      barplotACVF_ci(data_list[[i]], ymin, ymax, x$k, i, x$max_lag, x$n, plot_options=list(bar_colors=bar_color))
    }
  } else {
    for (i in 1:(x$k^2)){
      barplotACVF(data_list[[i]], ymin, ymax, x$k, i, x$max_lag, plot_options=list(bar_colors=bar_color))
    }
  }
  graphics::par(mfrow = c(1, 1))
  
  # print the overall labels
  graphics::mtext('Lag', side = 1, outer = TRUE, line = 2)
  graphics::mtext("Difference in ACVFs", side = 2, outer = TRUE, line = 2)
  at <- seq(1 / (2 * x$k), 1 - 1 / (2 * x$k), by = 1 / x$k)
  if(x$k>1) {
    graphics::mtext(var_names, side = 3, outer = TRUE, line = 0, at = at)
    graphics::mtext(rev(var_names), side = 4, outer = TRUE, line = 0, at = at)
  }
  
  # Reset graphics
  graphics::par(oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), xpd = FALSE)
  
}



#' Plot Difference of Autocovariance Functions in ggplot2
#' 
#' @description Will built a grid of \eqn{m x m} plotting the difference in the autocovariance function of two \eqn{m}-dimensional time series based on the output of autocovariance_test()
#' 
#' @param x acvfTest object, result of autocovariance_test
#' @param bar_color The color for the bars in the plot of differences, defaults at \code{"gray10"}.
#' @param var_names An optional vector of length matching the dimension of the time series. Only used when greater than one to label the grid plot.
#' @param include_ci A flag to indicate whether to include 95\% Simultaneous Confidence Intervals on the differences of autocovariances. See details.
#' 
#' @return A \code{gg} or \code{ggplot} object is supplied.
#' 
#' @details Confidence bands are plotted when \code{include_ci=TRUE} based on 95\% simultaneous confidence intervals. 
#' 
#' Specifically, a vector \eqn{\boldsymbol{\Delta}} is calculated as the difference of the autocovariance function between two series. Based on Bartlett's formula we calculate a covariance matrix for \eqn{\boldsymbol{\Delta}}, called \eqn{\mathbf{W}}.
#' The quadratic form \eqn{n\boldsymbol{\Delta}'\mathbf{W}^{-1}\boldsymbol{\Delta}} is asymptotically a Chi-square random variable. 95\% critical values are calculated from the corresponding chi square distribution along with marginal variance terms from \code{diag}\eqn{(\mathbf{W})} to construct a 95\% simultaneous confidence interval for each autocovarince difference in \eqn{\boldsymbol{\Delta}}. 
#' 
#' It is worth noting the multivariate nature of the proposed method. It is possible to have a significant result in \code{autocovariance_test} but for all differences to be within their respective simultaneous intervals.
#' 
#' @export 
#' 
#' @example 
#' 
#' set.seed(12345)
#' Sigma11 <- matrix(c(1, 0.6, 0.6, 1), 2, 2) # covariance of errors
#' Sigma12 <- matrix(c(0.1, 0.2, 0.2, 0.1), 2, 2) # Cross-covariance
#' Sigma <- rbind(cbind(Sigma11,Sigma12), cbind(Sigma12,Sigma11))
#' n <- 1000 # sample size
#' Errors <- mvtnorm::rmvnorm(n, mean = c(0, 0, 0, 0), sigma = Sigma)
#' 
#' # Simultate correlated AR(1) processes under the null hypothesis
#' X1 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 1]) 
#' X2 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 2]) 
#' Y1 <- arima.sim(model = list(ar  = 0.5), n, innov = Errors[ , 3])
#' Y2 <- arima.sim(model = list(ar  = 0.8), n, innov = Errors[ , 4])
#' X <- cbind(X1, X2)
#' Y <- cbind(Y1, Y2)
#' 
#' out <- autocovariance_test(X, Y, max_lag=5, test="Fixed_Lag")
#' autoplot(out)
#' autoplot(out, bar_color="blue")
#' autoplot(out, bar_color="blue", var_names=c("Dim One", "Dim 2"))
#' autoplot(out, bar_color="blue", var_names=c("Dim 1", "Dim two"), include_ci=TRUE)

autoplot.acvfTest <- function(x, 
                              bar_color="gray20", 
                              var_names=NULL,
                              include_ci=FALSE) {
  ## ggplot2 check
  pack_attached <- search()
  ggp_attached <- sum(ifelse(pack_attached == "package:ggplot2", 1, 0) )
  if(!ggp_attached)
    stop("The package ggplot2 must be installed and attached to the current R session")
  
  ## Error check on supplied color
  if(!tryCatch(is.matrix(grDevices::col2rgb(bar_color)), 
               error = function(e) FALSE) ) {
    warning("Not a valid color choice, using default")
    bar_color = c("gray20")
  }
  
  if(x$k >= 2) {
    k_sum <- x$k*(x$k+1)/2
    
    lag0_mat <- matrix(nrow=x$k, ncol=x$k, 0)
    cov0_mat <- matrix(nrow=x$k, ncol=x$k, 0)
    cnt <- 0
    for(i in 1:x$k) {
      for(j in i:x$k) {
        cnt <- cnt+1
        lag0_mat[i,j] <- x$delta[ cnt,1 ]
        cov0_mat[i,j] <- x$dep_cov[cnt,cnt]
      }
    }
    lag0_mat <- lag0_mat + t(lag0_mat) - diag(diag(lag0_mat))
    cov0_mat <- cov0_mat + t(cov0_mat) - diag(diag(cov0_mat))
    
    acvf_diff <- c(lag0_mat, x$delta[(k_sum+1):length(x$delta)] )
    acvf_diff_cov <- c(cov0_mat, diag(x$dep_cov)[(k_sum+1):length(diag(x$dep_cov))] )
    
    ### Quick error check on var_names
    if(is.null(var_names) ) {
      if(is.null(x$var_names)) {
        var_names <- paste0("Var", 1:x$k)
      } else {
        var_names <- x$var_names
      }
    }

    data_acvf <- data.frame(
      lag= rep(0:x$max_lag, each=x$k*x$k),
      parameter1 = rep(var_names, x$k*(x$max_lag+1)),
      parameter2 = rep(rep(var_names, each= x$k),(x$max_lag+1)),
      ACVF = acvf_diff,
      ACVF_cov = acvf_diff_cov)
    data_acvf$lag <- data_acvf$lag - 0.1
    data_acvf$lag_end <- data_acvf$lag + 0.2
    data_acvf$parameter1 <- factor(data_acvf$parameter1, levels=var_names)
    data_acvf$parameter2 <- factor(data_acvf$parameter2, levels=var_names)
    
    acvf_plot <- ggplot2::ggplot(data = data_acvf) +
      ggplot2::geom_rect(ggplot2::aes_string(xmin="lag", xmax="lag_end", ymin="0", ymax="ACVF"),
                color=bar_color, fill=bar_color) +
      ggplot2::facet_grid(data_acvf$parameter1 ~ data_acvf$parameter2) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = "0"), colour="gray40") + 
      ggplot2::labs(y="Difference in ACVFs", x="Lag")
  } else {
    data_acvf <- data.frame(
      lag= 0:x$max_lag,
      ACVF = x$delta,
      ACVF_cov <- diag(x$dep_cov) )
    data_acvf$lag <- data_acvf$lag - 0.05
    data_acvf$lag_end <- data_acvf$lag + 0.1
    
    acvf_plot <- ggplot2::ggplot(data = data_acvf) +
      ggplot2::geom_rect(ggplot2::aes_string(xmin="lag", xmax="lag_end", ymin="0", ymax="ACVF"),
                color=bar_color, fill=bar_color) +
      ggplot2::geom_hline(yintercept = 0, colour="gray40") + 
      ggplot2::labs(y="Difference in ACVFs", x="Lag")
  }
  if(include_ci){
    DF <- x$k*x$k*x$max_lag + x$k*(x$k+1)/2
    data_acvf$bound1 <- -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(data_acvf$ACVF_cov/(x$n) )
    data_acvf$bound2 <- sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(data_acvf$ACVF_cov/(x$n) )
    acvf_plot <- acvf_plot + 
      ggplot2::geom_line(ggplot2::aes_string(x="lag", y="bound1"),
              linetype=2, color="gray50", data=data_acvf) +
      ggplot2::geom_line(ggplot2::aes_string(x="lag", y="bound2"),
              linetype=2, color="gray50", data=data_acvf)
  }
  acvf_plot
}
