# Plot helper functions
#
# Side-by-side barplots (or solo) without confidence bands
# Used by the plot_autocovariances() or plot.acvfTest() functions

barplotACVF <- function(dataset, ymin, ymax, k, time, max_lag, plot_options=list(series_names=c("Series 1", "Series 2"),
                                                                           bar_colors = c("gray10", "gray40"))){
  M <- length(unique(dataset$timeseries))
  if(M==1) {
    line_len <- 2.4*max_lag   ## works well with only 1 set of bars
  }
  else
    line_len <- 3.5*max_lag   ## Works well with side-by-side bars
  
  if (k == 1){
    graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), data = dataset)
    graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
  } else {
    if (!(time %in% 1:k) & !(time %% k == 0)){
      graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", yaxt = "n", data = dataset)
      graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
    } else {
      if (time %in% 1:k & !(time %% k == 0)){
        graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", data = dataset)
        graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
      } else {
        if (!(time %in% 1:k) & time %% k == 0){
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), yaxt = "n", data = dataset)
          graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
        } else {
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), data = dataset)
          graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
        }
      }
    }
  }
}

# Side-by-side barplots (or solo) WITH confidence bands
# Used by the plot.acvfTest() functions if the user
# chooses to include 95% simultaneous chi-square based intervals
barplotACVF_ci <- function(dataset, ymin, ymax, k, time, max_lag, N, plot_options=list(series_names=c("Series 1", "Series 2"),
                                                                                 bar_colors = c("gray10", "gray40"))){
  M <- length(unique(dataset$timeseries))
  DF <- max_lag*k*k + k*(k+1)/2
  if(M==1) {
    line_len <- 2.4*max_lag   ## works well with only 1 set of bars
  }
  else 
    line_len <- 3.5*max_lag   ## Works well with side-by-side bars
  
  x_ci <- seq(1.5, line_len, line_len/length(unique(dataset$lags) ) )[1:length(unique(dataset$lags))]  
  if (k == 1){
    graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), data = dataset)
    graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
    graphics::lines(x = x_ci, y = sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
    graphics::lines(x = x_ci, y = -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
  } else {
    if (!(time %in% 1:k) & !(time %% k == 0)){
      graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", yaxt = "n", data = dataset)
      graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
      graphics::lines(x = x_ci, y = sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
      graphics::lines(x = x_ci, y = -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
    } else {
      if (time %in% 1:k & !(time %% k == 0)){
        graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), xaxt = "n", data = dataset)
        graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
        graphics::lines(x = x_ci, y = sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
        graphics::lines(x = x_ci, y = -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
      } else {
        if (!(time %in% 1:k) & time %% k == 0){
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), yaxt = "n", data = dataset)
          graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
          graphics::lines(x = x_ci, y = sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
          graphics::lines(x = x_ci, y = -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
        } else {
          graphics::barplot(acvf ~ timeseries + lags, xlab= "", ylab = "", col = plot_options$bar_colors, beside = TRUE, ylim = c(ymin, ymax), data = dataset)
          graphics::lines(x = c(0, line_len ), y = c(0, 0), col="gray40")
          graphics::lines(x = x_ci, y = sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
          graphics::lines(x = x_ci, y = -sqrt(stats::qchisq(0.95, df=(DF)) )*sqrt(dataset$acvf_cov/(N) ), lty=2, col="gray50")
        }
      }
    }
  }
}


# Function to create a autocovariance difference plot
# using Base R graphics
acvfPlot_base <- function(acvf, N, k, max_lag, plot_options){
  
  # Shape data for plotting
  acvfX <- data.frame(acvf = c(acvf[ , 1:k, 1:k]), lags = 0:max_lag, dim1 = rep(paste0("Dim", 1:k), each = (max_lag + 1) * k ), dim2 = rep(rep(paste0("Dim", 1:k), each = max_lag + 1), k), timeseries = "X")
  acvfY <- data.frame(acvf = c(acvf[ , (k + 1):(2 * k), (k + 1):(2 * k)]), lags = 0:max_lag, dim1 = rep(paste0("Dim", 1:k), each = (max_lag + 1) * k ), dim2 = rep(rep(paste0("Dim", 1:k), each = max_lag + 1), k), timeseries = "Y")
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
    barplotACVF(data_list[[i]], ymin, ymax, k, i, max_lag, plot_options=plot_options)
  }
  graphics::par(mfrow = c(1, 1))
  
  # print the overall labels
  graphics::mtext('Lag', side = 1, outer = TRUE, line = 2)
  graphics::mtext('Autocovariance Function', side = 2, outer = TRUE, line = 2)
  at <- seq(1 / (2 * k), 1 - 1 / (2 * k), by = 1 / k)
  if(k>1) {
    graphics::mtext(plot_options$var_names, side = 3, outer = TRUE, line = 0, at = at)
    graphics::mtext(rev(plot_options$var_names), side = 4, outer = TRUE, line = 0, at = at)
  }
  
  # add legend
  graphics::legend("topright", legend = plot_options$series_names, fill = plot_options$bar_colors, cex = 1/k)
  # Reset graphics
  graphics::par(oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0), xpd = FALSE)
  
}

# Function to create a autocovariance difference plot
# or the side-by-side autocovariances plot
# using ggplot2 based graphics
acvfPlot_ggplot <- function(acvf, N, k, max_lag, plot_options) {
  if(k==1) {
    data_acvf <- data.frame(
      lag=rep(0:max_lag, 2),
      series=c(rep(plot_options$series_names[1], max_lag+1), 
               rep(plot_options$series_names[2], max_lag+1)), 
      ACF=c(acvf[,1,1], acvf[,2,2]))
    data_acvf$lag <- data_acvf$lag + 0.005*c(rep(-1,(max_lag+1)*k*k), rep(1,(max_lag+1)*k*k))
    data_acvf$lag_end <- data_acvf$lag + 0.16*c(rep(-1,(max_lag+1)*k*k), rep(1,(max_lag+1)*k*k))
    data_acvf$series <- factor(data_acvf$series, levels=plot_options$series_names)
    
    acvf_plot <- ggplot2::ggplot() +
      ggplot2::geom_rect(ggplot2::aes(xmin=data_acvf$lag, xmax=data_acvf$lag_end, ymin=0, ymax=data_acvf$ACF,
                    color=data_acvf$series, fill=data_acvf$series)) +
      ggplot2::scale_color_manual(values = plot_options$bar_colors) +
      ggplot2::scale_fill_manual(values = plot_options$bar_colors) +
      ggplot2::geom_hline(yintercept = 0, colour="gray40") + 
      ggplot2::labs(y="Autocovariance Function", x="Lag", color = "Series", fill = "Series") + 
      ggplot2::theme(legend.position="bottom")
  } else {
    data_acvf <- data.frame(
      lag=rep(0:max_lag, 2*k*k),
      series=rep(plot_options$series_names, each=(max_lag+1)*k*k),
      parameter2 = rep(rep(plot_options$var_names, each= k*(max_lag+1)),2),
      parameter1 = rep(rep(plot_options$var_names, each= (max_lag+1)),2*k),
      ACF = c(as.vector(acvf[,1:k,1:k]), as.vector(acvf[,(k+1):(2*k),(k+1):(2*k)]) ) )
    data_acvf$lag <- data_acvf$lag + 0.015*c(rep(-1,(max_lag+1)*k*k), rep(1,(max_lag+1)*k*k))
    data_acvf$lag_end <- data_acvf$lag + 0.21*c(rep(-1,(max_lag+1)*k*k), rep(1,(max_lag+1)*k*k))
    data_acvf$series <- factor(data_acvf$series, levels=plot_options$series_names)
    data_acvf$parameter1 <- factor(data_acvf$parameter1, levels=plot_options$var_names)
    data_acvf$parameter2 <- factor(data_acvf$parameter2, levels=plot_options$var_names)
    
    acvf_plot <- ggplot2::ggplot(data = data_acvf) +
      ggplot2::geom_rect(ggplot2::aes_string(xmin="lag", xmax="lag_end", ymin="0", ymax="ACF",
                                      color="series", fill="series")) +
      ggplot2::geom_hline(ggplot2::aes_string(yintercept = "0"), colour="gray40") + 
      ggplot2::facet_grid(data_acvf$parameter1 ~ data_acvf$parameter2) +
      ggplot2::scale_color_manual(values = plot_options$bar_colors) +
      ggplot2::scale_fill_manual(values = plot_options$bar_colors) +
      ggplot2::labs(y="Autocovariance Function", x="Lag", color = "Series", fill = "Series") + 
      ggplot2::theme(legend.position="bottom")
    
  }
  acvf_plot
}

