## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
#' @useDynLib autocovarianceTesting, .registration=TRUE
#' 

# Function to run compatibility checks
compatibilityChecks <- function(X, Y, L, test, B, prewhiten){
  
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
#' @param B 
#' @param prewhiten 
#'
#' @return
#' @export
#'
#' @examples
autocovarianceTest <- function(X, Y, L, test = "Dependent", B = 500, prewhiten = TRUE){
  
  # Compatibility Tests
  compatibilityChecks(X, Y, L, test, B, prewhiten)
  
  # Get some info
  n <- nrow(X)
  k <- ncol(Y)
  
  # Initialize output list
  out <- list()
  
  # Compute independent and dependent tests if need be
  if ("Independent" %in% test | "Dependent" %in% test){
    # Compute dependent covariance
    deltaCovar <- calculateCovariance(X, Y, L)
    out <- c(out, deltaCovar)
    # Compute test statistics
    if ("Independent" %in% test){
      indTest <- calculateTestStat(deltaCovar$delta, deltaCovar$ind_cov, n, L, k)
      # Compute p-values
      indTest$pval <- pchisq(indTest$stat, df = indTest$df, lower.tail = FALSE)
      indTest$weight_pval <- pgamma(indTest$weight_stat, shape = indTest$alpha, scale = indTest$beta, lower.tail = FALSE)
      names(indTest) <- paste0("ind_",names(indTest))
      out <- c(out, indTest)
    }
    if ("Dependent" %in% test){
      depTest <- calculateTestStat(deltaCovar$delta, deltaCovar$dep_cov, n, L, k)
      # Compute p-values
      depTest$pval <- pchisq(depTest$stat, df = depTest$df, lower.tail = FALSE)
      depTest$weight_pval <- pgamma(depTest$weight_stat, shape = depTest$alpha, scale = depTest$beta, lower.tail = FALSE)
      names(depTest) <- paste0("dep_",names(depTest))
      out <- c(out, depTest)
    }
  }
  
  # Compute bootstrapped tests if need be
  if ("bootDependent" %in% test | "bootBartlett" %in% test){
    bootTest <- calculateBootTestStat(X, Y, L, B, prewhiten)
    out <- c(out, bootTest)
  }
  
  # Create a nice table to output
  
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
    out_table <- rbind(out_table, c("Independent", round(out$ind_stat, 3), out$ind_df, round(out$ind_pval, 3)))
    colnames(out_table) <- c("Test", "Chi-Sq", "df", "p-value")
    out_table_weight <- rbind(out_table_weight, c("Weighted Independent", round(out$ind_weight_stat, 3), round(out$ind_alpha, 3), round(out$ind_beta, 3), round(out$ind_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Gamma", "alpha", "beta", "p-value")
  }
  
  if ("Dependent" %in% test){
    # Add dependent tests
    out_table <- rbind(out_table, c("Dependent", round(out$dep_stat, 3), out$dep_df, round(out$dep_pval, 3)))
    colnames(out_table) <- c("Test", "Chi-Sq", "df", "p-value")
    out_table_weight <- rbind(out_table_weight, c("Weighted Dependent", round(out$dep_weight_stat, 3), round(out$dep_alpha, 3), round(out$dep_beta, 3), round(out$dep_weight_pval, 3)))
    colnames(out_table_weight) <- c("Test", "Gamma", "alpha", "beta", "p-value")
  }
  
  # Initialize bootstrapped tables
  if ("bootDependent" %in% test | "bootBartlett" %in% test){
    out_table_boot <- data.frame(0, 1, 2, 3)
    out_table_boot <- out_table_boot[-1, ]
  }
  
  if ("bootDependent" %in% test){
    # Add jin tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Jin", round(out$jin_stat, 3), round(out$jin_L, 3), round(out$jin_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statitic", "L hat", "p-value")
  }
  
  if ("bootBartlett" %in% test){
    # Add boot Bartlett tests
    out_table_boot <- rbind(out_table_boot, c("Bootstrap-Bartlett", round(out$bart_stat, 3), round(out$bart_L, 3), round(out$bart_pval, 3)))
    colnames(out_table_boot) <- c("Test", "Statitic", "L hat", "p-value")
  }
  
  # Set up knitr output
  table1 <- if ("Independent" %in% test | "Dependent" %in% test){ knitr::kable(out_table, row.names = FALSE, caption = "Fixed Lag Tests")} else {NULL}
  table2 <- if ("Independent" %in% test | "Dependent" %in% test){ knitr::kable(out_table_weight, row.names = FALSE, caption = "Weighted Fixed Lag Tests")} else {NULL}
  table3 <- if ("bootDependent" %in% test | "bootBartlett" %in% test){ knitr::kable(out_table_boot, row.names = FALSE, caption = "Automatic Lag Selection Tests")} else {NULL}
  
  print(knitr::kables(list(table1, table2, table3)))
  
  return(invisible(out))
  
}