test_that("Functions work", {
  
  ## Functions that were previously made in R. Now used for testing C++ versions.
  
  library(tidyverse)
  library(Matrix)
  library(mvtnorm)
  library(matrixcalc)
  library(parallel)
  library(testthat)
  library(forecast)
  library(Rcpp)
  library(autocovarianceTesting)
  library(microbenchmark)
  
  ################################
  ## get.covar.matrix()
  ##
  ## Function to compute autocovariance covariance matricies
  ## This function is called via the fixed lag code
  ## and the automated lag selection code.
  ##
  
  
  get.covar.matrix <- function(data, lag_L=c(5, 10), corr_only=FALSE){
    
    # Let's do the correlated version first
    L <- max(lag_L)
    n <- nrow(data)
    k <- ncol(data)/2
    
    # Naming data columns for convenience 
    colnames(data) <- sort(levels(interaction(c("X", "Y"), 1:k, sep="")))
    
    # Get all covariances (seems silly but I need the negative cross-covariance)
    ccf_index <- expand.grid(1:(2*k), 1:(2*k))
    ccfs <- mapply(function(one, two){(ccf(data[,one], data[,two], 
                                           lag.max = n-1, type = "covariance", 
                                           plot=FALSE, demean=TRUE)$acf)},
                   ccf_index[,1], ccf_index[,2]
    )
    
    # Give ccfs appropriate names
    colnames(ccfs) <- levels(interaction(levels(interaction(c("X", "Y"), 1:k, sep="")), 
                                         levels(interaction(c("X", "Y"), 1:k, sep="")), sep=""))
    # ccfs contains acf/ccf from -(n-1) to n-1
    
    # get vector of differences
    diff_vec <- c(t(as.matrix(ccfs[n + 0:L,grepl("X.X", colnames(ccfs))])[,order(colnames(ccfs)[grepl("X.X", colnames(ccfs))])] -
                      as.matrix(ccfs[n + 0:L,grepl("Y.Y", colnames(ccfs))])[,order(colnames(ccfs)[grepl("Y.Y", colnames(ccfs))])]))
    
    # replace by H0
    ccfs_H0 <- ccfs
    ccfs_H0[,sort(colnames(ccfs)[grepl("Y.Y", colnames(ccfs))])] <- ccfs_H0[,sort(colnames(ccfs)[grepl("X.X", colnames(ccfs))])] <- (ccfs_H0[,sort(colnames(ccfs)[grepl("X.X", colnames(ccfs))])] + ccfs_H0[,sort(colnames(ccfs)[grepl("Y.Y", colnames(ccfs))])])/2
    
    # Cute function to look up ccf
    # Ex: ts1 = "X1", ts2="Y2", lags=c(0,1,2) 
    ccf_lookup <- function(lags, ts1=NULL, ts2=NULL){
      if (is.null(ts1)){
        ccfs_H0[n + lags,]
      } else {
        ccfs_H0[n + lags, colnames(ccfs) %in% paste0(ts1, ts2)] 
      }
    }
    
    # truncation
    r <- seq(-floor(n^(1/3) + 0.00000001), floor(n^(1/3) + 0.00000001), 1)
    #r <- seq(-floor(sqrt(n)), floor(sqrt(n)), 1)
    # parameters
    # Find lower triangle and don't compute those terms
    loop_over <- expand.grid(colnames(ccfs), 0:L,  colnames(ccfs), 0:L)[,c(1,3,2,4)][c(t(upper.tri(matrix(NA, (L+1)*4*(k^2), (L+1)*4*(k^2)), diag=TRUE))),]
    
    # covariance matrix terms
    formula <- function(cov1, cov2, p, q){
      i <- substr(cov1, 1, 2)
      j <- substr(cov1, 3, 4)
      k <- substr(cov2, 1, 2)
      l <- substr(cov2, 3, 4)
      sum(ccf_lookup(r, i, k) * ccf_lookup(r - p + q, j, l) + ccf_lookup(r + q, i, l) * ccf_lookup(r - p, j, k))
    }
    
    terms <- mcmapply(formula, loop_over[,1], loop_over[,2], loop_over[,3], loop_over[,4])
    
    # Create dependent covariance matrix
    cov_mat <- matrix(0, nrow=ncol(ccfs)*(L+1) , ncol= ncol(ccfs)*(L+1) , byrow=T)
    cov_mat[lower.tri(cov_mat, diag=TRUE)] <- terms
    cov_mat <- cov_mat + t(cov_mat) - diag(diag(cov_mat))
    
    # create A matrix
    A <- bdiag(rep(list(cbind(diag(2*k^2 - k), matrix(0,ncol=2*k, nrow=(2*k^2 - k)), -diag(2*k^2 - k))[rep(rep(c(T,F), k), each=k)[1:(2*k^2 - k)],]), L+1))
    if(k==1) A <- t(A)
    eta <- as.matrix(c(t(ccfs[n + 0:L,])))
    
    if (corr_only==TRUE){
      list(tcrossprod(A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))], eta[1:(ncol(ccfs)*(L+1))]), 
           tcrossprod(tcrossprod(A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))], cov_mat[1:(ncol(ccfs)*(L+1)),1:(ncol(ccfs)*(L+1))]),A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))])) #independent
    } else {
      # Create independent covariance matrix
      cov_mat2 <- matrix(0, nrow=ncol(ccfs)*(L+1) , ncol= ncol(ccfs)*(L+1) , byrow=T)
      cov_mat2[lower.tri(cov_mat2, diag=TRUE)] <- (cbind(loop_over, terms) %>% mutate(terms = ifelse(grepl("X.X", Var1) & grepl("X.X", Var3), terms,
                                                                                                     ifelse(grepl("Y.Y", Var1) & grepl("Y.Y", Var3), terms, 0))))$terms
      cov_mat2 <- cov_mat2 + t(cov_mat2) - diag(diag(cov_mat2))
      
      list(tcrossprod(A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))], eta[1:(ncol(ccfs)*(L+1))]), 
           tcrossprod(tcrossprod(A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))], cov_mat[1:(ncol(ccfs)*(L+1)),1:(ncol(ccfs)*(L+1))]),A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))]), #dependent
           tcrossprod(tcrossprod(A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))], cov_mat2[1:(ncol(ccfs)*(L+1)),1:(ncol(ccfs)*(L+1))]),A[1:(k^2*(L+1)), 1:(ncol(ccfs)*(L+1))])) #independent
    }
    
  }
  
  #########################################
  ##
  ## Code for the fixed lag testing
  ##
  ## Essentially comparing Lund (2009) to 
  ##  the linearly dependent version
  ## 
  ## This code will handle both univariate
  ##  and multivariate comparison
  ##
  ## Written by Dan Cirkovic
  ##  edited by Tom Fisher (tjf)
  
  ##################################
  ## get.one.fixed.lag.test()
  ##
  ## This function will perform multiple test at specified
  ##   lags, lag_L, performing the Lund test (assuming independence)
  ##   and our new version that models linear dependence
  ## The data is passed in as one big matrix where the first
  ##   k columns are time series one, and the last k are the second
  ##
  get.one.fixed.lag.test <- function(data, lag_L=c(5, 10)){
    # Let's do the correlated version first
    L <- max(lag_L)
    n <- nrow(data)
    k <- ncol(data)/2
    
    statistics <- get.covar.matrix(data, lag_L)
    
    corr_results <- sapply(lag_L, function(x){
      # Get test stat
      if (k == 1) {
        base_stat <- (n)*tcrossprod(crossprod(statistics[[1]][1:(k^2*(x+1)),], solve(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))])), t(statistics[[1]][1:(k^2*(x+1)),]))
        weight_stat <- (n)*sum(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)*(statistics[[1]][1:(k^2*(x+1)),])^2)
        K1 <- sum(diag(statistics[[2]][1:(k^2*(x+1)), 1:(k^2*(x+1))])*rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2))
        K2 <- 2*sum(diag(tcrossprod(crossprod(tcrossprod(statistics[[2]][1:(k^2*(x+1)), 1:(k^2*(x+1))],diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2))),statistics[[2]][1:(k^2*(x+1)), 1:(k^2*(x+1))]),diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)))))
      } else {
        base_stat <- (n)*tcrossprod(crossprod(matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),], solve(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),])
        weight_stat <- (n)*sum(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]*(matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),])^2)
        K1 <- sum(diag(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])*rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])
        K2 <- 2*sum(diag(tcrossprod(crossprod(tcrossprod(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)],diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]),diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))))
      }
      
      c(as.numeric(base_stat), weight_stat, (K1^2)/K2, (K2)/K1)
    })
    
    cor_base_pvals <- pchisq(corr_results[1,], df=(k*(k+1)/2 + lag_L*k^2), lower.tail=FALSE)
    cor_weight_pvals <- pgamma(corr_results[2,], shape=corr_results[3,], scale=corr_results[4,], lower.tail=FALSE)
    
    ############### LUND #######################
    # Do the same stuff except cross cors are 0
    
    weight_results <- sapply(lag_L, function(x){
      # Get test stat
      if (k == 1) {
        base_stat <- (n)*tcrossprod(crossprod(statistics[[1]][1:(k^2*(x+1)),], solve(statistics[[3]][1:(k^2*(x+1)),1:(k^2*(x+1))])), t(statistics[[1]][1:(k^2*(x+1)),]))
        weight_stat <- (n)*sum(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)*(statistics[[1]][1:(k^2*(x+1)),])^2)
        K1 <- sum(diag(statistics[[3]][1:(k^2*(x+1)), 1:(k^2*(x+1))])*rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2))
        K2 <- 2*sum(diag(tcrossprod(crossprod(tcrossprod(statistics[[3]][1:(k^2*(x+1)), 1:(k^2*(x+1))],diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2))),statistics[[3]][1:(k^2*(x+1)), 1:(k^2*(x+1))]),diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)))))
      } else {
        base_stat <- (n)*tcrossprod(crossprod(matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),], solve(statistics[[3]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),])
        weight_stat <- (n)*sum(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]*(matrix(statistics[[1]][1:(k^2*(x+1)),])[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),])^2)
        K1 <- sum(diag(statistics[[3]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])*rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])
        K2 <- 2*sum(diag(tcrossprod(crossprod(tcrossprod(statistics[[3]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)],diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), statistics[[3]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]),diag(rep(rev(seq(0,1, by=1/(x+1))[-1]), each=k^2)[-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))))
      }
      
      c(as.numeric(base_stat), weight_stat, (K1^2)/K2, (K2)/K1)
    })
    
    lund_base_pvals <- pchisq(weight_results[1,], df=(k*(k+1)/2 + lag_L*k^2), lower.tail=FALSE)
    lund_weight_pvals <- pgamma(weight_results[2,], shape=weight_results[3,], scale=weight_results[4,], lower.tail=FALSE)
    
    
    pvals <- data.frame(lag = rep(lag_L, 4), method=c(rep("Lund Base", length(lag_L)), rep("Lund Weight", length(lag_L)), rep("Corr Base", length(lag_L)), rep("Corr Weight", length(lag_L))),
                        pval = c(lund_base_pvals, lund_weight_pvals, cor_base_pvals, cor_weight_pvals), 
                        stats = c(weight_results[1,], weight_results[2,], corr_results[1,], corr_results[2,]), stringsAsFactors = FALSE)
    
    pvals
  }
  
  # order selection test (independent and dependent)
  get.one.order.test <- function(data, L=NULL, alpha=0.05, B=500, prewhiten=TRUE){
    n <- nrow(data)
    k <- ncol(data)/2
    
    if(prewhiten == TRUE){
      # prewhiten
      p <- ceiling(min(10*log(n), 0.8*sqrt(n)))
      # Get "null" beta coefficients
      get.v.vec <- function(t){matrix(c(1,c(t(data[t:(t-p+1),]))))}
      V <- t(sapply(p:(n-1), get.v.vec))
      betas <- unname(solve(t(V)%*%V)%*%t(V)%*%data[-(1:p),])
      average.betas <- function(r){
        rbind(
          cbind((betas[r:(r+k-1),1:(1+k-1)] + betas[(r+k):(r+2*k-1),(k+1):(2*k)])/2,
                (betas[r:(r+k-1),(k+1):(2*k)] + betas[(r+k):(r+2*k-1),1:(1+k-1)])/2),
          cbind((betas[r:(r+k-1),(k+1):(2*k)] + betas[(r+k):(r+2*k-1),1:(1+k-1)])/2,
                (betas[r:(r+k-1),1:(1+k-1)] + betas[(r+k):(r+2*k-1),(k+1):(2*k)])/2)
        )
      }
      betas <- rbind(rep((betas[1,1:k] + betas[1,(k+1):(2*k)])/2,2), do.call(rbind,sapply(seq(2, nrow(betas), by=k*2), average.betas, simplify = F)))
      data <- data[-(1:p),] - V%*%betas
      n <- nrow(data)
    }
    
    # blocksize
    b <- floor(max(0.5*(n^(1/3)), 2))
    # no blocks
    K <- ceiling(n/b)
    L <- ifelse(is.null(L), ceiling(log2(n)^0.999), L)
    
    statistics <- get.covar.matrix(data, L, corr_only = TRUE)
    
    corr_results <- sapply(0:L, function(x){
      # Get test stat
      if (k == 1) {
        as.vector((n)*tcrossprod(crossprod(statistics[[1]][1:(k^2*(x+1)),], solve(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))])), t(statistics[[1]][1:(k^2*(x+1)),]))) - (2*(x+1))
      } else {
        as.vector((n)*tcrossprod(crossprod(statistics[[1]][1:(k^2*(x+1)),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)], solve(statistics[[2]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), t(statistics[[1]][1:(k^2*(x+1)),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))) - (2*(k*(k+1)/2 + x*k^2))
      }
    })
    
    ## Get bootstrap critical values
    
    # Basically sample acvfs if we can mean 0 data
    Delta <- do.call("rbind", replicate(L+1, t(matrix(data[,1:k], ncol=k)[,rep(1:k,each=k)]), simplify = FALSE)) * t(do.call("cbind",lapply(0:(L), function(x){lead(matrix(data[,1:k], ncol=k)[,rep(1:k,k)], n=x, default=0)}))) - do.call("rbind", replicate(L+1, t(matrix(data[,(k+1):(2*k)], ncol=k)[,rep(1:k,each=k)]), simplify = FALSE)) * t(do.call("cbind",lapply(0:(L), function(x){lead(matrix(data[,(k+1):(2*k)], ncol=k)[,rep(1:k,k)], n=x, default=0)})))
    # Get blocks
    getblocks = function(x){x:(x + b - 1)}
    blocks <- t(sapply(sample(1:(n - max(L, b)), size=B*K, replace=T),  getblocks))
    # Bootstrap ts indicies
    concatBlocks = function(x){c(t(blocks[((x - 1)*K + 1):(x*K),]))}
    boot_ind <- sapply(1:B, concatBlocks)[1:n,] # each column is a bs sample
    
    # Compute boot Deltas
    # each column is a bootstrap sample of the Delta vector
    Delta_l <- sapply(1:B, function(x){apply(Delta[,boot_ind[,x]], MARGIN = 1, mean)})
    Delta_bar <- matrix(apply(Delta_l, MARGIN=1, mean))
    # Compute bootstrapped covariance matricies
    boot_sigmas <- lapply(1:(L+1), function(r){lapply(1:B, function(l){(1/(B-1))*matrix(Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),])%*%t((Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),]))})})
    boot_cov <- lapply(1:B, function(x){get.covar.matrix(data[boot_ind[,x],], L)[[2]]})
    add <- function(x) Reduce("+", x)
    # Bootstrap estimated covariance matrix
    Sigma_r <- lapply(1:(L+1), function(x){n*add(boot_sigmas[[x]])})
    # Bootstrap estimated covariance matrix (our)
    Cov_r <- add(boot_cov)/(B-1)
    
    # Bootstrap distributions 
    # Us
    if (k==1){
      S_r_l_dep <- apply(sapply(1:(L+1), function(r){sapply(1:B, function(l){as.vector(n*t((Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),]))%*%solve(Cov_r[1:(k^2*(r)),1:(k^2*(r))])%*%t(t((Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),]))))})}, simplify = "array") - t(matrix(2*(0:L + 1)))[rep(1, B),], MARGIN=1, max)
      # jin 2019
      S_r_l <- apply(sapply(1:(L+1), function(r){sapply(1:B, function(l){n*t((Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),]))%*%solve(Sigma_r[[r]])%*%matrix(Delta_l[1:(k^2*r),l] - Delta_bar[1:(k^2*r),])})}, simplify = "array") - t(matrix(2*(0:L + 1)))[rep(1, B),], MARGIN=1, max)
      Sample_U_r_l <- lapply(1:(L+1), function(y){n*t(statistics[[1]][1:y,])%*%solve(Sigma_r[[y]])%*%statistics[[1]][1:y,] - 2*(y)})
    } else {
      S_r_l_dep <- apply(sapply(1:(L+1), function(r){sapply(1:B, function(l){as.vector(n*t((Delta_l[1:(k^2*r),l][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)] - Delta_bar[1:(k^2*r),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))%*%solve(Cov_r[1:(k^2*(r)),1:(k^2*(r))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE), -which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])%*%t(t((Delta_l[1:(k^2*r),l][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)] - Delta_bar[1:(k^2*r),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))))})}, simplify = "array") - t(matrix(2*(seq(3, (k*(k+1)/2 + L*k^2), by=2*k))))[rep(1, B),], MARGIN=1, max)
      # jin 2019
      S_r_l <- apply(sapply(1:(L+1), function(r){sapply(1:B, function(l){as.vector(n*t((Delta_l[1:(k^2*r),l][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)] - Delta_bar[1:(k^2*r),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))%*%solve(Sigma_r[[r]][1:(k^2*(r)),1:(k^2*(r))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE), -which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])%*%t(t((Delta_l[1:(k^2*r),l][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)] - Delta_bar[1:(k^2*r),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)]))))})}, simplify = "array") - t(matrix(2*(seq(3, (k*(k+1)/2 + L*k^2), by=2*k))))[rep(1, B),], MARGIN=1, max)
      Sample_U_r_l <- sapply(0:L, function(x){ as.vector((n)*tcrossprod(crossprod(statistics[[1]][1:(k^2*(x+1)),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)], solve(Sigma_r[[x + 1]][1:(k^2*(x+1)),1:(k^2*(x+1))][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE),-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])), t(statistics[[1]][1:(k^2*(x+1)),][-which(!c(upper.tri(matrix(0, nrow=k, ncol=k), diag=TRUE)) == TRUE)])))}) - 2*(seq(3, (k*(k+1)/2 + L*k^2), by=2*k))
      
    }
    
    data.frame(test = c("Corr Jin 2016", "Jin 2019"),
               stat = c(max(corr_results),  max(unlist(Sample_U_r_l))),
               crit_vals = c(quantile(S_r_l_dep, 1-alpha), quantile(S_r_l, 1-alpha)),
               pval = c((sum(S_r_l_dep > max(unlist(corr_results))) + 1)/(B+1), (sum(S_r_l > max(unlist(Sample_U_r_l))) + 1)/(B+1)),
               L.hat = c(which.max(corr_results)-1, which.max(unlist(Sample_U_r_l))-1) )
    
  }
  
  
  ## Tests for calculateAutocovariance function
  
  ## generate two univariate time series
  set.seed(1234)
  n <- 100
  X1 <- arima.sim(model = list(ar = 0.8), n)
  Y1 <- arima.sim(model = list(ma = 0.2), n)
  
  # My function output
  out1 <- calculateAutocovariance(as.matrix(X1), as.matrix(Y1), 6)
  
  # Correct base R output
  acfXX <- ccf(X1, X1, plot = FALSE, type = "covariance", demean = TRUE, lag.max = 5)$acf
  ccfXY <- ccf(X1, Y1, plot = FALSE, type = "covariance", demean = TRUE, lag.max = 5)$acf
  ccfYX <- ccf(Y1, X1, plot = FALSE, type = "covariance", demean = TRUE, lag.max = 5)$acf
  acfYY <- ccf(Y1, Y1, plot = FALSE, type = "covariance", demean = TRUE, lag.max = 5)$acf
  
  ## Try a multivariate example
  X <- cbind(rnorm(n), rnorm(n), rnorm(n))
  Y <- cbind(rnorm(n), rnorm(n), rnorm(n))
  
  # My function output
  out2 <- calculateAutocovariance(X, Y, 6)[ , , 6:11]
  out2 <- aperm(out2, c(3, 1, 2))
  
  # Correct base R output
  acfZ <- acf(cbind(X, Y), plot = FALSE, type = "covariance", demean = TRUE, lag.max = 5)$acf
  
  
  # Univariate time series example
  expect_equal(matrix(out1[1, 1, ]), matrix(acfXX))
  expect_equal(matrix(out1[1, 2, ]), matrix(ccfXY))
  expect_equal(matrix(out1[2, 1, ]), matrix(ccfYX))
  expect_equal(matrix(out1[2, 2, ]), matrix(acfYY))
  # Multivariate example
  expect_equal(out2, acfZ)
  
  
  ## Tests for dependentCovariance function
  set.seed(1234)
  n <- 1000
  X1 <- arima.sim(list(ar = 0.8), n)
  X2 <- arima.sim(list(ar = 0.8), n)
  Y1 <- arima.sim(list(ar = 0.8), n)
  Y2 <- arima.sim(list(ar = 0.8), n)
  
  # Try a multivariate example
  out1 <- get.covar.matrix(cbind(X1, X2, Y1, Y2), 5, corr_only = TRUE)
  out2 <- dependentCovariance(cbind(X1, X2), cbind(Y1, Y2), 5)
  
  # Try a univariate example
  out3 <- get.covar.matrix(cbind(X1, Y1), 5, corr_only = TRUE)
  out4 <- dependentCovariance(matrix(X1), matrix(Y1), 5)
  
  # Multivariate
  expect_equal(out1[[1]][-2], out2[[1]], check.attributes = FALSE, tolerance = 0.0000001)
  expect_equal(as.numeric(out1[[2]][-2, -2]), as.numeric(out2[[2]]), check.attributes = FALSE, tolerance = 0.00000001)
  # Univariate
  expect_equal(as.numeric(out3[[1]]), as.numeric(out4[[1]]), check.attributes = FALSE)
  expect_equal(as.numeric(out3[[2]]), as.numeric(out4[[2]]), check.attributes = FALSE)
  
  
  microbenchmark(
    get.covar.matrix(cbind(X1, X2, Y1, Y2), 5, corr_only = TRUE),
    dependentCovariance(cbind(X1, X2), cbind(Y1, Y2), 5),
    times = 10
  )
  # 80 times faster
  
})
