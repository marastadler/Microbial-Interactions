
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("R/slc_int.R")
source("R/utils.R")
library("trac")
## Viet's implementation:
classo_slc <- function(x, y, p_comp = NULL, q, ...) {
  #' wrapper around slc fitting to find stable features, attention select q
  #'  higher because we select the highest values across the lambda path and
  #'  the value q corresponds to the q largest values not the q covariates
  #'  with the highest coefs
  #' @param x (log transformed compositional) input
  #'    (if additional covariates concatinate and specify parameter p.
  #'    The first p values are assumed to be compositional)
  #' @param y outcome
  #' @param p how many of the parameters are compositional. Default NULL
  #' @param q maximum number of non zero coefficient to consider
  # fit model
  if (!is.null(p_comp)) {
    X_add <- x[, (p_comp + 1):ncol(x)]
    Z <- x[, 1:p_comp]
  } else {
    X_add <- NULL
    Z <- x
  }
  fit <- trac::sparse_log_contrast(Z = Z, y = y,
                                   additional_covariates = X_add,
                                   ...)
  # which coefficients are non-zero?
  # compute selection paths
  cf <- fit$beta

  #
  which_largest <- arrayInd(order(abs(cf), decreasing = TRUE)[1:(q)], .dim = dim(cf))
  which_largest <- unique(which_largest[, 1])
  selected <- logical(ncol(x))
  selected[which_largest] <- TRUE
  names(selected) <- colnames(x)

  ## return both
  return(list(selected = selected, path = cf))
}
# # # 11:41
# ### test slc
# set.seed(123)
# comp_x <- matrix(rnorm(100000), nrow = 100, ncol = 1000)
# colnames(comp_x) <- paste0("test", 1:1000)
# comp_x <- abs(comp_x) + 1
# comp_x_log <- log(comp_x)
# comp_x_clr <- compositions::clr(comp_x)
# additional_x <- matrix(rnorm(1000), nrow = 100, ncol = 10)
# colnames(additional_x) <- paste0("Add", 1:10)
# x <- cbind(comp_x_log, additional_x)
# fit <- classo_slc(x = x,
#                   y =  100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + 100 * additional_x[, 10],
#                   p_comp = 1000,
#                   q = 50)
# stab.glmnet <- stabs::stabsel(x = x,
#                        y = 100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + 100 * additional_x[, 10],
#                        fitfun = classo_slc, cutoff = 0.75,
#                        PFER = 1, args.fitfun = list(p_comp = 1000))
# stab.glmnet$selected

# My own extension to SLC + interactions

classo_slc_int <- function(x, y, q, ...) {
  #' wrapper around slc + int fitting to find stable features, attention select q
  #'  higher because we select the highest values across the lambda path and
  #'  the value q corresponds to the q largest values not the q covariates
  #'  with the highest coefs
  #' @param x (log transformed compositional) input
  #'    (if additional covariates concatinate and specify parameter p.
  #'    The first p values are assumed to be compositional)
  #' @param y outcome
  #' @param q maximum number of non zero coefficient to consider
  # fit model



  fit <- trac::sparse_log_contrast(Z = Z, y = y,
                                   additional_covariates = X_add,
                                   ...)
  # which coefficients are non-zero?
  # compute selection paths
  cf <- fit$beta

  #
  which_largest <- arrayInd(order(abs(cf), decreasing = TRUE)[1:(q)], .dim = dim(cf))
  which_largest <- unique(which_largest[, 1])
  selected <- logical(ncol(x))
  selected[which_largest] <- TRUE
  names(selected) <- colnames(x)

  ## return both
  return(list(selected = selected, path = cf))
}






slc.stabsel <- function(x, y, method = "regr", output = "raw",
                        quantile_cutoff = NULL, k_largest = NULL, q_max = 20){


  if(!is.null(quantile_cutoff) & !is.null(k_largest)){
    stop("Either provide 'quantile_cutoff' or 'k_largest'! Not both.")
  }
  if(is.null(quantile_cutoff) & is.null(k_largest)){
    stop("Either provide 'quantile_cutoff' or 'k_largest'!")
  }
  if(!is.null(k_largest)){
    if(k_largest > ncol(x)){
      stop("k_largest can't be greater than p!")
    }
  }


  beta_est_refit <- c()

  if(method == "regr") {
    refit_selected_model <- trac::refit_sparse_log_contrast_reg
  }
  if(method == "classif") {
    refit_selected_model <- trac::refit_sparse_log_contrast_classif
  }


  z <- log(x)

  fit <- trac::sparse_log_contrast(z, y, min_frac = 1e-3, nlam = 30, method = method)
  fit$refit <- FALSE#TRUE


  # which coefficients are non-zero?
  # compute selection paths
  cf <- fit$beta
  cf_rowsum <- rowSums(cf)
  if(!is.null(quantile_cutoff)){
    quant <- quantile(abs(cf_rowsum), probs = quantile_cutoff)
    cf_max <- cf_rowsum
    cf_max[abs(cf_max) < quant] <- 0
    cf_max[abs(cf_max) >= quant] <- 1
    ind_selected <- which(cf_max == 1)
    # which_largest <- arrayInd(order(abs(cf), decreasing = TRUE)[1:(q_max)], .dim = dim(cf))
    # which_largest <- unique(which_largest[, 1])
    selected <- cf_max
  }
  if(!is.null(k_largest)){
    cf_max <- cf_rowsum
    cf_max[order(cf_rowsum, decreasing = T)[1:k_largest]] <- 1
    cf_max[order(cf_rowsum, decreasing = T)[(k_largest + 1):length(cf_rowsum)]] <- 0
    ind_selected <- which(cf_max == 1)
    selected <- cf_max

  }

  names(selected) <- colnames(x)

  ## return both
  return(list(selected = selected, path = cf))




}






#' slc_int sparse log-contrast model with interactions (no train-test split included!) for Stabsel
#'  @normalize.int logical; whether interactions should be normalized
#'
slc_int.stabsel <- function(x, y, q_max = 60, method = "regr", output = "raw",
                            normalize.int = TRUE, p_main, quantile_cutoff = NULL, k_largest = NULL, ...){

  if(!is.null(quantile_cutoff) & !is.null(k_largest)){
    stop("Either provide 'quantile_cutoff' or 'k_largest'! Not both.")
  }
  if(is.null(quantile_cutoff) & is.null(k_largest)){
    stop("Either provide 'quantile_cutoff' or 'k_largest'!")
  }

  if(!is.null(k_largest)){
    if(k_largest > ncol(x)){
      stop("k_largest can't be greater than p!")
    }
  }

  X <- x
  # apply slc with interactions on 10 random train-test splits...
  x <- X[, 1:p_main]
  ## calculate interactions based on x with pseudo counts
  x_int <- X[, (p_main + 1):ncol(X)]


  beta_int_est_refit <- c()

  if(method == "regr") {
    refit_selected_model <- trac::refit_sparse_log_contrast_reg
  }
  if(method == "classif") {
    refit_selected_model <- trac::refit_sparse_log_contrast_classif
  }


  z <- log(x)




  if(normalize.int) {
    ### normalize interaction effects
    ## 1. clr transform main effects
    x_clr <- clr.trafo(x)
    ## 2. calculate norm of clr transformed main effects and interactions
    x_clr_norm <- apply(x_clr, 2, norm_vec)
    x_int_norm <- apply(x_int, 2, norm_vec)

    ## normalization
    avg_norm_clr_x <- mean(x_clr_norm)
    x_int_norm <- apply(x_int, 2, norm_vec)
    x_int_normalization <- t(t(x_int)/(x_int_norm)) *  avg_norm_clr_x
    # x_int_normalization_norm <- apply(x_int_normalization, 2 , norm_vec)
    x_int_normalization[is.na(x_int_normalization)] <- 0
  }
  else {
    x_int_normalization <- x_int
  }

  fit_int <- trac::sparse_log_contrast(z, y, additional_covariates = x_int_normalization,
                                 normalized = F, ## normalization already performed
                                 min_frac = 1e-3, nlam = 30, method = method)
  fit_int$refit <- FALSE#TRUE

  if(normalize.int) { ## back calculation after normalization
    fit_int$beta[(ncol(z) + 1): nrow(fit_int$beta), ] <-  fit_int$beta[(ncol(z) + 1): nrow(fit_int$beta), ] * avg_norm_clr_x / x_int_norm
    fit_int$beta[is.na(fit_int$beta)] <- 0
  }

  # which coefficients are non-zero?
  # compute selection paths
  cf <- fit_int$beta
  # matplot(x = log(fit_int$fraclist), t(as.matrix(cf)), type = "l")

  cf_rowsum <- rowSums(cf)
  if(!is.null(quantile_cutoff)){
    quant <- quantile(abs(cf_rowsum), probs = quantile_cutoff)
    cf_max <- cf_rowsum
    cf_max[abs(cf_max) < quant] <- 0
    cf_max[abs(cf_max) >= quant] <- 1
    ind_selected <- which(cf_max == 1)
    # which_largest <- arrayInd(order(abs(cf), decreasing = TRUE)[1:(q_max)], .dim = dim(cf))
    # which_largest <- unique(which_largest[, 1])
    selected <- cf_max
  }
  if(!is.null(k_largest)){
    cf_max <- cf_rowsum
    cf_max[order(cf_rowsum, decreasing = T)[1:k_largest]] <- 1
    cf_max[order(cf_rowsum, decreasing = T)[(k_largest + 1):length(cf_rowsum)]] <- 0
    ind_selected <- which(cf_max == 1)
    selected <- cf_max

  }
  names(selected) <- c(colnames(x), colnames(x_int))

  ## return both
  return(list(selected = selected, path = cf))




}


### test slc Stabsel
set.seed(123)
comp_x <- matrix(rnorm(2000), nrow = 100, ncol = 20)
colnames(comp_x) <- paste0("test", 1:20)
comp_x <- abs(comp_x) + 1
comp_x <- comp_x/sum(comp_x)
comp_x_log <- log(comp_x)
comp_x_clr <- compositions::clr(comp_x)
x <- comp_x



fit <- slc.stabsel(x = x[1:30,],
                  y =  100 * comp_x_clr[1:30, 1] - 100 * comp_x_clr[1:30, 2] + rnorm(30) * 0.3,
                  quantile_cutoff = .9, q_max =1)
## q_max doesn't have an influence, but stabsel needs it to not return an error...
which(fit$selected == TRUE)

stab.glmnet <- stabs::stabsel(x = x,
                              y =  100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + rnorm(100) * 0.3,
                              fitfun =  slc.stabsel, cutoff = 0.8, PFER = 1, args.fitfun = list(quantile_cutoff = .9),
                              sampling.type = "SS", B = 10)

stab.glmnet$selected
barplot(stab.glmnet$max[stab.glmnet$max > 0])


stab.glmnet <- stabs::stabsel(x = x,
                              y =  100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + rnorm(100) * 0.3,
                              fitfun =  slc.stabsel, cutoff = 0.6, PFER = 1,
                              args.fitfun = list(k_largest = 3, quantile_cutoff = NULL, method = "regr"),
                              sampling.type = "SS", B = 10)

stab.glmnet$selected

# ### test slc + int Stabsel
set.seed(123)
comp_x <- matrix(rnorm(2000), nrow = 100, ncol = 20)
colnames(comp_x) <- paste0("test", 1:20)
comp_x <- abs(comp_x) + 1
comp_x <- comp_x/sum(comp_x)
comp_x_log <- log(comp_x)
comp_x_clr <- compositions::clr(comp_x)
x <- comp_x

x_interact <- compute.interactions.aitchison(comp_x)


fit <- slc_int.stabsel(x = cbind(x, x_interact)[1:30,],
                  y =  100 * comp_x_clr[1:30, 1] - 100 * comp_x_clr[1:30, 2] + 100 * x_interact[1:30, "test1:test5"] + rnorm(30) * 0.3,
                  q_max = 70, p_main = 20, quantile_cutoff = .98)

which(fit$selected == TRUE)





stab.glmnet <- stabs::stabsel(x = cbind(x, x_interact),
                       y =  100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + 100 * x_interact[, "test1:test5"]  + rnorm(100) * 0.3,
                       fitfun =  slc_int.stabsel, cutoff = 0.7, PFER = 1, args.fitfun = list(p_main = 20, quantile_cutoff = .85),
                       sampling.type = "SS", B = 20)

stab.glmnet$selected
barplot(stab.glmnet$max[stab.glmnet$max > 0])


stab.glmnet2 <- stabs::stabsel(x = cbind(x, x_interact),
                              y =  100 * comp_x_clr[, 1] - 100 * comp_x_clr[, 2] + 100 * x_interact[, "test1:test5"]  + rnorm(100) * 0.3,
                              fitfun =  slc_int.stabsel, cutoff = 0.7, PFER = 1, args.fitfun = list(p_main = 20, k_largest = 3,
                                                                                                    quantile_cutoff = NULL, method = "regr"),
                              sampling.type = "SS", B = 20)

stab.glmnet2$selected
