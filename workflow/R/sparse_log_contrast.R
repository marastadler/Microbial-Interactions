## adapted functions for sparse log contrast regression
 


## make sure to set slc$refit <- TRUE


## This function is taken from https://github.com/viettr/trac/blob/master/R/cv_log_contrast.R
## and is extended to perform refitting in each fold


#' Perform cross validation for tuning parameter selection for sparse log contrast
#'
#' This function is to be called after calling \code{\link{sparse_log_contrast}}.
#' It performs \code{nfold}-fold cross validation.
#'
#' @param fit output of \code{\link{sparse_log_contrast}} function.
#' @param Z,y,additional_covariates same arguments as passed to
#'   \code{\link{sparse_log_contrast}}. C will be taken from fit object.
#' @param folds a partition of \code{1:nrow(Z)}.
#' @param nfolds number of folds for cross-validation
#' @param summary_function how to combine the errors calculated on each
#' observation within a fold (e.g. mean or median)
#' @export
#' 
cv_sparse_log_contrast_refit <- function(fit, Z, y, folds = NULL, nfolds = 5,
                                   summary_function = stats::median,
                                   additional_covariates = NULL) {
  
  n <- nrow(Z)
  p <- ncol(Z)
  stopifnot(length(y) == n)
  if (!is.null(additional_covariates) & !is.data.frame(additional_covariates)) {
    additional_covariates <- data.frame(additional_covariates)
  }
  if (is.null(folds)) folds <- ggb:::make_folds(n, nfolds)
  else
    nfolds <- length(folds)
  cv <- list()
  fit_folds <- list() # save this to reuse by log-ratio's cv function
  errs <- matrix(NA, ncol(fit$beta), nfolds)
  for (i in seq(nfolds)) {
    cat("fold", i, fill = TRUE)
    # add for backward compatibility
    if (is.null(fit$method)) fit$method <- "regr"
    if (is.null(fit$w_additional_covariates)) {
      fit$w_additional_covariates <- NULL
    }
    if (is.null(fit$rho)) fit$rho <- 0
    if (is.null(fit$normalized)) fit$normalized <- FALSE
    # train on all but i-th fold (and use settings from fit):
    fit_folds[[i]] <- sparse_log_contrast(Z[-folds[[i]], ],
                                          y[-folds[[i]]],
                                          additional_covariates[-folds[[i]], ],
                                          fit$C,
                                          fraclist = fit$fraclist,
                                          w_additional_covariates =
                                            fit$w_additional_covariates,
                                          method = fit$method,
                                          rho = fit$rho,
                                          normalized = fit$normalized)
    # if (fit$refit) stop("Not yet supported.")
    
    if (fit$refit) {
      fit_folds[[i]]$refit <- FALSE#TRUE
      ## for refitting the input matrix contains the compositional and also the additional covariates
      
      if(!is.null(additional_covariates)){
        tryCatch({
        fit_folds[[i]] <- trac::refit_sparse_log_contrast_classo(fit = fit_folds[[i]], Z = Z[-folds[[i]], ],
                                                    y = y[-folds[[i]]],
                                                    additional_covariates = additional_covariates[-folds[[i]], ])
      }, error=function(e){})
        # fit_folds[[i]] <- refit_sparse_log_contrast_old(fit = fit_folds[[i]], Z = Z[-folds[[i]], ],
        #                                                          y = y[-folds[[i]]],
        #                                                          additional_covariates = additional_covariates[-folds[[i]], ])

        
      } else{
        tryCatch({
        fit_folds[[i]] <- trac::refit_sparse_log_contrast_classo(fit_folds[[i]], Z[-folds[[i]], ],
                                                    y[-folds[[i]]],
                                                    additional_covariates = NULL)
        }, error=function(e){})
        # fit_folds[[i]] <- refit_sparse_log_contrast_old(fit = fit_folds[[i]], Z = Z[-folds[[i]], ],
        #                                                          y = y[-folds[[i]]],
        #                                                          additional_covariates = NULL)
        }
      
      
    }
    
    
    
    if (fit$method == "regr" | is.null(fit$method)) {
      errs[, i] <- apply((predict_trac(
        list(fit_folds[[i]]),
        Z[folds[[i]], ],
        additional_covariates[folds[[i]], ])[[1]] - y[folds[[i]]])^2,
        2, summary_function
      )
    }
    
    if (fit$method == "classif" |
        fit$method == "classif_huber") {
      # loss: max(0, 1 - y_hat * y)^2
      er <- sign(predict_trac(list(fit_folds[[i]]),
                              Z[folds[[i]],],
                              additional_covariates[folds[[i]],])[[1]]) !=
        c(y[folds[[i]]])
      errs[, i] <- colMeans(er)
    }
  }
  m <- rowMeans(errs)
  se <- apply(errs, 1, stats::sd) / sqrt(nfolds)
  ibest <- which.min(m)
  i1se <- min(which(m < m[ibest] + se[ibest]))
  cv <- list(errs = errs, m = m, se = se,
             lambda_best = fit$fraclist[ibest], ibest = ibest,
             lambda_1se = fit$fraclist[i1se], i1se = i1se,
             fraclist = fit$fraclist,
             nonzeros = colSums(abs(fit$beta) > 1e-5),
             fit_folds = fit_folds)
  list(cv = cv, folds = folds)
}


### This function is copied from https://github.com/viettr/trac/blob/prediction_during_cv/R/refit_log_contrast.R

#' Refit log-contrast regression to sparsity constraint
#'
#' Given output of \code{\link{sparse_log_contrast}}, solves the regression
#' problem with compositional constraint on the features selected by
#' sparse_log_contrast. In contrast to \code{\link{refit_sparse_log_contrast}}
#' this function only does the refit for one model specified by i_selected
#' or component_selected.
#' This i_selected usually comes from the i1se from the cross-validation output.
#'
#'
#' @param fit output of sparse_log_contrast
#' @param i_selected indicator which lambda is selected based on for example the
#'   cross-validation procedure
#' @param Z,y,additional_covariates same arguments as passed to
#'   \code{\link{sparse_log_contrast}}
#' @param tol tolerance for deciding whether a beta value is zero
#' @param component_selected vector with indices which component to include
#' @export
refit_sparse_log_contrast_reg <- function(fit, i_selected = NULL, Z, y,
                                          additional_covariates = NULL,
                                          tol = 1e-5,
                                          component_selected = NULL) {
  # check which components have non-zero coefficients
  if (is.null(component_selected)) {
    selected_variables <- abs(fit$beta[, i_selected]) > tol
    names_variables <- names(selected_variables)[which(selected_variables)]
    if (sum(selected_variables) < 1) stop("There are no selected components")
  } else {
    p_test <- ifelse(is.null(additional_covariates),  ncol(Z),
                     (ncol(Z) + ncol(additional_covariates)))
    if (max(component_selected) > p_test) {
      stop("The values of component_selected must be smaller than ncol(Z)
           + ncol(additional_covariates")
    }
    if (min(component_selected) < 1) {
      stop("The values of component_selected must be greater than 1")
    }
    selected_variables <- component_selected
    names_variables <- rownames(fit$beta)[which(selected_variables)]
  }
  if (is.null(names_variables)) {
    names_variables <- 1:length(which(selected_variables))
  }
  # for covariates check which of them need to satisfy the zero sum constraint
  C <- fit$C[, selected_variables]
  C <- matrix(C, nrow = 1)
  # extract meta information from the fit object
  rho <- fit$rho
  normalized <- fit$normalized
  intercept <- TRUE
  
  # Transformation of the additional covariates
  if (!is.null(additional_covariates)) {
    # transform additional covariates to data.frame --> easier to work with
    # different types of input
    if (!is.data.frame(additional_covariates)) {
      additional_covariates <- data.frame(additional_covariates)
    }
    # define the number of additional covariates
    p_x <- ncol(additional_covariates)
  }
  
  # normalize the non-compositional data if wanted
  if (!is.null(additional_covariates)) {
    if (normalized) {
      # call the normalization helper function
      normalized_values <-
        normalization_additional_covariates(additional_covariates =
                                              additional_covariates,
                                            p_x = p_x,
                                            intercept = intercept)
      additional_covariates <- normalized_values$X
    } else {
      # get the number of categorical variables if no normalization is applied
      categorical_list <- get_categorical_variables(additional_covariates)
      categorical <- categorical_list[["categorical"]]
      n_categorical <- categorical_list[["n_categorical"]]
      if (n_categorical > 0) {
        additional_covariates[, categorical] <-
          transform_categorical_variables(additional_covariates, categorical)
      }
    }
  }
  
  
  # Concenate the non-compositional data if available
  if (!is.null(additional_covariates)) {
    # add the non-compositional covariates
    X_classo <- as.matrix(cbind(Z, additional_covariates))
  } else {
    X_classo <- as.matrix(Z)
  }
  #  ybar <- mean(y)
  #  yt <- y - ybar
  X_classo <- X_classo[, selected_variables]
  
  n <- nrow(X_classo)
  p <- ncol(X_classo)
  
  # set up CLASSO problem:
  prob <- classo$classo_problem(X = X_classo,
                                C = C,
                                y = array(y))
  prob$formulation$classification <- FALSE
  prob$formulation$concomitant <- FALSE
  prob$formulation$intercept <- TRUE
  prob$formulation$huber <- FALSE
  prob$model_selection$PATH <- FALSE
  prob$model_selection$CV <- FALSE
  prob$model_selection$StabSel <- FALSE
  prob$model_selection$LAMfixed <- TRUE
  prob$model_selection$LAMfixedparameters$rescaled_lam <- TRUE
  prob$model_selection$LAMfixedparameters$lam <- 0.0
  
  
  
  # solve it
  prob$solve()
  # extract outputs
  beta <- (prob$solution$LAMfixed$beta)
  if (intercept) {
    # c-lasso can estimate beta0 --> select first column (estimated beta0)
    # delete the first column afterwards
    beta0 <- beta[1]
    beta <- beta[-1]
  }
  rownames(beta) <- names_variables
  beta <- t(beta)
  
  list(beta0 = beta0,
       beta = beta,
       C = C,
       refit = TRUE,
       method = fit$method,
       intercept = intercept,
       rho = rho,
       normalized = normalized)
}



#' Refit subject to sparsity constraints
#'
#' Given output of \code{\link{sparse_log_contrast}}, solves the least squares problem with
#' compositional constraint on the features selected by sparse_log_contrast.
#'
#' minimize_{beta, beta0} 1/(2n) || y - beta0 1_n - Z beta ||^2
#' subject to beta_{nonselected} = 0, 1_p^T beta = 0
#'
#' @param fit output of sparse_log_contrast
#' @param Z,y same arguments as passed to \code{\link{sparse_log_contrast}}
#' @param tol tolerance for deciding whether a beta value is zero
#' @export
refit_sparse_log_contrast_old <- function(fit, Z, y, additional_covariates = NULL, tol = 1e-4) {
  
  n <- nrow(Z)
  p1 <- ncol(Z)
  
  nlam <- length(fit$fraclist)
  
  if(!is.null(additional_covariates)){
    Z <- as.matrix(cbind(Z, additional_covariates))
    p2 <- ncol(additional_covariates)
    p <- p1 + p2
  }
  else { p <- p1 }
  
  for (i in seq(nlam)) {
    # let nz be nonzero elements from sparse log contrast and z be the zero ones:
    z <- which(abs(fit$beta[1:p1, i]) <= tol)
    if(!is.null(additional_covariates)){
      zmain <- z
      zint <- which(abs(fit$beta[, i]) <= tol)
      zint <- zint[!(zint %in% zmain)]
      z <- c(zmain, zint)
    }
    
    if (length(z) == p) {
      # all zero
      fit$beta0[i] <- fit$beta0[i]
      fit$beta[, i] <- fit$beta[, i]
      next
    }
    
    if(length(z) < p){
      # minimize_{beta0, beta} 1/(2n) || y - beta0 1_n - Z beta ||^2
      # subject to beta_z = 0, 1_p^T beta = 0
      # can be written as
      # minimize_{beta0, beta} 1/(2n) || y - beta0 1_n - Z_nz beta_nz ||^2
      # subject to beta_z = 0, 1_nz^T beta_nz = 0
      
      if(!is.null(additional_covariates)) {
        v <- c(rep(1, p1 - length(zmain)), numeric(p2 - length(zint)))
      } else { v <- rep(1, p1 - length(z))}
      
      # the above is equiv to solving
      # minimize_{beta0, beta_nz} 1/(2n) || y - beta0 1_n - Z_nz (I - P)beta_nz ||^2
      # where P = vv^T / ||v||^2
      # or
      # minimize_{beta0, beta_nz} 1/(2n) || y - M [beta0; beta_nz] ||^2
      # where M = [1_n, Z_nz (I - P)]
      # or [beta0; beta_nz] = M^{+} y
      if (length(z) == 0)
        Zz <- Z
      else
        Zz <- Z[, -z]
      
      ## special case where v consists out of zeros only (all main effects are selected): 
      if(all(v == 0)) {  M <- cbind(1, Zz) }
      else {
        ZzP <- (Zz %*% v) %*% matrix(v, nrow = 1) / sum(v^2)
        M <- cbind(1, Zz - ZzP)
      }
      
    
      sv <- svd(M)
      # M^{+} y = V D^{+} U^T y
      solution <- sv$v %*% (c(1 / sv$d[-length(sv$d)], 0) * crossprod(sv$u, y))
      fit$beta[z, i] <-  0
      if (length(z) == 0)
        fit$beta[, i] <- solution[-1]
      else
        fit$beta[-z, i] <- solution[-1]
      fit$beta0[i] <- solution[1]
    }
  }
  fit$tol <- tol
  fit$refit <- TRUE
  fit
}
