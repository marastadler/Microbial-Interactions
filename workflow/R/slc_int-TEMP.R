

compute.interactions.aitchison2 <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x2 <- matrix(nrow = n, ncol = p*(p + 1)/2 - p)
  ind <- 0
  col.names <- c()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      ind <- ind + 1
      x2[, ind] <- (x[, i] - x[, j])^2#(log(x[, i]) - log(x[, j]))^2
      col.names[ind] <- paste0(colnames(x)[i], ":", colnames(x)[j])
    }
  }
  colnames(x2) <- col.names
  return(x2)
}



#'  @normalize.int logical; whether interactions should be normalized
#'
sparse_log_contrast_interaction <- function(X, y, method = "classif",
                                            normalize.int = TRUE){


  beta_int_est_refit <- c()
  beta_int_est_norefit <- c()



  if(method == "regr") {
    refit_selected_model <- trac::refit_sparse_log_contrast_reg
  }
  if(method == "classif") {
    refit_selected_model <- trac::refit_sparse_log_contrast_classif
  }



  z <- log(X)


  ## calculate interactions based on X with pseudo counts

  x_int <- compute.interactions.aitchison(X)



  if(normalize.int) {
    ### normalize interaction effects
    ## 1. clr transform main effects
    x_clr <- clr.trafo(X)

    ## 2. calculate norm of clr transformed main effects and interactions
    x_clr_norm <- apply(x_clr, 2, norm_vec)
    x_int_norm <- apply(x_int, 2, norm_vec)


    ## normalization
    avg_norm_clr_X <- mean(x_clr_norm)
    x_int_norm <- apply(x_int, 2, norm_vec)
    x_int_normalization <- t(t(x_int)/(x_int_norm)) *  avg_norm_clr_X
    # x_int_normalization_norm <- apply(x_int_normalization, 2 , norm_vec)
    x_int_normalization[is.na(x_int_normalization)] <- 0

  }
  else {
    x_int_normalization <- x_int
  }

  fit_int <- sparse_log_contrast(z, y,
                                 additional_covariates =
                                   x_int_normalization,
                                 normalized = F,
                                 min_frac = 1e-3, nlam = 30, method = method)
  fit_int$refit <- FALSE#TRUE


  ## back calculation interaction coefficients beta
  if(normalize.int) {
    fit_int$beta[(ncol(z) + 1):nrow(beta_int_est_refit), ] <- fit_int$beta[(ncol(z) + 1):nrow(beta_int_est_refit), ] * avg_norm_clr_X / x_int_norm
    fit_int$beta[is.na(fit_int$beta)] <- 0
  }


  return(fit_int)

}



#'  @normalize.int logical; whether interactions should be normalized
#'
slc_slc_int_all_splits <- function(X, y, method = "classif", output = "class", nbsplit = 10,
                                   normalize.int = TRUE){


  # apply slc with interactions on 10 random train-test splits...
  set.seed(123)
  nsplit <- nbsplit
  ntot <- length(y)
  p <- ncol(X)
  n <- round(2/3 * ntot)

  tr <- list()
  fit <- list()
  cvfit <- list()
  yhat_tr <- list()
  yhat_te <- list()
  trainerr <- list()
  testerr <- list()
  nnz <- list()
  beta_est_refit <- matrix(nrow = ncol(X), ncol = nsplit)

  fit_int <- list()
  cvfit_int <- list()
  yhat_tr_int <- list()
  yhat_te_int <- list()
  trainerr_int <- list()
  testerr_int <- list()
  nnz_int <- list()
  beta_int_est_refit <- matrix(nrow = ncol(X) * (ncol(X) + 1)/2, ncol = nsplit)
  beta_int_est_norefit <- matrix(nrow = ncol(X) * (ncol(X) + 1)/2, ncol = nsplit)


  # ## CLR transformation of X
  # Z <- log(X)
  # Z_bar <- Matrix::rowMeans(Z)
  # X_CLR <- Z - Z_bar
  # colnames(X_CLR) <- 1:ncol(X_CLR)
  #

  for (r in seq(nsplit)) {
    cat("split ", r, fill = TRUE)

    # if(method == "classif") {
    #   ## make sure to draw equally often from both groups
    #   g1 <- which(y == 1)
    #   g2 <- which(y == -1)
    #   n1 <- round(2/3 * length(g1))
    #   n2 <- round(2/3 * length(g2))
    #   set.seed(r)
    #   tr1 <- sample(g1, n1)
    #   set.seed(r)
    #   tr2 <- sample(g2, n2)
    #   tr[[r]] <- c(tr1, tr2)
    # }
    # else {
    set.seed(r)
    tr[[r]] <- sample(ntot, n)
    # }


    fit[[r]] <- list()
    cvfit[[r]] <- list()
    yhat_tr[[r]] <- list()
    yhat_te[[r]] <- list()
    trainerr[[r]] <- list()
    testerr[[r]] <- list()
    nnz[[r]] <- list()

    if(method == "regr") {
      refit_selected_model <- trac::refit_sparse_log_contrast_reg
    }
    if(method == "classif") {
      refit_selected_model <- trac::refit_sparse_log_contrast_classif
    }

    ytr <- y[tr[[r]]]
    yte <- y[-tr[[r]]]

    ztr <- log(X[tr[[r]], ])
    zte <- log(X[-tr[[r]], ])

    fit[[r]] <- sparse_log_contrast(ztr, ytr, min_frac = 1e-3, nlam = 30, method = method)
    fit[[r]]$refit <- FALSE#TRUE
    cvfit[[r]] <- #cv_sparse_log_contrast_refit(fit[[r]], Z = ztr, y = ytr)
      cv_sparse_log_contrast(fit[[r]], Z = ztr, y = ytr)
    yhat_tr[[r]] <- predict_trac(list(fit[[r]]), new_Z = ztr, output = output)[[1]]
    yhat_te[[r]] <- predict_trac(list(fit[[r]]), new_Z = zte, output = output)[[1]]
    trainerr[[r]] <- sqrt(colMeans((yhat_tr[[r]] - ytr)^2))
    testerr[[r]] <- sqrt(colMeans((yhat_te[[r]] - yte)^2))
    nnz[[r]] <- colSums(abs(fit[[r]]$beta) > 1e-5 )



    ## refit if model is not empty
    if(nnz[[r]][cvfit[[r]]$cv$i1se] != 0){
      refit_lc <- refit_selected_model(fit[[r]], i_selected = cvfit[[r]]$cv$i1se,
                                                      Z = ztr, y = ytr)
      rownames(beta_est_refit) <- colnames(X)[1:p]
      beta_est_refit[colnames(refit_lc$beta), r] <- refit_lc$beta
    }
    beta_est_refit[is.na(beta_est_refit)] <- 0





    ## Model with interactions
    fit_int[[r]] <- list()
    cvfit_int[[r]] <- list()
    yhat_tr_int[[r]] <- list()
    yhat_te_int[[r]] <- list()
    trainerr_int[[r]] <- list()
    testerr_int[[r]] <- list()
    nnz_int[[r]] <- list()



    ## calculate interactions based on X with pseudo counts
    # xtr_int <- compute.interactions.aitchison2(X_CLR[tr[[r]], ])
    # xte_int <- compute.interactions.aitchison2(X_CLR[-tr[[r]], ])
    xtr_int <- compute.interactions.aitchison(X[tr[[r]], ])
    xte_int <- compute.interactions.aitchison(X[-tr[[r]], ])



    if(normalize.int) {
      ### normalize interaction effects
      ## 1. clr transform main effects
      xtr_clr <- clr.trafo(X[tr[[r]], ])
      xte_clr <- clr.trafo(X[-tr[[r]], ])

      ## 2. calculate norm of clr transformed main effects and interactions
      xtr_clr_norm <- apply(xtr_clr, 2, norm_vec)
      xtr_int_norm <- apply(xtr_int, 2, norm_vec)

      xte_clr_norm <- apply(xte_clr, 2, norm_vec)
      xte_int_norm <- apply(xte_int, 2, norm_vec)

      ## training set normalization
      avg_norm_clr_Xtr <- mean(xtr_clr_norm)
      xtr_int_norm <- apply(xtr_int, 2, norm_vec)
      xtr_int_normalization <- t(t(xtr_int)/(xtr_int_norm)) *  avg_norm_clr_Xtr
      # xtr_int_normalization_norm <- apply(xtr_int_normalization, 2 , norm_vec)
      xtr_int_normalization[is.na(xtr_int_normalization)] <- 0
      ## test set normalization wrt. main effects in training data
      # avg_norm_clr_Xte <- mean(xte_clr_norm)
      xte_int_norm <- apply(xte_int, 2, norm_vec)
      xte_int_normalization <- t(t(xte_int)/(xte_int_norm)) *  avg_norm_clr_Xtr
      # xte_int_normalization_norm <- apply(xte_int_normalization, 2 , norm_vec)
      xte_int_normalization[is.na(xte_int_normalization)] <- 0
    }
    else {
      xtr_int_normalization <- xtr_int
      xte_int_normalization <- xte_int
    }

    fit_int[[r]] <- sparse_log_contrast(ztr, ytr,
                                        additional_covariates =
                                          xtr_int_normalization,
                                        normalized = F,
                                        min_frac = 1e-3, nlam = 30, method = method)
    fit_int[[r]]$refit <- FALSE#TRUE


    cvfit_int[[r]] <- #cv_sparse_log_contrast_refit
      cv_sparse_log_contrast(fit_int[[r]], Z = ztr, y = ytr,
                                                   additional_covariates =
                                                     xtr_int_normalization,
                                                   folds = cvfit[[r]]$folds)


    yhat_tr_int[[r]] <- predict_trac(list(fit_int[[r]]), new_Z = ztr,
                                     new_additional_covariates =
                                       xtr_int_normalization, output = output)[[1]]
    yhat_te_int[[r]] <- predict_trac(list(fit_int[[r]]), new_Z = zte,
                                     new_additional_covariates =
                                       xte_int_normalization, output = output)[[1]]
    trainerr_int[[r]] <- sqrt(colMeans((yhat_tr_int[[r]] - ytr)^2))
    # plot(yhat_tr[[j]][, cvfit[[j]]$cv$ibest], ytr)
    testerr_int[[r]] <- sqrt(colMeans((yhat_te_int[[r]] - yte)^2))
    #plot(yhat_te[[j]][, cvfit[[j]]$cv$ibest], yte)

    ## back calculation interaction coefficients beta
    if(normalize.int) {
    fit_int[[r]]$beta[(ncol(ztr) + 1):nrow(beta_int_est_refit), ] <- fit_int[[r]]$beta[(ncol(ztr) + 1):nrow(beta_int_est_refit), ] * avg_norm_clr_Xtr / xtr_int_norm
    fit_int[[r]]$beta[is.na(fit_int[[r]]$beta)] <- 0
    }

    nnz_int[[r]] <- colSums(abs(fit_int[[r]]$beta) > 1e-5)

    print(nnz_int[[r]][cvfit_int[[r]]$cv$i1se])
    if(nnz_int[[r]][cvfit_int[[r]]$cv$i1se] > 1){#!= 0){

      refit_lc_int <- refit_selected_model(fit_int[[r]], i_selected = cvfit_int[[r]]$cv$i1se,
                                                          Z = ztr, y = ytr, additional_covariates =
                                                            xtr_int_normalization)
      #colnames(refit_lc_int$beta) <- gsub("X", "", colnames(refit_lc_int$beta))
      colnames(refit_lc_int$beta) <- gsub('[.]', ':', colnames(refit_lc_int$beta))#, fixed = T)
      rownames(beta_int_est_refit) <- c(colnames(X)[1:p], colnames(xtr_int_normalization))
      beta_int_est_refit[colnames(refit_lc_int$beta), r] <- refit_lc_int$beta

      # ## back calculation of the interaction coefficients after scaling
      if(normalize.int){
      beta_int_est_refit[(ncol(ztr) + 1):nrow(beta_int_est_refit), r] <-
        beta_int_est_refit[(ncol(ztr) + 1):nrow(beta_int_est_refit), r] * avg_norm_clr_Xtr / xtr_int_norm ## tr or te??? should be tr
      }
      ## no refitting
      beta_int_est_norefit[, r] <- fit_int[[r]]$beta[, cvfit_int[[r]]$cv$i1se]
    }
    beta_int_est_refit[is.na(beta_int_est_refit)] <- 0


  }

  list_res <- list("tr" = tr,
                   "fit" = fit,
                   "cvfit" = cvfit,
                   "yhat_tr" = yhat_tr,
                   "yhat_te" = yhat_te,
                   "trainerr" = trainerr,
                   "testerr" = testerr,
                   "nnz" = nnz,
                   "beta_est_refit" = beta_est_refit,

                   "fit_int" = fit_int,
                   "cvfit_int" = cvfit_int,
                   "yhat_tr_int" = yhat_tr_int,
                   "yhat_te_int" = yhat_te_int,
                   "trainerr_int" = trainerr_int,
                   "testerr_int" = testerr_int,
                   "nnz_int" = nnz_int,
                   "beta_int_est_refit" = beta_int_est_refit,
                   "beta_int_est_norefit" = beta_int_est_norefit)

  return(list_res)

}





