#' r: split
plot_misclass_oracle <- function(r, mod, yhat_te_ = yhat_te,
                                 yhat_te_int_ = yhat_te_int, y_ = y, tr_ = tr,
                                 nnz_ = nnz, nnz_int_ = nnz_int, ylim_ = c(0,1)){
  yhat_te_r <- yhat_te_[[r]]

  yhat_te_int_r <- yhat_te_int_[[r]]

  misclassif_te_slc_r <- colSums(y_[-tr_[[r]]] != yhat_te_r)/length(y_[-tr_[[r]]])
  misclassif_te_slc_int_r <- colSums(y_[-tr_[[r]]] != yhat_te_int_r)/length(y_[-tr_[[r]]])

  plot(nnz_int_[[r]], misclassif_te_slc_int_r, type = "l",
       col = "blue", ylab = "OOS misclassification rate", xlab = "Number of nonzero features",
       main = paste("Split", r), cex.main = .5, ylim = ylim_)
  points(nnz_[[r]], misclassif_te_slc_r, type = "l", col = "red")

  ind_lam_b <- mod$cvfit[[r]]$cv$ibest
  ind_lam_b2 <- mod$cvfit_int[[r]]$cv$ibest
  ind_lam_1se <- mod$cvfit[[r]]$cv$i1se
  ind_lam_1se2 <- mod$cvfit_int[[r]]$cv$i1se

  abline(v = mod$nnz[[r]][ind_lam_b], col = "red")
  abline(v = mod$nnz[[r]][ind_lam_1se], col = "red", lty = 2)
  abline(v = mod$nnz_int[[r]][ind_lam_b2], col = "blue")
  abline(v = mod$nnz_int[[r]][ind_lam_1se2], col = "blue", lty = 2)
 # legend("bottomright", legend = c("slc", "slc + int"), col = c("red", "blue"),
 #        lwd = 1)
}


#' @ii either "i1se", "ibest", or NULL, if solution is directly provided instead of entire path
#'
boxplot_misclassif <- function(y_ = y, tr_ = tr,
                               yhat_tr_ = yhat_tr, yhat_te_ = yhat_te,
                               cvfit_ = cvfit,
                               yhat_tr_int_ = yhat_tr_int,
                               yhat_te_int_ = yhat_te_int,
                               cvfit_int_ = cvfit_int,
                               n_ = n, ntot_ = ntot,
                               ii = "i1se", nsplit = 10
){
  misclassif_tr_slc <- c()
  misclassif_te_slc <- c()

  misclassif_tr_slc_int <- c()
  misclassif_te_slc_int <- c()

  for(r in seq(nsplit)){

    if(!is.null(ii)) {
      yhat_tr_sel <- yhat_tr_[[r]][,  cvfit_[[r]]$cv[[ii]]]
      yhat_te_sel <- yhat_te_[[r]][,  cvfit_[[r]]$cv[[ii]]]

      yhat_tr_int_sel <- yhat_tr_int_[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
      yhat_te_int_sel <- yhat_te_int_[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
    }
    else{ ## selected y_hat can also be passed directly
      yhat_tr_sel <- yhat_tr_[[r]]
      yhat_te_sel <- yhat_te_[[r]]

      yhat_tr_int_sel <- yhat_tr_int_[[r]]
      yhat_te_int_sel <- yhat_te_int_[[r]]
    }

    misclassif_tr_slc[r] <- sum(y_[tr_[[r]]] != yhat_tr_sel)/length(y_[tr_[[r]]])
    misclassif_te_slc[r] <- sum(y_[-tr_[[r]]] != yhat_te_sel)/length(y_[-tr_[[r]]])

    misclassif_tr_slc_int[r] <- sum(y_[tr_[[r]]] != yhat_tr_int_sel)/length(y_[tr_[[r]]])
    misclassif_te_slc_int[r] <- sum(y_[-tr_[[r]]] != yhat_te_int_sel)/length(y_[-tr_[[r]]])
  }

  boxplot(misclassif_tr_slc, misclassif_te_slc, misclassif_tr_slc_int, misclassif_te_slc_int,
          names = c(paste0("slc train (n=", n_, ")"),
                    paste0("slc test (n=", ntot_ - n_, ")"),
                    paste0("slc + int train (n=", n_, ")"),
                    paste0("slc + int test (n=", ntot_ - n_, ")")),
          cex.axis = 0.65,
          col = "white", ylab = "Misclassification rate")
}


library(pROC)
# raw = T: takes the score instead off class
# source("utils.R")
plot_auc_oracle <- function(r, mod, yhat_te_ = yhat_te,
                            yhat_te_int_ = yhat_te_int,
                            y_ = y, tr_ = tr,
                            nnz_ = nnz, nnz_int_ = nnz_int, raw = F,
                            fit_ = fit, fit_int_ = fit_int, X_ = X){



  ytr <- y_[tr_[[r]]]
  yte <- y_[-tr_[[r]]]




  auc_te_slc_r <- c()
  auc_te_slc_int_r <- c()

  if(raw){
    ztr <- log(X_[tr_[[r]], ])
    zte <- log(X_[-tr_[[r]], ])

    ## calculate interactions based on X with pseudo counts
    xtr_int <- compute.interactions.aitchison(X_[tr_[[r]], ])
    xte_int <- compute.interactions.aitchison(X_[-tr_[[r]], ])

    ### normalize interaction effects
    ## 1. clr transform main effects
    xtr_clr <- clr.trafo(X_[tr_[[r]], ])
    xte_clr <- clr.trafo(X_[-tr_[[r]], ])

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


    yhat_te_ <- list()
    yhat_te_int_ <- list()
    yhat_te_[[r]] <- predict_trac(list(fit_[[r]]), new_Z = zte, output = "raw")[[1]]
    yhat_te_int_[[r]] <- predict_trac(list(fit_int_[[r]]), new_Z = zte,
                                     new_additional_covariates =
                                       xte_int_normalization, output = "raw")[[1]]

  }

  yhat_te_r <- yhat_te_[[r]]
  yhat_te_int_r <- yhat_te_int_[[r]]


    for(i in seq(ncol(yhat_te_r))){
      auc_te_slc_r[i] <- as.numeric(auc(y_[-tr_[[r]]], yhat_te_r[, i]))
    }


    for(i in seq(ncol(yhat_te_int_r))){
      auc_te_slc_int_r[i] <- auc(y_[-tr_[[r]]], yhat_te_int_r[, i])
    }




  plot(nnz_int_[[r]],
    auc_te_slc_int_r, type = "l",
       col = "blue", ylab = "AUC", xlab = "Number of nonzero features",
       main = paste("Split", r), cex.main = .5, ylim = c(0,1))
  points( nnz_[[r]],
         auc_te_slc_r, type = "l", col = "red")

  ind_lam_b <- mod$cvfit[[r]]$cv$ibest
  ind_lam_b2 <- mod$cvfit_int[[r]]$cv$ibest
  ind_lam_1se <- mod$cvfit[[r]]$cv$i1se
  ind_lam_1se2 <- mod$cvfit_int[[r]]$cv$i1se

  abline(v = mod$nnz[[r]][ind_lam_b], col = "red")
  abline(v = mod$nnz[[r]][ind_lam_1se], col = "red", lty = 2)
  abline(v = mod$nnz_int[[r]][ind_lam_b2], col = "blue")
  abline(v = mod$nnz_int[[r]][ind_lam_1se2], col = "blue", lty = 2)
 # legend("bottomright", legend = c("slc", "slc + int"), col = c("red", "blue"),
  #       lwd = 1)
}

#' @ii either "i1se", "ibest", or NULL, if solution is directly provided instead of entire path
#'
boxplot_auc <- function(y_ = y, yhat_train = yhat_tr, yhat_test = yhat_te,
                        cvfit_ = cvfit,
                        yhat_train_int = yhat_tr_int,
                        yhat_test_int = yhat_te_int,
                        cvfit_int_ = cvfit_int, tr_ = tr, ii = "ibest", return= F){
  auc_tr_slc <- c()
  auc_te_slc <- c()

  auc_tr_slc_int <- c()
  auc_te_slc_int <- c()

  for(r in seq(nsplit)){

    if(!is.null(ii)) {
    yhat_tr_sel <- yhat_train[[r]][,  cvfit_[[r]]$cv[[ii]]]
    yhat_te_sel <- yhat_test[[r]][,  cvfit_[[r]]$cv[[ii]]]

    yhat_tr_int_sel <- yhat_train_int[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
    yhat_te_int_sel <- yhat_test_int[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
    }

    else{
      yhat_tr_sel <- as.vector(yhat_train[[r]])
      yhat_te_sel <- as.vector(yhat_test[[r]])

      yhat_tr_int_sel <- as.vector(yhat_train_int[[r]])
      yhat_te_int_sel <- as.vector(yhat_test_int[[r]])
    }

    if(!all(is.na(yhat_train[[r]]))){
    auc_tr_slc[r] <- auc(y_[tr_[[r]]], yhat_tr_sel)}
    else{auc_tr_slc[r] <- NA}
    if(!all(is.na(yhat_test[[r]]))){
    auc_te_slc[r] <- auc(y_[-tr_[[r]]], yhat_te_sel)}
    else{auc_te_slc[r] <- NA}
    if(!all(is.na(yhat_train_int[[r]]))){
    auc_tr_slc_int[r] <- auc(y_[tr_[[r]]], yhat_tr_int_sel)}
    else{auc_tr_slc_int[r] <- NA}
    if(!all(is.na(yhat_test_int[[r]]))){
    auc_te_slc_int[r] <- auc(y_[-tr_[[r]]], yhat_te_int_sel)}
    else{auc_te_slc_int[r] <- NA}
  }
  if(return){print(list("auc_tr_slc" = summary(auc_tr_slc), "auc_te_slc" = summary(auc_te_slc),
                        "auc_tr_slc_int" = summary(auc_tr_slc_int),
                        "auc_te_slc_int" = summary(auc_te_slc_int)))}

  boxplot(auc_tr_slc, auc_te_slc, auc_tr_slc_int, auc_te_slc_int,
          names = c(paste0("slc train (n=", n, ")"),
                    paste0("slc test (n=", ntot - n, ")"),
                    paste0("slc + int train (n=", n, ")"),
                    paste0("slc + int test (n=", ntot - n, ")")),
          cex.axis = 0.65,
          col = "white", ylab = "AUC" )
}


compute_auc_all_splits <- function(y_ = y, yhat_train = yhat_tr, yhat_test = yhat_te,
                        cvfit_ = cvfit,
                        yhat_train_int = yhat_tr_int,
                        yhat_test_int = yhat_te_int,
                        cvfit_int_ = cvfit_int, tr_ = tr, ii = "ibest"){
  auc_tr_slc <- c()
  auc_te_slc <- c()

  auc_tr_slc_int <- c()
  auc_te_slc_int <- c()

  for(r in seq(nsplit)){

    if(!is.null(ii)) {
      yhat_tr_sel <- yhat_train[[r]][,  cvfit_[[r]]$cv[[ii]]]
      yhat_te_sel <- yhat_test[[r]][,  cvfit_[[r]]$cv[[ii]]]

      yhat_tr_int_sel <- yhat_train_int[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
      yhat_te_int_sel <- yhat_test_int[[r]][,  cvfit_int_[[r]]$cv[[ii]]]
    }

    else{
      yhat_tr_sel <- as.vector(yhat_train[[r]])
      yhat_te_sel <- as.vector(yhat_test[[r]])

      yhat_tr_int_sel <- as.vector(yhat_train_int[[r]])
      yhat_te_int_sel <- as.vector(yhat_test_int[[r]])
    }

    if(!all(is.na(yhat_train[[r]]))){
      auc_tr_slc[r] <- auc(y_[tr_[[r]]], yhat_tr_sel)}
    else{auc_tr_slc[r] <- NA}
    if(!all(is.na(yhat_test[[r]]))){
      auc_te_slc[r] <- auc(y_[-tr_[[r]]], yhat_te_sel)}
    else{auc_te_slc[r] <- NA}
    if(!all(is.na(yhat_train_int[[r]]))){
      auc_tr_slc_int[r] <- auc(y_[tr_[[r]]], yhat_tr_int_sel)}
    else{auc_tr_slc_int[r] <- NA}
    if(!all(is.na(yhat_test_int[[r]]))){
      auc_te_slc_int[r] <- auc(y_[-tr_[[r]]], yhat_te_int_sel)}
    else{auc_te_slc_int[r] <- NA}
  }
  return(list("auc_tr_slc" = auc_tr_slc, "auc_te_slc" = auc_te_slc,
                        "auc_tr_slc_int" = auc_tr_slc_int,
                        "auc_te_slc_int" = auc_te_slc_int))


}

