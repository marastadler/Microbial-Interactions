library(glmnet)
norm_vec <- function(x) sqrt(sum(x^2))

#' This code is taken from https://github.com/viettr/trac.git and extended to interactions


#' Two stage regression
#'
#' Perform a two stage fitting procedure similar to the one proposed by
#' Bates, S., & Tibshirani, R. (2019). Log‐ratio lasso:
#' Scalable, sparse estimation for log‐ratio models. Biometrics, 75(2), 613-624.
#' This function uses the selected components by trac or log-ratio and
#' fits a sparse model based on all possible log-ratios. (Implementation: Viet Tran)
#'
#' This function is an extension to quadratic terms in the second stage.
#'
#' 1. Fit trac or sparse-log contrast and extract the selected components.
#' 2. Fit a second-stage based on the selected components with log ratios and
#'.   log-ratios^2.
#'
#' @param Z n by p matrix containing log(X) (see \code{Z} from
#'     \code{\link{trac}})
#' @param y n vector (response) (see \code{y} from \code{\link{trac}})
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root).
#'   Only needed for trac based models. (see \code{A} from \code{\link{trac}}).
#'   If A = NULL, a sparse-log contrast model is assumed.
#' @param additional_covariates n by p' matrix containing additional covariates,
#'   (see \code{additional_covariates} from \code{\link{trac}}).
#'   Non-compositional components are currently not penalized by the lasso
#' @param betas pre-screened coefficients
#'   (see output \code{gamma} from \code{\link{trac}} or
#'   output \code{beta} from \code{\link{sparse_log_contrast}})
#' @param topk maximum number of pre-screened coefficients to consider.
#'   Default NULL.
#' @param nfolds number of folds
#' @param method string which estimation method to use should be "regr" or
#'   "classif"
#' @param criterion which criterion should be used to select the coefficients
#'   "1se" or "min"
#' @param alpha nudge the model to select on a higher or lower level of the
#'   tree. Only relevant for trac based models.
#' @param foldid an optional vector of values between 1 and nfold identifying
#'   what fold each observation is in. If supplied, nfold can be missing.
#' @return list with: log_ratios: betas for log ratios; index: dataframe with
#'   index of the pre selected coefficients and the log ratio name; A: taxonomic
#'   tree information; method: regression or classification; cv_glmnet:
#'   output of glmnet, useful for prediction; criterion: which criterion to
#'   be used to select lambda based on cv (cross validation); foldid;
#'   expanded_z: input matrix log-ratios and log-ratios^2
#' @export
#' So far the extension to interactions only works for sparse_log_contrast, but not trac!



second_stage_interaction <- function(Z, A = NULL, y,
                                     additional_covariates = NULL, betas,
                                     topk = NULL, nfolds = 5,
                                     method = c("regr", "classif"),
                                     criterion = c("1se", "min"),
                                     alpha = 0,
                                     foldid = NULL, standardize = TRUE) {
  # partial matching of the input
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  # if no selected components --> stop
  stopifnot(alpha >= 0)
  # get the number of compositional inputs
  if (!is.null(A)) {
    n_alphas_compositional <- ncol(A)
  } else {
    n_alphas_compositional <- ncol(Z)
  }

  # get the betas of the compositional features
  pre_selected <- betas[1:n_alphas_compositional]
  # if no names provided --> create new names
  if (is.null(names(pre_selected))) {
    names(pre_selected) <- paste0("gamma_", 1:length(pre_selected))
  }
  # get the non zero components
  pre_selected_index <- abs(pre_selected) >= 1e-07 #>= 1e-3

  if (length(pre_selected[pre_selected_index]) <= 2) {
    stop("Only two or less preselected variables are passed.
       Therefore using a second step doesn't make sense.")
  }
  # if topk is specified take only the top k components
  if (!is.null(topk)) {
    if (!is.numeric(topk)) stop("Topk has to be numeric")
    if (topk <= 2) stop("Minimum number of topk is 3")
    if (length(pre_selected[pre_selected_index]) > topk) {
      # create df in order to get topk values and replace the rest with 0
      # reorder the df into original order
      id <- 1:length(pre_selected)
      topk_df <- data.frame(id = id, pre_selected = abs(pre_selected),
                            pre_selected_index = pre_selected_index)
      topk_df <- topk_df[order(topk_df$pre_selected, decreasing = TRUE) ,]
      topk_df$pre_selected[(topk + 1):length(pre_selected)] <- 0
      topk_df <- topk_df[order(topk_df$id, decreasing = FALSE) ,]
      pre_selected <- topk_df$pre_selected
      names(pre_selected) <- rownames(topk_df)
      pre_selected_index <- abs(pre_selected) > 0
    }
  }


  # Change to geometric mean --> see formula 2, only necessary for trac
  if (!is.null(A)) {
    Z_A <- as.matrix(Z %*% A)
    subtree_n <- colSums(as.matrix(A))
    Z_A <- t((1 / subtree_n) * t(Z_A))
  } else {
    Z_A <- Z
  }
  # create all logratios
  #  Z_A <- scale(Z_A, center = TRUE, scale = FALSE)
  # build pairs
  index <- expand.grid(which(pre_selected_index), which(pre_selected_index))
  index <- index[index$Var2 > index$Var1, ]
  # choose(length(pre_selected[pre_selected_index]), 2)
  index <- as.data.frame(index)
  if (is.null(colnames(Z_A))) {
    colnames(Z_A) <- paste0("comp_", 1:length(pre_selected))
  }
  Z_A_names <- colnames(Z_A)
  n <- nrow(Z_A)
  colnames(index) <- c("index1", "index2")

  index$variable1 <- Z_A_names[index$index1]
  index$variable2 <- Z_A_names[index$index2]
  # calculate the penalty factor for trac based models to nudge the model
  # to select pairs on a specific level, for sparse-log contrast set everything
  # to 1
  if (!is.null(A)) {
    penalty_weight_comp <- subtree_n[index$index1] * subtree_n[index$index2]
    penalty_weight_comp <- penalty_weight_comp^alpha
    penalty_weight_comp <- 1 / penalty_weight_comp
  } else {
    penalty_weight_comp <- rep(1, nrow(index))
  }


  # seperate name by ///
  index$log_name <- paste(index$variable1,
                          index$variable2, sep = "///")
  n_log_contrast <- nrow(index)

  # actually build the logratios
  expanded_z <- matrix(0, nrow = n, ncol = 2 * n_log_contrast)
  colnames(expanded_z) <- c(index$log_name, paste0(index$log_name,"_2"))

  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }

  expanded_z[, (n_log_contrast + 1):(2 * n_log_contrast)] <-
    expanded_z[, 1:n_log_contrast]^2

  ## Normalize expanded_z such that all columns are of norm 1
  # expanded_z_norm <- apply(expanded_z, 2, norm_vec)
  # expanded_z <- t(t(expanded_z)/(expanded_z_norm))

  n_comp <- ncol(expanded_z)

  # add non compositional covariates and transform them
  # non-compositional components are currently not penalized by the lasso
  if (!is.null(additional_covariates)) {
    categorical_list <- get_categorical_variables(additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      additional_covariates[, categorical] <-
        transform_categorical_variables(additional_covariates, categorical)
    }
    n_x <- ncol(additional_covariates)
    penalty_factor <- c(penalty_weight_comp, rep(0, n_x))
    expanded_z <- cbind(expanded_z, additional_covariates)
    expanded_z <- as.matrix(expanded_z)
  } else {
    penalty_factor <- rep(penalty_weight_comp, 2) ## same penalty for log-ratio and log-ratio^2
  }

  # change the label for glmnet
  if (method == "classif") {
    if (!is.factor(y)) {
      y_new <- as.factor(y)
    } else {
      y_new <- y
    }
  } else {
    y_new <- y
  }
  # run glmnet with cross-validation
  if (method == "regr") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new,
                                   alpha = 1,
                                   nfolds = nfolds,
                                   foldid = foldid,
                                   keep = TRUE,
                                   standardize = standardize,
                                   family = "gaussian",
                                   penalty.factor =
                                     penalty_factor)
  } else if (method == "classif") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new,
                                   alpha = 1,
                                   nfolds = nfolds,
                                   foldid = foldid,
                                   keep = TRUE,
                                   standardize = standardize,
                                   family = "binomial",
                                   type.measure = "class",
                                   penalty.factor = penalty_factor
                                   )
  }
  # select features based on one standard error rule or lambda that
  # minimizes the error
  if (criterion == "1se") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.1se"))[, 1]
  }
  if (criterion == "min") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.min"))[, 1]
  }

  # return results in list
  list(log_ratios = log_ratios[2:length(log_ratios)],
       log_ratios_path =  cv_glmnet$glmnet.fit$beta,
       beta0 = log_ratios[1],
       index = index,
       A = A,
       method = method,
       cv_glmnet = cv_glmnet,
       criterion = criterion,
       foldid = cv_glmnet$foldid,
       expanded_z = expanded_z)
}




#' Make predictions based on a second_stage fit
#'
#' @param new_Z a new data matrix (see \code{Z} from \code{\link{second_stage}})
#' @param new_additional_covariates a new data matrix
#'   (see \code{additional_covariates} from \code{\link{second_stage}})
#' @param fit output of the function \code{\link{second_stage}}
#' @param output string  either "raw", "probability" or "class" only relevant
#'   classification tasks
#' @return a vector of \code{nrow(new_Z) + nrow(new_additional_covariates)}
#'   predictions.
#' @importFrom stats predict
#' @export


predict_second_stage_interaction <- function(new_Z, new_additional_covariates = NULL, fit,
                                 output = c("raw", "probability", "class")) {
  # partial matching for input
  output <- match.arg(output)
  # create geometric mean
  if (!is.null(fit$A)) {
    A <- as.matrix(fit$A)
    Z_A <- as.matrix(new_Z %*% A)
    subtree_n <- colSums(A)
    Z_A <- t((1 / subtree_n) * t(Z_A))
  } else {
    Z_A <- new_Z
  }
  n <- nrow(Z_A)
  # create logratios
  index <- fit$index
  n_log_contrast <- nrow(index)

  expanded_z <- matrix(0, nrow = n, ncol = 2 * n_log_contrast)
  colnames(expanded_z) <- c(index$log_name, paste0(index$log_name,"_2"))


  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }

  expanded_z[, (n_log_contrast + 1):(2 * n_log_contrast)] <-
    expanded_z[, 1:n_log_contrast]^2

  ## Normalize expanded_z such that all columns are of norm 1
  # expanded_z_norm <- apply(expanded_z, 2, norm_vec)
  # expanded_z <- t(t(expanded_z)/(expanded_z_norm))

  # add aditional covariates
  if (!is.null(new_additional_covariates)) {
    categorical_list <- get_categorical_variables(new_additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      new_additional_covariates[, categorical] <-
        transform_categorical_variables(new_additional_covariates, categorical)
    }
    expanded_z <- cbind(expanded_z, new_additional_covariates)
    expanded_z <- as.matrix(expanded_z)
  }
  # specify how to select lambda
  if (fit$criterion == "1se") s <- "lambda.1se"
  if (fit$criterion == "min") s <- "lambda.min"
  # define different output possibilites for classification
  if (fit$method == "classif") {
    if (output == "raw") type <- "link"
    if (output == "probability") type <- "response"
    if (output == "class") type <- "class"
  } else if (fit$method == "regr") {
    type <- "link"
  }
  # make the prediction
  yhat <- predict(object = fit$cv_glmnet,
                  newx = expanded_z, s = s, type = type)
  # return prediction
  yhat
}


#' Two stage regression
#'
#' This is the version of the second_stage() function in
#' https://github.com/viettr/trac.git with slide adaptions to make it more
#' flexible: foldid as input, keep = TRUE
#'
#' Perform a two stage fitting procedure similar to the one proposed by
#' Bates, S., & Tibshirani, R. (2019). Log‐ratio lasso:
#' Scalable, sparse estimation for log‐ratio models. Biometrics, 75(2), 613-624.
#' This function uses the selected components by trac or log-ratio and
#' fits a sparse model based on all possible log-ratios.
#'
#' 1. Fit trac or sparse-log contrast and extract the selected components.
#' 2. Fit a second-stage based on the selected components.
#'
#' @param Z n by p matrix containing log(X) (see \code{Z} from
#'     \code{\link{trac}})
#' @param y n vector (response) (see \code{y} from \code{\link{trac}})
#' @param A p by (t_size-1) binary matrix giving tree structure (t_size is the
#'   total number of nodes and the -1 is because we do not include the root).
#'   Only needed for trac based models. (see \code{A} from \code{\link{trac}}).
#'   If A = NULL, a sparse-log contrast model is assumed.
#' @param additional_covariates n by p' matrix containing additional covariates,
#'   (see \code{additional_covariates} from \code{\link{trac}}).
#'   Non-compositional components are currently not penalized by the lasso
#' @param betas pre-screened coefficients
#'   (see output \code{gamma} from \code{\link{trac}} or
#'   output \code{beta} from \code{\link{sparse_log_contrast}})
#' @param topk maximum number of pre-screened coefficients to consider.
#'   Default NULL.
#' @param nfolds number of folds
#' @param method string which estimation method to use should be "regr" or
#'   "classif"
#' @param criterion which criterion should be used to select the coefficients
#'   "1se" or "min"
#' @param alpha nudge the model to select on a higher or lower level of the
#'   tree. Only relevant for trac based models.
#' @param foldid an optional vector of values between 1 and nfold identifying
#'   what fold each observation is in. If supplied, nfold can be missing.
#'
#' @return list with: log_ratios: betas for log ratios; index: dataframe with
#'   index of the pre selected coefficients and the log ratio name; A: taxonomic
#'   tree information; method: regression or classification; cv_glmnet:
#'   output of glmnet, useful for prediction; criterion: which criterion to
#'   be used to select lambda based on cv (cross validation); foldid
#' @export

second_stage_adapt <- function(Z, A = NULL, y, additional_covariates = NULL, betas,
                         topk = NULL, nfolds = 5,
                         method = c("regr", "classif"),
                         criterion = c("1se", "min"),
                         alpha = 0, foldid = NULL, standardize = TRUE) {
  # partial matching of the input
  criterion <- match.arg(criterion)
  method <- match.arg(method)
  # if no selected components --> stop
  stopifnot(alpha >= 0)
  # get the number of compositional inputs
  if (!is.null(A)) {
    n_alphas_compositional <- ncol(A)
  } else {
    n_alphas_compositional <- ncol(Z)
  }

  # get the betas of the compositional features
  pre_selected <- betas[1:n_alphas_compositional]
  # if no names provided --> create new names
  if (is.null(names(pre_selected))) {
    names(pre_selected) <- paste0("gamma_", 1:length(pre_selected))
  }
  # get the non zero components
  pre_selected_index <- abs(pre_selected) >= 1e-3

  if (length(pre_selected[pre_selected_index]) <= 2) {
    stop("Only two or less preselected variables are passed.
       Therefore using a second step doesn't make sense.")
  }
  # if topk is specified take only the top k components
  if (!is.null(topk)) {
    if (!is.numeric(topk)) stop("Topk has to be numeric")
    if (topk <= 2) stop("Minimum number of topk is 3")
    if (length(pre_selected[pre_selected_index]) > topk) {
      # create df in order to get topk values and replace the rest with 0
      # reorder the df into original order
      id <- 1:length(pre_selected)
      topk_df <- data.frame(id = id, pre_selected = abs(pre_selected),
                            pre_selected_index = pre_selected_index)
      topk_df <- topk_df[order(topk_df$pre_selected, decreasing = TRUE) ,]
      topk_df$pre_selected[(topk + 1):length(pre_selected)] <- 0
      topk_df <- topk_df[order(topk_df$id, decreasing = FALSE) ,]
      pre_selected <- topk_df$pre_selected
      names(pre_selected) <- rownames(topk_df)
      pre_selected_index <- abs(pre_selected) > 0
    }
  }


  # Change to geometric mean --> see formula 2, only necessary for trac
  if (!is.null(A)) {
    Z_A <- as.matrix(Z %*% A)
    subtree_n <- colSums(as.matrix(A))
    Z_A <- t((1 / subtree_n) * t(Z_A))
  } else {
    Z_A <- Z
  }
  # create all logratios
  #  Z_A <- scale(Z_A, center = TRUE, scale = FALSE)
  # build pairs
  index <- expand.grid(which(pre_selected_index), which(pre_selected_index))
  index <- index[index$Var2 > index$Var1, ]
  # choose(length(pre_selected[pre_selected_index]), 2)
  index <- as.data.frame(index)
  if (is.null(colnames(Z_A))) {
    colnames(Z_A) <- paste0("comp_", 1:length(pre_selected))
  }
  Z_A_names <- colnames(Z_A)
  n <- nrow(Z_A)
  colnames(index) <- c("index1", "index2")

  index$variable1 <- Z_A_names[index$index1]
  index$variable2 <- Z_A_names[index$index2]
  # calculate the penalty factor for trac based models to nudge the model
  # to select pairs on a specific level, for sparse-log contrast set everything
  # to 1
  if (!is.null(A)) {
    penalty_weight_comp <- subtree_n[index$index1] * subtree_n[index$index2]
    penalty_weight_comp <- penalty_weight_comp^alpha
    penalty_weight_comp <- 1 / penalty_weight_comp
  } else {
    penalty_weight_comp <- rep(1, nrow(index))
  }


  # seperate name by ///
  index$log_name <- paste(index$variable1,
                          index$variable2, sep = "///")
  n_log_contrast <- nrow(index)

  # actually build the logratios
  expanded_z <- matrix(0, nrow = n, ncol = n_log_contrast)
  colnames(expanded_z) <- index$log_name

  for (i in 1:n_log_contrast) {
    expanded_z[, i] <- Z_A[, index[i, "index1"]] - Z_A[, index[i, "index2"]]
  }
  n_comp <- ncol(expanded_z)

  # add non compositional covariates and transform them
  # non-compositional components are currently not penalized by the lasso
  if (!is.null(additional_covariates)) {
    categorical_list <- get_categorical_variables(additional_covariates)
    categorical <- categorical_list[["categorical"]]
    n_categorical <- categorical_list[["n_categorical"]]
    if (n_categorical > 0) {
      additional_covariates[, categorical] <-
        transform_categorical_variables(additional_covariates, categorical)
    }
    n_x <- ncol(additional_covariates)
    penalty_factor <- c(penalty_weight_comp, rep(0, n_x))
    expanded_z <- cbind(expanded_z, additional_covariates)
    expanded_z <- as.matrix(expanded_z)
  } else {
    penalty_factor <- penalty_weight_comp
  }

  # change the label for glmnet
  if (method == "classif") {
    if (!is.factor(y)) {
      y_new <- as.factor(y)
    } else {
      y_new <- y
    }
  } else {
    y_new <- y
  }



  # run glmnet with cross-validation


  if (method == "regr") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new, alpha = 1,
                                   nfolds = nfolds, standardize = standardize,
                                   family = "gaussian", penalty.factor =
                                     penalty_factor, foldid = foldid,
                                   keep = TRUE)
  } else if (method == "classif") {
    cv_glmnet <- glmnet::cv.glmnet(x = expanded_z, y = y_new, alpha = 1,
                                   nfolds = nfolds, standardize = standardize,
                                   family = "binomial", type.measure = "class",
                                   penalty.factor = penalty_factor,
                                   foldid = foldid,
                                   keep = TRUE)
  }
  # select features based on one standard error rule or lambda that
  # minimizes the error
  if (criterion == "1se") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.1se"))[, 1]
  }
  if (criterion == "min") {
    log_ratios <- as.matrix(
      glmnet::coef.glmnet(cv_glmnet, s = "lambda.min"))[, 1]
  }
  # return results in list
  list(log_ratios = log_ratios[2:length(log_ratios)],
       beta0 = log_ratios[1],
       index = index,
       A = A,
       method = method,
       cv_glmnet = cv_glmnet,
       criterion = criterion,
       foldid = cv_glmnet$foldid)
}




## transform beta coefficients log-ratio lasso such that they don't correspond to log-ratios but to individual features
library(dplyr)
#' @reg_second_stage output from second_stage_interaction()
#'
transform_beta_all_splits <- function(reg_second_stage = reg_second_stage_int_nf1, nsplit = 10, p, feature_names){
  Beta_lin.trafo <- matrix(ncol = p, nrow = nsplit)
  colnames(Beta_lin.trafo) <- feature_names

  for(r in seq(nsplit)){
    b <- reg_second_stage[[r]]$log_ratios
    ## extract main effect log-ratios
    b_lin_ind <- which(substr(names(b), nchar(names(b)) - 1, nchar(names(b))) != "_2")
    b_lin <- b[b_lin_ind]

    ## extract feature names
    names(b_lin) <- gsub("///", "-", names(b_lin))
    split_log_ratios <- strsplit(names(b_lin), '-')
    #ind_logratio <- regmatches(names(b_lin), gregexpr("[[:digit:]]+", names(b_lin)))
    ## from list to data frame with two columns
    split_log_ratios_df <- do.call(rbind, split_log_ratios )
    numerator_logratio <- split_log_ratios_df[, 1]
    denominator_logratio <- split_log_ratios_df[, 2]

    ## match coefficients to features (positive sign numerator and negative sign denominator)
    b_lin_numerator_logratio <- data.frame("feature" = numerator_logratio , "beta" = b_lin)
    b_lin_denominator_logratio <- data.frame("feature" = denominator_logratio , "beta" = -b_lin)

    ## merge numerator and denominator part for summary
    b_lin_transformed <- dplyr::bind_rows(b_lin_numerator_logratio, b_lin_denominator_logratio)

    ## sum over betas (group by feature), receive 30 entries
    beta_lin_trafo = b_lin_transformed %>%
      group_by(feature) %>%
      summarise(beta = sum(beta))



    Beta_lin.trafo[r , beta_lin_trafo$feature] <- beta_lin_trafo$beta
  }
  Beta_lin.trafo[is.na(Beta_lin.trafo)] <- 0
  Beta_main <- Beta_lin.trafo
  return(Beta_main)

}


## function to extract interaction coefficients from model object for all splits
get_interactions_all_splits <- function(reg_second_stage = reg_second_stage_int_nf1, nsplit = 10, p = 50, feature_names){
  ii = 0
  int_names <- c()
  for(i in 1:(p - 1)) {
    for(j in (i + 1): p){
      ii <- ii + 1
      int_names[ii] = paste0(feature_names[i], "///", feature_names[j], "_2")
    }
  }
  ## empty matrix with all possible combis of log-ratios
  Bint <- matrix(nrow = nsplit, ncol = length(int_names))
  colnames(Bint) <- int_names


  for(r in seq(nsplit)){
    b <- reg_second_stage[[r]]$log_ratios
    b_int_ind <- which(substr(names(b), nchar(names(b)) - 1, nchar(names(b))) == "_2")
    Bint[r, names(b)[b_int_ind]] <- b[b_int_ind]
  }
  return(Bint)
}

#'
transform_beta_path <- function(beta_path, p, feature_names){

  Beta_lin.trafo <- matrix(nrow = nrow(beta_path), ncol = p)
  colnames(Beta_lin.trafo) <- feature_names

  for(l in seq(nrow(beta_path))){
    b <- beta_path[l, ]
    ## extract main effect log-ratios
    b_lin_ind <- which(substr(names(b), nchar(names(b)) - 1, nchar(names(b))) != "_2")
    b_lin <- b[b_lin_ind]

    ## extract feature names
    names(b_lin) <- gsub("///", "-", names(b_lin))
    split_log_ratios <- strsplit(names(b_lin), '-')
    #ind_logratio <- regmatches(names(b_lin), gregexpr("[[:digit:]]+", names(b_lin)))
    ## from list to data frame with two columns
    split_log_ratios_df <- do.call(rbind, split_log_ratios )
    numerator_logratio <- split_log_ratios_df[, 1]
    denominator_logratio <- split_log_ratios_df[, 2]

    ## match coefficients to features (positive sign numerator and negative sign denominator)
    b_lin_numerator_logratio <- data.frame("feature" = numerator_logratio , "beta" = b_lin)
    b_lin_denominator_logratio <- data.frame("feature" = denominator_logratio , "beta" = -b_lin)

    ## merge numerator and denominator part for summary
    b_lin_transformed <- dplyr::bind_rows(b_lin_numerator_logratio, b_lin_denominator_logratio)

    ## sum over betas (group by feature), receive 30 entries
    beta_lin_trafo = b_lin_transformed %>%
      group_by(feature) %>%
      summarise(beta = sum(beta))



    Beta_lin.trafo[l , beta_lin_trafo$feature] <- beta_lin_trafo$beta
  }
  Beta_lin.trafo[is.na(Beta_lin.trafo)] <- 0
  Beta_main <- Beta_lin.trafo
  return(Beta_main)

}

get_interactions_path <- function(beta_path, p, feature_names){
  ii = 0
  int_names <- c()
  for(i in 1:(p - 1)) {
    for(j in (i + 1): p){
      ii <- ii + 1
      int_names[ii] = paste0(feature_names[i], "///", feature_names[j], "_2")
    }
  }
  ## empty matrix with all possible combis of log-ratios
  Bint <- matrix(nrow = nrow(beta_path), ncol = length(int_names))
  colnames(Bint) <- int_names


  for(l in seq(nrow(beta_path))){
    b <- beta_path[l, ]
    b_int_ind <- which(substr(names(b), nchar(names(b)) - 1, nchar(names(b))) == "_2")
    Bint[l, names(b)[b_int_ind]] <- b[b_int_ind]
  }
  return(Bint)
}
