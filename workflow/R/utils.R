compute.interactions.aitchison <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x2 <- matrix(nrow = n, ncol = p*(p + 1)/2 - p)
  ind <- 0
  col.names <- c()
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      ind <- ind + 1
      x2[, ind] <- (log(x[, i]) - log(x[, j]))^2
      col.names[ind] <- paste0(colnames(x)[i], ":", colnames(x)[j])
    }
  }
  colnames(x2) <- col.names
  return(x2)
}

norm_vec <- function(x) sqrt(sum(x^2))

clr.trafo <- function(X) (log(X) - rowMeans(log(X)))



plot_cv_slc <- function(cvfit_slc, superimpose = TRUE, main = NULL) {
  num_w <- length(cvfit_slc$cv)
  if (!superimpose) {
    r <- floor(sqrt(num_w))
    #graphics::par(mfrow = c(r, ceiling(length(iw) / r)))
    lapply(1:num_w, function(iw) plot_cv_trac_single_w(cvfit_slc$cv[[iw]]))
  } else {
    x <- cvfit_slc$cv
    graphics::par(mar = c(5, 5, 5, 1))
    xrang <- range(x$nonzeros)
    yrang <- range(c(x$m - x$se, x$m + x$se), na.rm = TRUE)
    graphics::plot(0, 0, xlab = "Number of nonzero coefficients",
                   ylab = "Cross-validation Error",
                   type = "n", xlim = xrang, ylim = yrang, main = main)
    ggb:::error_bars(x$nonzeros,
                     x$m - x$se,
                     x$m + x$se, width = 0.01,
                     col = "darkgrey")}
  graphics::lines(x$nonzeros, x$m, col = "darkgrey", pch = 19)
  graphics::points(x$nonzeros, x$m, pch = 19)
  graphics::abline(v = x$nonzeros[x$ibest], lty = 3)
  graphics::abline(v = x$nonzeros[x$i1se], lty = 3)


  invisible()
}


library(ggplot2)

plot_cv_slc_ggplot <- function(cvfit_slc, main = NULL) {
  x <- cvfit_slc$cv
  
  data <- data.frame(
    nonzeros = x$nonzeros,
    m = x$m,
    se = x$se
  )
  
  plt <- ggplot(data, aes(x = nonzeros, y = m)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.01) +
    labs(x = "Number of nonzero coefficients", y = "Cross-validation Error") +
    theme_minimal() +
    geom_vline(xintercept = x$nonzeros[x$ibest], linetype = "dashed") +
    geom_vline(xintercept = x$nonzeros[x$i1se], linetype = "dotted") +
    scale_linetype_manual(values = c("dashed", "dotted"), name = "Legend") +
    annotate("text", x = x$nonzeros[x$ibest], y = max(data$m), label = "Min.", vjust = -0.5, hjust = -0.5, color = "black") +
    annotate("text", x = x$nonzeros[x$i1se], y = max(data$m), label = "1se", vjust = -0.5, hjust = -0.5, color = "black")
  
  plt
}


plot_path <- function(
    slc_int = slc_int_nf1,
    cvslc_int = cvslc_int_nf1,
    p = 10,
    r = 1,
    colnames_main_nz = as.character(1:6),
    colnames_int_nz = c("1:2", "3:4", "8:9"),
    main = "Noise level", feature_names = names(beta_main_int_all),
    only_main = FALSE, only_int = FALSE
){
  
  
  
  B1 <- slc_int[[r]]$beta
  rownames(B1) = feature_names
  
  
  
  colplt <- rep("gray", p * (p + 1)/2)
  names(colplt) <- rownames(B1)
  colplt[colnames_main_nz] <- brewer.pal(length(colnames_main_nz), "Set2")
  colplt[colnames_int_nz] <- brewer.pal(length(colnames_int_nz), "Set1")
  
  if(only_main){
    
    par(oma=c(0, 0, 0, 5))
    matplot(x = log(slc_int[[r]]$fraclist), t(B1[1:p, ]), type = "l", col = colplt[1:p],
            ylab = expression(hat(beta)), xlab = expression(log(lambda)), lwd = 1.3,
            main = main, lty = 1)
    # if(ii == 1) {
    
    abline(v = log(cvslc_int[[r]]$cv$lambda_1se), lty = 2)
    abline(v = log(cvslc_int[[r]]$cv$lambda_best))
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend = c(colnames_main_nz, "spurious"), bg = "white",
           col = c(brewer.pal(length(colnames_main_nz), "Set2"), "gray"), lwd = 2)
    # }
    
    
    
  }

  if(only_int){

      par(oma=c(0, 0, 0, 5))
      matplot(x = log(slc_int[[r]]$fraclist), t(B1[(p + 1):nrow(B1), ]), type = "l", col = colplt[(p + 1):length(colplt)],
              ylab = expression(hat(beta)), xlab = expression(log(lambda)), lwd = 1.3,
              main = main, lty = 1)
      # if(ii == 1) {

      abline(v = log(cvslc_int[[r]]$cv$lambda_1se), lty = 2)
      abline(v = log(cvslc_int[[r]]$cv$lambda_best))
      legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend = c(colnames_int_nz, "spurious"), bg = "white",
             col = c(brewer.pal(length(colnames_int_nz), "Set2"), "gray"), lwd = 2)
      # }



  }
  if(!(only_main | only_int)){

    par(oma=c(0, 0, 0, 5))
    matplot(x = log(slc_int[[r]]$fraclist), t(B1), type = "l", col = colplt,
            ylab = expression(hat(beta)), xlab = expression(log(lambda)), lwd = 1.3,
            main = main, lty = 1)
    # if(ii == 1) {

    abline(v = log(cvslc_int[[r]]$cv$lambda_1se), lty = 2)
    abline(v = log(cvslc_int[[r]]$cv$lambda_best))
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, legend = c(colnames_main_nz, colnames_int_nz, "spurious"), bg = "white",
           col = c(brewer.pal(length(colnames_main_nz), "Set2"),
                   brewer.pal(length(colnames_int_nz), "Set1"), "gray"), lwd = 2)
    # }



  }
}



ggplot_path <- function(
    slc_int = slc_int_nf1,
    cvslc_int = cvslc_int_nf1,
    p = 10,
    r = 1,
    colnames_main_nz = as.character(1:6),
    colnames_int_nz = c("1:2", "3:4", "8:9"),
    main = "Noise level", feature_names = names(beta_main_int_all),
    only_main = FALSE, only_int = FALSE, legend_title = "True non-zero features"
){

  B1 <- slc_int[[r]]$beta
  rownames(B1) = feature_names
  B1 <- t(B1)
  # Create a data frame for ggplot
  beta_df <- data.frame(
    log_lambda = log(slc_int[[r]]$fraclist + 0.001),
    B1
  )



  colplt <- rep("gray", p * (p + 1)/2)
  names(colplt) <- colnames(B1)
  colplt[colnames_main_nz] <- brewer.pal(9, "Blues")[(9 - length(colnames_main_nz)) : 9]
  colplt[colnames_int_nz] <- brewer.pal(9, "Oranges")[(9 - length(colnames_int_nz)) : 9]

  # Reshape the data to long format for ggplot
  beta_df_long <- tidyr::gather(beta_df, key = "feature", value = "beta", -log_lambda)
  beta_df_long$feature <- gsub("\\.", ":", beta_df_long$feature)
  # Plot the beta path with ggplot2
  plt <- ggplot(beta_df_long, aes(x = log_lambda, y = beta, color = feature)) +
    geom_line() +
    scale_color_manual(values = colplt,
                       breaks = c("gray", colnames_main_nz, colnames_int_nz),
                       labels = c("spurious", colnames_main_nz, colnames_int_nz)) +
    labs(title = main,
         x = expression(log(lambda)),
         y = expression(hat(beta))) +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(fill=guide_legend(title="True nonzero features")) +
    geom_vline(xintercept = log(cvslc_int[[r]]$cv$lambda_best  + 0.001), linetype = "solid") +
    geom_vline(xintercept = log(cvslc_int[[r]]$cv$lambda_1se  + 0.001), linetype = "dashed") +
    guides(color=guide_legend(title=legend_title))

  return(plt)
}

