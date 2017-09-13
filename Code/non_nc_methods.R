## non-control-gene based methods

## regularized CATE -------------------------------------------------------------
cate_rr <- function(Y, X, num_sv, calibrate = FALSE) {
  calibrate <- as.logical(calibrate)
  cate_rr <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
                            X.nuis = X[, -2, drop = FALSE],
                            r = num_sv, adj.method = "rr",
                            calibrate = calibrate, fa.method = "pc")
  betahat     <- c(cate_rr$beta)
  sebetahat   <- c(sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X)))
  pvalues     <- c(cate_rr$beta.p.value)
  df          <- Inf

  if (calibrate) {
    lambda <- stats::mad(x = betahat / sebetahat, center = 0)
    sebetahat <- sebetahat * lambda
  }

  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}

## SVA --------------------------------------------------------------------------
sva <- function(Y, X, num_sv) {
  trash     <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
  X.sv      <- cbind(X, sva_out$sv)
  limma_out <- limma::lmFit(object = t(Y), design = X.sv)
  betahat   <- limma_out$coefficients[, 2]
  sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
  df        <- limma_out$df.residual[1]
  tstats    <- betahat / sebetahat
  pvalues   <- 2 * pt(-abs(tstats), df = df)
  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}


sva_voom <- function(Y, X, num_sv) {
  trash      <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
  X.sv       <- cbind(X, sva_out$sv)
  voom_out   <- limma::voom(counts = t(Y), design = X.sv)
  limma_out  <- limma::lmFit(object = voom_out)
  ebayes_out <- limma::ebayes(fit = limma_out)
  betahat    <- limma_out$coefficients[, 2]
  sebetahat  <- sqrt(ebayes_out$s2.post) * limma_out$stdev.unscaled[, 2]
  df         <- ebayes_out$df.total[1]
  pvalues    <- ebayes_out$p.value[, 2]
  return(list(betahat = betahat, sebetahat = sebetahat, df = df,
              pvalues = pvalues))
}


mouthwash <- function(Y, X, num_sv, likelihood = c("normal", "t"), alpha = 0, scale_var = FALSE,
                      var_inflate_pen = 0, mixing_dist = NULL) {
  likelihood <- match.arg(likelihood)
  if (is.null(mixing_dist) & likelihood == "t") {
    mixing_dist <- "sym_uniform"
  } else if (is.null(mixing_dist)) {
    mixing_dist <- "normal"
  }
  mout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, likelihood = likelihood,
                           scale_var = scale_var, sprop = alpha, mixing_dist = mixing_dist,
                           var_inflate_pen = var_inflate_pen)
  return_list <- list()
  return_list$betahat <- mout$result$PosteriorMean
  return_list$lfdr    <- mout$result$lfdr
  return_list$pi0hat  <- mout$pi0
  return(return_list)
}


backwash <- function(Y, X, num_sv, alpha = 0, scale_var = FALSE, var_inflate_pen = 0) {
  mout <- vicar::backwash(Y = Y, X = X, k = num_sv, scale_var = scale_var, sprop = alpha,
                          var_inflate_pen = var_inflate_pen)
  return_list <- list()
  return_list$betahat <- mout$result$PosteriorMean
  return_list$lfdr    <- mout$result$lfdr
  return_list$pi0hat  <- mout$pi0
  return(return_list)
}
