# I always assume that the variable of interest is the second one.

## First, re-define cate::fa.pc because it doesn't work when r = 1.
library(cate)
cate_fa <- function (Y, r)
{
  assertthat::assert_that(is.matrix(Y))
  assertthat::are_equal(length(r), 1)
  assertthat::assert_that(r >= 0 & r < min(dim(Y)))
  if (r == 0) {
    Gamma <- NULL
    Z <- NULL
    Sigma <- apply(Y, 2, function(x) mean(x^2))
  }
  else {
    svd_Y <- svd(Y)
    Gamma <- svd_Y$v[, 1:r, drop = FALSE] %*% diag(svd_Y$d[1:r], r, r)/sqrt(nrow(Y))
    Z <- sqrt(nrow(Y)) * svd_Y$u[, 1:r, drop = FALSE]
    Sigma <- apply(Y - Z %*% t(Gamma), 2, function(x) sum(x^2))/(nrow(Y) - r)
  }
  return(list(Gamma = Gamma, Z = Z, Sigma = Sigma))
}
R.utils::reassignInPackage(name = "fa.pc", pkgName = "cate", value = cate_fa)


## Basic OLS -----------------------------------------------------------------
ols <- function(Y, X, quant = 0.95) {
  limma_out <- limma::lmFit(object = t(Y), design = X)
  betahat   <- limma_out$coefficients[, 2]
  sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
  df        <- limma_out$df.residual[1]
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

## Specific version of RUVB I use -------------------------------------------
ruvb_bfa_gs_linked <- function(Y, X, control_genes, num_sv) {
  ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                         fa_func = vicar::bfa_gs_linked,
                         fa_args = list(use_code = "r", nsamp = 10000), ## nsamp should be 10000
                         cov_of_interest = 2)

  return(list(betahat = ruvbout$means,
              pvalues = ruvbout$lfsr2,
              lower = ruvbout$lower,
              upper = ruvbout$upper))
}


## Return betahat, sebetahat, and df
ruvb_bfa_gs_linked_se <- function(Y, X, control_genes, num_sv) {
  ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                         fa_func = vicar::bfa_gs_linked,
                         fa_args = list(use_code = "r", nsamp = 10000),
                         cov_of_interest = 2)

  return(list(betahat = ruvbout$means,
              sebetahat = ruvbout$sd,
              df = nrow(Y) - ncol(X) - num_sv,
              lower = ruvbout$lower,
              upper = ruvbout$upper))
}




## Methods to look at ---------------------------------------------------
cate_simp_nc_correction <- function(Y, X, num_sv, control_genes, calibrate = FALSE) {
  cate_nc <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
                            X.nuis = X[, -2, drop = FALSE],
                            r = num_sv, adj.method = "nc",
                            fa.method = "pc",
                            nc = as.logical(control_genes),
                            calibrate = calibrate,
                            nc.var.correction = TRUE)

  betahat   <- c(cate_nc$beta)
  sebetahat <- c(sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) /
                   sqrt(nrow(X)))
  df        <- nrow(Y) - ncol(X) - num_sv

  ## sanity check
  ## plot(c(cate_nc$beta.p.value), stats::pnorm(-abs(betahat / sebetahat)) * 2)

  return_list <- list(betahat = betahat, sebetahat = sebetahat, df = df)

  if (calibrate) {
    lambda <- stats::mad(x = betahat / sebetahat, center = 0)
    return_list$sebetahat <- sebetahat * lambda
    return_list$pvalues <- c(cate_nc$beta.p.value)
  }

  return(return_list)
}

ruv4_simp <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::vruv4(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                       likelihood = "normal", limmashrink = FALSE, gls = FALSE,
                       include_intercept = FALSE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_ols
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

cate_simp <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::vruv4(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                       likelihood = "normal", limmashrink = FALSE, include_intercept = FALSE,
                       gls = TRUE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_ols
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}



ruv2_simp <- function(Y, X, num_sv, control_genes) {
  vout <- ruv::RUV2(Y = Y, X = X[, 2, drop = FALSE],
                    ctl = control_genes, k = num_sv, Z = X[, -2, drop = FALSE])
  betahat   <- vout$betahat
  sebetahat <- sqrt(vout$sigma2 * vout$multiplier)
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

## Limma before gls for RUV4 ----------------------------------------
cate_limma <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::vruv4(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                       likelihood = "normal", limmashrink = TRUE, gls = TRUE,
                       include_intercept = FALSE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_ols
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

## RUV3 needs special attention since return 0's for betahats in ctl--

ruv3_simp <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = FALSE, include_intercept = FALSE, gls = TRUE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_unadjusted
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

ruv3_ctl_adjust <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = FALSE, include_intercept = FALSE, gls = TRUE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_adjusted
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

ruv3_limma_pre <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = TRUE, include_intercept = FALSE, gls = TRUE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_unadjusted
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

ruv3_limma_pre_adjust <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = TRUE, include_intercept = FALSE, gls = TRUE)
  betahat   <- vout$betahat
  sebetahat <- vout$sebetahat_adjusted
  df        <- nrow(Y) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

## Limma after gls for RUV3 --- since sebetahat is NA for control genes
ruv3_limma_post <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = FALSE, include_intercept = FALSE, gls = TRUE)

  lmout <- limma::squeezeVar(vout$simga2_unadjusted, df = nrow(Y) - ncol(X) - num_sv)
  betahat   <- vout$betahat
  sebetahat <- sqrt(vout$mult_matrix * lmout$var.post)
  sebetahat[control_genes] <- NA

  # Check
  # sqrt(vout$mult_matrix * vout$simga2_unadjusted)
  # vout$sebetahat_unadjusted
  df        <- lmout$df.prior + nrow(X) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

ruv3_limma_post_adjust <- function(Y, X, num_sv, control_genes) {
  vout <- vicar::ruv3(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2,
                      limmashrink = FALSE, include_intercept = FALSE, gls = TRUE)

  lmout <- limma::squeezeVar(vout$simga2_unadjusted, df = nrow(Y) - ncol(X) - num_sv)
  betahat   <- vout$betahat
  multiplier <- mean(t(vout$resid_mat) ^ 2 / lmout$var.post[control_genes])
  sebetahat <- sqrt(vout$mult_matrix * lmout$var.post * multiplier)
  sebetahat[control_genes] <- NA
  df        <- lmout$df.prior + nrow(X) - ncol(X) - num_sv
  return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}



## Limma shrinking variances (after gls for ruv3 and ruv4)----------------------------------
## Apply or not apply to all "simp" functions

#' @param obj A list whose first element is betahat, whose second element is sebetahat and
#'      whose third element is df.
limma_adjust <- function(obj) {
  betahat   <- obj[[1]]
  sebetahat <- obj[[2]]
  df        <- obj[[3]]
  lmout     <- limma::squeezeVar(var = sebetahat ^ 2, df = df)
  return(list(betahat = betahat, sebetahat = sqrt(lmout$var.post), df = df + lmout$df.prior))
}


## Adjustment of variances ------------------------------------
#' @param obj A list whose first element is betahat, whose second element is sebetahat
#'      and whose third element is df
#' @param control_genes The control genes
ctl_adjust <- function(obj, control_genes) {
  stopifnot(is.logical(control_genes))
  betahat   <- obj[[1]]
  sebetahat <- obj[[2]]
  df        <- obj[[3]]
  stopifnot(length(control_genes) == length(betahat))
  stopifnot(length(sebetahat) == length(betahat))
  stopifnot(all(sebetahat >= 0, na.rm = TRUE))
  mult_val      <- mean(betahat[control_genes] ^ 2 / sebetahat[control_genes] ^ 2)
  sebetahat_adjusted <- sqrt(mult_val) * sebetahat
  return(list(betahat = betahat, sebetahat = sebetahat_adjusted, df = df))
}

#' @param obj A list whose first element is betahat and whose second element is sebetahat,
#'     and whose third element is df
mad_adjust <- function(obj) {
  betahat   <- obj[[1]]
  sebetahat <- obj[[2]]
  df        <- obj[[3]]
  stopifnot(length(sebetahat) == length(betahat))
  stopifnot(all(sebetahat >= 0, na.rm = TRUE))
  mult_val <- stats::mad(betahat / sebetahat, center = 0, na.rm = TRUE)
  sebetahat_adjusted <- mult_val * sebetahat
  return(list(betahat = betahat, sebetahat = sebetahat_adjusted, df = df))
}

## Calculate CI and p-values
#' @param obj A list whose first element is betahat, whose second element is sebetahat,
#'     and whose third element is df
#' @param alpha We calculate 1 - alpha confidence intervals.
calc_ci_p <- function(obj, alpha = 0.05) {
  betahat   <- obj[[1]]
  sebetahat <- obj[[2]]
  df        <- obj[[3]]

  stopifnot(length(sebetahat) == length(betahat))
  stopifnot(all(sebetahat >= 0, na.rm = TRUE))

  tstats    <- betahat / sebetahat
  pvalues   <- 2 * stats::pt(-abs(tstats), df = df)
  tval  <- stats::qt(p = 1 - alpha / 2, df = df)
  lower <- betahat - tval * sebetahat
  upper <- betahat + tval * sebetahat
  return(list(betahat = betahat, pvalues = pvalues, lower = lower, upper = upper))
}

