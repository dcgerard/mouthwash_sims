#' @param calibrate a logical for whether to use MAD calibrated t-stats.
cate_nc <- function(Y, X, num_sv, control_genes, calibrate = FALSE, quant = 0.95) {
    calibrate <- as.logical(calibrate)
    cate_nc <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
                              X.nuis = X[, -2, drop = FALSE],
                              r = num_sv, adj.method = "nc",
                              nc = as.logical(control_genes),
                              calibrate = calibrate)

    betahat   <- c(cate_nc$beta)
    sebetahat <- c(sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) /
                   sqrt(nrow(X)))
    df        <- Inf
    pvalues   <- c(cate_nc$beta.p.value)

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat

    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}

#' @param calibrate a logical for whether to use MAD calibrated t-stats.
cate_rr <- function(Y, X, num_sv, calibrate = FALSE) {
    calibrate <- as.logical(calibrate)
    cate_rr <- cate::cate.fit(Y = Y, X.primary = X[, 2, drop = FALSE],
                              X.nuis = X[, -2, drop = FALSE],
                              r = num_sv, adj.method = "rr", calibrate = calibrate)
    betahat     <- c(cate_rr$beta)
    sebetahat   <- c(sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X)))
    pvalues     <- c(cate_rr$beta.p.value)
    df          <- Inf
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

leapp <- function(Y, X, num_sv) {
    ## LEAPP Sparse-------------------------------------------------------
    ## Sometimes stalls when sparse = TRUE
    trash2 <- tryCatch({
        trash <- R.utils::evalWithTimeout({
            leapp_sparse <- leapp::leapp(dat = t(Y), pred.prim = X[, 2, drop = FALSE],
                                         pred.covar = X[, -2, drop = FALSE],
                                         num.fac = num_sv)
            leapp_sparse_out         <- list()
            leapp_sparse_out$betahat <- leapp_sparse$gamma
            leapp_sparse_out$pvalues  <- leapp_sparse$p
            TRUE
        },
        timeout = 120,
        onTimeout = "silent")
    },
    error = function(e){NULL})
    if(is.null(trash2)) {
        leapp_sparse_out <- list()
        leapp_sparse_out$betahat <- rep(NA, length = ncol(Y))
        leapp_sparse_out$pvalues <- rep(NA, length = ncol(Y))
    } else if (is.null(trash)) {
        leapp_sparse_out <- list()
        leapp_sparse_out$betahat <- rep(NA, length = ncol(Y))
        leapp_sparse_out$pvalues <- rep(NA, length = ncol(Y))
    }
    return(leapp_sparse_out)
}

ols <- function(Y, X, quant = 0.95) {
    limma_out <- limma::lmFit(object = t(Y), design = X)
    betahat   <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df        <- limma_out$df.residual[1]
    tstats    <- betahat / sebetahat
    pvalues   <- 2 * pt(-abs(tstats), df = df)

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, upper = upper, lower = lower))
}

ruv2 <- function(Y, X, num_sv, control_genes, quant = 0.95) {
    ruv_ruv2 <- ruv::RUV2(Y = Y, X = X[, 2, drop = FALSE],
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = X[, -2, drop = FALSE])
    sebetahat <- sqrt(ruv_ruv2$sigma2 * ruv_ruv2$multiplier)
    betahat   <- c(ruv_ruv2$betahat)
    pvalues   <- c(ruv_ruv2$p)
    df        <- ruv_ruv2$df

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat

    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, upper = upper, lower = lower))
}

ruv4 <- function(Y, X, num_sv, control_genes, quant = 0.95) {
    ruv_ruv4 <- ruv::RUV4(Y = Y, X = X[, 2, drop = FALSE],
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = X[, -2, drop = FALSE])
    sebetahat <- sqrt(ruv_ruv4$sigma2 * ruv_ruv4$multiplier)
    betahat   <- c(ruv_ruv4$betahat)
    pvalues   <- c(ruv_ruv4$p)
    df        <- ruv_ruv4$df

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}

ruv4_rsvar_ebayes <- function(Y, X, num_sv, control_genes, quant = 0.95) {
    ruv_ruv4 <- ruv::RUV4(Y = Y, X = X[, 2, drop = FALSE],
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = X[, -2, drop = FALSE])
    ruv_ruv4inflate <- ruv::variance_adjust(fit = ruv_ruv4,
                                            evar = FALSE)
    sebetahat <- sqrt(ruv_ruv4inflate$varbetahat.rsvar.ebayes)
    betahat   <- c(ruv_ruv4inflate$betahat)
    pvalues   <- c(ruv_ruv4inflate$p.rsvar.ebayes)
    df        <- ruv_ruv4inflate$df.ebayes

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}


ruv2_rsvar_ebayes <- function(Y, X, num_sv, control_genes, quant = 0.95) {
    ruv_ruv2 <- ruv::RUV2(Y = Y, X = X[, 2, drop = FALSE],
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = X[, -2, drop = FALSE])
    ruv_ruv2inflate <- ruv::variance_adjust(fit = ruv_ruv2,
                                            evar = FALSE)
    sebetahat <- sqrt(ruv_ruv2inflate$varbetahat.rsvar.ebayes)
    betahat   <- c(ruv_ruv2inflate$betahat)
    pvalues   <- c(ruv_ruv2inflate$p.rsvar.ebayes)
    df        <- ruv_ruv2inflate$df.ebayes

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}

ruv4_rsvar <- function(Y, X, num_sv, control_genes) {
    ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, 2]),
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = as.matrix(X[, -2]))
    ruv_ruv4inflate <- ruv::variance_adjust(fit = ruv_ruv4,
                                            ebayes = FALSE,
                                            evar = FALSE)
    sebetahat <- sqrt(ruv_ruv4inflate$varbetahat.rsvar)
    betahat   <- c(ruv_ruv4inflate$betahat)
    pvalues   <- c(ruv_ruv4inflate$p.rsvar)
    df        <- rep(ruv_ruv4inflate$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ruv4_evar <- function(Y, X, num_sv, control_genes) {
    ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, 2]),
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = as.matrix(X[, -2]))
    ruv_ruv4inflate <- ruv::variance_adjust(fit = ruv_ruv4,
                                            ebayes = FALSE,
                                            rsvar = FALSE)
    sebetahat <- sqrt(ruv_ruv4inflate$varbetahat.evar)
    betahat   <- c(ruv_ruv4inflate$betahat)
    pvalues   <- c(ruv_ruv4inflate$p.evar)
    df        <- rep(ruv_ruv4inflate$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ruvinv <- function(Y, X, control_genes) {
    ruv_ruvinv <- ruv::RUVinv(Y = Y, X = as.matrix(X[, 2]),
                              ctl = as.logical(control_genes),
                              Z = as.matrix(X[, -2]), randomization = TRUE)
    sebetahat <- sqrt(ruv_ruvinv$sigma2 * ruv_ruvinv$multiplier)
    betahat   <- c(ruv_ruvinv$betahat)
    pvalues   <- c(ruv_ruvinv$p)
    df        <- rep(ruv_ruvinv$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

succotash <- function(Y, X, num_sv) {
    succ_out <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                      fa_method = "pca",
                                      optmethod = "em",
                                      two_step = TRUE,
                                      var_scale_init_type = "null_mle",
                                      z_init_type = "null_mle",
                                      likelihood = "normal")
    betahat <- succ_out$betahat
    lfdr    <- succ_out$lfdr
    pi0hat  <- succ_out$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

mouthwash <- function(Y, X, num_sv, likelihood = c("normal", "t"), scale_var = TRUE) {
    likelihood <- match.arg(likelihood)
    if (likelihood == "normal") {
        mout <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                 cov_of_interest = 2,
                                 include_intercept = FALSE,
                                 likelihood = "normal",
                                 mixing_dist = "normal",
                                 scale_var = scale_var)
    } else if (likelihood == "t") {
        mout <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                 cov_of_interest = 2,
                                 include_intercept = FALSE,
                                 likelihood = "t",
                                 mixing_dist = "sym_uniform",
                                 scale_var = scale_var)
    }
    betahat <- mout$result$PosteriorMean
    lfdr    <- mout$result$lfdr
    pi0hat  <- mout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

mouthwash_noscale <- function(Y, X, num_sv) {
    vout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, include_intercept = FALSE,
                             limmashrink = TRUE, likelihood = "t", mixing_dist = "uniform",
                             scale_var = FALSE)
    betahat <- vout$result$PosteriorMean
    lfdr    <- vout$result$lfdr
    pi0hat  <- vout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

mouthwash_scale <- function(Y, X, num_sv) {
    vout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, include_intercept = FALSE,
                             limmashrink = TRUE, likelihood = "t", mixing_dist = "uniform",
                             scale_var = TRUE)
    betahat <- vout$result$PosteriorMean
    lfdr    <- vout$result$lfdr
    pi0hat  <- vout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

mouthwash_normal_noscale <- function(Y, X, num_sv) {
    vout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, include_intercept = FALSE,
                             limmashrink = TRUE, likelihood = "normal", mixing_dist = "normal",
                             scale_var = FALSE)
    betahat <- vout$result$PosteriorMean
    lfdr    <- vout$result$lfdr
    pi0hat  <- vout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

mouthwash_normal_scale <- function(Y, X, num_sv) {
    vout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, include_intercept = FALSE,
                             limmashrink = TRUE, likelihood = "normal", mixing_dist = "normal",
                             scale_var = TRUE)
    betahat <- vout$result$PosteriorMean
    lfdr    <- vout$result$lfdr
    pi0hat  <- vout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}


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

vruv4 <- function(Y, X, num_sv, control_genes, adjust_bias = FALSE, quant = 0.95,
                  likelihood = c("normal", "t")) {
    likelihood <- match.arg(likelihood)
    vout <- vicar::vruv4(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes),
                         limmashrink = TRUE, cov_of_interest = 2, adjust_bias = adjust_bias,
                         likelihood = likelihood)
    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    pvalues   <- c(vout$pvalues)
    df        <- vout$degrees_freedom
    alpha     <- 1 - quant
    tval      <- qt(p = 1 - alpha / 2, df = df)
    lower     <- betahat - tval * sebetahat
    upper     <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}

vruv2 <- function(Y, X, num_sv, control_genes) {
    vout <- vicar::vruv2(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes),
                         limmashrink = TRUE, cov_of_interest = 2, likelihood = "t")
    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    pvalues   <- c(vout$pvalues)
    df        <- vout$degrees_freedom[1]
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}


vruvinv <- function(Y, X, control_genes) {
    vout <- vicar::vruvinv(Y = Y, X = X, ctl = as.logical(control_genes),
                           cov_of_interest = 2, likelihood = "t")
    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    pvalues   <- c(vout$pvalues)
    df        <- vout$degrees_freedom[1]
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}


ruv3 <- function(Y, X, num_sv, control_genes, multiplier = TRUE, quant = 0.95) {
    vout <- vicar::ruv3(Y = Y, X = X, ctl = as.logical(control_genes),
                        k = num_sv, cov_of_interest = 2)
    betahat   <- c(vout$betahat)
    if (multiplier) {
        sebetahat <- c(vout$sebetahat_adjusted)
        pvalues   <- c(vout$pvalues_adjusted)
    } else {
        sebetahat <- c(vout$sebetahat_unadjusted)
        pvalues   <- c(vout$pvalues_unadjusted)
    }
    df        <- vout$degrees_freedom

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues, lower = lower, upper = upper))
}


flashr_wrapper <- function(Y, max_rank) {
    if (!requireNamespace("flashr", quietly = TRUE)) {
        stop("Sorry, flashr needs to be installed to use flashr_wrapper.")
    }
    trash <- utils::capture.output(gout <- flashr::greedy(Y = Y, K = max_rank,
                                                   flash_para = list(partype = "var_col")))
    Yhat <- gout$l %*% t(gout$f)
    return(Yhat)
}


ruvimpute <- function(Y, X, control_genes, num_sv) {
    vout <- vicar::ruvimpute(Y = Y, X = X, k = num_sv,
                             ctl = as.logical(control_genes),
                             cov_of_interest = 2, do_variance = FALSE)
    betahat <- vout$betahat_long
    return(list(betahat = betahat))
}


ruvflash <- function(Y, X, control_genes) {
    max_rank <- min(sum(control_genes), nrow(X) - ncol(X)) - 1
    vout <- vicar::ruvimpute(Y = Y, X = X,
                             ctl = as.logical(control_genes),
                             impute_func = flashr_wrapper,
                             impute_args = list(max_rank = max_rank),
                             cov_of_interest = 2, do_variance = FALSE)
    betahat <- vout$betahat_long
    return(list(betahat = betahat))
}


ruvem <- function(Y, X, control_genes, num_sv, quant = 0.95) {
    ruvemout <- vicar::ruvem(Y = Y, X = X, ctl = control_genes, k = num_sv, cov_of_interest = 2)
    betahat <- ruvemout$betahat_long
    sebetahat <- ruvemout$sebetahat_limma
    df <- ruvemout$df_limma

    alpha <- 1 - quant
    tval  <- qt(p = 1 - alpha / 2, df = df)
    lower <- betahat - tval * sebetahat
    upper <- betahat + tval * sebetahat
    return(list(betahat = betahat,
                sebetahat = sebetahat,
                df = df,
                pvalues = ruvemout$pvalues_limma,
                lower = lower, upper = upper))
}

#' only allow for 95% credible intervals
ruvb_bfa_gs <- function(Y, X, control_genes, num_sv) {
    ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                           fa_func = vicar::bfa_gs, cov_of_interest = 2)
    return(list(betahat = ruvbout$means,
                pvalues = ruvbout$lfsr2,
                lower = ruvbout$lower,
                upper = ruvbout$upper))
}

#' only allow for 95% credible intervals
ruvb_bfa_gs_linked <- function(Y, X, control_genes, num_sv) {
    ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                           fa_func = vicar::bfa_gs_linked,
                           fa_args = list(use_code = "r", nsamp = 10000),
                           cov_of_interest = 2)
    return(list(betahat = ruvbout$means,
                pvalues = ruvbout$lfsr2,
                lower = ruvbout$lower,
                upper = ruvbout$upper))
}

#' only allow for 95% credible intervals
ruvb_bfa_gs_linked_norm <- function(Y, X, control_genes, num_sv) {
    ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                           fa_func = vicar::bfa_gs_linked, fa_args = list(use_code = "r"),
                           prior_fun = vicar::normal_prior, return_log = TRUE,
                           prior_args = list(prior_variance = 100),
                           cov_of_interest = 2)
    return(list(betahat = ruvbout$means,
                pvalues = ruvbout$lfsr2,
                lower = ruvbout$lower,
                upper = ruvbout$upper))
}

#' only allow for 95% credible intervals
ruvb_bfa_gs_linked_hier <- function(Y, X, control_genes, num_sv) {
    ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv,
                           fa_func = vicar::bfa_gs_linked, fa_args = list(use_code = "r"),
                           prior_fun = vicar::hier_fun, return_log = TRUE, cov_of_interest = 2)
    return(list(betahat = ruvbout$means,
                pvalues = ruvbout$lfsr2,
                lower = ruvbout$lower,
                upper = ruvbout$upper))
}


#' only allow for 95% credible intervals
ruvb_bfl <- function(Y, X, control_genes, num_sv) {
    ruvbout <- vicar::ruvb(Y = Y, X = X, ctl = control_genes, k = num_sv, fa_func = vicar::bfl,
                           cov_of_interest = 2)
    return(list(betahat = ruvbout$means,
                pvalues = ruvbout$lfsr2,
                lower = ruvbout$lower,
                upper = ruvbout$upper,
                svalues = ruvbout$svalues))
}

backwash <- function(Y, X, num_sv, scale_var = TRUE) {
    bout <- vicar::backwash(Y = Y, X = X, k = num_sv,
                            cov_of_interest = 2,
                            include_intercept = FALSE,
                            scale_var = scale_var)
    betahat <- bout$result$PosteriorMean
    lfdr    <- bout$result$lfdr
    pi0hat  <- bout$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}
