get_auc <- function(lfdr, which_null) {
    if (length(unique(which_null)) == 1) {
        auc <- NA
    } else {
        auc <- pROC::roc(predictor = lfdr, response = which_null)$auc
    }
    return(auc)
}

get_fdr <- function(pvalues, which_null, fdr_level = 0.1) {
    assertthat::assert_that(is.logical(which_null))
    padjusted <- stats::p.adjust(p = pvalues, method = "BH")
    if (sum(padjusted < fdr_level) == 0) {
        fdp <- 0
    } else {
        fdp <- mean(which_null[padjusted < fdr_level])
    }
    return(fdp)
}

get_mse <- function(betahat, beta_true) {
    mse <- mean((betahat - beta_true) ^ 2)
    return(mse)
}


get_ks <- function(pvalues, which_null) {
    nullp <- pvalues[which_null]
    ks_pvalue <- stats::ks.test(nullp, qunif)$p.value
    return(ks_pvalue)
}
