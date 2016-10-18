fit_ash <- function(betahat, sebetahat, df) {
    outputlevel <- 3
    ashout <- ashr::ash(betahat = betahat, sebetahat = sebetahat,
                        df = df, outputlevel = outputlevel)

    lfdr    <- c(ashout$result$lfdr)
    pi0hat  <- c(ashr::get_pi0(ashout))
    if (outputlevel > 2) {
        betahat <- c(ashout$result$PosteriorMean)
        return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
    }
    return(list(lfdr = lfdr, pi0hat = pi0hat))
}


fit_locfdr <- function(pvalues) {
    pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])
    pvalues[pvalues == 1] <- max(pvalues[pvalues != 1])
    zvalues <- qnorm(p = pvalues)
    locout <- locfdr::locfdr(zvalues)
    lfdr <- locout$fdr
    return(lfdr)
}

fit_qvalue <- function(pvalues) {
    qout    <- qvalue::qvalue(p = pvalues)
    lfdr    <- c(qout$lfdr)
    pi0hat  <- c(qout$pi0)
    return(list(lfdr = lfdr, pi0hat = pi0hat))
}
