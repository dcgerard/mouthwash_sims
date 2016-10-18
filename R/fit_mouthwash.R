library(reshape2)
library(ggplot2)
library(stringr)

proc_wrapper <- function(predictor, response) {
    pROC::roc(predictor = predictor, response = response)$auc
}

topk <- function(predictor, response, num_look = 100) {
    sum(response[order(predictor)[1:num_look]])
}

tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")
num_sv_seq <- readRDS("../output/ruvbout/num_sv.Rds")

nseq <- rep(NA, length = length(tissue_vec))

source("../Code/adjustment_methods.R")

do_ash <- function(args) {
    if (!is.null(args$sebetahat) & !is.null(args$betahat)) {
        if (!is.null(args$df)) {
            if (args$df == Inf) {
                args$df <- NULL
            }
        }
        ashout <- ashr::ash.workhorse(betahat = c(args$betahat),
                                      sebetahat = c(args$sebetahat),
                                      df = args$df,
                                      outputlevel = 2)
        return_list <- list()
        return_list$betahat <- ashout$result$PosteriorMean
        return_list$lfdr    <- ashout$result$lfdr
        return_list$pi0hat  <- ashr::get_pi0(ashout)
        return(return_list)
    }
}

do_qvalue <- function(args) {
    if (!is.null(args$pvalues)) {
        trash <- tryCatch({
            qout <- qvalue::qvalue(p = args$pvalues)
            return_list        <- list()
            return_list$lfdr   <- qout$lfdr
            return_list$pi0hat <- qout$pi0
            TRUE
        }, error = function(e){NULL})
        if (is.null(trash)) {
            return_list        <- list()
            return_list$lfdr   <- rep(NA, length(args$pvalues))
            return_list$pi0hat <- NA
        }
        return(return_list)
    }
}

get_lfdr <- function(args) {
    args$lfdr
}

get_pvalues <- function(args) {
    args$pvalues
}

get_pi0hat <- function(args) {
    args$pi0hat
}

get_betahat <- function(args) {
    args$betahat
}


plist <- list()
betahat_list <- list()
for(tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    num_sv <- num_sv_seq[tissue_index]

    dat <- readRDS(paste0("../Output/cleaned_gtex_data/", current_tissue, ".Rds"))
    onsex <- dat$chrom == "X" | dat$chrom == "Y"
    onsex[is.na(onsex)] <- FALSE
    dat$ctl[onsex] <- FALSE
    nseq[tissue_index] <- ncol(dat$Y)

    cat(tissue_index, "\n")

    Y <- t(dat$Y)
    X <- dat$X
    control_genes <- dat$ctl

    method_list            <- list()
    method_list$ols        <- ols(Y = Y, X = X)

    ## control gene methods --------------------------------------------------
    method_list$ruv2         <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes)
    method_list$ruv3_nomult  <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes,
                                     multiplier = FALSE)
    ## method_list$ruv3_mult    <- ruv3(Y = Y, X = X, num_sv = num_sv,
    ##                                  control_genes = control_genes,
    ##                                  multiplier = TRUE)
    method_list$ruv4         <- ruv4(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes)
    method_list$ruv4_rsvar   <- ruv4_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                  control_genes = control_genes)
    method_list$catenc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                        control_genes = control_genes,
                                        calibrate = FALSE)
    method_list$catenc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                        control_genes = control_genes,
                                        calibrate = TRUE)
    method_list$ruv4v_norm   <- vruv4(Y = Y, X = X,
                                      num_sv = num_sv,
                                      control_genes = control_genes,
                                      likelihood = "normal")
    method_list$ruv4v_t      <- vruv4(Y = Y, X = X,
                                      num_sv = num_sv,
                                      control_genes = control_genes,
                                      likelihood = "t")

    ## non control gene methods ----------------------------------------------
    method_list$sva          <- sva(Y = Y, X = X, num_sv = num_sv)
    method_list$caterr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
    method_list$caterr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)

    ash_list <- lapply(method_list, FUN = do_ash)
    ash_list$leapp <- NULL
    ash_list$catenc_cal <- NULL
    ash_list$caterr_cal <- NULL
    ash_list$mouthwash_norm <- mouthwash(Y = Y, X = X,
                                         num_sv = num_sv,
                                         likelihood = "normal",
                                         scale_var = TRUE)
    ash_list$mouthwash_t <- mouthwash(Y = Y, X = X,
                                      num_sv = num_sv,
                                      likelihood = "t")



    qvalue_list <- lapply(method_list, FUN = do_qvalue)




    ash_lfdr <- sapply(ash_list, FUN = get_lfdr)
    pvalmat  <- sapply(method_list, FUN = get_pvalues)
    colnames(ash_lfdr) <- paste0("ash_", colnames(ash_lfdr))
    colnames(pvalmat) <- paste0("pval_", colnames(pvalmat))
    plist[[tissue_index]] <- cbind(ash_lfdr, pvalmat)

    ash_betahat <- sapply(ash_list, FUN = get_betahat)
    reg_betahat <- sapply(method_list, FUN = get_betahat)
    betahat_list[[tissue_index]] <- cbind(ash_betahat, reg_betahat)
    colnames(betahat_list[[tissue_index]]) <- colnames(plist[[tissue_index]])

    num_methods <- length(ash_list) + length(qvalue_list)
    if (tissue_index == 1) {
        pi0mat <- matrix(NA, nrow = length(tissue_vec), ncol = num_methods)
        colnames(pi0mat) <- colnames(plist[[1]])
    }

    pi0mat[tissue_index, ] <- c(sapply(ash_list, get_pi0hat),
                                sapply(qvalue_list, get_pi0hat))

}

saveRDS(object = plist, file = "../Output/gtex_fits/plist.Rds")
saveRDS(object = betahat_list, file = "../Output/gtex_fits/betahat_list.Rds")
saveRDS(object = pi0mat, file = "../Output/gtex_fits/pi0mat.Rds")
