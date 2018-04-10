## This is the same as fit_mouthwash.R but using the control genes from Lin et al.
## I also don't run methods that don't require control genes here.

suppressMessages(library(tidyverse))
library(stringr)

# I'm assuming that the current working directory is the parent
# directory of this file.
source("./Code/nc_adjustment_methods.R")
source("./Code/non_nc_methods.R")
source("./Code/mosek.R")

# Although not strictly necessary, Rmosek dramatically speeds up the
# mouthwash and backwash methods. This will give an error if the
# Rmosek package is not available.
# library(Rmosek)
# library(REBayes)
# test_mosek()

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

nseq <- rep(NA, length = length(tissue_vec))

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
        return_list$betahat <- ashr::get_pm(ashout)
        return_list$lfdr    <- ashr::get_lfdr(ashout)
        return_list$pi0hat  <- ashr::get_pi0(ashout)
        return(return_list)
    }
}

do_qvalue <- function(args) {
  qout <- qvalue::qvalue(p = args)
  return_list        <- list()
  return_list$lfdr   <- qout$lfdr
  return_list$pi0hat <- qout$pi0
  return(return_list)
}

get_lfdr <- function(args) {
    args$lfdr
}

get_pvalues <- function(args) {
  if(!is.null(args$pvalues)) {
    return(list(pvalues = c(args$pvalues)))
  } else {
    pvalues <- 2 * stats::pt(-abs(args$betahat / args$sebetahat), df = args$df)
    return(list(pvalues = c(pvalues)))
  }
}

get_pi0hat <- function(args) {
    args$pi0hat
}

get_betahat <- function(args) {
    args$betahat
}


plist <- list()
betahat_list <- list()
num_sv_seq <- rep(NA, length = length(tissue_vec))
lin_ctl <- readRDS("./Output/cleaned_gtex_data/lin_ctl.Rds")
num_sv_seq <- readRDS("./Output/gtex_fits/num_sv.Rds")
for(tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    cat("Tissue",tissue_index,"=",current_tissue,"\n")

    dat <- readRDS(paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))
    dat$ctl <- lin_ctl[[tissue_index]]

    onsex <- dat$chrom == "X" | dat$chrom == "Y"
    onsex[is.na(onsex)] <- FALSE
    dat$ctl[onsex] <- FALSE
    nseq[tissue_index] <- ncol(dat$Y)

    num_sv <- num_sv_seq[tissue_index]

    Y <- t(dat$Y)
    X <- dat$X
    control_genes <- dat$ctl

    method_list      <- list()

    ## control gene methods --------------------------------------------------
    ## only do control gene methods here because non-control gene methods are done in fit_mouthwash.R
    cat(" - Running ruv2_o.\n")
    out <- system.time(ruv2_o <- ruv2_simp(Y = Y, X = X, num_sv = num_sv,
                                           control_genes = control_genes))
    cat("   Computation took",out["elapsed"],"seconds.\n")
    cat(" - Running ruv2.\n")
    out <- system.time(method_list$ruv2  <- limma_adjust(ruv2_o))
    cat("   Computation took",out["elapsed"],"seconds.\n")
    cat(" - Running ruv3.\n")
    out <- system.time(
      method_list$ruv3 <- ruv3_limma_pre(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes))
    cat("   Computation took",out["elapsed"],"seconds.\n")
    cat(" - Running cate.\n")
    out <- system.time(
      method_list$cate  <- cate_simp_nc_correction(Y = Y,X = X,
                             num_sv = num_sv,control_genes = control_genes))
    cat("   Computation took",out["elapsed"],"seconds.\n")
    cat(" - Running cate_madcal.\n")
    out <- system.time(method_list$cate_madcal <-
      cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv,
                              control_genes = control_genes,
                              calibrate = TRUE))
    cat("   Computation took",out["elapsed"],"seconds.\n")
    cat(" - Running cate_nccal.\n")
    out <- system.time(method_list$cate_nccal <-
                         ctl_adjust(obj = method_list$cate,
                                    control_genes = control_genes))
    cat("   Computation took",out["elapsed"],"seconds.\n")


    ash_list <- lapply(method_list, FUN = do_ash)
    ash_lfdr <- sapply(ash_list, FUN = get_lfdr)
    pvalmat  <- sapply(sapply(method_list, FUN = get_pvalues), c)
    colnames(pvalmat) <- stringr::str_replace(colnames(pvalmat), ".pvalues", "")
    qvalue_list <- apply(pvalmat, 2, FUN = do_qvalue)
    qvalue_lfdr <- sapply(qvalue_list, FUN = get_lfdr)

    colnames(qvalue_lfdr) <- paste0("qvalue_", colnames(qvalue_lfdr))
    colnames(ash_lfdr) <- paste0("ash_", colnames(ash_lfdr))
    colnames(pvalmat) <- paste0("pval_", colnames(pvalmat))
    plist[[tissue_index]] <- cbind(ash_lfdr, pvalmat, qvalue_lfdr)

    ash_betahat <- sapply(ash_list, FUN = get_betahat)
    reg_betahat <- sapply(method_list, FUN = get_betahat)
    betahat_list[[tissue_index]] <- cbind(ash_betahat, reg_betahat)
    colnames(betahat_list[[tissue_index]]) <- colnames(plist[[tissue_index]])[1:ncol(betahat_list[[tissue_index]])]

    num_methods <- length(ash_list) + length(qvalue_list)
    if (tissue_index == 1) {
        pi0mat <- matrix(NA, nrow = length(tissue_vec), ncol = num_methods)
        colnames(pi0mat) <- colnames(betahat_list[[1]])
    }

    pi0mat[tissue_index, ] <- c(sapply(ash_list, get_pi0hat),
                                sapply(qvalue_list, get_pi0hat))

}

saveRDS(object = plist, file = "./Output/gtex_fits_lin/plist.Rds")
saveRDS(object = betahat_list, file = "./Output/gtex_fits_lin/betahat_list.Rds")
saveRDS(object = pi0mat, file = "./Output/gtex_fits_lin/pi0mat.Rds")
saveRDS(object = num_sv_seq, "./Output/gtex_fits_lin/num_sv.Rds")

