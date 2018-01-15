# Hopefully these will be the final sims. I've run these too many
# times and I'm really really bored.

one_rep <- function(new_params, current_params) {

  # I'm assuming that the current working directory is the parent
  # directory of this file.
  source("./Code/nc_adjustment_methods.R")
  source("./Code/non_nc_methods.R")
  args_val <- append(current_params, new_params)
  set.seed(new_params$current_seed)

  ## Choose all of the genes because already got top expressed
  stopifnot(args_val$Ngene == ncol(args_val$mat))
  d_out <- seqgendiff::poisthin(mat = args_val$mat,
                                nsamp = args_val$Nsamp,
                                ngene = args_val$Ngene,
                                skip_gene = args_val$skip_gene,
                                signal_params = list(mean = 0, sd = args_val$log2foldsd),
                                gvec = rep(TRUE, length(args_val$Ngene)),
                                gselect = "custom",
                                prop_null = args_val$nullpi,
                                alpha = 0)

  which_null <- abs(d_out$beta) < 10 ^ -6
  nnull         <- sum(which_null)
  control_genes <- which_null
  control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

  beta_true <- d_out$beta

  X <- d_out$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(d_out$Y + 1)

  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

  method_list      <- list()
  method_list$ols  <- ols(Y = Y, X = X)

  ## control gene methods --------------------------------------------------
  ruv2_o <- ruv2_simp(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
  method_list$ruv2  <- limma_adjust(ruv2_o)
  method_list$ruv3 <- ruv3_limma_pre(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes)
  method_list$cate  <- cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv,
                                               control_genes = control_genes,
                                               calibrate = FALSE)
  method_list$cate_madcal <- cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv,
                                                  control_genes = control_genes,
                                                  calibrate = TRUE)
  method_list$cate_nccal <- ctl_adjust(obj = method_list$cate, control_genes = control_genes)

  ## non control gene methods ----------------------------------------------
  method_list$sva    <- sva_voom(Y = Y, X = X, num_sv = num_sv)
  method_list$caterr <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
  method_list$caterr_cal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)

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
                                    outputlevel = 3, alpha = 0)
      return_list <- list()
      return_list$betahat <- ashout$result$PosteriorMean
      return_list$lfdr    <- ashout$result$lfdr
      return_list$pi0hat  <- ashr::get_pi0(ashout)
      return(return_list)
    }
  }

  do_qvalue <- function(args) {
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

  get_pvalue <- function(args) {
    if(!is.null(args$pvalues)) {
      return(list(pvalues = args$pvalues))
    } else {
      pvalues <- 2 * stats::pt(-abs(args$betahat / args$sebetahat), df = args$df)
      return(list(pvalues = c(pvalues)))
    }
  }

  plist <- lapply(method_list, get_pvalue)

  ## Ash-like procedures --------------------------------------------------------
  ash_list <- lapply(method_list, FUN = do_ash)
  ash_list$mouthwash_norm <- mouthwash(Y = Y, X = X, num_sv = num_sv,
                                       likelihood = "normal", alpha = 0, scale_var = TRUE)
  ash_list$backwash <- backwash(Y = Y, X = X, num_sv = num_sv, alpha = 0, scale_var = TRUE)

  ## Qvalue procedures ---------------------------------------------------------
  qvalue_list <- lapply(plist, FUN = do_qvalue)


  get_mse <- function(args, beta_true) {
    mean((args$betahat - beta_true) ^ 2, na.rm = TRUE)
  }

  get_auc <- function(args, which_null) {
    ## Choose p-values or lfdr's for predictors
    if (!is.null(args$pvalues)) {
      pred <- args$pvalues
    } else if (!is.null(args$lfdr)) {
      pred <- args$lfdr
    }

    ## in all-null setting
    if (sum(which_null) == length(which_null)) {
      return(NA)
    }

    pROC::roc(predictor = pred, response = which_null)$auc
  }

  ## pi0hat ----------------------------------------------------------
  ash_pi0 <- sapply(ash_list, FUN = function(args) args$pi0hat)
  names(ash_pi0) <- paste0("pi0_ash_", names(ash_pi0))
  qvalue_pi0 <- sapply(qvalue_list, FUN = function(args) args$pi0hat)
  names(qvalue_pi0) <- paste0("pi0_qvalue_", names(qvalue_pi0))
  pi0_vec <- c(ash_pi0, qvalue_pi0)

  ## auc ------------------------------------------------------------
  ash_auc <- sapply(ash_list, FUN = get_auc, which_null = which_null)
  names(ash_auc) <- paste0("auc_ash_", names(ash_auc))
  pvalue_auc <- sapply(plist, FUN = get_auc, which_null = which_null)
  names(pvalue_auc) <- paste0("auc_pvalue_", names(pvalue_auc))

  auc_vec <- c(ash_auc, pvalue_auc)

  ## mse ------------------------------------------------------------
  ash_mse <- sapply(ash_list, FUN = get_mse, beta_true = beta_true)
  names(ash_mse) <- paste0("mse_ash_", names(ash_mse))
  reg_mse <- sapply(method_list, FUN = get_mse, beta_true = beta_true)
  names(reg_mse) <- paste0("mse_reg_", names(reg_mse))

  mse_vec <- c(ash_mse, reg_mse)
  return_vec <- c(mse_vec, auc_vec, pi0_vec)
  return(return_vec)
}

itermax <- 500
seed_start <- 2222

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(6, 10, 20, 40)
ncontrol_seq <- c(10, 100)

par_vals <- expand.grid(list((1 + seed_start):(itermax + seed_start),
                             nullpi_seq, Nsamp_seq, ncontrol_seq))
colnames(par_vals) <- c("current_seed", "nullpi", "Nsamp", "ncontrols")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
  par_list[[list_index]] <- list()
  for (inner_list_index in 1:ncol(par_vals)) {
    par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
    names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
  }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 0.8
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## Create muscle_mat with most expressed genes
mat <- t(as.matrix(read.csv("./Output/gtex_tissue_gene_reads_v6p/muscle.csv",
                            header = TRUE)[, -c(1,2)]))
args_val$mat <- mat[, order(apply(mat, 2, median), decreasing = TRUE)[1:args_val$Ngene]]
rm(mat)

# Number of threads to use.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  nc <- 1
} else {
  nc <- args[1]
}

# If on your own computer, use this.
library(parallel)
cl   <- makeCluster(nc)
sout <- t(parallel::parSapply(cl = cl, par_list, FUN = one_rep,
                              current_params = args_val))
stopCluster(cl)

saveRDS(cbind(par_vals, sout), "./Output/sims_out/sims_out.RDS")





