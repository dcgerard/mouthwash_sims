library(snow)

one_rep <- function(new_params, current_params) {
    return_vec <- tryCatch(expr = {
        source("./Code/data_generators.R")
        source("./Code/adjustment_methods.R")
        args_val <- append(current_params, new_params)
        set.seed(new_params$current_seed)
        d_out <- datamaker_counts_only(args_val)
        which_null <- d_out$meta$null
        control_genes <- as.logical(which_null)
        nnull         <- sum(control_genes)
        control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

        beta_true <- rep(0, length = args_val$Ngene)
        beta_true[!which_null] <- d_out$meta$true_log2foldchange

        X <- as.matrix(model.matrix(~d_out$input$condition))
        colnames(X) <- c("Intercept", "Treatment")
        Y <- t(log2(as.matrix(d_out$input$counts + 1)))


        num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

        method_list            <- list()
        method_list$ols        <- ols(Y = Y, X = X)

        ## control gene methods --------------------------------------------------
        method_list$ruv2         <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes)
        method_list$ruv3_nomult  <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes,
                                         multiplier = FALSE)
        method_list$ruv3_mult    <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes,
                                         multiplier = TRUE)
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

        ## LEAPP SUCKS! Slow and always with the BUGS!
        ## method_list$leapp <- leapp(Y = Y, X = X, num_sv = num_sv)


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
                                              outputlevel = 3)
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
        ash_list$backwash    <- backwash(Y = Y, X = X, num_sv = num_sv)


        qvalue_list <- lapply(method_list, FUN = do_qvalue)


        get_mse <- function(args, beta_true) {
            mean((args$betahat - beta_true) ^ 2)
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
        pvalue_auc <- sapply(method_list, FUN = get_auc, which_null = which_null)
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

    }, error = function(e){rep(NA, 78)})
    return(return_vec)
}

itermax <- 500
seed_start <- 2222

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(3, 5, 10, 20)
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
args_val$tissue       <- "muscle"
args_val$path         <- "./Output/gtex_tissue_gene_reads_v6p/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## one_rep(par_list[[3]], args_val)

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores() - 2)
sout <- t(snow::parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)


## save(sout, file = "general_sims2.Rd")
mse_mat <- cbind(par_vals, sout[, 1:27])
auc_mat <- cbind(par_vals, sout[, 28:54])
pi0_mat <- cbind(par_vals, sout[, 55:81])
write.csv(mse_mat, file = "./Output/sims_out/mse_mat2.csv", row.names = FALSE)
write.csv(auc_mat, file = "./Output/sims_out/auc_mat2.csv", row.names = FALSE)
write.csv(pi0_mat, file = "./Output/sims_out/pi0_mat2.csv", row.names = FALSE)
