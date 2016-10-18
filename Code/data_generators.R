################################
## Filename: big_sims.R
## Created by: David Gerard
## Created on: 02/19/2016
## Synopsis: SUCCOTASH (MOUTHWASH) with FLASH.
################################

pois_thin <- function(Nsamp, nullpi, path, tissue = "muscle", ncontrol = 100, Ngene = 1000) {
    ## these do not change
    args_val              <- list()
    args_val$log2foldsd   <- 1
    args_val$tissue       <- tissue
    args_val$path         <- path
    args_val$Ngene        <- Ngene
    args_val$log2foldmean <- 0
    args_val$skip_gene    <- 0
    args_val$Nsamp        <- Nsamp
    args_val$nullpi       <- nullpi

    if (nullpi != 1) {
        args_val$poisthin <- TRUE
    }

    d_out                  <- datamaker_counts_only(args_val)
    which_null             <- d_out$meta$null
    control_genes          <- as.logical(which_null)
    nnull                  <- sum(control_genes)
    control_genes[control_genes][sample(1:nnull, size = nnull - ncontrol)] <- FALSE
    beta_true              <- rep(0, length = args_val$Ngene)
    if (nullpi != 1) {
        beta_true[!which_null] <- d_out$meta$true_log2foldchange
    }
    X                      <- as.matrix(model.matrix(~d_out$input$condition))
    colnames(X)            <- c("Intercept", "Treatment")
    Y                      <- t(log2(as.matrix(d_out$input$counts + 1)))
    num_sv                 <- sva::num.sv(t(Y), mod = X, method = "be")
    return(list(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes,
                which_null = which_null, beta_true = beta_true))
}


rep_pois_thin <- function(nsim, Nsamp, nullpi, path, tissue = "muscle",
                          ncontrol = 100, Ngene = 1000) {
    return_list <- list()
    for (index in 1:nsim) {
        return_list[[index]] <- pois_thin(Nsamp = Nsamp, nullpi = nullpi, path = path,
                                          tissue = tissue, ncontrol = ncontrol, Ngene = Ngene)
    }
    return(return_list)
}




datamaker_counts_only <- function(args) {
    dfargs <- default_datamaker_args(args)

    ## rawdata1 <- readtissue(dfargs$path, dfargs$tissue[1])
    rawdata1 <- read.csv(paste0(dfargs$path, dfargs$tissue[1], ".csv"), header = TRUE)[, -c(1,2)]
    if (length(dfargs$tissue) > 1) {
        ## rawdata2 <- readtissue(dfargs$path, dfargs$tissue[2])
        rawdata2 <- read.csv(paste0(dfargs$path, dfargs$tissue[2], ".csv"), header = TRUE)[, -c(1,2)]

        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp <- min(dim(rawdata1)[2], dim(rawdata2)[2])
        }
        if (dim(rawdata1)[2] < dfargs$Nsamp | dim(rawdata2)[2] < dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }

        if (dfargs$nullpi == 0) {
            ## All genes are alternatives
            temp1 <- selectsample(rawdata1, dfargs$Nsamp, dfargs$breaksample)
            counts1 <- temp1$counts
            subsample1 <- temp1$subsample
            rm(temp1)
            temp2 <- selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 <- temp2$counts
            subsample2 <- temp2$subsample
            rm(temp2)

            counts <- cbind(counts1, counts2)
            subsample <- cbind(subsample1, subsample2)
        } else {
            ## Some genes are nulls, some are alternatives
            temp1 <- selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
            counts1 <- temp1$counts
            subsample1 <- temp1$subsample
            rm(temp1)
            temp2 <- selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts2 <- temp2$counts
            subsample2 <- temp2$subsample
            rm(temp2)
            counts <- cbind(counts1, counts2)
            subsample <- cbind(subsample1, subsample2)
        }
    } else {
        if (is.null(dfargs$Nsamp)) {
            dfargs$Nsamp <- floor(dim(rawdata1)[2] / 2)
        }
        if (dim(rawdata1)[2] < 2 * dfargs$Nsamp) {
            stop("Not enough samples in the raw dataset!")
        }


        temp <- selectsample(rawdata1, 2 * dfargs$Nsamp, dfargs$breaksample)
        counts <- temp$counts
        subsample <- temp$subsample
        rm(temp)
    }

    ## Remove genes without any reads
    subsample <- subsample[apply(counts, 1, sum) > 0, ]
    counts <- counts[apply(counts, 1, sum) > 0, ]

    ## Take the top Ngene high-expressed genes
    if (!is.null(dfargs$Ngene)) {
        dfargs$Ngene <- min(dfargs$Ngene, dim(counts)[1])
        subsample <- subsample[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
        counts <- counts[sort(order(rowSums(counts), decreasing = TRUE)[(dfargs$skip_genes + 1):(dfargs$Ngene + dfargs$skip_genes)]), ]
    }
    dfargs$Ngene <- dim(counts)[1]

    ## Model's design: Nsamp samples for group A and Nsamp samples for group B
    condition <- factor(rep(1:2, each = dfargs$Nsamp))
    design <- stats::model.matrix(~condition)

    ## Ground truth of null hypotheses: beta_g=0
    null <- rep(0, dfargs$Ngene)
    null[sample(dfargs$Ngene, round(dfargs$Ngene * dfargs$nullpi))] <- 1
    ## default is that dfargs$nullpi is 1, so all genes are null genes -- dcg

    ## Poisson thinning (optional)
    pois_out <- pois_thinning(counts, dfargs, null) ## does nothing if args$poisthin == FALSE -- dcg
    counts <- pois_out$counts
    true_log2foldchange <- pois_out$log2foldchanges

    ## Mix null and alternative genes from different samples (optional)
    mix <- mix_sample(counts, dfargs, null, subsample) ## only matters if not all genes are null genes -- dcg
    counts <- mix$counts
    subsample <- mix$subsample



    input <- list(counts = counts, condition = condition)
    meta <- list(null = null, true_log2foldchange = true_log2foldchange, dfargs = dfargs, subsample = subsample)

    return(list(meta = meta,input = input))
}

default_datamaker_args <- function(args) {
    ## poisthin: flag of Poisson thinning
    if (is.null(args$poisthin)) {
        args$poisthin <- FALSE
    }

    ## number of top genes to skip
    if (is.null(args$skip_genes)) {
        args$skip_genes <- 0
    }

    ## log2foldmean, log2foldsd: Poisson thinning params
    if (args$poisthin == TRUE) {
        if (is.null(args$log2foldmean)) {
            args$log2foldmean <- 0
        }
        if (is.null(args$log2foldsd)) {
            args$log2foldsd <- 1
        }
    }

    ## breaksample: flag of each gene randomly select samples
    if (is.null(args$breaksample)) {
        args$breaksample <- FALSE
    }


    ## nullpi: proportion of null genes
    if (is.null(args$nullpi)) {
        if (args$poisthin == TRUE) {
            args$nullpi <- 0.9
        } else if (length(args$tissue) == 1) {
            args$nullpi <- 1
        } else if (length(args$tissue) > 1) {
            args$nullpi <- 0
        } else if (is.null(args$tissue)) {
            args$nullpi <- 1
        }
    }

    ## pseudocounts: add pseudocounts to count matrix
    if (is.null(args$pseudocounts)) {
        args$pseudocounts <- 1
    }

    if(is.null(args$sig_diag)) {
        args$sig_diag <- rep(1, args$Ngene)
    }

    if(is.null(args$sig_alpha)) {
        args$sig_alpha <- 1
    }

    if(is.null(args$beta0)) {
         args$beta0 <- 10
    }

    if(is.null(args$get_null)) {
        args$get_null <- TRUE
    }

    if (is.null(args$alt_type)) {
        args$alt_type <- "normal"
    }

    if (is.null(args$log2fold_inflate_beta)) {
        args$log2fold_inflate_beta <- 1
    }

    return(args)
}

selectsample <- function(counts, Nsamp, breaksample) {
    if (breaksample == FALSE) {
        subsample <- sample(1:dim(counts)[2], Nsamp)
        counts <- counts[, subsample]
        subsample <- t(matrix(rep(subsample, dim(counts)[1]), ncol = dim(counts)[1]))
    } else {
        temp <- t(apply(counts, 1, sampleingene, Nsamp = Nsamp))
        counts <- temp[, 1:Nsamp]
        subsample <- temp[, (Nsamp + 1):(2 * Nsamp)]
    }
    return(list(counts = counts, subsample = subsample))
}

pois_thinning <- function(counts, args, null) {
    if (args$poisthin == TRUE) {
        if (args$alt_type == "normal") {
            log2foldchanges <- rnorm(sum(!null), mean = args$log2foldmean, sd = args$log2foldsd)
        } else if (args$alt_type == "mixnorm") {
            if (is.null(args$pi_vals)) {
                stop("args$alt_type = \"mixnorm\" but pi_vals not specified")
            } else if (is.null(args$tau_seq)) {
                stop("args$alt_type = \"mixnorm\" but tau_seq not specified")
            } else if (is.null(args$mu_seq)) {
                stop("args$alt_type = \"mixnorm\" but mu_seq not specified")
            }
            log2foldchanges <- rmixnorm(pi_vals = args$pi_vals,
                                        mu_seq = args$mu_seq,
                                        tau_seq = args$tau_seq,
                                        p = sum(!null))
            cat("yay!\n")
        }
        log2foldchanges <- log2foldchanges * args$log2fold_inflate_beta ## inflation defaults to 1.

        foldchanges <- 2 ^ log2foldchanges

        ## thin group A
        counts[which(!null)[log2foldchanges > 0], 1:args$Nsamp] <-
            matrix(rbinom(sum(log2foldchanges >
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges >
            0], 1:args$Nsamp])), prob = rep(1 / foldchanges[log2foldchanges > 0], args$Nsamp)),
            ncol = args$Nsamp)
        ## thin group B
        counts[which(!null)[log2foldchanges < 0], (args$Nsamp + 1):(2 * args$Nsamp)] <-
            matrix(rbinom(sum(log2foldchanges <
            0) * args$Nsamp, size = c(as.matrix(counts[which(!null)[log2foldchanges <
            0], (args$Nsamp + 1):(2 * args$Nsamp)])), prob = rep(foldchanges[log2foldchanges <
            0], args$Nsamp)), ncol = args$Nsamp)

    } else {
        log2foldchanges <- rep(0, length = args$Ngene)
    }
    return(list(counts = counts, log2foldchanges = log2foldchanges))
}

mix_sample <- function(counts, args, null, subsample) {
    if (args$nullpi < 1 & args$nullpi > 0 & args$breaksample == TRUE) {
        newcounts <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newcounts[as.logical(null), ] <- counts[as.logical(null), 1:(2 * args$Nsamp)]
        newcounts[!null, ] <- counts[!null, c(1:args$Nsamp, (2 * args$Nsamp + 1):(3 *
            args$Nsamp))]
        counts <- newcounts
        newsubsample <- matrix(rep(0, args$Ngene * 2 * args$Nsamp), nrow = args$Ngene)
        newsubsample[as.logical(null), ] <- subsample[as.logical(null), 1:(2 * args$Nsamp)]
        newsubsample[!null, ] <- subsample[!null, c(1:args$Nsamp, (2 * args$Nsamp +
            1):(3 * args$Nsamp))]
        subsample <- newsubsample
        rm(newcounts)
        rm(newsubsample)
    }
    return(list(counts = counts, subsample = subsample))
}

fit_ash <- function(out_obj) {
    if(is.null(out_obj$sebetahat)) { ## some methods do not return sebetahat.
        return()
    } else {
        ash_out <- ashr::ash(betahat = out_obj$betahat, sebetahat = out_obj$sebetahat,
                             df = out_obj$df)
        return(ash_out)
    }
}

#' Get wald p-values and adjust by benjamini and hochberg.
fit_freq_methods <- function(out_obj) {
    if (is.null(out_obj$pvalue)) { ## if given pvalues, use those
        if (is.null(out_obj$df)) {
            p_values <- 2 * (1 - pnorm(abs(out_obj$betahat / out_obj$sebetahat)))
        } else {
            p_values <- 2 * (1 - pt(abs(out_obj$betahat / out_obj$sebetahat), df = out_obj$df))
        }
    } else {
        p_values <- out_obj$pvalue
    }

    q_bh <- stats::p.adjust(p_values, method = "BH")

    if (all(is.na(p_values))) {
        q_storey <- NULL
    } else {
        q_storey <- qvalue::qvalue(p_values)
    }
    return(list(p_values = p_values, q_bh = q_bh, q_storey = q_storey))
}

extract_ashpi0 <- function(ash_obj) {
    return(ash_obj$fitted.g$pi[1])
}

extract_ashlfdr <- function(ash_obj) {
    return(ash_obj$lfdr)
}

extract_qvaluepi0 <- function(freq_obj) {
    return(freq_obj$q_storey$pi0)
}

extract_pvalues <- function(freq_obj) {
    return(freq_obj$p_values)
}


##' Simple ordinary least squares.
##'
##' @param log_counts A matrix of numerics. The responses.
##' @param condition A matrix of numerics. The predictors.
get_ols <- function(log_counts, condition) {
    limma_out <- limma::lmFit(log_counts, model.matrix(~condition))
    betahat <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df <- limma_out$df.residual
    return(list(betahat = betahat, sebetahat = sebetahat, df = df))
}

#' Generate a random sample from a mixture of normals.
#'
#' @param pi_vals The mixing proportions.
#' @param mu_seq The mixing means.
#' @param tau_seq The mixing standard deviations.
#' @param p The number of samples.
rmixnorm <- function (pi_vals, mu_seq, tau_seq, p)
{
    M <- length(pi_vals)
    beta <- rep(NA, length = p)
    which.mix <- sample(1:M, size = p, replace = TRUE, prob = pi_vals)
    for (index in 1:M) {
        current.ind <- which.mix == index
        n_m <- sum(current.ind)
        if (n_m > 0) {
            beta[current.ind] <- stats::rnorm(n = n_m, mean = mu_seq[index],
                sd = tau_seq[index])
        }
    }
    return(beta)
}



datamaker_theo <- function(Ngene, Nsamp, nullpi, ncontrols,
                           alpha_sd = 1, sigsd = 1, num_sv = 5) {
    Ntrt <- round(Nsamp / 2)
    ## X <- cbind(rep(1, length = Nsamp), c(rep(1, Ntrt), rep(0, Nsamp - Ntrt)))
    X <- matrix(rnorm(2 * Nsamp), nrow = Nsamp)
    ## Z <- matrix(NA, nrow = Nsamp, ncol = num_sv)
    ## for(index in 1:num_sv) {
    ##     Z[, index] <- sample(c(0, 1), size = Nsamp, replace = TRUE, prob = c(1/2, 1/2))
    ## }
    Z <- matrix(rnorm(num_sv * Nsamp), nrow = Nsamp)
    alpha <- matrix(rnorm(num_sv * Ngene, sd = alpha_sd), nrow = num_sv, ncol = Ngene)
    beta <- matrix(rnorm(2 * Ngene), nrow = 2, ncol = Ngene)

    which_null <- rep(FALSE, length = Ngene)
    which_null[sample(1:Ngene, round(nullpi * Ngene))] <- TRUE
    beta[2, which_null] <- 0

    E <- matrix(rnorm(Ngene * Nsamp, sd = sigsd), nrow = Nsamp)

    Y <- X %*% beta + Z %*% alpha + E

    retlist <- list()
    retlist$Y <- round(2 ^ Y * 256)
    retlist$X <- X
    retlist$beta <- beta
    retlist$Z <- Z
    retlist$alpha <- alpha
    retlist$E <- E
    retlist$which_null <- which_null
    return(retlist)
}


datamaker_change_ncovs <- function(Ngene, Nsamp, nullpi, ncontrols, ncovs,
                                   alpha_sd = 1, sigsd = 1, num_sv = 5) {
    Ntrt <- round(Nsamp / 2)
    ## X <- cbind(rep(1, length = Nsamp), c(rep(1, Ntrt), rep(0, Nsamp - Ntrt)))
    X <- matrix(rnorm(ncovs * Nsamp), nrow = Nsamp)
    ## Z <- matrix(NA, nrow = Nsamp, ncol = num_sv)
    ## for(index in 1:num_sv) {
    ##     Z[, index] <- sample(c(0, 1), size = Nsamp, replace = TRUE, prob = c(1/2, 1/2))
    ## }
    Z <- matrix(rnorm(num_sv * Nsamp), nrow = Nsamp)
    alpha <- matrix(rnorm(num_sv * Ngene, sd = alpha_sd), nrow = num_sv, ncol = Ngene)
    beta <- matrix(rnorm(2 * Ngene), nrow = ncovs, ncol = Ngene)

    which_null <- rep(FALSE, length = Ngene)
    which_null[sample(1:Ngene, round(nullpi * Ngene))] <- TRUE
    beta[, which_null] <- 0

    E <- matrix(rnorm(Ngene * Nsamp, sd = sigsd), nrow = Nsamp)

    Y <- X %*% beta + Z %*% alpha + E

    retlist <- list()
    retlist$Y <- round(2 ^ Y * 256)
    retlist$X <- X
    retlist$beta <- beta
    retlist$Z <- Z
    retlist$alpha <- alpha
    retlist$E <- E
    retlist$which_null <- which_null
    return(retlist)
}
