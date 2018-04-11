## Measure computation time on 100 repetitions

source("./Code/nc_adjustment_methods.R")
source("./Code/non_nc_methods.R")

log2foldsd   <- 0.8
Ngene        <- 1000
log2foldmean <- 0
skip_gene    <- 0
Nsamp_seq    <- c(6, 10, 20, 40)
nullpi       <- 0.9
ncontrol     <- 100
itermax      <- 100

mat <- t(as.matrix(read.csv("./Output/gtex_tissue_gene_reads_v6p/muscle.csv",
                                   header = TRUE)[, -c(1,2)]))
mat <- mat[, order(apply(mat, 2, median), decreasing = TRUE)[1:Ngene]]

simdat <- expand.grid(seed = 1:itermax, Nsamp = Nsamp_seq)
simdat$ols_time       <- NA
simdat$ruv2_time      <- NA
simdat$ruv3_time      <- NA
simdat$catenc_time    <- NA
simdat$sva_time       <- NA
simdat$caterr_time    <- NA
simdat$mouthwash_time <- NA
simdat$backwash_time  <- NA

for (index in 1:nrow(simdat)) {
  set.seed(simdat$seed)
  cat("Index:", index, "\n")
  Nsamp <- simdat$Nsamp[index]
  ## Choose all of the genes because already got top expressed
  d_out <- seqgendiff::poisthin(mat = mat,
                                nsamp = Nsamp,
                                ngene = Ngene,
                                skip_gene = 0,
                                signal_params = list(mean = 0, sd = log2foldsd),
                                gvec = rep(TRUE, length(Ngene)),
                                gselect = "custom",
                                prop_null = nullpi,
                                alpha = 0)

  which_null <- abs(d_out$beta) < 10 ^ -6
  nnull         <- sum(which_null)
  control_genes <- which_null
  control_genes[control_genes][sample(1:nnull, size = nnull - ncontrol)] <- FALSE
  beta_true <- d_out$beta
  X <- d_out$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(d_out$Y + 1)

  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

  simdat$ols_time[index] <- system.time(ols(Y = Y, X = X))["elapsed"]

  ## control gene methods --------------------------------------------------
  simdat$ruv2_time[index] <- system.time(ruv2_simp(Y = Y, X = X, num_sv = num_sv,
                                                   control_genes = control_genes))["elapsed"]

  simdat$ruv3_time[index] <- system.time(ruv3_limma_pre(Y = Y, X = X, num_sv = num_sv,
                                                        control_genes = control_genes))["elapsed"]

  simdat$catenc_time[index] <- system.time(cate_simp_nc_correction(Y = Y, X = X, num_sv = num_sv,
                                                                   control_genes = control_genes,
                                                                   calibrate = FALSE))["elapsed"]

  ## non control gene methods ----------------------------------------------
  simdat$sva_time[index] <- system.time(sva_voom(Y = Y, X = X, num_sv = num_sv))["elapsed"]
  simdat$caterr_time[index] <- system.time(cate_rr(Y = Y, X = X, num_sv = num_sv,
                                                   calibrate = FALSE))["elapsed"]

  ## Ash-like procedures --------------------------------------------------------
  simdat$mouthwash_time[index] <- system.time(mouthwash(Y = Y, X = X, num_sv = num_sv,
                                                        likelihood = "normal", alpha = 0,
                                                        scale_var = TRUE))["elapsed"]
  simdat$backwash_time[index] <- system.time(backwash(Y = Y, X = X, num_sv = num_sv, alpha = 0,
                                                      scale_var = TRUE))["elapsed"]
}

saveRDS(simdat, "./Output/computation/comp_time.Rds")






