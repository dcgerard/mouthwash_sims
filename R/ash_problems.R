library(ggplot2)

source("../Code/data_generators.R")
source("../Code/adjustment_methods.R")
## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "../Output/gtex_tissue_gene_reads_v6p/"
args_val$Ngene        <- 10000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0
args_val$nullpi       <- 1
args_val$Nsamp        <- 3
args_val$poisthin     <- FALSE


set.seed(249)

d_out <- datamaker_counts_only(args_val)

X <- as.matrix(model.matrix(~d_out$input$condition))
Y <- t(log2(as.matrix(d_out$input$counts + 1)))

vout <- limma::voom(d_out$input$counts)
lout <- limma::lmFit(vout, design = X)
eout <- limma::ebayes(lout)
betahat <- lout$coefficients[, 2]
sebetahat <- sqrt(eout$s2.post) * lout$stdev.unscaled[, 2]
df <- eout$df.total

ols_out <- ols(Y = Y, X = X)

ash_limma <- ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = df)
ash_ols   <- ashr::ash.workhorse(betahat = ols_out$betahat, sebetahat = ols_out$sebetahat,
                                 df = ols_out$df)

ashr::get_pi0(ash_ols)
ashr::get_pi0(ash_limma)

mouth_out <- vicar::mouthwash(Y = Y, X = X, cov_of_interest = 2)
mouth_out$pi0


lfdr_df <- data.frame(OLSASH = ash_ols$result$lfdr,
                      VLEASH = ash_limma$result$lfdr,
                      MOUTHWASH = mouth_out$result$lfdr)
longdat <- reshape2::melt(lfdr_df, id.vars = NULL)
longdat$variable

longdat$variable <- stringr::str_replace(longdat$variable, "OLSASH", "OLS + ASH")
longdat$variable <- stringr::str_replace(longdat$variable, "VLEASH", "voom-limma-ebayes + ASH")
longdat$variable <- factor(longdat$variable, levels = c("OLS + ASH", "voom-limma-ebayes + ASH",
                                                        "MOUTHWASH"))
names(longdat) <- c("Method", "lfdr")
pdf(file = "../Output/figures/ash_fail.pdf", family = "Times", width = 6.5, height = 2.5, colormodel = "cmyk")
ggplot(data = longdat, mapping = aes(x = lfdr, color = I("black"), fill = I("white"))) +
    facet_wrap(~Method) +
    geom_histogram(bins = 20) +
    theme_bw() + ylab("Count") +
    theme(strip.background = element_rect(fill="white"))
dev.off()
