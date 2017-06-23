library(ggplot2)
library(dplyr)

source("./Code/data_generators.R")
source("./Code/adjustment_methods.R")
## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "./Output/gtex_tissue_gene_reads_v6p/"
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

ash_limma <- ashr::ash.workhorse(betahat = betahat, sebetahat = sebetahat, df = df[1])
ash_ols   <- ashr::ash.workhorse(betahat = ols_out$betahat, sebetahat = ols_out$sebetahat,
                                 df = ols_out$df)

ashr::get_pi0(ash_ols)
ashr::get_pi0(ash_limma)

library(sva)
num_sv <- sva::num.sv(dat = t(Y), mod = X)

cat("num_sv = ", num_sv, "\n")

mouth_out <- vicar::mouthwash(Y = Y, X = X, k = num_sv, cov_of_interest = 2)
mouth_out$pi0

back_out <- vicar::backwash(Y = Y, X = X, k = num_sv, cov_of_interest = 2)
back_out$pi0

lfdr_df <- data.frame(OLSASH = ash_ols$result$lfdr,
                      VLEASH = ash_limma$result$lfdr,
                      MOUTHWASH = mouth_out$result$lfdr,
                      BACKWASH = back_out$result$lfdr)
longdat <- reshape2::melt(lfdr_df, id.vars = NULL)

longdat$variable <- stringr::str_replace(longdat$variable, "OLSASH", "OLS + ASH")
longdat$variable <- stringr::str_replace(longdat$variable, "VLEASH", "VOOM + ASH")
longdat$variable <- factor(longdat$variable, levels = c("OLS + ASH", "VOOM + ASH",
                                                        "MOUTHWASH", "BACKWASH"))
names(longdat) <- c("Method", "lfdr")
pdf(file = "./Output/figures/ash_fail.pdf", family = "Times", width = 6.5, height = 2.5, colormodel = "cmyk")
ggplot(data = longdat, mapping = aes(x = lfdr, color = I("black"), fill = I("white"))) +
    facet_grid(.~Method) +
    geom_histogram(bins = 20) +
    theme_bw() + ylab("Count") +
    theme(strip.background = element_rect(fill="white")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
dev.off()

smalldat <- filter(longdat, Method == "OLS + ASH" | Method  == "VOOM + ASH")

pdf(file = "./Output/figures/ash_fail_small.pdf", family = "Times", width = 4, height = 2.5, colormodel = "cmyk")
ggplot(data = smalldat, mapping = aes(x = lfdr, color = I("black"), fill = I("lightgrey"))) +
  facet_grid(.~Method) +
  geom_histogram(bins = 20) +
  theme_bw() + ylab("Count") +
  theme(strip.background = element_rect(fill="white")) +
  xlim(0, 1)
dev.off()

smalldat <- filter(longdat, Method == "MOUTHWASH")
pdf(file = "./Output/figures/ash_fail_small_mouth.pdf", family = "Times", width = 4, height = 2.5, colormodel = "cmyk")
ggplot(data = smalldat, mapping = aes(x = lfdr, color = I("black"), fill = I("lightgrey"))) +
  facet_grid(.~Method) +
  geom_histogram(bins = 20) +
  theme_bw() + ylab("Count") +
  theme(strip.background = element_rect(fill="white")) +
  xlim(0, 1)
dev.off()

