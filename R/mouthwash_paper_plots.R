library(stringr)
library(reshape2)
library(ggplot2)
library(dplyr)

## AUC --------------------------------------------------------------------------------
auc <- read.csv("../Output/sims_out/auc_mat2.csv")
auc <- select(auc, -auc_ash_ruv3_mult)
auc <- select(auc, -auc_pvalue_ruv3_mult)
auc_names <- names(auc)
auc_names <- str_replace(auc_names, "auc_", "")
auc_names <- str_replace(auc_names, "ash_mouthwash", "MOUTHWASH")
auc_names <- str_replace(auc_names, "pvalue_", "")
auc_names <- str_replace(auc_names, "ash_(.+)", "\\1+ASH")
auc_names <- str_replace(auc_names, "ruv", "RUV")
auc_names <- str_replace(auc_names, "RUV4v", "CATEv")
auc_names <- str_replace(auc_names, "_rsvar", "v")
auc_names <- str_replace(auc_names, "sva", "SVA")
auc_names <- str_replace(auc_names, "cate", "CATE")
auc_names <- str_replace(auc_names, "_nocal", "")
auc_names <- str_replace(auc_names, "_nomult", "")
auc_names <- str_replace(auc_names, "_cal", "+Cal")
auc_names <- str_replace(auc_names, "_mult", "+Cal")
auc_names <- str_replace(auc_names, "^ols", "OLS")
auc_names <- str_replace(auc_names, "_norm", "(Norm)")
auc_names <- str_replace(auc_names, "_t", "(t)")
names(auc) <- auc_names

longdat <- reshape2::melt(auc, id.vars = 1:5)
longdat <- filter(longdat, nullpi != 1)
longdat$Nsamp <- longdat$Nsamp * 2
auc_mean <- longdat %>% group_by(nullpi, Nsamp, ncontrols, variable) %>%
    summarise(avg = mean(value))

## get ordering
subsamp <- filter(auc_mean, ncontrols == 10 & nullpi == 0.9 & Nsamp == 40)
fa_levels <- subsamp$variable[order(subsamp$avg, decreasing = TRUE)]
auc_mean$variable <- factor(auc_mean$variable, levels = fa_levels)
auc_mean$Nsamp <- factor(auc_mean$Nsamp)
colnames(auc_mean)[2] <- "Sample Size"

auc_mean$`ASH-like` <- str_detect(auc_mean$variable, "ASH")

## Set up lines for max auc
dummy_dat <- auc_mean %>% group_by(nullpi, `Sample Size`, ncontrols) %>% summarize(Max = max(avg))
names(dummy_dat)[2] <- "Sample Size "

## plot!
pdf(file = "../Output/figures/auc_ave.pdf", family = "Times", width = 6.5, height = 7, colormodel = "cmyk")
ggplot(data = auc_mean, mapping = aes(pch = `Sample Size`, y = avg, x = variable)) +
    facet_grid(nullpi ~ ncontrols) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
    theme(strip.background = element_rect(fill="white")) +
    theme(text = element_text(size= 10)) +
    ylab("Mean AUC") + xlab("Method") +
    geom_hline(data = dummy_dat, mapping = aes(yintercept = Max, lty = `Sample Size `),
               alpha = 0.5)
dev.off()


## Pi0hat ----------------------------------------------------------------------------
rm(list = ls())
pi0 <- read.csv("../Output/sims_out/pi0_mat2.csv")
pi0 <- select(pi0, -pi0_ash_ruv3_mult)
pi0 <- select(pi0, -pi0_qvalue_ruv3_mult)
pi0_names <- names(pi0)
pi0_names <- str_replace(pi0_names, "pi0_", "")
pi0_names <- str_replace(pi0_names, "ash_mouthwash", "MOUTHWASH")
pi0_names <- str_replace(pi0_names, "qvalue_(.+)", "\\1+qvalue")
pi0_names <- str_replace(pi0_names, "ash_(.+)", "\\1+ASH")
pi0_names <- str_replace(pi0_names, "ruv", "RUV")
pi0_names <- str_replace(pi0_names, "RUV4v", "CATEv")
pi0_names <- str_replace(pi0_names, "_rsvar", "v")
pi0_names <- str_replace(pi0_names, "sva", "SVA")
pi0_names <- str_replace(pi0_names, "cate", "CATE")
pi0_names <- str_replace(pi0_names, "_nocal", "")
pi0_names <- str_replace(pi0_names, "_nomult", "")
pi0_names <- str_replace(pi0_names, "_cal", "+Cal")
pi0_names <- str_replace(pi0_names, "_mult", "+Cal")
pi0_names <- str_replace(pi0_names, "^ols", "OLS")
pi0_names <- str_replace(pi0_names, "_norm", "(Norm)")
pi0_names <- str_replace(pi0_names, "_t", "(t)")
names(pi0) <- pi0_names

longdat <- reshape2::melt(pi0, id.vars = 1:5)
longdat$Nsamp <- longdat$Nsamp * 2
pi0_median <- longdat %>% group_by(nullpi, Nsamp, ncontrols, variable) %>%
    summarise(med = median(value))
pi0_median$Nsamp <- factor(pi0_median$Nsamp)

## set up ordering
subsamp <- filter(pi0_median, nullpi == 0.9 & Nsamp == 40 & ncontrols == 10)
flevels <- levels(pi0_median$variable)[order(abs(subsamp$med - 0.9))]
pi0_median$variable <- factor(pi0_median$variable, levels = flevels)

pdf(file = "../Output/figures/pi0_med.pdf", family = "Times", colormodel = "cmyk", height = 7, width = 6.5)
ggplot(data = pi0_median, mapping = aes(x = variable, y = med, pch = Nsamp)) +
    facet_grid(nullpi ~ ncontrols) +
    geom_point() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
    theme(strip.background = element_rect(fill="white")) +
    theme(text = element_text(size= 10)) +
    ylab("Median Estimate of Pi0") + xlab("Method")
dev.off()


## Original pi0hat 0.9 ----------------------------------
subdat <- filter(longdat, nullpi == 0.9)

## find ordering
dist_from_9 <- function(value){mean((value - 0.9) ^ 2)}
subsubdat <- filter(subdat, Nsamp == 40 & ncontrols == 10)
mse_dat <- subsubdat %>% group_by(variable) %>% summarize(mse = dist_from_9(value))
fa_labels <- mse_dat$variable[order(mse_dat$mse)]

subdat$variable <- factor(subdat$variable, levels = fa_labels)

pdf(file = "../Output/figures/pi0_box_09.pdf", family = "Times", colormodel = "cmyk", height = 7, width = 6.5)
ggplot(data = subdat, mapping = aes(y = value, x = variable)) +
    facet_grid(Nsamp ~ ncontrols) +
    geom_boxplot(outlier.size = 0.1, lwd = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
    theme(strip.background = element_rect(fill="white")) +
    theme(text = element_text(size= 10)) +
    ylab("Estimate of Pi0") + xlab("Method") +
    geom_hline(yintercept = 0.9, alpha = 0.5, lty = 2)
dev.off()



## Original pi0hat 0.5 -------------------------------------------------
subdat <- filter(longdat, nullpi == 0.5)

## find ordering
dist_from_5 <- function(value){mean((value - 0.5) ^ 2)}
subsubdat <- filter(subdat, Nsamp == 40 & ncontrols == 10)
mse_dat <- subsubdat %>% group_by(variable) %>% summarize(mse = dist_from_5(value))
fa_labels <- mse_dat$variable[order(mse_dat$mse)]

subdat$variable <- factor(subdat$variable, levels = fa_labels)

pdf(file = "../Output/figures/pi0_box_05.pdf", family = "Times", colormodel = "cmyk", height = 7, width = 6.5)
ggplot(data = subdat, mapping = aes(y = value, x = variable)) +
    facet_grid(Nsamp ~ ncontrols) +
    geom_boxplot(outlier.size = 0.1, lwd = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
    theme(strip.background = element_rect(fill="white")) +
    theme(text = element_text(size= 10)) +
    ylab("Estimate of Pi0") + xlab("Method") +
    geom_hline(yintercept = 0.5, alpha = 0.5, lty = 2)
dev.off()

## Original pi0hat 1 -------------------------------------------------
subdat <- filter(longdat, nullpi == 1)

## find ordering
dist_from_1 <- function(value){mean((value - 1) ^ 2)}
subsubdat <- filter(subdat, Nsamp == 40 & ncontrols == 10)
mse_dat <- subsubdat %>% group_by(variable) %>% summarize(mse = dist_from_1(value))
fa_labels <- mse_dat$variable[order(mse_dat$mse)]

subdat$variable <- factor(subdat$variable, levels = fa_labels)

pdf(file = "../Output/figures/pi0_box_1.pdf", family = "Times", colormodel = "cmyk", height = 7, width = 6.5)
ggplot(data = subdat, mapping = aes(y = value, x = variable)) +
    facet_grid(Nsamp ~ ncontrols) +
    geom_boxplot(outlier.size = 0.1, lwd = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
    theme(strip.background = element_rect(fill="white")) +
    theme(text = element_text(size= 10)) +
    ylab("Estimate of Pi0") + xlab("Method") +
    geom_hline(yintercept = 1, alpha = 0.5, lty = 2)
dev.off()
