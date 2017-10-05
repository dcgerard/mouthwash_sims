library(stringr)
library(tidyverse)

do_topk <- function(lfdr, onsex, k) {
    return(sum(onsex[order(lfdr)][1:k], na.rm = TRUE))
}

tot_sig <- function(lfdr, onsex, alpha = 0.1) {
    sig <- onsex[lfdr < alpha]
    return(c(numsex = sum(sig, na.rm = TRUE), sumsig = length(sig)))
}

prop_sig <- function(lfdr, onsex, alpha = 0.1) {
    sig <- onsex[lfdr < alpha]
    propsig <- length(sig) / length(lfdr)
    propgen <- sum(sig, na.rm = TRUE) / length(sig)
    return(c(propsig = propsig, propgen = propgen))
}

replace_names <- function(x) {
  x <- stringr::str_replace(x, "pi0_", "")
  x <- stringr::str_replace(x, "auc_", "")
  x <- stringr::str_replace(x, "mse_", "")
  x <- stringr::str_replace(x, "ash_mouthwash", "MOUTHWASH")
  x <- stringr::str_replace(x, "ash_backwash", "BACKWASH")
  x <- stringr::str_replace(x, "ash_(.+)", "\\1+ASH")
  x <- stringr::str_replace(x, "qvalue_(.+)", "\\1+qvalue")
  x <- stringr::str_replace(x, "ruv", "RUV")
  x <- stringr::str_replace(x, "caterr_cal", "CATErr+MAD")
  x <- stringr::str_replace(x, "cate_nccal", "CATEnc+Cal")
  x <- stringr::str_replace(x, "caterr", "CATErr")
  x <- stringr::str_replace(x, "cate", "CATEnc")
  x <- stringr::str_replace(x, "_madcal", "+MAD")
  x <- stringr::str_replace(x, "pvalue_", "")
  x <- stringr::str_replace(x, "ols", "OLS")
  x <- stringr::str_replace(x, "sva", "SVA")
  x <- stringr::str_replace(x, "_norm", "")
}


tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")

good_tissue_labels <- c("Adipose Tissue", "Bladder", "Blood Vessel", "Breast",
                        "Colon", "Kidney", "Lung", "Nerve", "Pancreas",
                        "Skin", "Spleen", "Adrenal Gland", "Blood", "Brain",
                        "Esophagus", "Heart", "Liver", "Muscle", "Pituitary",
                        "Salivary Gland", "Small Intestine", "Stomach", "Thyroid")


## Change 20 if change number of methods!
top100 <- matrix(NA, nrow = length(tissue_vec), ncol = 20)

plist <- readRDS("./Output/gtex_fits/plist.Rds")
pi0mat <- readRDS("./Output/gtex_fits/pi0mat.Rds")
betahat_list <- readRDS("./Output/gtex_fits/betahat_list.Rds")
nseq <- rep(NA, length(tissue_vec))

for (tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    dat <- readRDS(paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))

    onsex <- dat$chrom == "X" | dat$chrom == "Y"
    control_genes <- dat$ctl

    pmat <- plist[[tissue_index]]
    pmat[control_genes, ] <- NA

    numzero <- apply(pmat, 2, function(x) sum(x == 0, na.rm = TRUE))
    if (any(numzero > 100)) {
        cat("LOOK HERE:", tissue_index, "\n")
    }

    betamat <- betahat_list[[tissue_index]]

    top100[tissue_index, ] <- apply(pmat[, 1:20], 2, do_topk, onsex = onsex, k = 200)
    nseq[tissue_index] <- ncol(dat$Y)

    countmat <- as.data.frame(t(apply(pmat, 2, tot_sig, onsex = onsex, alpha = 0.05)))
    propmat  <- as.data.frame(t(apply(pmat, 2, prop_sig, onsex = onsex, alpha = 0.05)))

    countmat <- countmat[!stringr::str_detect(rownames(countmat), "pval"), ]
    propmat  <- propmat[!stringr::str_detect(rownames(propmat), "pval"), ]

    namevec <- replace_names(rownames(countmat))

    countmat$method <- namevec
    propmat$method  <- namevec

    countmat$tissue <- tissue_vec[tissue_index]
    propmat$tissue  <- tissue_vec[tissue_index]

    if (tissue_index == 1) {
      colnames(top100) <- namevec
      allcountmat <- countmat
      allpropmat  <- propmat
    } else {
      allcountmat <- rbind(allcountmat, countmat)
      allpropmat  <- rbind(allpropmat, propmat)
    }
    cat(tissue_index, "\n")
}


## top 100 analysis
rownames(top100) <- paste0(good_tissue_labels, " (", nseq, ")")
top100p <- as_data_frame(top100 / apply(top100, 1, max)) %>%
  select(-`OLS+ASH`, -`OLS+qvalue`)
faorder <- order(apply(top100p, 2, median), decreasing = TRUE)
top100p$Tissue <- rownames(top100)




longdat <- gather(top100p, key = "Method", value = "Proportion", `RUV2+ASH`:`CATErr+qvalue`)
longdat$Tissue <- factor(longdat$Tissue, levels = rownames(top100)[order(nseq, decreasing = TRUE)])
longdat$Method <- factor(longdat$Method, levels = unique(longdat$Method)[faorder])


pdf(file = "./Output/figures/prop_max.pdf", height = 7.3, width = 6.5, family = "Times", color = "cmyk")
pl <- ggplot(data = longdat, mapping = aes(x = Method, y = Tissue, fill = Proportion)) +
    geom_raster() +
    scale_fill_gradient(high = "#ffffff",
                        low = "#000000",
                        na.value = "black",
                        guide = guide_colourbar(title = "Proportion\nfrom\nMaximum")) +
    theme_bw() +
    scale_size_manual(values=c(dot=2, no_dot=NA), guide="none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(pl)
dev.off()

setEPS()
postscript(file = "./Output/figures/prop_max.eps", height = 7.3, width = 6.5, family = "Times", color = "cmyk")
print(pl)
dev.off()

medpi0 <- data.frame(pi0hat = apply(pi0mat, 2, median))
rownames(medpi0) <- NULL
medpi0$Method <- namevec
medpi0 <- medpi0[, 2:1]
medpi0 <- medpi0[order(medpi0$pi0hat), ]

library(xtable)
print(xtable(medpi0), include.rownames = FALSE)


### Individual pmats -------------------------------------------------------------
rm(list = ls())
plist <- readRDS("./Output/gtex_fits/plist.Rds")
pi0mat <- readRDS("./Output/gtex_fits/pi0mat.Rds")
betahat_list <- readRDS("./Output/gtex_fits/betahat_list.Rds")

tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")

return_topk <- function(lfdr, onsex, k = 500) {
    olf <- order(lfdr)[1:k]
    return(cbind(lfdr[olf], onsex[olf]))
}

replace_names <- function(x) {
  x <- stringr::str_replace(x, "pi0_", "")
  x <- stringr::str_replace(x, "auc_", "")
  x <- stringr::str_replace(x, "mse_", "")
  x <- stringr::str_replace(x, "ash_mouthwash", "MOUTHWASH")
  x <- stringr::str_replace(x, "ash_backwash", "BACKWASH")
  x <- stringr::str_replace(x, "ash_(.+)", "\\1+ASH")
  x <- stringr::str_replace(x, "qvalue_(.+)", "\\1+qvalue")
  x <- stringr::str_replace(x, "ruv", "RUV")
  x <- stringr::str_replace(x, "caterr_cal", "CATErr+MAD")
  x <- stringr::str_replace(x, "cate_nccal", "CATEnc+Cal")
  x <- stringr::str_replace(x, "caterr", "CATErr")
  x <- stringr::str_replace(x, "cate", "CATEnc")
  x <- stringr::str_replace(x, "_madcal", "+MAD")
  x <- stringr::str_replace(x, "pvalue_", "")
  x <- stringr::str_replace(x, "ols", "OLS")
  x <- stringr::str_replace(x, "sva", "SVA")
  x <- stringr::str_replace(x, "_norm", "")
}

plot_now <- FALSE
for (tissue_index in 1:length(tissue_vec)) {

    current_tissue <- tissue_vec[tissue_index]
    dat <- readRDS(paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))
    onsex <- dat$chrom == "X" | dat$chrom == "Y"

    pmat <- as_data_frame(plist[[tissue_index]]) %>% select(contains("ash_"), contains("qvalue_"))

    aout <- lapply(pmat, return_topk, onsex = onsex)

    for (index in 1:length(aout)) {
        aout[[index]] <- as_data_frame(aout[[index]])
        aout[[index]]$method <- names(aout)[index]
        aout[[index]]$index  <- 1:nrow(aout[[index]])
    }

    longdat <- do.call(rbind, aout)
    longdat$tissue <- current_tissue
    names(longdat) <- c("lfdr", "onsex", "method", "index")

    if (tissue_index == 1) {
        biglongdat <- longdat
    } else {
        biglongdat <- rbind(biglongdat, longdat)
    }

    ## ggplot(data = longdat, mapping = aes(x = index, y = lfdr, color = method)) +
    ##    geom_line()

    if (plot_now) {
        par(ask = TRUE)
        p1 <- ggplot(data = longdat, mapping = aes(x = index,
                                                   y = lfdr,
                                                   color = as.factor(onsex), group = method)) +
            geom_point() +
            facet_wrap(~method)
        print(p1)
    }
}

colnames(biglongdat) <- c("lfdr", "onsex", "method", "index", "tissue")

dummydat <- biglongdat %>% group_by(method, index) %>% summarize(lfdr = median(lfdr), onsex = mean(onsex))



dummydat$method <- replace_names(dummydat$method)

names(dummydat) <- c("Method", "Rank", "lfdr", "Proportion")

biglongdat2 <- biglongdat
names(biglongdat2) <- c("lfdr", "onsex", "Method", "Rank", "Tissue")

biglongdat2$Method <- replace_names(biglongdat2$Method)

pdf(file = "./Output/figures/proponsex.pdf", height = 5.5, width = 6.5, family = "Times", color = "cmyk")
ggplot() + geom_point(data = dummydat, mapping = aes(x = Rank, y = lfdr,
                                                     color = Proportion), size = 0.5) +
    facet_wrap(~Method) +
    theme_bw() +
    ggthemes::scale_color_gradient_tableau(guide = guide_colourbar(title = "Proportion\non Sex\nChromosome"))+
                        #high = "black",
                        #low = "gray80",
    theme(strip.background = element_rect(fill="white"), text = element_text(size = 9)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Median lfdr")
dev.off()

pdf(file = "./Output/figures/lfdr_rank.pdf", height = 5.5, width = 6.5, family = "Times", color = "cmyk")
ggplot(data = biglongdat2, mapping = aes(y = lfdr, x = Rank, group = Tissue)) +
    geom_line(alpha = 1/4) +
    facet_wrap(~Method) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"), text = element_text(size = 10)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
