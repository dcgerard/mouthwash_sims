library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)

do_topk <- function(lfdr, onsex, k) {
    return(sum(onsex[order(lfdr)][1:k]))
}

tot_sig <- function(lfdr, onsex, alpha = 0.1) {
    sig <- onsex[lfdr < alpha]
    return(c(numsex = sum(sig), sumsig = length(sig)))
}

prop_sig <- function(lfdr, onsex, alpha = 0.1) {
    sig <- onsex[lfdr < alpha]
    propsig <- length(sig) / length(lfdr)
    propgen <- sum(sig) / length(sig)
    return(c(propsig = propsig, propgen = propgen))
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

#############################################################################
## CATErr+Cal and CATErr were misslabeld in simulations. When you re-run, switch them back!
#############################################################################


good_method_labels <- c("OLS+ASH", "RUV2+ASH", "RUV3+ASH", "RUV4+ASH",
                        "RUV4v+ASH", "CATEnc+ASH",
                        "CATEv(normal)+ASH", "CATEv(t)+ASH",
                        "SVA+ASH", "CATErr+ASH", "MOUTHWASH(normal)",
                        "MOUTHWASH(t)", "OLS", "RUV2", "RUV3", "RUV4",
                        "RUV4v", "CATEnc", "CATEnc+Cal",
                        "CATEv(normal)", "CATEv(t)", "SVA", "CATErr",
                        "CATErr+Cal")


top100 <- matrix(NA, nrow = length(tissue_vec), ncol = 24)

plist <- readRDS("./Output/gtex_fits/plist.Rds")
pi0mat <- readRDS("./Output/gtex_fits/pi0mat.Rds")
betahat_list <- readRDS("./Output/gtex_fits/betahat_list.Rds")
nseq <- rep(NA, length(tissue_vec))

for (tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    dat <- readRDS(paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))

    onsex <- dat$chrom == "X" | dat$chrom == "Y"

    pmat <- plist[[tissue_index]]

    numzero <- apply(pmat, 2, function(x) sum(x == 0, na.rm = TRUE))
    if (any(numzero > 100)) {
        cat("LOOK HERE:", tissue_index, "\n")
    }

    betamat <- betahat_list[[tissue_index]]

    top100[tissue_index, ] <- apply(pmat, 2, do_topk, onsex = onsex, k = 100)
    nseq[tissue_index] <- ncol(dat$Y)

    countmat <- as.data.frame(t(apply(pmat[, 1:12], 2, tot_sig, onsex = onsex, alpha = 0.05)))
    propmat  <- as.data.frame(t(apply(pmat[, 1:12], 2, prop_sig, onsex = onsex, alpha = 0.05)))

    namevec <- stringr::str_replace(rownames(countmat), "ash_", "")

    countmat$method <- namevec
    propmat$method  <- namevec

    countmat$tissue <- tissue_vec[tissue_index]
    propmat$tissue  <- tissue_vec[tissue_index]

    if (tissue_index == 1) {
        allcountmat <- countmat
        allpropmat  <- propmat
    } else {
        allcountmat <- rbind(allcountmat, countmat)
        allpropmat  <- rbind(allpropmat, propmat)
    }
    cat(tissue_index, "\n")
}


## top 100 analysis
colnames(top100) <- good_method_labels
rownames(top100) <- paste0(good_tissue_labels, " (", nseq, ")")
top100p <- as.data.frame(top100 / apply(top100, 1, max))
top100p <- dplyr::select(top100p, -OLS, -`OLS+ASH`)
faorder <- order(apply(top100p, 2, median), decreasing = TRUE)
top100p$Tissue <- rownames(top100p)




longdat <- melt(top100p, id.vars = "Tissue")
names(longdat) <- c("Tissue", "Method", "Proportion")
which_bad <- longdat$Method == "CATEnc" & longdat$Tissue == "Breast (214)"
longdat$Proportion[which_bad] <- NA
longdat$na <- FALSE
longdat$na[which_bad] <- TRUE

longdat$Tissue <- factor(longdat$Tissue, levels = rownames(top100)[order(nseq, decreasing = TRUE)])

longdat$Method <- factor(longdat$Method, levels = levels(longdat$Method)[faorder])


pdf(file = "./Output/figures/prop_max.pdf", height = 7.5, width = 6.5, family = "Times", color = "cmyk")
ggplot(data = longdat, mapping = aes(x = Method, y = Tissue, fill = Proportion)) +
    geom_raster() +
    scale_fill_gradient(high = "#ffffff",
                        low = "#000000",
                        na.value = "black",
                        guide = guide_colourbar(title = "Proportion\nfrom\nMaximum")) +
    theme_bw() +
    geom_point(aes(size=ifelse(na, "dot", "no_dot")), color = "white") +
    scale_size_manual(values=c(dot=2, no_dot=NA), guide="none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()

medpi0 <- data.frame(pi0hat = apply(pi0mat, 2, median))
rownames(medpi0) <- NULL
medpi0$Method <- good_method_labels
medpi0 <- medpi0[, 2:1]
medpi0 <- medpi0[order(medpi0$pi0hat), ]

medpi0$Method[!str_detect(medpi0$Method, "ASH")] <-
    paste0(medpi0$Method[!str_detect(medpi0$Method, "ASH")], "+qvalue")

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

plot_now <- FALSE
for (tissue_index in 1:length(tissue_vec)) {

    current_tissue <- tissue_vec[tissue_index]
    dat <- readRDS(paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))
    onsex <- dat$chrom == "X" | dat$chrom == "Y"

    pmat <- plist[[tissue_index]][, 1:12]

    aout <- lapply(as.data.frame(pmat), return_topk, onsex = onsex)

    for (index in 1:length(aout)) {
        aout[[index]] <- as.data.frame(aout[[index]])
        aout[[index]]$method <- names(aout)[index]
        aout[[index]]$index  <- 1:nrow(aout[[index]])
    }

    longdat <- do.call(rbind, aout)
    longdat$tissue <- current_tissue

    if (tissue_index == 1) {
        biglongdat <- longdat
    } else {
        biglongdat <- rbind(biglongdat, longdat)
    }

    names(longdat) <- c("lfdr", "onsex", "method", "index")

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

## ggplot(data = biglongdat, mapping = aes(y = lfdr, x = index)) +
##     facet_wrap(~method) +
##     geom_point(size = 0.2, mapping = aes(color = as.factor(onsex))) +
##     geom_line(alpha = 0.2, mapping = aes(group = tissue))


change_method_names <- function(method_names) {
    method_names <- str_replace(method_names, "ash_", "")
    method_names <- str_replace(method_names, "_nocal", "")
    method_names <- str_replace(method_names, "ruv4v", "CATEv")
    method_names <- str_replace(method_names, "_nomult", "")
    method_names <- str_replace(method_names, "_rsvar", "v")
    method_names <- str_to_upper(method_names)
    method_names <- str_replace(method_names, "_NORM", "(normal)")
    method_names <- str_replace(method_names, "_T", "(t)")
    method_names <- str_replace(method_names, "CATENC", "CATEnc")
    method_names <- str_replace(method_names, "CATERR", "CATErr")
    method_names <- str_replace(method_names, "RUV4V", "RUV4v")
    method_names <- str_replace(method_names, "CATEV", "CATEv")
    return(method_names)
}

dummydat$method <- change_method_names(dummydat$method)

names(dummydat) <- c("Method", "Rank", "lfdr", "Proportion")

biglongdat2 <- biglongdat
names(biglongdat2) <- c("lfdr", "onsex", "Method", "Rank", "Tissue")

biglongdat2$Method <- change_method_names(biglongdat2$Method)

pdf(file = "./Output/figures/proponsex.pdf", height = 5.5, width = 6.5, family = "Times", color = "cmyk")
ggplot() + geom_point(data = dummydat, mapping = aes(x = Rank, y = lfdr,
                                                    color = Proportion), size = 0.5) +
    facet_wrap(~Method) +
    theme_bw() +
    scale_color_gradient(high = "black",
                        low = "grey90",
                        guide = guide_colourbar(title = "Proportion\non Sex\nChromosome")) +
    theme(strip.background = element_rect(fill="white"), text = element_text(size = 10)) +
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
