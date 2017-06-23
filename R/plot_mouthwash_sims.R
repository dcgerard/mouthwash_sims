#####################
## Plot MOUTHWASH results
#####################

replace_names <- function(x) {
  x <- stringr::str_replace(x, "pi0_", "")
  x <- stringr::str_replace(x, "auc_", "")
  x <- stringr::str_replace(x, "mse_", "")
  x <- stringr::str_replace(x, "ash_(.+)", "\\1+ASH")
  x <- stringr::str_replace(x, "qvalue_(.+)", "\\1+qvalue")
  x <- stringr::str_replace(x, "ruv", "RUV")
  x <- stringr::str_replace(x, "cate", "CATE")
  x <- stringr::str_replace(x, "pvalue_", "")
  x <- stringr::str_replace(x, "ols", "OLS")
  x <- stringr::str_replace(x, "mouthwash", "MOUTHWASH")
  x <- stringr::str_replace(x, "backwash", "BACKWASH")
  x <- stringr::str_replace(x, "sva", "SVA")
  x <- stringr::str_replace(x, "_norm", "")
}

## pi0 first ------------------------------------------------------------
library(tidyverse)
dat <- as_data_frame(readRDS(file = "./Output/sims_out/sims_out.RDS"))
longdat <- select(dat, nullpi, Nsamp, ncontrols, contains("pi0")) %>%
  gather(key = "Method", value = "pi0hat", pi0_ash_ols:pi0_qvalue_caterr)
longdat$Method <- replace_names(longdat$Method)

for (nullpi_current in c(0.5, 0.9, 1)) {
  sublongdat <- filter(longdat, nullpi == nullpi_current)

  mse_dat <- sublongdat %>% filter(Nsamp == 40, ncontrols == 10) %>%
    group_by(Method) %>%
    summarize(mse = sum((pi0hat - nullpi_current) ^ 2))
  factor_levels <- mse_dat$Method[order(mse_dat$mse, decreasing = FALSE)]
  sublongdat$Method <- factor(sublongdat$Method, levels = factor_levels)

  pl <- ggplot(data = sublongdat, mapping = aes(y = pi0hat, x = Method)) +
    geom_hline(yintercept = nullpi_current, lty = 2) +
    facet_grid(Nsamp ~ ncontrols) +
    geom_boxplot() +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Estimate of Pi0")

  pdf(file = paste0("./Output/figures/pi0_box_", nullpi_current * 100, ".pdf"), colormodel = "cmyk",
      family = "Times", width = 6.5, height = 7)
  print(pl)
  dev.off()
}



## Now AUC --------------------------------------------------------------
dat <- as_data_frame(readRDS(file = "./Output/sims_out/sims_out.RDS"))
longdat <- select(dat, nullpi, Nsamp, ncontrols, contains("auc")) %>%
  gather(key = "Method", value = "auc", auc_ash_ols:auc_pvalue_caterr) %>%
  filter(nullpi != 1)

meddat <- longdat %>% group_by(nullpi, Nsamp, ncontrols, Method) %>%
  summarize(mean = mean(auc)) %>%
  ungroup()

meddat$Method <- replace_names(meddat$Method)

## Get ordering
subsamp       <- filter(meddat, Nsamp == 40, nullpi == 0.9, ncontrols == 10)
falevels      <- subsamp$Method[order(subsamp$mean, decreasing = TRUE)]
meddat$Method <- factor(meddat$Method, levels = falevels)
meddat$Nsamp  <- as.factor(meddat$Nsamp)

maxdat <- meddat %>%
  group_by(nullpi, Nsamp, ncontrols) %>%
  summarize(max = max(mean)) %>%
  ungroup()

pl <- ggplot(data = meddat, mapping = aes(x = Method, y = mean, pch = as.factor(Nsamp))) +
  facet_grid(nullpi ~ ncontrols)  +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_point() +
  geom_hline(data = maxdat, mapping = aes(yintercept = max, lty = Nsamp), alpha = 1/2) +
  scale_linetype_discrete(name = "Max Mean at\nSample Size") +
  scale_shape_discrete(name = "Sample Size") +
  ylab("Mean AUC")

pdf(file = "./Output/figures/auc_ave.pdf", colormodel = "cmyk",
    family = "Times", width = 6.5, height = 7)
print(pl)
dev.off()



