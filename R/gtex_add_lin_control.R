## Add new control genes to cleaned_gtex_data
suppressMessages(library(tidyverse))
tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")

hkvec <- read.csv("./Data/lin_hk_genes.csv")[, 1]
lin_ctl_list <- list()
for (tissue_index in 1:length(tissue_vec)) {
  current_tissue <- tissue_vec[tissue_index]
  cat("tissue",tissue_index,"=",current_tissue,"\n")

  final_tiss <- readRDS(file = paste0("./Output/cleaned_gtex_data/",
                                      current_tissue, ".Rds"))
  enz_vec <- stringr::str_replace_all(row.names(final_tiss$Y), "\\.\\d+$", "")
  lin_ctl_list[[tissue_index]] <- lin_ctl <- !is.na(match(enz_vec, hkvec))
}
names(lin_ctl_list) <- tissue_vec
saveRDS(object = lin_ctl_list, "./Output/cleaned_gtex_data/lin_ctl.Rds")
