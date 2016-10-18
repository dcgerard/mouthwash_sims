################
## Extract the the tissues and output to new csv files
################

sample_attributes <- "./Data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
gene_expression   <- "./Data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"

samp_df     <- read.csv(sample_attributes, sep = "\t")
gene_exp_df <- read.csv(gene_expression, sep = "\t", skip = 2)

samp_names <- samp_df$SAMPID
gene_exp_names <- sapply(X = colnames(gene_exp_df), FUN = gsub, pattern = "\\.", replacement = "-")

match_id <- match(gene_exp_names, samp_names)
samp_df_less <- samp_df[match_id, ]

samp_names <- samp_df_less$SAMPID

if(!all(samp_names == gene_exp_names, na.rm = TRUE)) {
    stop("no match")
}

tissues <- samp_df_less$SMTS
unique_tissues <- unique(tissues)
unique_tissues <- unique_tissues[!is.na(unique_tissues)]

output_filenames <- paste0(sapply(X = sapply(unique_tissues, tolower),
                                  FUN = gsub, pattern = " ", replacement = ""),
                           ".csv")

for (tiss_index in 1:length(unique_tissues)) {
    tiss_indicator <- tissues == unique_tissues[tiss_index]
    tiss_indicator[c(1, 2)] <- TRUE
    subdf <- gene_exp_df[, tiss_indicator]
    write.csv(subdf, file = paste0("./Output/gtex_tissue_gene_reads_v6p/",
                                   output_filenames[tiss_index]), row.names = FALSE)
}


## Implementation Check
rm(list = ls())
sample_attributes <- "./Data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
samp_df     <- read.csv(sample_attributes, sep = "\t")

tiss <- read.csv("./Output/gtex_tissue_gene_reads_v6p/heart.csv")
tiss_names <- sapply(X = colnames(tiss), FUN = gsub, pattern = "\\.", replacement = "-")
match_id <- match(tiss_names, samp_df$SAMPID)
samp_df$SMTS[match_id]
