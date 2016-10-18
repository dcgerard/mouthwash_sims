library(dplyr)
library(stringr)


tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")


## This is from gencode here http://useast.ensembl.org/info/data/ftp/index.html
## The GTF is what you want
## exon2 <- read.delim(file = "/home/david/Data/gtex/gene_reference/Homo_sapiens.GRCh38.85.gtf", skip = 5)



gtout <- read.table(file = "./Data/gencode.v19.genes.V6p_model.patched_contigs.gtf",
                    skip = 6, sep = "\t")
temp1 <- stringr::str_split(gtout[, 9], pattern = ";")

geneout <- stringr::str_replace(sapply(X = temp1, FUN = stringr::str_subset, pattern = "gene_id"),
                            pattern = "gene_id ", replacement = "")
transout <- stringr::str_replace(sapply(X = temp1, FUN = stringr::str_subset,
                                        pattern = "transcript_id"),
                                 pattern = " transcript_id ", replacement = "")

exon <- cbind(gtout[, 1:8], gene_id = geneout, transcript_id = transout)
colnames(exon)[1:8] <- c("CHR", "annotation_source", "feature_type", "gene_start",
                         "gene_end", "score", "strand", "phase")

## this file is from the GTEX portal
write.csv(exon, file = "./Output/gencode.v19.genes.V6p_model.patched_contigs.txt")

## Associate genes with chromosome position
gene_names <- rep(NA, length = length(exon$gene_id))
for (index in 1:length(gene_names)) {
    gene_names[index] <- sub("([[:alnum:]]+\\.[[:digit:]]+)\\_[[:digit:]]+", "\\1",
                             exon$gene_id[index])
}

for (tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    ## Read in Data
    tissuedat <- read.csv(paste0("./Output/gtex_tissue_gene_reads_v6p/", current_tissue,
                                 ".csv"))



    match_out <- match(tissuedat$Name, gene_names)
    tissuedat$CHR <- exon$CHR[match_out]

    ## ----------------------------------------------------------------------------

    ## Select all genes that have on average at least 10 reads per gene.
    expr <- select(tissuedat, -c(Name, Description, CHR))
    gene_levels <- rowSums(expr)
    ngene <- sum(gene_levels > 10 * ncol(expr))
    genes_order <- order(gene_levels, decreasing = TRUE)
    subtissue <- tissuedat[genes_order[1:ngene], ]

    rm(tissuedat)
    rm(expr)
    ## ----------------------------------------------------------------------------
    ngene ## Number of genes selected.

    ## Average technical replicates
    colname_vec <- rep(NA, ncol(subtissue))
    for (index in 3:(ncol(subtissue) - 1)) {
        colname_vec[index] <- sub(pattern = "^([[:alnum:]]+)\\.([[:alnum:]]+)\\..+",
                                  replacement = "\\1-\\2", x = names(subtissue)[index])
    }
    unique_subjects <- unique(colname_vec)[-1] ## to remove NA

    Y <- matrix(NA, ncol = length(unique_subjects), nrow = ngene)
    colnames(Y) <- unique_subjects
    rownames(Y) <- subtissue$Name

    for (index in 1:length(unique_subjects)) {
        current_cols <- colname_vec == unique_subjects[index]
        current_cols[is.na(current_cols)] <- FALSE
        Y[, index] <- rowMeans(log2(subtissue[, current_cols, drop = FALSE] + 1))
    }

    ## Now get the X matrix
    pheno <- read.table("./Data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",
                        header = TRUE)
    match_out <- match(colnames(Y), pheno$SUBJID)
    X <- model.matrix(~as.factor(pheno$GENDER[match_out]))
    colnames(X) <- c("Intercept", "Gender")

    ## Find housekeeping genes. Genes are from E. Eisenberg and E.Y. Levanon,
    ## Trends in Genetics, 29 (2013) and can be found at
    ## [http://www.tau.ac.il/~elieis/HKG/](http://www.tau.ac.il/~elieis/HKG/).

    ## Associating NCBI reference sequences with Ensembl annotations is done
    ## with the file gene2ensemble.gz which can be found at
    ## [ftp://ftp.ncbi.nih.gov/gene/DATA/](ftp://ftp.ncbi.nih.gov/gene/DATA/).

    hkgenes <- read.table("./Data/HK_genes.txt")
    gene2ensembl <- read.table("./Data/gene2ensembl.gz")
    ## V4 is NCBI, V3 is Ensemble

    ctl <- rep(FALSE, length = nrow(Y))
    for (index in 1:nrow(hkgenes)) {
        pattern_current <- paste0(hkgenes[index, 2], "\\.[[:digit:]]+")
        pos_pat <- grepl(pattern_current, gene2ensembl[, 4])
        if (sum(pos_pat) > 1) {
            stop("not correct match in gene2ensembl")
        } else if (sum(pos_pat) == 1) {

            index_pat <- which(pos_pat)

            enspat <- paste0(gene2ensembl[index_pat, 3], "\\.[[:digit:]]+")
            pos_ens <- grepl(enspat, rownames(Y))
            if (sum(pos_ens) == 1) {
                index_ens <- which(pos_ens)
                ctl[index_ens] <- TRUE
            } else if (sum(pos_ens) > 1) {
                stop ("to many matches in Y")
            }
        }
    }

    ## Save file
    final_tiss <- list(Y = Y, X = X, ctl = ctl, chrom = subtissue$CHR)
    saveRDS(object = final_tiss,
            file = paste0("./Output/cleaned_gtex_data/", current_tissue, ".Rds"))
}
