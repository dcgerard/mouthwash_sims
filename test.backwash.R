# Small demo of MOUTHWASH and BACKWASH on GTEx data set.
library(sva)
library(vicar)

cat("Loading data.\n")
dat   <- readRDS("Output/cleaned_gtex_data/muscle.Rds")
onsex <- dat$chrom == "X" | dat$chrom == "Y"
onsex[is.na(onsex)] <- FALSE
dat$ctl[onsex]      <- FALSE

cat("Running SVA.\n")
num_sv <- sva::num.sv(dat = dat$Y,mod = dat$X)
Y      <- t(dat$Y)
X      <- dat$X

cat("Running mouthwash.\n")
r <- system.time(mout <- mouthwash(Y = Y,X = X,k = num_sv))
cat("mouthwash took",r["elapsed"],"seconds.\n")

cat("Running backwash.\n")
r <- system.time(bout <- backwash(Y = Y,X = X,k = num_sv))
cat("Computation took",r["elapsed"],"seconds.\n")

