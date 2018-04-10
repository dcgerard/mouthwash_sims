# ADJUST THESE VARIABLES AS NEEDED TO SUIT YOUR COMPUTING ENVIRONMENT
# -------------------------------------------------------------------
# This variable specifies the number of threads to use for the
# parSapply calls. This could also be specified automatically using
# environment variables. For example, in SLURM, SLURM_CPUS_PER_TASK
# specifies the number of CPUs allocated for each task.
nc = 6

# R scripting front-end. Note that makeCluster sometimes fails to
# connect to a socker when using Rscript, so we are using the "R CMD
# BATCH" interface instead.
rexec = R CMD BATCH --no-save --no-restore

# AVOID EDITING ANYTHING BELOW THIS LINE
# --------------------------------------
# These variables specify the CSV files created by gtex_extract_tissues_v6p.R.
tissue_dir = ./Output/gtex_tissue_gene_reads_v6p
tissue_dat = $(tissue_dir)/adiposetissue.csv \
	     $(tissue_dir)/bladder.csv \
             $(tissue_dir)/bloodvessel.csv \
             $(tissue_dir)/breast.csv \
             $(tissue_dir)/colon.csv \
             $(tissue_dir)/fallopiantube.csv \
             $(tissue_dir)/kidney.csv \
             $(tissue_dir)/lung.csv \
             $(tissue_dir)/nerve.csv \
             $(tissue_dir)/pancreas.csv \
             $(tissue_dir)/prostate.csv \
             $(tissue_dir)/skin.csv \
             $(tissue_dir)/spleen.csv \
             $(tissue_dir)/testis.csv \
             $(tissue_dir)/uterus.csv \
             $(tissue_dir)/adrenalgland.csv \
             $(tissue_dir)/blood.csv \
             $(tissue_dir)/brain.csv \
             $(tissue_dir)/cervixuteri.csv \
             $(tissue_dir)/esophagus.csv \
             $(tissue_dir)/heart.csv \
             $(tissue_dir)/liver.csv \
             $(tissue_dir)/muscle.csv \
             $(tissue_dir)/ovary.csv \
             $(tissue_dir)/pituitary.csv \
             $(tissue_dir)/salivarygland.csv \
             $(tissue_dir)/smallintestine.csv \
             $(tissue_dir)/stomach.csv \
             $(tissue_dir)/thyroid.csv \
             $(tissue_dir)/vagina.csv
tissue_targets = $(addsuffix .%,$(basename $(tissue_dat)))

# Contains data in the format used as input for MOUTHWASH
cleaned_dir = ./Output/cleaned_gtex_data
cleaned_dat = $(cleaned_dir)/adiposetissue.Rds \
	      $(cleaned_dir)/bladder.Rds \
	      $(cleaned_dir)/bloodvessel.Rds \
	      $(cleaned_dir)/breast.Rds \
	      $(cleaned_dir)/esophagus.Rds \
	      $(cleaned_dir)/kidney.Rds \
	      $(cleaned_dir)/lung.Rds \
	      $(cleaned_dir)/nerve.Rds \
	      $(cleaned_dir)/pituitary.Rds \
	      $(cleaned_dir)/skin.Rds \
	      $(cleaned_dir)/spleen.Rds \
	      $(cleaned_dir)/thyroid.Rds \
	      $(cleaned_dir)/adrenalgland.Rds \
	      $(cleaned_dir)/blood.Rds \
	      $(cleaned_dir)/brain.Rds \
	      $(cleaned_dir)/colon.Rds \
	      $(cleaned_dir)/heart.Rds \
	      $(cleaned_dir)/liver.Rds \
	      $(cleaned_dir)/muscle.Rds \
	      $(cleaned_dir)/pancreas.Rds \
	      $(cleaned_dir)/salivarygland.Rds \
	      $(cleaned_dir)/smallintestine.Rds \
	      $(cleaned_dir)/stomach.Rds

# TO DO: Add description of these targets here. See above description
# of tissue_dat for an example.
gtex_dir  = ./Output/gtex_fits
gtex_fits = $(gtex_dir)/betahat_list.Rds \
            $(gtex_dir)/pi0mat.Rds \
            $(gtex_dir)/plist.Rds

# TO DO: Add description of these targets here. See above description
# of tissue_dat for an example.
sims_out = ./Output/sims_out/sims_out.Rds

all: sims gtex_analysis one_data

# Extract tissue data.
$(tissue_targets) : ./R/gtex_extract_tissues_v6p.R
	mkdir -p Output/gtex_tissue_gene_reads_v6p
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Clean tissue data for GTEx analysis.
$(cleaned_dat) : ./R/gtex_clean.R $(tissue_dat)
	mkdir -p Output/cleaned_gtex_data
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Get Lin et al control genes
$(cleaned_dir)/lin_ctl.Rds : ./R/gtex_add_lin_control.R ./Data/lin_hk_genes.csv $(cleaned_dat)
	mkdir -p Output/cleaned_gtex_data
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Run GTEx analysis.
$(gtex_fits) : ./R/fit_mouthwash.R $(cleaned_dat)
	mkdir -p Output/gtex_fits
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Create plots from results of GTEx analysis.
.PHONY : gtex_analysis
gtex_analysis : ./R/gtex_plots.R $(gtex_fits)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Run simulations.
$(sims_out) : ./Code/mouthwash_sims.R $(tissue_dat)
	mkdir -p Output/sims_out
	$(rexec) '--args nc=$(nc)' $< Output/$(basename $(notdir $<)).Rout

# Create plots from simulation experiments.
.PHONY : sims
sims : ./R/plot_mouthwash_sims.R $(sims_out)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Provide example figure of ash's poor performance in the presence of
# unobserved confounding.
.PHONY : one_data
one_data : ./R/ash_problems.R $(tissue_dat)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

clean:
	rm -f $(tissue_dat)
