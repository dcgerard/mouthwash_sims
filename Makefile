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

# Results of fitting the various confounder adjustment methods to all GTEx tissues.
# This uses the eisenburg control genes.
gtex_dir  = ./Output/gtex_fits
gtex_fits = $(gtex_dir)/betahat_list.Rds \
            $(gtex_dir)/pi0mat.Rds \
            $(gtex_dir)/plist.Rds

gtex_lin_dir  = ./Output/gtex_fits_lin
gtex_lin_fits = $(gtex_lin_dir)/betahat_list.Rds \
                $(gtex_lin_dir)/pi0mat.Rds \
                $(gtex_lin_dir)/plist.Rds

# Figures output from ash_problems.R
ash_prob_out_dir = ./Output/figures
ash_prob_out = $(ash_prob_out_dir)/ash_fail.eps \
               $(ash_prob_out_dir)/ash_fail_small.eps \
               $(ash_prob_out_dir)/ash_fail_small_mouth.eps

# Output from summarize_computation.R
sum_comp_out_dir = ./Output/computation
sum_comp_out = $(sum_comp_out_dir)/comp_tab.txt

# Output from the main simulation study
sims_out = ./Output/sims_out/sims_out.Rds

# Output from the simulation study using just control-gene methods
sims_out_control = ./Output/sims_out_control/sims_out_control.Rds

# Figures output from plot_mouthwash_sims.R
sim_fig_dir = ./Output/figures
sim_fig = $(sim_fig_dir)/pi0_box_50.eps \
          $(sim_fig_dir)/pi0_box_90.eps \
          $(sim_fig_dir)/pi0_box_100.eps \
          $(sim_fig_dir)/auc_ave.eps

# Figures output from gtex_plots.R
gtex_fig_dir = ./Output/figures
gtex_fig = $(gtex_fig_dir)/prop_max.eps \
           $(gtex_fig_dir)/proponsex.eps \
           $(gtex_fig_dir)/proponsex_bw.eps \
           $(gtex_fig_dir)/lfdr_rank.eps

# Figures output from gtex_plots_lin.R
gtex_fig_dir_lin = ./Output/figures
gtex_fig_lin = $(gtex_fig_dir_lin)/prop_max_lin.eps \
               $(gtex_fig_dir_lin)/proponsex_lin.eps \
               $(gtex_fig_dir_lin)/lfdr_rank_lin.eps

# Data frame of computation time for all methods
comp_time = ./Output/computation/comp_time.Rds

all: sims gtex_analysis gtex_analysis_lin one_data computation

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

# Re-run GTEx analysis with control genes from Lin et al.
$(gtex_lin_fits) : ./R/fit_gtex_lin.R $(cleaned_dat)
	mkdir -p Output/gtex_fits_lin
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Create plots from results of GTEx analysis.
$(gtex_fig) : ./R/gtex_plots.R $(gtex_fits)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

.PHONY : gtex_analysis
gtex_analysis : $(gtex_fig)

# Create same plots from results of GTEx analysis where we use
# the control genes from Lin et al.
$(gtex_fig_lin) : ./R/gtex_plots_lin.R $(gtex_lin_fits)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

.PHONY : gtex_analysis_lin
gtex_analysis_lin : $(gtex_fig_lin)

# Run simulations.
$(sims_out) : ./Code/mouthwash_sims.R $(tissue_dat)
	mkdir -p Output/sims_out
	$(rexec) '--args nc=$(nc)' $< Output/$(basename $(notdir $<)).Rout

# Create plots from simulation experiments.
$(sim_fig) : ./R/plot_mouthwash_sims.R $(sims_out)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

.PHONY : sims
sims : $(sim_fig)

# Provide example figure of ash's poor performance in the presence of
# unobserved confounding.
$(ash_prob_out) : ./R/ash_problems.R $(tissue_dat)
	mkdir -p Output/figures
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

.PHONY : one_data
one_data : $(ash_prob_out)

# Compuatation time simulations
$(comp_time) : ./R/computation_time.R $(tissue_dat)
	mkdir -p Output/computation
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

# Summarize computation time results
$(sum_comp_out) : ./R/summarize_computation.R Output/computation/comp_time.Rds
	mkdir -p Output/computation
	$(rexec) $< Output/$(basename $(notdir $<)).Rout

.PHONY : computation
computation : $(sum_comp_out)

clean:
	rm -f $(tissue_dat)
