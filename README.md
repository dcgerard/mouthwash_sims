Reproducing MOUTHWASH results in Gerard and Stephens (2016)
================

Introduction
============

This repository contains code to reproduce the empirical evaluations of Gerard and Stephens (2016). The new methods can be found in the [vicar](https://github.com/dcgerard/vicar) package.

If you are having trouble reproducing these results, it might be that you need to update some of your R packages. These are the versions that I used:

``` r
sessionInfo()
```

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 14.04.5 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] Rmpi_0.6-6        vicar_0.1.4       limma_3.26.3     
    ##  [4] sva_3.18.0        genefilter_1.52.0 mgcv_1.8-15      
    ##  [7] nlme_3.1-128      devtools_1.12.0   snow_0.4-1       
    ## [10] gridExtra_2.2.1   cate_1.0.4        ruv_0.9.6        
    ## [13] pROC_1.8          reshape2_1.4.1    ggplot2_2.1.0    
    ## [16] stringr_1.1.0     dplyr_0.5.0      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] splines_3.3.1        lattice_0.20-34      colorspace_1.2-6    
    ##  [4] htmltools_0.3.5      stats4_3.3.1         yaml_2.1.13         
    ##  [7] XML_3.98-1.4         survival_2.39-5      withr_1.0.2         
    ## [10] DBI_0.5-1            BiocGenerics_0.16.1  plyr_1.8.4          
    ## [13] munsell_0.4.3        leapp_1.2            gtable_0.2.0        
    ## [16] svd_0.4              evaluate_0.9         memoise_1.0.0       
    ## [19] Biobase_2.30.0       knitr_1.14           IRanges_2.4.6       
    ## [22] parallel_3.3.1       AnnotationDbi_1.32.3 esaBcv_1.2.1        
    ## [25] Rcpp_0.12.7          xtable_1.8-2         corpcor_1.6.8       
    ## [28] scales_0.4.0         formatR_1.4          S4Vectors_0.8.6     
    ## [31] annotate_1.48.0      digest_0.6.10        stringi_1.1.1       
    ## [34] grid_3.3.1           tools_3.3.1          magrittr_1.5        
    ## [37] tibble_1.2           RSQLite_1.0.0        MASS_7.3-45         
    ## [40] Matrix_1.2-7.1       assertthat_0.1       rmarkdown_1.0       
    ## [43] R6_2.1.3             compiler_3.3.1

As you can see above, I've also only tried this out on Ubuntu.

If you find a bug, please create an [issue](https://github.com/dcgerard/ruvb_sims/issues).

Instructions
============

To reproduce the results of Gerard and Stephens (2016), you need to (1) install the appropriate R packages, (2) obtain the appropriate data, (3) run `make` and (4) get some coffee while you wait.

Install R Packages
------------------

To install the needed R packages, run the following in R

``` r
install.packages(c("dplyr", "stringr", "ggplot2", "reshape2", "pROC",
                   "ruv", "cate", "gridExtra", "snow", "devtools", "Rmpi"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("sva", "limma"))
devtools::install_github("dcgerard/vicar")
```

Get Data
--------

Place the following files in the Data folder:

1.  [GTEx\_Data\_V6\_Annotations\_SampleAttributesDS.txt](http://www.gtexportal.org/home/datasets#filesetFilesDiv21)
2.  [GTEx\_Analysis\_v6p\_RNA-seq\_RNA-SeQCv1.1.8\_gene\_reads.gct.gz](http://www.gtexportal.org/home/datasets#filesetFilesDiv11)
3.  [gencode.v19.genes.V6p\_model.patched\_contigs.gtf](http://www.gtexportal.org/home/datasets#filesetFilesDiv14)
4.  [GTEx\_Data\_V6\_Annotations\_SubjectPhenotypesDS.txt](http://www.gtexportal.org/home/datasets#datasetDiv2)
5.  [HK\_genes.txt](http://www.tau.ac.il/~elieis/HKG/HK_genes.txt)
6.  gene2ensembl.gz at <ftp://ftp.ncbi.nih.gov/gene/DATA/>

1 through 4 of the above are only available if you are a registered user of the GTEx Portal. I don't think I'm allowed to release these data.

Run Make
--------

To reproduce all of the results in Gerard and Stephens (2016), simply run `make` from the terminal.

If you want to reproduce just the results from Section 5.1, run

``` bash
make sims
```

If you want to reproduce just the results from Section 5.2, run

``` bash
make gtex_analysis
```

If you want to reproduce the figure in the introduction, run

``` bash
make one_data
```

Get Coffee
----------

All of these runs (except the last one) should take a very long time (a day to a couple of days). You should get some coffee. If you're in the Chicago area, I would recommend [Plein Air Cafe](http://www.pleinaircafe.co/).

References
==========

Gerard, David, and Matthew Stephens. 2016. “MOUTHWASH: Maximizing Over Unobservables To Help With Adaptive SHrinkage.” *Github Preprint*. <https://github.com/stephenslab/RUV-B>.
