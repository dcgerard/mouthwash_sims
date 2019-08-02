
# Source to reproduce results from Gerard and Stephens (2018)

[![DOI](https://zenodo.org/badge/71280495.svg)](https://zenodo.org/badge/latestdoi/71280495)

This repository contains source code to reproduce the empirical
evaluations of Gerard and Stephens (2018). The new methods can be found
in the [vicar](https://github.com/dcgerard/vicar) package.

If you find a bug, please create an
[issue](https://github.com/dcgerard/ruvb_sims/issues).

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our paper:

> Gerard, David, and Matthew Stephens. 2018. “Empirical Bayes Shrinkage
> and False Discovery Rate Estimation, Allowing for Unwanted Variation.”
> *Biostatistics*, kxy029.
> <https://doi.org/10.1093/biostatistics/kxy029>.

## License

Copyright (c) 2017-2018, David Gerard.

All source code and software in this repository are made available under
the terms of the [GNU General Public
License](http://www.gnu.org/licenses/gpl.html). See the
[LICENSE](LICENSE) file for the full text of the license.

## Instructions

To reproduce the results of Gerard and Stephens (2018), you need to (i)
install the appropriate R packages, (ii) obtain the appropriate data,
(iii) run `make` and (iv) get some coffee while you wait. Please read
below for details on each of these steps.

### Install software and R packages

1.  Install [R](https://cran.r-project.org).
2.  Install [GNU Make](https://www.gnu.org/software/make).
3.  Install the required R packages by running the following commands in
    the R interactive environment. (Note the order of these commands is
    important—the Bioconductor packages should be installed before the
    CRAN packages.)

<!-- end list -->

``` r
install.packages(c("tidyverse", "stringr", "reshape2", "pROC",
                   "ruv", "cate", "devtools", "ashr", "bfa",
                   "xtable", "dplyr", "ggthemes",
                   "assertthat", "R.utils", "BiocManager"))
BiocManager::install(c("sva", "limma", "qvalue", "biomaRt"))
devtools::install_github("dcgerard/vicar", 
                         ref = "41fc3df8faed8d59af1f29410619908385970124")
devtools::install_github("dcgerard/seqgendiff", 
                         ref = "680e088c1a37879b82b7db86b24fbd657b5f6869")
```

Newer versions of `seqgendiff` should not work as the `poisthin()`
function has been replaced by `select_counts()` and `thin_2group()`.

See below for the versions of the packages that I used for the paper.

### Get data

Download the following files from the [GTEx
Portal](https://www.gtexportal.org) and copy them to the [Data](data)
subdirectory. In order to replicate my results, it is important to use
version “6p” of the GTEx data. To access these files, you will need to
register for a GTEx Portal account if you have not done so already.

1.  `GTEx_Data_V6_Annotations_SampleAttributesDS.txt`

2.  `GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz`

3.  `gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`

4.  `GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt`
    
    Next, download the list of human housekeeping genes from Eisenberg
    and Levanon (2013) and the NCBI-to-Ensembl gene mapping file, and
    copy these files to the [Data](data) directory:

5.  <http://www.tau.ac.il/~elieis/HKG/HK_genes.txt>

6.  <ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz>
    
    Finally, download the list of human housekeeping genes from Lin et
    al. (2017) by navigating to their [Shiny R
    application](http://shiny.maths.usyd.edu.au/scHK/), clicking on
    “Default Values”, then immediately clicking “Download”. Then place
    this file, labeled “h-scHKgenes.csv”, in the [Data](Data) directory.
    
    (Note: This shiny app is no longer being supported, so I placed the
    housekeeping genes that I downloaded in the [Data](Data) directory.)

7.  `h-scHKgenes.csv`
    
    I used the file [get\_ensembl.R](./R/get_ensembl.R) to correspond
    HGNC gene names to their ensembl annotation, but I don’t run this in
    make because I don’t want a dependency on the biomaRt package, so
    I’ve included the resulting file
    [lin\_hk\_genes.csv](./Data/lin_hk_genes.csv) in the [Data](Data)
    directory.

After completing these steps, the Data folder should look like this:

``` bash
cd Data
ls -1
```

    ## gencode.v19.genes.v6p_model.patched_contigs.gtf.gz
    ## gene2ensembl.gz
    ## GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz
    ## GTEx_Data_V6_Annotations_SampleAttributesDS.txt
    ## GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt
    ## HK_genes.txt
    ## h-scHKgenes.csv
    ## lin_hk_genes.csv

### Run make

To reproduce all of the results in Gerard and Stephens (2018), move to
the root directory in your local copy of this repository, and run `make`
from the terminal. This will run all the steps in the data processing
and analysis pipeline. Note that you may need to adjust some of the
settings in the [Makefile](Makefile) to suit your computing environment.

Since some of the steps take a long time to run, you may not want to run
everything at once. If you want to reproduce the simulation results, run

``` bash
make sims
```

If you want to reproduce just the results of the GTEx analysis using the
control genes from Eisenberg and Levanon (2013), run

``` bash
make gtex_analysis
```

If you want to reproduce just the results of the GTEx analysis using the
control genes from Lin et al. (2017), run

``` bash
make gtex_analysis_lin
```

If you want to reproduce the computation time calculations, run

``` bash
make computation
```

If you want to reproduce the figure in the introduction, run

``` bash
make one_data
```

At any point during or after the commands are running, it may be helpful
to inspect the `.Rout` files saved in the [Output](Output) subdirectory.

### Get Coffee

All of these runs (except the last one) should take a very long time (a
day to a couple of days). You should get some coffee. Here is a list of
some of my favorite places:

  - Chicago
      - [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
      - [Plein Air
        Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
  - Seattle
      - [Bauhaus
        Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
      - [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
  - Columbus
      - [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
      - [Stauf’s Coffee
        Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
      - [Caffe
        Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

If you are having trouble reproducing these results, it might be that
you need to update some of your R packages. These are the versions that
I used (including some versions of packages that are not actually needed
to run the code):

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
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
    ##  [1] vicar_0.1-9         seqgendiff_1.0.0    biomaRt_2.40.3     
    ##  [4] qvalue_2.16.0       limma_3.40.6        sva_3.32.1         
    ##  [7] BiocParallel_1.18.0 genefilter_1.66.0   mgcv_1.8-28        
    ## [10] nlme_3.1-141        R.utils_2.9.0       R.oo_1.22.0        
    ## [13] R.methodsS3_1.7.1   assertthat_0.2.1    ggthemes_4.2.0     
    ## [16] xtable_1.8-4        bfa_0.4             ashr_2.2-38        
    ## [19] devtools_2.1.0      usethis_1.5.1       cate_1.0.4         
    ## [22] ruv_0.9.6           pROC_1.15.3         reshape2_1.4.3     
    ## [25] forcats_0.4.0       stringr_1.4.0       dplyr_0.8.3        
    ## [28] purrr_0.3.2         readr_1.3.1         tidyr_0.8.3        
    ## [31] tibble_2.1.3        ggplot2_3.2.0       tidyverse_1.2.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] colorspace_1.4-1     rprojroot_1.3-2      corpcor_1.6.9       
    ##  [4] fs_1.3.1             rstudioapi_0.10      remotes_2.1.0       
    ##  [7] leapp_1.2            bit64_0.9-7          AnnotationDbi_1.46.0
    ## [10] lubridate_1.7.4.9000 xml2_1.2.1           codetools_0.2-16    
    ## [13] splines_3.6.1        pscl_1.5.2           doParallel_1.0.15   
    ## [16] knitr_1.23           pkgload_1.0.2        zeallot_0.1.0       
    ## [19] jsonlite_1.6         broom_0.5.2          annotate_1.62.0     
    ## [22] compiler_3.6.1       httr_1.4.0           backports_1.1.4     
    ## [25] Matrix_1.2-17        lazyeval_0.2.2       cli_1.1.0           
    ## [28] htmltools_0.3.6      prettyunits_1.0.2    tools_3.6.1         
    ## [31] coda_0.19-3          gtable_0.3.0         glue_1.3.1          
    ## [34] Rcpp_1.0.2           Biobase_2.44.0       cellranger_1.1.0    
    ## [37] vctrs_0.2.0          iterators_1.0.12     xfun_0.8            
    ## [40] ps_1.3.0             testthat_2.2.1       rvest_0.3.4         
    ## [43] XML_3.98-1.20        MASS_7.3-51.4        scales_1.0.0        
    ## [46] hms_0.5.0            parallel_3.6.1       yaml_2.2.0          
    ## [49] memoise_1.1.0        stringi_1.4.3        RSQLite_2.1.2       
    ## [52] SQUAREM_2017.10-1    S4Vectors_0.22.0     desc_1.2.0          
    ## [55] foreach_1.4.7        BiocGenerics_0.30.0  pkgbuild_1.0.3      
    ## [58] truncnorm_1.0-8      rlang_0.4.0          pkgconfig_2.0.2     
    ## [61] matrixStats_0.54.0   bitops_1.0-6         evaluate_0.14       
    ## [64] lattice_0.20-38      esaBcv_1.2.1         bit_1.1-14          
    ## [67] tidyselect_0.2.5     processx_3.4.1       plyr_1.8.4          
    ## [70] magrittr_1.5         R6_2.4.0             IRanges_2.18.1      
    ## [73] generics_0.0.2       DBI_1.0.0            pillar_1.4.2        
    ## [76] haven_2.1.1          withr_2.1.2          survival_2.44-1.1   
    ## [79] RCurl_1.95-4.12      mixsqp_0.1-97        modelr_0.1.4        
    ## [82] crayon_1.3.4         rmarkdown_1.14       progress_1.2.2      
    ## [85] grid_3.6.1           readxl_1.3.1         blob_1.2.0          
    ## [88] callr_3.3.1          digest_0.6.20        svd_0.4.3           
    ## [91] stats4_3.6.1         munsell_0.5.0        sessioninfo_1.1.1

I’ve also only tried this out on Ubuntu.

## Credits

This project was developed by [David Gerard](https://dcgerard.github.io)
at the University of Chicago.

[Peter Carbonetto](https://pcarbo.github.io/) made many fantastic
contributions and suggestions to increase the reproducibility of this
work.

Thanks to [Matthew Stephens](http://stephenslab.uchicago.edu) for his
support and mentorship.

## References

<div id="refs" class="references">

<div id="ref-eisenberg2013human">

Eisenberg, Eli, and Erez Y Levanon. 2013. “Human Housekeeping Genes,
Revisited.” *Trends in Genetics* 29 (10). Elsevier: 569–74.
<https://doi.org/10.1016/j.tig.2013.05.010>.

</div>

<div id="ref-gerard2018empirical">

Gerard, David, and Matthew Stephens. 2018. “Empirical Bayes Shrinkage
and False Discovery Rate Estimation, Allowing for Unwanted Variation.”
*Biostatistics*, kxy029. <https://doi.org/10.1093/biostatistics/kxy029>.

</div>

<div id="ref-lin2017housekeeping">

Lin, Yingxin, Shila Ghazanfar, Dario Strbenac, Andy Wang, Ellis Patrick,
Terence Speed, Jean Yang, and Pengyi Yang. 2017. “Housekeeping Genes,
Revisited at the Single-Cell Level.” *bioRxiv*. Cold Spring Harbor
Laboratory. <https://doi.org/10.1101/229815>.

</div>

</div>
