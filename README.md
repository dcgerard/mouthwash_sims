
# Source to reproduce results from Gerard and Stephens (2020)

[![DOI](https://zenodo.org/badge/71280495.svg)](https://zenodo.org/badge/latestdoi/71280495)

This repository contains source code to reproduce the empirical
evaluations of Gerard and Stephens (2020). The new methods can be found
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

To reproduce the results of Gerard and Stephens (2020), you need to (i)
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

To reproduce all of the results in Gerard and Stephens (2020), move to
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

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
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
    ##  [1] vicar_0.1-9         seqgendiff_1.1.1    biomaRt_2.42.0     
    ##  [4] qvalue_2.18.0       limma_3.42.0        sva_3.34.0         
    ##  [7] BiocParallel_1.20.1 genefilter_1.68.0   mgcv_1.8-31        
    ## [10] nlme_3.1-143        R.utils_2.9.2       R.oo_1.23.0        
    ## [13] R.methodsS3_1.7.1   assertthat_0.2.1    ggthemes_4.2.0     
    ## [16] xtable_1.8-4        bfa_0.4             ashr_2.2-39        
    ## [19] devtools_2.2.1      usethis_1.5.1       cate_1.1           
    ## [22] ruv_0.9.7.1         pROC_1.16.1         reshape2_1.4.3     
    ## [25] forcats_0.4.0       stringr_1.4.0       dplyr_0.8.3        
    ## [28] purrr_0.3.3         readr_1.3.1         tidyr_1.0.0        
    ## [31] tibble_2.1.3        ggplot2_3.2.1       tidyverse_1.3.0    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_1.4-1     ellipsis_0.3.0       rprojroot_1.3-2     
    ##   [4] corpcor_1.6.9        fs_1.3.1             rstudioapi_0.10     
    ##   [7] remotes_2.1.0        leapp_1.2            bit64_0.9-7         
    ##  [10] AnnotationDbi_1.48.0 fansi_0.4.1          lubridate_1.7.4     
    ##  [13] xml2_1.2.2           codetools_0.2-16     splines_3.6.2       
    ##  [16] pscl_1.5.2           doParallel_1.0.15    knitr_1.27          
    ##  [19] pkgload_1.0.2        zeallot_0.1.0        jsonlite_1.6        
    ##  [22] broom_0.5.3          annotate_1.64.0      dbplyr_1.4.2        
    ##  [25] compiler_3.6.2       httr_1.4.1           backports_1.1.5     
    ##  [28] Matrix_1.2-18        lazyeval_0.2.2       cli_2.0.1           
    ##  [31] htmltools_0.4.0      prettyunits_1.1.0    tools_3.6.2         
    ##  [34] coda_0.19-3          gtable_0.3.0         glue_1.3.1          
    ##  [37] rappdirs_0.3.1       Rcpp_1.0.3           Biobase_2.46.0      
    ##  [40] cellranger_1.1.0     vctrs_0.2.1          iterators_1.0.12    
    ##  [43] xfun_0.12            ps_1.3.0             testthat_2.3.1      
    ##  [46] rvest_0.3.5          lifecycle_0.1.0      XML_3.99-0.3        
    ##  [49] MASS_7.3-51.5        scales_1.1.0         hms_0.5.3           
    ##  [52] parallel_3.6.2       curl_4.3             yaml_2.2.0          
    ##  [55] memoise_1.1.0        gridExtra_2.3        stringi_1.4.5       
    ##  [58] RSQLite_2.2.0        SQUAREM_2020.1       S4Vectors_0.24.3    
    ##  [61] desc_1.2.0           foreach_1.4.7        BiocGenerics_0.32.0 
    ##  [64] pkgbuild_1.0.6       truncnorm_1.0-8      rlang_0.4.2         
    ##  [67] pkgconfig_2.0.3      matrixStats_0.55.0   bitops_1.0-6        
    ##  [70] evaluate_0.14        lattice_0.20-38      esaBcv_1.2.1        
    ##  [73] bit_1.1-15.1         tidyselect_0.2.5     processx_3.4.1      
    ##  [76] plyr_1.8.5           magrittr_1.5         R6_2.4.1            
    ##  [79] IRanges_2.20.2       generics_0.0.2       DBI_1.1.0           
    ##  [82] pillar_1.4.3         haven_2.2.0          withr_2.1.2         
    ##  [85] mixsqp_0.2-2         survival_3.1-8       RCurl_1.98-1.1      
    ##  [88] modelr_0.1.5         crayon_1.3.4         BiocFileCache_1.10.2
    ##  [91] rmarkdown_2.1        progress_1.2.2       grid_3.6.2          
    ##  [94] readxl_1.3.1         blob_1.2.1           callr_3.4.0         
    ##  [97] reprex_0.3.0         digest_0.6.23        svd_0.5             
    ## [100] openssl_1.4.1        stats4_3.6.2         munsell_0.5.0       
    ## [103] askpass_1.1          sessioninfo_1.1.1

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

Gerard, David, and Matthew Stephens. 2020. “Empirical Bayes shrinkage
and false discovery rate estimation, allowing for unwanted variation.”
*Biostatistics* 21 (1): 15–32.
<https://doi.org/10.1093/biostatistics/kxy029>.

</div>

<div id="ref-lin2017housekeeping">

Lin, Yingxin, Shila Ghazanfar, Dario Strbenac, Andy Wang, Ellis Patrick,
Terence Speed, Jean Yang, and Pengyi Yang. 2017. “Housekeeping Genes,
Revisited at the Single-Cell Level.” *bioRxiv*. Cold Spring Harbor
Laboratory. <https://doi.org/10.1101/229815>.

</div>

</div>
