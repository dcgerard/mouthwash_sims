
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

> Gerard, David, and Matthew Stephens. 2020. “Empirical Bayes Shrinkage
> and False Discovery Rate Estimation, Allowing for Unwanted Variation.”
> *Biostatistics*, 21(1), 15-32. doi:
> [10.1093/biostatistics/kxy029](https://doi.org/10.1093/biostatistics/kxy029).

## License

Copyright (c) 2017-2020, David Gerard.

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
    “Default Values,” then immediately clicking “Download.” Then place
    this file, labeled “h-scHKgenes.csv,” in the [Data](Data) directory.

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

-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air
        Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus
        Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Stauf’s Coffee
        Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe
        Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

If you are having trouble reproducing these results, it might be that
you need to update some of your R packages. These are the versions that
I used (including some versions of packages that are not actually needed
to run the code):

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.2 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
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
    ##  [1] vicar_0.1-11        seqgendiff_1.2.2    biomaRt_2.46.3     
    ##  [4] qvalue_2.22.0       limma_3.46.0        sva_3.38.0         
    ##  [7] BiocParallel_1.24.1 genefilter_1.72.1   mgcv_1.8-35        
    ## [10] nlme_3.1-152        R.utils_2.10.1      R.oo_1.24.0        
    ## [13] R.methodsS3_1.8.1   assertthat_0.2.1    ggthemes_4.2.4     
    ## [16] xtable_1.8-4        bfa_0.4             ashr_2.2-47        
    ## [19] devtools_2.4.0      usethis_2.0.1       cate_1.1.1         
    ## [22] ruv_0.9.7.1         pROC_1.17.0.1       reshape2_1.4.4     
    ## [25] forcats_0.5.1       stringr_1.4.0       dplyr_1.0.5        
    ## [28] purrr_0.3.4         readr_1.4.0         tidyr_1.1.3        
    ## [31] tibble_3.1.1        ggplot2_3.3.3       tidyverse_1.3.1    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] colorspace_2.0-0     ellipsis_0.3.2       rprojroot_2.0.2     
    ##   [4] corpcor_1.6.9        fs_1.5.0             rstudioapi_0.13     
    ##   [7] remotes_2.3.0        leapp_1.2            bit64_4.0.5         
    ##  [10] AnnotationDbi_1.52.0 fansi_0.4.2          lubridate_1.7.10    
    ##  [13] xml2_1.3.2           splines_4.0.5        cachem_1.0.4        
    ##  [16] knitr_1.33           pkgload_1.2.1        jsonlite_1.7.2      
    ##  [19] broom_0.7.6          annotate_1.68.0      dbplyr_2.1.1        
    ##  [22] compiler_4.0.5       httr_1.4.2           backports_1.2.1     
    ##  [25] Matrix_1.3-2         fastmap_1.1.0        cli_2.5.0           
    ##  [28] htmltools_0.5.1.1    prettyunits_1.1.1    tools_4.0.5         
    ##  [31] coda_0.19-4          gtable_0.3.0         glue_1.4.2          
    ##  [34] rappdirs_0.3.3       Rcpp_1.0.6           Biobase_2.50.0      
    ##  [37] cellranger_1.1.0     vctrs_0.3.7          xfun_0.22           
    ##  [40] ps_1.6.0             testthat_3.0.2       rvest_1.0.0         
    ##  [43] irlba_2.3.3          lifecycle_1.0.0      XML_3.99-0.6        
    ##  [46] edgeR_3.32.1         MASS_7.3-53.1        scales_1.1.1        
    ##  [49] hms_1.0.0            parallel_4.0.5       curl_4.3            
    ##  [52] yaml_2.2.1           memoise_2.0.0        gridExtra_2.3       
    ##  [55] stringi_1.5.3        RSQLite_2.2.7        SQUAREM_2021.1      
    ##  [58] S4Vectors_0.28.1     desc_1.3.0           BiocGenerics_0.36.1 
    ##  [61] pkgbuild_1.2.0       truncnorm_1.0-8      rlang_0.4.10        
    ##  [64] pkgconfig_2.0.3      matrixStats_0.58.0   invgamma_1.1        
    ##  [67] evaluate_0.14        lattice_0.20-41      esaBcv_1.2.1        
    ##  [70] bit_4.0.4            tidyselect_1.1.0     processx_3.5.1      
    ##  [73] plyr_1.8.6           magrittr_2.0.1       R6_2.5.0            
    ##  [76] IRanges_2.24.1       generics_0.1.0       DBI_1.1.1           
    ##  [79] pillar_1.6.0         haven_2.4.1          withr_2.4.2         
    ##  [82] mixsqp_0.3-43        survival_3.2-10      modelr_0.1.8        
    ##  [85] crayon_1.4.1         utf8_1.2.1           BiocFileCache_1.14.0
    ##  [88] rmarkdown_2.7        progress_1.2.2       locfit_1.5-9.4      
    ##  [91] grid_4.0.5           readxl_1.3.1         blob_1.2.1          
    ##  [94] callr_3.7.0          reprex_2.0.0         digest_0.6.27       
    ##  [97] svd_0.5              openssl_1.4.3        stats4_4.0.5        
    ## [100] munsell_0.5.0        askpass_1.1          sessioninfo_1.1.1

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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-eisenberg2013human" class="csl-entry">

Eisenberg, Eli, and Erez Y Levanon. 2013. “Human Housekeeping Genes,
Revisited.” *Trends in Genetics* 29 (10): 569–74.
<https://doi.org/10.1016/j.tig.2013.05.010>.

</div>

<div id="ref-gerard2020empirical" class="csl-entry">

Gerard, David, and Matthew Stephens. 2020. “<span
class="nocase">Empirical Bayes shrinkage and false discovery rate
estimation, allowing for unwanted variation</span>.” *Biostatistics* 21
(1): 15–32. <https://doi.org/10.1093/biostatistics/kxy029>.

</div>

<div id="ref-lin2017housekeeping" class="csl-entry">

Lin, Yingxin, Shila Ghazanfar, Dario Strbenac, Andy Wang, Ellis Patrick,
Terence Speed, Jean Yang, and Pengyi Yang. 2017. “Housekeeping Genes,
Revisited at the Single-Cell Level.” *bioRxiv*.
<https://doi.org/10.1101/229815>.

</div>

</div>
