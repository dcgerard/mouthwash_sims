# Source to reproduce results from Gerard & Stephens (2017)

This repository contains source code to reproduce the empirical
evaluations of [Gerard & Stephens (2017)](https://arxiv.org/abs/1709.10066).
The new methods can be found in the
[vicar](https://github.com/dcgerard/vicar) package.

If you find a bug, please create an
[issue](https://github.com/dcgerard/ruvb_sims/issues).

This code has been tested in ...

## Citing this work

If you find any of the source code in this repository useful for your
work, please cite our paper:

> Gerard, David, and Matthew Stephens. 2017. "Empirical Bayes Shrinkage
> and False Discovery Rate Estimation, Allowing for Unwanted Variation."
> *arXiv Preprint arXiv:1709.10066*. <https://arxiv.org/abs/1709.10066>.

## License

Copyright (c) 2017-2018, David Gerard.

All source code and software in this repository are made available
under the terms of the [GNU General Public
License](http://www.gnu.org/licenses/gpl.html). See the
[LICENSE](LICENSE) file for the full text of the license.

## Instructions

To reproduce the results of Gerard & Stephens (2017), you need to (1)
install the appropriate R packages, (2) obtain the appropriate data,
(3) run `make` and (4) get some coffee while you wait. Please read
below for details on each of these steps.

### Install software and R packages

1. Install [R](https://cran.r-project.org).

2. Install [GNU Make](https://www.gnu.org/software/make).

3. Install the required R packages by running the following commands
in the R interactive environment. (Note the order of these commands is
important---the Bioconductor packages should be installed before the
CRAN packages.)

```R
source("https://bioconductor.org/biocLite.R")
biocLite(c("sva","limma","qvalue"),suppressUpdates = TRUE)
install.packages(c("tidyverse", "stringr", "reshape2", "pROC",
                   "ruv", "cate", "devtools", "ashr", "bfa",
				   "xtable", "dplyr", "ggthemes", "qvalue",
				   "assertthat", "R.utils"))
devtools::install_github("dcgerard/seqgendiff")
devtools::install_github("dcgerard/vicar")
```

See below for the versions of the packages that I used for the paper.

### Get data

Download the following files from the
[GTEx Portal](https://www.gtexportal.org) and copy them to the
[Data](data) subdirectory. In order to replicate my results, it is
important to use version "6p" of the GTEx data. To access these files,
you will need to register for a GTEx Portal account if you have not
done so already.

1. `GTEx_Data_V6_Annotations_SampleAttributesDS.txt`
2. `GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz`
3. `gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`
4. `GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt`

Next, download the list of human housekeeping genes and the NCBI
NCBI-to-Ensembl gene mapping file, and copy these files to the
[Data](data) directory:

5. http://www.tau.ac.il/~elieis/HKG/HK_genes.txt
6. ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz

After completing these steps, the Data folder should look like this:

```bash
cd Data
ls -1
# GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz
# GTEx_Data_V6_Annotations_SampleAttributesDS.txt
# GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt
# HK_genes.txt
# gencode.v19.genes.v6p_model.patched_contigs.gtf.gz
# gene2ensembl.gz
```

### Run make

To reproduce all of the results in Gerard and Stephens (2017), move to
the root directory in your local copy of this repository, and run
`make` from the terminal. This will run all the steps in the data
processing and analysis pipeline. Note that you may need to adjust
some of the settings in the [Makefile](Makefile) to suit your
computing environment.

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

At any point during or after the commands are running, you can assess
progress and troubleshoot errors by inspecting the `.Rout` files saved
in the [Output](Output) subdirectory.

### Get Coffee

All of these runs (except the last one) should take a very long time
(a day to a couple of days). You should get some coffee. Here is a
list of some of my favorite places:

-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Stauf's Coffee Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

If you are having trouble reproducing these results, it might be that
you need to update some of your R packages. These are the versions
that I used (including some versions of packages that are not actually
needed to run the code):

```R
sessionInfo()
```

    ## R version 3.4.0 (2017-04-21)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 15063)
    ##
    ## Matrix products: default
    ##
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252
    ## [2] LC_CTYPE=English_United States.1252
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C             
    ## [5] LC_TIME=English_United States.1252
    ##
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base
    ##
    ## other attached packages:
    ##  [1] ashr_2.0.5          seqgendiff_0.1.0    Rmpi_0.6-6
    ##  [4] vicar_0.1.6         limma_3.32.10       sva_3.24.4
    ##  [7] BiocParallel_1.10.1 genefilter_1.58.0   mgcv_1.8-17
    ## [10] nlme_3.1-131        devtools_1.13.4     snow_0.4-2
    ## [13] gridExtra_2.3       cate_1.0.4          ruv_0.9.6
    ## [16] pROC_1.10.0         reshape2_1.4.3      forcats_0.2.0
    ## [19] stringr_1.2.0       dplyr_0.7.4         purrr_0.2.4
    ## [22] readr_1.1.1         tidyr_0.7.2         tibble_1.3.4
    ## [25] ggplot2_2.2.1       tidyverse_1.2.1
    ##
    ## loaded via a namespace (and not attached):
    ##  [1] svd_0.4              bitops_1.0-6         matrixStats_0.52.2
    ##  [4] lubridate_1.7.1      doParallel_1.0.10    httr_1.3.1
    ##  [7] rprojroot_1.2        tools_3.4.0          backports_1.1.1
    ## [10] R6_2.2.0             DBI_0.6-1            lazyeval_0.2.0
    ## [13] BiocGenerics_0.22.0  colorspace_1.3-2     withr_1.0.2
    ## [16] mnormt_1.5-5         compiler_3.4.0       cli_1.0.0
    ## [19] rvest_0.3.2          Biobase_2.36.0       xml2_1.1.1
    ## [22] scales_0.4.1         SQUAREM_2017.10-1    psych_1.7.3.21
    ## [25] esaBcv_1.2.1         digest_0.6.12        foreign_0.8-67
    ## [28] rmarkdown_1.6        pscl_1.4.9           pkgconfig_2.0.1
    ## [31] htmltools_0.3.6      rlang_0.1.4          readxl_1.0.0
    ## [34] rstudioapi_0.7       RSQLite_1.1-2        bindr_0.1
    ## [37] jsonlite_1.5         leapp_1.2            RCurl_1.95-4.8
    ## [40] magrittr_1.5         Matrix_1.2-9         Rcpp_0.12.14
    ## [43] munsell_0.4.3        S4Vectors_0.14.0     stringi_1.1.5
    ## [46] yaml_2.1.14          MASS_7.3-47          plyr_1.8.4
    ## [49] grid_3.4.0           parallel_3.4.0       crayon_1.3.4
    ## [52] lattice_0.20-35      haven_1.1.0          splines_3.4.0
    ## [55] annotate_1.54.0      hms_0.3              knitr_1.17
    ## [58] corpcor_1.6.9        codetools_0.2-15     stats4_3.4.0
    ## [61] XML_3.98-1.6         glue_1.1.1           evaluate_0.10
    ## [64] modelr_0.1.1         foreach_1.4.3        cellranger_1.1.0
    ## [67] gtable_0.2.0         assertthat_0.2.0     xtable_1.8-2
    ## [70] broom_0.4.2          survival_2.41-3      truncnorm_1.0-7
    ## [73] iterators_1.0.8      AnnotationDbi_1.38.0 memoise_1.1.0
    ## [76] IRanges_2.10.0       bindrcpp_0.2

As you can see, I've also only tried this out on Ubuntu.

## Credits

This project was developed by
[David Gerard](https://dcgerard.github.io) at the University of
Chicago.

Thanks to [Matthew Stephens](stephenslab.uchicago.edu) for his support
and mentorship.

