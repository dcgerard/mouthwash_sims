
Source to reproduce results from Gerard and Stephens (2017)
===========================================================

This repository contains source code to reproduce the empirical evaluations of Gerard and Stephens (2017). The new methods can be found in the [vicar](https://github.com/dcgerard/vicar) package.

If you find a bug, please create an [issue](https://github.com/dcgerard/ruvb_sims/issues).

Citing this work
----------------

If you find any of the source code in this repository useful for your work, please cite our paper:

> Gerard, David, and Matthew Stephens. 2017. "Empirical Bayes Shrinkage and False Discovery Rate Estimation, Allowing for Unwanted Variation." *arXiv Preprint arXiv:1709.10066*. <https://arxiv.org/abs/1709.10066>.

License
-------

Copyright (c) 2017-2018, David Gerard.

All source code and software in this repository are made available under the terms of the [GNU General Public License](http://www.gnu.org/licenses/gpl.html). See the [LICENSE](LICENSE) file for the full text of the license.

Instructions
------------

To reproduce the results of Gerard and Stephens (2017), you need to (i) install the appropriate R packages, (ii) obtain the appropriate data, (iii) run `make` and (iv) get some coffee while you wait. Please read below for details on each of these steps.

### Install software and R packages

1.  Install [R](https://cran.r-project.org).
2.  Install [GNU Make](https://www.gnu.org/software/make).
3.  Install the required R packages by running the following commands in the R interactive environment. (Note the order of these commands is important---the Bioconductor packages should be installed before the CRAN packages.)

``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("sva", "limma", "qvalue", "biomaRt"), suppressUpdates = TRUE)
install.packages(c("tidyverse", "stringr", "reshape2", "pROC",
                   "ruv", "cate", "devtools", "ashr", "bfa",
                   "xtable", "dplyr", "ggthemes",
                   "assertthat", "R.utils"))
devtools::install_github("dcgerard/seqgendiff")
devtools::install_github("dcgerard/vicar")
```

See below for the versions of the packages that I used for the paper.

### Get data

Download the following files from the [GTEx Portal](https://www.gtexportal.org) and copy them to the [Data](data) subdirectory. In order to replicate my results, it is important to use version "6p" of the GTEx data. To access these files, you will need to register for a GTEx Portal account if you have not done so already.

1.  `GTEx_Data_V6_Annotations_SampleAttributesDS.txt`
2.  `GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz`
3.  `gencode.v19.genes.v6p_model.patched_contigs.gtf.gz`
4.  `GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt`

    Next, download the list of human housekeeping genes from Eisenberg and Levanon (2013) and the NCBI-to-Ensembl gene mapping file, and copy these files to the [Data](data) directory:

5.  <http://www.tau.ac.il/~elieis/HKG/HK_genes.txt>
6.  <ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz>

    Finally, download the list of human housekeeping genes from Lin et al. (2017) by navigating to their [Shiny R application](http://shiny.maths.usyd.edu.au/scHK/), clicking on "Default Values", then immediately clicking "Download". Then place this file, labeled "h-scHKgenes.csv", in the [Data](Data) directory.

7.  `h-scHKgenes.csv`

    I used the file [get\_ensembl.R](./R/get_ensembl.R) to correspond HGNC gene names to their ensembl annotation, but I don't run this in make because I don't want a dependency on the biomaRt package, so I've included the resulting file [lin\_hk\_genes.csv](./Data/lin_hk_genes.csv) in the [Data](Data) directory.

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
    ## hsc.txt~
    ## lin_hk_genes.csv

### Run make

To reproduce all of the results in Gerard and Stephens (2017), move to the root directory in your local copy of this repository, and run `make` from the terminal. This will run all the steps in the data processing and analysis pipeline. Note that you may need to adjust some of the settings in the [Makefile](Makefile) to suit your computing environment.

Since some of the steps take a long time to run, you may not want to run everything at once. If you want to reproduce the simulation results, run

``` bash
make sims
```

If you want to reproduce just the results of the GTEx analysis using the control genes from Eisenberg and Levanon (2013), run

``` bash
make gtex_analysis
```

If you want to reproduce just the results of the GTEx analysis using the control genes from Lin et al. (2017), run

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

At any point during or after the commands are running, it may be helpful to inspect the `.Rout` files saved in the [Output](Output) subdirectory.

### Get Coffee

All of these runs (except the last one) should take a very long time (a day to a couple of days). You should get some coffee. Here is a list of some of my favorite places:

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

If you are having trouble reproducing these results, it might be that you need to update some of your R packages. These are the versions that I used (including some versions of packages that are not actually needed to run the code):

``` r
sessionInfo()
```

    ## R version 3.4.4 (2018-03-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/atlas-base/atlas/libblas.so.3.0
    ## LAPACK: /usr/lib/atlas-base/atlas/liblapack.so.3.0
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
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.4.4  backports_1.0.5 magrittr_1.5    rprojroot_1.3-2
    ##  [5] tools_3.4.4     htmltools_0.3.5 yaml_2.1.14     Rcpp_0.12.16   
    ##  [9] stringi_1.1.7   rmarkdown_1.9   knitr_1.20      stringr_1.3.0  
    ## [13] digest_0.6.12   evaluate_0.10

As you can see, I've also only tried this out on Ubuntu.

Credits
-------

This project was developed by [David Gerard](https://dcgerard.github.io) at the University of Chicago.

[Peter Carbonetto](https://pcarbo.github.io/) made many fantastic contributions and suggestions to increase the reproducibility of this work.

Thanks to [Matthew Stephens](http://stephenslab.uchicago.edu) for his support and mentorship.

References
----------

Eisenberg, Eli, and Erez Y Levanon. 2013. “Human Housekeeping Genes, Revisited.” *Trends in Genetics* 29 (10). Elsevier: 569–74. doi:[10.1016/j.tig.2013.05.010](https://doi.org/10.1016/j.tig.2013.05.010).

Gerard, David, and Matthew Stephens. 2017. “Unifying and Generalizing Methods for Removing Unwanted Variation Based on Negative Controls.” *arXiv Preprint arXiv:1705.08393*. <https://arxiv.org/abs/1705.08393>.

Lin, Yingxin, Shila Ghazanfar, Dario Strbenac, Andy Wang, Ellis Patrick, Terence Speed, Jean Yang, and Pengyi Yang. 2017. “Housekeeping Genes, Revisited at the Single-Cell Level.” *bioRxiv*. Cold Spring Harbor Laboratory. doi:[10.1101/229815](https://doi.org/10.1101/229815).
