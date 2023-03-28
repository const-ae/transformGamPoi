
# transformGamPoi

<!-- badges: start -->
<!-- badges: end -->

R package that accompanies our paper ‘Transformation and Preprocessing
of Single-Cell RNA-Seq Data’
(<https://www.biorxiv.org/content/10.1101/2021.06.24.449781v1>).

`transformGamPoi` provides methods to stabilize the variance of single
cell count data:

- acosh transformation based on the delta method
- shifted logarithm (log(x + c)) with a pseudo-count c, so that it
  approximates the acosh transformation
- randomized quantile and Pearson residuals

## Installation

You can install the current development version of `transformGamPoi` by
typing the following into the [*R*](https://cloud.r-project.org/)
console:

``` r
# install.packages("devtools")
devtools::install_github("const-ae/transformGamPoi")
```

The installation should only take a few seconds and work across all
major operating systems (MacOS, Linux, Windows).

## Example

Let’s compare the different variance-stabilizing transformations.

We start by loading the `transformGamPoi` package and setting a seed to
make sure the results are reproducible.

``` r
library(transformGamPoi)
set.seed(1)
```

We then load some example data, which we subset to 1000 genes and 500
cells

``` r
sce <- TENxPBMCData::TENxPBMCData("pbmc4k")
#> snapshotDate(): 2022-10-31
#> Warning: package 'GenomicRanges' was built under R version 4.2.2
#> Warning: package 'S4Vectors' was built under R version 4.2.2
#> Warning: package 'GenomeInfoDb' was built under R version 4.2.2
#> see ?TENxPBMCData and browseVignettes('TENxPBMCData') for documentation
#> loading from cache
sce_red <- sce[sample(which(rowSums2(counts(sce)) > 0), 1000),
               sample(ncol(sce), 500)]
```

We calculate the different variance-stabilizing transformations. We can
either use the generic `transformGamPoi()` method and specify the
`transformation`, or we use the low-level functions `acosh_transform()`,
`shifted_log_transform()`, and `residual_transform()` which provide more
settings. All functions return a matrix, which we can for example insert
back into the `SingleCellExperiment` object:

``` r
assay(sce_red, "acosh") <- transformGamPoi(sce_red, transformation = "acosh")
assay(sce_red, "shifted_log") <- shifted_log_transform(sce_red, overdispersion = 0.1)
# For large datasets, we can also do the processing without 
# loading the full dataset into memory (on_disk = TRUE)
assay(sce_red, "rand_quant") <- residual_transform(sce_red, "randomized_quantile", on_disk = FALSE)
assay(sce_red, "pearson") <- residual_transform(sce_red, "pearson", clipping = TRUE, on_disk = FALSE)
```

Finally, we compare the variance of the genes after transformation using
a scatter plot

``` r
par(pch = 20, cex = 1.15)
mus <- rowMeans2(counts(sce_red))
plot(mus, rowVars(assay(sce_red, "acosh")), log = "x", col = "#1b9e77aa", cex = 0.6,
     xlab =  "Log Gene Means", ylab = "Variance after transformation")
points(mus, rowVars(assay(sce_red, "shifted_log")), col = "#d95f02aa", cex = 0.6)
points(mus, rowVars(assay(sce_red, "pearson")), col = "#7570b3aa", cex = 0.6)
points(mus, rowVars(assay(sce_red, "rand_quant")), col = "#e7298aaa", cex = 0.6)
legend("topleft", legend = c("acosh", "shifted log", "Pearson Resid.", "Rand. Quantile Resid."),
       col = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"), pch = 16)
```

![](man/figures/README-plotMeanVar-1.png)<!-- -->

# See also

There are a number of preprocessing methods and packages out there. Of
particular interests are

- [sctransform](https://github.com/ChristophH/sctransform) by Christoph
  Hafemeister and the [Satija lab](https://satijalab.org/). The original
  developers of the Pearson residual variance-stabilizing transformation
  approach for single cell data.
- [scuttle::logNormCounts()](https://bioconductor.org/packages/release/bioc/html/scuttle.html)
  by Aaron Lun. This is an alternative to the `shifted_log_transform()`
  and plays nicely together with the Bioconductor universe. For more
  information, I highly recommend to take a look at the
  [normalization](https://bioconductor.org/books/release/OSCA/normalization.html)
  section of the [OSCA
  book](https://bioconductor.org/books/release/OSCA/).
- [Sanity](https://github.com/jmbreda/Sanity) by Jérémie Breda *et al.*.
  This method is not directly concerned with variance stabilization, but
  still provides an interesting approach for single cell data
  preprocessing.

# Session Info

``` r
sessionInfo()
#> R version 4.2.1 RC (2022-06-17 r82503)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur ... 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] TENxPBMCData_1.16.0         HDF5Array_1.26.0           
#>  [3] rhdf5_2.42.0                DelayedArray_0.24.0        
#>  [5] Matrix_1.5-3                SingleCellExperiment_1.20.0
#>  [7] SummarizedExperiment_1.28.0 Biobase_2.58.0             
#>  [9] GenomicRanges_1.50.2        GenomeInfoDb_1.34.9        
#> [11] IRanges_2.32.0              S4Vectors_0.36.2           
#> [13] BiocGenerics_0.44.0         MatrixGenerics_1.10.0      
#> [15] matrixStats_0.63.0          transformGamPoi_1.4.0      
#> 
#> loaded via a namespace (and not attached):
#>  [1] httr_1.4.5                    bit64_4.0.5                  
#>  [3] AnnotationHub_3.6.0           DelayedMatrixStats_1.20.0    
#>  [5] shiny_1.7.4                   interactiveDisplayBase_1.36.0
#>  [7] highr_0.10                    BiocManager_1.30.20          
#>  [9] BiocFileCache_2.6.1           blob_1.2.3                   
#> [11] GenomeInfoDbData_1.2.9        yaml_2.3.7                   
#> [13] BiocVersion_3.16.0            pillar_1.8.1                 
#> [15] RSQLite_2.3.0                 lattice_0.20-45              
#> [17] glue_1.6.2                    digest_0.6.31                
#> [19] promises_1.2.0.1              XVector_0.38.0               
#> [21] htmltools_0.5.4               httpuv_1.6.9                 
#> [23] pkgconfig_2.0.3               zlibbioc_1.44.0              
#> [25] purrr_1.0.1                   xtable_1.8-4                 
#> [27] later_1.3.0                   tibble_3.1.8                 
#> [29] KEGGREST_1.38.0               generics_0.1.3               
#> [31] ellipsis_0.3.2                withr_2.5.0                  
#> [33] cachem_1.0.7                  cli_3.6.0                    
#> [35] magrittr_2.0.3                crayon_1.5.2                 
#> [37] mime_0.12                     memoise_2.0.1                
#> [39] evaluate_0.20                 fansi_1.0.4                  
#> [41] tools_4.2.1                   lifecycle_1.0.3              
#> [43] Rhdf5lib_1.20.0               AnnotationDbi_1.60.0         
#> [45] Biostrings_2.66.0             compiler_4.2.1               
#> [47] rlang_1.0.6                   grid_4.2.1                   
#> [49] RCurl_1.98-1.10               rhdf5filters_1.10.0          
#> [51] rstudioapi_0.14               rappdirs_0.3.3               
#> [53] glmGamPoi_1.11.7              bitops_1.0-7                 
#> [55] rmarkdown_2.20                ExperimentHub_2.6.0          
#> [57] DBI_1.1.3                     curl_5.0.0                   
#> [59] R6_2.5.1                      knitr_1.42                   
#> [61] dplyr_1.1.0                   fastmap_1.1.1                
#> [63] bit_4.0.5                     utf8_1.2.3                   
#> [65] filelock_1.0.2                Rcpp_1.0.10                  
#> [67] vctrs_0.5.2                   png_0.1-8                    
#> [69] sparseMatrixStats_1.10.0      dbplyr_2.3.1                 
#> [71] tidyselect_1.2.0              xfun_0.37
```
