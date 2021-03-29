# cacoa

<!-- badges: start -->
<!-- badges: end -->

Case-Control Analysis of scRNA-seq experiments

## Installation

Prior to installing the package, dependencies have to be installed:

```r
BiocManager::install(c("clusterProfiler", "DESeq2", "DOSE", "EnhancedVolcano", "enrichplot", "fabia", "GOfuncR", "Rgraphviz"))
```

Also make sure to install the latest version of sccore (not the one from CRAN):

``` r
devtools::install_github("kharchenkolab/sccore", ref="dev")
```

And if you're going to use Conos as a data object, the latest version is also recommended:

``` r
devtools::install_github("kharchenkolab/conos", ref="dev")
```

To install the package use:

``` r
devtools::install_github("kharchenkolab/cacoa")
```

## Example

See [the vignette](http://pklab.med.harvard.edu/viktor/cacoa/ep.html) for an example.
