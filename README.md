# cacoa

<!-- badges: start -->
<!-- badges: end -->

Case-Control Analysis of scRNA-seq experiments

As the method is not published yet, there is no text description available. To get some idea of how it works you may want to check [the slides](https://slides.com/vpetukhov/cacoa-scs-sept-2021) from the method presentation.

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
## Usage

Cacoa currently supports input in several formats (see below). Most of them require the following metadata:

- `sample.groups`: vector with condition labels per sample named with sample ids
- `cell.groups`: cell type annotation vector named by cell ids
- `sample.per.cell`: vector with sample labels per cell named with cell ids
- `ref.level`: id of the condition, corresponding to the reference (i.e. control)
- `target.level`: id of the condition, corresponding to the target (i.e. case)

### Raw or normalized joint count matrix `cm`

```r
cao <- Cacoa$new(cm, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
                 ref.level=ref.level, target.level=target.level, embedding=embedding)
```

Parameter `embedding` contains a matrix or data.frame with a cell embedding. Rownames should match to the cell ids. It is used for visualization and some cluster-free analysis.

### Seurat object `so`

```r
cao <- Cacoa$new(so, sample.groups=sample.groups, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
                 ref.level=ref.level, target.level=target.level, graph.name=graph.name)
```

Parameter `graph.name` is required for cluster-free analysis, and must contain a name of joint graph in Seurat object. For that, the Seurat object must have a joint graph estimated (see [FindNeighbors](https://satijalab.org/seurat/reference/findneighbors)). For visualization purpuses, Seurat also must have cell embedding estimated or the embedding data frame must be provided in the `embedding` parameter.

### Conos object `co`

```r
cao <- Cacoa$new(co, sample.groups=sample.groups, cell.groups=cell.groups, 
                 ref.level=ref.level, target.level=target.level)
```

For visualization purpuses, Conos must have cell embedding estimated or the embedding data frame must be provided in the `embedding` parameter. And for cluster-free analysis it should have a joint graph (see [Conos$buildGraph](https://cran.r-project.org/web/packages/conos/conos.pdf) method).

## Example

See [the vignette](http://pklab.med.harvard.edu/viktor/cacoa/ep.html) for an example.
