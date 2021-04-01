
<!-- README.md is generated from README.Rmd. Please edit that file -->

# noderank

<!-- badges: start -->

<!-- badges: end -->

noderank maps differential gene expression results onto an established
network of protein-protein associations and employs network analysis
methods to select a prioritized list of important nodes. For data in
which the mechanism of perturbation is known, each ranked list result is
scored by the position of the causal gene.

## Installation

You can install the released version of noderank from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("noderank")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(noderank)
## basic example code
```

## A roadmap for noderank:

Starting material: differential gene expression analysis results,
i.e. the output of topTable() from *limma* or results() from *DESeq2*.

1.  
<!-- end list -->

``` r
# summary(cars)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
