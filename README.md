
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

## A roadmap for noderank:

Input: Differential gene expression analysis results including log
fold-changes and p-values, i.e. the output of topTable() from *limma* or
results() from *DESeq2*.

1.  Generate protein-protein interaction network.
2.  Map DEA results onto the PPI.
3.  Prune the network by user-specified thresholds for fold-change and
    p-value.
4.  Implement user-specified network analysis algorithm to rank the
    nodes of the network.
5.  Calculate overall network score and uncertainty metrics.

Output:

  - The final network composed of the important genes.

  - A dataframe of important genes with their associated ranking
    metrics.

  - The mean score and p-value for the overall network.
