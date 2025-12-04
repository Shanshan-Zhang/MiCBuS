
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MiCBuS

Marker Gene Mining for Unknown Cell Types Using Bulk and Single Cell
RNA-Seq Data

![](man/figures/Fig1.png)

A framework of MiCBuS to identify marker genes for unknown cell types.
The unknown cell types, highlighted in orange, are present in bulk
RNA-seq data but remain undetected in incomplete scRNA-seq data. MiCBuS
takes both bulk and single-cell RNA-seq data as input and follows three
major steps. MiCBuS estimates cell type proportions in bulk RNA-seq data
using scRNA-seq as a reference. It then proceeds to generate
Dirichlet-pseudo-bulk RNA-seq samples. By comparing these simulated
Dirichlet-pseudo-bulk RNA-seq samples with the original bulk RNA-seq
data, MiCBuS can identify pseudo marker genes of unknown cell types,
referred to as psMarker.

## Installation

You can install the development version of MiCBuS like so:

``` r
# install devtools if necessary
install.packages('devtools')

# install the MiCBuS package
devtools::install_github('Shanshan-Zhang/MiCBuS')

# load
library(MiCBuS)
```

## Example

For an example of how to use **MiCBuS**, please see the vignette:

- **[Intro to
  MiCBuS](https://shanshan-zhang.github.io/MiCBuS/doc/Intro_to_MiCBuS.html)**

You can also access it locally after installing MiCBuS:

\`\`\`r vignette(“Intro_to_MiCBuS”)

For an example how to use MiCBuS Example datasets are given in the
Package
