
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CrosstalkeR

<!-- badges: start -->
<!-- badges: end -->

The goal of `CrosstalkeR` is to deduce possible pathway involved in a
certain cellular phenomenon from RNA-Sequence data, and then indicate
possible target downstream or upstream gene for research.

## Description

`CrosstalkR` is a tool that analyzes RNA-Seq data to identify potential
pathways involved in a cellular phenomenon. It predicts key upstream and
downstream genes, helping researchers uncover biological processes,
explore disease mechanisms, and identify therapeutic targets. By
integrating gene expression with pathway analysis, `CrosstalkR` enables
more informed hypothesis generation and experimental design.

## Installation

To install the latest version of the package:

``` r
# install.packages("pak")
install.packages("devtools")
library("devtools")
devtools::install_github("ziyue-li317/CrosstalkeR", build_vignettes = TRUE)
library("CrosstalkeR")
```

To run the Shiny app:

## Overview

``` r
ls("package:<CrosstalkeR>")
data(package = "<CrosstalkeR>") 
browseVignettes("<CrosstalkeR>")
```

`CrosstalkeR` contains three functions:

1.**prepare_interaction_data**: This function extracts significant GO
and KEGG pathways, filters them by q-value, and maps each pathway to
associated genes. It then intersects the gene list with differentially
expressed genes (DEGs) to create a final pathway-gene mapping.

2.**plot_pathway_gene_interactions**: This function visualizes the
interaction between pathways and genes by creating a network graph. It
uses igraph to plot pathways as nodes and DEGs as connected genes,
providing insights into gene-pathway relationships with a visually
appealing layout.

3.**rank_pathways**: This function ranks pathways by the number of DEGs
associated with each. It creates a bar plot to visualize the ranking,
helping identify the most relevant pathways based on DEG involvement.
Pathways are sorted by gene count and plotted in descending order.

## Contribution

The author of package is Ziyue Li. The author write
**prepare_interaction_data** to processes the Gene Ontology (GO) and
Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment results,
filtering for pathways with a q-value less than 0.1. It then maps the
significant pathways to the corresponding genes, focusing only on
differentially expressed genes (DEGs) by filtering genes based on their
Entrez gene IDs. The function returns a list of pathways with their
associated DEG gene sets, which can be used for further analysis.

Developed by the author, **plot_pathway_gene_interactions** visualizes
the interaction network between pathways and genes. It creates a graph
using the igraph package, where the nodes represent genes and pathways,
and edges represent the association between them. The function uses the
ggraph package to generate a visually appealing plot, which helps in
understanding the complex relationships between DEGs and biological
pathways.

**rank_pathways** ranks pathways based on the number of associated DEGs
and generates a bar plot to visualize the top pathways by gene
involvement. The function provides a clear ranking of pathways, aiding
in the identification of the most relevant pathways in the context of
DEG analysis.

## Example

This is a basic example which shows you how to solve a common problem:

``` r

## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
