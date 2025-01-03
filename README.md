
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

``` r
runCrosstalkeR()
```

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

Refer to package vignettes for more details.

An overview of the package is illustrated below.

![](./inst/extdata/Ziyue_Li.png)

## Contribution

The author of package is Ziyue Li.

The author write **prepare_interaction_data** to processes the Gene
Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG)
enrichment results, filtering for pathways with a q-value less than 0.1.
It then maps the significant pathways to the corresponding genes,
focusing only on differentially expressed genes (DEGs) by filtering
genes based on their Entrez gene IDs. The function returns a list of
pathways with their associated DEG gene sets, which can be used for
further analysis.

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

Also, OpenAI’s ChatGPT, were used to help structure the Roxygen
documentation and provide suggestions for optimizing the code logic.

## Reference

Silva A. (2019). GitHub - anjalisilva/TestingPackage: R Package
Illustrating Components of an R package for BCB410H - Applied
Bioinformatics (2019-2023), University of Toronto, Canada. GitHub.
<https://github.com/anjalisilva/TestingPackage>

OpenAI. (2024). ChatGPT (Version 3.5). Retrieved from
<https://chat.openai.com/chat>

Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change
and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550.
<doi:10.1186/s13059-014-0550-8>.

Thomas PD, Ebert D, Muruganujan A, Mushayahama T, Albou LP, Mi H (2022).
“PANTHER: Making genome-scale phylogenetics accessible to all.” Protein
Science, 31(1), 8-22. <doi:10.1002/pro.4218>.

Kanehisa M, Goto S, Sato Y, et al. (2012). “KEGG for integration and
interpretation of large-scale molecular data sets.” Nucleic Acids
Research, 40(D1), D109-D114. <doi:10.1093/nar/gkr988>.

## Acknowledgement

This package was developed as part of an assessment for 2024 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `CrosstalkeR` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.

## Example

This is a basic example which shows you how to utilize the functions:


    # Example Usage: Assuming DESeq2 and enrichment data is in the 'data' file
    # Load data (Replace 'data' with your actual file names)

    diff_genes <- read.csv("data/diff_gene_deseq2_lactate.csv")  # DESeq2 results with annotated DEGs
    go_enrich <- readRDS("data/go_enrich_results.rds")           # Precomputed GO enrichment results
    kegg_enrich <- readRDS("data/kegg_enrich_results.rds")       # Precomputed KEGG enrichment results

    # Prepare interaction data
    pathway_gene_map <- prepare_interaction_data(go_enrich, kegg_enrich, diff_genes)

    # Plot pathway-gene interaction network
    plot_pathway_gene_interactions(pathway_gene_map)

    # Plot ranked pathways
    rank_pathways(pathway_gene_map)
