








library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(igraph)
library(ggplot2)


# Prepare_interaction_data function
prepare_interaction_data <- function(go_enrich, kegg_enrich, diff_genes) {
  # Ensure they are data.frame and filter by qvalue
  significant_go <- go_enrich[go_enrich$qvalue < 0.1, ]
  significant_kegg <- kegg_enrich[kegg_enrich$qvalue < 0.1, ]

  # Merge pathways with DEG data
  pathway_gene_map <- list()

  # GO Pathway Mapping
  for (i in 1:nrow(significant_go)) {
    genes <- strsplit(significant_go[i, "geneID"], "/")[[1]]
    pathway_gene_map[[significant_go[i, "Description"]]] <- genes
  }

  # KEGG Pathway Mapping
  for (i in 1:nrow(significant_kegg)) {
    genes <- strsplit(significant_kegg[i, "geneID"], "/")[[1]]
    pathway_gene_map[[significant_kegg[i, "Description"]]] <- genes
  }

  # Filter to DEGs only
  deg_entrez <- diff_genes$entrezgene_id
  pathway_gene_map <- lapply(pathway_gene_map, function(genes) {
    intersect(genes, deg_entrez)
  })

  return(pathway_gene_map)
}




# Visualize Pathway-Gene Interactions as a Graph
plot_pathway_gene_interactions <- function(pathway_gene_map) {
  # Check if pathway_gene_map is empty or not
  if(length(pathway_gene_map) == 0) {
    stop("Pathway-Gene map is empty. Please check the enrichment results and DEG data.")
  }

  # Print the pathway_gene_map to inspect its structure
  print("Pathway-Gene Map:")
  print(pathway_gene_map)

  # Create edges for the graph
  edges <- data.frame(
    from = rep(names(pathway_gene_map), sapply(pathway_gene_map, length)),
    to = unlist(pathway_gene_map)
  )

  # Check if edges have at least two columns
  if (ncol(edges) < 2) {
    stop("The edges data frame should contain at least two columns: 'from' (pathway) and 'to' (gene).")
  }

  # Print edges to check its structure
  print("Edges DataFrame:")
  print(head(edges))

  # Create an igraph object
  g <- graph_from_data_frame(edges)

  # Plot the graph
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(alpha = 0.5), show.legend = FALSE) +
    geom_node_point(size = 5, color = "steelblue") +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_minimal() +
    labs(title = "Pathway-Gene Interaction Network")
}

# Rank Pathways by DEG Involvement
rank_pathways <- function(pathway_gene_map) {
  # Count the number of genes associated with each pathway
  pathway_ranks <- data.frame(
    Pathway = names(pathway_gene_map),
    GeneCount = sapply(pathway_gene_map, length)
  )

  # Sort by GeneCount
  pathway_ranks <- pathway_ranks[order(-pathway_ranks$GeneCount), ]

  # Plot barplot
  ggplot(pathway_ranks, aes(x = reorder(Pathway, -GeneCount), y = GeneCount)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Pathway Rankings by DEG Involvement",
         x = "Pathway",
         y = "Number of DEGs")
}

# Example Usage: Assuming DESeq2 and enrichment data is in the 'data' file
# Load data (Replace 'data' with your actual file names)
diff_genes <- read.csv("data/diff_gene_deseq2_lactate.csv")  # DESeq2 results with annotated DEGs
go_enrich <- readRDS("data/go_enrich_results.rds")           # Precomputed GO enrichment results
kegg_enrich <- readRDS("data/kegg_enrich_results.rds")       # Precomputed KEGG enrichment results

class(go_enrich)
class(kegg_enrich)



# Prepare interaction data
pathway_gene_map <- prepare_interaction_data(go_enrich, kegg_enrich, diff_genes)

# Plot pathway-gene interaction network
plot_pathway_gene_interactions(pathway_gene_map)

# Plot ranked pathways
rank_pathways(pathway_gene_map)
