#' Plot Pathway-Gene Interactions
#'
#' Visualizes the relationship between pathways and DEGs as a network graph.
#'
#' @param pathway_gene_map A named list created by `prepare_interaction_data`, where each name is a pathway
#'   and its value is a vector of associated gene IDs.
#' @return A ggplot object representing the interaction network.
#' @examples
#' # Example pathway-gene map:
#' pathway_gene_map <- list("Pathway1" = c("Gene1", "Gene2"), "Pathway2" = c("Gene3"))
#' plot_pathway_gene_interactions(pathway_gene_map)
#'
#' @export
plot_pathway_gene_interactions <- function(pathway_gene_map) {
  # Create edges for the graph
  edges <- data.frame(
    from = rep(names(pathway_gene_map), sapply(pathway_gene_map, length)),
    to = unlist(pathway_gene_map)
  )

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
