#' Rank Pathways by DEG Involvement
#'
#' Ranks pathways based on the number of associated DEGs and visualizes the results as a bar plot.
#'
#' @param pathway_gene_map A named list created by `prepare_interaction_data`, where each name is a pathway
#'   and its value is a vector of associated gene IDs.
#' @return A ggplot object displaying a bar plot of pathways ranked by DEG count.
#' @examples
#' # Example pathway-gene map:
#' pathway_gene_map <- list("Pathway1" = c("Gene1", "Gene2"), "Pathway2" = c("Gene3"))
#' rank_pathways(pathway_gene_map)
#'
#' @export
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
