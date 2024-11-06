# Load required packages
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(dplyr)
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(crosstalkR)  # Ensure this is installed
library(igraph)

# Load DESeq2 and Pathway Results (from previous steps)
diff_gene_deseq2_lactate <- read.csv("diff_gene_deseq2_lactate.csv", row.names = 1)
kegg_enrich <- read.csv("kegg_enrichment_results.csv", row.names = 1)
go_enrich <- read.csv("go_enrichment_results.csv", row.names = 1)







#[1]Prepare Data for Interaction Visualization ---

# Select significant genes from DESeq2 results
sig_genes <- diff_gene_deseq2_lactate %>%
  filter(padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1)) %>%
  select(entrezgene_id, log2FoldChange, padj) %>%
  arrange(padj)

# Get the pathway names from KEGG enrichment (top 10 pathways)
top_kegg_pathways <- head(kegg_enrich, 10)

# Merge KEGG pathways with significant genes
# Assuming that the KEGG pathway genes are available and match entrezgene_id
kegg_pathway_genes <- merge(sig_genes, top_kegg_pathways, by = "entrezgene_id", all = FALSE)






# [2]Create Gene-Pathway Interaction Network ---
# Define edges (interaction between gene and pathway)
edges <- data.frame(
  from = kegg_pathway_genes$entrezgene_id,
  to = kegg_pathway_genes$pathway_id
)

# Create an igraph object for interaction network
graph <- graph_from_data_frame(edges, directed = FALSE)

#Plot Gene-Pathway Interaction Network
# Plot the interaction network
plot(graph, vertex.size = 10, vertex.label.cex = 0.7,
     vertex.color = "lightblue", edge.color = "gray",
     main = "Gene-Pathway Interaction Network (Top KEGG Pathways)")







# [3] Rank Pathways and Visualize ---
# Rank pathways based on enrichment p-values and visualize the top pathways
ranked_pathways <- top_kegg_pathways %>%
  arrange(pvalue) %>%
  head(10)

# Create a ranked barplot for pathways
barplot(-log10(ranked_pathways$pvalue), names.arg = ranked_pathways$pathway,
        col = "blue", main = "Top Ranked KEGG Pathways", las = 2)








# [4]Generate Schematic Figure Showing Gene-Pathway Interaction ---
# visualize the top 5 KEGG pathways
top_5_pathways <- head(ranked_pathways, 5)



# [5] Create Interactive Visualization with crosstalkR ---
# Create an interactive network using crosstalkR
# For the crosstalkR network, create a list of selected genes and pathways

network_data <- data.frame(
  gene = rep(sig_genes$entrezgene_id, each = 5),
  pathway = rep(top_5_pathways$pathway, times = nrow(sig_genes))
)

# Create an interactive network using crosstalkR (simplified)
crosstalk_plot <- crosstalkR::createNetwork(
  data = network_data, source = "gene", target = "pathway",
  type = "interaction", directed = TRUE
)





# [6] Display and Save Interactive Plot ---
# Display the interactive plot
crosstalk_plot

# [7]save the plot for later use
ggsave("gene_pathway_interaction_network.png", plot = crosstalk_plot)
write.csv(network_data, "gene_pathway_interaction_data.csv")

