
# load packages
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(dplyr)
library(DESeq2)

# data preparation
mycounts_unf <-final_data
head(mycounts_unf)
rownames(mycounts_unf) <-mycounts_unf[, 1]
mycounts_unf <-mycounts_unf[, -1]
head(mycounts_unf)

# filtering lowly expressed genes
cpm_value <- cpm(mycounts_unf)
keep_genes <- rowSums(cpm_value >=0.5) >= 7
mycounts <- mycounts_unf[keep_genes, ]

# condition
condition <- factor(c(rep("lactate_control",2), rep("lactate_treat",2),
                      rep("pregnant_control",2), rep("pregnant_treat", 2),
                      rep("virgin_control", 2), rep("virgin_treat", 2)),
                    level = c("lactate_control", "lactate_treat",
                              "pregnant_control", "pregnant_treat",
                              "virgin_control", "virgin_treat"))

# colData
colData <- data.frame(row.names = colnames(mycounts), condition)
colData

# create dds
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~ condition)
# dds normalization
dds_norm <- DESeq(dds)
# gain data after normalization
normalized_counts <- counts(dds_norm, normalized = TRUE)
head(normalized_counts)

# order data by difference in expression : mad value
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad,
                                             decreasing = T), ]

# log
rld <- rlog(dds_norm, blind = FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing = T), ]

# data output after normalization
write.table(normalized_counts, file = "dds_normalized_counts.xls",
            quote = F, sep = "\t", row.names = T, col.names = T)
write.table(rlogMat, file = "dds_normalized_counts_rlog.xls",
            quote = F, sep = "\t", row.names = T, col.names = T)


# analyse differences in expression
# lactate
res_lactate <- results(dds_norm, contrast = c("condition",
                                      "lactate_control", "lactate_treat" ))
res_lactate = res_lactate[order(res_lactate$pvalue), ]
head(res_lactate)
summary(res_lactate)
write.csv(res_lactate, file = "lactate_result.csv")

diff_gene_deseq2_lactate <- subset(res_lactate,padj < 0.05 &
                             (log2FoldChange >1 | log2FoldChange < -1))
dim(diff_gene_deseq2_lactate) # 6321*6
head(diff_gene_deseq2_lactate)
write.csv(diff_gene_deseq2_lactate, file = "diff_gene_deseq2_lactate.csv")

#annotation
library("biomaRt")
library("curl")

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
entrez_ids <- row.names(diff_gene_deseq2_lactate)

my_gene_ids <- getBM(
  attributes = c("entrezgene_id", "ensembl_gene_id"),
  filters = "entrezgene_id",
  values = entrez_ids,
  mart = mart
)

head(my_gene_ids)
my_ensembl_gene_id <- my_gene_ids[, 2]

hg_symbols <-getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name',"description" ),
  filters = "ensembl_gene_id",
  values = my_ensembl_gene_id,
  mart = mart
)
head(hg_symbols)


all_symbols <-merge(hg_symbols, my_gene_ids, by="ensembl_gene_id", all = FALSE)


entrezgene_id <- rownames(diff_gene_deseq2_lactate)
diff_gene_deseq2_lactate <-cbind(entrezgene_id, diff_gene_deseq2_lactate)
colnames(diff_gene_deseq2_lactate)[1] <- c("entrezgene_id")
diff_name <- merge(diff_gene_deseq2_lactate, all_symbols, by="entrezgene_id",
                   all = FALSE)
head(diff_gene_deseq2_lactate)
head(diff_name)












# pathway enrichment analysis
# Load additional required packages for enrichment analysis
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
install.packages("ggplot2")
library(ggplot2)

# Ensure that the entrezgene_id is in the right format
diff_gene_deseq2_lactate$entrezgene_id <- as.character(diff_gene_deseq2_lactate$entrezgene_id)

# Perform pathway enrichment analysis using KEGG
kegg_enrich <- enrichKEGG(gene = diff_gene_deseq2_lactate$entrezgene_id,
                          organism = "mmu", # mouse
                          pvalueCutoff = 0.05)

# View the enrichment results for KEGG
head(kegg_enrich)

# Save the KEGG enrichment results to a CSV file
write.csv(as.data.frame(kegg_enrich), file = "kegg_enrichment_results.csv", row.names = FALSE)

# Visualize the KEGG enrichment results
dotplot(kegg_enrich) + ggtitle("KEGG Pathway Enrichment")
ggsave("kegg_enrichment_dotplot.png")

# Perform Gene Ontology (GO) enrichment analysis using enrichGO instead of Reactome
go_enrich <- enrichGO(gene = diff_gene_deseq2_lactate$entrezgene_id,
                      OrgDb = org.Mm.eg.db,  # Mouse annotation
                      ont = "BP",  # Biological Process (can also be "CC" for Cellular Component or "MF" for Molecular Function)
                      pvalueCutoff = 0.05)

# View the GO enrichment results
head(go_enrich)

# Save the GO enrichment results to a CSV file
write.csv(as.data.frame(go_enrich), file = "go_enrichment_results.csv", row.names = FALSE)

# Visualize the GO enrichment results
dotplot(go_enrich) + ggtitle("GO Pathway Enrichment")
ggsave("go_enrichment_dotplot.png")

# Save the key results to files after enrichment analysis
write.csv(diff_gene_deseq2_lactate, file = "diff_gene_deseq2_lactate.csv", row.names = TRUE)
write.csv(as.data.frame(kegg_enrich), file = "kegg_enrichment_results.csv", row.names = FALSE)
write.csv(as.data.frame(reactome_enrich), file = "reactome_enrichment_results.csv", row.names = FALSE)


