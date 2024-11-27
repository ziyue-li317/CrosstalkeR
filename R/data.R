

#install basic packages
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db",
                       "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)
library(dplyr)

install.packages("dplyr", dependencies = TRUE)
library(dplyr)

#data input and arrange
input_data_un <-read.delim("data/GSE60450_Lactation-GenewiseCounts.txt",
                           stringsAsFactors = FALSE)
sampleinfo_un <- read.delim("data/SampleInfo.txt", stringsAsFactors = TRUE)
rownames(sampleinfo_un) <-sampleinfo_un$FileName
sampleinfo <- sampleinfo_un[, -1]


transposed_data <- t(input_data_un)
transposed_data <- as.data.frame(transposed_data)

dim(sampleinfo_un) # 12 4
dim(transposed_data) # 14 27179

fixed_rows <-transposed_data[1:2, ]
remianing_rows <- transposed_data[-(1:2), ]

dim(sampleinfo) # 12 4
dim(remianing_rows) # 12 27179

combined_rows <- merge(sampleinfo, remianing_rows,
                       by = "row.names", all = TRUE)
dim(combined_rows) # 12 27183
ordered_rows <-combined_rows %>% arrange(Status)
row_to_move <- ordered_rows[9, ]
ordered_rows <- ordered_rows[-9, ]
ordered_rows <- rbind(ordered_rows, row_to_move)

new_first_column <- c("lactate control1",
                      "lactate control2", "lactate treat1",
                      "lactate treat2", "pregnant control1",
                      "pregnant control2", "pregnant treat1",
                      "pregnant treat2", "virgin control1",
                      "virgin control2", "virgin treat1", "virgin treat2")
ordered_rows <-cbind(NewColumn = new_first_column, ordered_rows)
rownames(ordered_rows) <- ordered_rows$Row.names
ordered_rows <- ordered_rows[, -(1:2)]
ordered_rows <-cbind(NewColumn = new_first_column, ordered_rows)

rownames(ordered_rows) <- ordered_rows$NewColumn
ordered_rows <- ordered_rows[, -(1:3)]
ordered_rows <- ordered_rows[, -1]




final_data <-rbind(fixed_rows, ordered_rows)






final_data <-t(final_data)
final_data <-as.data.frame(final_data)
final_data <-final_data[, -(2)]
rownames(final_data) <- final_data[, 1]
my_data <- final_data[, -1]
# here we generate right formal for DESeq2
write.csv(my_data, file = "my_data.csv", row.names = TRUE)

