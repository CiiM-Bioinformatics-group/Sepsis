########################################################################################################################################
#### SCRNASEQ SEPSIS PROJECT ##################################################################################################################
########################################################################################################################################
library(dplyr)
library(Seurat)
library(ggplot2)

### Import data
setwd('path/Sepsis project')
load('12.02.2021 Sepsis normalized PBMC seurat object.RData')

### Select monocytes only
Mono <- subset(mydata, idents = c("Monocytes"))

### marker genes enrichment analysis
BiocManager::install("AUCell")
require(AUCell)
Mono$cluster.factor <- factor(Mono$seurat_clusters)

###calculate cell rankings#####
tmp<-as.data.frame(as.matrix(Mono@assays$RNA@data))
filtered_matrix<-rowSums(tmp > 0)
filtered_matrix<-filtered_matrix[filtered_matrix>10]
tmp <- tmp[row.names(tmp) %in% names(filtered_matrix),]
gc()
memory.limit()
memory.limit(20000)
cells_rankings <- AUCell_buildRankings(as.matrix(tmp))

#######  read a list of markers from file #####
markers <- openxlsx::read.xlsx(xlsxFile = 'your file.xlsx')
markers_list<-as.list(markers)

cells_AUC <- AUCell_calcAUC(markers_list, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)
tmp <- getAUC(cells_AUC)
tmp <- as.data.frame(t(tmp))

df <- data.frame(row.names = row.names(Mono@meta.data), cluster = Mono$cluster.factor, stringsAsFactors = F)

df <- as.data.frame(merge(df, tmp, by = 0))
df$geneset<-df[,"X1"]   ## this can be cell type name when input a list
df %>% group_by(cluster) %>% summarize(Mean = mean(geneset, na.rm=TRUE)) -> df.tmp
df.tmp$Mean <- scale(df.tmp$Mean)
df <- merge(df, df.tmp, key = "cluster", all = T)

df <- df[grep("mix", df$cluster, invert = T),]
#df$geneset<-df[,i]
dev.off()
gg <- ggplot(df, aes(y = geneset, x = cluster, group = cluster, fill = Mean)) + 
  geom_violin(scale = "area", trim = T, adjust = 3, draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30")+
  scale_fill_gradient2(low="blue", mid="white",high="red") +
  ylab("AUC score") + xlab("cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3))#+ggtitle(paste(i))

gg