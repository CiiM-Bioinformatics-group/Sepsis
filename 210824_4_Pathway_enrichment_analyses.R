########################################################################################################################################
#### SCRNASEQ SEPSIS PROJECT ##################################################################################################################
########################################################################################################################################
library(dplyr)
library(Seurat)
library(ggplot2)

### Import data
setwd('path/Sepsis project')
load('12.02.2021 Sepsis normalized PBMC seurat object.RData')

#######################################################################################################
### Remove Dendritic cells, platelets and undefined (not enough cells)
mydata <- subset(mydata, subset = celltype != 'Dendritic cells')
mydata <- subset(mydata, subset = celltype != 'Platelets')
mydata <- subset(mydata, subset = celltype != 'Undefined')
mydata$celltype <- Idents(mydata)

####################################################################################################
### Pathway enrichment analyses - Sepsis vs. healthy
library(ggsignif)
library(reshape2)
library(clusterProfiler)
library(enrichplot)
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(tidyverse)

Idents(mydata) <- "Disease"

DEG.up <- list()
DEG.down <- list()
DEG.all <- list()
DEG.celltypes <- data.frame(Gene=character(), Celltype=character(), P=double(), LogFC=double(), stringsAsFactors = F)
for(celltype in unique(mydata$celltype)){
  data.celltype <- FetchData(mydata, vars = "celltype")
  data.celltype <- mydata[, which(x = data.celltype == celltype)]
  DEG.celltype <- FindMarkers(data.celltype, logfc.threshold = 0, ident.1 = 'Sepsis', ident.2 = 'Healthy', min.cells.feature = 1, min.cells.group = 1)
  sig.Deg <- DEG.celltype[which(DEG.celltype$p_val_adj <= 0.05),]
  id <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(sig.Deg), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  row.names(id) <- id$SYMBOL
  up <- id[rownames(sig.Deg[which(sig.Deg$avg_logFC > 0),]),2]
  down <- id[rownames(sig.Deg[which(sig.Deg$avg_logFC < 0),]),2]
  both <- id[,2]
  DEG.up[[celltype]] <- up
  DEG.down[[celltype]] <- down
  DEG.all[[celltype]] <- both
  if(nrow(sig.Deg)>0){
    DEG <- data.frame(Gene=row.names(sig.Deg), Celltype=celltype, P=sig.Deg$p_val_adj, LogFC=sig.Deg$avg_logFC, stringsAsFactors = F)
    DEG.celltypes <- rbind(DEG.celltypes, DEG)
  }
}
KEGG.up <- compareCluster(geneClusters = DEG.up, fun = "enrichKEGG")
KEGG.down <- compareCluster(geneClusters = DEG.down, fun = "enrichKEGG")
KEGG.all <- compareCluster(geneClusters = DEG.all, fun = "enrichKEGG")
GO.up <- compareCluster(geneClusters = DEG.up, fun = "enrichGO", OrgDb = organism)
GO.down <- compareCluster(geneClusters = DEG.down, fun = "enrichGO", OrgDb = organism)
GO.all <- compareCluster(geneClusters = DEG.all, fun = "enrichGO", OrgDb = organism)

print(enrichplot::dotplot(KEGG.up) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: Up regulated DEGs"))
print(enrichplot::dotplot(KEGG.down) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: Down regulated DEGs"))
print(enrichplot::dotplot(KEGG.all) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: All DEGs"))
print(enrichplot::dotplot(GO.up)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: Up regulated DEGs"))
print(enrichplot::dotplot(GO.down)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: Down regulated DEGs"))
print(enrichplot::dotplot(GO.all)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: All DEGs"))

####################################################################################################
### Pathway analyses - MALS vs. Immune paralysis
Idents(mydata) <- "Immune_class"

DEG.up <- list()
DEG.down <- list()
DEG.all <- list()
DEG.celltypes <- data.frame(Gene=character(), Celltype=character(), P=double(), LogFC=double(), stringsAsFactors = F)
for(celltype in unique(mydata$celltype)){
  data.celltype <- FetchData(mydata, vars = "celltype")
  data.celltype <- mydata[, which(x = data.celltype == celltype)]
  DEG.celltype <- FindMarkers(data.celltype, logfc.threshold = 0, ident.1 = 'MALS', ident.2 = 'Immune paralysis ', min.cells.feature = 1, min.cells.group = 1)
  sig.Deg <- DEG.celltype[which(DEG.celltype$p_val_adj <= 0.05),]
  id <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(sig.Deg), columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  row.names(id) <- id$SYMBOL
  up <- id[rownames(sig.Deg[which(sig.Deg$avg_logFC > 0),]),2]
  down <- id[rownames(sig.Deg[which(sig.Deg$avg_logFC < 0),]),2]
  both <- id[,2]
  DEG.up[[celltype]] <- up
  DEG.down[[celltype]] <- down
  DEG.all[[celltype]] <- both
  if(nrow(sig.Deg)>0){
    DEG <- data.frame(Gene=row.names(sig.Deg), Celltype=celltype, P=sig.Deg$p_val_adj, LogFC=sig.Deg$avg_logFC, stringsAsFactors = F)
    DEG.celltypes <- rbind(DEG.celltypes, DEG)
  }
}
KEGG.up <- compareCluster(geneClusters = DEG.up, fun = "enrichKEGG")
KEGG.down <- compareCluster(geneClusters = DEG.down, fun = "enrichKEGG")
KEGG.all <- compareCluster(geneClusters = DEG.all, fun = "enrichKEGG")
GO.up <- compareCluster(geneClusters = DEG.up, fun = "enrichGO", OrgDb = organism)
GO.down <- compareCluster(geneClusters = DEG.down, fun = "enrichGO", OrgDb = organism)
GO.all <- compareCluster(geneClusters = DEG.all, fun = "enrichGO", OrgDb = organism)

print(enrichplot::dotplot(KEGG.up) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: Up regulated DEGs"))
print(enrichplot::dotplot(KEGG.down) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: Down regulated DEGs"))
print(enrichplot::dotplot(KEGG.all) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + ggtitle("KEGG pathways: All DEGs"))
print(enrichplot::dotplot(GO.up)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: Up regulated DEGs"))
print(enrichplot::dotplot(GO.down)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: Down regulated DEGs"))
print(enrichplot::dotplot(GO.all)+ theme(axis.text.x = element_text(angle = 45, hjust=1))+ ggtitle("GO terms: All DEGs"))