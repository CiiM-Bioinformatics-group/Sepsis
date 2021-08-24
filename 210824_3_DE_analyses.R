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

#####################################################################################################
significance <- function(x, treshold = 1, pval = 0.05){
  x$sig <- "nonsig"
  x[x$p_val_adj <= pval, "sig"] <- "sig"
  x <- data.frame(gene = rownames(x), x, stringsAsFactors = FALSE)
  return(x)
}  

# Volcano Plot Function
volcano <- function(x, title, top = 10){
  x <- x[order(abs(x$p_val_adj), decreasing = FALSE),]
  volcano.labels <- ifelse(x$sig %in% c("sig"), as.character(x$gene), "")
  volcano.labels[(top+1):length(volcano.labels)] <- ""
  
  ggplot(x, aes(x = avg_logFC, y = -log(p_val_adj, 10), label = gene)) +
    geom_point(aes(color = sig), size = 6) +
    geom_label_repel(aes(label = volcano.labels), 
                     hjust = 1.25, vjust = 0, size = 8) +
    scale_colour_manual(name = "Significance", 
                        values = c("sig" = "black", "nonsig" = "gray"),
                        breaks = c("sig", "nonsig"),
                        labels = c("sig", "nonsig")) +
    ggtitle(label = title) +
    ylab("-log10(p_val_adj)") +
    xlab("avg_logFC") +
    geom_hline(yintercept = -log(0.05,10)) +
    geom_vline(xintercept = c(0)) +
    theme(axis.text = element_text(size = 30),
          axis.title = element_text(size = 30),
          title = element_text(size = 30),
          legend.text = element_text(size = 30))
}

#####################################################################################################
### DE analyses Sepsis vs HC
mydata$celltype.disease <- paste(Idents(mydata), mydata$Disease, sep = "_")
Idents(mydata) <- "celltype.disease"

mono_DE_dis <- FindMarkers(mydata, ident.1 = "Monocytes_Sepsis", ident.2 = "Monocytes_Healthy", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(mono_DE_dis)
volcano(x = volcano.data, title = "DEGs between sepsis and control - monocytes")

CD4_DE_dis <- FindMarkers(mydata, ident.1 = "CD4 T cells_Sepsis", ident.2 = "CD4 T cells_Healthy", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(CD4_DE_dis)
volcano(x = volcano.data, title = "DEGs between sepsis and control - CD4 T cells")

CD8_DE_dis <- FindMarkers(mydata, ident.1 = "CD8 T cells_Sepsis", ident.2 = "CD8 T cells_Healthy", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(CD8_DE_dis)
volcano(x = volcano.data, title = "DEGs between sepsis and control - CD8 T cells")

NK_DE_dis <- FindMarkers(mydata, ident.1 = "NK cells_Sepsis", ident.2 = "NK cells_Healthy", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(NK_DE_dis)
volcano(x = volcano.data, title = "DEGs between sepsis and control - NK cells")

B_DE_dis <- FindMarkers(mydata, ident.1 = "B cells_Sepsis", ident.2 = "B cells_Healthy", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(B_DE_dis)
volcano(x = volcano.data, title = "DEGs between sepsis and control - B cells")

#####################################################################################################
### DE analyses MALS vs Immune paralysis
mydata$celltype.immune.class <- paste(mydata$celltype, mydata$Immune_class, sep = "_")
Idents(mydata) <- "celltype.immune.class"

mono_DE_imm <- FindMarkers(mydata, ident.1 = "Monocytes_MALS", ident.2 = "Monocytes_Immune paralysis ", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(mono_DE_imm)
volcano(x = volcano.data, title = "DEGs between MALS and immune paralysis - monocytes")

CD4_DE_imm <- FindMarkers(mydata, ident.1 = "CD4 T cells_MALS", ident.2 = "CD4 T cells_Immune paralysis ", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(CD4_DE_imm)
volcano(x = volcano.data, title = "DEGs between MALS and immune paralysis - CD4 T cells")

B_DE_imm <- FindMarkers(mydata, ident.1 = "B cells_MALS", ident.2 = "B cells_Immune paralysis ", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(B_DE_imm)
volcano(x = volcano.data, title = "DEGs between MALS and immune paralysis - B cells")

NK_DE_imm <- FindMarkers(mydata, ident.1 = "NK cells_MALS", ident.2 = "NK cells_Immune paralysis ", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(NK_DE_imm)
volcano(x = volcano.data, title = "DEGs between MALS and immune paralysis - NK cells")

CD8_DE_imm <- FindMarkers(mydata, ident.1 = "CD8 T cells_MALS", ident.2 = "CD8 T cells_Immune paralysis ", verbose = FALSE, logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1)
volcano.data <- significance(CD8_DE_imm)
volcano(x = volcano.data, title = "DEGs between MALS and immune paralysis - CD8 T cells")