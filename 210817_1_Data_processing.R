### Loading packages
library(Seurat)
library(dplyr)
library(openxlsx)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(SingleR)
library(Matrix)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(DirichletReg)

### Import 10X data
pool2 <- Read10X(data.dir = "path/Data/Pool2")
pool3 <- Read10X(data.dir = "path/Data/Pool3")
pool4 <- Read10X(data.dir = "path/Data/Pool4")
pool5 <- Read10X(data.dir = "path/Data/Pool5")
pool6 <- Read10X(data.dir = "path/Data/Pool6")

### Import demultiplexed data
clusters2 <- read.table('path/Data/pool_2.demultiplex_results.txt', header = TRUE, stringsAsFactors = FALSE)
clusters3 <- read.table('path/Data/pool_3.demultiplex_results.txt', header = TRUE, stringsAsFactors = FALSE)
clusters4 <- read.table('path/Data/pool_4.demultiplex_results.txt', header = TRUE, stringsAsFactors = FALSE)
clusters5 <- read.table('path/Data/pool_5.demultiplex_results.txt', header = TRUE, stringsAsFactors = FALSE)
clusters6 <- read.table('path/Data/pool_6.demultiplex_results.txt', header = TRUE, stringsAsFactors = FALSE)

### Merge demultiplexed data
clusters2$pool <- 'pool 2'
clusters3$pool <- 'pool 3'
clusters4$pool <- 'pool 4'
clusters5$pool <- 'pool 5'
clusters6$pool <- 'pool 6'
clusters <- rbind(clusters2, clusters3, clusters4, clusters5, clusters6)
clusters <- subset(clusters, status == "Singlet" | status == "singlet")

### Remove all data that are not in singlets files
pool2 <-	pool2[, colnames(pool2) %in% clusters$barcode]
pool3 <-	pool3[, colnames(pool3) %in% clusters$barcode]
pool4 <-	pool4[, colnames(pool4) %in% clusters$barcode]
pool5 <-	pool5[, colnames(pool5) %in% clusters$barcode]
pool6 <-	pool6[, colnames(pool6) %in% clusters$barcode]

### Adding the cluster assignment	and metadata
cluster.assignment <- clusters$assignment
names(cluster.assignment) <- clusters$barcode
metadata <- openxlsx::read.xlsx(xlsxFile = 'path/Data/Demographics.xlsx', colNames = TRUE, rowNames = FALSE)
clusters <- merge(clusters, metadata, by.x = 2, by.y = 1, all.x = T)
clusters <- subset(clusters, !duplicated(clusters$barcode))
rownames(clusters) <- clusters$barcode

### Adding the metadata and create seurat object and merge
pool2 <- CreateSeuratObject(pool2, min.cells = 2, meta.data = clusters, project = "Pool 2")
pool3 <- CreateSeuratObject(pool3, min.cells = 2, meta.data = clusters, project = 'Pool 3')
pool4 <- CreateSeuratObject(pool4, min.cells = 2, meta.data = clusters, project = 'Pool 4')
pool5 <- CreateSeuratObject(pool5, min.cells = 2, meta.data = clusters, project = 'Pool 5')
pool6 <- CreateSeuratObject(pool6, min.cells = 2, meta.data = clusters, project = 'Pool 6')
mydata23 <- merge(pool2, y = pool3, project = 'sepsis')
mydata45 <- merge(pool4, y = pool5, project = 'sepsis')
mydata456 <- merge(mydata45, y = pool6, project = 'sepsis')
mydata <- merge(mydata23, y = mydata456, project = 'sepsis')
rm(mydata23, mydata45, mydata456, metadata, pool2, pool3, pool4, pool5, pool6)

### Remove sepsis_13 sample (clear outlier on UMAP)
mydata <- subset(mydata, subset = assignment != 'Sepsis_13')

###########################################################################################

### Calculate mitochondrial percentage
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")

### Remove mitochondrial genes
counts <- GetAssayData(mydata, assay = "RNA")
counts <- counts[-grep('^MT-',rownames(counts)),]  
mydata <- subset(mydata, features = rownames(counts))

### Calculate ribosomal percentage
mydata[["percent.rb"]] <- PercentageFeatureSet(mydata, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

### Remove ribosomal genes
counts <- GetAssayData(mydata, assay = "RNA")
counts <- counts[-grep('^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA',rownames(counts)),]  
mydata <- subset(mydata, features = rownames(counts))

### Filter based on quality measures
mydata <- subset(mydata, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)

### Normalize data using global-scaling normalization method and scale data
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(mydata)
rm(counts, clusters, cluster.assignment)
mydata <- ScaleData(mydata, features = all.genes, vars.to.regress = c("percent.mt"))

###########################################################################################

### Dimensionality reduction
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
mydata <- RunPCA(mydata, verbose = FALSE)
ElbowPlot(mydata, ndims = 40)
mydata <- FindNeighbors(mydata, reduction = "pca", dims = 1:20)
mydata <- FindClusters(mydata, resolution = 0.4)
mydata <- RunUMAP(mydata, dims = 1:20)
DimPlot(mydata, reduction = "umap", group.by = 'seurat_clusters')

###########################################################################################
### Visualize marker genes

FeaturePlot(mydata, features = c("MS4A1", 'CD79A'))
FeaturePlot(mydata, features = c("GNLY", 'NKG7', 'GZMB'))
FeaturePlot(mydata, features = c("CD8A", 'GZMK'))
FeaturePlot(mydata, features = c("PPBP"))
FeaturePlot(mydata, features = c('IL7R', 'CCR7', 'TCF7', 'S100A4'))
FeaturePlot(mydata, features = c("CD14", "LYZ", "FCGR3A", 'MS4A7'))
FeaturePlot(mydata, features = c("FCER1A", 'CST3', 'CD86'))

VlnPlot(mydata, features = c("MS4A1", "CD79A"))
VlnPlot(mydata, features = c("GNLY", 'NKG7', 'GZMB'))
VlnPlot(mydata, features = c("CD8A", 'GZMK'))
VlnPlot(mydata, features = c("PPBP"))
VlnPlot(mydata, features = c('IL7R', 'CCR7', 'TCF7', 'S100A4'), ncol =2)
VlnPlot(mydata, features = c("CD14", "LYZ", "FCGR3A", 'MS4A7'), ncol = 2)
VlnPlot(mydata, features = c("FCER1A", 'CST3', 'CD86'))

###################################################################################################
### Annotate cell types to clusters
new.cluster.ids <- c('CD4 T cells', 'Monocytes','CD8 T cells','Monocytes', 'B cells', 'NK cells',
                     'Monocytes','Undefined','CD4 T cells','Monocytes','CD4 T cells','Undefined',
                     'Platelets','Dendritic cells')

names(new.cluster.ids) <- levels(mydata)
mydata <- RenameIdents(mydata, new.cluster.ids)
markers_pbmc <- c("IL7R","CCR7","TCF7", 'S100A4', "CD8A","GZMK", 'GNLY', "NKG7","GZMB","CD79A","MS4A1","CD14","LYZ","FCGR3A",'MS4A7','FCER1A','CST3','CD86',"PPBP")
mydata$celltype <- Idents(mydata)

### Visualize marker genes in dotplot
DotPlot(mydata,features=markers_pbmc,cols="RdBu")+coord_flip() +
theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
ggtitle(label = "PBMC Marker Genes")

### Create R object for further analyses
save(mydata, file = "12.02.2021 Sepsis normalized PBMC seurat object.RData")