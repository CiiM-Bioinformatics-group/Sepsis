########################################################################################################################################
#### SCRNASEQ SEPSIS PROJECT ##################################################################################################################
########################################################################################################################################

### Import data
setwd('path/Sepsis project')
load('12.02.2021 Sepsis normalized PBMC seurat object.RData')

### TF motif enrichment analysis
library(Matrix)
library(GENIE3)
library(readr)
library(DT)
library(reshape2)
library(visNetwork)
library(RcisTarget)

##### download the reference dataset: 
motifRankings <- importRankings("path/Sepsis project/hg38_ref/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")  #### your path to the downloaded reference dataset
data(motifAnnotations_hgnc)

### Marker genes
Mono <- subset(mydata, idents = c("Monocytes"))
Idents(Mono) <- "Immune_class"
Mono <- subset(Mono, subset = Immune_class != "Healthy control")
Mono.markers <- FindAllMarkers(Mono, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
Unclassified <- subset(Mono.markers, cluster == 'Unclassified')
Immune_paralysis <- subset(Mono.markers, cluster == 'Immune paralysis ')
MALS <- subset(Mono.markers, cluster == 'MALS')

### TF enrichement analysis
targets <- row.names(Unclassified)  ##### your gene list

motifEnrichmentTable_wGenes <- cisTarget(targets, 
                                         motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)

motifs_AUC <- calcAUC(targets, motifRankings, nCores=1)
auc <- getAUC(motifs_AUC)

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_hgnc)

motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable,
                                                   rankings=motifRankings, 
                                                   geneSets=targets)

motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]

datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf","geneSet","rankAtMax"), with=FALSE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$TF_highConf,
                            motifEnrichmentTable$geneSet),
                      function(x) {
                        genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(genes, "; ")))
                        return(genesSplit)
                      })

##### the result of TF ######
tf.list <- anotatedTfs$geneSet