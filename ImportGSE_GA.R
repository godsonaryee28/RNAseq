#### To start go to session -> set working directory -> to source file location #####
#Set working directory#

#### Major code provided by https://satijalab.org/seurat/ ##### Seurat V3 ####
#### Pre-processing workflow modifications created for specific data (MGH36, MGH53 and MGH54 tumors) ####
#Load required libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

# Loading GES data provided from 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/#SD1'
GSE <- read.table(file = 'GSE70630.txt', 
                  as.is = TRUE, 
                  header = TRUE
)
#### Using his a a standard pre-processing workflow ####
# Extract MGH36, MGH53 and MGH54 tumors
coln = colnames(GSE)
MGH36_vals = coln[grepl( "MGH36" , coln )]
MGH53_vals = coln[grepl( "MGH53" , coln )]
MGH54_vals = coln[grepl( "MGH54" , coln )]
combined = c(MGH36_vals, MGH54_vals, MGH53_vals)
three_mgh <- GSE[, combined]
View(combined)
View(three_mgh)
# Create a Seurat object
gse.object <- CreateSeuratObject(counts = three_mgh, project = "cancertumor")
#get mito percent, ribosomal proteins and Visualize QC metrics as a violin plot, featureScatter 
gse.object[["percent.mt"]] <- PercentageFeatureSet(gse.object, pattern = "^MT-")
gse.object[["percent.rb"]] <- PercentageFeatureSet(gse.object, pattern = "^RP[SL]")

# Add number of genes per UMI for each cell
#nUMI = nCount_RNA, nGene = nFeature_RNA
gse.object[["GenesPerUMI"]] <- log10(gse.object$nFeature_RNA) / log10(gse.object$nCount_RNA)
# Add cell IDs to metadata
#gse.object$cells <- rownames(gse.object@meta.data)
#head(gse.object)

#filter 
#seu.filtered <- subset(gse.object, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & nCount_RNA >800 & percent.rb <5)
#plot VlnPlot
VlnPlot(gse.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.1) & theme(plot.title = element_text(size = 10))
VlnPlot(gse.object, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
plot1 <- FeatureScatter(gse.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gse.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(gse.object, feature1 =  "nCount_RNA", feature2 = "percent.rb")
plot4 <- FeatureScatter(gse.object, feature1 =  "nCount_RNA", feature2 = "GenesPerUMI")
plot1 + plot2 +plot3+ plot4
## Normalizing the data ##
#using scaling normalization method “LogNormalize” 
gse.object <- NormalizeData(gse.object, normalization.method = "LogNormalize", scale.factor = 10000)
gse.object <- FindVariableFeatures(gse.object, selection.method = "vst", nfeatures = 2000)

# Identify the 10 & 20 most highly variable genes
top10 <- head(VariableFeatures(gse.object), 10)
top20<-head(VariableFeatures(gse.object), 20)
#Identify the 10 least variable genes
down10<-tail(VariableFeatures(gse.object),10)
# plot variable features with and without labels
plot5 <- VariableFeaturePlot(gse.object)
plot6 <- LabelPoints(plot = plot5, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot7 <- LabelPoints(plot = plot5, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
plot8 <- LabelPoints(plot = plot5, points = down10, repel = TRUE, xnudge = 0, ynudge = 0)
plot5 + plot7
plot5 +plot8
plot5 + plot6 +plot7+ plot8
##Scaling the data: either using SCTransform() or ScaleData()
#The ScaleData() function:Shifts the expression of each gene, so that the mean expression across cells is 0
all.genes <- rownames(gse.object)
gse.object<- ScaleData(gse.object, features = all.genes)
gse.object <- ScaleData(gse.object)
gse.object <- ScaleData(gse.object, vars.to.regress = "percent.mt")
gse.object <- ScaleData(gse.object, vars.to.regress = "percent.rb")
gse.object <- ScaleData(gse.object, vars.to.regress = "GenesPerUMI")
gse.object <- RunPCA(gse.object, features = VariableFeatures(object = gse.object))
VizDimLoadings(gse.object, dims = 1:2, reduction = "pca")
DimPlot(gse.object, reduction = "pca")
#DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, 
DimHeatmap(gse.object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(gse.object, dims = 1:15, cells = 500, balanced = TRUE)
#Determine the ‘dimensionality’ of the dataset#
gse.object <- JackStraw(gse.object, num.replicate = 100)
gse.object <- ScoreJackStraw(gse.object, dims = 1:20)
JackStrawPlot(gse.object, dims = 1:20)
# reduce computation time by using ElbowPlot()
ElbowPlot(gse.object)
## Perform Clustering of cells
gse.object <- FindNeighbors(gse.object, dims = 1:20)
gse.object <- FindClusters(gse.object, resolution = 0.5)
# Checkout the cluster IDs of the first 10 cells
head(Idents(gse.object), 10)
##Run non-linear dimensional reduction (UMAP/tSNE)
gse.object <- RunUMAP(gse.object, dims = 1:10)
DimPlot(gse.object, reduction = "umap")
#Saveout the output
saveRDS(gse.object, file = "/Volumes/GoogleDrive/My Drive/GSE70630.rds")

#### Finding differentialy expressed features (cluster biomarkers) ####
# find all markers of cluster 2
cluster2.markers <- FindMarkers(gse.object, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)
# find all markers distinguishing cluster 1 from clusters 3 and 5
cluster1.markers <- FindMarkers(gse.object, ident.1 = 1, ident.2 = c(3, 5), min.pct = 0.25)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(gse.object, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster51.markers, n = 10)
head(cluster5.markers, n = 10)
# find markers for every cluster compared to all remaining cells, report only the positive ones
gse.object.markers <- FindAllMarkers(gse.object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gse.object.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
#Seurat has several tests for differential expression which can be set with the test.use parameter  
cluster0.markers <- FindMarkers(gse.object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
print (cluster0.markers)
VlnPlot(gse.object, features = c("FOS", "CCL3"))
RidgePlot(gse.object, features = c("FOS", "CCL3"))

#plots of counts#
VlnPlot(gse.object, features = c("FOS", "CCL3"), slot = "counts", log = TRUE)
#Featureplots of the top10 genes
FeaturePlot(gse.object, features = c("FOS","CCL3","CCL4","RGS1", "CCL3L1","CCL4L1","CCL4L2","CD74","SPP1","CCL3L3"))
#DoHeatmap() generates an expression heatmap for given cells and features. 
gse.object.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(gse.object, features = top10$gene) + NoLegend()
#Saveout the output
saveRDS(gse.object, file = "/Volumes/GoogleDrive/My Drive/GSE70630final.rds")

