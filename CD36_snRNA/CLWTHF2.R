####################################
## Single Nuclei RNA Seq Analysis ##
####################################

##### PROJECT: WT-2

# Packages we are going to use for the analysis
library(dplyr)
library(readr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(Biobase)
library(plotly)
library(patchwork) 

###################################################################
## 1) Load a collection of 10X data-sets produced by cell ranger ##
###################################################################

# Read single-nuclei data from sample WT-2
sc <- Read10X(data.dir = "CLWTHF2/")

############################
## 2. Build Seurat object ##
############################

## We are now ready to build our Seurat R object in order to start analyzing UMI distributions of single cells.

seurat_object <- CreateSeuratObject(counts = sc, project = "CLWTHF2", min.cells = 3, min.features = 200)
dim(seurat_object)

## [1] 15758  2544

#####################################################
## 3. Quality control, Filtering and Normalization ##
#####################################################

## genes == nFeatures_RNA
## UMIs == nCount_RNA

## We will incorporate VlnPlot function (violin plot) to detect cells with unusual number of UMIs
## or those with extreme features

seurat_object %>% VlnPlot(c("nCount_RNA","nFeature_RNA"))
## These cells with high UMIs do most likely have high number of features as well. 
## For that, we will use the FeatureScatter function.
seurat_object %>% FeatureScatter(feature1 = "nCount_RNA",feature2 = "nFeature_RNA")

## Some cells seem to have extreme number of UMIs which may indicate that these barcodes actually
## captured RNA from multiple cells as these can often happen in droplet-based scRNA sequencing methods.

## On the other hand, we may also have barcodes with low number of total UMIs which may indicate that those
## barcodes may associated with (i) empty droplets with no cells captured (cell-free RNA is captured instead), 
## (ii) cells with low RNA content that damages cell type annotation and clustering.

## Additional quality control measurements

## Apoptotic (dying) cells express mitochondrial genes and export these transcripts 
## (MGI gene symbol starting with "mt-" (mouse and rat)) to the cytoplasm in mammalian cells.

## Looking for mt- genes captured in the sequencing
rownames_seurat_object <- rownames(seurat_object)
rownames_seurat_object[grepl(pattern = "mt-",x = rownames_seurat_object)]

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
## We can now visualize all three quality control measures and determine filtering thresholds.

QC_Plot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path = "../results/CLWTHF2/", filename = "CLWTHF2_QC_Violin1.png", height = 10, width = 13, units = "in", dpi = 300)
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
QC_FeaturePlot <- plot1 + plot2
ggsave(path = "../results/CLWTHF2/", filename = "CLWTHF2_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)


## We can hence remove those cells with lower than 200 genes, higher than some total of 3000 UMIs and 
## those cells with higher than total MT percentage of 10. We use the subset function to delete those cells 
## with those extreme properties.
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
seurat_object

dim(seurat_object)
## [1] 15758  1190

## Normalization
## Seurat incorporates NormalizeData function to rescale UMI counts to a global-scale by dividing each count 
## with the UMI sum of cells, then it multiplies this count with a default scale.factor of 10000.
## Finally, we log transform scaled UMI counts to normalized skew data such as RNA counts.


seurat_object <- NormalizeData(object = seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

## The resulting normalized counts are stored in a separate UMI table matrix, thus raw RNA counts are preserved
## for additional analysis.

#######################################################
## 5) Feature Selection and Dimensionality Reduction ##
#######################################################

## We first select a number of genes more variable (changing gene expression more frequently
## across cells). We use the FindVariableFeatures function to select those cells, the default value
## for the number of most variable features (or genes) is 2000.

seurat_object <- FindVariableFeatures(object = seurat_object, nfeatures = 2000,selection.method = "vst")
seurat_object

## 5.1. Cell Cycle Regression
#############################

# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html

## Cell cycle variation is a common source of uninteresting variation in single-cell RNA-seq data. To examine 
## cell cycle variation in our data, we assign each cell a score, based on its expression of G2/M and S phase 
## markers, and regress these out of the data during pre-processing.

cc_genes <- readLines(con = "../additional_tables/cell_cycle_genes_mouse.txt")
s_genes <- cc_genes[1:43]
g2m_genes <- cc_genes[44:97]
seurat_object <- CellCycleScoring(object = seurat_object, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
head(x = seurat_object@meta.data)

seurat_object@meta.data$CC.Difference <- seurat_object@meta.data$S.Score - seurat_object@meta.data$G2M.Score
colnames(seurat_object@meta.data)


## 5.2. Scaling the data 
########################

## Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to 
## dimensional reduction techniques like PCA. The ScaleData() function:
##  a) Shifts the expression of each gene, so that the mean expression across cells is 0
##  b) Scales the expression of each gene, so that the variance across cells is 1 
##    (This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)
##  c)The results of this are stored in seurat_object@assays$RNA@scale.data

## First, we remove the unwanted sources of variation:
##  we also use the ScaleData() function to remove unwanted sources of variation from a single-cell dataset.
## For example, we could 'regress out' heterogeneity associated with (for example) cell cycle stage, or 
## mitochondrial contamination.

all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, vars.to.regress=c('CC.Difference', 'percent.mt'), verbose = TRUE, features=all.genes)
seurat_object <- ScaleData(seurat_object)

## 5.3. Perform linear dimension reduction
##########################################

## Next we perform PCA on the scaled data. 
seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(object = seurat_object))

## Seurat provides several useful ways of visualizing both cells and features that define the PCA, 
## including VizDimReduction(), DimPlot(), and DimHeatmap()

print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
DimPlot(seurat_object, reduction = "pca")
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)

png('../results/CLWTHF2/CLWTHF2_PCA_heatmap.png', width = 15, height = 15, units = "in", res = 300)
DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


## 5.3..1 Determine the 'dimensionality' of the dataset.

## To overcome the extensive technical noise in any single feature for scRNA-seq data, 
## Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metafeature' that combines
## information across a correlated feature set. The top principal components therefore represent a robust compression of
## the dataset. However, how many components should we choose to include? 10? 20? 100?

## We used an heuristic method that generates an 'Elbow plot': a ranking of principle components based on the percentage of
## variance explained by each one (ElbowPlot() function)

ElbowPlot(object = seurat_object, ndims = 50)
ggsave('../results/CLWTHF2/elbow_plot.png', height = 15, 
       width = 15, dpi = 300, units = "in")

## 5.4 Run non-linear dimensional reduction (UMAP/tSNE) ##
##########################################################

## We can incorporate the method t-distributed stochastic 
## neighbor embedding (t-SNE), that reduce high dimensional datasets to 2 dimensional datasets by capturing 
## locally similar samples (or cells).

seurat_object <- RunTSNE(seurat_object, dims = 1:30)
DimPlot(seurat_object, reduction = "tsne", pt.size = 2)
ggsave(filename = "CLWTHF2_tsne01.png",path = "../results/CLWTHF2/",width = 15, height = 15, units = "in")

## We can alternatively incorporate another dimensionality reduction method, which is called 
## Uniform Manifold Approximation and Projection or UMAP.

seurat_object <- RunUMAP(object = seurat_object, dims = 1:30)
DimPlot(seurat_object, reduction = "umap", pt.size = 1)
ggsave(filename = "CLWTHF2_umap01.png",path = "../results/CLWTHF2/",width = 15, height = 15, units = "in")

## 5.5 Clustering
#################

## We use the shared nearest neighbor to cluster cells.

seurat_object <- FindNeighbors(object = seurat_object, dims = 1:30)

## We will now move on to clustering with the computed nearest neighbors. 
## We can visualize the results of multiple resolutions in the same time to determine the best number of clusters.

seurat_object <- FindClusters(object = seurat_object, resolution = c(0.2,0.5,0.7,1,1.3))

colnames(seurat_object@meta.data)

g1 <- DimPlot(object = seurat_object, reduction = "umap", label = T, group.by = "RNA_snn_res.0.2")
g2 <- DimPlot(object = seurat_object, reduction = "umap", label = T, group.by = "RNA_snn_res.0.5")
g3 <- DimPlot(object = seurat_object, reduction = "umap", label = T, group.by = "RNA_snn_res.0.7")
g4 <- DimPlot(object = seurat_object, reduction = "umap", label = T, group.by = "RNA_snn_res.1")
g5 <- DimPlot(object = seurat_object, reduction = "umap", label = T, group.by = "RNA_snn_res.1.3")

g1 + g2 + g3 + g4 + g5
ggsave(filename = "CLWTHF2_umap_resolutions.png", path = "../results/CLWTHF2/", height = 15, 
       width = 15, dpi = 300, units = "in")

## The number of cluster to use will depend on your analysis. In our case resolution of 0.7 fitted correctly

Idents(seurat_object) <- "RNA_snn_res.0.7"


## 5.6. Save the Seurat Object

saveRDS(object = seurat_object,file = '../results/CLWTHF2/CLWTHF2_Seuratv4.rds')
