##################################
## Single cell RNA Seq Analysis ##
##################################

## Sample Healthy 4

## Packages we are going to use for the analysis
library(Seurat)
library(ggplot2)

#####################################################################
## 1) Reading human liver scRNA-seq dataset from patient Healthy 4 ##
#####################################################################

## Loading data from Healthy 4 which is subset into 3 different files. The data 
## is available in the GEO database: GSE136103.

sc_CD45pos <- Read10X("../data/Healthy_4/CD45+/")
sc_CD45neg <- Read10X("../data/Healthy_4/CD45-/")

############################
## 2. Build Seurat object ##
############################

## We are now ready to build our Seurat R objects in order to start analyzing UMI distributions of single cells.

seurat_object_pos <- CreateSeuratObject(counts = sc_CD45pos, project = "Healthy_4", min.cells = 10, min.features = 200)
dim(seurat_object_pos)
# [1] 14501  4932

seurat_object_neg <- CreateSeuratObject(counts = sc_CD45neg, project = "Healthy_4", min.cells = 10, min.features = 200)
dim(seurat_object_neg)
# [1] 15502  3306

######################################
## 3. Quality control and filtering ##
######################################

## The column percent.mt is added to the object metadata with the function PercentageFeatureSet(),
## which calculates the percentage of counts from a set of features. It contains 
## the mitochondrial QC metrics.
seurat_object_pos[["percent.mt"]] <- PercentageFeatureSet(seurat_object_pos, pattern = "^MT-")
seurat_object_neg[["percent.mt"]] <- PercentageFeatureSet(seurat_object_neg, pattern = "^MT-")

## We can now visualize all three quality control measures and determine filtering thresholds.
## Here, we are using the function FeatureScatter to display feature-feature relationship.
plot1 <- FeatureScatter(seurat_object_pos, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object_pos, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Healthy_4", filename = "Healthy_4_pos_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

plot1 <- FeatureScatter(seurat_object_neg, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object_neg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Healthy_4", filename = "Healthy_4_neg_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

##################################
## 4. Sub-setting the data sets ##
##################################

## We can hence remove those cells with lower than 200 genes. We also filter them 
## by the total amount of UMIs and the total MT percentage according to Quality Control
## measures from each sample. We use the subset function to delete those cells 
## that pass the thresholds chosen for each sample.
seurat_object_pos <- subset(seurat_object_pos, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10)
#We can see the reduction in dimension after sub-setting
dim(seurat_object_pos)
## [1] 14501  4835

seurat_object_neg <- subset(seurat_object_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 12.5)
dim(seurat_object_neg)
## [1] 15502  3191

######################################
## 5. Analysis of unintegrated data ##
######################################

## We prepare the three seurat objects before merging them into a single seurat object
## by adding information in the metadata about its origin
CD45pos <- AddMetaData(object = seurat_object_pos, metadata = "CD45+", col.name = "cell_population")
CD45neg <- AddMetaData(object = seurat_object_neg, metadata = "CD45-", col.name = "cell_population")

# Merge Seurat objects adding to each of them an identification/sample name
data_combined <-  merge(x = CD45pos, y = CD45neg, 
                        add.cell.ids = c("Healthy_cd45pos", "Healthy_cd45neg"))

data_combined[["RNA"]] <- split(data_combined[["RNA"]], f = data_combined$cell_population)

## We run standard analysis workflow without integration

data_combined <- NormalizeData(data_combined)
data_combined <- FindVariableFeatures(data_combined)
data_combined <- ScaleData(data_combined)
data_combined <- RunPCA(data_combined)

# Visualization of unintegrated samples
DimPlot(data_combined)

####################################
## 6. Perform integrated analysis ##
####################################
data_combined <- IntegrateLayers(object = data_combined, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)

# re-join layers after integration
data_combined[["RNA"]] <- JoinLayers(data_combined[["RNA"]])

## We are going to determine the 'dimensionality' of the dataset. Here, we use an
## heuristic method that generates a ranking with the top principal components, from PCA,
## ordered according to the percentage of variance explained by each one.

ElbowPlot(data_combined, ndims = 50)

data_combined <- FindNeighbors(data_combined, reduction = "integrated.cca", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.5)

## Run non-linear dimensional reduction (UMAP/tSNE)

data_combined <- RunUMAP(data_combined, dims = 1:30, reduction = "integrated.cca")
data_combined <- RunTSNE(data_combined, dims = 1:30, reduction = "integrated.cca")

### Visualization of non-linear dimensional reduction. the resulting clusters are both
### defined by sample origin and cell-type (unbiased clusters)

DimPlot(data_combined, reduction = "umap", group.by = "cell_population")
DimPlot(data_combined, reduction = "tsne", group.by = "cell_population")

## Save the combined seurat object

saveRDS(object = data_combined,file = '../results/Healthy_4/Healthy_4_data_combined.rds')

