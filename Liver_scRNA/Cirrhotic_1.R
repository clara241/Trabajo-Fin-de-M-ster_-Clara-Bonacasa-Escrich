##################################
## Single cell RNA Seq Analysis ##
##################################

## Sample Cirrhotic 1:

## Packages we are going to use for the analysis
library(Seurat)
library(ggplot2) 

#######################################################################
## 1) Reading human liver scRNA-seq dataset from patient Cirrhotic 1 ##
#######################################################################

## Loading data from Cirrhotic 1 which is subset into 3 different files. The data 
## is available in the GEO database: GSE136103.

sc_CD45pos <- Read10X("../data/Cirrhotic_1/CD45+/")
sc_CD45negA <- Read10X("../data/Cirrhotic_1/CD45-A/")
sc_CD45negB <- Read10X("../data/Cirrhotic_1/CD45-B/")

############################
## 2. Build Seurat object ##
############################
## We are now ready to build our Seurat R objects in order to start analyzing UMI distributions of single cells.

seurat_object_pos <- CreateSeuratObject(counts = sc_CD45pos, project = "Cirrhotic_1", min.cells = 10, min.features = 200)
dim(seurat_object_pos)
# [1] 15008  2344

seurat_object_negA <- CreateSeuratObject(counts = sc_CD45negA, project = "Cirrhotic_1", min.cells = 10, min.features = 200)
dim(seurat_object_negA)
# [1] 15516  1652

seurat_object_negB <- CreateSeuratObject(counts = sc_CD45negB, project = "Cirrhotic_1", min.cells = 10, min.features = 200)
dim(seurat_object_negB)
# [1] 15661  2004

######################################
## 3. Quality control and filtering ##
######################################

## genes == nFeatures_RNA
## UMIs == nCount_RNA

## Some cells seem to have extreme number of UMIs which may indicate that these barcodes actually
## captured RNA from multiple cells as these can often happen in droplet-based scRNA sequencing methods.

## On the other hand, we may also have barcodes with low number of total UMIs which may indicate that those
## barcodes may associated with (i) empty droplets with no cells captured (cell-free RNA is captured instead), 
## (ii) cells with low RNA content that damages cell type annotation and clustering.

## Apoptotic (dying) cells express mitochondrial genes and export these transcripts 
## (MGI gene symbol starting with "mt-" (mouse and rat)) to the cytoplasm in mammalian cells.
## Looking for MT- genes captured in the sequencing
rownames_seurat_object <- rownames(seurat_object)
rownames_seurat_object[grepl(pattern = "mt-",x = rownames_seurat_object)]

## The column percent.mt is added to the object metadata with the function PercentageFeatureSet(),
## which calculates the percentage of counts from a set of features. It contains 
## the mitochondrial QC metrics.
seurat_object_pos[["percent.mt"]] <- PercentageFeatureSet(seurat_object_pos, pattern = "^MT-")
seurat_object_negA[["percent.mt"]] <- PercentageFeatureSet(seurat_object_negA, pattern = "^MT-")
seurat_object_negB[["percent.mt"]] <- PercentageFeatureSet(seurat_object_negB, pattern = "^MT-")

## We can now visualize all three quality control measures and determine filtering thresholds.
## Here, we are using the function FeatureScatter to display feature-feature relationship.
plot1 <- FeatureScatter(seurat_object_pos, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object_pos, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Cirrhotic_1", filename = "Cirrhotic_1_pos_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

plot1 <- FeatureScatter(seurat_object_negA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object_negA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Cirrhotic_1", filename = "Cirrhotic_1_negA_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

plot1 <- FeatureScatter(seurat_object_negB, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object_negB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Cirrhotic_1", filename = "Cirrhotic_1_negB_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

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
## [1] 15008  2065

seurat_object_negA <- subset(seurat_object_negA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 12.5)
dim(seurat_object_negA)
## [1] 15516  1393

seurat_object_negB <- subset(seurat_object_negB, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 12.5)
dim(seurat_object_negB)
## [1] 15661  1542

######################################
## 5. Analysis of unintegrated data ##
######################################

## We prepare the three seurat objects before merging them into a single seurat object
## by adding information in the metadata about its origin
CD45pos <- AddMetaData(object = seurat_object_pos, metadata = "CD45+", col.name = "cell_population")
CD45negA <- AddMetaData(object = seurat_object_negA, metadata = "CD45-A", col.name = "cell_population")
CD45negB <- AddMetaData(object = seurat_object_negB, metadata = "CD45-B", col.name = "cell_population")

# Merge Seurat objects adding to each of them an identification/sample name

data_combined <-  merge(x = CD45pos, y = c(CD45negA,CD45negB), 
                        add.cell.ids = c("Cirrhotic_cd45pos", "Cirrhotic_cd45nega","Cirrhotic_cd45negb"))


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

## We used an heuristic method that generates an 'Elbow plot': a ranking of principle components based on the percentage of
## variance explained by each one (ElbowPlot() function)

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

saveRDS(object = data_combined, file = '../results/Cirrhotic_1/Cirrhotic_1_data_combined.rds')

