##################################
## Single cell RNA Seq Analysis ##
##################################

## Sample Cirrhotic 4

## Packages we are going to use for the analysis
library(Seurat)
library(ggplot2)

#######################################################################
## 1) Reading human liver scRNA-seq dataset from patient Cirrhotic 4 ##
#######################################################################

## Loading data from Cirrhotic 4 which is a single file. The data 
## is available in the GEO database: GSE136103.

sc <- Read10X("../data/Cirrhotic_4/")

############################
## 2. Build Seurat object ##
############################

## We are now ready to build our Seurat R object in order to start analyzing UMI distributions of single cells.

seurat_object <- CreateSeuratObject(counts = sc, project = "Cirrhotic_4", min.cells = 10, min.features = 200)
dim(seurat_object)
# [1] 14800  4819

######################################
## 3. Quality control and filtering ##
######################################

## The column percent.mt is added to the object metadata with the function PercentageFeatureSet(),
## which calculates the percentage of counts from a set of features. It contains 
## the mitochondrial QC metrics.
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

## We can now visualize all three quality control measures and determine filtering thresholds.
## Here, we are using the function FeatureScatter to display feature-feature relationship.
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(path = "../results/Cirrhotic_4", filename = "Cirrhotic_4_QC_Feature_plot.png", height = 10, width = 13, units = "in", dpi = 300)

#################################
## 4. Sub-setting the data set ##
#################################

## We can hence remove those cells with lower than 200 genes, higher than some total of 4000 UMIs and 
## those cells with higher than total MT percentage of 10. We use the subset function to delete those cells 
## with those extreme properties.
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

#We can see the reduction in dimension after sub-setting
dim(seurat_object)
## [1] 14800  4776

saveRDS(object = seurat_object,file = '../results/Cirrhotic_4/Cirrhotic_4_data_combined.rds')
