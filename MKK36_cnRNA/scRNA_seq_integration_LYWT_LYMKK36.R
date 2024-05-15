###########################
## scRNA-seq INTEGRATION ##
###########################

##### PROJECT: Integration MKK3/6 data

# 0) Load all the required Libraries
library(Seurat)
library(readxl)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(Nebulosa)
library(dplyr)
library(Biobase)
library(plotly)
library(cowplot)
library(SingleR)
library(scater)
library(DESeq2)
library(DescTools)
library(BiocParallel)
library(grid)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(data.table)
library(readr)
library(tibble)
library(stringr)

###############################
## 1)Setup the Seurat object ##
###############################

# 1.1 Read Gene expression matrices for all populations

LYMKK361 <- readRDS(file = "../results/LYMKK361/LYMKK361_Seuratv4.rds")
LYMKK362 <- readRDS(file = "../results/LYMKK362/LYMKK362_Seuratv4.rds")
LYWT1 <- readRDS(file = "../results/LYWT1/LYWT1_Seuratv4.rds")
LYWT2 <- readRDS(file = "../results/LYWT2/LYWT2_Seuratv4.rds")

# 1.2 Add extra column to the Seurat objects with the condition

LYWT1 <- AddMetaData(object = LYWT1, metadata = "Control", col.name = "Condition")
LYWT2 <- AddMetaData(object = LYWT2, metadata = "Control", col.name = "Condition")
LYMKK361 <- AddMetaData(object = LYMKK361, metadata = "KO", col.name = "Condition")
LYMKK362 <- AddMetaData(object = LYMKK362, metadata = "KO", col.name = "Condition")


# 1.3 Merge Seurat objects

merged_seurat <- merge(x = LYWT1, y = c(LYWT2, LYMKK361, LYMKK362), 
                       add.cell.ids = c("LYWT1", "LYWT2","LYMKK361","LYMKK362"))
merged_seurat
head(colnames(merged_seurat))
tail(colnames(merged_seurat ))
dim(merged_seurat)
# [1] 17538 10111

# 1.4 split the dataset into a list of two seurat objects

list_seurat <- SplitObject(merged_seurat, split.by = "orig.ident")

### Let us remove the objects that we do not need any more to save memory

rm(LYMKK361,LYMKK362, LYWT1, LYWT2)

############################
## 2) Perform Integration ##
############################

data_combined <- merged_seurat

data_combined[["RNA"]] <- split(data_combined[["RNA"]], f = data_combined$Condition)

## We run standard analysis workflow without integration

data_combined <- NormalizeData(data_combined)
data_combined <- FindVariableFeatures(data_combined)
data_combined <- ScaleData(data_combined)
data_combined <- RunPCA(data_combined)

DimPlot(data_combined)


ElbowPlot(data_combined, ndims = 50, reduction = "pca")

data_combined <- FindNeighbors(data_combined, dims = 1:30, reduction = "pca")
data_combined <- FindClusters(data_combined, resolution = 0.5, cluster.name = "unintegrated_clusters")

DimPlot(data_combined)

data_combined <- RunUMAP(data_combined, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(data_combined, reduction = "umap.unintegrated", group.by = c("Condition", "seurat_clusters"))
DimPlot(data_combined, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))


## Integration

data_combined <- IntegrateLayers(object = data_combined, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)


####################################
## 3) Perform integrated analysis ##
####################################


# re-join layers after integration
data_combined[["RNA"]] <- JoinLayers(data_combined[["RNA"]])

ElbowPlot(data_combined, ndims = 50, reduction = "pca")


data_combined <- FindNeighbors(data_combined, reduction = "integrated.cca", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.5)

DimPlot(data_combined)

data_combined <- RunUMAP(data_combined, dims = 1:30, reduction = "integrated.cca")
data_combined <- RunTSNE(data_combined, dims = 1:30, reduction = "integrated.cca")

#  Visualization

## Combined by groups
DimPlot(data_combined, reduction = "umap", group.by = "Condition", pt.size = 1, cols = alpha(c("gray50","firebrick2"),0.66))
ggsave(filename = "umap_clusters_groups.png", path = "../results/Integrated_Analysis/", width = 13, height = 8,
       units = "in", dpi = 300)

## Combined by cluster
DimPlot(data_combined, reduction = "umap", label = TRUE, pt.size = 1)
ggsave(filename = "umap_clusters.png", path = "../results/Integrated_Analysis/", width = 13, height = 8,
       units = "in", dpi = 300)

p1 <- DimPlot(data_combined, reduction = "tsne", group.by = "Condition", pt.size = 1, cols = alpha(c("gray50","firebrick2"),0.66))
p2 <- DimPlot(data_combined, reduction = "tsne", label = TRUE, pt.size = 1)
p1 + p2
ggsave(filename = "tsne_clusters_groups.png", path = "../results/Integrated_Analysis/", width = 13, height = 8,
       units = "in", dpi = 300)

## To visualize the two conditions side-by-side, we can use the split.by argument to 
## show each condition colored by cluster

## Split by Condition
DimPlot(data_combined, reduction = "umap", split.by = "Condition",pt.size = 1, label = TRUE, repel = T)
ggsave(filename = "umap_split_by_condition.png", path = "../results/Integrated_Analysis/", width = 15, height = 8,
       units = "in", dpi = 300)

###Split by Cluster (Run this when you have classified your clusters; it will allow you the plot sample by sample)
DimPlot(data_combined, reduction = "umap", split.by = "orig.ident", label = FALSE, label.size=2, pt.size = 1.0)
ggsave(filename = "umap_split_by_cluster.png", path = "../results/Integrated_Analysis/", width = 20, height = 8,
       units = "in", dpi = 300)

## Save the combined seurat object

saveRDS(object = data_combined, file = '../results/Integrated_Analysis/seurat_integration_Condition.rds')
data_combined <- readRDS(file = "../results/Integrated_Analysis/seurat_integration_Condition.rds")

#############################################
## 4) Identify conserved cell type markers ##
#############################################

## To identify canonical cell type marker genes that are conserved across conditions, we provide 
## the FindConservedMarkers() function. This function performs differential gene expression testing 
## for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R 
## package. 

Idents(data_combined) <- "seurat_clusters"

for (i in levels(Idents(data_combined))) {
  markers <- FindConservedMarkers(data_combined, ident.1 = i, grouping.var = "Condition", verbose = F)
  markers <- rownames_to_column(markers, var = "gene") %>% as_tibble()
  write_delim(x = markers, file = paste("../results/Integrated_Analysis/markers_cluster_", i, ".tsv", sep = ""), delim = "\t")
}

# For performing differential expression after integration, make sure to switch back to the original
# data

DefaultAssay(data_combined) <- "RNA"

## Hepatic CD45+ cells markers table from Blériot C., et al., 2021 

KCs_markers <- read_excel(path = "../additional_tables/KCs_Cells_Markers_literature.xlsx")

## Violin plots of markers

for (x in c(1:11)) {
  my_features <- KCs_markers$Markers[x] %>% str_split(pattern = ",") %>% unlist()
  file_name <- paste0(KCs_markers$Cell_type[x],"_VlnPlot.png")
  Marker_plot <- VlnPlot(object = data_combined, features = my_features)
  ggsave(filename = file_name, plot = Marker_plot, device = "png", 
         path = "../results/Integrated_Analysis/VlnPlots_Markers/",
         width = 20, height = 10,dpi = 300, units = "in")
}

# Feature plots allows to see the distribution of the expression of each marker in our data

for (x in c(1:10)) {
  my_features <- KCs_markers$Markers[x] %>% str_split(pattern = ",") %>% unlist()
  Marker_plot <- FeaturePlot(object = data_combined,reduction = "umap", features = my_features,label = TRUE, repel = TRUE)
  file_name <- paste0(KCs_markers$Cell_type[x],"_FeaturePlot.png")
  ggsave(filename = file_name, plot = Marker_plot, device = "png", 
         path = "../results/Integrated_Analysis/FeaturePlots_Markers/",
         width = 30, height = 20,dpi = 300)
}


############################
## 5) Assign clusters IDs ##
############################

# Manual annotation of cell types
KCs_markers_clusters <- read_excel(path = "../additional_tables/KCs_Markers_clusters.xlsx")

# Save in a new variable the annotated cell types
cluster_ids <- KCs_markers_clusters$Name

names(cluster_ids) <- levels(data_combined)

data_combined$type <- data_combined@active.ident

# Assign the cell type labels to each cell in the data
data_combined <- RenameIdents(data_combined, cluster_ids)

Idents(data_combined)

## First, let's reorder the clusters so it is easier to interpret in the visualization

clusters_order <- c("B cells","T cells","NK cells", "NKT cells","Monocytes + cDCs", "LAMs","pDCs", "Unknown")

data_combined@active.ident <- factor(x = data_combined@active.ident, levels = clusters_order)
data_combined$type <- Idents(data_combined)

# Visualization
plot1 <- DimPlot(data_combined, reduction = "umap") + NoLegend()
LabelClusters(plot1, id = "ident", color = unique(ggplot_build(plot1)$data[[1]]$colour), size = 9, repel = T,  box.padding = 1)
ggsave(filename = "Integrated_clusters_IDs_umap_TFM.png", path = "../results/Integrated_Analysis/", width = 10,
       height = 9, units = "in", dpi = 300)

plot2 <- DimPlot(data_combined, reduction = "tsne") + NoLegend()
LabelClusters(plot2, id = "ident", color = unique(ggplot_build(plot1)$data[[1]]$colour), size = 10, repel = T,  box.padding = 1)
ggsave(filename = "Integrated_clusters_IDs_tsne_TFM.png", path = "../results/Integrated_Analysis/", width = 10,
       height = 9, units = "in", dpi = 300)

###Split by Cluster (Run this when you have classified your clusters; it will allow you the plot sample by sample)
DimPlot(data_combined, reduction = "umap", split.by = "orig.ident", label = TRUE, label.size=2, pt.size = 1.0)
ggsave(filename = "umap_split_by_cluster_with_IDs.png", path = "../results/Integrated_Analysis/", width = 20, height = 8,
       units = "in", dpi = 300)

DimPlot(data_combined, reduction = "umap", split.by = "Condition",pt.size = 1, label = TRUE, repel = T)
ggsave(filename = "umap_split_by_condition_with_IDs.png", path = "../results/Integrated_Analysis/", width = 20, height = 8,
       units = "in", dpi = 300)
DimPlot(data_combined, reduction = "tsne", split.by = "Condition",pt.size = 1, label = TRUE, repel = T)
ggsave(filename = "tsne_split_by_condition_with_IDs.png", path = "../results/Integrated_Analysis/", width = 20, height = 8,
       units = "in", dpi = 300)

DimPlot(data_combined, reduction = "tsne", pt.size = 1, label = TRUE, repel = T)
ggsave(filename = "tsne_with_IDs.png", path = "../results/Integrated_Analysis/", width = 20, height = 8,
       units = "in", dpi = 300)

VlnPlot(data_combined, features = c("Cx3cr1", "Adgre1", "Ccr2", "Ccl5", "Ctla4", "Cd3g"), ncol = 3, 
        pt.size = 0)
ggsave(filename = "../results/Integrated_Analysis/Violin plots.png",
       width = 10, height = 6, units = "in", dpi = 300)

# Heatmap visualization
marcadores <- FindAllMarkers(data_combined, min.pct = 0.2,logfc.threshold = 1)
## log2order_markers <- marcadores %>% group_by(cluster) %>% slice_max(n=length(marcadores$gene),order_by = avg_log2FC)

top10 <- marcadores %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
DoHeatmap(object = data_combined, features = top10$gene, label = T, size = 3, hjust = T) + NoLegend()
# + FontSize(y.text = 10) +  theme(text = element_text(size = 20)) 
ggsave(filename = "Heatmap_clusters_IDs_top10_sinlabels.png", path = "../results/Integrated_Analysis/", width = 10, height = 15)


## Find strong markers for each cluster

KCs_markers_unique <- KCs_markers_clusters$Marker %>% str_split(pattern = ",") %>% unlist() %>% unique()

VlnPlot(object = data_combined, features = KCs_markers_unique,pt.size = 0)& 
  theme( plot.title = element_text( face = "italic") )
ggsave(filename = "Main_markers_Vlnplot.png",device = "png", 
       path = "../results/Integrated_Analysis/",
       width = 30, height = 30,dpi = 300)

## Identify the most interesting markers to represent each cluster to avoid excesive information

## The DotPlot() function with the split.by parameter can be useful for viewing conserved cell type markers across 
## conditions, showing both the expression level and the percentage of cells in a cluster expressing any given gene. 

markers_to_plot <- c("Cd19","Cd79a","Cd3g","Trac","Ccl5","Ctla4", "Gzmk",
                     "Ccr2","Itgam","Dpp4","Nr1h3","Siglech")

DotPlot(data_combined, features = markers_to_plot,dot.scale = 10, 
        cols = c("white","firebrick2"), assay = "RNA", col.min = 0) +
  RotatedAxis() + theme(axis.text = element_text(size=12), axis.text.x = element_text(face = "italic"))
ggsave(filename = "DotPlot_markers_clusters.png", path = "../results/Integrated_Analysis/",width = 12, height = 6,dpi = 300)

VlnPlot(object = data_combined, features = markers_to_plot, pt.size = 0, ncol = 1,same.y.lims = T)& 
  theme( plot.title = element_text( face = "italic"), axis.title.x = element_blank(), axis.text.x = element_blank())
ggsave(filename = "Selected_markers_Vlnplot.png",device = "png", 
       path = "../results/Integrated_Analysis/",
       width = 15, height = 35,dpi = 300)
ggsave(filename = "Selected_markers_Vlnplot.pdf",device = "pdf", 
       path = "../results/Integrated_Analysis/",
       width = 15, height = 35)

#####################################################
## 6) Calculating cell proportions in each cluster ##
#####################################################

# We see how  many cells are in each cluster
clusters_counts <- table(Idents(data_combined))

# How does cluster membership vary by condition?
conditions_cluster_counts <- table(Idents(data_combined), data_combined$Condition)

conditions_cluster_perc <- prop.table(table(Idents(data_combined), 
                                            data_combined$Condition), margin = 2)*100 ## Multiply by 100 to have percentages

difference_percentage <- conditions_cluster_perc[,2] - conditions_cluster_perc[,1]

cell_percentages_table <- as.data.frame(
  cbind(conditions_cluster_perc,difference_percentage,conditions_cluster_counts,clusters_counts))
colnames(cell_percentages_table) <- c("CLWT(%)","CLCD36(%)","Percentage_difference",
                                      "CLWT","CLCD36","Total_cells")

column_totals <- c(sum(conditions_cluster_perc[,1]),sum(conditions_cluster_perc[,2]),
                   NA,
                   sum(conditions_cluster_counts[,1]),sum(conditions_cluster_counts[,2]),
                   sum(clusters_counts))

cell_percentages_table <- rbind(cell_percentages_table,column_totals)
rownames(cell_percentages_table) <- c(rownames(cell_percentages_table[-8,]),"Total")

write.table(x = cell_percentages_table, file = "../results/Integrated_Analysis/cell_percentages_clusters.txt",
            sep = "\t",row.names = TRUE, col.names = NA)

## Let's create a table with counts and percentages per sample


# We see how does cluster membership vary by sample
conditions_cluster_counts_sample <- table(Idents(data_combined), data_combined$orig.ident)
conditions_cluster_counts_sample <- cbind(conditions_cluster_counts_sample[,3:4],conditions_cluster_counts_sample[,1:2])

conditions_cluster_perc_sample <- prop.table(table(Idents(data_combined), 
                                                   data_combined$orig.ident), margin = 2)*100 ## Multiply by 100 to have percentages
conditions_cluster_perc_sample <- cbind(conditions_cluster_perc_sample[,3:4],conditions_cluster_perc_sample[,1:2])


cell_percentages_table_sample <- as.data.frame(
  cbind(conditions_cluster_perc_sample,conditions_cluster_counts_sample,clusters_counts))

colnames(cell_percentages_table_sample) <- c("LYWT1_(%)","LYWT2_(%)","LYMKK361_(%)","LYMKK362_(%)",
                                             "LYWT1","LYWT2","LYMKK361","LYMKK362","Total_cells")

column_totals_sample <- c(sum(conditions_cluster_perc_sample[,1]),sum(conditions_cluster_perc_sample[,2]),
                          sum(conditions_cluster_perc_sample[,3]),sum(conditions_cluster_perc_sample[,4]),
                          sum(conditions_cluster_counts_sample[,1]),sum(conditions_cluster_counts_sample[,2]),
                          sum(conditions_cluster_counts_sample[,3]),sum(conditions_cluster_counts_sample[,4]),
                          sum(clusters_counts))

cell_percentages_table_sample <- rbind(cell_percentages_table_sample,column_totals_sample)
rownames(cell_percentages_table_sample) <- c(rownames(cell_percentages_table_sample[-8,]),"Total")


write.table(x = cell_percentages_table_sample, file = "../results/Integrated_Analysis/cell_percentages_clusters_sample.txt",
            sep = "\t",row.names = TRUE, col.names = NA)


## Bar Plots of the cell proportions

ggplotColours <- function(n = 6, h = c(0, 360) + 15){ ## Function to extract colors from ggplot (the scheme that DimPlot uses)
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=7) ## Asign as many colors as clusters
color_list_0.4 <- adjustcolor(color_list, alpha.f = 0.4)
color_list_0.8 <- adjustcolor(color_list, alpha.f = 0.8)



bplot <- barplot(conditions_cluster_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Kupffer cells", ylim = c(0,70),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))
# barplot(conditions_cluster_perc,  legend.text = FALSE, col = color_list_0.4,border = color_list, beside = TRUE, 
        # ylab =  quote(bold("Percentage (%)")), 
        # main = "Kupffer cells", ylim = c(0,60), )
text(bplot[,1], (conditions_cluster_perc[,1]+2), paste(round(conditions_cluster_perc[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc[,2]+2), paste(round(conditions_cluster_perc[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
#legend(x="topright", legend = rownames(conditions_cluster_perc), col = color_list_0.4, 
       #inset =c(-0.05,0),cex = 1, bty="n", pch=19, pt.cex = 1.5)
points(x = bplot[,1], y = conditions_cluster_perc_sample[,1], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,1], y = conditions_cluster_perc_sample[,2], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,2], y = conditions_cluster_perc_sample[,3], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,2], y = conditions_cluster_perc_sample[,4], col=color_list_0.8, pch=19, cex=1.1)



png(filename = "../results/Integrated_Analysis/cell_percentages_cluster.png", width = 14, height = 8, units = "in", res = 300)
barplot(conditions_cluster_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
        ylab = quote(bold("Percentage (%)")), 
        main = "Kupffer cells", ylim = c(0,70))
text(bplot[,1], (conditions_cluster_perc[,1]+4), paste(round(conditions_cluster_perc[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc[,2]+4), paste(round(conditions_cluster_perc[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
# legend(x="topright", legend = rownames(conditions_cluster_perc), col = color_list_0.8,cex = 1.5, bty="n", pch=19, pt.cex = 1.5)
points(x = bplot[,1], y = conditions_cluster_perc_sample[,1], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,1], y = conditions_cluster_perc_sample[,2], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,2], y = conditions_cluster_perc_sample[,3], col=color_list_0.8, pch=19, cex=1.1)
points(x = bplot[,2], y = conditions_cluster_perc_sample[,4], col=color_list_0.8, pch=19, cex=1.1)
dev.off()


## Each sample separately

bplot <- barplot(conditions_cluster_perc_sample,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Kupffer cells", ylim = c(0,75),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))
# barplot(conditions_cluster_perc_sample,  legend.text = FALSE, col = color_list_0.4,border = color_list, beside = TRUE, 
       # ylab = quote(bold("Percentage (%)")), 
      #  main = "Kupffer cells", ylim = c(0,55))
# legend(x="topright", legend = rownames(conditions_cluster_perc), col = color_list,cex = 1, bty="n", pch=19, pt.cex = 1.5)

png(filename = "../results/Integrated_Analysis/cell_percentages_cluster_sample.png", width =14, height = 8, units = "in", res = 300)
barplot(conditions_cluster_perc_sample,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
        ylab = quote(bold("Percentage (%)")), 
        main = "Kupffer cells", ylim = c(0,75), args.legend=list(x = "topright", bty = "n", inset=c(-0.03, -0.1)))
#legend(x="topright", legend = rownames(conditions_cluster_perc), col = color_list,cex = 1, bty="n", pch=19, pt.cex = 1.5)
dev.off()

################################################################
## 7) Identify differential expressed genes across conditions ##
################################################################

## Now that we've aligned the KO and control cells, we can start to do comparative 
## analyses and look at the differences induced by HFD

data_combined$celltype_cond <- paste(Idents(data_combined), data_combined$Condition, sep = "_")
data_combined$celltype <- Idents(data_combined)
Idents(data_combined) <- "celltype_cond"

# Select the different cell types we are using in our analysis
clusters_KCs <- unique(cluster_ids)

for (x in 1:length(clusters_KCs)) {
  ## Find DEGs
  print(clusters_KCs[x])
  DEGs <- FindMarkers(object = data_combined, 
                      ident.1 = paste(clusters_KCs[x],"_KO",sep = ""),
                      ident.2 = paste(clusters_KCs[x],"_Control",sep = ""), logfc.threshold = 0)
  ## Save table with DEGs
  write.table(x = DEGs, file = paste("../results/Integrated_Analysis/","DEGS_",clusters_KCs[x],".txt",sep = ""),
              sep = "\t", quote = F)
  
  keyvals <- ifelse(
    DEGs$avg_log2FC < -1 & DEGs$p_val_adj < 0.05, "navy",
    ifelse(DEGs$avg_log2FC  > 1 & DEGs$p_val_adj < 0.05,"#D82632",
           "grey"))
  names(keyvals)[keyvals == "navy"] <- "Downregulated"
  names(keyvals)[keyvals == "#D82632"] <- "Upregulated"
  names(keyvals)[keyvals == "grey"] <- "NS"
  DEGs$gene <- rownames(DEGs)
  # selected_labs <- c(slice_min(.data = DEGs, order_by = p_val_adj, n = 10)$gene)
  volcano_plot <- EnhancedVolcano(toptable = DEGs,lab = rownames(DEGs),
                  x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 1,
                  ylim = , xlim = range(DEGs$avg_log2FC), # selectLab = selected_labs,
                  drawConnectors = TRUE,widthConnectors = 0.5,typeConnectors = ,
                  endsConnectors = 'first',labSize = 4,gridlines.minor = F,gridlines.major = F
                  ,pointSize = 3,colAlpha = 0.5, title = clusters_KCs[x],
                  subtitle = "Adj p-value cutoff (dashed line): p<0.05
  Log2 FC cutoff (dashed line): 1",colCustom = keyvals,
                  labFace = "italic",boxedLabels = TRUE, max.overlaps = Inf)
  ggsave(filename = paste("../results/Integrated_Analysis/","Volcano_DEGs_",clusters_KCs[x],".png",sep = ""),
         plot = volcano_plot,device = png,
         width = 15, height = 10,dpi = 300)
}

## Barplot of DEGS counts in each cluster

DEGs_count <- c()

for (i in 1:length(clusters_KCs)) {
  DEGs <- read.delim(paste("../results/Integrated_Analysis/DEGS_",clusters_KCs[i],".txt",sep = ""),sep = "\t")
  DEGs_count[i*2-1] <- length((DEGs %>% filter(p_val_adj<0.05 & avg_log2FC>0))$p_val)
  DEGs_count[i*2] <- length((DEGs %>% filter(p_val_adj<0.05 & avg_log2FC<0))$p_val)
}

DEGs_df <- data.frame(Cluster=rep(clusters_KCs[1:9],each=2),Expression=rep(c("Up","Down"),9),
                      Count=DEGs_count)

ggplot(data = DEGs_df, aes(x = Cluster, y = Count, fill=Expression)) + 
  geom_bar(stat="identity", position=position_dodge()) + theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "DEGs Count") +
  scale_fill_manual(values = alpha(c("slateblue2","firebrick1"),0.80))
ggsave(filename = "../results/Integrated_Analysis/DEGs_count_barplot.png",width = 7.5, height = 5, units = "in",
       dpi = 300)


###########################
## 8) Automated Anotation##
###########################

library(celldex)
library(SingleR)

data_combined <- readRDS(file = "../results/Integrated_Analysis/seurat_integration_Condition.rds")

# We save in the variable referencia the built-in reference we are going to use for our dataset
referencia <- celldex::MouseRNAseqData()

# Perform automated annotation with SingleR
pbmc.data <- GetAssayData(data_combined, layer = "counts")
pred <- SingleR(test = pbmc.data,
                ref = referencia,
                labels = referencia$label.main)

data_combined$labels <- pred$labels[match(rownames(data_combined@meta.data),rownames(pred))]
esto <- DimPlot(data_combined, reduction = "umap", group.by = "labels")
LabelClusters(esto, id = "labels", color = unique(ggplot_build(esto)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)

# We erase the cell types found as a consequence of contamination that in fact are found in very small quantities
Idents(data_combined) <- data_combined$labels
data_combined <- subset(data_combined, ident = c("B cells", "Granulocytes", "Macrophages", "Monocytes","NK cells", "T cells"))

esto <- DimPlot(data_combined, reduction = "umap", group.by = "labels")
LabelClusters(esto, id = "labels", color = unique(ggplot_build(esto)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)
ggsave("../results/Integrated_Analysis/Unsupervised_annotation.png", width = 14, height = 10, units = "in", dpi = 300)

############################
## 9) Study of cluster 12 ##
############################

Idents(data_combined) <- data_combined$seurat_clusters
cluster12 <- subset(data_combined, ident = "12")
Idents(cluster12) <- cluster12$type
counts <- table(Idents(cluster12))
counts_perc <- prop.table(table(Idents(cluster12)))*100
counts_matrix <- cbind(counts, counts_perc)

## Bar Plots of the cell proportions
FeaturePlot(cluster12, features = c("Cx3cr1", "Adgre1", "Ccr2", "Ccl5", "Ctla4", "Cd3g"), reduction = "umap", ncol = 3)
ggsave(filename = "../results/Integrated_Analysis/Feature plots.png",
       width = 12, height = 8, units = "in", dpi = 300)

# Visualization of cell proportion in cluster 12 at the different conditions
ggplotColours <- function(n = 6, h = c(0, 360) + 15){ ## Function to extract colors from ggplot (the scheme that DimPlot uses)
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=4) ## Asign as many colors as clusters
color_list_0.4 <- adjustcolor(color_list, alpha.f = 0.4)
color_list_0.8 <- adjustcolor(color_list, alpha.f = 0.8)

bplot <- barplot(counts_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Unknown (cluster 12)", ylim = c(0,70),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot, (counts_perc +2), paste(round(counts_perc,digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
points(x = bplot, y = counts_perc, col=color_list_0.8, pch=19, cex=1.1)

png(filename = "../results/Integrated_Analysis/cell_percentages_cluster.png", width = 14, height = 8, units = "in", dpi = 300)
bplot <- barplot(counts_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Unknown (cluster 12)", ylim = c(0,70),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot, (counts_perc +2), paste(round(counts_perc,digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
points(x = bplot, y = counts_perc, col=color_list_0.8, pch=19, cex=1.1)

ggsave("../results/Integrated_Analysis/cell_percentages_cluster.png", width = 14, height = 8, units = "in", dpi = 300)

