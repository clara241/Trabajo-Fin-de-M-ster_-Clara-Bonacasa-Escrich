###########################
## scRNA-seq INTEGRATION ##
###########################

# 0) Load all the required Libraries
library(Seurat)
library(readr)
library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)
library(Nebulosa)
library(SingleR)
library(Biobase)
library(plotly)
library(patchwork) 

###################################################################
## 1) Read a collection of data sets for single cell integration ##
###################################################################

## 1.1 Read filtered seurat objects of each sample for integration 

H1 <- read_rds("../results/Healthy_1/Healthy_1_data_combined.rds")
H4 <- read_rds("../results/Healthy_4/Healthy_4_data_combined.rds")
C1 <- read_rds("../results/Cirrhotic_1/Cirrhotic_1_data_combined.rds")
C4 <- read_rds("../results/Cirrhotic_4/Cirrhotic_4_data_combined.rds")

## 1.2 Add extra column to the Seurat objects with the condition
H1 <- AddMetaData(object = H1, metadata = "Healthy", col.name = "Condition")
H4 <- AddMetaData(object = H4, metadata = "Healthy", col.name = "Condition")
C1 <- AddMetaData(object = C1, metadata = "Cirrhosis", col.name = "Condition")
C4 <- AddMetaData(object = C4, metadata = "Cirrhosis", col.name = "Condition")

# Merge seurat objects adding to each of them an identification/sample name
data_combined <- merge(x = H1, y = c(H4,C1,C4),
                       add.cell.ids = c("Healthy_1","Healthy_4","Cirrhosis_1", "Cirrhosis_4"))

## We erase these variables as they are not going to be used in the following analysis
rm(H1,H4,C1,C4)

## 1.3. Unintegrated analysis

data_combined[["RNA"]] <- split(data_combined[["RNA"]], f = data_combined$Condition)

# We run standard analysis workflow
data_combined <- NormalizeData(data_combined)
data_combined <- FindVariableFeatures(data_combined)
data_combined <- ScaleData(data_combined)
data_combined <- RunPCA(data_combined)

## 1.4. Integration
## The aim of the integration is to link the datasets from the two different conditions
## so the cell types or subpopulations will cluster together. In order to achieve that,
## seurat v5 integration gives an unified dimensional reduction that captures the 
## common sources of variation across several layers. This enables cells in similar
## biological states to form distinct clusters
data_combined <- IntegrateLayers(object = data_combined, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = FALSE)

## we re-join layer after integration
data_combined[["RNA"]] <- JoinLayers(data_combined[["RNA"]])

## We are going to determine the 'dimensionality' of the dataset. Here, we use an
## heuristic method that generates a ranking with the top principal components, from PCA,
## ordered according to the percentage of variance explained by each one.
ElbowPlot(data_combined, ndims = 50, reduction = "pca")

data_combined <- FindNeighbors(data_combined, reduction = "integrated.cca", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.5)

data_combined <- RunUMAP(data_combined, dims = 1:30, reduction = "integrated.cca")
data_combined <- RunTSNE(data_combined, dims = 1:30, reduction = "integrated.cca")

## Visualization
DimPlot(data_combined, reduction = "umap", label = T)
ggsave("../results/Integrated_Analysis/umap_seurat_clusters.png", height = 10, width = 12, units = "in", dpi = 300)

DimPlot(data_combined, reduction = "tsne", label = T)
ggsave("../results/Integrated_Analysis/tsne_seurat_clusters.png", height = 10, width = 12, units = "in", dpi = 300)

DimPlot(data_combined, reduction = "umap", group.by = c("orig.ident", "Condition"))
ggsave("../results/Integrated_Analysis/umap_condition_origident.png", height = 10, width = 12, units = "in", dpi = 300)

## Save the combined seurat object
saveRDS(object = data_combined, file = '../results/Integrated_Analysis/seurat_integration_Condition.rds')
data_combined <- readRDS(file = "../results/Integrated_Analysis/seurat_integration_Condition.rds")

#################################
## 2. Identify cluster markers ##
#################################
## We read an excel file containing cell type and markers specific to each cell lineage
## based on the paper where these data sets were retrieve from https://doi.org/10.1038/s41598-021-98806-y
liver_markers <- read_excel(path = "../additional_tables/Liver_cells_markers_literature.xlsx")

## Violin plots of markers

for (x in c(1:14)) {
  my_features <- liver_markers$Markers[x] %>% str_split(pattern = ",") %>% unlist()
  file_name <- paste0(liver_markers$Cell_type[x],"_VlnPlot.png")
  Marker_plot <- VlnPlot(object = data_combined, features = my_features)
  ggsave(filename = file_name, plot = Marker_plot, device = "png", 
         path = "../results/Integrated_Analysis/VlnPlots_Markers/",
         width = 20, height = 10,dpi = 300, units = "in")
}

############################
## 3. Assign clusters IDs ##
############################

liver_marker_cells <- read_excel(path = "../additional_tables/Liver_cluster_markers.xlsx")
cluster_ids <- liver_marker_cells$Name

names(cluster_ids) <- levels(data_combined)

data_combined$type <- data_combined@active.ident

data_combined <- RenameIdents(data_combined, cluster_ids)

Idents(data_combined)

## First, let's reorder the clusters so it is easier to interpret in the visualization

clusters_order <- c("HSCs","LSECs 1","LSECs 2", "Cholangiocytes","T cells","NK cells", "B cells", 
                    "Plasma cells", "Dendritic cells","MDMs", "Kupffer cells", "Cycling cells")

data_combined@active.ident <- factor(x = data_combined@active.ident, levels = clusters_order)
data_combined$cell_type <- Idents(data_combined)

### Visualization
plot1 <- DimPlot(data_combined, reduction = "umap") + NoLegend()
LabelClusters(plot1, id = "ident", color = unique(ggplot_build(plot1)$data[[1]]$colour), size = 9, repel = T,  box.padding = 2.5)
ggsave(filename = "Integrated_clusters_IDs_umap.png", path = "../results/Integrated_Analysis/", width =11,
       height = 9, units = "in", dpi = 300)

plot2 <- DimPlot(data_combined, reduction = "tsne") + NoLegend()
LabelClusters(plot2, id = "ident", color = unique(ggplot_build(plot1)$data[[1]]$colour), size = 10, repel = T,  box.padding = 1)
ggsave(filename = "Integrated_clusters_IDs_tsne.png", path = "../results/Integrated_Analysis/", width = 10,
       height = 9, units = "in", dpi = 300)

plot3 <- DimPlot(data_combined, reduction = "umap", split.by = "Condition") + NoLegend()
LabelClusters(plot3, id = "ident", color = rep(unique(ggplot_build(plot3)$data[[1]]$colour),2), size = 6, repel = T,  box.padding = 3)
ggsave(filename = "Integrated_clusters_splitby_condition_umap.png", path = "../results/Integrated_Analysis/", width = 14,
       height = 9, units = "in", dpi = 300)

DimPlot(data_combined, reduction = "umap", group.by = "Condition")
ggsave(filename = "Integrated_clusters_condition.png", path = "../results/Integrated_Analysis/", width = 10,
       height = 9, units = "in", dpi = 300)


DotPlot(data_combined, features = c("COL1A1","DCN","SPARCL1","CLEC14A","STAB1",
                                    "FCN2","SOX9","EPCAM","CD3D","TRAC","NKG7","GNLY","MS4A1",
                                    "CD79A","IGLC2","IGHG1","IRF8","LILRA4","MNDA","CD68",
                                    "MARCO","TIMD4","UBE2C","TOP2A"), dot.scale = 10,
        cols = c("white","firebrick2"), assay = "RNA", col.min = 0) +
        RotatedAxis() + theme(axis.text = element_text(size=12), axis.text.x = element_text(face = "italic"))
ggsave("../results/Integrated_Analysis/DotPlot_gene_expression.png",width = 12,height = 10, units = "in", dpi = 600)
#####################################################
## 6) Calculating cell proportions in each cluster ##
#####################################################

# How many cells are in each cluster
clusters_counts <- table(Idents(data_combined))

# How does cluster membership vary by condition?
conditions_cluster_counts <- table(Idents(data_combined), data_combined$Condition)

conditions_cluster_counts <- cbind(conditions_cluster_counts[,2], conditions_cluster_counts[,1])

conditions_cluster_perc <- prop.table(table(Idents(data_combined), 
                                            data_combined$Condition), margin = 2)*100 ## Multiply by 100 to have percentages

conditions_cluster_perc <- cbind(conditions_cluster_perc[,2],conditions_cluster_perc[,1])

difference_percentage <- conditions_cluster_perc[,2] - conditions_cluster_perc[,1]

cell_percentages_table <- as.data.frame(
  cbind(conditions_cluster_perc,difference_percentage,conditions_cluster_counts,clusters_counts))
colnames(cell_percentages_table) <- c("Healthy(%)","Cirrhosis(%)","Percentage_difference",
                                      "Healthy","Cirrhosis","Total_cells")

column_totals <- c(sum(conditions_cluster_perc[,1]),sum(conditions_cluster_perc[,2]),
                   NA,
                   sum(conditions_cluster_counts[,1]),sum(conditions_cluster_counts[,2]),
                   sum(clusters_counts))

cell_percentages_table <- rbind(cell_percentages_table,column_totals)
rownames(cell_percentages_table) <- c(rownames(cell_percentages_table[-13,]),"Total")

write.table(x = cell_percentages_table, file = "../results/Integrated_Analysis/cell_percentages_clusters.txt",
            sep = "\t",row.names = TRUE, col.names = NA)

## Let's create a table with counts and percentages per sample

# How does cluster membership vary by sample?
conditions_cluster_counts_sample <- table(Idents(data_combined), data_combined$orig.ident)
conditions_cluster_counts_sample <- cbind(conditions_cluster_counts_sample[,3:4],conditions_cluster_counts_sample[,1:2])

conditions_cluster_perc_sample <- prop.table(table(Idents(data_combined), 
                                                   data_combined$orig.ident), margin = 2)*100 ## Multiply by 100 to have percentages
conditions_cluster_perc_sample <- cbind(conditions_cluster_perc_sample[,3:4],conditions_cluster_perc_sample[,1:2])


cell_percentages_table_sample <- as.data.frame(
  cbind(conditions_cluster_perc_sample,conditions_cluster_counts_sample,clusters_counts))

colnames(cell_percentages_table_sample) <- c("Healthy_1(%)","Healthy_4(%)","cirrhotic_1(%)","cirrhotic_4(%)",
                                             "Healthy_1","Healthy_4","cirrhotic_1","cirrhotic_4",
                                             "Total_cells")

column_totals_sample <- c(sum(conditions_cluster_perc_sample[,1]),sum(conditions_cluster_perc_sample[,2]),
                          sum(conditions_cluster_perc_sample[,3]),sum(conditions_cluster_perc_sample[,4]),
                          sum(conditions_cluster_counts_sample[,1]),sum(conditions_cluster_counts_sample[,2]),
                          sum(conditions_cluster_counts_sample[,3]),sum(conditions_cluster_counts_sample[,4]),
                          sum(clusters_counts))

cell_percentages_table_sample <- rbind(cell_percentages_table_sample,column_totals_sample)
rownames(cell_percentages_table_sample) <- c(rownames(cell_percentages_table_sample[-13,]),"Total")


write.table(x = cell_percentages_table_sample, file = "../results/Integrated_Analysis/cell_percentages_clusters_sample.txt",
            sep = "\t",row.names = TRUE, col.names = NA)


## Bar Plots of the cell proportions

ggplotColours <- function(n = 6, h = c(0, 360) + 15){ ## Function to extract colors from ggplot (the scheme that DimPlot uses)
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=12) ## Asign as many colors as clusters
color_list_0.4 <- adjustcolor(color_list, alpha.f = 0.4)
color_list_0.8 <- adjustcolor(color_list, alpha.f = 0.8)

colnames(conditions_cluster_perc) <- c("Healthy", "Cirrhosis")
bplot <- barplot(conditions_cluster_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Liver cells", ylim = c(0,45),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot[,1], (conditions_cluster_perc[,1]+2), paste(round(conditions_cluster_perc[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc[,2]+2), paste(round(conditions_cluster_perc[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
#points(x = bplot[,1], y = conditions_cluster_perc_sample[,1], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,1], y = conditions_cluster_perc_sample[,2], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,3], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,4], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,5], col=color_list_0.8, pch=19, cex=1.1)

png(filename = "../results/Integrated_Analysis/cell_percentages_cluster.png", width = 17, height = 11, units = "in", res = 600)
barplot(conditions_cluster_perc,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
        ylab = quote(bold("Percentage (%)")), 
        main = "Liver cells", ylim = c(0,45))
text(bplot[,1], (conditions_cluster_perc[,1]+2), paste(round(conditions_cluster_perc[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc[,2]+2), paste(round(conditions_cluster_perc[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)

#points(x = bplot[,1], y = conditions_cluster_perc_sample[,1], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,1], y = conditions_cluster_perc_sample[,2], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,3], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,4], col=color_list_0.8, pch=19, cex=1.1)
#points(x = bplot[,2], y = conditions_cluster_perc_sample[,5], col=color_list_0.8, pch=19, cex=1.1)
dev.off()

## Plotting cell proportion of Kupffer cell across conditions
conditions_cluster_perc <- conditions_cluster_perc[11,]

bplot <- barplot(conditions_cluster_perc,  legend.text = F, col = "#FFB5C5",border = "#FF69B4", beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Kupffer cells", ylim = c(0,10),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot, (conditions_cluster_perc), paste(round(conditions_cluster_perc,digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)

png(filename = "../results/Integrated_Analysis/Kupffer_cell_percentages.png", width = 6, height = 7, units = "in", res = 300)

bplot <- barplot(conditions_cluster_perc,  legend.text = F, col = "#FFB5C5",border = "#FF69B4", beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Kupffer cells", ylim = c(0,10),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot, (conditions_cluster_perc), paste(round(conditions_cluster_perc,digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
dev.off()


## Plotting cell proportions of each sample
bplot <- barplot(conditions_cluster_perc_sample,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Liver cells", ylim = c(0,65),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot[,1], (conditions_cluster_perc_sample[,1]+2), paste(round(conditions_cluster_perc_sample[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc_sample[,2]+2), paste(round(conditions_cluster_perc_sample[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,3], (conditions_cluster_perc_sample[,3]+2), paste(round(conditions_cluster_perc_sample[,3],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,4], (conditions_cluster_perc_sample[,4]+2), paste(round(conditions_cluster_perc_sample[,4],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)

png(filename = "../results/Integrated_Analysis/cell_percentages_cluster_samples.png", width = 18, height = 13, units = "in", res = 300)
bplot <- barplot(conditions_cluster_perc_sample,  legend.text = TRUE, col = color_list_0.4,border = color_list, beside = TRUE, 
                 ylab = quote(bold("Percentage (%)")), 
                 main = "Liver cells", ylim = c(0,65),
                 args.legend = list(x="topright", inset =c(-0.05,0),cex = 0.8, bty="n"))

text(bplot[,1], (conditions_cluster_perc_sample[,1]+2), paste(round(conditions_cluster_perc_sample[,1],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,2], (conditions_cluster_perc_sample[,2]+2), paste(round(conditions_cluster_perc_sample[,2],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,3], (conditions_cluster_perc_sample[,3]+2), paste(round(conditions_cluster_perc_sample[,3],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
text(bplot[,4], (conditions_cluster_perc_sample[,4]+2), paste(round(conditions_cluster_perc_sample[,4],digits = 2)," %",sep = ""),
     cex = 1, srt=33, adj = 0.5)
dev.off()

################################################################
## 7) Identify differential expressed genes across conditions ##
################################################################

## Now that we've aligned the stimulated and control cells, we can start to do comparative 
## analyses and look at the differences induced by stimulation.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(EnhancedVolcano)



data_combined$celltype_cond <- paste(Idents(data_combined), data_combined$Condition, sep = "_")
data_combined$celltype <- Idents(data_combined)
Idents(data_combined) <- "celltype_cond"


clusters_liver <- unique(cluster_ids)

# I am using CLWT as the control
for (x in 1:length(clusters_liver)) {
  ## Find DEGs
  print(clusters_liver[x])
  DEGs <- FindMarkers(object = data_combined, 
                      ident.1 = paste(clusters_liver[x],"_Cirrhosis",sep = ""),
                      ident.2 = paste(clusters_liver[x],"_Healthy",sep = ""), logfc.threshold = 0.1)
  ## Save table with DEGs
  write.table(x = DEGs, file = paste("../results/Integrated_Analysis/","DEGS_",clusters_liver[x],".txt",sep = ""),
              sep = "\t", quote = F)
  
  keyvals <- ifelse(
    DEGs$avg_log2FC < -1 & DEGs$p_val_adj < 0.05, "navy",
    ifelse(DEGs$avg_log2FC  > 1 & DEGs$p_val_adj < 0.05,"#D82632",
           "grey"))
  names(keyvals)[keyvals == "navy"] <- "Downregulated"
  names(keyvals)[keyvals == "#D82632"] <- "Upregulated"
  names(keyvals)[keyvals == "grey"] <- "NS"
  DEGs$gene <- rownames(DEGs)
  
  selected_labs <- DEGs %>% filter(p_val_adj < 0.001 & abs(avg_log2FC) > 2) 
  selected_labs <- c(slice_min(.data = selected_labs, order_by = p_val_adj, n = 30)$gene)
  volcano_plot <- EnhancedVolcano(toptable = DEGs,lab = rownames(DEGs),
                                  x = "avg_log2FC",y = "p_val_adj",pCutoff = 0.05,FCcutoff = 1,
                                  ylim = , xlim = range(DEGs$avg_log2FC), selectLab = selected_labs,
                                  drawConnectors = TRUE,widthConnectors = 0.5,typeConnectors = ,
                                  endsConnectors = 'first',labSize = 4,gridlines.minor = F,gridlines.major = F
                                  ,pointSize = 3,colAlpha = 0.5, title = clusters_liver[x],
                                  subtitle = "Adj p-value cutoff (dashed line): p<0.05
  Log2 FC cutoff (dashed line): 1",colCustom = keyvals,
                                  labFace = "italic",boxedLabels = TRUE, max.overlaps = Inf)
  ggsave(filename = paste("../results/Integrated_Analysis/Volcano_plots/","Volcano_DEGs_",clusters_liver[x],".png",sep = ""),
         plot = volcano_plot,device = png,
         width = 15, height = 10,dpi = 300)
}

#############################
## 8) Automated annotation ##
#############################

library(celldex)
library(SingleR)
referencia <- celldex::HumanPrimaryCellAtlasData()

pbmc.data <- GetAssayData(data_combined, layer = "counts")
pred <- SingleR(test = pbmc.data,
                ref = referencia,
                labels = referencia$label.main)

data_combined$labels <- pred$labels[match(rownames(data_combined@meta.data),rownames(pred))]
esto <- DimPlot(data_combined, reduction = "umap", group.by = "labels")
LabelClusters(esto, id = "labels", color = unique(ggplot_build(esto)$data[[1]]$colour), size = 6, repel = T,  box.padding = 1)
 
###########################
## 9) IPA representation ##
###########################
can_path <- read_excel("../results/Integrated_Analysis_CD45/IPA/Canonical_Kupffer.xlsx")
colnames(can_path)[1] <- "Pathway"

can_path[1:8,] %>% arrange(desc(`-log(p-value)`)) %>% 
  mutate(Pathway = factor(x = Pathway, levels = unique(Pathway))) %>% # Reorder factors (genes) by AvA_AAVgd
  ggplot(aes(x = `-log(p-value)`, y = Pathway, size = Ratio, fill = `z-score`)) + # Create plot, choose fill = z-score to color the points later
  geom_point(alpha = 0.8, shape = 21, color = "gray50") + # select alpha, shape and color of the points
  #geom_text(aes(x = AvA_AAVgd, y = Gene, label = Phosphosite, size =1)) + # to add the phosphosites
  scale_fill_gradient2(low = "navy", mid = "white", high = "orange", midpoint = 0,  name = "Activation\nz-score") + # color the points according to p-value
  scale_size(range = c(1, 10), name="Ratio") + # size the points according to number of proteins
  scale_y_discrete(limits = rev) + # Reorder labels in axis
  theme_cowplot() + geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  xlim(1, 7) +
  ylab("Canonical Pathways") + labs(title = "Lipid Metabolism Pathways")
ggsave(filename = "../results/Integrated_Analysis_CD45/IPA/results/Kupffer_canonical.png", width = 8, height = 7, dpi = 600, units = "in")  


