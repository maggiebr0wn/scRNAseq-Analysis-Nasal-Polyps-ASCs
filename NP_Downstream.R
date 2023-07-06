#!/usr/bin/env Rscript

library(pheatmap)
library(readxl)
library(Seurat)

############################################################# 
#               9-15-2020                                   # 
# 1) Isotype labels:                                        #
#      i.) VDJ labels if available                          #
#      ii.) If no VDJ available, LogNorm predicted labels   #
# 2) Redo clustering analysis:                              #
#     i.) Unintegrated                                      #
#     ii.) Integrated                                       # 
#     iii.) DEG testing                                     #
############################################################# 

best_ascs <- readRDS("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")
# 18953 features across 6449 samples within 1 assay 

############################
# Redo Clustering Analysis #
############################


### UNINTEGRATED ###

ascs_unint <- best_ascs
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Unintegrated")

### Recluster
ascs_unint <- NormalizeData(ascs_unint, normalization.method = "LogNormalize", scale.factor = 10000)
ascs_unint <- FindVariableFeatures(ascs_unint, selection.method = "vst", nfeatures = 2750)
ascs_unint <- ScaleData(ascs_unint)
ascs_unint <- RunPCA(ascs_unint, npcs = 30)

# PCA Plots:
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$orig.ident)
pca_orig <- DimPlot(ascs_unint, reduction = "pca")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
pca_iso <- DimPlot(ascs_unint, reduction = "pca")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Label_Source)
pca_lab <- DimPlot(ascs_unint, reduction = "pca")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$mtGroup)
pca_mt <- DimPlot(ascs_unint, reduction = "pca")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$GNGroup)
pca_gn <- DimPlot(ascs_unint, reduction = "pca")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$IgGroup)
pca_ig <- DimPlot(ascs_unint, reduction = "pca")

pdf(file = "PCA_plots.pdf")
pca_orig
pca_iso
pca_lab
pca_mt
pca_gn
pca_ig
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(ascs_unint) # Choose 2:13 PCs
dev.off()

# Cluster
ascs_unint <- FindNeighbors(ascs_unint, dims = 2:13)
ascs_unint <- FindClusters(ascs_unint, resolution = 0.5)
ascs_unint <- RunUMAP(ascs_unint, dims = 2:13)

# UMAP Plots
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$seurat_clusters)
umap_clusts <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$orig.ident)
umap_orig <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
umap_iso <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Label_Source)
umap_lab <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$mtGroup)
umap_mt <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$GNGroup)
umap_gn <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$IgGroup)
umap_ig <- DimPlot(ascs_unint, reduction = "umap", label = TRUE)

pdf(file = "UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
umap_lab
umap_mt
umap_gn
umap_ig
dev.off()

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$seurat_clusters)
umap_split <- DimPlot(ascs_unint, reduction = "umap", split.by = "orig.ident")
umap_isosplit <- DimPlot(ascs_unint, reduction = "umap", split.by = "Isotype")

ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
umap_split_Iso <- DimPlot(ascs_unint, reduction = "umap", split.by = "orig.ident")
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$orig.ident)
umap_isosplit_Iso <- DimPlot(ascs_unint, reduction = "umap", split.by = "Isotype")

pdf(file = "UMAPs_split.pdf", width = 12, height = 3)
umap_split
umap_isosplit
umap_split_Iso
umap_isosplit_Iso
dev.off()

# Plot Features
markerGenes <- c(
  "MS4A1", # non-ASC B Cells
  "JCHAIN", # ASCs
  "SDC1", # CD138 - for ASC subpops ID
  "CD19", # for ASC subpops ID
  "CD38", # for ASC subpops ID
  "IRF4", "PRDM1", "BCL6", "IL10",
  "XBP1", "ZBTB20", "PAX5", "CD40", 
  "IL4", "MITF", "MTA3", "BACH2",
  "CR2", "FCER2", "IL6R", "CD27", "CD24",
  "TBX21", "FAS", "CD86", "ITGAX"
)

HomingMarkers <- c(
  "CXCR3", "CXCR4", "CXCR5", "CXCR6", # Homing Markers
  "ITGA1", "ITGA4", "ITGAL", "ITGB1", "ITGB2", "ITGB7", # Homing Markers
  "CCR7", "CCR10", # Homing Markers
  "SELL", # SELL = CD62L, GPR28 = CCR9
  "S1PR1", "GPR183" # GPR183 = EBI2
)


ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
pdf(file = "ASCs_subpop_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs_unint, reduction = "umap", features = markerGenes)
DotPlot(ascs_unint, features = markerGenes, assay = "RNA") + RotatedAxis()
FeaturePlot(ascs_unint, reduction = "umap", features = HomingMarkers)
DotPlot(ascs_unint, features = HomingMarkers, assay = "RNA") + RotatedAxis()
dev.off()


isotype_genes <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4",
                   "IGHA1", "IGHA2",
                   "IGHD",
                   "IGHM",
                   "IGHE")

pdf(file = "Isotype_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs_unint, features = isotype_genes)
DotPlot(ascs_unint, features = isotype_genes, assay = "RNA") + RotatedAxis()
dev.off()

pdf(file = "VlnPlots.pdf", width = 10, height = 10)
VlnPlot(ascs_unint, features = isotype_genes, sort = TRUE, pt.size = 0.5)
VlnPlot(ascs_unint, features = markerGenes, sort = TRUE, pt.size = 0.5)
VlnPlot(ascs_unint, features = HomingMarkers, sort = TRUE, pt.size = 0.5)
dev.off()

# Heatmap 
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
pdf(file = "heatmap_isotype.pdf", width = 15, height = 5)
DoHeatmap(ascs_unint, draw.lines = FALSE, features = isotype_genes, angle = 90, size = 2)
dev.off()

saveRDS(ascs_unint, "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Unintegrated.rds")

### INTEGRATED ###

ascs_int <- best_ascs
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Integrated")

integrate.list <- SplitObject(ascs_int, split.by = "orig.ident")
for (i in 1:length(integrate.list)) {
  integrate.list[[i]] <- NormalizeData(integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(integrate.list[[i]], selection.method = "vst", 
                                              nfeatures = 2750, verbose = FALSE)
}
reference.list <- integrate.list[c("cls063_ascs", "cls063_pop5", "cls071_ascs", "cls071_pop5", "cls072_ascs", "cls072_pop5")]
ascs.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
ascs_int <- IntegrateData(anchorset = ascs.anchors, dims = 1:30)
DefaultAssay(ascs_int) <- "integrated"

### Recluster
ascs_int <- ScaleData(ascs_int)
ascs_int <- RunPCA(ascs_int, npcs = 30)

# PCA Plots:
ascs_int <- SetIdent(ascs_int, value = ascs_int$orig.ident)
pca_orig <- DimPlot(ascs_int, reduction = "pca")

ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
pca_iso <- DimPlot(ascs_int, reduction = "pca")

ascs_int <- SetIdent(ascs_int, value = ascs_int$Label_Source)
pca_lab <- DimPlot(ascs_int, reduction = "pca")

ascs_int <- SetIdent(ascs_int, value = ascs_int$mtGroup)
pca_mt <- DimPlot(ascs_int, reduction = "pca")

ascs_int <- SetIdent(ascs_int, value = ascs_int$GNGroup)
pca_gn <- DimPlot(ascs_int, reduction = "pca")

ascs_int <- SetIdent(ascs_int, value = ascs_int$IgGroup)
pca_ig <- DimPlot(ascs_int, reduction = "pca")

pdf(file = "PCA_plots.pdf")
pca_orig
pca_iso
pca_lab
pca_mt
pca_gn
pca_ig
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(ascs_int) # Choose 2:12 PCs
dev.off()

# Cluster
ascs_int <- FindNeighbors(ascs_int, dims = 2:12)
ascs_int <- FindClusters(ascs_int, resolution = 0.5)
ascs_int <- RunUMAP(ascs_int, dims = 2:12)

# UMAP Plots
ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)
umap_clusts <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$orig.ident)
umap_orig <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
umap_iso <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$Label_Source)
umap_lab <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$mtGroup)
umap_mt <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$GNGroup)
umap_gn <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

ascs_int <- SetIdent(ascs_int, value = ascs_int$IgGroup)
umap_ig <- DimPlot(ascs_int, reduction = "umap", label = TRUE)

pdf(file = "UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
umap_lab
umap_mt
umap_gn
umap_ig
dev.off()

ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)
umap_split <- DimPlot(ascs_int, reduction = "umap", split.by = "orig.ident")
umap_isosplit <- DimPlot(ascs_int, reduction = "umap", split.by = "Isotype")

ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
umap_split_Iso <- DimPlot(ascs_int, reduction = "umap", split.by = "orig.ident")
ascs_int <- SetIdent(ascs_int, value = ascs_int$orig.ident)
umap_isosplit_Iso <- DimPlot(ascs_int, reduction = "umap", split.by = "Isotype")

pdf(file = "UMAPs_split.pdf", width = 12, height = 3)
umap_split
umap_isosplit
umap_split_Iso
umap_isosplit_Iso
dev.off()

# Plot Features
ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
pdf(file = "ASCs_subpop_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs_int, reduction = "umap", features = markerGenes)
DotPlot(ascs_int, features = markerGenes, assay = "RNA") + RotatedAxis()
FeaturePlot(ascs_int, reduction = "umap", features = HomingMarkers)
DotPlot(ascs_int, features = HomingMarkers, assay = "RNA") + RotatedAxis()
dev.off()

pdf(file = "Isotype_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs_int, features = isotype_genes)
DotPlot(ascs_int, features = isotype_genes, assay = "RNA") + RotatedAxis()
dev.off()

pdf(file = "VlnPlots.pdf", width = 10, height = 10)
VlnPlot(ascs_int, assay = "RNA", features = isotype_genes, sort = TRUE, pt.size = 0.5)
VlnPlot(ascs_int, assay = "RNA", features = markerGenes, sort = TRUE, pt.size = 0.5)
VlnPlot(ascs_int, assay = "RNA", features = HomingMarkers, sort = TRUE, pt.size = 0.5)
dev.off()

# Heatmap 
ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
pdf(file = "heatmap_isotype.pdf", width = 15, height = 5)
DoHeatmap(ascs_int, draw.lines = FALSE, features = isotype_genes, angle = 90, size = 2)
dev.off()

saveRDS(ascs_int, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

#########################
#      DEG TESTING:     #
# Use Unintegrated data #
#########################

### DEGs Testing: Each Isotype vs All; CD19+ and Pop5 ###
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/DEGs/CD19_and_Pop5")

## DEGs Testing: Each Isotype vs All; CD19+ and Pop5 Cells ##
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
for (i in levels(ascs_unint)){
  print(i)
  if(i != "IgD"){ # Only 2 IgD cells
    markers <- FindMarkers(ascs_unint, ident.1 = i, ident.2 = NULL, assay = "RNA")
    filenme <- paste(i, "_vs_ALL_DEGs.csv", sep = "")
    write.csv(markers, file = filenme)}
}

## DEGs Testing: Each Isotype vs All; CD19+ and Pop5 Cells ##
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$Isotype)
for (i in levels(ascs_unint)){
  print(i)
  if(i != "IgD"){ # Only 2 IgD cells
    markers <- FindMarkers(ascs_unint, ident.1 = i, ident.2 = NULL, assay = "RNA")
    filenme <- paste(i, "_vs_ALL_DEGs.csv", sep = "")
    write.csv(markers, file = filenme)}
}

## DEGs Testing: Pairwise Isotypes; CD19+ and Pop5 Cells ##
for (i in levels(ascs_unint)){
  if(i != "IgD"){
    for (j in levels(ascs_unint)){
      if (i != j){
        if(j != "IgD"){
          markers <- FindMarkers(ascs_unint, ident.1 = i, ident.2 = j, assay = "RNA")
          filenme <- paste(i, "_vs_", j, "_DEGs.csv", sep = "")
          write.csv(markers, file = filenme) }}}}
}

## Dotplots for IgE vs ALL DEGs; CD19+ and Pop5 Cells##
ige_vs_all_1to100 <- read_excel("IgE_vs_ALL_DEGs_Significant_Only.xls", sheet = "1-100")
ige_vs_all_101to201 <- read_excel("IgE_vs_ALL_DEGs_Significant_Only.xls", sheet = "101-201")

levels(ascs_cd19_only) <- c("IgM", "IgG4","IgG3",  "IgG2", "IgG1", "IgE", "IgD", "IgA2", "IgA1")
pdf(file = "Plots/IgE_vs_All_CD19_Only.pdf", width = 20, height = 5)
DotPlot(ascs_cd19_only, features = ige_vs_all_1to80$gene, assay = "RNA") + RotatedAxis()
DotPlot(ascs_cd19_only, features = ige_vs_all_81to161$gene, assay = "RNA") + RotatedAxis()
dev.off()


### DEGs Testing: Each Isotype vs All; CD19+ Only ###
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/DEGs/CD19_Only")
ascs_unint <- SetIdent(ascs_unint, value = ascs_unint$orig.ident)
levels(ascs_final)

cd19 <- c("cls063_ascs", "cls071_ascs", "cls072_ascs")
ascs_cd19_only <- subset(ascs_unint, idents = cd19)

## DEGs Testing: Each Isotype vs All; CD19+ Only Cells ##
ascs_cd19_only <- SetIdent(ascs_cd19_only, value = ascs_cd19_only$Isotype)
for (i in levels(ascs_cd19_only)){
  print(i)
  if(i != "IgD"){ # Only 2 IgD cells
    markers <- FindMarkers(ascs_cd19_only, ident.1 = i, ident.2 = NULL, assay = "RNA")
    filenme <- paste(i, "_vs_ALL_DEGs.csv", sep = "")
    write.csv(markers, file = filenme)}
}

## DEGs Testing: Each Isotype vs All; CD19+ Only Cells ##
ascs_cd19_only <- SetIdent(ascs_cd19_only, value = ascs_cd19_only$Isotype)
for (i in levels(ascs_cd19_only)){
  print(i)
  if(i != "IgD"){ # Only 2 IgD cells
    markers <- FindMarkers(ascs_cd19_only, ident.1 = i, ident.2 = NULL, assay = "RNA")
    filenme <- paste(i, "_vs_ALL_DEGs.csv", sep = "")
    write.csv(markers, file = filenme)}
}

## DEGs Testing: Pairwise Isotypes; CD19+ Only Cells ##
for (i in levels(ascs_cd19_only)){
  if(i != "IgD"){
    for (j in levels(ascs_cd19_only)){
      if (i != j){
        if(j != "IgD"){
          markers <- FindMarkers(ascs_cd19_only, ident.1 = i, ident.2 = j, assay = "RNA")
          filenme <- paste(i, "_vs_", j, "_DEGs.csv", sep = "")
          write.csv(markers, file = filenme) }}}}
}

## Dotplots for IgE vs ALL DEGs; CD19+ Only ##
ige_vs_all_up <- read_excel("IgE_vs_ALL_DEGs_Significant_Only.xlsx", sheet = "Upregulated")
ige_vs_all_down <- read_excel("IgE_vs_ALL_DEGs_Significant_Only.xlsx", sheet = "Downregulated")

pdf(file = "Plots/IgE_vs_All_CD19_Only.pdf", width = 25, height = 5)
DotPlot(ascs_cd19_only, features = ige_vs_all_up$gene, assay = "RNA") + RotatedAxis()
DotPlot(ascs_cd19_only, features = ige_vs_all_down$gene, assay = "RNA") + RotatedAxis()
dev.off()

### Other genes of interest ###
genes_tfs <- c("IRF4", "PRDM1", "BCL6", "IL10",
               "XBP1", "ZBTB20", "PAX5", "CD40", 
               "IL4", "MITF", "MTA3", "BACH2", 
               "CXCR5", "CD19","CR2", "FCER2", 
               "IL6R", "CXCR6", "CXCR4", "CXCR3",
               "CD38", "SDC1") # PRDM1 = BLIMP, IL10 = AID, CR2 = CD21, FCER2 = CD23

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Plots")
pdf(file = "CD19_pos_Dotplot.pdf", width = 12, height = 10)
DotPlot(ascs_cd19_only, features = HomingMarkers, group.by = "Isotype") + RotatedAxis()
DotPlot(ascs_cd19_only, features = markerGenes, group.by = "Isotype") + RotatedAxis()
DotPlot(ascs_cd19_only, features = genes_tfs, group.by = "Isotype") + RotatedAxis()
dev.off()

pdf(file = "CD19andPop5_pos_Dotplot.pdf", width = 12, height = 10)
DotPlot(ascs_unint, features = HomingMarkers, group.by = "Isotype") + RotatedAxis()
DotPlot(ascs_unint, features = markerGenes, group.by = "Isotype") + RotatedAxis()
DotPlot(ascs_unint, features = genes_tfs, group.by = "Isotype") + RotatedAxis()
dev.off()
