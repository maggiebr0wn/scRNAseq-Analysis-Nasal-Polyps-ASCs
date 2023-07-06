#!/usr/bin/env Rscript

library(pheatmap)
library("rfUtilities")
library(Seurat)
library("zoo")

### 8-19-2020 Redo Nasal Polyps ###
# 1. Look at each patient individually
# 2. Assess mitochondrial content, gene count, % IgGenes
# 3. Integrated vs Not Integrated
# 4. Assess whether GEX can be used for missing VDJ isotype labels:
#     - LogNorm vs Scaled GEX

##############
##### QC #####
##############

### CLS063 CD19+ + Pop5; exclude VDJ for now
# load, create seurat object, merge
cls063ascs.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS063_ASCs_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls063pop5.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS063_pop5_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls063ascs <- CreateSeuratObject(counts = cls063ascs.data, project = "cls063_ascs", min.cells = 3, min.features = 0)
cls063pop5 <- CreateSeuratObject(counts = cls063pop5.data, project = "cls063_pop5", min.cells = 3, min.features = 0)

### CLS071 CD19+ + Pop5
# load, create seurat object
cls071ascs.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS071_ASCs_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls071pop5.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS071_pop5_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls071ascs <- CreateSeuratObject(counts = cls071ascs.data, project = "cls071_ascs", min.cells = 3, min.features = 0)
cls071pop5 <- CreateSeuratObject(counts = cls071pop5.data, project = "cls071_pop5", min.cells = 3, min.features = 0)

### CLS072 CD19+ + Pop5
cls072ascs.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS072_ASCs_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls072pop5.data <- Read10X_h5(filename = "/Users/maggiebrown/scProject/nasal_polyps/nasal_polyps/Nasal/CLS072_pop5_GEX/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
cls072ascs <- CreateSeuratObject(counts = cls072ascs.data, project = "cls072_ascs", min.cells = 3, min.features = 0)
cls072pop5 <- CreateSeuratObject(counts = cls072pop5.data, project = "cls072_pop5", min.cells = 3, min.features = 0)

### Merge all donors, assess QC together
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/QC/Plots")

merge_all <- merge(cls063ascs, y = c(cls063pop5, cls071ascs, cls071pop5, cls072ascs, cls072pop5), 
                   add.cell.ids = c("cls063ascs", "cls063pop5", "cls071ascs", "cls071pop5", "cls072ascs", "cls072pop5"), 
                   project = "nasal")
# assess
merge_all[["percent.mt"]] <- PercentageFeatureSet(merge_all, pattern = "^MT-")
merge_all[["percent.IgGenes"]] <- PercentageFeatureSet(merge_all, pattern = "^IGH|^IGL|^IGK")

pdf(file = "violin_plot.pdf")
VlnPlot(merge_all, features = "nFeature_RNA")
VlnPlot(merge_all, features = "nCount_RNA")
VlnPlot(merge_all, features = "percent.mt")
VlnPlot(merge_all, features = "percent.IgGenes")
dev.off()

plot1 <- FeatureScatter(merge_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merge_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(merge_all, feature1 = "nCount_RNA", feature2 = "percent.IgGenes")
pdf(file = "scatter_plots.pdf")
plot1
plot2
plot3
dev.off()

###################################################################################
### Assign cells as: low/mid/hi %mt, low/mid/hi gene_count, low/mid/hi %IgGenes ###
###################################################################################

# Group %mt
merge_all[['mtGroup']] = "NA"
merge_all$mtGroup[which(merge_all$percent.mt <= 10)] = "LowMT"
merge_all$mtGroup[which(merge_all$percent.mt > 10 & merge_all$percent.mt <= 50)] = "MidMT"
merge_all$mtGroup[which(merge_all$percent.mt > 50)] = "HighMT"
table(merge_all$mtGroup)
# HighMT  LowMT  MidMT 
# 227  14611   1104

# Group genecount
merge_all[['GNGroup']] = "NA"
merge_all$GNGroup[which(merge_all$nFeature_RNA <= 1000)] = "LowGN"
merge_all$GNGroup[which(merge_all$nFeature_RNA > 1000 & merge_all$nFeature_RNA <= 4000)] = "MidGN"
merge_all$GNGroup[which(merge_all$nFeature_RNA > 4000)] = "HighGN"
table(merge_all$GNGroup)
# HighGN  LowGN  MidGN 
# 993   8035   6914

# Group %IgGenes
merge_all[['IgGroup']] = "NA"
merge_all$IgGroup[which(merge_all$percent.IgGenes <= 25)] = "LowIg"
merge_all$IgGroup[which(merge_all$percent.IgGenes > 25 & merge_all$percent.IgGenes <= 75)] = "MidIg"
merge_all$IgGroup[which(merge_all$percent.IgGenes > 75)] = "HighIg"
table(merge_all$IgGroup)
# HighIg  LowIg  MidIg 
# 1517   2305  12120 

saveRDS(merge_all, "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Merged_Unfiltered.rds")

#########################################
### Detect/Remove Contaminating Cells ###
#########################################

merge_all <- readRDS("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Merged_Unfiltered.rds")

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Contaminating_Cells")
markerGenes <- c(
  "MS4A1", # non-ASC B Cells
  "JCHAIN", # ASCs
  "SDC1", # CD138 - for ASC subpops ID
  "CD19", # for ASC subpops ID
  "CD38", # for ASC subpops ID
  "IL7R", "S100A4", # Memory CD4+ T Cells
  "CD8A", # CD8+ T Cells
  "FCER1A", "CST3", # Dendritic Cells
  "FCGR3A", "MS4A7", # FCGR3A+ Monocytes
  "CD14", "LYZ", # CD14+ Monocytes
  "GNLY", "NKG7", # NK Cells
  "PPBP" # Platelets
)

remove_contam <- merge_all
# Quick and Dirt Clustering
remove_contam <- NormalizeData(remove_contam, normalization.method = "LogNormalize", scale.factor = 10000)
remove_contam <- FindVariableFeatures(remove_contam, selection.method = "vst", nfeatures = 2750)
all.genes <- rownames(remove_contam)
remove_contam <- ScaleData(remove_contam, features = all.genes)
remove_contam <- RunPCA(remove_contam, features = VariableFeatures(object = remove_contam))
# PCA Plots
remove_contam <- SetIdent(remove_contam, value = remove_contam$orig.ident)
pca_orig <- DimPlot(remove_contam, reduction = "pca")

remove_contam <- SetIdent(remove_contam, value = remove_contam$IgGroup)
pca_Ig <- DimPlot(remove_contam, reduction = "pca")

remove_contam <- SetIdent(remove_contam, value = remove_contam$mtGroup)
pca_mt <- DimPlot(remove_contam, reduction = "pca")

remove_contam <- SetIdent(remove_contam, value = remove_contam$GNGroup)
pca_GN <- DimPlot(remove_contam, reduction = "pca")

pdf(file = "PCA_plots.pdf")
pca_orig
pca_Ig
pca_mt
pca_GN
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(remove_contam) # Choose 15 PCs
dev.off()

# Cluster
remove_contam <- FindNeighbors(remove_contam, dims = 1:15)
remove_contam <- FindClusters(remove_contam, resolution = 0.5)
remove_contam <- RunUMAP(remove_contam, dims = 1:15)

# UMAP Plots
remove_contam <- SetIdent(remove_contam, value = remove_contam$seurat_clusters)
umap_clusts <- DimPlot(remove_contam, reduction = "umap", label = TRUE)

remove_contam <- SetIdent(remove_contam, value = remove_contam$orig.ident)
umap_orig <- DimPlot(remove_contam, reduction = "umap", label = TRUE)

remove_contam <- SetIdent(remove_contam, value = remove_contam$IgGroup)
umap_Ig <- DimPlot(remove_contam, reduction = "umap", label = TRUE)

remove_contam <- SetIdent(remove_contam, value = remove_contam$mtGroup)
umap_mt <- DimPlot(remove_contam, reduction = "umap", label = TRUE)

remove_contam <- SetIdent(remove_contam, value = remove_contam$GNGroup)
umap_GN <- DimPlot(remove_contam, reduction = "umap", label = TRUE)

pdf(file = "UMAP_plots_QCs.pdf")
umap_clusts
umap_orig
umap_Ig
umap_mt
umap_GN
dev.off()

# ID Contaminating Cells
remove_contam <- SetIdent(remove_contam, value = remove_contam$seurat_clusters)
marker_genes_featplot <- FeaturePlot(remove_contam, features = markerGenes, )
marker_genes_dotplot <- DotPlot(remove_contam, features = markerGenes) + RotatedAxis()
pdf(file = "markerGene_plots.pdf", width = 10, height = 10)
marker_genes_featplot
marker_genes_dotplot
dev.off()

### Contaminating Cells to Remove:
# Clusters: 9, 12, 13, 14; T Cells, Dendritic Cells, NK Cells
# 9 = 333 cells, 12 = 123 cells, 13 = 109 cells, 14 = 63 cells
remove_contam <- SetIdent(remove_contam, value = remove_contam$seurat_clusters)
# > remove_contam
# An object of class Seurat 
# 18953 features across 15942 samples within 1 assay 
# Active assay: RNA (18953 features, 2750 variable features)
# 2 dimensional reductions calculated: pca, umap

pure_ascs <- subset(remove_contam, 
                    idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "11"))

saveRDS(pure_ascs, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_ASCs.rds")
#########################
### PURE ASCs Metrics ###
#########################
# > pure_ascs
# An object of class Seurat 
# 18953 features across 15314 samples within 1 assay 
# Active assay: RNA (18953 features, 2750 variable features)
# 2 dimensional reductions calculated: pca, umap
#
# > table(pure_ascs$orig.ident)
# cls063_ascs cls063_pop5 cls071_ascs cls071_pop5 cls072_ascs cls072_pop5 
# 6105        1501        3815        2149        1392         352 
# 
# > table(pure_ascs$IgGroup)
# HighIg  LowIg  MidIg 
# 1517   1723  12074 
# 
# > table(pure_ascs$mtGroup)
# HighMT  LowMT  MidMT 
# 221  14057   1036 
# 
# > table(pure_ascs$GNGroup)
# HighGN  LowGN  MidGN 
# 956   7765   6593 

#######################################################################
### Redo Clustering Analysis for Pure ASCs; Still Have bad QC Cells ### 
#######################################################################

pure_ascs <- readRDS("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_ASCs.rds")

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Unfiltered_ASCs")
# Quick and Dirt Clustering
pure_ascs <- NormalizeData(pure_ascs, normalization.method = "LogNormalize", scale.factor = 10000)
pure_ascs <- FindVariableFeatures(pure_ascs, selection.method = "vst", nfeatures = 2750)
all.genes <- rownames(pure_ascs)
pure_ascs <- ScaleData(pure_ascs, features = all.genes)
pure_ascs <- RunPCA(pure_ascs, features = VariableFeatures(object = pure_ascs))
# PCA Plots
pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$orig.ident)
pca_orig <- DimPlot(pure_ascs, reduction = "pca")

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$IgGroup)
pca_Ig <- DimPlot(pure_ascs, reduction = "pca")

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$mtGroup)
pca_mt <- DimPlot(pure_ascs, reduction = "pca")

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$GNGroup)
pca_GN <- DimPlot(pure_ascs, reduction = "pca")

pdf(file = "PCA_plots.pdf")
pca_orig
pca_Ig
pca_mt
pca_GN
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(pure_ascs) # Choose 15 PCs
dev.off()

# Cluster
pure_ascs <- FindNeighbors(pure_ascs, dims = 1:12)
pure_ascs <- FindClusters(pure_ascs, resolution = 0.5)
pure_ascs <- RunUMAP(pure_ascs, dims = 1:12)

# UMAP Plots
pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$seurat_clusters)
umap_clusts <- DimPlot(pure_ascs, reduction = "umap", label = TRUE)

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$orig.ident)
umap_orig <- DimPlot(pure_ascs, reduction = "umap", label = TRUE)

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$IgGroup)
umap_Ig <- DimPlot(pure_ascs, reduction = "umap", label = TRUE)

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$mtGroup)
umap_mt <- DimPlot(pure_ascs, reduction = "umap", label = TRUE)

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$GNGroup)
umap_GN <- DimPlot(pure_ascs, reduction = "umap", label = TRUE)

pdf(file = "UMAP_plots_QCs.pdf")
umap_clusts
umap_orig
umap_Ig
umap_mt
umap_GN
dev.off()

markerGenes <- c(
  "MS4A1", # non-ASC B Cells
  "JCHAIN", # ASCs
  "SDC1", # CD138 - for ASC subpops ID
  "CD19", # for ASC subpops ID
  "CD38" # for ASC subpops ID
)

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$seurat_clusters)
pdf(file = "markerGene_plots.pdf", width = 10, height = 10)
FeaturePlot(pure_ascs, features = markerGenes)
DotPlot(pure_ascs, features = markerGenes) + RotatedAxis()
dev.off()

################################
### Remove Low Quality Cells ###
################################

# Confirm apoptosis in Mid-Hight MT Cells
pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$mtGroup)

anti_apop_genes <- c("BCL2", "MCL1", "BCL2L1", "BCL2L2", "BCL2A1")
anti_featplot <- FeaturePlot(pure_ascs, features = anti_apop_genes)
anti_dotplot <- DotPlot(pure_ascs, features = anti_apop_genes)

pro_apop_genes <- c("P53", "BAX", "BAK1", "BOK", "BID", "BMF")
pro_featplot <- FeaturePlot(pure_ascs, features = pro_apop_genes)
pro_dotplot <- DotPlot(pure_ascs, features = pro_apop_genes)

mt_umap_split <- DimPlot(pure_ascs, reduction = "umap", split.by = "mtGroup")

pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$GNGroup)
gn_umap_split <- DimPlot(pure_ascs, reduction = "umap", split.by = "GNGroup")

pdf(file = "mt_Assess_plots.pdf", width = 10, height = 10)
mt_umap_split
anti_featplot
anti_dotplot
pro_featplot
pro_dotplot
gn_umap_split
dev.off()

# Hight MT cells not showing evidence for anti-apoptosis or pro-apoptosis; low GN

#         HighGN LowGN MidGN
# HighMT      0   221     0
# LowMT     956  6588  6513
# MidMT       0   956    80

### Remove Bad Quality Cells ###
pure_ascs <- SetIdent(pure_ascs, value = pure_ascs$orig.ident)
# assess
pure_ascs[["percent.mt"]] <- PercentageFeatureSet(pure_ascs, pattern = "^MT-")
pure_ascs[["percent.IgGenes"]] <- PercentageFeatureSet(pure_ascs, pattern = "^IGH|^IGL|^IGK")

pdf(file = "violin_plot.pdf")
VlnPlot(pure_ascs, features = "nFeature_RNA")
VlnPlot(pure_ascs, features = "nCount_RNA")
VlnPlot(pure_ascs, features = "percent.mt")
VlnPlot(pure_ascs, features = "percent.IgGenes")
dev.off()

plot1 <- FeatureScatter(pure_ascs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pure_ascs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(pure_ascs, feature1 = "nCount_RNA", feature2 = "percent.IgGenes")
pdf(file = "scatter_plots.pdf")
plot1
plot2
plot3
dev.off()

### Cross Check Tables ###

# 1. GN vs MT
table(pure_ascs$GNGroup, pure_ascs$mtGroup)
#         HighMT LowMT MidMT
# HighGN      0   956     0
# LowGN     221  6588   956   
# MidGN       0  6513    80

# 2. GN vs IgGenes
table(pure_ascs$GNGroup, pure_ascs$IgGroup)
#         HighIg LowIg MidIg
# HighGN      0   266   690
# LowGN     810  1089  5866
# MidGN     707   368  5518

# 3. MT vs IgGenes
table(pure_ascs$mtGroup, pure_ascs$IgGroup)
#         HighIg LowIg MidIg
# HighMT      0   210    11
# LowMT    1515  1098 11444
# MidMT       2   415   619

### Groups to remove:
# Hight MT, Low GN, Low IgGenes

## New Thresholds:
# Mit: x < 25
# GN: 1000 < x < 5000
# IgGenes: x > 25

filtered_ascs <- subset(pure_ascs, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 25 & percent.IgGenes > 25)
# 18953 features across 6450 samples within 1 assay 

saveRDS(filtered_ascs, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")

##########################################
### Redo Analysis: Filtered, Pure ASCs ###
##########################################

filtered_ascs <- readRDS("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Unintegrated")
# Re-Normalize, Scale, and PCA
filtered_ascs <- NormalizeData(filtered_ascs, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_ascs <- FindVariableFeatures(filtered_ascs, selection.method = "vst", nfeatures = 2750)
all.genes <- rownames(filtered_ascs)
filtered_ascs <- ScaleData(filtered_ascs, features = all.genes)
filtered_ascs <- RunPCA(filtered_ascs, features = VariableFeatures(object = filtered_ascs))
# PCA Plots
filtered_ascs <- SetIdent(filtered_ascs, value = filtered_ascs$orig.ident)
pca_orig <- DimPlot(filtered_ascs, reduction = "pca")
pdf(file = "PCA_plots.pdf")
pca_orig
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(filtered_ascs) # Choose 13 PCs
dev.off()

# Cluster
filtered_ascs <- FindNeighbors(filtered_ascs, dims = 2:15)
filtered_ascs <- FindClusters(filtered_ascs, resolution = 0.5)
filtered_ascs <- RunUMAP(filtered_ascs, dims = 2:15)

# UMAP Plots
filtered_ascs <- SetIdent(filtered_ascs, value = filtered_ascs$seurat_clusters)
umap_clusts <- DimPlot(filtered_ascs, reduction = "umap", label = TRUE)

filtered_ascs <- SetIdent(filtered_ascs, value = filtered_ascs$orig.ident)
umap_orig <- DimPlot(filtered_ascs, reduction = "umap", label = TRUE)

umap_split <- DimPlot(filtered_ascs, reduction = "umap", split.by = "orig.ident")

pdf(file = "UMAP_plots.pdf")
umap_clusts
umap_orig
dev.off()

pdf(file = "UMAP_plot_split.pdf", width = 15, height = 5)
umap_split
dev.off()

markerGenes <- c(
  "MS4A1", # non-ASC B Cells
  "JCHAIN", # ASCs
  "SDC1", # CD138 - for ASC subpops ID
  "CD19", # for ASC subpops ID
  "CD38" # for ASC subpops ID
)

pdf(file = "ASCs_subpop_markers.pdf", width = 10, height = 10)
FeaturePlot(filtered_ascs, features = markerGenes)
DotPlot(filtered_ascs, features = markerGenes) + RotatedAxis()
dev.off()

# ID Isotypes w/o removing cont cells
filtered_ascs <- SetIdent(filtered_ascs, value = filtered_ascs$seurat_clusters)
isotype_genes <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4",
                   "IGHA1", "IGHA2",
                   "IGHD",
                   "IGHM",
                   "IGHE")

pdf(file = "Isotype_markers.pdf", width = 10, height = 10)
FeaturePlot(filtered_ascs, features = isotype_genes)
DotPlot(filtered_ascs, features = isotype_genes) + RotatedAxis()
dev.off()

# Group genecount
merge_all[['GNGroup']] = "NA"
merge_all$GNGroup[which(merge_all$nFeature_RNA <= 1000)] = "LowGN"
merge_all$GNGroup[which(merge_all$nFeature_RNA > 1000 & merge_all$nFeature_RNA <= 4000)] = "MidGN"
merge_all$GNGroup[which(merge_all$nFeature_RNA > 4000)] = "HighGN"
table(merge_all$GNGroup)

filtered_ascs <- SetIdent(filtered_ascs, value = filtered_ascs$GNGroup)
DimPlot(filtered_ascs, reduction = "umap")

# See how many cells express min for gene
# subset(filtered_ascs, subset = IGHE > 1)
# subset(filtered_ascs, subset = IGHA1 > 1)
# subset(filtered_ascs, subset = IGHA2 > 1)
# subset(filtered_ascs, subset = IGHG1 > 1)
# subset(filtered_ascs, subset = IGHG2 > 1)
# subset(filtered_ascs, subset = IGHG3 > 1)
# subset(filtered_ascs, subset = IGHG4 > 1)
# subset(filtered_ascs, subset = IGHD > 1)
# subset(filtered_ascs, subset = IGHM > 1)

### Look at IGH Genes using LogNormalized counts ###
filtered_ascs_lognorm <- filtered_ascs
# Get isotype gene expressions
for (i in isotype_genes){
  print(i)
  filtered_ascs_lognorm[[i]] = GetAssayData(object = filtered_ascs[["RNA"]], slot = "data")[i,]
}

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Unintegrated/Isotype_Assignments_GE/LogNormalized_counts")

df <- data.frame(filtered_ascs_lognorm$IGHA1, filtered_ascs_lognorm$IGHA2,
                 filtered_ascs_lognorm$IGHD,
                 filtered_ascs_lognorm$IGHE,
                 filtered_ascs_lognorm$IGHG1, filtered_ascs_lognorm$IGHG2, filtered_ascs_lognorm$IGHG3, filtered_ascs_lognorm$IGHG4,
                 filtered_ascs_lognorm$IGHM)

write.table(df, "log_normalized_isotype_genes.txt", quote = FALSE, sep = "\t")

# use isotype_genes.txt with Isotype_Assignments_GE/assign_isotypes.py to get assignments.txt

lognorm_assignments <- read.delim("assignments.txt", header=FALSE)

filtered_ascs_lognorm[["Isotype"]] = lognorm_assignments$V2

filtered_ascs_lognorm <- SetIdent(filtered_ascs_lognorm, value = filtered_ascs_lognorm$Isotype)

pdf(file = "Isotype_assignments_PCA_UMAP.pdf", width = 7, height = 7)
DimPlot(filtered_ascs_lognorm, reduction = "pca")
DimPlot(filtered_ascs_lognorm, reduction = "umap")
dev.off()

pdf(file = "Isotype_assignments_splitUMAP.pdf", width = 10, height = 3)
DimPlot(filtered_ascs_lognorm, reduction = "umap", split.by = "Isotype")
dev.off()

pdf(file = "heatmap_isotype.pdf", width = 15, height = 4)
DoHeatmap(filtered_ascs_lognorm, draw.lines = FALSE, features = isotype_genes, angle = 90, size = 2)
dev.off()

# LogNorm GEX Predictions
# > table(filtered_ascs_lognorm$Isotype)
# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4 
# 1687 2035    3  123 1060 1358   80  104

genes_tfs <- c("IRF4", "PRDM1", "BCL6", "IL10",
               "XBP1", "ZBTB20", "PAX5", "CD40", 
               "IL4", "MITF", "MTA3", "BACH2") # PRDM1 = BLIMP, IL10 = AID
pdf(file = "Bcell_genes_dotplot_isotype.pdf")
DotPlot(filtered_ascs_lognorm, features = genes_tfs, assay = "RNA") + RotatedAxis()
DotPlot(filtered_ascs_lognorm, features = isotype_genes, assay = "RNA") + RotatedAxis()
dev.off()

### Look at IGH Genes using LogNormalized and Scaled counts ###
filtered_ascs_scaled <- filtered_ascs
# Get isotype gene expressions
for (i in isotype_genes){
  print(i)
  filtered_ascs_scaled[[i]] = GetAssayData(object = filtered_ascs[["RNA"]], slot = "scale.data")[i,]
}

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Unintegrated/Isotype_Assignments_GE/Scaled_Counts")

df <- data.frame(filtered_ascs_scaled$IGHA1, filtered_ascs_scaled$IGHA2,
                 filtered_ascs_scaled$IGHD,
                 filtered_ascs_scaled$IGHE,
                 filtered_ascs_scaled$IGHG1, filtered_ascs_scaled$IGHG2, filtered_ascs_scaled$IGHG3, filtered_ascs_scaled$IGHG4,
                 filtered_ascs_scaled$IGHM)

write.table(df, "scaled_isotype_genes.txt", quote = FALSE, sep = "\t")

# use isotype_genes.txt with assign_isotypes.py to get assignments.txt

scaled_assignments <- read.delim("assignments.txt", header=FALSE)

filtered_ascs_scaled[["Isotype"]] = scaled_assignments$V2

filtered_ascs_scaled <- SetIdent(filtered_ascs_scaled, value = filtered_ascs_scaled$Isotype)

pdf(file = "Isotype_assignments_PCA_UMAP.pdf", width = 7, height = 7)
DimPlot(filtered_ascs_scaled, reduction = "pca")
DimPlot(filtered_ascs_scaled, reduction = "umap")
dev.off()

pdf(file = "Isotype_assignments_splitUMAP.pdf", width = 10, height = 3)
DimPlot(filtered_ascs_scaled, reduction = "umap", split.by = "Isotype")
dev.off()

pdf(file = "heatmap_isotype.pdf", width = 15, height = 4)
DoHeatmap(filtered_ascs_scaled, draw.lines = FALSE, features = isotype_genes, angle = 90, size = 2)
dev.off()

# table(filtered_ascs_scaled$Isotype)
# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4 
# 1295 1867   73  332  866 1389  434  194

# Heatmap: Dotplot: Scaled Isotypes
filtered_ascs_scaled <- SetIdent(filtered_ascs_scaled, value = filtered_ascs_scaled$Isotype)
pdf(file = "heatmap_isotype.pdf", width = 15, height = 4)
DoHeatmap(filtered_ascs_scaled, draw.lines = FALSE, features = isotype_genes, angle = 90, size = 2)
dev.off()

genes_tfs <- c("IRF4", "PRDM1", "BCL6", "IL10",
               "XBP1", "ZBTB20", "PAX5", "CD40", 
               "IL4", "MITF", "MTA3", "BACH2") # PRDM1 = BLIMP, IL10 = AID
pdf(file = "Bcell_genes_dotplot_isotype.pdf")
DotPlot(filtered_ascs_scaled, features = genes_tfs, assay = "RNA") + RotatedAxis()
DotPlot(filtered_ascs_scaled, features = isotype_genes, assay = "RNA") + RotatedAxis()
dev.off()

#######################################################
### Redo Analysis: Integration, Filtered, Pure ASCs ###
#######################################################

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Filtered_ASCs/Integrated/")
integrate <- filtered_ascs
integrate[["LogNormIsotype"]] = lognorm_assignments$V2
integrate[["ScaledIsotype"]] = scaled_assignments$V2

# Integrate
integrate.list <- SplitObject(integrate, split.by = "orig.ident")
for (i in 1:length(integrate.list)) {
  integrate.list[[i]] <- NormalizeData(integrate.list[[i]], verbose = FALSE)
  integrate.list[[i]] <- FindVariableFeatures(integrate.list[[i]], selection.method = "vst", 
                                              nfeatures = 2750, verbose = FALSE)
}
reference.list <- integrate.list[c("cls063_ascs", "cls063_pop5", "cls071_ascs", "cls071_pop5", "cls072_ascs", "cls072_pop5")]
ascs.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
ascs.integrated <- IntegrateData(anchorset = ascs.anchors, dims = 1:30)
DefaultAssay(ascs.integrated) <- "integrated"

### Recluster
ascs.integrated <- ScaleData(ascs.integrated)
ascs.integrated <- RunPCA(ascs.integrated, npcs = 30)

ascs.integrated[["LogNormIsotype"]] = lognorm_assignments$V2
ascs.integrated[["ScaledIsotype"]] = scaled_assignments$V2

# PCA Plots: LogNormalised Isotypes
ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$orig.ident)
pca_orig <- DimPlot(ascs.integrated, reduction = "pca")

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$LogNormIsotype)
pca_iso <- DimPlot(ascs.integrated, reduction = "pca")

pdf(file = "LogNorm_Isotype_Labels_PCA_plots.pdf")
pca_orig
pca_iso
dev.off()

# PCA Plots: Scaled Isotypes
ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$orig.ident)
pca_orig <- DimPlot(ascs.integrated, reduction = "pca")

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$ScaledIsotype)
pca_iso <- DimPlot(ascs.integrated, reduction = "pca")

pdf(file = "Scaled_Isotype_LabelsPCA_plots.pdf")
pca_orig
pca_iso
dev.off()

# Elbow Plot
pdf(file = "elbow_plot.pdf")
ElbowPlot(ascs.integrated) # Choose 2:12 PCs
dev.off()

# Cluster
ascs.integrated <- FindNeighbors(ascs.integrated, dims = 2:12)
ascs.integrated <- FindClusters(ascs.integrated, resolution = 0.5)
ascs.integrated <- RunUMAP(ascs.integrated, dims = 2:12)

# UMAP Plots; PCA Plots: LogNormalised Isotypes
ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$seurat_clusters)
umap_clusts <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$orig.ident)
umap_orig <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

umap_split <- DimPlot(ascs.integrated, reduction = "umap", split.by = "orig.ident")

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$LogNormIsotype)
umap_iso <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

umap_isosplit <- DimPlot(ascs.integrated, reduction = "umap", split.by = "LogNormIsotype")

pdf(file = "LogNorm_Isotype_Labels_UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
dev.off()

pdf(file = "LogNorm_Isotype_Labels_UMAPs_split.pdf", width = 12, height = 3)
umap_split
umap_isosplit
dev.off()

# UMAP Plots; PCA Plots: Scaled Isotypes
ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$seurat_clusters)
umap_clusts <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$orig.ident)
umap_orig <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

umap_split <- DimPlot(ascs.integrated, reduction = "umap", split.by = "orig.ident")

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$ScaledIsotype)
umap_iso <- DimPlot(ascs.integrated, reduction = "umap", label = TRUE)

umap_isosplit <- DimPlot(ascs.integrated, reduction = "umap", split.by = "ScaledIsotype")

pdf(file = "Scaled_Isotype_Labels_UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
dev.off()

pdf(file = "Scaled_Isotype_Labels_UMAPs_split.pdf", width = 12, height = 3)
umap_split
umap_isosplit
dev.off()

# Plot Features
ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$seurat_clusters)
markerGenes <- c(
  "MS4A1", # non-ASC B Cells
  "JCHAIN", # ASCs
  "SDC1", # CD138 - for ASC subpops ID
  "CD19", # for ASC subpops ID
  "CD38" # for ASC subpops ID
)

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$seurat_clusters)
pdf(file = "ASCs_subpop_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs.integrated, reduction = "umap", features = markerGenes, slot = "data")
DotPlot(ascs.integrated, features = markerGenes, assay = "RNA") + RotatedAxis()
dev.off()

ascs.integrated <- SetIdent(ascs.integrated, value = ascs.integrated$seurat_clusters)
isotype_genes <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4",
                   "IGHA1", "IGHA2",
                   "IGHD",
                   "IGHM",
                   "IGHE")

pdf(file = "Isotype_markers.pdf", width = 10, height = 10)
FeaturePlot(ascs.integrated, features = isotype_genes)
DotPlot(ascs.integrated, features = isotype_genes, assay = "RNA") + RotatedAxis()
dev.off()

###########################
### Save Seurat Objects ###
###########################

# Unintegrated
#saveRDS(filtered_ascs_lognorm, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/filtered_ascs_lognorm.rds")
#saveRDS(filtered_ascs_scaled, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/filtered_ascs_scaled.rds")
unintegrated <- filtered_ascs
unintegrated[["LogNormIsotype"]] = lognorm_assignments$V2
unintegrated[["ScaledIsotype"]] = scaled_assignments$V2

saveRDS(unintegrated, "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")

# Integrated
saveRDS(ascs.integrated, file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Integrated_ASCs.rds")

###################################
### Compare Isotype Assignments ###
###################################

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Compare_Isotype_Labels/LogNormalized/")

write.table(ascs.integrated$LogNormIsotype, file = "LogNorm_Isotypes.txt", quote = FALSE, col.names = FALSE)
vdj_labels <- read.table("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Compare_Isotype_Labels/FIXED_vdj_labels_ordered_for_Seurat.txt", quote="\"", comment.char="")

### Compare VDJ Isotype Labels vs LogNormalized Expression Labels ###

table(filtered_ascs_lognorm$Isotype)

# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4 
# 1687 2035    3  123 1060 1358   80  104

# 18953 features across 6450 samples within 1 assay 

filtered_ascs_lognorm <- AddMetaData(
  object = filtered_ascs_lognorm,
  metadata = vdj_labels$V1,
  col.name = "VDJ_labels"
)

#             IgA1 IgA2 IgD IgE IgG1 IgG2 IgG3 IgG4 IgM IgNONE
# cls063_ascs  270  487   0  47  234  410   35    0  27      1
# cls063_pop5   71  393   0   0   33  153    1    0   3      0
# cls071_ascs  242  162   0  18  198  215   38    1  33      0
# cls071_pop5  129  348   2   0   94  171   17    0  34      0
# cls072_ascs  259  122   0  27  139   51   18    1   9      0
# cls072_pop5   61  129   0   0    5   14    1    0   3      0

cM <- confusionMatrix(paste0(filtered_ascs_lognorm$VDJ_labels), paste0(filtered_ascs_lognorm$Isotype))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

pdf("vdj_vs_expression.pdf")
p
dev.off()
write.table(cM, file = "vdj_vs_expression.txt")

cM <- confusionMatrix(paste0(filtered_ascs_lognorm$Isotype), paste0(filtered_ascs_lognorm$VDJ_labels))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

pdf("expression_vs_vdj.pdf")
p
dev.off()
write.table(cM, file = "expression_vs_vdj.txt")

temp_lognorm <- SetIdent(filtered_ascs_lognorm, value = filtered_ascs_lognorm$VDJ_labels)
# 18953 features across 6450 samples within 1 assay 
temp_lognorm <- subset(temp_lognorm, idents = c("IgA1", "IgA2", "IgD", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))
# 18953 features across 4705 samples within 1 assay 

accuracy(paste0(temp_lognorm$Isotype), paste0(temp_lognorm$VDJ_labels))

# Accuracy (PCC): 94.5377258235919% 
# 
# Cohen's Kappa: 0.9278 
# 
# Users accuracy: 
#  IgA1  IgA2   IgD   IgE  IgG1  IgG2  IgG3  IgG4   IgM 
#  97.8  98.5 100.0  98.9  98.2  98.3  37.3 100.0   0.0 
#  
# 
# Producers accuracy: 
# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4  IgM 
# 92.0 97.9 66.7 89.2 96.5 97.8 83.7  2.9  NaN 
#  
# 
# Confusion matrix 
#       y
# x      IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4  IgM
#   IgA1 1009   13    0    0    4    5    1    0   65
#   IgA2    4 1616    0    0    7    3    0    0   21
#   IgD     0    0    2    0    0    0    0    0    1
#   IgE     0    0    0   91    0    0    0    0   11
#   IgG1    6    4    0    0  690    9    2    0    4
#   IgG2    6    6    0    0    2  997    2    0    6
#   IgG3    6    1    0    0    0    0   41    0    1
#   IgG4    1    1    0    1    0    0   64    2    0
  
### Compare VDJ Isotype Labels vs Scaled Expression Labels ###

table(ascs.integrated$ScaledIsotype)
table(filtered_ascs_scaled$Isotype)

# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4 
# 1295 1867   73  332  866 1389  434  194

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Compare_Isotype_Labels/Scaled")

filtered_ascs_scaled <- AddMetaData(
  object = filtered_ascs_scaled,
  metadata = vdj_labels$V1,
  col.name = "VDJ_labels"
)

cM <- confusionMatrix(paste0(filtered_ascs_scaled$VDJ_labels), paste0(filtered_ascs_scaled$Isotype))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

pdf("vdj_vs_expression.pdf")
p
dev.off()
write.table(cM, file = "vdj_vs_expression.txt")

cM <- confusionMatrix(paste0(filtered_ascs_scaled$Isotype), paste0(filtered_ascs_scaled$VDJ_labels))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
)
p

pdf("expression_vs_vdj.pdf")
p
dev.off()
write.table(cM, file = "expression_vs_vdj.txt")

temp_scaled <- SetIdent(filtered_ascs_scaled, value = filtered_ascs_scaled$VDJ_labels)
# 18953 features across 6450 samples within 1 assay 
temp_scaled <- subset(temp_scaled, idents = c("IgA1", "IgA2", "IgD", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))
# 18953 features across 4705 samples within 1 assay 

accuracy(paste0(temp_scaled$Isotype), paste0(temp_scaled$VDJ_labels))

# Accuracy (PCC): 87.3326248671626% 
# 
# Cohen's Kappa: 0.8355 
# 
# Users accuracy: 
#  IgA1  IgA2   IgD   IgE  IgG1  IgG2  IgG3  IgG4   IgM 
#  91.4  95.6 100.0 100.0  73.1  93.2  40.0  50.0   0.0 
#  
# 
# Producers accuracy: 
# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4  IgM 
# 95.9 97.1  3.6 39.5 94.0 95.9 21.6  1.2  NaN 
#  
# 
# Confusion matrix 
#       y
# x      IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4  IgM
#   IgA1  943   11    0    0    1    2    0    0   26
#   IgA2    6 1568    0    0    3    3    0    0   34
#   IgD    15    7    2    0   11    4    0    0   17
#   IgE    15   17    0   92   55   20   11    1   22
#   IgG1   12   11    0    0  514    7    0    0    3
#   IgG2   15   17    0    0    3  945    2    0    3
#   IgG3   21    8    0    0  116   13   44    0    2
#   IgG4    5    2    0    0    0   20   53    1    2
  
### Compare lognorm vs scaled
accuracy(paste0(filtered_ascs_lognorm$Isotype), paste0(filtered_ascs_scaled$Isotype))
# Accuracy (PCC): 82.6046511627907% 
# 
# Cohen's Kappa: 0.7788 
# 
# Users accuracy: 
# IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4 
# 98.7 97.2  4.1 37.0 83.5 87.9 18.0 44.8 
#  
# 
# Producers accuracy: 
#  IgA1  IgA2   IgD   IgE  IgG1  IgG2  IgG3  IgG4 
#  75.8  89.2 100.0 100.0  68.2  89.9  97.5  83.7 
#  
# 
# Confusion matrix 
#       y
# x      IgA1 IgA2  IgD  IgE IgG1 IgG2 IgG3 IgG4
#   IgA1 1278   49   27   57   99   87   58   32
#   IgA2   13 1815   18   29   44   75   32    9
#   IgD     0    0    3    0    0    0    0    0
#   IgE     0    0    0  123    0    0    0    0
#   IgG1    2    2   15   73  723    6  224   15
#   IgG2    2    1    9   32    0 1221   42   51
#   IgG3    0    0    1    1    0    0   78    0
#   IgG4    0    0    0   17    0    0    0   87

# LogNorm prediction is more accurate, use those labels for missing VDJ data

#####################################
# Add labels to Unintegrated Object #
#####################################

unintegrated <- readRDS("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")

vdj_labels$V1[is.na(vdj_labels$V1)] <- "None"

unintegrated <- AddMetaData(
  object = unintegrated,
  metadata = vdj_labels$V1,
  col.name = "VDJ_labels"
)

### VDJ + High Confidence ###

# Get cells with VDJ Labels

# > table(unintegrated$orig.ident, unintegrated$VDJ_labels)
#             IgA1 IgA2 IgD IgE IgG1 IgG2 IgG3 IgG4 IgM IgNONE None
# cls063_ascs  270  487   0  47  234  410   35    0  27      1  627
# cls063_pop5   71  393   0   0   33  153    1    0   3      0  128
# cls071_ascs  242  162   0  18  198  215   38    1  33      0  454
# cls071_pop5  129  348   2   0   94  171   17    0  34      0  258
# cls072_ascs  259  122   0  27  139   51   18    1   9      0  241
# cls072_pop5   61  129   0   0    5   14    1    0   3      0   36

unintegrated <- SetIdent(unintegrated, value = "VDJ_labels")
# 18953 features across 6450 samples within 1 assay 
unintegrated_high_confidence <- subset(unintegrated, idents = c("IgG2", "IgA2", "IgA1",
                                                                 "IgG1", "IgG3", "IgE",   
                                                                 "IgM", "IgG4", "IgD"))

# 18953 features across 4705 samples within 1 assay 

# For VDJ labeled cells, keep this label
unintegrated_high_confidence <- AddMetaData(
  object = unintegrated_high_confidence,
  metadata = unintegrated_high_confidence$VDJ_labels,
  col.name = "Isotype"
)

unintegrated_high_confidence <- AddMetaData(
  object = unintegrated_high_confidence,
  metadata = "VDJ",
  col.name = "Label_Source"
)

isotype_genes <- c("IGHA1", "IGHA2", "IGHD", "IGHE","IGHG1", "IGHG2","IGHG3","IGHG4","IGHM")

new_object <- unintegrated
table(Idents(new_object))
# IgG2 IgA2 IgA1 IgG1 IgG3  IgE  IgM IgG4  IgD 
# 1014 1641 1032  703  110   92  109    2    2 

new_object <- subset(new_object, ident = "None") # only get cells w/o VDJ labels
# 18953 features across 1744 samples within 1 assay 

# For cells w/o VDJ label, use LogNorm label
new_object <- AddMetaData(
  object = new_object,
  metadata = new_object$LogNormIsotype,
  col.name = "Isotype"
)

new_object <- AddMetaData(
  object = new_object,
  metadata = "RNA",
  col.name = "Label_Source"
)

# merge cells (with VDJ Labels) and (w/o VDJ + High Confidence)
best_ascs <- merge(unintegrated_high_confidence, y = new_object, add.cell.ids = c("vdj_labels", "predicted_labels"), project = "nasal_best")
# 18953 features across 6449 samples within 1 assay; dropped IgNone 

# table(best_ascs$orig.ident, best_ascs$Isotype)
#             IgA1 IgA2 IgD IgE IgG1 IgG2 IgG3 IgG4 IgM
# cls063_ascs  417  625   0  56  362  580   49   21  27
# cls063_pop5   92  443   0   0   53  188    2    1   3
# cls071_ascs  434  233   0  26  293  285   46   11  33
# cls071_pop5  215  424   2   0  140  215   22    1  34
# cls072_ascs  388  153   0  31  195   67   21    3   9
# cls072_pop5   76  147   0   0    5   17    1    0   3

# table(best_ascs$Isotype, best_ascs$Label_Source)
#       RNA  VDJ
# IgA1  590 1032
# IgA2  384 1641
# IgD     0    2
# IgE    21   92
# IgG1  345  703
# IgG2  338 1014
# IgG3   31  110
# IgG4   35    2
# IgM     0  109


best_ascs <- SetIdent(best_ascs, value = best_ascs$Isotype)
pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/Pure_ASCs/Compare_Isotype_Labels/Final_Isotype_Labels.pdf")
levels(best_ascs) <- c("IgM", "IgG4","IgG3",  "IgG2", "IgG1", "IgE", "IgD", "IgA2", "IgA1")
DotPlot(best_ascs, features = genes_tfs, assay = "RNA") + RotatedAxis()
DotPlot(best_ascs, features = isotype_genes, assay = "RNA") + RotatedAxis()
VlnPlot(best_ascs, features = isotype_genes, sort = TRUE)
dev.off()

saveRDS(best_ascs, "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Pure_Filtered_ASCs.rds")

### Redo Clustering with Nasal_Polyps_Downstream.R
