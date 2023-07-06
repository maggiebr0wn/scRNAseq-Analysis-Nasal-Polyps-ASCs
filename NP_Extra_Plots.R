#!/usr/bin/env Rscript

library(pheatmap)
library("rfUtilities")
library(Seurat)
library("zoo")
library(ggplot2)
library(readr)
library(dplyr)

# Extra plots:

# 3-14-2021: Percent Stacked Bar Chart: Isotypes per Cluster (integrated clusters)
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
library(ggplot2)

# create a dataframe
data <- data.frame(ascs_int$Isotype,ascs_int$seurat_clusters,1)

# Stacked + percent
p <- ggplot(data, aes(fill=ascs_int.Isotype, y=X1, x=ascs_int.seurat_clusters)) + 
  ggtitle("Proportion of Isotypes per Cluster") +
  xlab("Cluster Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual("legend", values = c("IgA1" = "gold", "IgA2" = "orange1", 
                                         "IgD" = "indianred", "IgE" = "palegreen4",
                                         "IgG1" = "lightskyblue", "IgG2" = "steelblue3", 
                                         "IgG3" = "royalblue","IgG4" = "royalblue4", 
                                         "IgM" = "purple3"))

p <- p + guides(fill=guide_legend(title="Isotypes"))

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
ggsave(filename = "Isotype_Proportions_Integrated_Clusters.jpeg",
       plot = print(p))

# 3-15-2021: subset IgE cells; PCA & UMAP
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

ige_only <- subset(ascs_int, idents = "IgE")

### Recluster

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/IgE_Subset/")

ige_only <- ScaleData(ige_only)
ige_only <- RunPCA(ige_only, npcs = 30)

# PCA Plots:
ige_only <- SetIdent(ige_only, value = ige_only$orig.ident)
pca_orig <- DimPlot(ige_only, reduction = "pca")

ige_only <- SetIdent(ige_only, value = ige_only$Isotype)
pca_iso <- DimPlot(ige_only, reduction = "pca")

ige_only <- SetIdent(ige_only, value = ige_only$Label_Source)
pca_lab <- DimPlot(ige_only, reduction = "pca")

ige_only <- SetIdent(ige_only, value = ige_only$mtGroup)
pca_mt <- DimPlot(ige_only, reduction = "pca")

ige_only <- SetIdent(ige_only, value = ige_only$GNGroup)
pca_gn <- DimPlot(ige_only, reduction = "pca")

ige_only <- SetIdent(ige_only, value = ige_only$IgGroup)
pca_ig <- DimPlot(ige_only, reduction = "pca")

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
ElbowPlot(ige_only) # Choose 2:12 PCs
dev.off()

# Cluster 0.75 res
ige_only <- FindNeighbors(ige_only, dims = 2:12)
ige_only <- FindClusters(ige_only, resolution = 0.75)
ige_only <- RunUMAP(ige_only, dims = 2:12)

# UMAP Plots
ige_only <- SetIdent(ige_only, value = ige_only$seurat_clusters)
umap_clusts <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$orig.ident)
umap_orig <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$Isotype)
umap_iso <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$Label_Source)
umap_lab <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$mtGroup)
umap_mt <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$GNGroup)
umap_gn <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$IgGroup)
umap_ig <- DimPlot(ige_only, reduction = "umap", label = TRUE)

pdf(file = "res0.75_UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
umap_lab
umap_mt
umap_gn
umap_ig
dev.off()

# Cluster 0.5 res
ige_only <- FindNeighbors(ige_only, dims = 2:12)
ige_only <- FindClusters(ige_only, resolution = 0.5)
ige_only <- RunUMAP(ige_only, dims = 2:12)

# UMAP Plots
ige_only <- SetIdent(ige_only, value = ige_only$seurat_clusters)
umap_clusts <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$orig.ident)
umap_orig <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$Isotype)
umap_iso <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$Label_Source)
umap_lab <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$mtGroup)
umap_mt <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$GNGroup)
umap_gn <- DimPlot(ige_only, reduction = "umap", label = TRUE)

ige_only <- SetIdent(ige_only, value = ige_only$IgGroup)
umap_ig <- DimPlot(ige_only, reduction = "umap", label = TRUE)

pdf(file = "res0.5_UMAPs.pdf")
umap_clusts
umap_orig
umap_iso
umap_lab
umap_mt
umap_gn
umap_ig
dev.off()

# Confirm ASCs:
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

ASC_markers = c("JCHAIN", "CD38", "SDC1", "IRF4", "PRDM1", "XBP1", "CD19", "PAX5", "MS4A1")

DefaultAssay(ascs_int) <- "RNA"

pdf(file = "ASC_Markers.pdf")
DotPlot(ascs_int, features = ASC_markers, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(ascs_int, features = ASC_markers, assay = "RNA", group.by = "Isotype") + RotatedAxis()
FeaturePlot(ascs_int, features = ASC_markers)
dev.off()

# Volcano Plot for IgE vs CD19+ Only DEGs

library(EnhancedVolcano)
ascs_int <- SetIdent(ascs_int, value = "Isotype")
degs <- FindMarkers(ascs_int, ident.1 = "IgE")

pdf(file = "DEGs_Volcano.pdf")
EnhancedVolcano(degs,
                lab = rownames(degs),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 0.25,
                pCutoff = 0.05,
                subtitle = "IgE ASCs",
                legendPosition = "bottom"
                )
dev.off()

### 3-17-2021 
# Assess Cytokine Genes
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")

cytokines_geneset <- read.table("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/cytokines_geneset.txt", quote="\"", comment.char="")
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

pdf(file = "cytokine_genes.pdf", width = 20, height = 4)
DotPlot(ascs_int, features = cytokines_geneset$V1, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis() + theme(axis.text.x = element_text(size = 5))
DotPlot(ascs_int, features = cytokines_geneset$V1,, assay = "RNA", group.by = "Isotype") + RotatedAxis() + theme(axis.text.x = element_text(size = 5))
dev.off()

# Assess Chemokine Genes
chemokine_geneset <- read.table("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/chemokine_geneset.txt", quote="\"", comment.char="")

pdf(file = "chemokine_genes.pdf", width = 20, height = 4)
DotPlot(ascs_int, features = chemokine_geneset$V1, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis() + theme(axis.text.x = element_text(size = 5))
DotPlot(ascs_int, features = chemokine_geneset$V1,, assay = "RNA", group.by = "Isotype") + RotatedAxis() + theme(axis.text.x = element_text(size = 5))
dev.off()

# 3-19-2020 Confirm ASCs:
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

DefaultAssay(ascs_int) <- "RNA"

#pdf(file = "ASC_Markers.pdf")
DotPlot(ascs_int, features = ASC_markers, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(ascs_int, features = ASC_markers, assay = "RNA", group.by = "Isotype") + RotatedAxis()
FeaturePlot(ascs_int, features = "CCR2")
#dev.off()

x <- FindMarkers(
  ascs_int,
  ident.1 = "IgE",
  logfc.threshold = 0.05,
  min.pct = 0.05
)

# 3-20-2021

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

DefaultAssay(ascs_int) <- "RNA"

dn2_prosurvival <- c("LAPTM5", "B2M", "HLA-DRA", "HLA-DRB1", "CD74",
        "PECAM1", "CXCR4", "TNFRSF13C")

#pdf(file = "ASC_Markers.pdf")
DotPlot(ascs_int, features = dn2_prosurvival, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(ascs_int, features = dn2_prosurvival, assay = "RNA", group.by = "Isotype") + RotatedAxis()

# 3-29-2021 Antigen presenting

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

DefaultAssay(ascs_int) <- "RNA"

ap_genes <- c("HLA-A", "HLA-B", "B2M", "HLA-E", "HLA-F", # MHC I
             "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "AP1S2", "CD74", "CTSS" # MHC II
             )

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
pdf(file = "MHC_degs.pdf")
DotPlot(ascs_int, features = ap_genes, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(ascs_int, features = ap_genes, assay = "RNA", group.by = "Isotype") + RotatedAxis()
dev.off()

# IFN-gamma + LAPTM5 and SARAF

ifn_misc <- c("HLA-A", "HLA-B", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1",
              "HLA-DRA", "HLA-DRB1", "B2M", "HLA-E", "HLA-F", "ISG20", "PTPN6",
              "LAPTM5", "SARAF")

pdf(file = "IFN-gama_misc_degs.pdf")
DotPlot(ascs_int, features = ifn_misc, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(ascs_int, features = ifn_misc, assay = "RNA", group.by = "Isotype") + RotatedAxis()
dev.off()

# Highlight IgE cells in UMAP

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")

IgE <- WhichCells(ascs_int, idents = "IgE")

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/IgE_highlighted_umap.pdf")
DimPlot(ascs_int, reduction = "umap", group.by = "Isotype", cells.highlight = IgE) + 
  scale_color_manual(labels = c("Other", "IgE"), values = c("grey", "red"))
dev.off()

# highlight all isotypes on UMAP
pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Isotypes_highlighted_umap.pdf")
for (i in levels(ascs_int)){
  isotype <- WhichCells(ascs_int, idents = i)
  print(DimPlot(ascs_int, reduction = "umap", group.by = "Isotype", cells.highlight = isotype) + 
    scale_color_manual(labels = c("Other", i), values = c("grey", "red")))
}
dev.off()

### Survival Factors; Remove IgD cells ### 

subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

DefaultAssay(subascs_int) <- "RNA"

surv_facts <- c("TNFRSF13C", # BAFF-R
                "TNFRSF13B", # TACI
                "TNFRSF17" # BCMA
                )

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Survival_Factors.pdf")
DotPlot(subascs_int, features = surv_facts, assay = "RNA", group.by = "Isotype") + RotatedAxis()
VlnPlot(subascs_int, features = surv_facts, assay = "RNA", group.by = "Isotype", sort = TRUE,  pt.size = 0)
FeaturePlot(subascs_int, features = surv_facts)
dev.off()

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Survival_Factors_VlnOnly.pdf", height = 5, width= 8)
VlnPlot(subascs_int, features = surv_facts, assay = "RNA", group.by = "Isotype", sort = TRUE,  pt.size = 0)
dev.off()

### Redo Some Plots: Isotype Plots - no IgD, no points ###

isotype_genes <- c("IGHA1", "IGHA2",
                   "IGHD",
                   "IGHE",
                   "IGHG1", "IGHG2", "IGHG3", "IGHG4",
                   "IGHM"
                   )

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/No_IgD_IGHgenes.pdf")
DotPlot(subascs_int, features = isotype_genes, assay = "RNA", group.by = "Isotype") + RotatedAxis()
VlnPlot(subascs_int, features = isotype_genes, assay = "RNA", group.by = "Isotype", pt.size = 0, sort = TRUE)
dev.off()

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/No_IgD_IGHgenes_heamtap.pdf", width = 20, height = 10)
DoHeatmap(subascs_int, features = isotype_genes, slot = "data")
dev.off()

### IGE DETFs ###

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")

subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

tfs <- c("ATF5", "ARID5B", "KLF2", "ZBTB38", "MEF2C", "KLF3",
         "NR4A1", "JUNB", "XBP1", "EGR1")

DefaultAssay(subascs_int) <- "RNA"

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/IgE_DETFs.pdf")
DotPlot(subascs_int, features = tfs, assay = "RNA", group.by = "Isotype") + RotatedAxis()
VlnPlot(subascs_int, features = tfs, assay = "RNA", group.by = "Isotype", pt.size = 0, sort = TRUE)
FeaturePlot(subascs_int, features = tfs)
dev.off()


### Cluster specific Pops 1-5 markers ###

pops <- c("CD19", "SDC1", "CD38")

DefaultAssay(ascs_int) <- "RNA"

ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)
pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Cluster_specific_PopMarks.pdf")
DotPlot(ascs_int, features = pops, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
VlnPlot(ascs_int, features = pops, assay = "RNA", group.by = "seurat_clusters", pt.size = 0, sort = FALSE)
VlnPlot(ascs_int, features = pops, assay = "RNA", group.by = "seurat_clusters", pt.size = 0, sort = TRUE)
FeaturePlot(ascs_int, features = tfs)
dev.off()

### Replot UMAP with Isotype Labels w/o IgD cells 4-13-2021 ###

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Isotype_Labels_UMAP_noIgD.pdf")
DimPlot(subascs_int, reduction = "umap", label = FALSE)
DimPlot(subascs_int, reduction = "umap", label = TRUE)
dev.off()

# Redo isotype proportions

# create a dataframe
data <- data.frame(subascs_int$Isotype,subascs_int$seurat_clusters,1)

# Stacked + percent
p <- ggplot(data, aes(fill=subascs_int.Isotype, y=X1, x=subascs_int.seurat_clusters)) + 
  ggtitle("Proportion of Isotypes per Cluster") +
  xlab("Cluster Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual("legend", values = c("IgA1" = "gold", "IgA2" = "orange1", 
                                         "IgE" = "palegreen4",
                                         "IgG1" = "lightskyblue", "IgG2" = "steelblue3", 
                                         "IgG3" = "royalblue","IgG4" = "royalblue4", 
                                         "IgM" = "purple3")) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black"))

p <- p + guides(fill=guide_legend(title="Isotypes"))

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
ggsave(filename = "Isotype_Proportions_Integrated_Clusters_NoIgD.jpeg",
       plot = print(p))


### Redo MHC I and II w/o IgD

DefaultAssay(subascs_int) <- "RNA"

ap_genes <- c("HLA-A", "HLA-B", "B2M", "HLA-E", "HLA-F", # MHC I
              "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "AP1S2", "CD74", "CTSS" # MHC II
)

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
pdf(file = "MHC_degs_noIgD.pdf")
DotPlot(subascs_int, features = ap_genes, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(subascs_int, features = ap_genes, assay = "RNA", group.by = "Isotype") + RotatedAxis()
dev.off()

# IFN-gamma + LAPTM5 and SARAF

ifn_misc <- c("HLA-A", "HLA-B", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1",
              "HLA-DRA", "HLA-DRB1", "B2M", "HLA-E", "HLA-F", "ISG20", "PTPN6",
              "LAPTM5", "SARAF")

pdf(file = "IFN-gama_misc_degs_noIgD.pdf")
DotPlot(subascs_int, features = ifn_misc, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(subascs_int, features = ifn_misc, assay = "RNA", group.by = "Isotype") + RotatedAxis()
dev.off()

# Canonical TFs

tfs <- c("IRF4", "PRDM1", "XBP1", "CD19", "CD38", "SDC1", "PAX5")

pdf(file = "Canonical_TFs_noIgD.pdf")
DotPlot(subascs_int, features = tfs, assay = "RNA", group.by = "seurat_clusters") + RotatedAxis()
DotPlot(subascs_int, features = tfs, assay = "RNA", group.by = "Isotype") + RotatedAxis()
FeaturePlot(subascs_int, features = tfs)
dev.off()

#### 4-28-2021 Volcano Plot for All Isotype's DEGs, CD19+ and CD19- ####

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")

subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

library(EnhancedVolcano)
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/Volcano_Plots/")

DefaultAssay(subascs_int) <- "RNA"

for (i in levels(subascs_int)){
  print(i)
  
  degs <- FindMarkers(subascs_int, ident.1 = i)
  filename <- paste("Isotype_", i, "_DEGs_Volcano.pdf", sep = "")
  plot_subtitle <- paste(i, "ASCs", sep = " ")
  
  pdf(file = filename)
  print(
  EnhancedVolcano(degs,
                  lab = rownames(degs),
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  FCcutoff = 0.25,
                  pCutoff = 0.05,
                  subtitle = plot_subtitle,
                  legendPosition = "bottom"
  ))
  dev.off()
}

DefaultAssay(sle_int) <- "RNA"

FeaturePlot(sle_int, features = "CR2")
DotPlot(sle_int, features = c("TNFRSF13C", "TNFRSF17", "TNFRSF11A", "TNFRSF3", "TNFRSF25", "TNFRSF6B"))

### B-cell like Genes 5-06-2021 ###

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

subascs_int <- SetIdent(subascs_int, value = "seurat_clusters")

DefaultAssay(subascs_int) <- "RNA"

genes <- c(
  "ATF6",
  "CD27", "BACH1", "PAX5", "CD40",
  "CXCR4", "SELL", "PRDM1", "XBP1",
  "MS4A1", "BCL2",
  "ZBTB32", "ZBTB20", "ZBTB24",
  "LTB", "ERN1", "STAT3",
  "SETBP1", # TBET targets
  "CD86", "CCRL2", "IL18RAP" # TBET targets
)

genes <- c(
  "CXCR4", "IL4R", "FGR", "TRAF5", "LYN",
  "TNFSR13B", "TNFSF4",
  "HMGB1", "HMGB2", "EZH2"
)

genes <- c(
  "BACH2", "PAX5", "MS4A1", "BCL2", "BCL6", "XBP1"
)

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
pdf("Immaturity_dotplotv2.pdf")
DotPlot(subascs_int, features = genes, 
        group.by = "seurat_clusters") + RotatedAxis()

DotPlot(subascs_int, features = genes, 
        group.by = "Isotype") + RotatedAxis()
dev.off()

DotPlot(subascs_int, features = c("FCER2", "BCL6", "EGR1", "IRF4", "PAX5", "STAT6"), 
        group.by = "Isotype") + RotatedAxis()

# GC genes
gc_genes <- c(
  "XBP1", "BACH2", "HMGB2", "EZH2", "MCM2", 
  "MCM7", "TCF19", "FOXM1", "NR3C1", "RB1", 
  "UHRF1", "E2F1", "SMARCA4", "CD19", "CD38", "SDC1"
)

# Interestingly, several of the transcription factors that were 
# expressed during the earliest stages of B cell development were 
# also found to be highly expressed in germinal center (GC) cells 
# (GC B cells or CD19+ peripheral blood B cells, data not shown), 
# including XBP1, BACH2, HMGB2, EZH2, MCM2, MCM7, TCF19, FOXM1, NR3C1, 
# RB1, UHRF1, E2F1, and SMARCA4. https://www.jimmunol.org/content/179/6/3662.long

gc_genesv2 <- c(
  "BACH2", "HMGB2", "EZH2", "MCM2", 
  "MCM7", "TCF19", "FOXM1", 
  "UHRF1", "E2F1", "SMARCA4", "DNMT1",
  "CD19", "CD38", "SDC1"
)

pdf(file = "gc_genes_C89.pdf")
DotPlot(subascs_int, features = gc_genesv2, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10))
dev.off()

DotPlot(subascs_int, features = gc_genes, 
        group.by = "Isotype") + RotatedAxis()

mitogenes <- c(
  "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4L",
  "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB",
  "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6",
  "MT-ATP8", "MT-RNR2", "ZBTB20"
  )

DotPlot(ascs_int, features = mitogenes, 
        group.by = "seurat_clusters") + RotatedAxis()

VlnPlot(ascs_int, features = ascs_int$percent.mt, group.by = "seurat_clusters")

## Immune Eff. Response

genes <- c("INB6","PKM","TOLLIP","PLCG2","SCIMP","HSPA8","CNN2","PMAIP1","DPP7",
           "WIPF1","NOP53","TBX21","GMFG","POU2AF1","POU2F2","IFI16",
           "APPL1","LSM14A","CTSH","IL4R","CTSS","CTSZ","CYBA","CYBB","PSAP")

DotPlot(ascs_int, features = genes, 
        group.by = "Isotype") + RotatedAxis()

DN2_vs_non_DN2_DEGs_4_14 <- read_csv("~/Dropbox (GaTech)/For_Erin/DARs_vs_DEGs/DN2/DN2_vs_non-DN2_DEGs_4-14.csv")

### 5-12-2021: Check B Cell DEGs (from Erin's ASC vs Non-ASC DEGs) ###

ASC_vs_B_DEGs_updated_3_15 <- read_csv("~/Dropbox (GaTech)/For_Erin/DARs_vs_DEGs/ASCs/ASC_vs_B_DEGs_updated-3-15.csv")

# only keep down regulated degs
B_cell_genes <- ASC_vs_B_DEGs_updated_3_15[ASC_vs_B_DEGs_updated_3_15$avg_log2FC <= 0,]

# plot genes
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

DefaultAssay(subascs_int) <- "RNA"
DotPlot(subascs_int, features = B_cell_genes$X1) + RotatedAxis()

bcell1 <- B_cell_genes[1:200, ]
p1 <- DotPlot(subascs_int, features = bcell1$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell2 <- B_cell_genes[201:400, ]
p2 <- DotPlot(subascs_int, features = bcell2$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell3 <- B_cell_genes[401:600, ]
p3 <- DotPlot(subascs_int, features = bcell3$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell4 <- B_cell_genes[601:800, ]
p4 <- DotPlot(subascs_int, features = bcell4$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell5 <- B_cell_genes[801:1000, ]
p5 <- DotPlot(subascs_int, features = bcell5$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell6 <- B_cell_genes[1001:1200, ]
p6 <- DotPlot(subascs_int, features = bcell6$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell7 <- B_cell_genes[1201:1400, ]
p7 <- DotPlot(subascs_int, features = bcell7$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell8 <- B_cell_genes[1401:1571, ]
p8 <- DotPlot(subascs_int, features = bcell8$X1, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/B_Cell_Genes")
pdf(file = "B_Cell_DEGs_by_Cluster.pdf",
    width = 20, height = 5)
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()

### group by isotype

bcell1 <- B_cell_genes[1:200, ]
p1 <- DotPlot(subascs_int, features = bcell1$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell2 <- B_cell_genes[201:400, ]
p2 <- DotPlot(subascs_int, features = bcell2$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell3 <- B_cell_genes[401:600, ]
p3 <- DotPlot(subascs_int, features = bcell3$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell4 <- B_cell_genes[601:800, ]
p4 <- DotPlot(subascs_int, features = bcell4$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell5 <- B_cell_genes[801:1000, ]
p5 <- DotPlot(subascs_int, features = bcell5$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell6 <- B_cell_genes[1001:1200, ]
p6 <- DotPlot(subascs_int, features = bcell6$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell7 <- B_cell_genes[1201:1400, ]
p7 <- DotPlot(subascs_int, features = bcell7$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

bcell8 <- B_cell_genes[1401:1571, ]
p8 <- DotPlot(subascs_int, features = bcell8$X1, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/B_Cell_Genes")
pdf(file = "B_Cell_DEGs_by_Isotype.pdf",
    width = 20, height = 5)
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()

# 5-13-2021 Subset B cell like genes - IP #

genes <- c("LAPTM5", "SELL", "SP140", "PXK", "STAP1", "MCUB", " HLA−F",
           "HLA−B", "NCR3", "DDAH2", "C6orf48", "GPSM3", "HLA−DRA", "HLA−DRB5",
           "HLA−DRB1", "HLA−DQA1", "HLA−DQB1", "HLA−DOB", "PSMB8−AS1", "PSMB9",
           "HLA−DMB", "HLA−DMA", "HLA−DOA", "HLA−DPA1", "HLA−DPB1", "CRYBG1", 
           "ITGB2", "AFF4", "SECISBP2L", "CLEC2D")

DefaultAssay(ascs_int) <- "RNA"
DotPlot(ascs_int, features = c("CLEC2D","CD83", "CXCR4", "TNFSF8"), 
        group.by = "seurat_clusters")


# Genes for naïve/memory/germinal center B cells/DNs from Rich:

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")

gc <- c("BCL6", "BACH2", "PAX5", "S1PR1", "EBI2", "STAT6", "MYC", 
        "MEF2B", "MEF2C", "EBF1")
mem <- c("HHEX", "TLE3", "SKI", "ZBTB32", "CD80", "PDL2")
nv_dn_mem <- c("BACH2", "PAX5", "IRF8")

DefaultAssay(subascs_int) <- "RNA"

pdf(file = "gc_mem_bcells.pdf")
DotPlot(subascs_int, features = gc, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
DotPlot(subascs_int, features = mem, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
DotPlot(subascs_int, features = nv_dn_mem, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
DotPlot(subascs_int, features = gc, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
DotPlot(subascs_int, features = mem, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
DotPlot(subascs_int, features = nv_dn_mem, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
dev.off()


### CLEC2D

subascs_int <- SetIdent(subascs_int, value = "seurat_clusters")

c5 <- subset(subascs_int, idents = "5")

DefaultAssay(c5) <- "RNA"

DotPlot(c5, features = c("CLEC2D", "CD83", "CXCR4"), group.by = "Isotype")+ 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

DotPlot(subascs_int, features = c("CLEC2D", "CD83", "CXCR4"), group.by = "seurat_clusters")+ 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))

DefaultAssay(subascs_int) <- "integrated"

FeaturePlot(subascs_int, features = c("CXCR4", "CLEC2D"), blend = TRUE)



DefaultAssay(subascs_int) <- "RNA"
VlnPlot(subascs_int, features = c("CLEC2D", "CD83", "CXCR4"), 
        sort = TRUE, group.by = "seurat_clusters")

VlnPlot(subascs_int, features = c("PAX5", "MS4A1"), 
        sort = TRUE, group.by = "seurat_clusters")

DotPlot(subascs_int, features = c("TBX21", "ITGAX"), 
  group.by = "seurat_clusters")

VlnPlot(subascs_int, features = c("TBX21", "ITGAX"), 
        sort = TRUE, group.by = "Isotype")

VlnPlot(c5, features = c("CLEC2D", "CD83", "CXCR4"), 
        sort = TRUE, group.by = "Isotype")

c5_CLEC2D <- WhichCells(c5, expression = CLEC2D > 0)
c5$ids <- colnames(c5)
c5 <- SetIdent(c5, value = c5$ids)
c5_sub <- subset(c5, idents = c5_CLEC2D)

table(c5_sub$Isotype)
# IgA1 IgA2  IgE IgG1 IgG2 IgG3  IgM 
# 16    2    6   25    9    8    2

c5_sub <- SetIdent(c5_sub, value = c5_sub$Isotype)
c5_IgE1 <- subset(c5_sub, idents = "IgE")

c5_CXCR4 <- WhichCells(c5, expression = CXCR4 > 0)
c5_sub2 <- subset(c5, idents = c5_CXCR4)

table(c5_sub2$Isotype)
# IgA1 IgA2  IgE IgG1 IgG2 IgG3 IgG4  IgM 
# 62   15   69  146   62   26    2    6 

c5_sub2 <- SetIdent(c5_sub2, value = c5_sub2$Isotype)

c5_IgE <- subset(c5_sub2, idents = "IgE")

### 5-20-2021

# DotPlot

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))


genes <- c("CD83", "CD86", "CXCR4", "CLEC2D", "POLH")

DefaultAssay(subascs_int) <- "RNA"

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/dotplot_5-20-2021_meeting.pdf")
DotPlot(subascs_int, features = genes, group.by = "seurat_clusters")+ 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10))
dev.off()

# Cell Cycle

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/dotplot_Sphase_meeting.pdf")
DotPlot(subascs_int, features = s.genes, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=6))
dev.off()

pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/dotplot_G2Mphase_meeting.pdf")
DotPlot(subascs_int, features = g2m.genes, group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=5))
dev.off()

subascs_int <- CellCycleScoring(subascs_int, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

VlnPlot(subascs_int, features = subascs_int$S.Score)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=11)

df <- data.frame(subascs_int$seurat_clusters, subascs_int$S.Score, subascs_int$G2M.Score)

p1 <- ggplot(df, aes(x=subascs_int.seurat_clusters, y=subascs_int.S.Score, fill = subascs_int.seurat_clusters)) + 
  geom_violin() +
  scale_fill_manual(values=color_list) +
  NoLegend() +
  labs(x="Cluster", y = "S Phase Score")
  
p2 <- ggplot(df, aes(x=subascs_int.seurat_clusters, y=subascs_int.G2M.Score, fill = subascs_int.seurat_clusters)) + 
  geom_violin() +
  scale_fill_manual(values=color_list) +
  NoLegend() +
  labs(x="Cluster", y = "G2M Phase Score")


pdf(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/cell_cycle_VlnPlots.pdf")
p1
p2
dev.off()


### Check CD19 expression 5-26-2021 ###

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

subascs_int <- SetIdent(subascs_int, value = "seurat_clusters")

DefaultAssay(ascs_int) <- "RNA"

# combine isotype & sort labels

subascs_int <- SetIdent(subascs_int, value = subascs_int$orig.ident)

new.cluster.ids <- c("CD19+", "CD19-", "CD19+", "CD19-", "CD19+", "CD19-")
names(new.cluster.ids) <- levels(subascs_int)
subascs_int <- RenameIdents(subascs_int, new.cluster.ids)

subascs_int$FACS <- Idents(subascs_int)

subascs_int$Iso_Sort <- paste(subascs_int$FACS, subascs_int$Isotype, sep = " ")

DefaultAssay(subascs_int) <- "RNA"

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/")
pdf(file = "cd19_assess.pdf")
DotPlot(subascs_int, features = "CD19", group.by = "Iso_Sort") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10))
dev.off()

### Check CD19 expression 5-26-2021 ###

mhc_genes <- c("HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "AP1S2", "CD74", "CTSS" # MHC II
)

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
pdf(file = "MHC_degsv2.pdf")
DotPlot(subascs_int, features = mhc_genes, group.by = "Isotype") + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10)) +
  xlab("MHC II Genes") + ylab("Isotypes")
dev.off()

### MISC check for Eun 12-07-2021 ###

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

genes <-c("TNFAIP3", "TNFAIP8", "TNFRSF10A", "TNFRSF13B",
          "TNFRSF13C", "TNFRSF14", "TNFRSF17", 'TNFRSF18',
          "TNFRSF4", "TNFSF10")

genes <-c("TNFAIP3", "TNFRSF13B",
          "TNFRSF13C", "TNFRSF14", "TNFRSF17", 'TNFRSF18',
          "TNFRSF4", "TNFSF10")

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")
pdf(file = "TNF-per-Eun-1206-2021.pdf")
DotPlot(subascs_int, features = genes, group.by = "Isotype", assay = "RNA", scale = FALSE) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10)) +
  xlab("TNF Genes") + ylab("Isotypes")

DotPlot(subascs_int, features = genes, group.by = "seurat_clusters", assay = "RNA", scale = FALSE) + 
  RotatedAxis() +
  scale_colour_gradient2(low = "blue", mid = "gray88", high = "firebrick1") +
  theme(axis.text.x = element_text(size=10)) +
  xlab("TNF Genes") + ylab("Isotypes")
dev.off()

# 06-16-2022
# Assess:
# 1.0 Confusion matrix of Ig vs Non-Ig clustering
# 2.0 RandIndex of Ig vs Non-Ig clustering

library(ArchR)
library(pheatmap)
library(Seurat)
library("zoo")
library(ggplot2)
library(readr)
library(dplyr)
library(Dune)
library(RColorBrewer)
library(tidyr)
library(purrr)

data("nuclei", package = "Dune")
theme_set(theme_classic())

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = "Isotype")
subascs_int <- subset(ascs_int, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

ascs_noIG <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Assess_Ig_Clustering/No_IgGenes/No_IgGenes.rds")
ascs_noIG <- SetIdent(ascs_noIG, value = "Isotype")
subascs_noIG <- subset(ascs_noIG, idents = c("IgA1", "IgA2", "IgE", "IgG1", "IgG2", "IgG3", "IgG4", "IgM"))

# Cluster Confusion Matrix

cM <- confusionMatrix(paste0(subascs_int$seurat_clusters), paste0(subascs_noIG$seurat_clusters))
cM

cM <- cM / Matrix::rowSums(cM)
# cluste 5: 85% retention
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots/")

pdf(file = "ClusterConfusionMatrix_pheatmap.pdf")
p
dev.off()

# Rand Index

cMRI <- data.frame(subascs_int$seurat_clusters, 
                    subascs_noIG$seurat_clusters)

colnames(cMRI) <- c("With_Ig", "Without_Ig")

cMRI <- apply(as.matrix.noquote(cMRI),
              2,
              as.numeric)

merger <- Dune(clusMat = cMRI, verbose = TRUE)
# [1] "Without_Ig" "1"          "0"         
# [1] "With_Ig" "1"       "0"      
# [1] "With_Ig" "3"       "0"      
# [1] "With_Ig" "7"       "0"      
# [1] "With_Ig" "0"       "10"     
# [1] "Without_Ig" "6"          "5" 
names(merger)
# "initialMat" "currentMat" "merges"     "ImpMetric"  "metric"


pdf(file = "Ig_RandIndex.pdf")
plotARIs(clusMat = merger$currentMat)
dev.off()






