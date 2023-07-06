#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(pheatmap)
library("rfUtilities")
library(Seurat)
library("zoo")
library(ggplot2)
library(tidyverse) 
library(ReactomeGSA)
library(ReactomeGSA.data)
library(pathfindR)

# NP CLuster Specific Analysis: Reactome

#BiocManager::install("ReactomeGSA.data")

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)
#gsva_result <- analyse_sc_clusters(ascs_int, verbose = TRUE) # PACE

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/Reactome")
gsva_result <- readRDS("gsva_result.RDS")

plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])

pathway_expression <- pathways(gsva_result)

# simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
pathway_expression[1:3,]

# find the maximum differently expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

# get top 50 pathways
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/Reactome/50_pathways")
for (i in 1:50){
  print(i)
  p <- plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[i])
  filename = paste(i, "_pathway.pdf", sep = "")
  pdf(file = filename)
  print(p)
  dev.off()
}

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/Reactome")
pdf(file = "gsva_heatmap.pdf", height = 10, width = 8)
plot_gsva_heatmap(gsva_result, max_pathways = 50, margins = c(4,20))
dev.off()

### 4-18-2021 ; use Immune Genes as Background set ###
DefaultAssay(ascs_int) <- "RNA"

ascs_int
# 20953 features across 6449 samples within 2 assays 

gex <- GetAssayData(ascs_int, slot = "data")
rowSums(gex)
quantile(rowSums(gex))
# > quantile(rowSums(gex))
# 0%          25%          50%          75%         100% 
# 0.000000     5.902509    80.902919   376.483716 33096.111800 

# > mean(rowSums(gex))
# [1] 531.1325
# > nrow(as.data.frame(gex_df[gex_df$`rowSums(gex)` >= 531,]))
# [1] 3693

gex_df <- as.data.frame(rowSums(gex))
gex_df$gene <- rownames(gex_df)
gex_df_new <- as.data.frame(gex_df[gex_df$`rowSums(gex)` >= 5.9,])
# > nrow(gex_df_new)
# [1] 14215

gex_df_new1 <- as.data.frame(gex_df[gex_df$`rowSums(gex)` >= 80.9,])
# > nrow(gex_df_new1)
# [1] 9477

# Will use 80.9 cut off for 9477 genes
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Immune_Background_Pathway_Analysis/")

write.table(gex_df_new1, file = "immune_genes_bgd.txt")

## Pathway analysis

# IDd by ToppFun: Background Gene List (9054 / 9477)
# Use PathfindR to plot results

toppfun_updegs <- read.delim("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Immune_Background_Pathway_Analysis/IgE_vs_CD19+_UpDEGs.txt")

pathways <- toppfun_updegs[toppfun_updegs$Category %like% "Pathway",]
sorted <- pathways[order(-log10(pathways$q.value.Bonferroni)),]

filtered <- as.data.frame(sorted[sorted$q.value.Bonferroni <= 0.05,])

p <- ggplot(filtered, aes(x=-log10(q.value.Bonferroni), y=Name)) +
  geom_point(size = filtered$Hit.Count.in.Query.List/2) +
  xlab("-log10(Bonferroni p-value)")
  
#geom_dotplot(binaxis = 'y', dotsize = (pathways$Hit.Count.in.Query.List)/100)

## Pathway analysis 5-02-2021 ##

### pathway analysis ###

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)

### Cluster Specific ##

for (i in levels(ascs_int)){
  
  print(i)
  setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/All_Clusters/")
  
  filename <- paste(i, "_vs_ALL_DEGs.csv", sep = "")
  ci <- read.csv(filename)
  ci_df <- data.frame(ci$X, ci$avg_log2FC, ci$p_val_adj)
  
  # Upregulation
  setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/All_Clusters/PathFindR/Upregulation/")
  
  colnames(ci_df) <- c("Gene.symbol", "logFC", "adj.P.Val")
  subset_data <- filter(ci_df, adj.P.Val <= 0.05)
  subset_data <- filter(subset_data, logFC >= 0)
  
  outputdir <- paste("C", i, sep = "")
  
  output_df <- run_pathfindR(subset_data,
                             gene_sets = "KEGG",
                             min_gset_size = 10,
                             output_dir = outputdir)
  
  setwd(outputdir)
  write.table(output_df, file = "pathfindR_out.txt")
  
  pdf("enrichment_plot.pdf", width = 7, height = 5)
  print(enrichment_chart(output_df, top_terms = 15))
  dev.off()

  # Downregulation
  setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/Add.Comparisons/All_Clusters/PathFindR/Downregulation/")
  
  colnames(ci_df) <- c("Gene.symbol", "logFC", "adj.P.Val")
  subset_data <- filter(ci_df, adj.P.Val <= 0.05)
  subset_data <- filter(subset_data, logFC <= 0)
  
  outputdir <- paste("C", i, sep = "")
  output_df <- run_pathfindR(subset_data,
                             gene_sets = "KEGG",
                             min_gset_size = 10,
                             output_dir = outputdir)
  
  setwd(outputdir)
  write.table(output_df, file = "pathfindR_out.txt")
  
  pdf("enrichment_plot.pdf", width = 7, height = 5)
  print(enrichment_chart(output_df, top_terms = 15)) 
  dev.off()

}

## Cell Cycle

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ascs_int <- CellCycleScoring(ascs_int, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(ascs_int, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

ascs_int <- RunPCA(ascs_int, features = c(s.genes, g2m.genes))
DimPlot(ascs_int)



# Exploratory
ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")
ascs_int <- SetIdent(ascs_int, value = ascs_int$seurat_clusters)

DefaultAssay(ascs_int) <- "RNA"

DotPlot(ascs_int, features = c("CD19", "CD38", "SDC1", "CST3", "GSTM2", "IGLL5", "CREB3L2", "MCL1", "ZBTB20", "TBX21"), group.by = "seurat_clusters") + RotatedAxis()

DotPlot(ascs_int, features = c("ITGAX", "CD86", "NR4A1", "CR2", "TFRC", "FAS"), group.by = "seurat_clusters")

DotPlot(ascs_int, features = c("BATF", "NR4A1", "BACH1", "PAX5", "IRF4", "IRF8", "FCRL4"), group.by = "seurat_clusters")

DotPlot(ascs_int, features = c("CD27", "CD24", "FCRL5", "CXCR5", "SLAMF7", "MS4A1", "ITGAX", "TBX21"), 
        group.by = "seurat_clusters") + RotatedAxis()

DotPlot(ascs_int, features = c("PRDM1", "ATF4", "TCF3", "MAZ", "RUNX2", "SNAI3", "TFDP1", "BATF", "ETV6", "MLX", "MIXL1", "BHLHE41", "IRF7", "RREB1", "XBP1", "CREB3", "IRF4", "MYBL2", "KLF13"), 
        group.by = "Isotype") + RotatedAxis()

DotPlot(ascs_int, features = c("NUP133"), 
        group.by = "Isotype") + RotatedAxis()

DotPlot(ascs_int, features = c("NUP133"), 
        group.by = "seurat_clusters") + RotatedAxis()

c5 <- subset(ascs_int, idents = "5")

DotPlot(c5, features = c("CD19", "CD38", "SDC1", "CST3", "GSTM2", "IGLL5", "CREB3L2", "MCL1", "ZBTB20", "TBX21"), group.by = "Isotype") + RotatedAxis()

DotPlot(c5, features = c("ITGAX", "CD86", "NR4A1", "CR2", "TFRC", "FAS"), group.by = "Isotype")

DotPlot(c5, features = c("BATF", "NR4A1", "BACH1", "PAX5", "IRF4", "IRF8", "FCRL4"), group.by = "Isotype")

DotPlot(c5, features = c("CD27", "CD24", "FCRL5", "CXCR5", "SLAMF7", "MS4A1", "ITGAX", "TBX21"), 
        group.by = "Isotype") + RotatedAxis()

DotPlot(c5, features = c("PRDM1", "ATF4", "TCF3", "MAZ", "RUNX2", "SNAI3", "TFDP1", "BATF", "ETV6", "MLX", "MIXL1", "BHLHE41", "IRF7", "RREB1", "XBP1", "CREB3", "IRF4", "MYBL2", "KLF13"), 
        group.by = "Isotype") + RotatedAxis()



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

c5 <- CellCycleScoring(c5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(c5, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)



c5 <- SetIdent(c5, value = ascs_int$Isotype)
C5_IgE_vs_C5<- FindMarkers(c5, ident.1 = "IgE")

table(c5$Phase, c5$Isotype)

#      IgA1 IgA2 IgE IgG1 IgG2 IgG3 IgG4 IgM
# G1    46   11  35   80   36   13    2   5
# G2M   18    8  17   46   21    7    1   0
# S     20    8  27   38   19    8    0   2

### 5-06-2021 fgsea analysis ###

library(fgsea)
library(data.table)
library(ggplot2)
library(tidyverse)

### PRACTICE on IgE vs CD19+ ###

# import DEGs
setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/CD19_Only")
ige_degs <- read.csv("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/CD19_Only/IgE_vs_ALL_DEGs.csv")

# filter significant DEGs:
sig_ige_degs <- filter(ige_degs, ige_degs$p_val_adj <= 0.05)
degs_df <- data.frame(sig_ige_degs$X, sig_ige_degs$avg_log2FC)

# Rank
ranks <- deframe(degs_df)
head(ranks, 20)

# Load the pathways into a named list
pathways.hallmark <- gmtPathways("/Users/maggiebrown/GSEA/gene_sets/h.all.v7.4.symbols.gmt")

# run the fgsea algorithm with 1000 permutations:
fgseaRes <- fgsea(pathways=pathways.all, stats=ranks)

# Tidy results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

sig_fgseaResTidy <- filter(fgseaResTidy, fgseaResTidy$padj <= 0.05)

# Plot the normalized enrichment scores. 
# Color the bar indicating whether or not the pathway was significant:
  
ggplot(sig_fgseaResTidy, aes(reorder(pathway, NES), NES)) +
geom_col(aes(fill=padj<0.05)) +
coord_flip() +
labs(x="Pathway", y="Normalized Enrichment Score",
     title="Hallmark pathways NES from GSEA") + 
theme_minimal()

plotEnrichment(pathways.all[["GSE22886_NAIVE_BCELL_VS_BLOOD_PLASMA_CELL_UP"]], ranks)

topPathways <- fgseaRes[head(order(pval), n=20)][order(NES), pathway]
plotGseaTable(pathways.all[topPathways], ranks,
              fgseaRes, gseaParam=0.5)

pathways <- reactomePathways(names(ranks))
fgseaRes <- fgsea(pathways, exampleRanks, maxSize=500)
head(fgseaRes)


sig_fgseaResTidy[sig_fgseaResTidy$pathway=="GSE22886_NAIVE_BCELL_VS_BM_PLASMA_CELL_UP",]

### 5-06-2021 Toppfun pathway plot ###

ige_up <- read_excel("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Topp_Fun_Results/IgE_vs_CD19pos_upDEGsfiltered.xlsx")
ige_degs <- read.csv("~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/DEGs/CD19_Only/IgE_vs_ALL_DEGs.csv")

# Get Average(Average log2FC) for each pathway
avgavg = c()
for (gene_list in ige_up$`Hit in Query List`){
  print(gene_list)
  splt_gene_list <- strsplit(gene_list, split = ",")
  
  add_log2fc = 0
  for (gene in 1:length(splt_gene_list[[1]])){
    gene_name <- splt_gene_list[[1]][[gene]]
    print(gene_name)
    generow <- filter(ige_degs, ige_degs$X == gene_name)
    add_log2fc <- add_log2fc + generow$avg_log2FC
  }
  # calc mean
  logfc_mean <- add_log2fc/length(splt_gene_list[[1]])
  avgavg <- c(avgavg, logfc_mean)
}

ige_up_new <- data.frame(ige_up, avgavg)

ige_top50 <- top_n(ige_up_new, 50, ige_up_new$q.value.Bonferroni)

ige_top50$Name <- factor(ige_top50$Name, levels = ige_top50$Name[order(-log10(ige_top50$q.value.Bonferroni))])

pdf(file = "~/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Topp_Fun_Results/IgE_Toppfun.pdf", height = 10, width = 10)
ggplot(ige_top50, aes(avgavg, Name)) +
  geom_point(aes(size = Hit.Count.in.Query.List, colour = -log10(q.value.Bonferroni))) +
  xlab("Average Log2FC of \n DEGs in Pathway") +
  ylab("Pathways") +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  scale_colour_gradient(low = "gray", high = "firebrick1") +
  labs(size = "# Genes", col = "-log10(adj. P-value)") +
  labs(title = "Pathways Enriched in IgE ASCs")
dev.off()


