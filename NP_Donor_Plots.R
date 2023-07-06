#!/usr/bin/env Rscript

library(pheatmap)
library("rfUtilities")
library(Seurat)
library("zoo")
library(ggplot2)

### 4-05-2021 "Tack on" v1 analyses for manuscript: ###

# 1.0 Assess donor contibutions to:
#   - each cluster
#   - each isotype
# 2.0 Characterize each cluster (briefly)
#   - DEGs: cluster vs all; maybe predict a trajectory of sorts?
#   - ID mito high and prolif. pops
# 3.0 Cluster 2 (mostly IgA2 cells) vs Cluster 10 (mostly IgA1 cells)
#   - DEGs, pathways
# 4.0 Clean up figs (color coordinate with Rich):
#   - Dotplots, UMAPs, barplots

ascs_int <- readRDS(file = "/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/Pure_Filtered_ASCs/3-11-2021_Final/data/Labeled_Integrated.rds")

# 1.0 Proportions of donors per cluster/isotype

setwd("/Users/maggiebrown/Dropbox (GaTech)/Nasal_Polyps/For_Rich/All_Nasal_Polyps_Data/3-11-2021_Final/Extra_Plots")

ascs_int <- SetIdent(ascs_int, value = ascs_int$orig.ident)
levels(ascs_int)

donor_ids <- c(
  "AFRS_1_CD19+", "AFRS_1_CD19-",
  "AFRS_2_CD19+", "AFRS_2_CD19-",
  "AFRS_3_CD19+", "AFRS_3_CD19-"
)

names(donor_ids) <- levels(ascs_int)

ascs_int <- RenameIdents(ascs_int, donor_ids)
ascs_int[["Donor_Labels"]] <- Idents(ascs_int)

ascs_int <- SetIdent(ascs_int, value = ascs_int$Donor_Labels)

# create a dataframe
data <- data.frame(ascs_int$Donor_Labels,ascs_int$seurat_clusters,1)

# Stacked + percent
p <- ggplot(data, aes(fill=ascs_int.Donor_Labels, y=X1, x=ascs_int.seurat_clusters)) + 
  ggtitle("Proportion of Donors per Cluster") +
  xlab("Cluster Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual("legend", values = c(  "AFRS_1_CD19+" = "deepskyblue1", 
                                           "AFRS_1_CD19-" = "dodgerblue2",
                                           "AFRS_2_CD19+" = "chartreuse",
                                           "AFRS_2_CD19-" = "chartreuse4",
                                           "AFRS_3_CD19+" = "darkgoldenrod1",
                                           "AFRS_3_CD19-" = "darkgoldenrod3"
                                         ))

p <- p + guides(fill=guide_legend(title="Donors"))

ggsave(filename = "Donor_Proportions_Integrated_Clusters.jpeg",
       plot = print(p))

# Donor per isotype

# create a dataframe
data <- data.frame(ascs_int$Isotype,ascs_int$Donor_Labels,1)

# Stacked + percent
p <- ggplot(data, aes(fill=ascs_int.Donor_Labels, y=X1, x=ascs_int.Isotype)) + 
  ggtitle("Proportion of Donors per Cluster") +
  xlab("Cluster Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual("legend", values = c(  "AFRS_1_CD19+" = "deepskyblue1", 
                                           "AFRS_1_CD19-" = "dodgerblue2",
                                           "AFRS_2_CD19+" = "chartreuse",
                                           "AFRS_2_CD19-" = "chartreuse4",
                                           "AFRS_3_CD19+" = "darkgoldenrod1",
                                           "AFRS_3_CD19-" = "darkgoldenrod3"
  ))

p <- p + guides(fill=guide_legend(title="Donors"))

ggsave(filename = "Donor_Proportions_Isotypes_with_IgD.jpeg",
       plot = print(p))

# repeat, remove IgD
ascs_int <- SetIdent(ascs_int, value = ascs_int$Isotype)
subascs_int <- subset(ascs_int, idents = c("IgG2", "IgA2", "IgA1", "IgG1", "IgG3", "IgE",  "IgM",  "IgG4"))

# create a dataframe
data <- data.frame(subascs_int$Isotype,subascs_int$Donor_Labels,1)

# Stacked + percent
p <- ggplot(data, aes(fill=subascs_int.Donor_Labels, y=X1, x=subascs_int.Isotype)) + 
  ggtitle("Proportion of Donors per Cluster") +
  xlab("Cluster Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual("legend", values = c(  "AFRS_1_CD19+" = "deepskyblue1", 
                                           "AFRS_1_CD19-" = "dodgerblue2",
                                           "AFRS_2_CD19+" = "chartreuse",
                                           "AFRS_2_CD19-" = "chartreuse4",
                                           "AFRS_3_CD19+" = "darkgoldenrod1",
                                           "AFRS_3_CD19-" = "darkgoldenrod3"
  ))

p <- p + guides(fill=guide_legend(title="Donors"))

ggsave(filename = "Donor_Proportions_Isotypes_withoutgD.jpeg",
       plot = print(p))



