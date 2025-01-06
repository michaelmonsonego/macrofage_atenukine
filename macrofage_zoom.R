

library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(sleepwalk)
library(gridExtra)
library(future)
library(reshape2)
library(reticulate)
library(ggplot2)
library(plyr)
library(stringr)
library(RColorBrewer)
library(grid)
library(cowplot)
library(ggsignif)
library(DT)
library(Cairo)
library(snakecase)
library(glue)
library(ggrepel)

setwd("D:/Michael/git_check/macrofage_atenukine")

x1=4.6
#T_cells = readRDS("objects/tils_all_.45_integrgate_annotated (2).rds")
myeloid = readRDS("objects/Attenukine_Saline_Isotype.Myeloid.subset.Merged.rds")

DimPlot(myeloid, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 10) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.title.x = element_text(size = 30),  # X-axis title size
    axis.title.y = element_text(size = 30),  # Y-axis title size
    axis.text.x = element_text(size = 24),   # X-axis tick labels size
    axis.text.y = element_text(size = 24),   # Y-axis tick labels size
    axis.ticks = element_line(size = 1)      # Adjust axis tick size
  )

markers_per_cluster <- FindAllMarkers(myeloid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



Top50Markers_T = markers_per_cluster %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers_T, "excels/ateno_DE genes.csv")



# subset for macrofage clusters-------------
macrofages <- subset(myeloid, idents = c("0" = "0_Cxcl16+ Macrophages", "1" = "1_Cxcl10+ Monocytes")) 


macrofages = NormalizeData(macrofages, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
macrofages = FindVariableFeatures(macrofages, assay = "RNA", selection.method = "vst", nfeatures = 2000) 
macrofages = ScaleData(macrofages, verbose = FALSE)
macrofages = RunPCA(macrofages, npcs = 30, verbose = FALSE)
DimPlot(macrofages, reduction = "pca", group.by="Sample", pt.size=1.1) #  cols = c("#FDCDAC", "#B3E2CD")
ElbowPlot(macrofages, ndims = 30) + theme_classic()
ggsave(file = "figures/macrofage_zoom/elbow_plot.png", dpi=300, width=10, height=10)

macrofages <- FindNeighbors(macrofages, reduction = "pca", dims = 1:14) #M# choose ?
macrofages = RunTSNE(macrofages, dims = 1:14)
#M# macrofages = RunUMAP(macrofages, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  macrofages.res_test <- FindClusters(macrofages, resolution = res)
  
  res_tSNE  <- DimPlot(macrofages.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/macrofage_zoom/macrofage_zoom_0_14_dim_res_test.png', dpi=300, height=10, width=16)

macrofages <- FindClusters(macrofages, resolution = 0.35)
# saveRDS(macrofages, file = "objects/macrofage_zoom_.35.rds") #M# change if adjusted
macrofages <- readRDS("objects/macrofage_zoom_.35.rds")

macrofages.markers = FindAllMarkers(macrofages, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
Top50Markers =
  macrofages.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers, "Excels/macrofage_zoom_DE genes.csv")


#M# DE genes between groups, not clusters
macrofages <- SetIdent(macrofages, value = "Sample") #M# this enables DE genes between each group to the rest, change back to seurat clusters after running:
#M# macrofages <- SetIdent(macrofages, value = "seurat_clusters")
DE_all <- FindAllMarkers(macrofages, only.pos = TRUE)
top_100_per_group <-  DE_all %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)

write_csv(top_100_per_group, "Excels/DE_between_groups.csv")


DimPlot(macrofages, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 10) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.title.x = element_text(size = 30),  # X-axis title size
    axis.title.y = element_text(size = 30),  # Y-axis title size
    axis.text.x = element_text(size = 24),   # X-axis tick labels size
    axis.text.y = element_text(size = 24),   # Y-axis tick labels size
    axis.ticks = element_line(size = 1)      # Adjust axis tick size
  )
ggsave(file = "figures/macrofage_zoom/tsne_no_annotation.png", dpi=300, width=10, height=10)


FeaturePlot(macrofages, features = c("Mki67"), order=TRUE,pt.size=0.5, reduction="tsne", ncol=1)
ggsave(file = "figures/macrofage_zoom/proliferation.png", dpi=300, width=4, height=4)






