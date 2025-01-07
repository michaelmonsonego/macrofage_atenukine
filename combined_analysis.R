

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

macrofages <- readRDS("objects/macrofage_zoom_.35.rds")
macrophage_memory <- readRDS("objects/macrophage_memory_.35.rds")

macrofages$Batch <- "one_shot"

macrophage_memory$Batch <- "memory"


combined <- merge(x = macrofages, 
                    y = macrophage_memory,
                    add.cell.ids = c("one_shot","memory"))  

combined = NormalizeData(combined, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
combined = FindVariableFeatures(combined, assay = "RNA", selection.method = "vst", nfeatures = 2000) 
combined = ScaleData(combined, verbose = FALSE)
combined = RunPCA(combined, npcs = 30, verbose = FALSE)
DimPlot(combined, reduction = "pca", group.by="Sample", pt.size=1.1) #  cols = c("#FDCDAC", "#B3E2CD")
ElbowPlot(combined, ndims = 30) + theme_classic()
ggsave(file = "figures/combined/elbow_plot.png", dpi=300, width=10, height=10)

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:14) #M# choose ?
combined = RunTSNE(combined, dims = 1:14)
#M# combined = RunUMAP(combined, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  combined.res_test <- FindClusters(combined, resolution = res)
  
  res_tSNE  <- DimPlot(combined.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/combined/combined_0_14_dim_res_test.png', dpi=300, height=10, width=16)

combined <- FindClusters(combined, resolution = 0.35)
# saveRDS(combined, file = "objects/combined_macrophage_.35.rds") #M# change if adjusted
combined <- readRDS("objects/combined_macrophage_.35.rds")


combined.markers = FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
Top50Markers =
  combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers, "Excels/combined_macrofage_zoom_DE genes.csv")





































