

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

macrofages$Batch <- "early"
macrophage_memory$Batch <- "intermediate"
colnames(macrofages@meta.data) #M# verification


combined <- merge(x = macrofages, 
                    y = macrophage_memory,
                    add.cell.ids = c("early","intermediate"))  

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
# combined <- readRDS("objects/combined_macrophage_.35.rds")

combined.markers = FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
Top50Markers =
  combined.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers, "Excels/combined_macrofage_zoom_DE genes.csv")

#M# infiltration signiture---------------
infiltration_list = list(c("Ccr2", "Cx3cr1", "Ccr5", "Ccr7", "Itgb2", "Itgam", "Mmp9", "Mmp12", "Spp1", "Cd11c"))
combined = AddModuleScore(object = combined, features = infiltration_list, name = "infiltration", assay = "RNA")
g <- FeaturePlot(object = combined, features = "infiltration1",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "infiltration", subtitle = "()")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/infiltration.png", dpi=300, width=5, height=5)


#M# proliferation signiture---------------
proliferation_list = list(c("Mki67", "Pcna", "Ccna2", "Ccnb1", "Ccnb2", "Ccne1", "Cdk1", "Aurka", "Aurkb", "Top2a", "Hist1h1", "Hist1h2", "Hist1h3"))
combined = AddModuleScore(object = combined, features = proliferation_list, name = "proliferation", assay = "RNA")
g <- FeaturePlot(object = combined, features = "proliferation1",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "proliferation", subtitle = "()")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/proliferation.png", dpi=300, width=5, height=5)


#M# M1 M2 tsne------------
M2_list = list(c("Il10", "Tgfb1", "Ccl17","Ccl22", "Cd163", "Mrc1", "Clec10a", "Arg1", "Stat3", "Pparg","C1qb", "C1qa","Cd72","Aif1","Trem2","C3ar1","F13a1","Gpnmb"))
combined = AddModuleScore(object = combined, features = M2_list, name = "M2", assay = "RNA")
g <- FeaturePlot(object = combined, features = "M21",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M2")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/M2_michael_ayelet.png", dpi=300, width=5, height=5)
M1_list = list(c("Tnf","Il1b", "Il6", "Cxcl9", "Cxcl10", "Cd86", "Cd80", "Hla-dra","Nos2", "Ido1", "Stat1", "Irf5","Ifit2","Gbp5","Msrb1","Il18","Hilpda"))
combined = AddModuleScore(object = combined, features = M1_list, name = "M1", assay = "RNA")
g <- FeaturePlot(object = combined, features = "M11",pt.size=0.75, reduction = "umap", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M1")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/M1_michael_ayelet.png", dpi=300, width=5, height=5)


#M# infiltration &  proliferation signiture umap ---------------
infiltration_list = list(c("Ccr2", "Cx3cr1", "Ccr5", "Ccr7", "Itgb2", "Itgam", "Mmp9", "Mmp12", "Spp1", "Cd11c"))
combined = AddModuleScore(object = combined, features = infiltration_list, name = "infiltration", assay = "RNA")
g <- FeaturePlot(object = combined, features = "infiltration1",pt.size=0.75, reduction = "umap", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "infiltration")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/infiltration_umap.png", dpi=300, width=5, height=5)

proliferation_list = list(c("Mki67", "Pcna", "Ccna2", "Ccnb1", "Ccnb2", "Ccne1", "Cdk1", "Aurka", "Aurkb", "Top2a", "Hist1h1", "Hist1h2", "Hist1h3"))
combined = AddModuleScore(object = combined, features = proliferation_list, name = "proliferation", assay = "RNA")
g <- FeaturePlot(object = combined, features = "proliferation1",pt.size=0.75, reduction = "umap", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "proliferation")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/proliferation_umap.png", dpi=300, width=5, height=5)




#M# batch based tsne ---------------
DimPlot(combined, reduction = "tsne", group.by = "Batch", label = TRUE)
ggsave(file = "figures/combined/tsne_combined_by_groups.png", dpi=300, width=6, height=6, limitsize=FALSE)

DimPlot(combined, reduction = "tsne", repel = T, label = TRUE, label.size = 5)
ggsave(file = "figures/combined/tsne_combined_by_clusters.png", dpi=300, width=6, height=6, limitsize=FALSE)


#M# umap ------------

combined = RunUMAP(combined, dims = 1:14)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 10) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.title.x = element_text(size = 30),  # X-axis title size
    axis.title.y = element_text(size = 30),  # Y-axis title size
    axis.text.x = element_text(size = 24),   # X-axis tick labels size
    axis.text.y = element_text(size = 24),   # Y-axis tick labels size
    axis.ticks = element_line(size = 1)      # Adjust axis tick size
  )
ggsave(file = "figures/combined/umap_no_annotation_dims_14.png", dpi=300, width=10, height=10)

DimPlot(combined, reduction = "umap", group.by = "Batch", label = TRUE)
ggsave(file = "figures/combined/umap_combined_by_clusters.png", dpi=300, width=6, height=6, limitsize=FALSE)



#M# M1 M2 umap------------
M2_list = list(c("Il10", "Tgfb1", "Ccl17","Ccl22", "Cd163", "Mrc1", "Clec10a", "Arg1", "Stat3", "Pparg","C1qb", "C1qa","Cd72","Aif1","Trem2","C3ar1","F13a1","Gpnmb"))
combined = AddModuleScore(object = combined, features = M2_list, name = "M2", assay = "RNA")
g <- FeaturePlot(object = combined, features = "M21",pt.size=0.75, reduction = "umap", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M2")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/umap_M2_sig.png", dpi=300, width=5, height=5)

M1_list = list(c("Tnf","Il1b", "Il6", "Cxcl9", "Cxcl10", "Cd86", "Cd80", "Hla-dra","Nos2", "Ido1", "Stat1", "Irf5","Ifit2","Gbp5","Msrb1","Il18","Hilpda"))
combined = AddModuleScore(object = combined, features = M1_list, name = "M1", assay = "RNA")
g <- FeaturePlot(object = combined, features = "M11",pt.size=0.75, reduction = "umap", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M1")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/combined/umap_M1_sig.png", dpi=300, width=5, height=5)




#M# vln between groups -------------
VlnPlot(combined, features = c("infiltration1"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Sample",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#A4DEF9","#5C7D9D", "#CFBAE1", "#BEE3DB", "#A8DCC1", "#93C9AA")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/combined/infiltration_signature.png", dpi=300, width=6, height=6, limitsize=FALSE)


VlnPlot(combined, features = c("proliferation1"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Sample",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#A4DEF9","#5C7D9D", "#CFBAE1", "#C4E7DA", "#8ADBAE", "#5CCBA3")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/combined/proliferation_signature.png", dpi=300, width=6, height=6, limitsize=FALSE)

VlnPlot(combined, features = c("M11"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Sample",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#A4DEF9","#5C7D9D", "#CFBAE1", "#C4E7DA", "#8ADBAE", "#5CCBA3")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/combined/M1_signature.png", dpi=300, width=6, height=6, limitsize=FALSE)


VlnPlot(combined, features = c("M21"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Sample",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#A4DEF9","#5C7D9D", "#CFBAE1", "#C4E7DA", "#8ADBAE", "#5CCBA3")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/combined/M2_signature.png", dpi=300, width=6, height=6, limitsize=FALSE)
























