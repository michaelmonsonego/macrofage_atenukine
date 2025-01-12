
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
myeloid_memory = readRDS("objects/Myeloid_object.rds")
unique(myeloid_memory$Treatment)
macrophage_memory <-  subset(myeloid_memory, idents = c("1" = "1_M2-like", "2" = "2_Non-classical monocytes/\nanti-inflammatory macrophages")) 

macrophage_memory = NormalizeData(macrophage_memory, assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000) 
macrophage_memory = FindVariableFeatures(macrophage_memory, assay = "RNA", selection.method = "vst", nfeatures = 2000) 
macrophage_memory = ScaleData(macrophage_memory, verbose = FALSE)
macrophage_memory = RunPCA(macrophage_memory, npcs = 30, verbose = FALSE)
DimPlot(macrophage_memory, reduction = "pca", group.by="Sample", pt.size=1.1) #  cols = c("#FDCDAC", "#B3E2CD")
ElbowPlot(macrophage_memory, ndims = 30) + theme_classic()
ggsave(file = "figures/macrophage_memory/elbow_plot.png", dpi=300, width=10, height=10)

macrophage_memory <- FindNeighbors(macrophage_memory, reduction = "pca", dims = 1:14) #M# choose ?
macrophage_memory = RunTSNE(macrophage_memory, dims = 1:14)
#M# macrophage_memory = RunUMAP(macrophage_memory, dims = 1:18)

res_seq <- c(.05,.3, .35, .4, .45, .5)
for(res in res_seq){
  macrophage_memory.res_test <- FindClusters(macrophage_memory, resolution = res)
  
  res_tSNE  <- DimPlot(macrophage_memory.res_test, reduction = "tsne",
                       repel = T, label = TRUE, label.size = 5) +
    theme(legend.position = "none") + 
    plot_annotation(title = paste("Res of", res))
  assign(paste0("tSNE_",res), res_tSNE)
}

tSNE_ls <- list(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5)
all_tSNE <- plot_grid(tSNE_0.05,tSNE_0.3,tSNE_0.35,tSNE_0.4,tSNE_0.45,tSNE_0.5, ncol = 3) 
ggsave(all_tSNE,filename = 'figures/macrophage_memory/macrophage_memory_0_14_dim_res_test.png', dpi=300, height=10, width=16)

macrophage_memory <- FindClusters(macrophage_memory, resolution = 0.35)
# saveRDS(macrophage_memory, file = "objects/macrophage_memory_.35.rds") #M# change if adjusted
macrophage_memory <- readRDS("objects/macrophage_memory_.35.rds")

macrophage_memory.markers = FindAllMarkers(macrophage_memory, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "RNA")
Top50Markers =
  macrophage_memory.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(Top50Markers, "Excels/macrophage_memory_DE genes.csv")


#M# DE genes between groups, not clusters
macrophage_memory <- SetIdent(macrophage_memory, value = "Sample") #M# this enables DE genes between each group to the rest, change Idents back to seurat clusters after running:
#M# macrophage_memory <- SetIdent(macrophage_memory, value = "seurat_clusters")
DE_all <- FindAllMarkers(macrophage_memory, only.pos = TRUE)
top_100_per_group <-  DE_all %>% 
  group_by(cluster) %>% 
  top_n(n = 100, wt = avg_log2FC) %>% 
  as.data.frame %>% 
  arrange(cluster,-avg_log2FC)
write_csv(top_100_per_group, "Excels/macrophage_memory_DE_between_groups.csv")

macrophage_memory <- SetIdent(macrophage_memory, value = "seurat_clusters") #M# change Idents back to seurat clusters after running

DimPlot(macrophage_memory, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 10) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.title.x = element_text(size = 30),  # X-axis title size
    axis.title.y = element_text(size = 30),  # Y-axis title size
    axis.text.x = element_text(size = 24),   # X-axis tick labels size
    axis.text.y = element_text(size = 24),   # Y-axis tick labels size
    axis.ticks = element_line(size = 1)      # Adjust axis tick size
  )
ggsave(file = "figures/macrophage_memory/tsne_no_annotation.png", dpi=300, width=10, height=10)


macrophage_memory = RunUMAP(macrophage_memory, dims = 1:14)
DimPlot(macrophage_memory, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 10) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.title.x = element_text(size = 30),  # X-axis title size
    axis.title.y = element_text(size = 30),  # Y-axis title size
    axis.text.x = element_text(size = 24),   # X-axis tick labels size
    axis.text.y = element_text(size = 24),   # Y-axis tick labels size
    axis.ticks = element_line(size = 1)      # Adjust axis tick size
  )
ggsave(file = "figures/macrophage_memory/umap_no_annotation_dims_14.png", dpi=300, width=10, height=10)



#M# checking ----------------
FeaturePlot(macrophage_memory, features = c("Mki67", "CD11b"), order=TRUE,pt.size=0.5, reduction="tsne", ncol=1)
ggsave(file = "figures/macrophage_memory/proliferation.png", dpi=300, width=4, height=4)

FeaturePlot(macrophage_memory, features = c("Itgb2", "Itgam"), order=TRUE,pt.size=0.5, reduction="tsne", ncol=2)



#M# infiltration signiture ------------------
infiltration_list = list(c("Ccr2", "Cx3cr1", "Ccr5", "Ccr7", "Itgb2", "Itgam", "Mmp9", "Mmp12", "Spp1", "Cd11c"))
macrophage_memory = AddModuleScore(object = macrophage_memory, features = infiltration_list, name = "infiltration", assay = "RNA")
g <- FeaturePlot(object = macrophage_memory, features = "infiltration1",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "infiltration", subtitle = "()")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/macrophage_memory/macrophage_memory_infiltration.png", dpi=300, width=5, height=5)


#M# proliferation signiture ---------------------
proliferation_list = list(c("Mki67", "Pcna", "Ccna2", "Ccnb1", "Ccnb2", "Ccne1", "Cdk1", "Aurka", "Aurkb", "Top2a", "Hist1h1", "Hist1h2", "Hist1h3"))
macrophage_memory = AddModuleScore(object = macrophage_memory, features = proliferation_list, name = "proliferation", assay = "RNA")
g <- FeaturePlot(object = macrophage_memory, features = "proliferation1",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#e46467", "#b33336", "#A73033"), order = TRUE)+labs(title = "proliferation", subtitle = "()")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/macrophage_memory/macrophage_memory_proliferation.png", dpi=300, width=5, height=5)


#M# M1 M2 tsne ----------------
M2_list = list(c("Il10", "Tgfb1", "Ccl17","Ccl22", "Cd163", "Mrc1", "Clec10a", "Arg1", "Stat3", "Pparg","C1qb", "C1qa","Cd72","Aif1","Trem2","C3ar1","F13a1","Gpnmb"))
macrophage_memory = AddModuleScore(object = macrophage_memory, features = M2_list, name = "M2", assay = "RNA")
g <- FeaturePlot(object = macrophage_memory, features = "M21",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M2")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/macrophage_memory/macrophage_memory_M2.png", dpi=300, width=5, height=5)

VlnPlot(macrophage_memory, features = c("M21"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Treatment",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#CFBAE1", "#5C7D9D", "#A4DEF9")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE) +labs(title = "M2")
ggsave(file = "figures/macrophage_memory/M2_signature_vln_intermidiate.png", dpi=300, width=6, height=6, limitsize=FALSE)

M1_list = list(c("Tnf","Il1b", "Il6", "Cxcl9", "Cxcl10", "Cd86", "Cd80", "Hla-dra","Nos2", "Ido1", "Stat1", "Irf5","Ifit2","Gbp5","Msrb1","Il18","Hilpda"))
macrophage_memory = AddModuleScore(object = macrophage_memory, features = M1_list, name = "M1", assay = "RNA")
g <- FeaturePlot(object = macrophage_memory, features = "M11",pt.size=0.75, reduction = "tsne", cols=c("grey","grey","#E46467", "#B33336", "#A73033"), order = TRUE)+labs(title = "M1")+ theme(plot.subtitle = element_text(hjust =0.5))
g
ggsave(file = "figures/macrophage_memory/macrophage_memory_M1.png", dpi=300, width=5, height=5)


VlnPlot(macrophage_memory, features = c("M11"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Treatment",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#CFBAE1", "#5C7D9D", "#A4DEF9")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)+labs(title = "M1")
ggsave(file = "figures/macrophage_memory/M1_signature_intermidiate.png", dpi=300, width=6, height=6, limitsize=FALSE)


unique(macrofages$Annotation)
unique(macrofages$Sample)
unique(macrophage_memory$Treatment)
unique(macrophage_memory@active.ident)


#M# boxplot for M1/M2 ------------

df <- macrophage_memory@meta.data %>%
  as.data.frame()

df <- df %>%
  mutate(M1_M2_ratio = M11 / M21)

df <- df %>%
  filter(!is.na(M1_M2_ratio) & M1_M2_ratio != Inf & M1_M2_ratio != 0)

df <- df %>%
  mutate(M1_M2_logRatio = log2(M1_M2_ratio))



ggplot(df, aes(x = Treatment, y = M1_M2_logRatio, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  theme_classic() +
  labs(
    title = "M1/M2 Ratio (log10) by Treatment",
    x = NULL,
    y = "log10(M1/M2 Ratio)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values= c("#CFBAE1", "#5C7D9D", "#A4DEF9")) +
  coord_cartesian(ylim = c(-3, 3))
ggsave(file = "figures/macrophage_memory/M1_M2_signature_ratio_intermidiate.png", dpi=300, width=6, height=6, limitsize=FALSE)

2












#M# vln between groups -------------
VlnPlot(macrophage_memory, features = c("infiltration1"),
        assay = "RNA",
        flip= TRUE,
        group.by = "Sample",
        pt.size = 0
)+ 
  theme_classic() + scale_fill_manual(values= c("#A4DEF9","#5C7D9D", "#CFBAE1")) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 24, face = "italic"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.ticks.y = element_line(size = 0.5),
    strip.text.y = element_text(angle = 0, size = 16, face = "bold")
  ) +
  geom_boxplot(alpha = 0.3, show.legend = FALSE)
ggsave(file = "figures/macrophage_memory/infiltration_signature.png", dpi=300, width=6, height=10, limitsize=FALSE)

























