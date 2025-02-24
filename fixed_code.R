setwd('C:/Users/Zz/Desktop/xupeng')

library(data.table)
library(Seurat)
library(plyr)
library(NMF)
library(scater)
library(tidyverse)
library(future)
library(reshape2)
library(tibble)
library(ggalluvial)
library(harmony)
library(cowplot)
library(ggpubr)
library(ggbeeswarm)
library(forcats)
library(pheatmap)
library(circlize)
library(viridis)
library(future.apply)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(factoextra)
library(gridExtra)
library(doParallel)
library(foreach)
library(GOSemSim)



mycol <- c("#0070B2", "#5CB3DA", "#B8E3EA", "#DA1735", "#F15E4C", "#FF9F99",
           "#A231A1", "#A37CB7", "#F2D7EE", "#B91372", "#E93B8C", "#ECB2C8",
           "#FF7149", "#F7AE24", "#FBDD7E", "#679436", "#8BBE53", "#CDE391",
           "#067D69", "#00A385", "#98D4C6", "#114B5F", "#028090", "#B2DBBF",
           "#A23E48", "#CD6981", "#FBD0C0", "#788585", "#9CAEA9", "#CCDAD1")

col <- c(colorRampPalette(c(magma(323, begin = 0.15)[1]))(10),
         magma(323, begin = 0.18))

col2 <- colorRampPalette(c("#0070b2","#5ec7dd",
                           "#f3f3f1","#f79676","#da1735"))(50)

col3 <- colorRampPalette(c("#f3f3f1","#b8e3ea",
                           "#5ec7dd","#009bc7","#0070b2"))(50)

col_rb <- c("#343391","#0064af","#0090cc",
            "#00b6db","#01b7c2","#53c0a3",
            "#8dcb8a","#bbd967","#fbd324",
            "#f6bd25","#f4a02e","#ed6f32",
            "#ea5c2e","#d5452f","#c02e2f","#8b2a21")


datafilt<-readRDS('data/fixed_main.rds')

datafilt <- subset(datafilt, cluster_standard =='Macrophages')


ndim = 15
neigh = 20
dist = 0.5
res = 0.5

datafilt <- FindNeighbors(datafilt, k.param = neigh,
                          dims = 1:ndim, reduction = "pca")

datafilt <- FindClusters(datafilt, resolution = res, n.iter = 50)


datafilt <- RunUMAP(datafilt, dims = 1:ndim,
                    n.neighbors = neigh, min.dist = dist, 
                    reduction = "pca", reduction.name = "umap")


plot<-DimPlot(datafilt, pt.size = 0.7, cols = mycol, label = T, repel = T,
        raster = FALSE, label.size = 5, reduction = "umap",
        group.by = "sample")

ggsave("momac/umap_batch.png", plot, dpi = 200,
       width = 8, height = 7, limitsize = FALSE)

plot <- DimPlot(datafilt, pt.size = 0.7, cols = mycol, label = T, repel = T,
                raster = FALSE, label.size = 5, reduction = "umap",
                group.by = "seurat_clusters")

ggsave("momac/umap_clusters.png", plot, dpi = 300, width = 10, height = 8)



select <- c("Ctsd","C1qc", "Fn1",'RT1-Da','Spp1', "Mki67")


plot <- FeaturePlot(datafilt, features = select, #cols = col, 
                    reduction = "umap", ncol = 3,
                    pt.size = 0.2, label = F, order = T)

height = round(length(select)/3 + 0.3, 0) * 5

ggsave("momac//umap_marker2.png", plot, dpi = 200,
       width = 17.5, height = height, limitsize = FALSE)


Idents(datafilt)<-'seurat_clusters'

diff <- FindAllMarkers(datafilt,#ident.1 = 'BM',ident.2 = 'NB',
                       group.by = 'seurat_clusters', #features = genes,
                       min.pct = 0.2,#test.use='MAST',
                       logfc.threshold = 0.3,
                       only.pos = T)
diff <- diff[,-1]

diff_flt <- diff[(abs(diff$avg_log2FC) > 0.5 & diff$p_val_adj < 0.05),]

write.csv(diff_flt,'momac/momac_diff.csv')


celltype=data.frame(ClusterID=0:7,
                    celltype='un')
celltype[celltype$ClusterID %in% c( '5'),2]='FN1'  
celltype[celltype$ClusterID %in% c( '3'),2]='RT1-Da'
celltype[celltype$ClusterID %in% c( '0','6'),2]='Ctsd'   
celltype[celltype$ClusterID %in% c( '1','2'),2]='C1qc'
celltype[celltype$ClusterID %in% c( '4' ),2]='Spp1' 
celltype[celltype$ClusterID %in% c( '7' ),2]='Top2a' 

table(celltype)

datafilt@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  datafilt@meta.data[which(datafilt@meta.data$RNA_snn_res.0.4 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(datafilt@meta.data$celltype)



plot <- DimPlot(datafilt, pt.size = 0.7, cols = mycol, label = T, repel = T,
                raster = FALSE, label.size = 5, reduction = "umap",
                group.by = "celltype")

ggsave("momac/umap_celltype.png", plot, dpi = 300, width = 10, height = 8)





