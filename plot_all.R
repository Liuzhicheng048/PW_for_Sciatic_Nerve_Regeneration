library(data.table)
library(Seurat)
library(plyr)
library(RColorBrewer)
library(scater)
library(SingleCellExperiment)
library(tidyverse)
library(future)
library(SingleR)
library(reshape2)
library(tibble)
library(ggalluvial)
library(cowplot)
library(ggpubr)
library(forcats)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)


# 设定颜色 ------------------------------

mycol <- c(brewer.pal(5,"Set1"), brewer.pal(8,"Set2"),
           brewer.pal(11,"Set3"), brewer.pal(12,"Paired"),
           brewer.pal(8,"Accent"), brewer.pal(11,"Spectral"),
           brewer.pal(11,"BrBG"), brewer.pal(11,"PiYG"),
           brewer.pal(11,"PuOr"),brewer.pal(11,"RdBu"))

col <- c(colorRampPalette(c(magma(323, begin = 0.15)[1]))(10),
         magma(323, begin = 0.18))

col2 <- colorRampPalette(c("#0070b2","#5ec7dd",
                           "#f3f3f1","#f79676","#da1735"))(50)

col3 <- colorRamp2(c(-2,-1,0,1,2),
                   c("#0070b2","#5ec7dd","#f3f3f1","#f79676","#da1735"))

#读取数据
datafilt <- readRDS("C:/Users/Zz/Desktop/xupeng/data/mac_dpt.rds")

plot<-DimPlot(datafilt, pt.size =1, cols = mycol,
              label=F, repel=T, raster=FALSE,split.by = 'group',
              label.size=5, reduction = "umap", 
              group.by = "celltype"
)
ggsave("result/umap_celltype.pdf", plot,
       width = 15, height = 7, limitsize = FALSE)

Idents(datafilt)<-'celltype'
diff<-FindAllMarkers(datafilt,min.pct = 0.3,only.pos = T)
top_diff<-diff%>%group_by(.,cluster)%>%slice_max(n = 10, order_by = avg_log2FC)
  
list_genes=list(C1qc=c("C1qc","Cd163",'Lyve1','Mrc1'),
                Ctsd=c('Ctsd',"Apoc4",'Fabp5','Dgat2'),
                FN1=c("S100a4","Fn1","Ccl9",'Mmp9'),
                RT1=c("RT1-Da","RT1-Db1","RT1-Bb",'Cd74'),
                Spp1=c("Spp1","Lgals1",'Apod'),
                Top2a=c("Cdkn3", "Top2a", "Mki67")
)


plot <-DotPlot(datafilt, scale = T, col.min = -1,
               group.by = "celltype", col.max = 1, scale.max = 100,features = list_genes) + 
  scale_colour_gradient2(low = "#0070b2", mid = "#f3f3f1", high = "#da1735") + 
  RotatedAxis() + # 来自Seurat
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank()
  )
ggsave("plot/dotplot_clusters.pdf", plot,  width = 15, height = 7)



# A1 sample time散点图与密度图 =================================================

# 提取坐标数据 ------------------------------

plotdata <- Embeddings(datafilt, reduction = 'umap') %>% data.frame()

plotdata$group <- datafilt$group
plotdata$celltype <- datafilt$celltype


# sample time的分群散点图 ------------------------------

plot <- ggplot(plotdata, aes(x = plotdata[,1],y = plotdata[,2])) + 
  geom_point(aes(color = celltype), size = 0.2) + 
  scale_color_manual(values=c(mycol)) + 
  
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  
  theme_bw() + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  facet_wrap(~group, nrow = 1)

# 保存图形

output_name <- paste0("allcell/group_celltype.png")
ggsave(output_name, plot, dpi = 300, width = 11, height = 5)


#密度图 ------------------------------

plot <- ggplot(plotdata, aes(x = plotdata[,1],y = plotdata[,2])) +
  
  stat_density_2d(aes(fill = ..density..), geom = "raster", h = 2, contour = FALSE) + 
  geom_density2d(size = 0.1, colour = "#FDAF9199", alpha = 0.35, bins = 15, h = 2) +
  scale_fill_viridis(option="magma") + 
  
  theme_bw() + 
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_wrap(~group, nrow = 1)

# 保存图形

output_name <- paste0("allcell/group_denseplot.pdf")
ggsave(output_name, plot,  width = 10, height = 5)


# B2 以sample time为主要分组的比例图 ===========================================

# 提取数据并预处理 ------------------------------

plotdata <- data.frame(cluster = datafilt$cluster_standard,
                       group = datafilt$group,
                       #time = datafilt$sample_time,
                       value = 1)

plotdata <- aggregate(plotdata[,3], by = list(plotdata$group,
                                              plotdata$cluster), FUN = sum)

names(plotdata)[1:3] <- c("group", "cluster",  "value")


# 把数量转成比例 ------------------------------

plotdata<-plotdata %>%
    group_by(group,cluster) %>%
    summarise(n = sum(value)) %>%
    mutate(prop = n / sum(n))


plotdata$label <- paste0(round(plotdata$prop * 100, 0), "%")

plotdata<-plotdata[!plotdata$label=="0%",]

# 开始画图 ------------------------------

plot <- ggplot(plotdata, aes(x = cluster, fill = group)) +
  
  geom_bar(stat = "identity", data = subset(plotdata, group == "Con"), aes(y = prop)) +
  geom_text(data = subset(plotdata, group == "Con"), 
            aes(y = prop, label = label), size = 3, hjust = -0.1) +
  
  geom_bar(stat = "identity", data = subset(plotdata, group == "PB"), aes(y = prop * (-1)) ) +
  geom_text(data = subset(plotdata, group == "PB"), 
            aes(y = prop * (-1), label = label), size = 3, hjust = 1.1) +
  
  scale_fill_manual(values=c(mycol)) + 
  theme_bw() + theme(axis.text = element_text(colour = "black"),
                     panel.grid = element_blank()) + 
  coord_flip() + 
  ylab("") + xlab("") + 
  
  ylim(-max(plotdata$prop)-0.035, max(plotdata$prop)+0.035)

# 保存图形

output_name <- paste0("allcell/all_prop.png")
ggsave(output_name, plot, dpi = 300, width = 6, height = 3)


##-----------------------------macrophage------------------------------------------


datafilt<-readRDS('ajustment/mac_dpt.rds')


# A1 sample time散点图与密度图 =================================================

# 提取坐标数据 ------------------------------

plotdata <- Embeddings(datafilt, reduction = 'umap') %>% data.frame()

plotdata$group <- datafilt$group
plotdata$celltype <- datafilt$celltype


# sample time的分群散点图 ------------------------------

plot <- ggplot(plotdata, aes(x = plotdata[,1],y = plotdata[,2])) + 
  geom_point(aes(color = celltype), size = 0.2) + 
  scale_color_manual(values=c(mycol)) + 
  
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  guides(colour = guide_legend(override.aes = list(size = 3))) + 
  
  theme_bw() + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  facet_wrap(~group, nrow = 1)

# 保存图形

output_name <- paste0("momac/group_celltype.png")
ggsave(output_name, plot, dpi = 300, width = 11, height = 5)


#密度图 ------------------------------

plot <- ggplot(plotdata, aes(x = plotdata[,1],y = plotdata[,2])) +
  
  stat_density_2d(aes(fill = ..density..), geom = "raster", h = 2, contour = FALSE) + 
  geom_density2d(size = 0.1, colour = "#FDAF9199", alpha = 0.35, bins = 15, h = 2) +
  scale_fill_viridis(option="magma") + 
  
  theme_bw() + 
  labs(x = 'UMAP1',y= 'UMAP2',title = '') + 
  theme(strip.background = element_rect(colour = NA, fill = 'grey90'),
        strip.text.x = element_text(size = 15), 
        panel.grid = element_blank(),
        strip.placement = 'outside',
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  facet_wrap(~group, nrow = 1)

# 保存图形

output_name <- paste0("momac/group_denseplot.png")
ggsave(output_name, plot, dpi = 300, width = 10, height = 5)


# B2 以sample time为主要分组的比例图 ===========================================

# 提取数据并预处理 ------------------------------

plotdata <- data.frame(cluster = datafilt$celltype,
                       group = datafilt$group,
                       #time = datafilt$sample_time,
                       value = 1)

plotdata <- aggregate(plotdata[,3], by = list(plotdata$group,
                                              plotdata$cluster), FUN = sum)

names(plotdata)[1:3] <- c("group", "cluster",  "value")


# 把数量转成比例 ------------------------------

plotdata<-plotdata %>%
  group_by(group,cluster) %>%
  summarise(n = sum(value)) %>%
  mutate(prop = n / sum(n))


plotdata$label <- paste0(round(plotdata$prop * 100, 0), "%")

plotdata<-plotdata[!plotdata$label=="0%",]

# 开始画图 ------------------------------

plot <- ggplot(plotdata, aes(x = cluster, fill = group)) +
  
  geom_bar(stat = "identity", data = subset(plotdata, group == "Con"), aes(y = prop)) +
  geom_text(data = subset(plotdata, group == "Con"), 
            aes(y = prop, label = label), size = 3, hjust = -0.1) +
  
  geom_bar(stat = "identity", data = subset(plotdata, group == "PB"), aes(y = prop * (-1)) ) +
  geom_text(data = subset(plotdata, group == "PB"), 
            aes(y = prop * (-1), label = label), size = 3, hjust = 1.1) +
  
  scale_fill_manual(values=c(mycol)) + 
  theme_bw() + theme(axis.text = element_text(colour = "black"),
                     panel.grid = element_blank()) + 
  coord_flip() + 
  ylab("") + xlab("") + 
  
  ylim(-max(plotdata$prop)-0.035, max(plotdata$prop)+0.035)

# 保存图形

output_name <- paste0("momac/all_prop.png")
ggsave(output_name, plot, dpi = 300, width = 6, height = 3)



 





