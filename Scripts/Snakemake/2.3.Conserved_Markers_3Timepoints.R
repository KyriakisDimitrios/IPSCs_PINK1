## ----setup, include=FALSE--------------------------------------------------------------------
set.seed(123)
# library(reticulate)
options(future.globals.maxSize= 2122317824)
library(sctransform)
library(Seurat)
library(RColorBrewer)
library(tictoc)
library(crayon)
library(stringr)
library(Routliers)
library(jcolors)
library(cluster)
library(NMF)
library(ggplot2)
library(ggpubr)
library(cowplot)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
olor_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)


dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/2.Differential_Epression/Conserved_Markers_3Timepoints")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/2.Differential_Epression/Conserved_Markers_3Timepoints")


Combined <- readRDS("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess/IPSCs_Combined.rds")


Conserved_M <- subset(Combined,subset= condition %in% c("Control_D06","Control_D15","Control_D21","PINK1_D06","PINK1_D15","PINK1_D21"))
DefaultAssay(Conserved_M) <- "RNA"
Conserved_M <-NormalizeData(Conserved_M)
Conserved_M <-ScaleData(Conserved_M)
Conserved_M$Timepoints <- as.vector(Conserved_M$condition)
Conserved_M$Timepoints[grep("D06",Conserved_M$condition)] <-"Day06"
Conserved_M$Timepoints[grep("D15",Conserved_M$condition)] <-"Day15"
Conserved_M$Timepoints[grep("D21",Conserved_M$condition)] <-"Day21"
Idents(Conserved_M) <- "Treatment"
markers <- FindConservedMarkers(Conserved_M,ident.1 = "Control",ident.2 = "PINK",grouping.var = "Timepoints",test.use="MAST",latent.vars="nCount_RNA",logfc.threshold=0.0)

index_fc <- c( sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
sub_markers2 <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(sub_markers2)
write.table(sub_markers2,"Conserved_3.txt")
