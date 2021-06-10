
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
## ----Venn------------------------------------------------------------------------------------
library(VennDiagram)
library(EnhancedVolcano)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
olor_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/")
load("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/Michi_Data/2020-04-17_seurat_elbow_TRUE_Mito-FALSE_Ribo-FALSE_SCT-TRUE_criteria_pass-3/IPSCs_PINK.RData")


# =============================== PAIRWISE DF ===============================================
Combined$condition <- as.factor(Combined$condition)
Idents(Combined) <- as.factor(Combined$condition)
cl_combinations <- combn(levels(Combined$condition),2)
cl_combinations <- cl_combinations[,c(5,13,25,30)]
DefaultAssay(Combined) <- "RNA"
Combined <- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))
Combined$Timepoints <- factor(unlist(lapply(as.vector(Combined$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21"))


#Combined_Sub1 <- subset(Combined,subset=condition%in%c("Control_D21","PINK1_D21"))


#Combined_Sub2 <- subset(Combined,subset=condition%in%c("Control_D21","PINK1_D21"))
#Combined_Sub2 <- subset(Combined,subset=condition%in%c("Control_D15","PINK1_D15"))
Combined_Sub2 <- subset(Combined,subset=condition%in%c("Control_D06","PINK1_D06"))
Combined_Sub2 <- subset(Combined,subset=condition%in%c("Control_IPSCs","PINK1_IPSCs"))
Combined_Sub2[["RNA"]]@counts<-as.matrix(Combined_Sub2[["RNA"]]@counts)+1
Combined_Sub2 <- NormalizeData(Combined_Sub2)


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = Combined_Sub2[["RNA"]]@counts,
                              colData = Combined_Sub2@meta.data,
                              design = ~ Treatment)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# dds$condition <- factor(dds$condition, levels = c("untreated","treated"))
dds$condition <- relevel(dds$Treatment, ref = "Control")

dds1 <- DESeq2::estimateSizeFactors(object = dds)
dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
dds1 <- DESeq2::nbinomWaldTest(object = dds1)
res <- DESeq2::results(dds1)

write.table(res,"D00_Deseq.tsv",sep="\t")
res00 <- res

min_vlc <- min(res$log2FoldChange)
max_vlc <- max(res$log2FoldChange)

pdf("D00_Deseq.pdf")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',subtitle = paste("Wilcox"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("Wilcox")
dev.off()







res <- DESeq2::results(
object = dds1,
contrast = c("group", "Group1", "Group2"),
alpha = 0.05,
...
)


pbmc.markers <- FindMarkers(object = Combined_Sub2,ident.1 = "Control_D21",ident.2 = "PINK1_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "DESeq2",pseudocount.use = 0)




b11 <- log(x = mean(x = expm1(x = b1)) +1)
b22 <- log(x = mean(x = expm1(x = b2)) +1)
b11-b22




pbmc.markers0 <- FindMarkers(object = Combined_Sub1,ident.1 = "Control_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,pseudocount.use = 1)
pbmc.markers01 <- FindMarkers(object = Combined_Sub2,ident.1 = "Control_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,pseudocount.use = 1)


pbmc.markers1 <- FindMarkers(object = Combined_Sub1,ident.1 = "Control_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "MAST",pseudocount.use = 1,latent.vars = c("nCount_RNA"))

pbmc.markers12 <- FindMarkers(object = Combined_Sub2,ident.1 = "Control_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "MAST",pseudocount.use = 0,latent.vars = c("nCount_RNA"))



pbmc.markers0[c("DLK1"),]
pbmc.markers1[c("DLK1"),]
pbmc.markers12[c("DLK1"),]
pbmc.markers2[c("DLK1"),]
pbmc.markers3[c("DLK1"),]
b1 <- Combined_Sub2[["RNA"]]@counts[c("DLK1"),Combined_Sub1$condition=="Control_D21"]
b2 <- Combined_Sub2[["RNA"]]@counts[c("DLK1"),Combined_Sub1$condition=="PINK1_D21"]
log2(mean(b1)/mean(b2))
log(mean(b1)/mean(b2))
b1 <- Combined_Sub1[["RNA"]]@data[c("DLK1"),Combined_Sub1$condition=="Control_D21"]
b2 <- Combined_Sub1[["RNA"]]@data[c("DLK1"),Combined_Sub1$condition=="PINK1_D21"]
log2((mean(b1)+1)/(mean(b2)+1))
log(mean(b1)/mean(b2))
