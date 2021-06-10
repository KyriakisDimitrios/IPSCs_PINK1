
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
library(viridis)
## ----Venn------------------------------------------------------------------------------------
library(VennDiagram)
library(EnhancedVolcano)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021")
load("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/Michi_Data/2020-04-17_seurat_elbow_TRUE_Mito-FALSE_Ribo-FALSE_SCT-TRUE_criteria_pass-3/IPSCs_PINK.RData")

Combined$Timepoints <- as.vector(Combined$condition)
Combined$Timepoints[grep("IPSC",Combined$condition)] <-"iPSCs"
Combined$Timepoints[grep("D06",Combined$condition)] <-"Day06"
Combined$Timepoints[grep("D10",Combined$condition)] <-"Day10"
Combined$Timepoints[grep("D15",Combined$condition)] <-"Day15"
Combined$Timepoints[grep("D21",Combined$condition)] <-"Day21"
Combined$Timepoints <- factor(Combined$Timepoints,c("iPSCs","Day06","Day10","Day15","Day21"))



DefaultAssay(Combined) <- "RNA"




Conserved_M <- subset(Combined,subset= condition !="Control_D10")
DefaultAssay(Conserved_M) <- "RNA"
Conserved_M <-NormalizeData(Conserved_M)
Conserved_M <-ScaleData(Conserved_M)
Conserved_M$Timepoints <- as.vector(Conserved_M$condition)
Conserved_M$Timepoints[grep("IPSC",Conserved_M$condition)] <-"IPSCs"
Conserved_M$Timepoints[grep("D06",Conserved_M$condition)] <-"Day06"
Conserved_M$Timepoints[grep("D15",Conserved_M$condition)] <-"Day15"
Conserved_M$Timepoints[grep("D21",Conserved_M$condition)] <-"Day21"
Idents(Conserved_M) <- "Treatment"




# ========================================================== Pairwise DEGs (Group A) =======================================================================
# ========================================================== Pairwise DEGs (Group A) =======================================================================
# ========================================================== Pairwise DEGs (Group A) =======================================================================
dirs_pairs <- list.dirs("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise",full.names = TRUE )[-1]
dirs_pairs <- grep('IPSC|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)
df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()
Timepoint <- c("D06","D15","D21","IPSCs")
for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    timp <- Timepoint[iter]
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
    print(file)
    l1 <- read.table(file,header=TRUE)
    
    l1$cluster <- l1$avg_logFC
    l1$cluster[ l1$avg_logFC2<0] <- "PINK"
    l1$cluster[ l1$avg_logFC2>0] <- "Control"
    l1$p.adj_BY <-p.adjust(l1$p_val,method="BY")
    l1$p.adj_bonferroni <-p.adjust(l1$p_val,method="bonferroni")
    l1$p.adj_BH <- p.adjust(l1$p_val,method="BH")
    ctrl_l1 <- l1[grep("Control",l1$cluster),]
    pink_l1 <- l1[grep("PINK",l1$cluster),]
	all_l1 <-  l1
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.01 & abs(ctrl_l1$avg_logFC ) >0.1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.01 & abs(pink_l1$avg_logFC ) >0.1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]]) 
}



# # =============  Intersect All Common Genes 
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
# 13
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)
# 14
# # ============= Intersect 3 Timepoints Common Genes
cntrl_intesect3 <- Reduce(intersect, df_return_nt_cntrl[1:3])
print(cntrl_intesect3)
# 27
pink_intesect3 <- Reduce(intersect, df_return_nt_pink[1:3])
print(pink_intesect3)
# 28

degs_group_a <- unique(c(cntrl_intesect,cntrl_intesect3,pink_intesect,pink_intesect3))
length(degs_group_a )
# [1] 55

saveRDS(degs_group_a,"Group_A_DEGS.rds")
# 64
# 86


# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------


# ====================================================== FIND CONSERVED MARKERS (Group B -C) ================================================================
# ====================================================== FIND CONSERVED MARKERS (Group B -C) ================================================================
# ====================================================== FIND CONSERVED MARKERS (Group B -C) ================================================================
markers <- FindConservedMarkers(Conserved_M,ident.1 = "Control",ident.2 = "PINK",grouping.var = "Timepoints",test.use="MAST",latent.vars="nCount_RNA",logfc.threshold=0.0)

# ====================================== Group B ===============================================
# ====================================== Group B ===============================================
# ================================== 151 Genes Group B =========================================
index_fc <- c(sign(markers$IPSCs_avg_logFC)==sign(markers$Day06_avg_logFC) & sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("IPSCs_avg_logFC","Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
degs_group_b <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(degs_group_b)
saveRDS(degs_group_b,"Group_B_DEGS.rds")
# 151

# ====================================== Group C ===============================================
# ====================================== Group C ===============================================
# ================================= 172 Genes Group C ===============================================================================
index_fc <- c( sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("IPSCs_avg_logFC","Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
degs_group_c <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(degs_group_c)
saveRDS(degs_group_c,"Group_C_DEGS.rds")
#172

# ====================================== Group D ===============================================
# ====================================== Group D ===============================================
# ================================= 285 Genes Group D ===============================================================================
Conserved_D <- subset(Combined,subset= condition %in% c("Control_D06","Control_D15","Control_D21","PINK1_D06","PINK1_D15","PINK1_D21"))
DefaultAssay(Conserved_D) <- "RNA"
Conserved_D <-NormalizeData(Conserved_D)
Conserved_D <-ScaleData(Conserved_D)
Conserved_D$Timepoints <- as.vector(Conserved_D$condition)
Conserved_D$Timepoints[grep("D06",Conserved_D$condition)] <-"Day06"
Conserved_D$Timepoints[grep("D15",Conserved_D$condition)] <-"Day15"
Conserved_D$Timepoints[grep("D21",Conserved_D$condition)] <-"Day21"
Idents(Conserved_D) <- "Treatment"
markers_d <- FindConservedMarkers(Conserved_D,ident.1 = "Control",ident.2 = "PINK",grouping.var = "Timepoints",test.use="MAST",latent.vars="nCount_RNA",logfc.threshold=0.0)

index_fc <- c( sign(markers_d$Day06_avg_logFC)==sign(markers_d$Day15_avg_logFC) & sign(markers_d$Day15_avg_logFC)==sign(markers_d$Day21_avg_logFC))
sub_markers <- markers_d[markers_d$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
degs_group_d <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(degs_group_d)
saveRDS(degs_group_d,"Group_D_DEGS.rds")
# 285


# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------


# ====================================================== ALL DEGS =================================================================
DEGs <- unique(c(degs_group_a,rownames(degs_group_b),rownames(degs_group_c),rownames(degs_group_d)))
length(DEGs)
saveRDS(DEGs,"DEGS.rds")
# 292 