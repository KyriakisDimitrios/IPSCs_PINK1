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
library(future)
# library(clusterProfiler)
# library(enrichplot)
# library(ReactomePA)
# library(org.Hs.eg.db)
list.of.packages <- c("multtest", "EnhancedVolcano")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)


colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

# ================ READ ARGS =================
args = commandArgs(trailingOnly=TRUE)
print(args)

source("workflow/scripts/Functions.R")
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
    print("Arguments Passed")
}
input_file <- args[1]
threads <- as.numeric(args[2])
output_group_B <- args[3]
output_group_C <- args[4]
# ---------------------------------------------




Combined <- readRDS(input_file)
Conserved_M <- subset(Combined,subset= condition !="Control_D10")
DefaultAssay(Conserved_M) <- "RNA"
Conserved_M <-NormalizeData(Conserved_M)
all.genes <- rownames(Conserved_M)
Conserved_M <- ScaleData(Conserved_M, features = all.genes)
Conserved_M$Timepoints <- as.vector(Conserved_M$date)
Conserved_M$Timepoints[grep("IPSC",Conserved_M$date)] <-"IPSCs"
Conserved_M$Timepoints[grep("D06",Conserved_M$date)] <-"Day06"
Conserved_M$Timepoints[grep("D15",Conserved_M$date)] <-"Day15"
Conserved_M$Timepoints[grep("D21",Conserved_M$date)] <-"Day21"
Idents(Conserved_M) <- "Treatment"

plan("multiprocess", workers = threads)
markers <- FindConservedMarkers(Conserved_M,ident.1 = "Control",ident.2 = "PINK1",grouping.var = "Timepoints",test.use="MAST",latent.vars="nCount_RNA",logfc.threshold=0.0)
index_fc <- c(sign(markers$IPSCs_avg_logFC)==sign(markers$Day06_avg_logFC) & sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("IPSCs_avg_logFC","Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
sub_markers2 <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(sub_markers2)
# 151
write.table(sub_markers2,output_group_B)


# # ======================= Enrichment ==================================
# control_m <- rownames(sub_markers2[sub_markers2$avg_FC < 0,] )
# pink_m <- sub_markers2[sub_markers2$avg_FC > 0,])
# control_m_e <- bitr(control_m, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# pink_m_e <- bitr(pink_m, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# n_c <- length(control_m_e$ENTREZID)
# n_p <- length(pink_m_e$ENTREZID)
# c_d <- data.frame(ENTREZID = control_m_e$ENTREZID , Condition ="Control")
# p_d <- data.frame(ENTREZID = pink_m_e$ENTREZID , Condition ="PINK1")
# all_d <- rbind(c_d,p_d)
# mydf <- as.data.frame(all_d)
# mydf$ENTREZID <- as.vector(mydf$ENTREZID)
# mydf$ENTREZID <-as.numeric(mydf$ENTREZID )
# GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
# Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
# Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


# pdf("result/DEG/GO_Enrichment_DEG_Group_B.pdf",width=22,height=11)
# dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dev.off()


index_fc <- c( sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("IPSCs_avg_logFC","Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
sub_markers2 <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(sub_markers2)
#172
write.table(sub_markers2,output_group_C)


# --------------------------------------------------------------------