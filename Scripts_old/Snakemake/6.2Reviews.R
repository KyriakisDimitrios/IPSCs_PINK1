options(future.globals.maxSize= 2122317824)
# py_config()

#library(ICSWrapper)
library(sctransform)
library(Seurat)
library( RColorBrewer)
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
set.seed(123)
tool="seurat"
project ="IPSCs_pink1"
dataset <- project

# DATADIR<- "home/users/dkyriakis/PhD/Projects/Michi_Data/DATA/"
DATADIR<-"C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/Michi_Data/DATA/"

list_of_files <-c(
    paste0(DATADIR,"DADD6_S3_DGE.txt"),
    paste0(DATADIR,"DADA8_S4_DGE.txt"),
    paste0(DATADIR,"DADA6_S2_DGE.txt"),
    paste0(DATADIR,"DADA5_S1_DGE.txt"),
    paste0(DATADIR,"DADD2_S2_DGE.txt"),
    paste0(DATADIR,"DADD4_S1_DGE.txt")
)
    

condition_names <- c(
    "IPSCs",
    "D06",
    "D15",
    "D21",
    "D26",
    "D34")

organism<- "human"


dir.create("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Reviewes")
setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Reviewes")

source("C:/Users/dimitrios.kyriakis/Desktop/PhD/Scripts/ipscs_pink1/Scripts/Functions.R")
source("C:/Users/dimitrios.kyriakis/Desktop/PhD/Scripts/ipscs_pink1/Scripts/Minors.R")
source("C:/Users/dimitrios.kyriakis/Desktop/PhD/Scripts/ipscs_pink1/Scripts/Networks.R")



colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")[c(5,1,2,3,4,9,6,7,8)]
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]

color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)
# color_cells <-primary.colors(15, steps = 3, no.white = TRUE)

setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Reviewes")


# ================================ SETTING UP ======================================== #
# Number of cells to use
imputation = FALSE
remove_mt=FALSE
remove_ribsomal=FALSE
n_cores=1
elbow = TRUE
SCT=TRUE
criteria_pass=3
min.cells <- 10
min.features <- 200
data_10x <- FALSE

## ----read------------------------------------------------------------------------------------
# options(future.globals.maxSize= 2122317824)
# ==============================================================================================
# ================================ Setup the Seurat objects ====================================
# ==============================================================================================
# ======== Perform an integrated analysis ====
#NewDir <- paste0(Sys.Date(),"_",tool,"_elbow_",elbow,"_Mito-",remove_mt,"_Ribo-",remove_ribsomal,"_SCT-",SCT,"_criteria_pass-",criteria_pass)
#dir.create(NewDir)
#setwd(NewDir)
dir.create("QC")
setwd("QC")
# debugonce(create_cds)


Return_fun <- ICSWrapper::create_cds2(list_of_files=list_of_files,
                                      condition_names=condition_names,
                                      min.features =min.features,min.cells=min.cells,
                                      remove_mt=remove_mt,data_10x=data_10x,
                                      elbow = elbow,tool=tool,n_cores=1,SCT=SCT,
                                      criteria_pass = criteria_pass,vars.to.regress=c("nCount_RNA"))

Combined  <- Return_fun$Combined
Data_List <- Return_fun$Data_List
setwd("../")



dir.create("Clusters")
setwd("Clusters")

set.seed(123)
library(RcppAnnoy)
# Combined <- ReduceDim(Combined,method="umap",project=project)$Combined
# debugonce(ICSWrapper::reduce_dim)
Combined <- ICSWrapper::reduce_dim(Combined,project=project,assay = "RNA")$Combined#,resolution=c(0.1))$Combined

Combined$condition <- factor(as.factor(Combined$condition), levels = c(
    "IPSCs",
    "D06",
    "D15",
    "D21",
    "D26",
    "D34"))


p3 <- DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
p4 <- DimPlot(object = Combined, reduction = "umap", label = TRUE,cols = color_clust)
pdf(paste(Sys.Date(),project,"umap","Seurat_PINK.pdf",sep="_"))
print(p3)
print(p4)
dev.off()

setwd("../")


DefaultAssay(Combined) <- "RNA"


Rgl1<- c(
"FABP7",
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"RSPO2",
"MSX1"
)


NProg <- c(
"ASCL1",
"NEUROD1",
"DDC",
"DCX")


DA0<- c(
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"TH"
)

DA1<-c(
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"TH",
"BNC2",
"SLC18A2",
"CALB1") 


split_list <- list(Rgl1=Rgl1,NProg=NProg,DA0=DA0,DA1=DA1)


Combined <-AddModuleScore(Combined,features = split_list,name = "DA_s")
colnames(Combined@meta.data)[21:24]<-names(split_list)

unlist_feat <- unlist(split_list)
Idents(Combined) <- "condition"
Combined <- ScaleData(Combined,unlist_feat)
pdf("Reviews2_PINK.pdf")
DimPlot(Combined,group.by="condition")
FeaturePlot(Combined,features = names(split_list),order=T)
#dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()

Combined$Timepoints <- factor(unlist(lapply(as.vector(Combined$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D22","D26","D35","D50"))



saveRDS(Combined,"All_PINKS.rds")
Controls <- readRDS("Reviews/All_PINKS.rds")
DotPlot(Controls,features = unique(unlist_feat),group.by = "condition")+CoordFlip()


library(Seurat)

my_genes <- unique(c(unlist_feat))


exp <- FetchData(Controls, my_genes)

mat <- matrix(ncol=6,nrow=19)
j=0
for (i in levels(Controls$condition)){
    j=j+1
    min_mat <- exp[Controls$condition==i,]
    mat[,j] <- as.matrix(colMeans(min_mat  > 0))*100
}
rownames(mat) <- rownames(as.matrix(colMeans(min_mat  > 0))*100)
colnames(mat) <- levels(Controls$condition)
print(mat)



ord_mat <- mat[rev(1:18),]
annot_map <- ord_mat[,1:2]
annot_map[,1]<- c("DA1","DA1","DA1","DA1","DA1","DA1","DA1","DA1",
"NProg","NProg","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1")
annot_map <- as.data.frame(annot_map)
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
annotation_row = as.data.frame(annot_map[,1]), display_numbers = round(ord_mat,2))


pdf("Reviews_DotPlots_PINK.pdf")
DotPlot(Controls,features = rev(names(split_list)),group.by = "condition")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
DotPlot(Controls,features = rev(unique(unlist_feat)),group.by = "condition")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
dev.off()

pdf("Reviews_Heatmaps_PINK.pdf")
DoHeatmap(Controls,unlist_feat,raster=FALSE)
DoHeatmap(Controls,unlist_feat)
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
annotation_row = as.data.frame(annot_map[,1]), display_numbers = round(ord_mat,2))
dev.off()
