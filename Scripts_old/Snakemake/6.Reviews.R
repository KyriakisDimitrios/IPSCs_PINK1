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
    paste0(DATADIR,"DADD5_S2_DGE.txt"),
    paste0(DATADIR,"DADA4_S4_DGE.txt"),
    paste0(DATADIR,"DADA3_S3_DGE.txt"),
    paste0(DATADIR,"DADA2_S2_DGE.txt"),
    paste0(DATADIR,"DADA1_S1_DGE.txt"),
    paste0(DATADIR,"DADB4_S4_DGE.txt"),
    paste0(DATADIR,"DADB7_S3_DGE.txt"),
    paste0(DATADIR,"DADC3_S3_DGE.txt")
)
    

condition_names <- c(
    "Control_IPSCs",
    "Control_D06",
    "Control_D10",
    "Control_D15",
    "Control_D21",
    "Control_D26",
    "Control_D35",
    "Control_D50")

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

Combined$condition <- factor(as.factor(Combined$condition), levels = c("Control_IPSCs", "Control_D06"  ,"Control_D10",   "Control_D15",   "Control_D21","Control_D22","Control_D26",
                                                                       "Control_D35"  ,"Control_D50"))


p3 <- DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
p4 <- DimPlot(object = Combined, reduction = "umap", label = TRUE,cols = color_clust)
pdf(paste(Sys.Date(),project,"umap","Seuratw.pdf",sep="_"))
print(p3)
print(p4)
dev.off()

setwd("../")


DefaultAssay(Combined) <- "RNA"


saveRDS(Combined,"All_controls.rds")




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
Controls <- readRDS("./Reviewes/All_controls.rds")


# ==== Fig 2 IPSCs Reviewers
ipsc_markers<-c("MYC",
"POU5F1",
"NANOG",
"TDGF-1",
"L1TD1",
"USP44",
"POLR3G",
"TERF1" ,
"IFITM1",
"DPPA4",
"PRDX1" ,
"LCK",
"DNMT3B",
"IDO1",
"LIN28A"
)
pdf("Figure2rev.pdf",width= 8)
DefaultAssay(Combined) <- "RNA"
Combined<-ScaleData(Combined,ipsc_markers)
DoHeatmap(Combined,ipsc_markers,hjust=0.5,angle=0,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
Controls_DM<-ScaleData(Controls_DM,ipsc_markers)
DoHeatmap(Controls_DM,ipsc_markers,hjust=0.5,angle=0,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()



mDA2 <- c(
"PTCH1"	,
"FZD7"	,
"FABP7",
"CTNNB1",
"SHH"	,
"SOX2"	,
"FOXA2"	,
"OTX2"	,
"PBX1"	,
"HES1"	,
"NTN1",
"SOX6"	,
"DCX"	,
"DDC"	,
"TCF12"	,
"ALCAM",
"SLIT2"	,
"LMO3",
"LMX1A"	,
"WNT5A"	,
"RSPO2"	,
"MSX1"	,
"NEUROG2",
"MSX1"	,
"DMRTA2",	
"PITX2",
"ASCL1",
"PITX3"	,
"NEUROD1",
"CORIN",
"SLIT1",
"TH",
"NR4A2"
)
pdf("Figure3b2.pdf",width= 8)
DefaultAssay(Controls) <- "RNA"
Controls<-ScaleData(Controls,mDA2)
DoHeatmap(Controls,mDA2,angle=0,hjust=0.5,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()



IPSC <- ipsc_markers

Rgl1<- c(
"FABP7",
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"RSPO2",
"MSX1",
"SOX6"
)

Rgl3<- c(
"FABP7",
"SOX2",
"WNT5A",
"CORIN",
"BNC2"
)

ProgM <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"CORIN",
"DCC"
)

ProgFPs <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"DCC"
)



NProg <- c(
"SOX2",
"FOXA2",
"OTX2",
"WNT5A",
"ASCL1",
"NEUROD1",
"NEUROG2",
"TBB3",
"DDC",
"DCX")

NbM <- c(
"FOXA2",
"NEUROD1",
"NEUROG2",
"TBB3",
"DCX")


DA0<- c(
"FOXA2",
"LMX1A",
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH"
)

DA1<- c(
"FOXA2",
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH",
"BNC2",
"SLC18A2",
"SLC6A3"
)

DA2<- c(
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH",
"BNC2",
"SLC18A2",
"SLC6A3",
"LM03",
"ALDH1A1",
"SOX6"
)


split_list <- list(IPSC =IPSC,Rgl=unique(c(Rgl3,Rgl1)),Prog=unique(c(ProgM,ProgFPs)),NProg=unique(c(NProg,NbM)),DA=unique(c(DA0,DA1,DA2)))
unlist_feat <- unlist(split_list)
colnames(Controls@meta.data)
Controls@meta.data <- Controls@meta.data[,c(1:20)]
Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D26","D35","D50"))
Controls <-AddModuleScore(Controls,features = split_list,name = "DA_s")
colnames(Controls@meta.data)[22:26]<-names(split_list)



DotPlot(Controls,features = unique(unlist_feat),group.by = "condition")+coord_flip()
library(Seurat)
my_genes <- unique(c(unlist_feat))
exp <- FetchData(Controls, my_genes)
mat <- matrix(ncol=8,nrow=35)
j=0
for (i in levels(Controls$Timepoints)){
    j=j+1
    min_mat <- exp[Controls$Timepoints==i,]
    mat[,j] <- as.matrix(colMeans(min_mat  > 0))*100
}
rownames(mat) <- rownames(as.matrix(colMeans(min_mat  > 0))*100)
colnames(mat) <- levels(Controls$Timepoints)
print(mat)


ord_mat <- mat[rev(1:35),]
annot_map <- ord_mat[,1:35]
# annot_map[,1]<- c("DA1","DA1","DA1","DA1","DA1","DA1","DA1","DA1",
# "NProg","NProg","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1")
# annot_map <- as.data.frame(annot_map)
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
annotation_row = as.data.frame(annot_map[,1]), display_numbers = round(ord_mat,2))




saveRDS(Combined,"All_controls.rds")

Controls <- readRDS("../Reviews/All_controls.rds")
split_list <- list(IPSC =IPSC,Rgl=unique(c(Rgl3,Rgl1)),Prog=unique(c(ProgM,ProgFPs,NProg,NbM)),DA=unique(c(DA0,DA1,DA2)))
unlist_feat <- unlist(split_list)

colnames(Controls@meta.data)
Controls@meta.data <- Controls@meta.data[,c(1:20)]
Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D26","D35","D50"))

Controls <-AddModuleScore(Controls,features = split_list,name = "DA_s")
colnames(Controls@meta.data)[22:25]<-names(split_list)



DotPlot(Controls,features = unique(unlist_feat),group.by = "condition")+coord_flip()


library(Seurat)

my_genes <- unique(c(unlist_feat))


exp <- FetchData(Controls, my_genes)

mat <- matrix(ncol=8,nrow=27)
j=0
for (i in levels(Controls$Timepoints)){
    j=j+1
    min_mat <- exp[Controls$Timepoints==i,]
    mat[,j] <- as.matrix(colMeans(min_mat  > 0))*100
}
rownames(mat) <- rownames(as.matrix(colMeans(min_mat  > 0))*100)
colnames(mat) <- levels(Controls$Timepoints)
print(mat)



ord_mat <- mat[rev(1:27),]
annot_map <- ord_mat[,1:27]
# annot_map[,1]<- c("DA1","DA1","DA1","DA1","DA1","DA1","DA1","DA1",
# "NProg","NProg","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1","Rgl1")
# annot_map <- as.data.frame(annot_map)
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
annotation_row = as.data.frame(annot_map[,1]), display_numbers = round(ord_mat,2))


pdf("Reviews_DotPlots.pdf")
DotPlot(Controls,features = rev(names(split_list)),group.by = "Timepoints")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
DotPlot(Controls,features = rev(unique(unlist_feat)),group.by = "Timepoints")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
dev.off()

pdf("Reviews_Heatmapss.pdf")
DoHeatmap(Controls,unlist_feat,raster=FALSE,group.by = "Timepoints")
DoHeatmap(Controls,unlist_feat,group.by = "Timepoints")
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
 display_numbers = round(ord_mat,2))#annotation_row = as.data.frame(annot_map[,1]),
dev.off()


pa <- DotPlot(Controls,features = c("Rgl3","Rgl1" ,   "ProgM", "NProg", "NbM",   "DA0",   "DA1"  ,"DA2"  ),group.by="Timepoints")
pa <- DotPlot(Controls,features = c("IPSC","Rgl","Prog","NProg","DA"  ),group.by="Timepoints")

toplot <- pa$data[,c(2,3,4)]
#toplot$features.plot <- factor(toplot$features.plot,levels=c("Rgl3","Rgl1" ,   "ProgM", "NProg", "NbM",   "DA0",   "DA1"  ,"DA2"  ))
toplot$features.plot <- factor(toplot$features.plot,levels=c("IPSC","Rgl","Prog","NProg","DA"  ))

pdf("Skypin_plot.pdf",width=12)
ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    geom_text(aes(label=round(pct.exp,0)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=2.5)+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")

ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")

dev.off()

pdf("Plot_reviews.pdf",width=12)
DotPlot(Controls,features = rev(names(split_list)),group.by = "Timepoints")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
dev.off()

DotPlot(Controls,features = rev(unique(unlist_feat)),group.by = "Timepoints")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    geom_text(aes(label=round(pct.exp,0)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=2.5)+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")
DoHeatmap(Controls,unlist_feat,raster=FALSE,group.by = "Timepoints")
pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
 display_numbers = round(ord_mat,2))#annotation_row = as.data.frame(annot_map[,1]),
dev.off()




scores_meta <-  Controls@meta.data[22:26]
scores_meta_assign <- apply(scores_meta, 1,function(x) {
    if(x[which.max(x)]<0){
        col_n <- "Undefined"
    }else{
        col_n <- colnames(scores_meta)[which.max(x)]
    }
    return(col_n)
})
Controls$Phase_Dif <- factor(scores_meta_assign,
    levels=c("Undefined","IPSC","Rgl","Prog","NProg","DA"))


table_pl <- rbind(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D06"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D06"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D10"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D10"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D15"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D15"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D21"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D21"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D26"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D26"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D35"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D35"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D50"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D50"]))
)
rownames(table_pl) <-  c("IPSCs","D06","D10","D15","D21","D26","D35"  ,"D50")
toplot <- reshape2::melt(table_pl)
colnames(toplot) <- c("TimePoint","Group","value")
pdf("Cell_Percentage.pdf",width=12)
ggplot(data=toplot, aes(x=TimePoint, y=value, fill=Group)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    geom_text(aes(label=round(value*100,0)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=2.5)+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")
dev.off()

















scores_meta[Combined$Timepoints=="IPSCs"&Combined$Phase_Dif=="DA2",]

scores_meta_assign <- apply(scores_meta, 1,function(x) {
    if(x[which.max(x)]<0){
        col_n <- "Negative_Score"
    }else{
        col_n <- colnames(scores_meta)[which.max(x)]
    }
    return(col_n)
})


scores_meta_assign <- apply(scores_meta, 1,function(x) {
    print(length(x))
    col_n <- colnames(scores_meta)[which.max(x)]
    return(col_n)
})


Combined$Phase_Dif <- factor(scores_meta_assign,levels=c("Negative_Score","Rgl3","Rgl1" ,   "ProgM", "NProg", "NbM",   "DA0",   "DA1"  ,"DA2"  ))

DimPlot(Combined,group.by=c("Phase_Dif","Timepoints"), pt.size=2 )

table_pl <- rbind(table(Combined$Phase_Dif[Combined$Timepoints=="IPSCs"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="IPSCs"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D06"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D06"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D10"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D10"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D15"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D15"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D21"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D21"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D26"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D26"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D35"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D35"])),
table(Combined$Phase_Dif[Combined$Timepoints=="D50"]) / sum(table(Combined$Phase_Dif[Combined$Timepoints=="D50"]))
)
rownames(table_pl) <-  c("IPSCs","D06","D10","D15","D21","D26","D35"  ,"D50")
toplot <- reshape2::melt(table_pl)
colnames(toplot) <- c("TimePoint","Group","value")
pdf("Cell_Percentage.pdf",width=12)
ggplot(data=toplot, aes(x=TimePoint, y=value, fill=Group)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    geom_text(aes(label=round(value*100,0)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=2.5)+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")
dev.off()



VlnPlot(Combined,c("TRIM22",
                   "LCK",
                   "IDO1",
                   "HESX1",
                   "DNMT3B",
                   "CXCL5"))

DotPlot(Combined,c("TRIM22",
                   "LCK",
                   "IDO1",
                   "HESX1",
                   "DNMT3B",
                   "CXCL5"))+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

















Rgl1<- c(
"FABP7",
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"RSPO2",
"MSX1",
"SOX6"
)

rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[Rgl1,]
sums_col <- colSums(da2_mat)
Combined$Rgl1 <- sums_col


Rgl3<- c(
"FABP7",
"SOX2",
"WNT5A",
"CORIN",
"BNC2"
)
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[Rgl3,]
sums_col <- colSums(da2_mat)
Combined$Rgl3 <- sums_col


ProgM <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"CORIN",
"DCC"
)
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[ProgM,]
sums_col <- colSums(da2_mat)
Combined$ProgM <- sums_col


ProgFPs <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"DCC"
)
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[ProgFPs,]
sums_col <- colSums(da2_mat)
Combined$ProgFPs <- sums_col



NProg <- c(
"SOX2",
"FOXA2",
"OTX2",
"WNT5A",
"ASCL1",
"NEUROD1",
"DDC",
"DCX")
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[NProg,]
sums_col <- colSums(da2_mat)
Combined$NProg <- sums_col


NbM <- c(
"FOXA2",
"NEUROD1",
"DCX")
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[NbM,]
sums_col <- colSums(da2_mat)
Combined$NbM <- sums_col


DA0<- c(
"FOXA2",
"LMX1A",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"TH"
)
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[DA0,]
sums_col <- colSums(da2_mat)
Combined$DA0 <- sums_col

DA1<- c(
"FOXA2",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"TH",
"BNC2",
"SLC18A2")
rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[DA1,]
sums_col <- colSums(da2_mat)
Combined$DA1 <- sums_col


DA2<- c(
"DDC",
"DCX",
"NR4A2",
"PBX1",
"TH",
"BNC2",
"SLC18A2",
"ALDH1A1",
"SOX6"
)

rna_mat <- Combined@assays$RNA@counts
da2_mat  <- rna_mat[DA2,]
sums_col <- colSums(da2_mat)

Combined$DA2 <- sums_col


table(sums_col[Combined$Timepoints=="D21"] > 0)
sum(table(sums_col[Combined$Timepoints=="D21"] > 0))


DotPlot(Combined,features = rev(names(split_list)),group.by = "Timepoints")+coord_flip()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

pa <- DotPlot(Combined,features = c("Rgl3","Rgl1" ,   "ProgM", "NProg", "NbM",   "DA0",   "DA1"  ,"DA2"  ))
toplot <- pa$data[,c(2,3,4)]

pdf("Skypin_plot.pdf",width=12)
ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    geom_text(aes(label=round(pct.exp,2)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=2.5)+
    theme_minimal()


























# ====================================== PROTEOMICS =========================================
#  ================= Figure 7 
# FIG 7A
proteomics_volc_D25 <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/Proteomics_D25_Volcano.csv")
proteomics_volc_D40 <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/Proteomics_D40_Volcano.csv")





# FIG 7B
proteomics <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/selected.proteins.short.csv")

# PINK1_D25_R1 	PINK1 	D25 	R1 	LFQ.intensity.D25.66.8
# PINK1_D25_R2 	PINK1 	D25 	R2 	LFQ.intensity.D25.66.7
# Ctrl_D25_R1 	Ctrl 	D25 	R1 	LFQ.intensity.D25.C8
# Ctrl_D25_R2 	Ctrl 	D25 	R2 	LFQ.intensity.D25.C7
# PINK1_D40_R1 	PINK1 	D40 	R1 	LFQ.intensity.D40.66.8
# PINK1_D40_R2 	PINK1 	D40 	R2 	LFQ.intensity.D40.66.7
# Ctrl_D40_R1 	Ctrl 	D40 	R1 	LFQ.intensity.D40.C.8
# Ctrl_D40_R2 	Ctrl 	D40 	R2 	LFQ.intensity.D40.C.7


col_keep <- c("Gene.names",
"LFQ.intensity.D25.66.7",
"LFQ.intensity.D25.66.8",       
"LFQ.intensity.D25.C7",      
"LFQ.intensity.D25.C8",        
"LFQ.intensity.D40.66.7",       
"LFQ.intensity.D40.66.8",       
"LFQ.intensity.D40.C.7",       
"LFQ.intensity.D40.C.8")

proteomics_filt <- proteomics[,col_keep]

# > grep("CPE",proteomics_filt$Gene.names)
# integer(0)
# > match("CPE",proteomics_filt$Gene.names)
# [1] NA
# > grep("ALCAM",proteomics_filt$Gene.names)
# integer(0)
proteomics_filt




