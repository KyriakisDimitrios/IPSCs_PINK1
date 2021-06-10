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

DATADIR<-"/home/DATA/"
list_of_files <-c(
    paste0(DATADIR,"DADA1_S1_DGE.txt"),
    paste0(DATADIR,"DADA2_S2_DGE.txt"),
    paste0(DATADIR,"DADA3_S3_DGE.txt"),
    paste0(DATADIR,"DADA4_S4_DGE.txt"),
    paste0(DATADIR,"DADA5_S1_DGE.txt"),
    paste0(DATADIR,"DADA6_S2_DGE.txt"),
    paste0(DATADIR,"DADA8_S4_DGE.txt"),
    paste0(DATADIR,"DADD5_S2_DGE.txt"),
    paste0(DATADIR,"DADD6_S3_DGE.txt"))

condition_names <- c(
    "Control_D21",
    "Control_D15",
    "Control_D10",
    "Control_D06",
    "PINK1_D21",
    "PINK1_D15",
    "PINK1_D06",
    "Control_IPSCs",
    "PINK1_IPSCs")

organism<- "human"


dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess")


setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/NEW")

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


## ----remapping-------------------------------------------------------------------------------
dir.create("Aligned_Cond_RegPhase")
setwd("Aligned_Cond_RegPhase")
# ================================== ALLIGN CONDITIONS =========================================
DefaultAssay(Combined) <- "RNA"

Combined$condition <- factor(as.factor(Combined$condition), levels = c("Control_IPSCs", "Control_D06"  ,"Control_D10",   "Control_D15",   "Control_D21",
"PINK1_IPSCs","PINK1_D06",     "PINK1_D15",     "PINK1_D21"))

Combined$Treatment <-as.vector(Combined$condition)
Combined$Treatment[grep("Control",Combined$Treatment)] <- "Control"
Combined$Treatment[grep("PINK",Combined$Treatment)] <- "PINK"
pink.list <-SplitObject(Combined,split.by = "Treatment")

for (i in 1:length(pink.list)) {
    pink.list[[i]] <- SCTransform(pink.list[[i]], verbose = FALSE,vars.to.regress=c("G2M.Score","S.Score"))
}
 # doi: https://doi.org/10.1101/576827
int.features <- SelectIntegrationFeatures(object.list = pink.list, nfeatures = 3000)
pink.list <- PrepSCTIntegration(object.list = pink.list, anchor.features = int.features,
                                     verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = pink.list, normalization.method = "SCT",
                                            anchor.features = int.features, verbose = FALSE)
Seurat.combined <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                      verbose = FALSE)
DefaultAssay(object = Seurat.combined) <- "integrated"
#Seurat.combined$condition <- Idents(object = Seurat.combined)

Combined <- Seurat.combined
Combined$condition <- factor(as.factor(Combined$condition), levels = c("Control_IPSCs", "Control_D06"  ,"Control_D10",   "Control_D15",   "Control_D21",
                                                                       "PINK1_IPSCs","PINK1_D06",     "PINK1_D15",     "PINK1_D21"))

setwd("../")

saveRDS(Combined,"test.rds")
Combined <- readRDS("test.rds")
## ----Clustering------------------------------------------------------------------------------
# ================================== Clustering =========================================
dir.create("Clusters")
setwd("Clusters")

set.seed(123)
library(RcppAnnoy)
# Combined <- ReduceDim(Combined,method="umap",project=project)$Combined
# debugonce(ICSWrapper::reduce_dim)
Combined <- ICSWrapper::reduce_dim(Combined,project=project,assay = "SCT")$Combined#,resolution=c(0.1))$Combined


pdf(paste(Sys.Date(),project,"tsne","projection.pdf",sep="_"))
ICSWrapper::plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,color_list = color_list)
ICSWrapper::plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,color_list = color_list)
dev.off()

ICSWrapper::plot_nFeatures(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
ICSWrapper::plot_tot_mRNA(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)


if(tolower(tool)=="seurat" & elbow){
    p3 <- DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
    p4 <- DimPlot(object = Combined, reduction = "umap", label = TRUE,cols = color_clust)
    pdf(paste(Sys.Date(),project,"umap","Seuratw.pdf",sep="_"))
    print(p3)
    print(p4)
    dev.off()
}
setwd("../")

dim(Combined@assays$RNA@data)

S1<-Combined@assays$RNA@counts
S2<-Combined@assays$RNA@data
Combined@meta.data
SS <- NormalizeData(Combined)
S3<-SS@assays$RNA@data[1:10,1:10]

meta <- Combined@meta.data[,c(1,2,3,4,5,6,9,10,11,12,13)]

write.csv(meta, file=gzfile("Metadata.csv.gz"))
write.csv(S1, file=gzfile("Raw_data.csv.gz"))
write.csv(S2, file=gzfile("Normalized_data.csv.gz"))

S3

saveRDS(Combined,paste0("Clustered_",NewDir,".rds"))

pdf(paste(Sys.Date(),project,"_projection_Aligned_Treatment.pdf",sep="_"))
ICSWrapper::plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list)
ICSWrapper::plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list)
ICSWrapper::plot_cells(Combined,target="Phase",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list)
ICSWrapper::plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list = color_list)
ICSWrapper::plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list = color_list)
ICSWrapper::plot_cells(Combined,target="Phase",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list = color_list)
dev.off()
# ---------------------------------------------------------------------------------------


pdf("Combined_QC.pdf")
res<-ICSWrapper::scatter_gene(Combined,features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.rb"),size=0.9)
print(res)
dev.off()






## ----Developmental_Markers-------------------------------------------------------------------
# ================================== Developmental Stages =========================================
dir.create("Developmental_Markers")
setwd("Developmental_Markers")
DefaultAssay(Combined) <- "RNA"
file <- paste0(WORKDIR,"/Gene_Lists/Paper_IPCS_genes.txt")
file <- "C:/Users/dimitrios.kyriakis/Desktop/PhD/Scripts/ipscs_pink1/Paper_IPCS_genes.txt"

genes_state <-read.table(file)

for(category in levels(as.factor(genes_state$V1))){
  category_genes <- toupper(as.vector(genes_state[genes_state$V1==category,2]))
  category_genes_l <- category_genes[category_genes%in%rownames(Combined)]
  Combined <- AddModuleScore(Combined,features = list(category_genes_l),name = category)

  pdf(paste0(category,"_2umap_projection_condition_regPhase.pdf"),width = 8,height = 8)
  res <- ICSWrapper::scatter_gene(Combined,features = category_genes_l,ncol = 2,nrow = 2,size=1.1)
  plot(res)
  dev.off()
  
}
# early_genes <- toupper(as.vector(genes_state[genes_state$V1=="Early",2]))
# mid_genes <- toupper(as.vector(genes_state[genes_state$V1=="Mid",2]))
# late_genes <- toupper(as.vector(genes_state[genes_state$V1=="Late",2]))
# early_genes_l <- early_genes[early_genes%in%rownames(Combined)]
# mid_genes_l <- mid_genes[mid_genes%in%rownames(Combined)]
# late_genes_l <- late_genes[late_genes%in%rownames(Combined)]
# 
# Combined <- AddModuleScore(Combined,features = list(early_genes_l),name = "Early")
# Combined <- AddModuleScore(Combined,features = list(mid_genes_l),name = "Mid")
# Combined <- AddModuleScore(Combined,features = list(late_genes_l),name = "Late")
# 
# pdf("Early_umap_projection_condition_regPhase.pdf",width = 8,height = 8)
# res <- ICSWrapper::scatter_gene(Combined,features = early_genes_l,ncol = 2,nrow = 2,size=1.1)
# ggarrange(plotlist=res,ncol = 2,nrow = 2)
# dev.off()
# pdf("Mid_umap_projection_condition_regPhase.pdf",width = 8,height = 8)
# res <- ICSWrapper::scatter_gene(Combined,features = mid_genes_l,ncol = 2,nrow = 2,size=1.1)
# ggarrange(plotlist=res,ncol = 2,nrow = 2)
# dev.off()
# pdf("Late_umap_projection_condition_regPhase.pdf",width = 8,height = 8)
# res <- ICSWrapper::scatter_gene(Combined,features = late_genes_l,ncol = 2,nrow = 2,size=1.1)
# ggarrange(plotlist=res,ncol = 2,nrow = 2)
# dev.off()

features <- c("iPSC_identity1","Mda_identity_stage11", "Mda_identity_stage21","Mda_identity_stage31","Mda_identity_stage41", "Non.Mda1")    

pdf("2Development_umap_projection_condition_regPhase.pdf",width = 12,height = 8)
res <- ICSWrapper::scatter_gene(Combined,features = features,ncol = 3,nrow = 2,size=1.1)
print(ggarrange(plotlist=res,ncol = 3,nrow = 2))
dev.off()


scan_dim <- function (object, group.by = "Cell_Type", features, assay = "RNA", 
    method = "heat", organism = "human", cellheight = 20, 
    cellwidth = 20, width = 10) 
{
    title <- paste0(Sys.Date(), "_", group.by)
    require(NMF)
    require(ggplot2)
    require(dplyr)
    require(viridis)
    require(Seurat)
    graphics.off()
    scaled_data <- t(as.matrix(object@assays[[assay]]@counts)[features, 
        ])
    df <- as.data.frame(scaled_data)
    if (organism != "human") {
        colnames(df) <- unlist(lapply(tolower(colnames(df)), 
            ICSWrapper::simpleCap))
        features <- unlist(lapply(tolower(colnames(df)), ICSWrapper::simpleCap))
    }
    df[[group.by]] <- object[[group.by]]
    heat_cl <- aggregate(df[, 1:dim(df)[2] - 1], list(df[, dim(df)[2]])[[1]], 
        mean)
    row.names(heat_cl) <- heat_cl[[group.by]]
    heat_cl[[group.by]] <- NULL
    print("The heatmap created with c1 scale")
    aheatmap(heat_cl, color = viridis(1000), scale = "column", 
        distfun = "correlation", cellwidth = cellwidth,Rowv = NA, 
        cellheight = cellheight, border_color = "gray", 
        filename = paste0(title, "_Mean_heat_extra.pdf"),width=10,height=8)
    library(reshape2)
    id <- df[[group.by]]
    df[[group.by]] <- NULL
    df <- cbind(id = id, df)
    melted_df <- melt(df)
    violin_df <- melted_df
    library(ggplot2)
    print("The Violin created with log1p counts")
    pdf(paste0(title, "_Violin_extra_log.pdf"), width = 20)
    plot(ggplot(violin_df, aes(x = variable, y = log1p(value), 
        fill = get(group.by))) + geom_violin(scale = "width", 
        width = 0.7) + facet_grid(get(group.by) ~ ., switch = "y", 
        space = "free") + cowplot::theme_cowplot() + theme(axis.text.x = element_text(angle = 45, 
        hjust = 1), strip.text.y = element_text(angle = 180), 
        strip.placement = "outside", strip.background = element_rect(colour = "white", 
            fill = "white")) + scale_x_discrete(limits = features) + 
        NoLegend() + ylab("") + xlab("") + scale_fill_manual(values = color_list[[group.by]], 
        name = group.by, na.value = "gray"))
    dev.off()
    pdf(paste0(title, "_Jitter_extra_log.pdf"), width = 20)
    plot(ggplot(violin_df, aes(x = variable, y = log1p(value), 
        fill = get(group.by))) + geom_jitter(aes(color = get(group.by))) + 
        scale_fill_manual(values = color_list[[group.by]], name = group.by, 
            na.value = "gray") + facet_grid(get(group.by) ~ 
        ., switch = "y", space = "free") + cowplot::theme_cowplot() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
            strip.text.y = element_text(angle = 180), strip.placement = "outside", 
            strip.background = element_rect(colour = "white", 
                fill = "white")) + scale_x_discrete(limits = features) + 
        NoLegend() + ylab("") + xlab(""))
    dev.off()
}

Combined <- ScaleData(Combined,rownames(Combined))
category_genes <- toupper(as.vector(genes_state[,2]))
category_genes_l <- category_genes[category_genes%in%rownames(Combined)]
ICSWrapper::annotated_heat(Combined,row_annotation = c(1),gene_list = category_genes_l,ordering = "condition",title="Development_Markers")
ICSWrapper::ics_scanpy(Combined,features = category_genes_l,group.by = "condition")#,Rowv = NA,scale="c1")
setwd("../")
# --------------------------------------------------------------------------------------------------


save.image("IPSCs_PINK.RData")

saveRDS(Combined,"IPSCs_Combined.rds")