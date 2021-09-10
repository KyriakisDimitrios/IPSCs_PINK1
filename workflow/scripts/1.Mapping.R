args = commandArgs(trailingOnly=TRUE)
print(args)

# ================= LIBRARIES ===================
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
set.seed(123)

project ="IPSCs_pink1"
colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")[c(5,1,2,3,4,9,6,7,8)]
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

source("workflow/scripts/Functions.R")
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
    print("Arguments Passed")
}




input <- args[1:9]
threads <- as.numeric(args[10])
output <- args[11]

print(input)
print(output)


tm_list <- mclapply(input,function(x){
    temp_s <- readRDS(x)
}, mc.cores = threads)



# select features that are repeatedly variable across datasets for integration
int.features <- SelectIntegrationFeatures(object.list = tm_list, nfeatures = 3000)
tm_list <- PrepSCTIntegration(object.list = tm_list, 
                anchor.features = int.features, verbose = FALSE)
map.anchors <- FindIntegrationAnchors(object.list = tm_list, 
    normalization.method = "SCT", anchor.features = int.features, 
    verbose = FALSE)
# this command creates an 'integrated' data assay
Combined <- IntegrateData(anchorset = map.anchors,
    normalization.method = "SCT", verbose = FALSE)
DefaultAssay(object = Combined) <- "integrated"

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Combined <- CellCycleScoring(Combined, s.features = s.genes,g2m.features = g2m.genes)

DefaultAssay(Combined) <- 'RNA'
pink.list <-SplitObject(Combined,split.by = "condition")
pink.list <- mclapply(pink.list,function(temp_ss){
    temp_ss <- SCTransform(temp_ss, verbose = FALSE,vars.to.regress=c("G2M.Score","S.Score"))
}, mc.cores = threads)
# doi: https://doi.org/10.1101/576827
int.features <- SelectIntegrationFeatures(object.list = pink.list, nfeatures = 3000)
pink.list <- PrepSCTIntegration(object.list = pink.list, anchor.features = int.features,
                                     verbose = FALSE)
int.anchors <- FindIntegrationAnchors(object.list = pink.list, normalization.method = "SCT",
                                            anchor.features = int.features, verbose = FALSE)
Seurat <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                                      verbose = FALSE)
DefaultAssay(object = Seurat) <- "integrated"
# all.genes <- rownames(Seurat)
# Seurat <- ScaleData(Seurat, features = all.genes)


print("Projection")
set.seed(123)
library(RcppAnnoy)
Seurat <- RunPCA(object= Seurat, verbose =FALSE)
pdf("result/Mapping/ElbowPlot.pdf")
ElbowPlot(Seurat)
dev.off()
num_dim <- ICSWrapper::calc_num_pc(object = Seurat, 0.95)
Seurat <- RunUMAP(object = Seurat, reduction = "pca", dims = 1:num_dim)

Seurat <- RunTSNE(object = Seurat, reduction = "pca", dims = 1:num_dim)
Seurat <- FindNeighbors(object = Seurat, reduction = "pca", dims = 1:num_dim)
optimal_output <- optimal_clusters(Seurat, k.max = 10, save = TRUE, resolution = FALSE)
Seurat <- optimal_output$object
opt_num <- optimal_output$opt_num
sil_scor <- optimal_output$sil_scor
Seurat$Treatment <- Seurat$condition
Seurat$condition <- Seurat$orig.ident

print("Plotting")
pdf(paste("result/Mapping/",Sys.Date(),"_umap.pdf",sep=""))
DimPlot(object = Seurat, reduction = "umap", group.by = "Treatment")
DimPlot(object = Seurat, reduction = "umap", group.by = "condition")
DimPlot(object = Seurat, reduction = "tsne", group.by = "condition")
DimPlot(object = Seurat, reduction = "umap", group.by = "date")
DimPlot(object = Seurat, reduction = "umap")
dev.off()


DefaultAssay(object = Seurat) <- "RNA"
all.genes <- rownames(Seurat)
Seurat <- ScaleData(Seurat, features = all.genes)

print("SAVE RDS")
saveRDS(Seurat,output)


# Seurat$sample <- factor(as.factor(Seurat$sample), 
# levels = c("Control_IPSCs", "Control_D06"  ,"Control_D10",   "Control_D15",   "Control_D21",
# "PINK1_IPSCs","PINK1_D06",   "PINK1_D15",   "PINK1_D21"))








# pdf(paste("result/Mapping/",Sys.Date(),"_tsne_projection.pdf",sep=""))
# ICSWrapper::plot_cells(Seurat,target="sample",leg_pos="right",save=FALSE,ncol=1,color_list = color_list)
# ICSWrapper::plot_cells(Seurat,target="Cluster",leg_pos="right",save=FALSE,ncol=1,color_list = color_list)
# dev.off()
# ICSWrapper::plot_nFeatures(Seurat,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
# ICSWrapper::plot_tot_mRNA(Seurat,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
# p3 <- DimPlot(object = Seurat, reduction = "umap", group.by = "sample",cols = color_cond)
# p4 <- DimPlot(object = Seurat, reduction = "umap", label = TRUE,cols = color_clust)
# pdf(paste("result/Mapping/",Sys.Date(),project,"umap","Seuratw.pdf",sep="_"))
# print(p3)
# print(p4)
# dev.off()





