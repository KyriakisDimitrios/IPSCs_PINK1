
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


dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/2.Differential_Epression/Pairwise")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/2.Differential_Epression/Pairwise")


Combined <- readRDS("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess/IPSCs_Combined.rds")

# =============================== PAIRWISE DF ===============================================
Combined$condition <- as.factor(Combined$condition)
Idents(Combined) <- as.factor(Combined$condition)
cl_combinations <- combn(levels(Combined$condition),2)
cl_combinations <- cl_combinations[,c(5,19,25,30)]
DefaultAssay(Combined) <- "RNA"
Combined <- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))

library(parallel)
pairwise_df <- function (comb,object,cl_combinations){
    DefaultAssay(object) <- "RNA"

    # for(comb in 1:dim(cl_combinations)[2]){
    title <- paste(cl_combinations[,comb],collapse = "_")
    dir.create(title)
    setwd(title)
    target <- "condition"
    idents <- as.vector(cl_combinations[,comb])
    ident.1 <- idents[1]
    print(ident.1)

    ident.2 <- idents[2]

    pbmc.markers <- FindMarkers(object = object,
                                    ident.1 = ident.1,
                                    ident.2 =ident.2,
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "MAST",latent.vars = c("nCount_RNA"))
    pbmc.markers$gene <- rownames(pbmc.markers)
    qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(object@assays$RNA@counts)[1])
    pbmc.markers$qvalue <- qvalue

    top <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]

    to_fc <- top[order(abs(top$avg_logFC),decreasing = TRUE),]
    to_fc_gene <- rownames(to_fc)[1:50]
    #top10 <- top %>% top_n(n = 50, wt = abs(avg_logFC))
    #top10_genes<- rownames(top10)

    temp <- object[,object$condition%in%c(ident.1,ident.2)]
    temp$condition <- as.factor(as.vector(temp$condition))

    # debugonce(annotated_heat)
    pdf("Volcano.pdf")
    plot(EnhancedVolcano(pbmc.markers,
                    lab = pbmc.markers$gene,
                    x = 'avg_logFC',
                    y = 'p_val_adj',subtitle = paste(ident.1,"vs",ident.2,"(FCcutoff=0.6)"),
                    xlim = c(-2, 2),FCcutoff = 0.6))
    dev.off()

    ICSWrapper::annotated_heat(object=temp,
                   row_annotation=c(1),
                   gene_list=to_fc_gene,
                   Rowv=TRUE,
                   gene_list_name="DF_genes",
                   title=title,
                   ordering="condition",One_annot = TRUE)

    DefaultAssay(temp) <- "integrated"
    write.table(pbmc.markers, file = paste0(Sys.Date(),"_TO_EXP_each_",target,"_",title,".tsv"),row.names=FALSE, na="", sep="\t")
    setwd("../")
}

mclapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,object=Combined,cl_combinations=cl_combinations,mc.cores=1)
setwd("../")
# ----------------------------------------------------------------------------------------------
