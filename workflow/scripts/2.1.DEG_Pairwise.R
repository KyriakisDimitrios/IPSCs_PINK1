dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path

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
list.of.packages <- c("VennDiagram", "EnhancedVolcano")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
library(VennDiagram)
library(EnhancedVolcano)

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
output_group_A <- args[3]
# ---------------------------------------------


Combined <- readRDS(input_file)

# =============================== PAIRWISE DF ===============================================
Combined$condition <- as.factor(Combined$condition)
Idents(Combined) <- as.factor(Combined$condition)
cl_combinations <- combn(levels(Combined$condition),2)

print(cl_combinations)

cl_combinations <- cl_combinations[,c(5,19,25,30)]


DefaultAssay(Combined) <- "RNA"
Combined <- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))
Combined$Timepoints <- factor(unlist(lapply(as.vector(Combined$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21"))

setwd("result/Pairwise/")
library(parallel)
pairwise_df <- function (comb,object,cl_combinations){
    DefaultAssay(object) <- "RNA"

    # for(comb in 1:dim(cl_combinations)[2]){
    title <- paste(cl_combinations[,comb],collapse = "_")
    dir.create(paste0(title))
    setwd(paste0(title))
    target <- "condition"
    idents <- as.vector(cl_combinations[,comb])
    ident.1 <- idents[1]
    print(ident.1)

    ident.2 <- idents[2]
    print(ident.1)
    print(ident.2)
    pbmc.markers <- FindMarkers(object = object,
                                    ident.1 = ident.1,
                                    ident.2 =ident.2,
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "MAST",latent.vars = c("nCount_RNA"))
    pbmc.markers$gene <- rownames(pbmc.markers)

    top <- pbmc.markers[pbmc.markers$p_val_adj<0.05,]

    to_fc <- top[order(abs(top$avg_logFC),decreasing = TRUE),]
    to_fc_gene <- rownames(to_fc)[1:50]
    #top10 <- top %>% top_n(n = 50, wt = abs(avg_logFC))
    #top10_genes<- rownames(top10)

    temp <- object[,object$condition%in%c(ident.1,ident.2)]
    temp$condition <- as.factor(as.vector(temp$condition))
    pbmc.markers$avg_logFC2 <-  pbmc.markers$avg_logFC/log(2)
    
    min_vlc <- min(pbmc.markers$avg_logFC2)

    max_vlc <- max(pbmc.markers$avg_logFC2)

    # debugonce(annotated_heat)
    pdf(paste0(title,"_Volcano.pdf"))
    plot(EnhancedVolcano(pbmc.markers,
                    lab = pbmc.markers$gene,
                    x = 'avg_logFC2',
                    y = 'p_val_adj',subtitle = paste(ident.1,"vs",ident.2,"(FCcutoff=0.6)"),
                    xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 0.6))
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
Idents(Combined) <- Combined$condition
mclapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,object=Combined,cl_combinations=cl_combinations,mc.cores=threads)
print(getwd())
# ----------------------------------------------------------------------------------------------

dirs_pairs <- list.dirs(full.names = TRUE )[-1]
dirs_pairs <- grep('D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)


df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()

for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
    print(file)
    l1 <- read.table(file,header=TRUE)
    l1$cluster <- l1$avg_logFC
    l1$cluster[ l1$avg_logFC<0] <- "PINK"
    l1$cluster[ l1$avg_logFC>0] <- "Control"
    ctrl_l1 <- l1[grep("Control",l1$cluster),]
    pink_l1 <- l1[grep("PINK",l1$cluster),]
	all_l1 <-  l1
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.01 & abs(ctrl_l1$avg_logFC) >0.4,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.01 & abs(pink_l1$avg_logFC) >0.4,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]])
}



# # ============= Intersect Common Genes
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)

# df_return_nt_cntrl[[1]]
# df_return_nt_pink[[1]]

all_intesect <- Reduce(intersect, df_return_nt_all)
print(all_intesect)

# pdf("All_venn_diagramm.pdf")
# myCol <- brewer.pal(3, "Pastel2")
# draw.triple.venn(area1 = 117, area2 = 110, area3 = 114, n12 = 23, n23 = 18, n13 = 10,
#     n123 = 7, category = c("Day06", "Day15", "Day21"),
#     fill = myCol)
# dev.off()



# =========== Venn Diagrams DF genes
myCol <- brewer.pal(3, "Pastel2")
pdf("Control_venn_diagramm.pdf")
day06 <- df_return_nt_cntrl[[1]]
day15 <- df_return_nt_cntrl[[2]]
day21 <- df_return_nt_cntrl[[3]]

# Generate plot
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
# have a look at the default plot
grid.newpage()
grid.draw(v)
# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

v[[11]]$label <- paste(intersect(intersect(day06, day15),day21), collapse="\n")  
# plot  
grid.newpage()
grid.draw(v)
dev.off()


pdf("Pink1_venn_diagramm.pdf")
day06 <- df_return_nt_pink[[1]]
day15 <- df_return_nt_pink[[2]]
day21 <- df_return_nt_pink[[3]]

# Generate plot
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
# have a look at the default plot
grid.newpage()
grid.draw(v)
# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

v[[11]]$label <- paste(intersect(intersect(day06, day15),day21), collapse="\n")  
# plot  
grid.newpage()
grid.draw(v)
dev.off()





























# # ====================================== COMMONS ===================================
# dirs_pairs <- list.dirs(full.names = TRUE )[-1]
# dirs_pairs <- grep('IPSCs|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)

# df_return_nt_cntrl <- list()
# df_return_nt_pink <- list()
# df_return_nt_all <- list()




# Timepoint <- c("D06","D15","D21","IPSCs")
# for (iter in 1:length(dirs_pairs)){
#     dirs_iter <- dirs_pairs[iter]
#     timp <- Timepoint[iter]
#     file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
#     print(file)
#     l1 <- read.table(file,header=TRUE)
    
#     l1$cluster <- l1$avg_logFC
#     l1$cluster[ l1$avg_logFC2<0] <- "PINK"
#     l1$cluster[ l1$avg_logFC2>0] <- "Control"
#     l1$p.adj_BY <-p.adjust(l1$p_val,method="BY")
#     l1$p.adj_bonferroni <-p.adjust(l1$p_val,method="bonferroni")
#     l1$p.adj_BH <- p.adjust(l1$p_val,method="BH")
#     ctrl_l1 <- l1[grep("Control",l1$cluster),]
#     pink_l1 <- l1[grep("PINK",l1$cluster),]
# 	all_l1 <-  l1
#     df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.05 & abs(ctrl_l1$avg_logFC2) >1,"gene"])
#     df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.05 & abs(pink_l1$avg_logFC2) >1,"gene"])
#     print(length(df_return_nt_cntrl[[iter]]))
#     print(length(df_return_nt_pink[[iter]]))
#     df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]]) 
    
#     pdf(paste0(timp,"_Volcanos.pdf"))
#     l1$avg_logFC2 <- -l1$avg_logFC2
#     l1$avg_logFC <-  -l1$avg_logFC   
#     min_vlc <- min(l1$avg_logFC2)
#     max_vlc <- max(l1$avg_logFC2)
#     plot(EnhancedVolcano(l1,
#                         lab = l1$gene,
#                         x = 'avg_logFC2',
#                         y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
#                         xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+ggtitle(timp))
#     min_vlc <- min(l1$avg_logFC)
#     max_vlc <- max(l1$avg_logFC)
#     plot(EnhancedVolcano(l1,
#                         lab = l1$gene,
#                         x = 'avg_logFC',
#                         y = 'p_val_adj',subtitle = "FC (FCcutoff=1)",
#                         xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change"))
#     min_vlc <- min(l1$avg_logFC2)
#     max_vlc <- max(l1$avg_logFC2)
#     plot(EnhancedVolcano(l1,
#                         lab = l1$gene,
#                         x = 'avg_logFC2',
#                         y = 'p.adj_BY',subtitle = "BY FC2 (FCcutoff=1)",
#                         xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change")+ggtitle(timp))

#     dev.off()

# }











# library(UpSetR )
# glmer_res <- readRDS("../glmer_pval_mat_10percent.rds")
# glmer_res$max_FC2 <- -glmer_res$max_FC/log(2)
# glmer_res$FC2 <- -glmer_res$FC/log(2)
# glmer_res$p_val_adj <- p.adjust(glmer_res$pvalues_l,method="bonferroni",n=dim(Combined@assays$RNA@counts)[1]) 

# pdf("GLMER_Volcanos.pdf")
# min_vlc <- min(glmer_res$max_FC2)
# max_vlc <- max(glmer_res$max_FC2)
# plot(EnhancedVolcano(glmer_res,
#                         lab = rownames(glmer_res),
#                         x = 'max_FC2',
#                         y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
#                         xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Log2 fold change")+ggtitle("GLMM"))

# dev.off()


# glmm_genes  <- rownames(glmer_res[glmer_res$p_val_adj<0.05 & abs(glmer_res$max_FC2)> 1,])
# df_return_nt_all[[5]] <- glmm_genes
# names(df_return_nt_all) <- c("D06","D15","D21","IPSCs","GLMM")
# pdf("Common.pdf",width=12)
# UpSetR::upset(fromList(df_return_nt_all),text.scale =1.9)
# dev.off()

# unlist(df_return_nt_cntrl)





# # =========================================== ENRICHMENT ANALYSIS ==================================
# r_d00 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_IPSCs_PINK1_IPSCs/2021-04-01_TO_EXP_each_condition_Control_IPSCs_PINK1_IPSCs.tsv",header=T)
# r_d06 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D06_PINK1_D06/2021-04-01_TO_EXP_each_condition_Control_D06_PINK1_D06.tsv",header=T)
# r_d15 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D15_PINK1_D15/2021-04-01_TO_EXP_each_condition_Control_D15_PINK1_D15.tsv",header=T)
# r_d21 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D21_PINK1_D21/2021-04-01_TO_EXP_each_condition_Control_D21_PINK1_D21.tsv",header=T)



# r_d00_f <- r_d00[r_d00$p_val_adj<0.05 & abs(r_d00$avg_logFC2) > 0.5,]
# r_d00_f_c <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2>0])
# r_d00_f_p <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2<0])

# r_d06_f <- r_d06[r_d06$p_val_adj<0.05 & abs(r_d06$avg_logFC2) > 0.5,]
# r_d06_f_c <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2>0])
# r_d06_f_p <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2<0])

# r_d15_f <- r_d15[r_d15$p_val_adj<0.05 & abs(r_d15$avg_logFC2) > 0.5,]
# r_d15_f_c <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2>0])
# r_d15_f_p <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2<0])

# r_d21_f <- r_d21[r_d21$p_val_adj<0.05 & abs(r_d21$avg_logFC2) > 0.5,]
# r_d21_f_c <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2>0])
# r_d21_f_p <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2<0])


# control_m<- unique(c(r_d00_f_c,r_d06_f_c,r_d15_f_c,r_d21_f_c))
# length(r_d00_f_c)
# length(r_d06_f_c)
# length(r_d15_f_c)
# length(r_d21_f_c)
# length(control_m)
# pink_m <- unique(c(r_d00_f_p,r_d06_f_p,r_d15_f_p,r_d21_f_p))
# length(r_d00_f_p)
# length(r_d06_f_p)
# length(r_d15_f_p)
# length(r_d21_f_p)
# length(pink_m)

# library(clusterProfiler)
# library(enrichplot)
# library(ReactomePA)
# library(org.Hs.eg.db)

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
# GOclusterplot <- 0
# Keggclusterplot<-0
# Reactomeclusterplot<-0
# GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
# Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
# Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


# pdf("GO_Enrichment_DF_Clusters_log2FC_0.5.pdf",width=22,height=11)
# dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dev.off()








# control_m_e <- bitr(r_d06_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# pink_m_e <- bitr(r_d06_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

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


# pdf("D06_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
# dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dev.off()



# control_m_e <- bitr(r_d15_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# pink_m_e <- bitr(r_d15_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


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


# pdf("D15_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
# dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dev.off()



# control_m_e <- bitr(r_d21_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# pink_m_e <- bitr(r_d21_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

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


# pdf("D21_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
# dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
# dev.off()

