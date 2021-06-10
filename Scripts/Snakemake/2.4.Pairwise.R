
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


Combined_Sub1 <- subset(Combined,subset=condition%in%c("Control_D21","PINK1_D21"))
Combined_Sub2 <- subset(Combined,subset=condition%in%c("Control_D21","PINK1_D21"))
Combined_Sub2[["RNA"]]@counts<-as.matrix(Combined_Sub2[["RNA"]]@counts)+1
Combined_Sub1 <- NormalizeData(Combined_Sub1)
Combined_Sub2 <- NormalizeData(Combined_Sub2)

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

pbmc.markers3 <- FindMarkers(object = Combined_Sub2,ident.1 = "Control_D21",ident.2 = "PINK1_D21",
                                   assay ="RNA",min.pct =0.1,
                                   logfc.threshold=0.0,
                                   only.pos = FALSE,
                                   test.use = "DESeq2",pseudocount.use = 0)

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



b11 <- log(x = mean(x = expm1(x = b1)) +1)
b22 <- log(x = mean(x = expm1(x = b2)) +1)
b11-b22

log(mean(expm1(b1)))-log(mean(expm1(b2)))
(b11-b22)/log(2)

b1 <- Combined_Sub1[["RNA"]]@data[c("DLK1"),Combined_Sub1$condition=="Control_D21"]
b2 <- Combined_Sub1[["RNA"]]@data[c("DLK1"),Combined_Sub1$condition=="PINK1_D21"]

log2(mean(b1)/mean(b2))
log(mean(b1)/mean(b2))

pdf("Volcanos.pdf")
pbmc.markers0$avg_logFC2 <-  pbmc.markers0$avg_logFC/log(2)
min_vlc <- min(pbmc.markers0$avg_logFC2)
max_vlc <- max(pbmc.markers0$avg_logFC2)
EnhancedVolcano(pbmc.markers0,
                lab = rownames(pbmc.markers0),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("Wilcox"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("Wilcox")


pbmc.markers01$avg_logFC2 <-  pbmc.markers01$avg_logFC/log(2)
min_vlc <- min(pbmc.markers01$avg_logFC2)
max_vlc <- max(pbmc.markers01$avg_logFC2)
EnhancedVolcano(pbmc.markers01,
                lab = rownames(pbmc.markers01),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("Wilcox (Pseud+1)"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("Wilcox (Pseud+1)")


pbmc.markers1$avg_logFC2 <-  pbmc.markers1$avg_logFC/log(2)
min_vlc <- min(pbmc.markers1$avg_logFC2)
max_vlc <- max(pbmc.markers1$avg_logFC2)
xlim = c(min_vlc-0.5,max_vlc +0.5)
EnhancedVolcano(pbmc.markers1,
                lab = rownames(pbmc.markers1),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("(FCcutoff=1)"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("MAST")

pbmc.markers12$avg_logFC2 <-  pbmc.markers12$avg_logFC/log(2)
min_vlc <- min(pbmc.markers12$avg_logFC2)
max_vlc <- max(pbmc.markers12$avg_logFC2)
EnhancedVolcano(pbmc.markers12,
                lab = rownames(pbmc.markers12),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("(FCcutoff=1)"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("MAST (Pseud+1)")

pbmc.markers3$avg_logFC2 <-  pbmc.markers3$avg_logFC/log(2)
min_vlc <- min(pbmc.markers3$avg_logFC2)
max_vlc <- max(pbmc.markers3$avg_logFC2)
EnhancedVolcano(pbmc.markers3,
                lab = rownames(pbmc.markers3),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("Not normalized (FCcutoff=1)"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 1)+ggtitle("DESeq2 (Pseud+1)")
dev.off()


comparisons <- list(c("Control_D21", "PINK1_D21"))

pdf("VLN.pdf")
gene_sig <- rownames(pbmc.markers0[order(pbmc.markers0$p_val_adj,decreasing = F),][1:4,])
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig1", test_sign = comparisons, y_max = 7)
gene_sig <- rownames(pbmc.markers01[order(pbmc.markers01$p_val_adj,decreasing = F),][1:4,])
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig2", test_sign = comparisons, y_max = 5)
gene_sig <- rownames(pbmc.markers1[order(pbmc.markers1$p_val_adj,decreasing = F),][1:4,])
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig3", test_sign = comparisons, y_max = 5)
gene_sig <- rownames(pbmc.markers12[order(pbmc.markers12$p_val_adj,decreasing = F),][1:4,])
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig4", test_sign = comparisons, y_max = 5)
gene_sig <- rownames(pbmc.markers3[order(pbmc.markers3$p_val_adj,decreasing = F),][1:4,])
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig5", test_sign = comparisons, y_max = 7)
dev.off()



pbmc.markers2
    pbmc.markers$gene <- rownames(pbmc.markers)



vp_case1 <- function(gene_signature, file_name, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(Combined_Sub1, features = signature,
            pt.size = 0.1, 
           group.by = "condition", 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
}

gene_sig <- c("MALAT1", "DLK1","CRABP1","MAP1B","NEFL")

gene_sig <- rownames(pbmc.markers0[order(pbmc.markers0$p_val_adj,decreasing = F),][1:4,])


pbmc.markers[rownames(pbmc.markers)%in% gene_sig,]
comparisons <- list(c("Control_D21", "PINK1_D21"))
vp_case1(gene_signature = gene_sig, file_name = "./gene_sig", test_sign = comparisons, y_max = 9)


FindMarkers(Combined_Sub, ..., test.use = "DESeq2", slot = "counts")


pbmc.markers2$avg_logFC2 <-  pbmc.markers2$avg_logFC/log(2)
min_vlc <- min(pbmc.markers2$avg_logFC2)
max_vlc <- max(pbmc.markers2$avg_logFC2)

pbmc.markers$avg_logFC2 <-  pbmc.markers$avg_logFC/log(2)
min_vlc <- min(pbmc.markers$avg_logFC2)
max_vlc <- max(pbmc.markers$avg_logFC2)
# debugonce(annotated_heat)
pdf("Volcano2.pdf")
# sub_pbm <- pbmc.markers[abs(pbmc.markers$avg_logFC2) > 0.1,]

plot(EnhancedVolcano(pbmc.markers,
                lab = rownames(pbmc.markers),
                x = 'avg_logFC2',
                y = 'p_val_adj',subtitle = paste("(FCcutoff=0.6)"),
                xlim = c(min_vlc-0.5,max_vlc +0.5),FCcutoff = 0.5))
dev.off()

pbmc.markers2[c("DLK1"),]
pbmc.markers[c("DLK1"),]


b1 <- Combined_Sub[["RNA"]]@counts[c("DLK1"),Combined_Sub$condition=="Control_D21"]
b2 <- Combined_Sub[["RNA"]]@counts[c("DLK1"),Combined_Sub$condition=="PINK1_D21"]
log2(mean(b1)/mean(b2))
log(mean(b1)/mean(b2))

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
mclapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,object=Combined,cl_combinations=cl_combinations,mc.cores=1)
setwd("../")
# ----------------------------------------------------------------------------------------------



# ====================================== COMMONS ===================================
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
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.05 & abs(ctrl_l1$avg_logFC2) >1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.05 & abs(pink_l1$avg_logFC2) >1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]]) 
    
    pdf(paste0(timp,"Volcanos.pdf"))
    l1$avg_logFC2 <- -l1$avg_logFC2
    l1$avg_logFC <-  -l1$avg_logFC   
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+ggtitle(timp))
    min_vlc <- min(l1$avg_logFC)
    max_vlc <- max(l1$avg_logFC)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC',
                        y = 'p_val_adj',subtitle = "FC (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change"))
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p.adj_BY',subtitle = "BY FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change")+ggtitle(timp))

    dev.off()

}


  
library(UpSetR )
glmer_res <- readRDS("../glmer_pval_mat_10percent.rds")
glmer_res$max_FC2 <- -glmer_res$max_FC/log(2)
glmer_res$FC2 <- -glmer_res$FC/log(2)
glmer_res$p_val_adj <- p.adjust(glmer_res$pvalues_l,method="bonferroni",n=dim(Combined@assays$RNA@counts)[1]) 

pdf("GLMER_Volcanos.pdf")
min_vlc <- min(glmer_res$max_FC2)
max_vlc <- max(glmer_res$max_FC2)
plot(EnhancedVolcano(glmer_res,
                        lab = rownames(glmer_res),
                        x = 'max_FC2',
                        y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Log2 fold change")+ggtitle("GLMM"))

dev.off()


glmm_genes  <- rownames(glmer_res[glmer_res$p_val_adj<0.05 & abs(glmer_res$max_FC2)> 1,])
df_return_nt_all[[5]] <- glmm_genes
names(df_return_nt_all) <- c("D06","D15","D21","IPSCs","GLMM")
pdf("Common.pdf",width=12)
UpSetR::upset(fromList(df_return_nt_all),text.scale =1.9)
dev.off()

unlist(df_return_nt_cntrl)

# # ============= Intersect Common Genes
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)

df_return_nt_cntrl[[1]]
df_return_nt_pink[[1]]


all_intesect <- Reduce(intersect, df_return_nt_all)
print(all_intesect)

pdf("All_venn_diagramm.pdf")
myCol <- brewer.pal(3, "Pastel2")
draw.triple.venn(area1 = 117, area2 = 110, area3 = 114, n12 = 23, n23 = 18, n13 = 10,
    n123 = 7, category = c("Day06", "Day15", "Day21"),
    fill = myCol)
dev.off()



# =========== Venn Diagrams DF genes
library(RColorBrewer)
library(VennDiagram)

myCol <- brewer.pal(4, "Pastel2")
pdf("Control_venn_diagramm.pdf")
day06 <- df_return_nt_cntrl[[1]]
day15 <- df_return_nt_cntrl[[2]]
day21 <- df_return_nt_cntrl[[3]]
IPSC <- df_return_nt_cntrl[[4]]



# Generate plot
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21,IPSC=IPSC),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
# have a look at the default plot
grid.newpage()
grid.draw(v)
# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

v[[14]]$label <- paste(intersect(intersect(intersect(day06, day15),day21),IPSC), collapse="\n")  
# plot  
grid.newpage()
grid.draw(v)
dev.off()


pdf("Pink1_venn_diagramm.pdf")
day06 <- df_return_nt_pink[[1]]
day15 <- df_return_nt_pink[[2]]
day21 <- df_return_nt_pink[[3]]
IPSC <- df_return_nt_pink[[4]]

# Generate plot
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21,IPSC=IPSC),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5,0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
# have a look at the default plot
grid.newpage()
grid.draw(v)
# have a look at the names in the plot object v
lapply(v,  names)
# We are interested in the labels
lapply(v, function(i) i$label)

v[[14]]$label <- paste(intersect(intersect(intersect(day06, day15),day21),IPSC), collapse="\n")  
# plot  
grid.newpage()
grid.draw(v)
dev.off()


# =========================================== ENRICHMENT ANALYSIS ==================================
r_d00 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_IPSCs_PINK1_IPSCs/2021-04-01_TO_EXP_each_condition_Control_IPSCs_PINK1_IPSCs.tsv",header=T)
r_d06 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D06_PINK1_D06/2021-04-01_TO_EXP_each_condition_Control_D06_PINK1_D06.tsv",header=T)
r_d15 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D15_PINK1_D15/2021-04-01_TO_EXP_each_condition_Control_D15_PINK1_D15.tsv",header=T)
r_d21 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/Pairwise/Control_D21_PINK1_D21/2021-04-01_TO_EXP_each_condition_Control_D21_PINK1_D21.tsv",header=T)



r_d00_f <- r_d00[r_d00$p_val_adj<0.05 & abs(r_d00$avg_logFC2) > 0.5,]
r_d00_f_c <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2>0])
r_d00_f_p <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2<0])

r_d06_f <- r_d06[r_d06$p_val_adj<0.05 & abs(r_d06$avg_logFC2) > 0.5,]
r_d06_f_c <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2>0])
r_d06_f_p <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2<0])

r_d15_f <- r_d15[r_d15$p_val_adj<0.05 & abs(r_d15$avg_logFC2) > 0.5,]
r_d15_f_c <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2>0])
r_d15_f_p <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2<0])

r_d21_f <- r_d21[r_d21$p_val_adj<0.05 & abs(r_d21$avg_logFC2) > 0.5,]
r_d21_f_c <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2>0])
r_d21_f_p <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2<0])


control_m<- unique(c(r_d00_f_c,r_d06_f_c,r_d15_f_c,r_d21_f_c))
length(r_d00_f_c)
length(r_d06_f_c)
length(r_d15_f_c)
length(r_d21_f_c)
length(control_m)
pink_m <- unique(c(r_d00_f_p,r_d06_f_p,r_d15_f_p,r_d21_f_p))
length(r_d00_f_p)
length(r_d06_f_p)
length(r_d15_f_p)
length(r_d21_f_p)
length(pink_m)

library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)

control_m_e <- bitr(control_m, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pink_m_e <- bitr(pink_m, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


n_c <- length(control_m_e$ENTREZID)
n_p <- length(pink_m_e$ENTREZID)
c_d <- data.frame(ENTREZID = control_m_e$ENTREZID , Condition ="Control")
p_d <- data.frame(ENTREZID = pink_m_e$ENTREZID , Condition ="PINK1")
all_d <- rbind(c_d,p_d)
mydf <- as.data.frame(all_d)
mydf$ENTREZID <- as.vector(mydf$ENTREZID)
mydf$ENTREZID <-as.numeric(mydf$ENTREZID )
GOclusterplot <- 0
Keggclusterplot<-0
Reactomeclusterplot<-0
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


pdf("GO_Enrichment_DF_Clusters_log2FC_0.5.pdf",width=22,height=11)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()








control_m_e <- bitr(r_d06_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pink_m_e <- bitr(r_d06_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

n_c <- length(control_m_e$ENTREZID)
n_p <- length(pink_m_e$ENTREZID)
c_d <- data.frame(ENTREZID = control_m_e$ENTREZID , Condition ="Control")
p_d <- data.frame(ENTREZID = pink_m_e$ENTREZID , Condition ="PINK1")
all_d <- rbind(c_d,p_d)
mydf <- as.data.frame(all_d)
mydf$ENTREZID <- as.vector(mydf$ENTREZID)
mydf$ENTREZID <-as.numeric(mydf$ENTREZID )
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


pdf("D06_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()



control_m_e <- bitr(r_d15_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pink_m_e <- bitr(r_d15_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


n_c <- length(control_m_e$ENTREZID)
n_p <- length(pink_m_e$ENTREZID)
c_d <- data.frame(ENTREZID = control_m_e$ENTREZID , Condition ="Control")
p_d <- data.frame(ENTREZID = pink_m_e$ENTREZID , Condition ="PINK1")
all_d <- rbind(c_d,p_d)
mydf <- as.data.frame(all_d)
mydf$ENTREZID <- as.vector(mydf$ENTREZID)
mydf$ENTREZID <-as.numeric(mydf$ENTREZID )
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


pdf("D15_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()



control_m_e <- bitr(r_d21_f_c, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
pink_m_e <- bitr(r_d21_f_p, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

n_c <- length(control_m_e$ENTREZID)
n_p <- length(pink_m_e$ENTREZID)
c_d <- data.frame(ENTREZID = control_m_e$ENTREZID , Condition ="Control")
p_d <- data.frame(ENTREZID = pink_m_e$ENTREZID , Condition ="PINK1")
all_d <- rbind(c_d,p_d)
mydf <- as.data.frame(all_d)
mydf$ENTREZID <- as.vector(mydf$ENTREZID)
mydf$ENTREZID <-as.numeric(mydf$ENTREZID )
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


pdf("D21_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()















# Combined21 <- subset(Combined,subset= condition %in% c("Control_D21","PINK1_D21"))
# DefaultAssay(Combined21)<-"RNA"
# pbmc <- NormalizeData(Combined21)
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# pbmc <- RunPCA(pbmc)

# pbmc <- FindNeighbors(pbmc, dims = 1:10)
# pbmc <- FindClusters(pbmc, resolution = 0.5)
# pbmc <- RunUMAP(pbmc, dims = 1:10)


# arnes <- list(
# c("FABP7","SOX2","WNT5A","CORIN","BNC2"),#Rgl3
# c("FABP7","SOX2","FOXA2","LMX1A","OTX2","WNT5A","RSPO2","MSX1","SOX6"),#Rgl1
# c("SOX2","FOXA2","LMX1A","OTX2","WNT5A","CORIN","DDC"),#ProgM
# c("SOX2","FOXA2","LMX1A","OTX2","WNT5A","DDC"),#ProgFPs
# c("SOX2","FOXA2","OTX2","WNT5A","ASCL1","NEUROG2","NEUROD1","TUBB3","DDC","DCX"),#Nprog
# c("FOXA2","NEUROG2","NEUROD1","TUBB3","DDC"),#NbM
# c("FOXA2","LMX1A","TUBB3","DDC","DCX","NR4A2","PBX1","PITX3","EN1","TH"),#DA0
# c("FOXA2","TUBB3","DDC","DCX","NR4A2","PBX1","PITX3","EN1","TH","BNC2","SLC18A2","SLC6A3","CALB1"),#DA1
# c("TUBB3","DDC","DCX","NR4A2","PBX1","PITX3","EN1","TH","BNC2","SLC18A2","SLC6A3","CALB1","LM03","ALDH1A1","SOX6")#DA2
# )

# pbmc <-AddModuleScore(pbmc,features = arnes,name = "arnes")
# pbmc@meta.data[,paste0("arnes",c(1:9))]


# pdf("arnes.pdf")
# DimPlot(pbmc, reduction = "umap")
# FeaturePlot(pbmc,features=paste0("arnes",c(1:9)),order=T)
# dev.off()






# skata <- read.csv("Table S2. Binarized Genes Across Cell Types, Manno et al 2016.csv",,header = T)
# skata[,1] <- as.vector(skata[,1])
# gene_l <-  list()
# for (i in 2:dim(skata)[2]){
#     gene_l[[i-1]] = skata[,1][skata[,i]==1]
# }


# types <- as.vector(colnames(skata)[2:27])

# Combined <-AddModuleScore(Combined,features = gene_l,name = "DA_s")
# scores <-  Combined@meta.data[,paste0("DA_s",c(1:26))]
# colnames(scores)<-types
# Combined_DA_category <- list()
# for(i in 1:dim(scores)[1]){
#     Combined_DA_category[[i]]  = types[which.max(scores[i,])]

#     scores[i,] <- scales::rescale(unlist(c(scores[i,])),c(0,1))
# }
# Combined$DA_category <- unlist(Combined_DA_category)


# library(ggradar)
# library(dplyr)
# library(scales)
# library(tibble)

# radar_data <- cbind(scores,Combined$DA_category)
# colnames(radar_data)[dim(radar_data)[2]] <- "DA_category"
# radar_data <- radar_data  %>% mutate_at(vars(-DA_category),rescale_max)

# radar_data_plot <- radar_data %>% group_by(DA_category) %>% summarise_all(funs(mean))
# radar_data_plot #<- radar_data_plot[,-dim(radar_data_plot)[2]]

# radar_data_plot %>%
#     mutate_at(vars(-"DA_category"), rescale)  %>%
#     ggradar()#font.radar = "Circular Air"

# #dim(Combined@meta.data)[2]-26
# #22
# # dim(Combined@meta.data)[2]-1

# colnames(Combined@meta.data)[22:47] <- types
# length(types)
# pdf("Scores.pdf")
# FeaturePlot(Combined,types[1:4],order=T)
# FeaturePlot(Combined,types[5:8],order=T)
# FeaturePlot(Combined,types[9:12],order=T)
# FeaturePlot(Combined,types[13:16],order=T)
# FeaturePlot(Combined,types[17:20],order=T)
# FeaturePlot(Combined,types[21:24],order=T)
# FeaturePlot(Combined,types[24:27],order=T)
# dev.off()

# FeaturePlot(Combined,paste0("DA_s",c(1:26)))
# names(gene_l) <- colnames(skata)

# ScaleData(Combined,types)
# DoHeatmap(Combined,types)


# skata[,1][skata[,2]==1]

# gene_l <- c("FABP7",
# "WNT5A",
# "RSPO2",
# "CORIN",
# "NEUROG2",
# "NEUROD1",
# "TUBB3",
# "DCX",
# "PBX1",
# "BNC2",
# "SLC6A3",
# "LM03",
# "SOX6")

# FeaturePlot(Combined,gene_l,order=T)
