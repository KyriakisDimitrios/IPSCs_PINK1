
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

Combined_Sub1 <- subset(Combined,subset=Timepoints%in%c("IPSCs","D06","D15","D21"))
Combined_Sub1 <- NormalizeData(Combined_Sub1)

glmm_mat <- as.data.frame(t(Combined_Sub1@assays$RNA@data))
glmm_mat_r <- as.data.frame(t(Combined_Sub1@assays$RNA@counts))
dim(glmm_mat)
glmm_mat_f <- glmm_mat[,colSums(glmm_mat_r>0) > (3893*10/100)]
m1 <- colnames(glmm_mat_f)

glmm_mat_f$Timepoints <- Combined_Sub1$Timepoints
glmm_mat_f$Timepoints <- Combined_Sub1$Timepoints
table(Combined_Sub1$Treatment)
Combined_Sub1$Treatment <- factor(Combined_Sub1$Treatment,levels=c("Control","PINK"))

glmm_mat_f$Treatment <- Combined_Sub1$Treatment

library(pbapply)
library(lme4)
library(parallel)
cl <- makeCluster(2L)
clusterExport(cl,c("glmm_mat_f","m1"))
## parallel with progress bar: snow type cluster
## (RNG is set in the main process to define the object bid)
system.time(res1pbcl <- pblapply(m1,function(feat) {
    formula <- paste("Treatment ~ " , feat , " + (1 | Timepoints)" )
    lmm <- lme4::glmer(formula, data = glmm_mat_f, family = binomial)
    sum_lmm <- summary(lmm)
    a <- as.data.frame(sum_lmm$coefficients)
    a$`Pr(>|z|)`[2]
}, cl = cl))

pvalues_l <- res1pbcl
names(pvalues_l) <- m1
adj_pval_BH <- unlist(p.adjust(pvalues_l,method="BH"))
adj_pval_bondferroni <- unlist(p.adjust(pvalues_l,method="bonferroni"))
adj_pval_BY <- unlist(p.adjust(pvalues_l,method="BY"))

pval_mat <- cbind(pvalues_l,adj_pval_BY,adj_pval_BH,adj_pval_bondferroni)
pval_mat <- as.data.frame(pval_mat)

pval_mat$adj_pval_BY <- unlist(pval_mat$adj_pval_BY)
pval_mat$adj_pval_bondferroni <- unlist(pval_mat$adj_pval_bondferroni)
pval_mat$adj_pval_BH <- unlist(pval_mat$adj_pval_BH)
rownames(pval_mat) <- m1



b <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$Treatment=="Control"]
bp <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$Treatment=="PINK"]
b_c <- apply(b,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
bp_c <- apply(bp,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})




# ======= Substract genes ======
b00 <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="Control_IPSCs"]
b00p <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="PINK1_IPSCs"]
b06 <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="Control_D06"]
b06p <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="PINK1_D06"]
b15 <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="Control_D15"]
b15p <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="PINK1_D15"]
b21 <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="Control_D21"]
b21p <- Combined_Sub1[["RNA"]]@data[m1,Combined_Sub1$condition=="PINK1_D21"]
# ---------------------------------

# ========= Calculate FC =========
b00_c <- apply(b00,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b00p_c <- apply(b00p,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b06_c <- apply(b06,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b06p_c <- apply(b06p,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b15_c <- apply(b15,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b15p_c <- apply(b15p,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b21_c <- apply(b21,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})
b21p_c <- apply(b21p,1,function(b_s) {log(x = mean(x = expm1(x = b_s)) +1)})

pval_mat$FC_D00 <- b00_c-b00p_c
pval_mat$FC_D06 <- b06_c-b06p_c
pval_mat$FC_D15 <- b15_c-b15p_c
pval_mat$FC_D21 <- b21_c-b21p_c
pval_mat$FC <- b_c-bp_c
# ---------------------------------

# ============ MAX FC =============
max_fc_l <- apply(pval_mat[,c("FC_D00","FC_D06","FC_D15","FC_D21")],1,function(gene) {
    max_fc <- gene[which.max(abs(gene))]
    max_fc
})
pval_mat$max_FC <- unlist(max_fc_l)
# ---------------------------------



saveRDS(pval_mat,"glmer_pval_mat_10percent.rds")


stopCluster(cl)
# -------------------------------------------------------------------------


# ============================= PLOT VOLCANOS ==========================
pdf("Volc_glmer.pdf")

min_vlc <- min(pval_mat$FC_D00)
max_vlc <- max(pval_mat$FC_D00)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC_D00',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("FC: IPSCs")

min_vlc <- min(pval_mat$FC_D06)
max_vlc <- max(pval_mat$FC_D06)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC_D06',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("FC: D06")

min_vlc <- min(pval_mat$FC_D15)
max_vlc <- max(pval_mat$FC_D15)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC_D15',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("FC: D15")

min_vlc <- min(pval_mat$FC_D21)
max_vlc <- max(pval_mat$FC_D21)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC_D21',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("FC: D21")

min_vlc <- min(pval_mat$max_FC)
max_vlc <- max(pval_mat$max_FC)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'max_FC',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("MAX FC")



min_vlc <- min(pval_mat$FC)
max_vlc <- max(pval_mat$FC)
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 1,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("Overall FC")
EnhancedVolcano(pval_mat,
                lab = rownames(pval_mat),
                x = 'FC',
                y = 'adj_pval_BY',
                xlim = c(min_vlc-0.5,max_vlc +0.5),
                FCcutoff = 0.5,subtitle ="Treatment ~ feat + (1 | Timepoints)")+ggtitle("Overall FC")
toplot <- pval_mat[pval_mat$adj_pval_BY<0.05 & abs(pval_mat$FC)>0.4 ,]
DoHeatmap(Combined_Sub1,rownames(toplot[order(toplot$FC,decreasing=T),]),size = 4.5) + 
    theme(text = element_text(size = 9))
dev.off()






# ======================= Enrichment ==================================
sig_genes <- pval_mat[pval_mat$adj_pval_BY < 0.05& abs(pval_mat$FC)>0.1,]
control_m <- rownames(sig_genes[sig_genes$FC < 0,] )
pink_m <-  rownames(sig_genes[sig_genes$FC > 0,])
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
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = mydf$ENTREZID,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

pdf("GLME_GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
# dotplot(GOclusterplot)
# dotplot(Keggclusterplot)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()

