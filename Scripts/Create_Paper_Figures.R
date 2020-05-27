# ======================== LOAD =========================
library("Seurat")
library("tidyverse")
library("viridis")
library("ggpubr")
library("gridExtra")
library("VennDiagram")
library("RColorBrewer")
Combined<- readRDS("Seurat.rds")
DefaultAssay(Combined) <- "SCT"
Combined$Timepoints <- as.vector(Combined$condition)
Combined$Timepoints[grep("IPSC",Combined$condition)] <-"iPSCs"
Combined$Timepoints[grep("D06",Combined$condition)] <-"Day06"
Combined$Timepoints[grep("D10",Combined$condition)] <-"Day10"
Combined$Timepoints[grep("D15",Combined$condition)] <-"Day15"
Combined$Timepoints[grep("D21",Combined$condition)] <-"Day21"
Combined$Timepoints <- factor(Combined$Timepoints,c("iPSCs","Day06","Day10","Day15","Day21"))

# -------------------------------------------------------


# ======================= FIGURE 3 ==========================
# === Figure a
dir.create("DF_Conditions")
setwd("DF_Conditions")
DefaultAssay(Combined) <-"RNA"
Combined <- NormalizeData(Combined)
Combined <- ScaleData(Combined)
Controls_DM <- subset(Combined,subset= Treatment !="PINK")
DefaultAssay(Controls_DM) <- "RNA"
Idents(Controls_DM)<-as.factor(Controls_DM$condition)
pbmc.markers <- FindAllMarkers(object = Controls_DM,
                            assay ="RNA",min.pct =0.1,
                            logfc.threshold=0.1,
                            only.pos = TRUE,
                            test.use = "MAST",latent.vars = c("nCount_RNA"))

pbmc.markers$gene <- rownames(pbmc.markers)
qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(Controls_DM@assays$RNA@counts)[1])
pbmc.markers$qvalue <- qvalue

top <- pbmc.markers[pbmc.markers$p_val_adj<0.01 & abs(pbmc.markers$avg_logFC) >0.1,]
to_heat <- top %>% group_by(cluster) %>% top_n(n = 15, wt = abs(avg_logFC)) 
to_heat <- to_heat[order(to_heat$cluster),]
Controls_DM$condition <- as.factor(Controls_DM$condition)
ICSWrapper::annotated_heat(object=Controls_DM,
                   row_annotation=c(1),
                   gene_list=to_heat$gene,
                   Rowv=NA,
                   gene_list_name="DF_genes",
                   title="Figure3a",
                   ordering="Treatment",One_annot = TRUE,color_list = color_list)
 

# === Figure b
Stemness_markers<- c("SOX2", "MYC", "POU5F1", "NANOG", "LIN28")
differentiation_path<- c("OTX2", "LMX1B", "LMX1A", "FOXA2", "PTCH1", "FZD7")
last <- c(Stemness_markers,differentiation_path)
pdf("Figure3b.pdf",width= 8)
DoHeatmap(Combined,last,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()
# -----------------------------------------------------------






# ======================= FIGURE 4 ==========================
mDA <- c("TCF12", "ALCAM", "PITX2", "DDC", "ASCL1")
mDA <- mDA[mDA %in% rownames(Combined@assays$RNA@counts)]
p0 <- DimPlot(Combined,group.by = "condition",cols=color_cond)
#FeaturePlot(Combined,features = c("TH","KCNJ6"),pt.size = 2,,order = TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"))
results <- ICSWrapper::scatter_gene(object = Combined,features = c("TH","KCNJ6"),ncol =1 ,nrow =2)+ theme_void()
plot4 <- DoHeatmap(Combined,mDA,group.by = "Timepoints", raster = FALSE)+  
scale_fill_viridis(option="inferno")
plot5 <- grid.arrange(p0, results, widths=c(2/3, 1/2), ncol=2, nrow=1)
pdf("Figure4.pdf",width=12,height=10)
grid.arrange(plot5,plot4,ncol=1, nrow=2,clip=TRUE)
dev.off()
# -----------------------------------------------------------


# ================== FIGURE 5 ==============================
color_list_three <- color_list
color_list_three$condition <- color_list_three$condition[-c(1,3,6)]
dirs_pairs <- list.dirs("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/Michi_Data/2020-04-17_seurat_elbow_TRUE_Mito-FALSE_Ribo-FALSE_SCT-TRUE_criteria_pass-3/DF_Pairwise",full.names = TRUE )[-1]
dirs_pairs <- grep('IPSC|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)
dirs_pairs <- dirs_pairs[-4]
df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()
for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
    l1 <- read.table(file,header=TRUE)
    l1$cluster <- l1$avg_logFC
    l1$cluster[ l1$avg_logFC<0] <- "PINK"
    l1$cluster[ l1$avg_logFC>0] <- "Control"
    ctrl_l1 <- l1[grep("Control",l1$cluster),]
    pink_l1 <- l1[grep("PINK",l1$cluster),]
	all_l1 <-  l1
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.01 & abs(ctrl_l1$avg_logFC) >0.1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.01 & abs(pink_l1$avg_logFC) >0.1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]])
}

# # ============= Intersect Common Genes
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)
length(cntrl_intesect)
length(pink_intesect)


Conserved_M <- subset(Combined,subset= Timepoints %in% c("Day06","Day15","Day21"))

Conserved_M$condition<- factor(Conserved_M$condition)
ICSWrapper::annotated_heat(object=Conserved_M,
                   row_annotation=c(1),
                   gene_list=c(cntrl_intesect,pink_intesect),
                   Rowv=NA,
                   gene_list_name="DF_genes",
                   title="Figure5a",
                   ordering="Treatment",One_annot = TRUE,color_list =  color_list_three)


# ============================= Figure 5b
df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()
for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
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
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
myCol <- brewer.pal(3, "Pastel2")
# === COntrol
pdf("Control_venn_diagramm.pdf")
day06 <- c(df_return_nt_cntrl[[1]])
day15 <- c(df_return_nt_cntrl[[2]])
day21 <- c(df_return_nt_cntrl[[3]])
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)
lapply(v,  names)
lapply(v, function(i) i$label)
v[[11]]$label <- paste(intersect(intersect(day06, day15),day21), collapse="\n")  
grid.newpage()
grid.draw(v)
dev.off()
v1 <-v
#------------

# === PINK
pdf("PINK_venn_diagramm.pdf")
day06 <- c(df_return_nt_pink[[1]])
day15 <- c(df_return_nt_pink[[2]])
day21 <- c(df_return_nt_pink[[3]])
v <- venn.diagram(list(Day06=day06, Day15=day15,Day21=day21),
                  fill = myCol,
                  alpha = c(0.5, 0.5, 0.5), cat.cex = 1.5, cex=1.5,
                  filename=NULL)
grid.newpage()
grid.draw(v)
lapply(v,  names)
lapply(v, function(i) i$label)

v[[11]]$label <- paste(intersect(intersect(day06, day15),day21), collapse="\n")  
grid.newpage()
p1 <-grid.draw(v)
dev.off()
v2 <-v
# ------------------ 

# ===== Plot
library(grid)
library(gridBase)
graphics.off()
pdf("Figure5b.pdf")
plot.new()
vp1 <- viewport(x=0,y=0.5,width=0.5, height=0.5, just = c("left", "bottom"))
vp2 <- viewport(x=0,y=0,width=0.5, height=0.5, just = c("left", "bottom"))
upViewport()
pushViewport(vp1)
grid.draw(v1)
par(new=TRUE, fig=gridFIG())

upViewport()
pushViewport(vp2)
grid.draw(v2)
par(new=TRUE, fig=gridFIG())
dev.off()

# -----------------------------------------------------------




# ==================== Correlation Network ==================
# -----------------------------------------------------------