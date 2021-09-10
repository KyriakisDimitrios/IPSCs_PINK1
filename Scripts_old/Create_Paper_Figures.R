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
Controls_DM <-ScaleData(Controls_DM,last)
pdf("Figure3b.pdf",width= 8)
DoHeatmap(Controls_DM,last,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()







genes_f <- c(
"SHH"
"FGF8"
"PAX2"
"WNT1"
"EN2"
"MSX1"
"FERD3L"
"GLI1"
"SNCA",
"RELN",
"CTNNB1",
"VIP",
"DMRTA2",
"SOX2",
"HES1",
"HES5",
"ADH1B",
"CALB1",
"EPHA5",
"NTN1",
"EBF1",
'HES1',
"NKX2-1",
"SLIT1",
"SLIT2",
"ROBO1",
"ROBO2",
"SEMA3A")
pdf("Figure3b2.pdf",width= 8)
Controls_DM<-ScaleData(Controls_DM,genes_f)
DoHeatmap(Controls_DM,genes_f,cells=Controls_DM$Timepoints!="Day10",group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()



# -----------------------------------------------------------

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
# "ROBO1"	,
# "ROBO2"	,
# "NTN1"	,
# "SLIT1"	,
# "OTX2"	,
# "KCNH6" ,	
"PITX2",
"ASCL1",
"PITX3"	,
"NEUROD1",
"CORIN",
"SLIT1",
"TH",
"NR4A2"
# "CORIN",
# "NEUROD1",
# "TUBB3",
# "NR4A2"	
)
pdf("Figure3b2.pdf",width= 8)
DefaultAssay(Controls) <- "RNA"
Controls<-ScaleData(Controls,mDA2)
DoHeatmap(Controls,mDA2,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")
dev.off()


mDA21 <- c("FABP7"	,
"SOX2"	,
"FOXA2"	,
"LMX1A"	,
"OTX2"	,
"WNT5A"	,
"RSPO2"	,
"MSX1"	,
"CORIN"	,
"ASCL1"	,
	      "NEUROG2",
	       "NEUROD1",
	       "TUBB3",
"DDC"	,
"DCX"	,
"NR4A2"	,
"PBX1"	,
"PITX3"	,
"EN1"	,
"TH"	,
'BNC2'	,
"SLC18A2"	,
	      "SLC6A3",
	       "CALB1",
	       "LMO3",
	      #"ALDH1A1",
"SOX6"	,
"SHH"	,
"MSX1"	,
"FERD3L",	
"CTNNB1",	
"DMRTA2",	
"HES1"	,
"SLIT2"	,
"ROBO1"	,
"ROBO2"	,
"NTN1"	,
"SLIT1"	,
"OTX2"	,
"PTCH1"	,
"FZD7"	,
"KCNH6" ,	
"TCF12"	,
"ALCAM",
"PITX2"
)
Controls_DM 

# ======================= FIGURE 4 ==========================
mDA <- c("TCF12", "ALCAM", "PITX2", "DDC", "ASCL1")
mDA <- mDA[mDA %in% rownames(Combined@assays$RNA@counts)]
p0 <- DimPlot(Combined,group.by = "condition",cols=color_cond)
#FeaturePlot(Combined,features = c("TH","KCNJ6"),pt.size = 2,,order = TRUE,cols = c("lightgrey","#FDBB84","#EF6548","#D7301F","#B30000","#7F0000"))
results <- ICSWrapper::scatter_gene(object = Combined,features = c("TH","KCNJ6"),ncol =1 ,nrow =2)+ theme_void()
plot4 <- DoHeatmap(Combined,mDA,group.by = "Timepoints", raster = FALSE)+  
scale_fill_viridis(option="inferno")
DefaultAssay(Combined) <- "SCT"
Combined <- ScaleData(Combined,mDA2)
plot4.1 <- DoHeatmap(Combined,mDA2,group.by = "Timepoints", raster = FALSE)+  
scale_fill_viridis(option="inferno")
DefaultAssay(Combined) <- "RNA"
Combined <- ScaleData(Combined,mDA2)
plot4.1  <- DoHeatmap(Combined,mDA2,group.by = "Timepoints", raster = FALSE)+  
scale_fill_viridis(option="inferno")

pdf("Figure4c.pdf",width=12,height=10)
plot4.1 
dev.off()


plot5 <- grid.arrange(p0, results, widths=c(2/3, 1/2), ncol=2, nrow=1)

plot5 <- ggarrange(plotlist=list(p0, results), widths=c(2/3, 1/2), ncol=2, nrow=1)

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