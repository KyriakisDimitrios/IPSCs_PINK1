
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

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
olor_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/3.Venn_Pairwise/")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/3.Venn_Pairwise")


Combined <- readRDS("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess/IPSCs_Combined.rds")


# ===================== OPEN FILES TAKE THE P.Adj- FC Genes
dirs_pairs <- list.dirs("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/2.Differential_Epression/Pairwise/",full.names = TRUE )[-1]
dirs_pairs <- list.dirs("C:/Users/dimitrios.kyriakis/Desktop/PINK1/NEW/2.1/Control_IPSCs_PINK1_IPSCs/Control_D10_PINK1_D06/",full.names = TRUE )[-1]

dirs_pairs <- grep('IPSC|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)

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





myCol <- brewer.pal(3, "Pastel2")
pdf("Summary_venn_diagramm.pdf")
day06 <- c(df_return_nt_cntrl[[1]],df_return_nt_pink[[1]])
day15 <- c(df_return_nt_cntrl[[2]],df_return_nt_pink[[2]])
day21 <- c(df_return_nt_cntrl[[3]],df_return_nt_pink[[3]])

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









# ======================= Enrichment ==================================
control_m <-unique(unlist(df_return_nt_cntrl))
pink_m <-  unique(unlist(df_return_nt_pink))

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




pdf("GO_Enrichment_DF_Clusters.pdf",width=22,height=11)
# dotplot(GOclusterplot)
# dotplot(Keggclusterplot)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()
# --------------------

unlist(df_return_nt_cntrl)

