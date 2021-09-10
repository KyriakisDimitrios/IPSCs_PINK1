## ----setup, include=FALSE--------------------------------------------------------------------
set.seed(123)
# library(reticulate)
options(future.globals.maxSize= 2122317824)
library(sctransform)
library(Seurat)
library(Routliers)
library(cluster)
library(NMF)

library(RColorBrewer)
library(jcolors)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(igraph)
library(gridExtra)
library(grid)
library(crayon)

library(tictoc)
library(stringr)
library(plyr)
library(jsonlite)
library(purrr)
library(data.table)

library("Hmisc")
require(tidyverse)
require(corrr)
require(qgraph)


colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

# ================ READ ARGS =================
args = commandArgs(trailingOnly=TRUE)
print(args)

source("workflow/scripts/Functions.R")
source("workflow/scripts/rstring.r")

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
    print("Arguments Passed")
}
input_file <- args[1]
csv1 <- args[2]
csv2<- args[3]
csv3 <- args[4]
csv4 <- args[5]
output_file <- args[6]
# ---------------------------------------------

Combined <- readRDS(input_file)
DefaultAssay(Combined) <- "RNA"
Combined<- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))
Control<-subset(Combined,subset = Treatment=="Control")
PINK1<-subset(Combined,Treatment=="PINK1") 

# ========================================= DF NETWORK ====================================================
# ========================================= DF NETWORK ====================================================
# ========================================= DF NETWORK ====================================================
dataset <- as.data.frame(Combined@assays$RNA@data)

net_analysis="extented"
net_analysis="net"
# DF
graph_annotation <- read.csv(csv1)
graph_annotation$group1 <- as.vector(graph_annotation$group1)
graph_annotation$group1[graph_annotation$group1==""] <- "NA"
first_graph <- read.csv(csv2)


# =============== ANNOTATION ========================
ann_genes <- c(as.vector(graph_annotation$Id))
ann_genes <- toupper(ann_genes)
ann_genes <- gsub("-", ".", ann_genes, fixed = TRUE)
ann_genes[ann_genes=="ENSG00000173575"] <- "CHD2" 
ann_genes[ann_genes=="ENSG00000279576"] <- "MALAT1" 
ann_genes[ann_genes=="GPR128"] <- "ADGRG7" 
ann_genes[ann_genes=="DJ.1"] <- "PARK7" 
ann_genes[ann_genes=="PARKIN"] <- "PARK2" 
ann_genes[ann_genes=="PD2"] <- "PAF1"
ann_genes[ann_genes=="ERV1"] <- "GFER"
ann_genes[ann_genes=="RM1"] <- "TIPARP"
# ----------------------------------------------------

# ================ NETWORK ========================
first_graph$Source <- toupper(gsub("-", ".", first_graph$Source , fixed = TRUE))
first_graph$Target <- toupper(gsub("-", ".", first_graph$Target , fixed = TRUE))
# ====== SOURCE
first_graph$Source[first_graph$Source=="ENSG00000173575"] <- "CHD2"
first_graph$Source[first_graph$Source=="ENSG00000279576"] <- "MALAT1"
first_graph$Source[first_graph$Source=="GPR128"] <- "ADGRG7"
first_graph$Source[first_graph$Source=="DJ.1"] <- "PARK7"
first_graph$Source[first_graph$Source=="PARKIN"] <- "PARK2"
first_graph$Source[first_graph$Source=="PD2"] <- "PAF1"
first_graph$Source[first_graph$Source=="ERV1"] <- "GFER"
first_graph$Source[first_graph$Source=="RM1"] <- "TIPARP"
#=== Target
first_graph$Target[first_graph$Target=="ENSG00000173575"] <- "CHD2"
first_graph$Target[first_graph$Target=="ENSG00000279576"] <- "MALAT1"
first_graph$Target[first_graph$Target=="GPR128"] <- "ADGRG7"
first_graph$Target[first_graph$Target=="DJ.1"] <- "PARK7"
first_graph$Target[first_graph$Target=="PARKIN"] <- "PARK2"
first_graph$Target[first_graph$Target=="PD2"] <- "PAF1"
first_graph$Target[first_graph$Target=="ERV1"] <- "GFER"
first_graph$Target[first_graph$Target=="RM1"] <- "TIPARP"
# ----------------------------------------------------

# first_graph <- first_graph[first_graph$Target%in%df151_genes & first_graph$Source %in% df151_genes,]


f_g_genes <- unique(c(as.vector(first_graph$Source),as.vector(first_graph$Target))) 
f_g_genes <- toupper(f_g_genes)
f_g_genes <- gsub("-", ".", f_g_genes, fixed = TRUE)
r_f_g_genes <- f_g_genes[f_g_genes%in% rownames(dataset)]
cat(paste("Genes not in dataset:",
          length(f_g_genes)-length(r_f_g_genes)))
missing <- setdiff(f_g_genes, r_f_g_genes)
cat(missing)
# # res_mis <- GeneSymbolThesarus(missing)
# match(missing,annot_29$Id)
# match(missing,annot_29$Label)


r_first_graph <- first_graph[first_graph$Source %in% r_f_g_genes,]
r2_first_graph <- r_first_graph[r_first_graph$Target %in% r_f_g_genes,]
dim(first_graph)
dim(r2_first_graph)
dataset_r <-dataset[r_f_g_genes,]

# === cOR NET
retur_graph <- cor_net(dataset = dataset_r, cor_thres=0.1)
g2 <- retur_graph
# retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.1)
# print(retur_graph)
# g2 <- retur_graph$graph
print(g2)

# retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.2)
# g_cor2 <- retur_graph$graph

retur_graph <- cor_net(dataset = dataset_r, cor_thres=0.2)
g_cor2 <- retur_graph
print(g_cor2)



# retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.3)
# g_cor3 <- retur_graph$graph
retur_graph <- cor_net(dataset = dataset_r, cor_thres=0.3)
g_cor3 <- retur_graph
print(g_cor3)



E(g2)$weight <- rep(1,length(as_ids(E(g2))))
edge_cor2 <- as_ids(E(g2))[as_ids(E(g2))%in%as_ids(E(g_cor2))]
E(g2)$weight[match(edge_cor2,as_ids(E(g2)))]<-2.5
edge_cor3 <- as_ids(E(g2))[as_ids(E(g2))%in%as_ids(E(g_cor3))]
E(g2)$weight[match(edge_cor3,as_ids(E(g2)))]<-5

# ----------

# ===================== COLOR PALETTE ==============================
color_clust<-categorical_pal(8)
color_clust[5] <- color_clust[7]
color_clust[1] <-"#0072B2"
color_clust[6] <- "red"
color_clust[2] <- "gray"

g1 <- graph_from_data_frame(first_graph,directed = FALSE)
color_facts <- as.factor(graph_annotation$group1[match(as_ids(V(g1)),ann_genes)])
Color_rest <-color_clust[color_facts]
Color_rest[Color_rest=="#D55E00"] <- "red"
V(g1)$cex <- 0.05
V(g1)$label.cex <-0.25
V(g1)$label.color="black"
V(g1)$color <- Color_rest

set.seed(1234)
g1<- igraph::simplify(g1)
# ============================


# ========================== PLOT ===================================
pdf("result/Network/Network_DE_only.pdf",width=15,height=15)
set.seed(123)
locs <- layout_on_sphere(g1)*0.7

e <- get.edgelist(g1,names=FALSE)
locs <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g1),
                                          area=50*(vcount(g1)^3),repulse.rad=(vcount(g1)^4.1))



col_added<- unlist(lapply(color_facts,FUN = function(x){ if(x %in% c("A","B","C","B-PD")){return(1)}else{return(2)}}))    

V(g1)$color <- Color_rest

plot(g1,  vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=unlist(lapply(V(g1)$name,nchar))*2.8,vertex.size2=5, vertex.label.cex=0.98,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F,cex=1.5)
FigSup_1<- recordPlot()
plot.new() ## clean up device


bet<-betweenness(g1)
n <- vcount(g1)
test_bet <- (2 * bet) / (n*n - 3*n + 2) * 500

test_bet<-scale(bet,center = 0)+4

plot(g1, vertex.label.color="black",vertex.label.font=11,
     vertex.size2=unlist(lapply(V(g1)$name,nchar))*2.5, vertex.label.cex=1,vertex.size=test_bet,
     layout=locs, 
     main="Original")
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F,cex=1.5)
dev.off()


# Remove a few vertices
g3_c <- induced_subgraph(g1, as_ids(V(g2)))
test1<-delete_edges(g3_c, E(g3_c))
test2<-add_edges(test1,get.edgelist(g2))
test3<- igraph::simplify(test2, remove.multiple = TRUE, remove.loops = TRUE)

# Add weight
E(test3)$weight <- rep(1.0,length(as_ids(E(test3))))
g4_edge_cor2 <- as_ids(E(test3))[as_ids(E(test3))%in%as_ids(E(g_cor2))]
E(test3)$weight[match(g4_edge_cor2,as_ids(E(test3)))]<-2.5
g4_edge_cor3 <- as_ids(E(test3))[as_ids(E(test3))%in%as_ids(E(g_cor3))]
E(test3)$weight[match(g4_edge_cor3,as_ids(E(test3)))]<-5.0
cor_r <- E(test3)$weight
cor_r[cor_r==1.0] <- 0.1
cor_r[cor_r==2.5] <- 0.2
cor_r[cor_r==5.0] <- 0.3


corr_edge_color <- cor_r
corr_edge_color[corr_edge_color==0.1] <- "gray"
corr_edge_color[corr_edge_color==0.2] <- "orange"
corr_edge_color[corr_edge_color==0.3] <- "red"

corr_edge_width <- cor_r
corr_edge_width[cor_r==0.1] <- 1
corr_edge_width[cor_r==0.2] <- 2.5
corr_edge_width[cor_r==0.3] <- 5


E(test3)$weight <- 1
E(test3)$color <- corr_edge_color
E(test3)$width <- corr_edge_width

## To check the result
pdf("result/Network/1.Correlation_Networks.pdf",width=15,height=15)
plot(test3, vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=25,
     vertex.size=unlist(lapply(V(test3)$name,nchar))*3.2,vertex.size2=6, vertex.label.cex=1, layout=locs[as_ids(V(g1)) %in% as_ids(V(test3)),])


legend("bottomleft",inset = c(0, 0),bty = "n",
       legend=levels(color_facts),
       fill= color_clust,
       horiz=F,title="Nodes",cex=1.5)

legend("bottomleft",inset = c(0.15, 0),
       bty = "n", title="Edges",
       legend=levels(as.factor(cor_r)),
       fill=  c("gray","orange","red"),
       horiz=T,cex=1.5)
FigSup_2<- recordPlot()
plot.new() ## clean up device

# ----------------------------------------------------------------------------------------------------


# =============================== COMMON EDGES =======================================================
common1 <-  g1%s%g2
common2 <-  graph_from_data_frame(as_edgelist(igraph::simplify(common1), names = TRUE),directed = FALSE)

common3 <- induced_subgraph(g1, as_ids(V(common2)))
common4<-delete_edges(common3, E(common3))
common5<-add_edges(common4,get.edgelist(common2))
common5<- igraph::simplify(common5, remove.multiple = TRUE, remove.loops = TRUE)

# Add weight
E(common5)$weight <- rep(1,length(as_ids(E(common5))))
g4_edge_cor2 <- as_ids(E(common5))[as_ids(E(common5))%in%as_ids(E(g_cor2))]
E(common5)$weight[match(g4_edge_cor2,as_ids(E(common5)))]<-2.5
g4_edge_cor3 <- as_ids(E(common5))[as_ids(E(common5))%in%as_ids(E(g_cor3))]
E(common5)$weight[match(g4_edge_cor3,as_ids(E(common5)))]<-5


cor_r <- E(common5)$weight
cor_r[cor_r==1.0] <- 0.1
cor_r[cor_r==2.5] <- 0.2
cor_r[cor_r==5.0] <- 0.3

common_edge_color <- cor_r
common_edge_color[common_edge_color==0.1] <- "gray"
common_edge_color[common_edge_color==0.2] <- "orange"
common_edge_color[common_edge_color==0.3] <- "red"

common_edge_width <- cor_r
common_edge_width[cor_r==0.1] <- 1
common_edge_width[cor_r==0.2] <- 2.5
common_edge_width[cor_r==0.3] <- 5


E(common5)$weight <- 1
E(common5)$color <- common_edge_color
E(common5)$width <- common_edge_width

plot(common5,vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=25,
     vertex.size=unlist(lapply(V(common5)$name,nchar))*2.5,vertex.size2=6, vertex.label.cex=1,
     layout=locs[as_ids(V(g1)) %in% as_ids(V(common5)),])
legend("bottomleft",inset = c(0.15, 0),
       bty = "n", title="Edges",
       legend=levels(as.factor(cor_r)),
       fill=  c("gray","orange","red"),
       horiz=F,cex=1.5)

legend("bottomleft",inset = c(0, 0),bty = "n",
       legend=levels(color_facts),
       fill= color_clust,
       horiz=F,title="Nodes",cex=1.5)
FigSup_3<- recordPlot()
plot.new() ## clean up device

dev.off()

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------



# ========================================= DF+ADDED NETWORK ====================================================
# ========================================= DF+ADDED NETWORK ====================================================
# ========================================= DF+ADDED NETWORK ====================================================
dataset <- as.data.frame(Combined@assays$RNA@data)
net_analysis="extented"
# Ubiq
graph_annotation <- read.csv(csv3)
# Net
first_graph <- read.csv(csv4)
graph_annotation$group1 <- as.vector(graph_annotation$group1)
graph_annotation$group1 <- unlist(lapply(graph_annotation$group1,FUN = function(x){ if(x %in% c("A","B","C","B-PD")){return("DE")}else{return(x)}}))

# =============== ANNOTATION ========================
ann_genes <- c(as.vector(graph_annotation$Id))
ann_genes <- toupper(ann_genes)
ann_genes <- gsub("-", ".", ann_genes, fixed = TRUE)
ann_genes[ann_genes=="ENSG00000173575"] <- "CHD2" 
ann_genes[ann_genes=="ENSG00000279576"] <- "MALAT1" 
ann_genes[ann_genes=="GPR128"] <- "ADGRG7" 
ann_genes[ann_genes=="DJ.1"] <- "PARK7" 
ann_genes[ann_genes=="PARKIN"] <- "PARK2" 
ann_genes[ann_genes=="PD2"] <- "PAF1"
ann_genes[ann_genes=="ERV1"] <- "GFER"
ann_genes[ann_genes=="RM1"] <- "TIPARP"
# ----------------------------------------------------

# ================ NETWORK ========================
first_graph$Source <- toupper(gsub("-", ".", first_graph$Source , fixed = TRUE))
first_graph$Target <- toupper(gsub("-", ".", first_graph$Target , fixed = TRUE))
# ====== SOURCE
first_graph$Source[first_graph$Source=="ENSG00000173575"] <- "CHD2"
first_graph$Source[first_graph$Source=="ENSG00000279576"] <- "MALAT1"
first_graph$Source[first_graph$Source=="GPR128"] <- "ADGRG7"
first_graph$Source[first_graph$Source=="DJ.1"] <- "PARK7"
first_graph$Source[first_graph$Source=="PARKIN"] <- "PARK2"
first_graph$Source[first_graph$Source=="PD2"] <- "PAF1"
first_graph$Source[first_graph$Source=="ERV1"] <- "GFER"
first_graph$Source[first_graph$Source=="RM1"] <- "TIPARP"
#=== Target
first_graph$Target[first_graph$Target=="ENSG00000173575"] <- "CHD2"
first_graph$Target[first_graph$Target=="ENSG00000279576"] <- "MALAT1"
first_graph$Target[first_graph$Target=="GPR128"] <- "ADGRG7"
first_graph$Target[first_graph$Target=="DJ.1"] <- "PARK7"
first_graph$Target[first_graph$Target=="PARKIN"] <- "PARK2"
first_graph$Target[first_graph$Target=="PD2"] <- "PAF1"
first_graph$Target[first_graph$Target=="ERV1"] <- "GFER"
first_graph$Target[first_graph$Target=="RM1"] <- "TIPARP"
# ----------------------------------------------------

# first_graph <- first_graph[first_graph$Target%in%df151_genes & first_graph$Source %in% df151_genes,]


f_g_genes <- unique(c(as.vector(first_graph$Source),as.vector(first_graph$Target))) 
f_g_genes <- toupper(f_g_genes)
f_g_genes <- gsub("-", ".", f_g_genes, fixed = TRUE)
r_f_g_genes <- f_g_genes[f_g_genes%in% rownames(dataset)]
cat(paste("Genes not in dataset:",
          length(f_g_genes)-length(r_f_g_genes)))
missing <- setdiff(f_g_genes, r_f_g_genes)
cat(missing)
# res_mis <- GeneSymbolThesarus(missing)
# match(missing,annot_29$Id)
# match(missing,annot_29$Label)


r_first_graph <- first_graph[first_graph$Source %in% r_f_g_genes,]
r2_first_graph <- r_first_graph[r_first_graph$Target %in% r_f_g_genes,]
dim(first_graph)
dim(r2_first_graph)
dataset_r <-dataset[r_f_g_genes,]

# === cOR NET
retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.1)
g2 <- retur_graph$graph

retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.2)
g_cor2 <- retur_graph$graph

retur_graph <- ics_net(dataset = dataset_r , method = "correlation",threshold = 0.3)
g_cor3 <- retur_graph$graph

E(g2)$weight <- rep(1,length(as_ids(E(g2))))
edge_cor2 <- as_ids(E(g2))[as_ids(E(g2))%in%as_ids(E(g_cor2))]
E(g2)$weight[match(edge_cor2,as_ids(E(g2)))]<-2.5
edge_cor3 <- as_ids(E(g2))[as_ids(E(g2))%in%as_ids(E(g_cor3))]
E(g2)$weight[match(edge_cor3,as_ids(E(g2)))]<-5

# ----------

# ===================== COLOR PALETTE ==============================
color_clust<-categorical_pal(8)

color_clust[5] <- color_clust[7]
color_clust[3] <- "red"

#color_clust<-sns.palplot(sns.color_palette("colorblind", 8))
color_clust<- c( "#E69F00","#0072B2", "red","#56B4E9","#009E73","#F0E442", "#D55E00", "#CC79A7","#999999")


g1 <- graph_from_data_frame(first_graph,directed = FALSE)
color_facts <- as.factor(graph_annotation$group1[match(as_ids(V(g1)),ann_genes)])
Color_rest <-color_clust[color_facts]
Color_rest[Color_rest=="#D55E00"] <- "red"
V(g1)$cex <- 0.05
V(g1)$label.cex <-0.25
V(g1)$label.color="black"
V(g1)$color <- Color_rest

set.seed(1234)
g1<- igraph::simplify(g1)
# ============================


# ========================== PLOT ===================================

pdf("result/Network/2.Network_extented.pdf",width=12,height=10)
set.seed(12345)
E(g1)$weight = 0.15
locs <- layout_with_fr(g1,dim=2, niter=1000)

col_added<- unlist(lapply(color_facts,FUN = function(x){ if(x %in% c("DE")){return(1)}else{return(2)}}))



col_added[col_added==1] <- "Differentially Expressed"
col_added[col_added==2] <- "Added"
col_added<-as.factor(col_added)
V(g1)$color <- col_added


plot(g1, vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=unlist(lapply(V(g1)$name,nchar))*2.5,vertex.size2=5, vertex.label.cex=0.7,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(col_added),fil=categorical_pal(8)[c(2,1)], horiz=F,cex=1.5)


V(g1)$color <- Color_rest
ela <- V(g1)$name
plot(g1,  vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=unlist(lapply(ela,nchar))*2.5,vertex.size2=4, vertex.label.cex=0.8,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F,cex=1.5)

FigSup_4<- recordPlot()
plot.new() ## clean up device


plot(g1,  vertex.label.color="black",vertex.shape="circle", vertex.size=5, vertex.label.cex=0.9,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F)

plot(g1,  vertex.label.color="black",vertex.shape="circle", vertex.size=5, vertex.label=NA,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F)


bet<-betweenness(g1)
n <- vcount(g1)
test_bet <- scale(bet,center = 0)+4

plot(g1, vertex.label.color="black",vertex.label.font=11,
     vertex.size2=5, vertex.label.cex=0.7,vertex.size=test_bet,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F)



# Remove a few vertices
g3_c <- induced_subgraph(g1, as_ids(V(g2)))
test1<-delete_edges(g3_c, E(g3_c))
test2<-add_edges(test1,get.edgelist(g2))
test3<- igraph::simplify(test2, remove.multiple = TRUE, remove.loops = TRUE)

# Add weight
E(test3)$weight <- rep(1.0,length(as_ids(E(test3))))
g4_edge_cor2 <- as_ids(E(test3))[as_ids(E(test3))%in%as_ids(E(g_cor2))]
E(test3)$weight[match(g4_edge_cor2,as_ids(E(test3)))]<-2.5
g4_edge_cor3 <- as_ids(E(test3))[as_ids(E(test3))%in%as_ids(E(g_cor3))]
E(test3)$weight[match(g4_edge_cor3,as_ids(E(test3)))]<-5.0
cor_r <- E(test3)$weight
cor_r[cor_r==1.0] <- 0.1
cor_r[cor_r==2.5] <- 0.2
cor_r[cor_r==5.0] <- 0.3


corr_edge_color <- cor_r
corr_edge_color[corr_edge_color==0.1] <- "gray"
corr_edge_color[corr_edge_color==0.2] <- "orange"
corr_edge_color[corr_edge_color==0.3] <- "red"

corr_edge_width <- cor_r
corr_edge_width[cor_r==0.1] <- 1
corr_edge_width[cor_r==0.2] <- 2.5
corr_edge_width[cor_r==0.3] <- 5


E(test3)$weight <- 1
E(test3)$color <- corr_edge_color
E(test3)$width <- corr_edge_width


## To check the result
plot(test3, vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=11,vertex.size2=5, vertex.label.cex=0.7, layout=locs[as_ids(V(g1)) %in% as_ids(V(test3)),])


legend("bottomleft",inset = c(0, 0),bty = "n",
       legend=levels(color_facts),
       fill= color_clust,
       horiz=F,title="Nodes")

legend("bottomleft",inset = c(0.15, 0),
       bty = "n", title="Edges",
       legend=levels(as.factor(cor_r)),
       fill=  c("gray","orange","red"),
       horiz=F)

# -------------------------------------------------------------------------



common1 <-  g1%s%g2
common2 <-  graph_from_data_frame(as_edgelist(igraph::simplify(common1), names = TRUE),directed = FALSE)

common3 <- induced_subgraph(g1, as_ids(V(common2)))
common4<-delete_edges(common3, E(common3))
common5<-add_edges(common4,get.edgelist(common2))
common5<- igraph::simplify(common5, remove.multiple = TRUE, remove.loops = TRUE)

# Add weight
E(common5)$weight <- rep(1,length(as_ids(E(common5))))
g4_edge_cor2 <- as_ids(E(common5))[as_ids(E(common5))%in%as_ids(E(g_cor2))]
E(common5)$weight[match(g4_edge_cor2,as_ids(E(common5)))]<-2.5
g4_edge_cor3 <- as_ids(E(common5))[as_ids(E(common5))%in%as_ids(E(g_cor3))]
E(common5)$weight[match(g4_edge_cor3,as_ids(E(common5)))]<-5


cor_r <- E(common5)$weight
cor_r[cor_r==1.0] <- 0.1
cor_r[cor_r==2.5] <- 0.2
cor_r[cor_r==5.0] <- 0.3

common_edge_color <- cor_r
common_edge_color[common_edge_color==0.1] <- "gray"
common_edge_color[common_edge_color==0.2] <- "orange"
common_edge_color[common_edge_color==0.3] <- "red"

common_edge_width <- cor_r
common_edge_width[cor_r==0.1] <- 1
common_edge_width[cor_r==0.2] <- 2.5
common_edge_width[cor_r==0.3] <- 5


E(common5)$weight <- 1
E(common5)$color <- common_edge_color
E(common5)$width <- common_edge_width




plot(common5,vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=11,vertex.size2=5, vertex.label.cex=0.7,
     layout=locs[as_ids(V(g1)) %in% as_ids(V(common5)),])
legend("bottomleft",inset = c(0.15, 0),
       bty = "n", title="Edges",
       legend=levels(as.factor(cor_r)),
       fill=  c("gray","orange","red"),
       horiz=F)
legend("bottomleft",inset = c(0, 0),bty = "n",
       legend=levels(color_facts),
       fill= color_clust,
       horiz=F,title="Nodes")


dev.off()


# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# pdf("result/Network/3.1.Network_Sup_Figures.pdf",width=75,height=15)
# plot(FigSup_1)
# plot(ggarrange(plotlist=list(FigSup_2,FigSup_3,FigSup_4),ncol=3))
# dev.off()

# pdf("result/Network/3.2.Network_Sup_Figures.pdf",width=20,height=10)
# FigSup_1
# FigSup_2
# FigSup_3
# FigSup_4
# dev.off()



saveRDS(list(FigSup_1,FigSup_2,FigSup_3,FigSup_4),output_file)


