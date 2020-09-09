
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
library("qgraph")
library("igraph")

## ----Venn------------------------------------------------------------------------------------
library(VennDiagram)
library(EnhancedVolcano)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
olor_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/5.Correlation_Net/")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/5.Correlation_Net/")


Combined <- readRDS("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess/IPSCs_Combined.rds")



# ========================================= FUNCTIONS ====================================================
# ========================================= FUNCTIONS ====================================================
# ========================================= FUNCTIONS ====================================================
#' Build Correlation Network
#'
#' @author Dimitrios Kyriakis
#' @export
#' @param dataset : A Matrix with columns the nodes
#' @param regulated : A vector with "UP Regulated" and "DOWN Regulated".
#' @param cor_thres : Correlation Threshold
#' @param corr.method : What kind of correlations should be computed? Default is "pearson", but "spearman" and "kendall" are also supported.
#' @return Network
#' @examples cor_net(heat_matrix)
#'
cor_net <- function (dataset,regulated=NULL,cor_thres=0.5,corr.method="pearson"){
    require(tidyverse)
    require(corrr)
    require(igraph)
    require(ggraph)

    net_mat_plot <- as.matrix(dataset)

    # # Create a tidy data frame of correlations
    # d <- as.data.frame(t(net_mat_plot))
    # tidy_cors <- d %>%
    #     correlate(method=corr.method) %>%
    #     stretch()


    tryCatch(
        {
            # Convert correlations stronger than some value
            # to an undirected graph object
            library("Hmisc")
            res2 <- rcorr(t(as.matrix(dataset_r)))
            flattenCorrMatrix <- function(cormat, pmat) {
                ut <- upper.tri(cormat)
                data.frame(
                    row = rownames(cormat)[row(cormat)[ut]],
                    column = rownames(cormat)[col(cormat)[ut]],
                    cor  =(cormat)[ut],
                    p = pmat[ut]
                    )
            }
            return_cor <-flattenCorrMatrix(res2$r, res2$P)
            useful::corner(return_cor)

            graph_cors <- return_cor %>%
                dplyr::filter(abs(cor) > cor_thres)%>%
                dplyr::filter(abs(p) < 0.05) %>%
                graph_from_data_frame(directed = FALSE)
                
            # graph_cors <- tidy_cors %>%
            #     dplyr::filter(abs(r) > cor_thres) %>%
            #     graph_from_data_frame(directed = FALSE)

            G2 <-graph_cors
            return(G2)
        },
        error=function(cond){
            return("No Edge Found. Try to decrease MMHC threshold")
        }
    )

}









#' Build Regulatory Network Using Genie3
#'
#' @author Dimitrios Kyriakis
#' @export
#' @param dataset : A Matrix with columns the nodes
#' @param regulated : A vector with "UP Regulated" and "DOWN Regulated".
#' @param weight.threshold : Weight Threshold (less adds more edges)
#' @param nCores : Number of cores to use for parallel computing. Default: 4.
#' @param nTrees : Number of trees in an ensemble for each target gene. Default: 2000.
#' @return Network
#' @examples genie3_net(dataset,regulated=NULL,weight.threshold=0.05,title="",nCores=4,nTrees=2000)
#'
genie3_net <- function (dataset,regulated=NULL,weight.threshold=0.05,title="",nCores=4,nTrees=2000){
    require(GENIE3)
    require(igraph)
    require(RCy3)
    require(Rgraphviz)
    net_mat_plot <- as.matrix(dataset)
    tryCatch(
        {

        weightMat <- GENIE3(net_mat_plot, nCores=nCores, verbose=TRUE,nTrees=nTrees)
        link.list <- getLinkList(weightMat,threshold=weight.threshold)
        edge_listsi <- link.list[!duplicated(link.list),]

        Gsi <- graph.data.frame(edge_listsi,directed = FALSE)
        Asi <- get.adjacency(Gsi,sparse = F,attr = "weight",type = "both")
        G2 <- graph.adjacency(Asi,mode = "undirected",weighted = T)

        return(G2)
    },
    error=function(cond){
        return("Debugonce()")
    }
    )
}





#' Build Network
#'
#' @author Dimitrios Kyriakis
#' @export
#' @param dataset : A Matrix with columns the nodes
#' @param regulated : A vector with "UP Regulated" and "DOWN Regulated".
#' @param threshold : Threshold (aplha parameter for mmhc/h2pc, pvalue for corellation, weeight for genie3)
#' @param nCores : Parallel (aplied only for genie3)
#' @param nTrees : Number of trees for Random forest (aplied only for genie3)
#' @param corr.method : What kind of correlations should be computed? Default is "pearson", but "spearman" and "kendall" are also supported.
#' @param title : Name of the pdf
#' @return Network
#' @examples ics_net(dataset,regulated=NULL,threshold=0.01,nCores=4,nTrees=2000,title="Day09",method="mmhc",corr.method="pearson")
#'
ics_net <- function (dataset,regulated=NULL,threshold=0.05,nCores=4,nTrees=2000,title="",method="mmhc",corr.method="pearson"){
    require(bnlearn)
    require(bnviewer)
    require(igraph)

    labelx <- paste0(Sys.Date(),"_",title)
    graphics.off()
    tryCatch(
        {
            # ===== MMPC ======================
            if(method%in%c("pc","aracne","mmhc")){
                net_mat_plot <- as.data.frame(t(as.matrix(dataset)))
                if(method=="pc"){
                    temp_network <- bnlearn::h2pc(net_mat_plot,restrict.args=list("alpha"=threshold))
                }else if (method=="aracne"){
                    temp_network <- bnlearn::aracne(net_mat_plot)
                }else{
                    temp_network <- bnlearn::mmhc(net_mat_plot,restrict.args=list("alpha"=threshold))
                }
                G1 <- bnviewer::bn.to.igraph(temp_network)
            }else{
                if(method=="correlation"){
                    G1 <- cor_net(dataset=dataset,regulated=regulated,cor_thres=threshold,corr.method=corr.method)
                }else if(method=="genie3"){
                    G1 <-genie3_net(dataset=dataset,regulated=regulated,weight.threshold=threshold,nCores=nCores,nTrees=nTrees)
                }else{
                    print("Not available method")
                    return()
                }
            }



            G2 <- igraph::simplify(G1)
            if(method=="aracne"){
                G2 <- as.undirected(G2)

            }
            if(!method%in%c("pc","mmhc")){
                clusterlouvain <- cluster_louvain(G2)
            }else{
                G3 <- as.undirected(G2)
                clusterlouvain <- cluster_louvain(G3)
            }

            V(G2)$cex <- 0.11
            V(G2)$label.cex <-0.3
            G3 <- G2
            Degree_net <- igraph::degree(G2, mode = "in")
            Between_net <- igraph::betweenness(G2)

            pdf(paste0(labelx,"_Net_",method,".pdf"))
            # =================================== UP AND DOWN REGULATION ==========================================
            if(!is.null(regulated)){
                regulated <- regulated[names(regulated)%in%vertex_attr(G2)$name]
                Color_rest <-c("skyblue","pink")[as.factor(regulated)]
                set.seed(1234)  # set seed to make the layout reproducible
                V(G2)$color <- Color_rest
                plot(G2, edge.arrow.size=0.25, main="Regulated",cex=0.1,main=labelx)
                legend("bottomleft",bty = "n",
                       legend=levels(as.factor(regulated)),
                       fill= c("skyblue","pink"), horiz=F)
            }else{
                Color_rest <-c("skyblue")
            }
            set.seed(1234)  # set seed to make the layout reproducible
            V(G2)$color <- Color_rest
            
            set.seed(1234)
            plot(G2, edge.arrow.size=0.25, main="Gene Network",cex=0.1,main=labelx)
            set.seed(1234)
            plot(clusterlouvain,G2, edge.arrow.size=0.25, main="Gene Network",cex=0.1,main=labelx )
            # -----------------------------------------------------------------------------------------------------


            # =================================== Betweenness CENTRALITY ==========================================
            G4 <- G2
            V(G4)$color <-Color_rest
            #Edge Options: Color
            E(G4)$color <- "grey"
            set.seed(1234)  # set seed to make the layout reproducible
            print(Between_net)
            if(length(unique(Between_net))==1){
                plot(G4,cex=0.1, edge.arrow.size=0.25,main="Betweeness Centrality")
            }else{
                V(G4)$size=ICSWrapper::rescale_ics(Between_net) #because we have wide range, I am dividing by 5 to keep the high in-degree nodes from overshadowing everything else.
                plot(G4, edge.arrow.size=0.25,main="Betweeness Centrality")

            }

            # -----------------------------------------------------------------------------------------------------




            ## ====================================== DEGREE ======================================================
            G5 <- G2
            #Node or Vetex Options: Size and Color
            V(G5)$color <- Color_rest
            #Edge Options: Color
            E(G5)$color <- "grey"
            #Plotting, Now Specifying an arrow size and getting rid of arrow heads
            #We are letting the color and the size of the node indicate the directed nature of the graph
            set.seed(1234)  # set seed to make the layout reproducible
            print(Degree_net)
            if(length(unique(Degree_net))==1){
                plot(G5,cex=0.1, edge.arrow.size=0.25,main="Degree In")
            }else{
                V(G5)$size=ICSWrapper::rescale_ics(Degree_net) #because we have wide range, I am dividing by 5 to keep the high in-degree nodes from overshadowing everything else.
                plot(G5, edge.arrow.size=0.25,main="Degree In")#, vertex.label = NA
            }
            dev.off()
            # -----------------------------------------------------------------------------------------------------
            graphics.off()

            # ====================================== DATA FRAME ===================================================
            if(!method%in%c("pc","mmhc")){
                if(!is.null(regulated)){
                    df_G2 <- data.frame(Gene=vertex_attr(G2)$name,Cluster=clusterlouvain$membership,Degree=Degree_net,BTC=Between_net,Regulated=regulated,Color=vertex_attr(G2)$color)
                }else{
                    df_G2 <- data.frame(Gene=vertex_attr(G2)$name,Cluster=clusterlouvain$membership,Degree=Degree_net,BTC=Between_net,Color=vertex_attr(G2)$color)
                }
            }else{
                if(!is.null(regulated)){
                    df_G2 <- data.frame(Gene=vertex_attr(G2)$name,Degree=Degree_net,BTC=Between_net,Regulated=regulated,Color=vertex_attr(G2)$color)
                }else{
                    df_G2 <- data.frame(Gene=vertex_attr(G2)$name,Degree=Degree_net,BTC=Between_net,Color=vertex_attr(G2)$color)
                }
            }
            write.table( df_G2,paste0(labelx,"_",method,".tsv"), row.names=FALSE,sep="\t")
            graphics.off()
            # -----------------------------------------------------------------------------------------------------
            return(list("graph"=G3,"df"=df_G2))
    },
    error=function(cond){
        return("Try exception gave an error. Check the treshold. Or use debugonce(ics_net) to check where the error comes from and report it")
    }
    )

}
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------










# ========================================= DF NETWORK ====================================================
# ========================================= DF NETWORK ====================================================
# ========================================= DF NETWORK ====================================================

DefaultAssay(Combined) <- "RNA"
Combined<- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))

Control<-subset(Combined,subset = Treatment=="Control")
PINK1<-subset(Combined,Treatment=="PINK") 


dataset <- as.data.frame(Combined@assays$RNA@data)


net_analysis="extented"
net_analysis="net"
# DF
graph_annotation <- read.csv("NODES_and_pathways_11.6.20_9.csv")
graph_annotation$group1 <- as.vector(graph_annotation$group1)
graph_annotation$group1[graph_annotation$group1==""] <- "NA"
first_graph <- read.csv("RUN_12_EDGES_ST_GM_REMERGE_used_11_6_20.csv")




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
match(missing,annot_29$Id)
match(missing,annot_29$Label)


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
pdf("Network_DE_only.pdf",width=12,height=10)
set.seed(123)
locs <- layout_on_sphere(g1)*0.4

e <- get.edgelist(g1,names=FALSE)
locs <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g1),
                                          area=50*(vcount(g1)^3),repulse.rad=(vcount(g1)^4.1))



col_added<- unlist(lapply(color_facts,FUN = function(x){ if(x %in% c("A","B","C","B-PD")){return(1)}else{return(2)}}))    

V(g1)$color <- Color_rest

plot(g1,  vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=11,vertex.size2=5, vertex.label.cex=0.7,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F)


bet<-betweenness(g1)
n <- vcount(g1)
test_bet <- (2 * bet) / (n*n - 3*n + 2) * 500

test_bet<-scale(bet,center = 0)+4

plot(g1, vertex.label.color="black",vertex.label.font=11,
     vertex.size2=5, vertex.label.cex=0.7,vertex.size=test_bet,
     layout=locs, 
     main="Original")
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




# ========================================= DF+ADDED NETWORK ====================================================
# ========================================= DF+ADDED NETWORK ====================================================
# ========================================= DF+ADDED NETWORK ====================================================
library("qgraph")
library("igraph")
library("Seurat")
setwd("C:\\Users\\dimitrios.kyriakis\\Desktop\\PhD\\Projects\\Michi_Data\\2020-04-17_seurat_elbow_TRUE_Mito-FALSE_Ribo-FALSE_SCT-TRUE_criteria_pass-3\\Gabriela_Networks\\Final_Network/")

DefaultAssay(Combined) <- "RNA"
Combined<- NormalizeData(Combined)
Combined <- ScaleData(Combined,rownames(Combined@assays$RNA@counts))

Control<-subset(Combined,subset = Treatment=="Control")
PINK1<-subset(Combined,Treatment=="PINK") 


dataset <- as.data.frame(Combined@assays$RNA@data)


net_analysis="extented"

# Ubiq
graph_annotation <- read.csv("NODES_May_29_Manual_Mito_Ubiq.csv")
# Net
first_graph <- read.csv("EDGES_part1_manual_May 29.csv")
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
match(missing,annot_29$Id)
match(missing,annot_29$Label)


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

color_clust<-sns.palplot(sns.color_palette("colorblind", 8))
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

pdf("Network_extented.pdf",width=12,height=10)
set.seed(12345)
E(g1)$weight = 0.15
locs <- layout_with_fr(g1,dim=2, niter=1000)

col_added<- unlist(lapply(color_facts,FUN = function(x){ if(x %in% c("DE")){return(1)}else{return(2)}}))



col_added[col_added==1] <- "Differentially Expressed"
col_added[col_added==2] <- "Added"
col_added<-as.factor(col_added)
V(g1)$color <- col_added


plot(g1, vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=11,vertex.size2=5, vertex.label.cex=0.7,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(col_added),fil=categorical_pal(8)[c(2,1)], horiz=F)

V(g1)$color <- Color_rest

plot(g1,  vertex.label.color="black",vertex.shape="rectangle",vertex.label.font=11,
     vertex.size=11,vertex.size2=5, vertex.label.cex=0.7,
     layout=locs)
legend("bottomleft",bty = "n",
       legend=levels(color_facts),
       fill= color_clust, horiz=F)


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
