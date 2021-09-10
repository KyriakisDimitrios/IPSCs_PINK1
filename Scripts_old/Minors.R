get_os <- function () {
    sysinf <- Sys.info()
    if (!is.null(sysinf)) {
        os <- sysinf["sysname"]
        if (os == "Darwin") 
            os <- "osx"
    }
    else {
        os <- .Platform$OS.type
        if (grepl("^darwin", R.version$os)) 
            os <- "osx"
        if (grepl("linux-gnu", R.version$os)) 
            os <- "linux"
    }
    tolower(os)
}


Vln_QC <- function (df_qc, condition, title, outliers = NULL) {
    library(gridExtra)
    if (is.null(outliers)) {
        outliers <- 1
    }
    else {
        outliers <- as.numeric(!outliers) + 1
    }
    p1 <- ggplot(df_qc, aes(factor(Cond), nFeatures, fill = Cond)) + 
        geom_violin(width = 0.8, show.legend = FALSE) + ggtitle("nFeatures") + 
        xlab("") + ylab("Expression Level") + geom_jitter(width = 0.3, 
        size = 1, show.legend = FALSE, colour = outliers)
    p2 <- ggplot(df_qc, aes(factor(Cond), Total_MRNA, fill = Cond), 
        width = 1) + geom_violin(width = 0.8, scale = "count", 
        show.legend = FALSE, adjust = 1/2) + ggtitle("nCounts") + 
        xlab("") + ylab("Expression Level") + geom_jitter(width = 0.3, 
        size = 1, show.legend = FALSE, colour = outliers)
    p3 <- ggplot(df_qc, aes(factor(Cond), Percent_Mit, fill = Cond), 
        width = 1) + geom_violin(width = 0.8, show.legend = FALSE) + 
        ggtitle("Percent Mit") + xlab("") + ylab("Expression Level") + 
        geom_jitter(width = 0.3, size = 1, show.legend = FALSE, 
            colour = outliers)
    pdf(paste(Sys.Date(), "QC", condition, title, ".pdf", 
        sep = "_"))
    grid.arrange(p1, p2, p3, nrow = 1)
    dev.off()
}


object_identifier<- function (object) {
    if (class(object) == "Seurat") {
        tool = "seurat"
    }
    else {
        tool = "monocle"
    }
    return(tool)
}



scatter_gene<-function (object, features, assay = "RNA", reduction = "umap", 
    size = 2, ncol = 2, nrow = 2) {
    DefaultAssay(object) <- "RNA"
    proj <- object@reductions[[reduction]]@cell.embeddings
    exp <- as.matrix(object@assays[[assay]]@data)
    df <- cbind(proj, object@meta.data, t(exp))
    library(dplyr)
    myPalette <- colorRampPalette(brewer.pal(9, "OrRd"))
    myPalette <- colorRampPalette(c("lightgrey", "#FDBB84", 
        "#EF6548", "#D7301F", "#B30000", "#7F0000"))
    sc <- scale_colour_gradientn(colours = myPalette(9))
    results <- lapply(features, function(x) {
        ggplot(df %>% arrange(get(x)), aes(x = get(colnames(df)[1]), 
            y = get(colnames(df)[2]), color = get(x))) + geom_point(size = size) + 
            ggtitle(x) + sc + theme_cowplot() + theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5), legend.title = element_blank()) + 
            xlab(colnames(df)[1]) + ylab(colnames(df)[2])
    })
    p1 <- ggarrange(plotlist = results, ncol = ncol, nrow = nrow)
    return(p1)
}



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
            corner(return_cor)

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
                    G1 <- ICSWrapper::cor_net(dataset=dataset,regulated=regulated,cor_thres=threshold,corr.method=corr.method)
                }else if(method=="genie3"){
                    G1 <- ICSWrapper::genie3_net(dataset=dataset,regulated=regulated,weight.threshold=threshold,nCores=nCores,nTrees=nTrees)
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









