#' Metrics Calculation Mit,RB,ERCC
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: Gene expression matrix.
#'
#' @return Gene expression Matrix
#' @examples metrics_calc(object)
metrics_calc <- function(object,remove_mt=TRUE,remove_rb=TRUE){
    # ============= ERCC Percentage ==============
    ERCC.index <- grep(pattern = "^ERCC", x = rownames(object), value = FALSE)
    percent.ERCC <- Matrix::colSums(object[ERCC.index, ])
    object <- object[-ERCC.index, ]
    # --------------------------------------------

    # ========= Mitochondria Percentage ==========
    Mit.index <- grep(pattern = "^MT-|^MT\\.", x = rownames(object), value = FALSE)
    percent.Mit <- Matrix::colSums(object[Mit.index, ])/Matrix::colSums(object)
    # --------------------------------------------


    if(remove_mt){
      object <- object[-Mit.index, ]
    }

    # ========= Ribosomal Percentage ==========
    RIB.index <- grep(pattern = "^RPL|^RPS", x = rownames(object), value = FALSE)
    percent.RIB <- Matrix::colSums(object[RIB.index, ])/Matrix::colSums(object)
    # -----------------------------------------

    if(remove_rb){
        object <- object[-RIB.index, ]
    }
    return(list("object"=object,"percent.ERCC"=percent.ERCC,"percent.mito"=percent.Mit,"percent.rb"=percent.RIB))
}




create_cds2 <-function (list_of_files, condition_names, min.features = 200, 
    min.cells = 5, remove_mt = TRUE, remove_rb = TRUE, outlier_detector = "MAD", 
    data_10x = FALSE, elbow = FALSE, imputation = FALSE, tool = "Seurat", 
    n_cores = 2, criteria_pass = 2, SCT = FALSE, vars.to.regress = NULL) 
{
    require(Seurat)
    require(monocle)
    tic()
    cat(cyan("Selected Parameteres :\n") %+% paste("\tRemove Mitochondrial Genes : ", 
        remove_mt, "\n\tRemove Ribosomal genes : ", remove_rb, 
        "\n\tImputation : ", imputation, "\n\t10xData : ", 
        data_10x, "\n\tElbow : ", elbow, "\n\tOutlier Detector : ", 
        outlier_detector))
    cat(cyan("\n\nLoading Data : ") %+% as.character(Sys.time()) %+% 
        "\n")
    tm_list <- list()
    init_num <- 1
    print("QC for every file")
    if (get_os() == "windows") {
        cat(yellow("Warning :") %+% "Your os is Windows. The Parallel file reading is implemented only for UNIX/MAC\n")
        n_cores = 1
    }
    tm_list <- mclapply(FUN = read_file2, c(1:length(list_of_files)), 
        list_of_files = list_of_files, condition_names = condition_names, 
        data_10x = data_10x, elbow = elbow, imputation = imputation, 
        outlier_detector = outlier_detector, min.features = min.features, 
        min.cells = min.cells, remove_mt = remove_mt, remove_rb = remove_rb, 
        mc.cores = n_cores, criteria_pass = criteria_pass, SCT = SCT, 
        vars.to.regress = vars.to.regress)
    if (length(list_of_files) > 1) {
        if (SCT == TRUE) {
            int.features <- SelectIntegrationFeatures(object.list = tm_list, 
                nfeatures = 3000)
            tm_list <- PrepSCTIntegration(object.list = tm_list, 
                anchor.features = int.features, verbose = FALSE)
            int.anchors <- FindIntegrationAnchors(object.list = tm_list, 
                normalization.method = "SCT", anchor.features = int.features, 
                verbose = FALSE)
            Seurat.combined <- IntegrateData(anchorset = int.anchors, 
                normalization.method = "SCT", verbose = FALSE)
        }
        else {
            Seurat.anchors <- FindIntegrationAnchors(object.list = tm_list, 
                dims = 1:20)
            Seurat.combined <- IntegrateData(anchorset = Seurat.anchors, 
                dims = 1:20)
        }
        DefaultAssay(object = Seurat.combined) <- "integrated"
        Seurat.combined$condition <- Idents(object = Seurat.combined)
    }
    else {
        Seurat.combined <- tm_list[[1]]
        Seurat.combined$condition <- Idents(object = Seurat.combined)
    }
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    Seurat.combined <- CellCycleScoring(Seurat.combined, s.features = s.genes, 
        g2m.features = g2m.genes)
    object <- Seurat.combined
    if (tolower(tool) == "monocle") {
        object <- seurat_to_monocle(object, pca = FALSE)
    }
    return(list(Combined = object, Data_List = tm_list))
}




read_file2<-function (iter_qc, list_of_files, condition_names, min.features = 200, 
    min.cells = 5, remove_mt = TRUE, remove_rb = TRUE, outlier_detector = "MAD", 
    data_10x = FALSE, elbow = FALSE, imputation = FALSE, tool = "Seurat", 
    criteria_pass = 2, SCT = FALSE, vars.to.regress = NULL) 
{
    init_num <- iter_qc
    file <- list_of_files[iter_qc]
    object <- load_files(file = file, data_10x = data_10x)
    if (elbow == TRUE) {
        Before_col_names <- colnames(object)
        object <- elbow_calc(object, condition_names, iter_qc)
        After_col_names_Elbow <- colnames(object)
        Rm_Cell_N_Elbow <- Before_col_names[!Before_col_names %in% 
            After_col_names_Elbow]
    }
    else {
        After_col_names_Elbow <- colnames(object)
        Rm_Cell_N_Elbow <- c()
    }
    Outlier_object <- Advanced_Outlier_detection(object, condition = condition_names[init_num], 
        criteria_pass = criteria_pass)
    object <- Outlier_object$object_filtered
    oultliers_index <- Outlier_object$oultliers_index
    After_col_names_MAD <- colnames(object)
    Rm_Cell_N_MAD <- After_col_names_Elbow[!After_col_names_Elbow %in% 
        After_col_names_MAD]
    Rm_Cell_N_Elbow <- c(Rm_Cell_N_Elbow, Rm_Cell_N_MAD)
    before_genes <- rownames(object)
    metrics_output <- metrics_calc(object = object, remove_mt = remove_mt, 
        remove_rb = remove_rb)
    object <- metrics_output$object
    percent.mito <- metrics_output$percent.mito
    percent.rb <- metrics_output$percent.rb
    Seurat <- CreateSeuratObject(counts = object, project = condition_names[init_num], 
        min.cells = min.cells, min.features = min.features, meta.data = data.frame(percent.mito = percent.mito, 
            percent.rb = percent.rb))
    Seurat$stim <- condition_names[init_num]
    After_col_names_SEURAT <- colnames(Seurat)
    Rm_Cell_N_SEURAT <- After_col_names_MAD[!After_col_names_MAD %in% 
        After_col_names_SEURAT]
    Rm_Cell_N_Elbow <- c(Rm_Cell_N_Elbow, Rm_Cell_N_SEURAT)
    write.table(Rm_Cell_N_Elbow, paste0(condition_names[init_num], 
        "Cells_Removed.tsv"), sep = "\t")
    after_genes <- rownames(Seurat@assays$RNA@counts)
    rm_genes <- before_genes[!before_genes %in% after_genes]
    write.table(rm_genes, paste0(condition_names[init_num], "Genes_Removed.tsv"), 
        sep = "\t")
    if (imputation == TRUE) {
        cat(green("\nImputation with SAVER\n"))
        allcells <- as.matrix(Seurat@assays$RNA@counts)
        library(SAVER)
        allcells[is.na(allcells)] <- 0
        cortex.saver <- saver(allcells, ncores = n_cores)
        allcells <- as.matrix(cortex.saver$estimate)
        allcells[is.na(allcells)] <- 0
        Seurat_imputed <- CreateSeuratObject(counts = allcells, 
            project = condition_names[init_num], min.cells = min.cells, 
            min.features = min.features, meta.data = data.frame(percent.mito = Seurat$percent.mito, 
                percent.rb = Seurat$percent.rb))
        Seurat_imputed$stim <- condition_names[init_num]
        Seurat <- Seurat_imputed
    }
    if (SCT == TRUE) {
        Seurat <- SCTransform(object = Seurat, verbose = FALSE, 
            vars.to.regress = vars.to.regress)
    }
    else {
        Seurat <- NormalizeData(object = Seurat, normalization.method = "LogNormalize", 
            scale.factor = 10000)
        Seurat <- FindVariableFeatures(object = Seurat, selection.method = "vst", 
            nfeatures = 5000)
        top10 <- head(x = VariableFeatures(object = Seurat), 
            10)
        plot1 <- VariableFeaturePlot(object = Seurat)
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        pdf(paste(Sys.Date(), "_top10_var_feat_", Seurat$stim, 
            ".pdf"))
        print(plot2)
        dev.off()
    }
    tm1 <- Seurat
    print(dim(tm1))
    return(tm1)
}


load_files<-function (file, data_10x) {
    if (data_10x) {
        print(file)
        barcode.path <- paste0(file, "barcodes.tsv")
        features.path <- paste0(file, "features.tsv")
        matrix.path <- paste0(file, "matrix.mtx")
        mat <- as.matrix(readMM(file = matrix.path))
        feature.names = read.delim(features.path, header = FALSE, 
            stringsAsFactors = FALSE)
        barcode.names = read.delim(barcode.path, header = FALSE, 
            stringsAsFactors = FALSE)
        colnames(mat) = barcode.names$V1
        rownames(mat) = feature.names$V2
        object <- mat[, order(colSums(as.matrix(mat)), decreasing = T)]
    }
    else {
        object <- read.csv(file, header = T, row.names = 1, sep = "\t")
        object <- object[, order(colSums(object), decreasing = T)]
    }
    rows_nms <- str_replace_all(toupper(rownames(object)), "-", 
        ".")
    rownames(object) <- rows_nms
    return(object)
}




Advanced_Outlier_detection <- function (object, filtering = TRUE, min.features = 200, min.cells = 5, 
    condition = "", criteria_pass = 2) 
{
    cat(green("Outlier Detection Started (At least 2 criteria)"))
    nFeatures <- colSums(object != 0)
    object <- object[rowSums(object != 0) >= min.cells, nFeatures >= 
        min.features]
    nCounts <- colSums(object)
    nFeatures <- colSums(object != 0)
    percent.mit <- colSums(object[grep("^MT-|^MT\\.", rownames(object)), 
        ])/nCounts
    percent.rb <- colSums(object[grep("^RPL|^RPS", rownames(object)), 
        ])/nCounts
    df_qc <- data.frame(Cond = rep(condition, dim(object)[2]), 
        nFeatures = nFeatures, Total_MRNA = nCounts, Percent_Mit = percent.mit)
    if (filtering) {
        print("MAD outlier Detection")
        false_list <- rep(TRUE, length(nFeatures))
        res_out_nFeat <- outliers_mad(nFeatures, threshold = 2)
        out_nFeat <- false_list
        out_nFeat[res_out_nFeat$outliers_pos] <- FALSE
        res_out_nCounts <- outliers_mad(nCounts, threshold = 2)
        out_nCounts <- res_out_nCounts$outliers_pos
        out_nCount <- false_list
        out_nCount[out_nCounts] <- FALSE
        if (data_10x) {
            doublerate <- round(dim(object)[2] * 0.9/1000)
        }
        res_out_percent.mit <- outliers_mad(percent.mit, threshold = 1.5)
        out_percent.mit <- percent.mit[res_out_percent.mit$outliers_pos] > 
            mean(percent.mit)
        out_percent.mits <- res_out_percent.mit$outliers_pos[out_percent.mit]
        out_percent.mit <- false_list
        out_percent.mit[out_percent.mits] <- FALSE
        QC_Matrix <- data.frame(nCounts = out_nCount, nFeatures = out_nFeat, 
            percent.mit = out_percent.mit)
        outlier_index <- rep(FALSE, length(nFeatures))
        QC_Matrix$Passed_Tests <- rowSums(QC_Matrix)
        outlier_index[QC_Matrix$Passed_Tests >= criteria_pass] <- TRUE
        df_filtered <- object[, outlier_index]
        if (res_out_nCounts$limit[1] < 0) {
            res_out_nCounts$limit[1] <- 0
        }
        res_out_nCounts <- res_out_nCounts$limit
        res_out_nFeat <- c(min.features, res_out_nFeat$limit[2])
        res_out_percent.mit <- res_out_percent.mit$limit[2]
    }
    else {
        df_filtered <- object
        oultliers_index <- NULL
        return(list(object_filtered = df_filtered, oultliers_index = oultliers_index))
    }
    pdf(paste("result/Preprocess/",sample,"/",Sys.Date(), "_Outliers_hist_",condition_name,".pdf", sep = ""), width = 7, height = 3)
    par(mfrow = c(1, 3))
    hist(nCounts, col = "lightblue", main = "Histogram for nCounts")
    abline(v = res_out_nCounts, col = "red", lwd = 3, lty = 2)
    hist(nFeatures, col = "lightblue", main = "Histogram for nFeatures")
    abline(v = res_out_nFeat, col = "red", lwd = 3, lty = 2)
    hist(percent.mit, col = "lightblue", main = "Histogram for percent.mit")
    abline(v = res_out_percent.mit, col = "red", lwd = 3, 
        lty = 2)
    dev.off()
    print("Dimensions after Filtering out Low Quality Cells")
    print(dim(df_filtered))
    cat(green("Outlier Detection Finished\n"))
    Vln_QC(df_qc, condition = condition, title = "before", 
        outliers = outlier_index)
    nCounts <- colSums(df_filtered)
    nFeatures <- colSums(df_filtered != 0)
    percent.mit <- colSums(df_filtered[grep("^MT-|^MT\\.", 
        rownames(df_filtered)), ])/nCounts
    percent.rb <- colSums(df_filtered[grep("^RPL|^RPS", 
        rownames(df_filtered)), ])/nCounts
    df_qc <- data.frame(Cond = rep(condition, dim(df_filtered)[2]), 
        nFeatures = nFeatures, Total_MRNA = nCounts, Percent_Mit = percent.mit)
    Vln_QC(df_qc, condition = condition, title = "after")
    return(list(object_filtered = df_filtered, oultliers_index = QC_Matrix))
}







reduce_dim <-function (object, assay = "RNA", project, resolution = FALSE) {
    tool <- object_identifier(object)
    if (tolower(tool) == "monocle") {
        object <- detectGenes(object, min_expr = 0.1)
        disp_table <- dispersionTable(object)
        unsup_clustering_genes <- subset(disp_table, mean_expression >= 
            0.1)
        object <- setOrderingFilter(object, unsup_clustering_genes$gene_id)
    }
    else {
        print(paste0("The projection will be on ", assay))
        if (assay == "RNA") {
            DefaultAssay(object) <- "RNA"
            object <- NormalizeData(object)
            object <- FindVariableFeatures(object = object, selection.method = "vst", 
                nfeatures = 5000)
            all.genes <- rownames(object)
            object <- ScaleData(object, features = all.genes)
        }
        else if (assay == "integrated") {
            DefaultAssay(object) <- "integrated"
            all.genes <- rownames(object)
            object <- ScaleData(object, features = all.genes)
        }
        else {
            DefaultAssay(object) <- "integrated"
        }
        object <- RunPCA(object = object, verbose = FALSE)
    }
    pca_elbow(object, project)
    num_dim <- calc_num_pc(object = object, 0.95)
    if (tolower(tool) == "monocle") {
        object <- reduceDimension(object, max_components = 2, 
            num_dim = num_dim, reduction_method = "tSNE", 
            residualModelFormulaStr = "~num_genes_expressed", 
            verbose = T, check_duplicates = TRUE)
    }
    else {
        object <- RunUMAP(object = object, reduction = "pca", 
            dims = 1:num_dim)
        object <- RunTSNE(object = object, reduction = "pca", 
            dims = 1:num_dim)
        object <- FindNeighbors(object = object, reduction = "pca", 
            dims = 1:num_dim)
    }
    optimal_output <- optimal_clusters(object, k.max = 10, save = TRUE, 
        resolution = resolution)
    object <- optimal_output$object
    opt_num <- optimal_output$opt_num
    sil_scor <- optimal_output$sil_scor
    return(list(Combined = object, Num_dim = num_dim))
}








optimal_clusters<-function (object, k.max = 10, save = TRUE, resolution = FALSE) {
    colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
    color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")[c(5,1,2,3,4,9,6,7,8)]
    color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
    color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
    color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
    color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

    mean_score <- c()
    Clusts_num <- c()
    sd_score <- c()
    
    if (!resolution) {
        resolution <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
    }
    for (i in c(1:length(resolution))) {
        res <- resolution[i]
        object <- FindClusters(object, resolution = res)
        object$Cluster <- Idents(object = object)
        projection_mat <- Embeddings(object = object$pca)
        sil_scor <- silhouette(as.numeric(object$Cluster), 
            dist(projection_mat))
        mean_score[i] = summary(sil_scor)$si.summary[4]
        wt_size <- summary(sil_scor)$clus.sizes/sum(summary(sil_scor)$clus.sizes)
        xm <- mean_score[i]
        v <- sum(wt_size * (summary(sil_scor)$clus.avg.widths - 
            xm)^2)/sum(wt_size)
        sd_score[i] = sqrt(v)
        Clusts_num[i] = length(levels(object$Cluster))
    }
    silhouette_df <- data.frame(resolution = resolution, 
        Mean_Silhouette = mean_score)
    opt_res <- resolution[which.max(mean_score)]
    p1 <- ggplot(data = silhouette_df, aes(x = resolution, 
        y = Mean_Silhouette)) + geom_point() + geom_line() + 
        scale_x_continuous() + geom_vline(xintercept = opt_res, 
        linetype = "dashed", color = "red") + 
        ggtitle("Optimal number of clusters") + ylab("Mean Silhouette score") + 
        xlab("Resolution")
    object <- FindClusters(object, resolution = opt_res)
    object$Cluster <- as.factor(as.numeric(object$seurat_clusters))
    names(object$Cluster) <- colnames(object)
    projection_mat <- Embeddings(object = object$umap)
    sil_scor <- silhouette(as.numeric(object$Cluster), dist(projection_mat))
    print_text <- paste("The optimal Resolution based on Mean Silhouette score is ", 
        opt_res)
    print(paste(rep("#", nchar(print_text)), collapse = ""))
    print(print_text)
    print(paste(rep("-", nchar(print_text)), collapse = ""))
    optimal <- opt_res
    
    data <- data.frame(cbind(factor(object$Cluster), paste(" ", 
        object$condition)))
    colnames(data) <- c("Cluster", "Condition")
    pdf(paste("result/Mapping/",Sys.Date(), "_Barplot_num-Cond_per_Cluster.pdf",sep = ""))
    print(ggplot(data, aes(Cluster)) + geom_bar(aes(fill = Condition)) + 
        guides(col = guide_legend(ncol = 2, )) + scale_fill_brewer(palette = "Dark2") + 
        theme(legend.position = "bottom"))
    dev.off()
    num_clusters <- length(unique(object$Cluster))
    if (save) {
        pdf(paste("result/Mapping/",Sys.Date(),"_pd_optimal_cluster.pdf",sep = ""))
        print(p1)
        dev.off()
        pdf(paste("result/Mapping/",Sys.Date(), "_pd_silhouette.pdf", sep = ""))
        plot(sil_scor, col = color_list$Cluster[1:num_clusters])
        dev.off()
    }
    else {
        print(p1)
    }
    return(list(object = object, opt_num = optimal, sil_scor = sil_scor))
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
    pdf(paste("result/Preprocess/",sample,"/",Sys.Date(), "_QC_", title, ".pdf", 
        sep = ""))
    grid.arrange(p1, p2, p3, nrow = 1)
    dev.off()
}



#' Elbow Calculation
#' @author Dimitrios Kyriakis
#'
#' @param object: Gene expression matrix.
#' @param sample: Condition names.
#' @param iter_qc: Iteration.
#'
#' @return Gene expression Matrix
#' @examples elbow_calc(object,condition_names,iter_qc)
elbow_calc <- function(object,sample,iter_qc){
    require(ecp)
    #====================== ELBOW ==========================================
    pdf(paste("result/Preprocess/",sample,"/",Sys.Date(),"_QC_Kneeplot_Zeros.pdf",sep=""))
    # Knee plot
    knee_data <- object
    # plot(colSums(is.na(knee_data)),main="Number of missing genes per cell",ylab="Number of missing genes",xlab="Cells")
    knee_data[is.na(knee_data)] <- 0
    plot(colSums(knee_data == 0),main=sample,ylab="Number of non expressed genes",xlab="Cells")
    # print("No memory error")
    # ============== Elbow Analysis ============= #

    y<- cumsum(sort(colSums(knee_data),TRUE))
    x<-c(1:length(y))
    y1<- as.matrix(y)
    y2<- as.matrix(y[500:length(y)])

    real_index <- e.divisive(diff(y1),k=1,min.size=2)
    real_index <- real_index$considered.last
    print(paste("The real elbow value is ",real_index))

    # pdf("knee.pdf")
    plot(x, y, pch=19,xlab="Cells",ylab="Cumulative Total mRNA")
    if (real_index<500){
        y2<- as.matrix(y[real_index:length(y)])
        indicies <- e.divisive(diff(y2),k=1,min.size=2)$considered.last +500
        print(paste("The elbow value adjusted to be between 500 and 1000. The new value is ",indicies))
        points(x[indicies], y[indicies], pch=19, col='red',main=paste0("Elbow ",real_index))
    }else{
        indicies <- real_index
    }

    points(x[real_index], y[real_index], pch=19, col='lightblue')
    abline(v=500,col="red")
    abline(v=1000,col="red")
    # dev.off()
    if(indicies>1000& real_index<1000){
      object<- object[,1:700]
    }else{
      object<- object[,1:indicies]
    }

    # ---------------------------------------------
    dev.off()
    return(object)
}







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

