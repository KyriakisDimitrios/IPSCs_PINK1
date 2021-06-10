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
    pdf(paste(Sys.Date(), "_Outliers_hist_", condition, 
        ".pdf", sep = ""), width = 7, height = 3)
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
    tool <- object_identifier(object)
    mean_score <- c()
    Clusts_num <- c()
    sd_score <- c()
    if (tolower(tool) == "monocle") {
        dist_mat <- dist(t(as.matrix(object@reducedDimA)))
        for (i in c(2:k.max)) {
            object <- clusterCells(object, num_clusters = i)
            if (length(levels(object$Cluster)) < i) {
                object <- clusterCells(object, num_clusters = i + 
                  1, verbose = F)
            }
            sil_scor <- silhouette(as.numeric(object$Cluster), 
                dist_mat)
            mean_score[i - 1] = summary(sil_scor)$si.summary[4]
            wt_size <- summary(sil_scor)$clus.sizes/sum(summary(sil_scor)$clus.sizes)
            xm <- mean_score[i - 1]
            v <- sum(wt_size * (summary(sil_scor)$clus.avg.widths - 
                xm)^2)/sum(wt_size)
            sd_score[i - 1] = sqrt(v)
            Clusts_num[i - 1] = length(levels(object$Cluster))
        }
        silhouette_df <- data.frame(Clusters = Clusts_num, Mean_Silhouette = mean_score)
        num_clusters <- which.max(mean_score) + 1
        p1 <- ggplot(data = silhouette_df, aes(x = Clusters, 
            y = Mean_Silhouette)) + geom_point() + geom_line() + 
            scale_x_continuous(breaks = seq(2, 15, 1)) + geom_vline(xintercept = num_clusters, 
            linetype = "dashed", color = "red") + 
            ggtitle("Optimal number of clusters") + ylab("Mean Silhouette score") + 
            xlab("Number of Clusters")
        cat(cyan("Clustering : ") %+% paste("The optimal number of clusters based on Mean Silhouette score is ", 
            num_clusters, ".\n"))
        object <- clusterCells(object, num_clusters = num_clusters)
        if (length(levels(object$Cluster)) != num_clusters) {
            object <- clusterCells(object, num_clusters = num_clusters + 
                1)
        }
        sil_scor <- silhouette(as.numeric(object$Cluster), dist_mat)
        optimal <- num_clusters
    }
    else {
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
    }
    data <- data.frame(cbind(factor(object$Cluster), paste(" ", 
        object$condition)))
    colnames(data) <- c("Cluster", "Condition")
    pdf(paste(Sys.Date(), "Barplot_num-Cond_per_Cluster.pdf", 
        sep = "_"))
    print(ggplot(data, aes(Cluster)) + geom_bar(aes(fill = Condition)) + 
        guides(col = guide_legend(ncol = 2, )) + scale_fill_brewer(palette = "Dark2") + 
        theme(legend.position = "bottom"))
    dev.off()
    num_clusters <- length(unique(object$Cluster))
    if (save) {
        pdf(paste(Sys.Date(), "pd_optimal_cluster.pdf", 
            sep = "_"))
        print(p1)
        dev.off()
        pdf(paste(Sys.Date(), "pd_silhouette.pdf", sep = "_"))
        plot(sil_scor, col = color_list$Cluster[1:num_clusters])
        dev.off()
    }
    else {
        print(p1)
    }
    return(list(object = object, opt_num = optimal, sil_scor = sil_scor))
}