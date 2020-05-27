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