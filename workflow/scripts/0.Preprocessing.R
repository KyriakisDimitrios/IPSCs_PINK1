# ================= LIBRARIES ===================
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
library(viridis)
set.seed(123)


# ================ READ ARGS =================
args = commandArgs(trailingOnly=TRUE)
print(args)

source("workflow/scripts/Functions.R")
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else {
    print("Arguments Passed")
}
# ---------------------------------------------



file <- args[1]
sample <- args[2]
condition_name <- strsplit(args[2],"_")[[1]][1]
date <- strsplit(args[2],"_")[[1]][2]
output_file <- args[3]


print(condition_name)
# condition_name <- args[2]
# sample <- args[3]
# date <- args[4]


imputation = FALSE
remove_mt=FALSE
remove_ribsomal=FALSE
remove_rb = TRUE

n_cores=1
elbow = TRUE
SCT=TRUE
criteria_pass=3
min.cells <- 10
min.features <- 200
data_10x <- FALSE
vars.to.regress=c("nCount_RNA")
outlier_detector = "MAD"

# # ------------------------------------------------
colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")[c(5,1,2,3,4,9,6,7,8)]
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

object <- load_files(file = file, data_10x = data_10x)

if (elbow == TRUE) {
    Before_col_names <- colnames(object)
    object <- elbow_calc(object, sample, iter_qc)
    After_col_names_Elbow <- colnames(object)
    Rm_Cell_N_Elbow <- Before_col_names[!Before_col_names %in% After_col_names_Elbow]
} else {
    After_col_names_Elbow <- colnames(object)
    Rm_Cell_N_Elbow <- c()
}


Outlier_object <- Advanced_Outlier_detection(object, 
        condition = sample, 
        criteria_pass = criteria_pass)

#saveRDS(Outlier_object,file=output_file)


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


Seurat <- CreateSeuratObject(counts = object, project = sample, 
    min.cells = min.cells, min.features = min.features, 
    meta.data = data.frame(percent.mito = percent.mito, 
        percent.rb = percent.rb,
        condition=condition_name,
        date=date))

Seurat$stim <- sample
After_col_names_SEURAT <- colnames(Seurat)
Rm_Cell_N_SEURAT <- After_col_names_MAD[!After_col_names_MAD %in% 
    After_col_names_SEURAT]
Rm_Cell_N_Elbow <- c(Rm_Cell_N_Elbow, Rm_Cell_N_SEURAT)

write.table(Rm_Cell_N_Elbow, paste0("result/Preprocess/",sample,"/Cells_Removed.tsv"), sep = "\t")

after_genes <- rownames(Seurat@assays$RNA@counts)
rm_genes <- before_genes[!before_genes %in% after_genes]
write.table(rm_genes, paste0("result/Preprocess/",sample,"/Genes_Removed.tsv"), sep = "\t")

Seurat <- SCTransform(object = Seurat, verbose = FALSE, vars.to.regress = vars.to.regress)

saveRDS(Seurat,file=output_file)
