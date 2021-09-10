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

library(igraph)
library(plyr)
library(jsonlite)
library(purrr)
library(data.table)

list.of.packages <- c("STRINGdb")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
library(STRINGdb)


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
group_B <- args[2]
group_C<- args[3]
group_D <- args[4]
output_file <- args[5]
# ---------------------------------------------

Combined <- readRDS(input_file)

# ============== READ DF GROUPS ================
# gA <- read.table(args[1])
gB <- read.table(group_B)
gC <- read.table(group_C)
gD <- read.table(group_D)
DEGs <- unique(c(rownames(gB),rownames(gC),rownames(gD)))
all_g <- rownames(Combined@assays$RNA@counts)
rest <- all_g[!all_g%in%DEGs]
length(all_g)
# [1] 18097
length(DEGs)
# [1] 292
length(rest)
# [1] 17805
# -----------------------------------------------


mean_l <- c()
median_l <- c()
degree_l <- list()
num_nodes <- c()
num_inter <- c()
string_db = STRINGdb$new(version="11",species=9606)

string_plot_net <- function(i,mean_l,median_l,degree_l){
    # n_genes <- 292
    # r_list1 <- sample(rest,n_genes)
    # write.table(r_list1,paste0('Random_Genes',i,"_",n_genes,".txt"),row.names = F,col.names = F)
    # Rerun 
    r_list1<- read.table(paste0('DATA/Random_Genes/Random_Genes',i,".txt"),col.names = F)[,1]
    Genes <- string_db$mp(r_list1)
    #string_db$plot_network(Genes)
    inter_g <- string_db$get_interactions(Genes)
    full.graph <- string_db$get_subnetwork(Genes)
    degrees_f <- igraph::degree(full.graph)
    mean_l <-c (mean_l,mean(degrees_f))
    median_l <-c (median_l,median(degrees_f))
    num_nodes <- c(num_nodes,length(degrees_f))
    degree_l[[i]] <-degrees_f
    num_inter <- c(num_inter,ecount(full.graph))
    return(list("mean_l"=mean_l,"median_l"=median_l,"degree_l"=degree_l,"num_nodes"=num_nodes,"num_inter"=num_inter))
}

for(i in 1:50){
    return_res <- string_plot_net(i,mean_l,median_l,degree_l)
    mean_l <- return_res$mean_l
    median_l <- return_res$median_l
    num_nodes <- return_res$num_nodes
    degree_l <- return_res$degree_l
    num_inter <- return_res$num_inter
}

lapply(list(mean_l,median_l,degree_l,num_inter),length)


# pdf("result/Networks/Ours.pdf")
#write.table(DEGs,paste0('DEGs',51,".txt"),row.names = F,col.names = F)
# DEGs <- read.table(paste0('DEGs',51,".txt"),col.names = F)[,1]
string_db = STRINGdb$new(version="11.0",species=9606)
Genes <- string_db$mp(DEGs)
print("OK1")
# string_db$plot_network(Genes)
# print("OK2")

full.graph <- string_db$get_subnetwork(Genes)
degrees_f <- igraph::degree(full.graph)
mean_l <-c (mean_l,mean(degrees_f))
median_l <-c (median_l,median(degrees_f))
num_nodes <- c(num_nodes,length(degrees_f))
num_inter <- c(num_inter,ecount(full.graph))
degree_l[[51]] <-degrees_f
# dev.off()

print("OK3")

rest_degree <- unlist(degree_l[1:50])
deg_degree <- unlist(degree_l[51])
df <- data.frame("degree"=unlist(degree_l))
length(deg_degree)
length(unlist(degree_l))
n_random <- length(unlist(degree_l)) - length(deg_degree)
df$ord <- "Random"
df$ord[n_random:length(unlist(degree_l))] <- "DEG"

rest_num_nodes <- unlist(num_nodes[1:50])
deg_num_nodes <- unlist(num_nodes[51])

rest_num_inter <- unlist(num_inter[1:50])
deg_num_inter <- unlist(num_inter[51])


d <- density(rest_num_nodes)


print("OK4")

mu <- ddply(df, "ord", summarise, grp.mean=mean(degree))


p1 <- df %>%
    ggplot( aes(x=degree, fill=ord)) +
    geom_density(alpha=0.8)+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=ord),size=2,
               linetype="dashed")+theme_cowplot()+
    theme(legend.position = "none")   + ylab("Degree distribution") +xlab("Degree")       



df <- data.frame("nodes"=rest_num_nodes)
df$ord <- "RANDOM"

p2 <- df %>%
    ggplot( aes(x=nodes, color=ord)) +
    geom_density(fill="#00BFC4", color="#00BFC4",alpha=0.8)+
    geom_vline(xintercept = deg_num_nodes,color="#F8766D",linetype="dashed", size = 2)+theme_cowplot()+
    theme(legend.position = "top") + ylab("Probability") + xlab("Number of nodes")+         
    annotate(geom="text", x=190, y=0.045, label="Random",color="#00BFC4",size=7)+
    annotate(geom="text", x=240, y=0.045, label="DEG",color="#F8766D",size=7)+xlim(170,260)


df <- data.frame("interactions"=rest_num_inter)
df$ord <- "RANDOM"

p6 <- df %>%
    ggplot( aes(x=interactions, color=ord)) +
    geom_density(fill="#00BFC4", color="#00BFC4",alpha=0.8)+
    geom_vline(xintercept = deg_num_inter,color="#F8766D",linetype="dashed",size = 2)+
    theme_cowplot()+
    theme(legend.position = "top") + ylab("Probability")  + xlab("Number of interactions")    +xlim(c(400,1650)) +
    annotate(geom="text", x=300, y=0.003, label="Random",color="#00BFC4",size=7)+
    annotate(geom="text", x=1300, y=0.003, label="DEG",color="#F8766D",size=7)+xlim(100,1600)
    
#  #69b3a2
# #e9ecef





dt_list <- map(degree_l, as.data.table)
dt <- rbindlist(dt_list, fill = TRUE, idcol = T)
colnames(dt)[2] <- "Degree"
dt$ord <- "RANDOM"
dt$ord[dt$.id==51]<- "DEG"
colnames(dt)[3] <- "Condition"



# p3 <- ggplot(dt, aes(x=.id, y=Degree,group=.id,fill=Condition)) + 
#     geom_boxplot()+theme_cowplot()


# ggboxplot(dt, x = ".id", y = "V1",group=".id")+
#     stat_compare_means() +                                         # Global p-value
#     stat_compare_means(ref.group = 51, label = "p.signif",
#                        label.y = c(22, 29))      

my_comparisons <- list( c("DEG","RANDOM"))

# p4 <- ggboxplot(dt, x = "Condition", y = "Degree",
#                 fill = "Condition")+ 
#     stat_compare_means(comparisons = my_comparisons)+ theme_cowplot()  + ylab("Degree")  # Add pairwise comparisons p-value
    # stat_compare_means(label.y = 50) 



dt2 <- dt[order(Condition),]
p4 <- ggboxplot(dt2, x = "Condition", y = "Degree",
                fill = "Condition")+ 
    stat_compare_means(comparisons = my_comparisons, ref.group = "RANDOM")  +
    #annotate(geom="text", x="DEG", y=90, label="***",color="black",size=5)+theme_cowplot()+
    theme(legend.position = "none")     + ylab("Degree") +xlab("Condition")     

pdf(output_file)
p2
p6
p4
dev.off()


# p5 <- ggplot(dt, aes(x=Condition, y=Degree,group=Condition,fill=Condition)) + 
#         geom_boxplot()+ 
#     stat_compare_means(method = "wilcox.test")+      # Add global p-value
#     stat_compare_means(label = "p.signif", method = "t.test",
#                        ref.group = "RANDOM")+theme_cowplot()+
#     theme(legend.position = "top")          




# cor_list<- readRDS("Correlation_Networks.rds")
# p26<-ggarrange(plotlist=list(p2,p6),nrow = 2)
# g14 <- p1 + annotation_custom(ggplotGrob(p4), xmin = 30, xmax = 80, 
#                        ymin = 0.04, ymax = 0.11)
# sup6_p2 <- ggarrange(plotlist=list(p26,g14),nrow = 1)

# FigSup_1 <- cor_list[[1]]
# FigSup_2 <- cor_list[[2]]
# FigSup_3 <- cor_list[[3]]
# FigSup_4 <- cor_list[[4]]


# sup6_p1 <- ggarrange(plotlist=list(FigSup_1,sup6_p2),ncol=2)

# sup6_p3  <- ggarrange(plotlist=list(FigSup_2,FigSup_3,FigSup_4),ncol=3)


# pdf("result/Network/Network_Sup_Figures.pdf",width=20,height=10)
# FigSup_1 
# FigSup_2 
# FigSup_3 
# FigSup_4 
# sup6_p2
# dev.off()






# pdf("Sup6")
# ggarrange(plotlist=list(p26,g14),nrow = 1)
# dev.off()













# p12 <- ggarrange(plotlist=list(p1,p2),nrow = 1)
# p45 <- ggarrange(plotlist=list(p4,p5),nrow = 1)

# pdf("QC_292.pdf",height=12,width=8)
# ggarrange(plotlist=list(p12,p3,p45),nrow = 3)
# dev.off()


# pdf("QC2_292.pdf",height=8,width=15)
# p125 <- ggarrange(plotlist=list(p1,p2,p6,p4),ncol = 4)
# ggarrange(plotlist=list(p125,p3),nrow = 2)
# dev.off()

# pdf("QC3_292.pdf",width=15,height=10)
# ggarrange(plotlist=list(p2,p1,p6,p4),nrow = 2,ncol=2)
# dev.off()


# pdf("QC4_292.pdf",width=8,height=4)
# p2
# p1
# p6
# p4
# dev.off()

# p6+ annotation_custom(ggplotGrob(p2), xmin = 1, xmax = 3, 
#                        ymin = -0.3, ymax = 0.6)


# pc <- ggarrange(plotlist=list(p6,p4),ncol=1)
# pb <- ggarrange(plotlist=list(p2,p1),ncol=1)

# ggarrange(plotlist=list(FigSup_1,FigSup_3,FigSup_4),c("d","e","f"),ncol=3)

# ggarrange(plotlist=list(FigSup_2,FigSup_3,FigSup_4),c("d","e","f"),ncol=3)



