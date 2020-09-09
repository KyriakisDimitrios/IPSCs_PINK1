
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

library(STRINGdb)
library(igraph)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
olor_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)


dir.create("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/4.Network_Validation/")
setwd("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/4.Network_Validation/")

Combined <- readRDS("/home/users/dkyriakis/PhD/Projects/IPSCs_pink1/1.Preprocess/IPSCs_Combined.rds")



# === READ DF GROUPS
g1 <- read.table("Conserved_all.txt")
g2 <- read.table("Conserved_all_alt.txt")
g3 <- read.table("Conserved_3.txt")

all_g <- rownames(Combined@assays$RNA@counts)
write.table(all_g,"All_genes.txt")
length(all_g)
all_df_g <- unique(c(rownames(g1),rownames(g2),rownames(g3)))
length(all_df_g)
rest <- all_g[!all_g%in%all_df_g]


mean_l <- c()
median_l <- c()
degree_l <- list()
num_nodes <- c()
num_inter <- c()
string_db = STRINGdb$new(version="10",species=9606)

string_plot_net <- function(i,mean_l,median_l,degree_l){
    r_list1 <- sample(rest,292)
    write.table(r_list1,paste0('Random_Genes',i,".txt"),row.names = F,col.names = F)
    Genes <- string_db$mp(r_list1)
    #string_db$plot_network(Genes)
    inter_g <- string_db$get_interactions(Genes)
    full.graph <- string_db$get_subnetwork(Genes)
    degrees_f <- degree(full.graph)
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


pdf("Ours.pdf")
write.table(all_df_g,paste0('all_df_g',51,".txt"),row.names = F,col.names = F)
string_db = STRINGdb$new(version="10",species=9606)
Genes <- string_db$mp(all_df_g)
string_db$plot_network(Genes)
full.graph <- string_db$get_subnetwork(Genes)
degrees_f <- degree(full.graph)
mean_l <-c (mean_l,mean(degrees_f))
median_l <-c (median_l,median(degrees_f))
num_nodes <- c(num_nodes,length(degrees_f))
num_inter <- c(num_inter,ecount(full.graph))
degree_l[[51]] <-degrees_f
dev.off()

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



library(plyr)
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
    theme(legend.position = "top") + ylab("Probability") + xlab("Number of nodes")         
#    annotate(geom="text", x=225, y=0.045, label="Random",color="#00BFC4",size=4)+
#    annotate(geom="text", x=245, y=0.045, label="DEG",color="#F8766D",size=4)+


df <- data.frame("interactions"=rest_num_inter)
df$ord <- "RANDOM"

p6 <- df %>%
    ggplot( aes(x=interactions, color=ord)) +
    geom_density(fill="#00BFC4", color="#00BFC4",alpha=0.8)+
    geom_vline(xintercept = deg_num_inter,color="#F8766D",linetype="dashed",size = 2)+
    theme_cowplot()+
    theme(legend.position = "top") + ylab("Probability")  + xlab("Number of interactions")    +xlim(c(400,1550)) 
    #annotate(geom="text", x=840, y=0.003, label="Random",color="#00BFC4",size=4)+
    #annotate(geom="text", x=1360, y=0.003, label="DEG",color="#F8766D",size=4)+
    
#  #69b3a2
# #e9ecef



library(jsonlite)
library(purrr)
library(data.table)

dt_list <- map(degree_l, as.data.table)
dt <- rbindlist(dt_list, fill = TRUE, idcol = T)
colnames(dt)[2] <- "Degree"
dt$ord <- "RANDOM"
dt$ord[dt$.id==51]<- "DEG"
colnames(dt)[3] <- "Condition"

pdf("292.pdf")
plot(mean_l,median_l)
plot(mean_l,num_nodes)
plot(median_l,num_nodes)
dev.off()

p3 <- ggplot(dt, aes(x=.id, y=Degree,group=.id,fill=Condition)) + 
    geom_boxplot()+theme_cowplot()


# ggboxplot(dt, x = ".id", y = "V1",group=".id")+
#     stat_compare_means() +                                         # Global p-value
#     stat_compare_means(ref.group = 51, label = "p.signif",
#                        label.y = c(22, 29))      


p4 <- ggboxplot(dt, x = "Condition", y = "Degree",
                fill = "Condition")+ 
    stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 50) + theme_cowplot()  + ylab("Degree") 


my_comparisons <- list( c("DEG","RANDOM"))
dt2 <- dt[order(Condition),]
p4 <- ggboxplot(dt2, x = "Condition", y = "Degree",
                fill = "Condition")+ 
    stat_compare_means(comparisons = my_comparisons, ref.group = "RANDOM")  +
    #annotate(geom="text", x="DEG", y=90, label="***",color="black",size=5)+theme_cowplot()+
    theme(legend.position = "none")     + ylab("Degree") +xlab("Condition")     

    
p5 <- ggplot(dt, aes(x=Condition, y=Degree,group=Condition,fill=Condition)) + 
        geom_boxplot()+ 
    stat_compare_means(method = "wilcox.test")+      # Add global p-value
    stat_compare_means(label = "p.signif", method = "t.test",
                       ref.group = "RANDOM")+theme_cowplot()+
    theme(legend.position = "top")          



p12 <- ggarrange(plotlist=list(p1,p2),nrow = 1)
p45 <- ggarrange(plotlist=list(p4,p5),nrow = 1)

pdf("QC_292.pdf",height=12,width=8)
ggarrange(plotlist=list(p12,p3,p45),nrow = 3)
dev.off()


pdf("QC2_292.pdf",height=8,width=15)
p125 <- ggarrange(plotlist=list(p1,p2,p6,p4),ncol = 4)
ggarrange(plotlist=list(p125,p3),nrow = 2)
dev.off()

pdf("QC3_292.pdf",width=15,height=10)
ggarrange(plotlist=list(p2,p1,p6,p4),nrow = 2,ncol=2)
dev.off()


pdf("QC4_292.pdf",width=8,height=4)
p2
p1
p6
p4
dev.off()

