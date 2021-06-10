
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
library(viridis)
## ----Venn------------------------------------------------------------------------------------
library(VennDiagram)
library(EnhancedVolcano)

colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

setwd("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021")
load("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/Michi_Data/2020-04-17_seurat_elbow_TRUE_Mito-FALSE_Ribo-FALSE_SCT-TRUE_criteria_pass-3/IPSCs_PINK.RData")

Combined$Timepoints <- as.vector(Combined$condition)
Combined$Timepoints[grep("IPSC",Combined$condition)] <-"iPSCs"
Combined$Timepoints[grep("D06",Combined$condition)] <-"Day06"
Combined$Timepoints[grep("D10",Combined$condition)] <-"Day10"
Combined$Timepoints[grep("D15",Combined$condition)] <-"Day15"
Combined$Timepoints[grep("D21",Combined$condition)] <-"Day21"
Combined$Timepoints <- factor(Combined$Timepoints,c("iPSCs","Day06","Day10","Day15","Day21"))

Controls <- readRDS("../Reviewes/All_controls.rds")

Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D26","D35","D50"))

Controls_DM <- subset(Combined,subset= Treatment !="PINK")
Controls_AD10 <- subset(Combined,subset= Timepoints !="Day10")
Controls_CD10 <- subset(Controls_DM,subset= Timepoints !="Day10")

DefaultAssay(Controls) <- "RNA"
DefaultAssay(Controls_DM) <- "RNA"
DefaultAssay(Combined) <- "RNA"
DefaultAssay(Controls_AD10) <- "RNA"
DefaultAssay(Controls_CD10) <- "RNA"


# TH EXPRESSION 
p1 <- DotPlot(Combined,
    features = "TH",
    group.by = "Timepoints")+coord_flip()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+ 
    theme(text = element_text(size = 12))
p2 <- DotPlot(Controls,
    features = "TH",
    group.by = "Timepoints")+coord_flip()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+ 
    theme(text = element_text(size = 12))

as.vector(p1$data$pct.exp)
p2$data$pct.exp


th_expressed_all <- split(Combined@assays$RNA@counts["TH",] > 0, Combined$Timepoints)
th_expressed_controls <- split(Controls@assays$RNA@counts["TH",] > 0, Controls$Timepoints)

unlist(lapply(th_expressed_all,mean))*100
unlist(lapply(th_expressed_controls,mean))*100



# ======= Figure 1 B
fig1_markers <- c("MYC","POU5F1",
"PTCH1","FZD7",
"HES1","OTX2",
"SLIT2","LMX1A",
"DCX","DDC")

pdf("Figure1B.pdf",width= 8)
DefaultAssay(Controls_AD10) <- "RNA"
Controls_AD10$Timepoints <- factor(Controls_AD10$Timepoints)
Controls_AD10<-ScaleData(Controls_AD10,fig1_markers)
DoHeatmap(Controls_AD10,fig1_markers,hjust=0.5,angle=0,
    group.by = "Timepoints", raster = FALSE)+
    scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))

Controls_CD10$Timepoints <- factor(Controls_CD10$Timepoints)
Controls_CD10<-ScaleData(Controls_CD10,fig1_markers)
DoHeatmap(Controls_CD10,fig1_markers,hjust=0.5,angle=0,
    group.by = "Timepoints", raster = FALSE)+
    scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))
dev.off()




# ==== Fig 2 IPSCs Reviewers
ipsc_markers<-c("MYC",
"POU5F1",
# "NANOG",
"L1TD1",
"USP44",
"TDGF1",
"POLR3G",
"TERF1" ,
"IFITM1",
"PRDX1" ,
# "LCK",
"DNMT3B",
"DPPA4",
# "IDO1",
"LIN28A"
)

grep("TDGF",rownames(Combined@assays$RNA@counts))
pdf("Fig2B.pdf",width= 8)
DefaultAssay(Combined) <- "RNA"
Combined<-ScaleData(Combined,ipsc_markers)
# DoHeatmap(Combined,ipsc_markers,hjust=0.5,angle=0,
#     group.by = "Timepoints", raster = FALSE)+
#     scale_fill_viridis(option="inferno")+ 
#     theme(text = element_text(size = 20))
Controls_DM<-ScaleData(Controls_DM,ipsc_markers)
DoHeatmap(Controls_DM,ipsc_markers,hjust=0.5,angle=0,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))
dev.off()




# ================================== MODULE SCORES 
# ================================== MODULE SCORES 
# ================================== MODULE SCORES 

IPSC <- ipsc_markers
Rgl1<- c(
"FABP7",
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"RSPO2",
"MSX1",
"SOX6"
)

Rgl3<- c(
"FABP7",
"SOX2",
"WNT5A",
"CORIN",
"BNC2"
)

ProgM <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"CORIN",
"DCC"
)

ProgFPs <- c(
"SOX2",
"FOXA2",
"LMX1A",
"OTX2",
"WNT5A",
"DCC"
)



NProg <- c(
"SOX2",
"FOXA2",
"OTX2",
"WNT5A",
"ASCL1",
"NEUROD1",
"NEUROG2",
"TBB3",
"DDC",
"DCX")

NbM <- c(
"FOXA2",
"NEUROD1",
"NEUROG2",
"TBB3",
"DCX")


DA0<- c(
"FOXA2",
"LMX1A",
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH"
)

DA1<- c(
"FOXA2",
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH",
"BNC2",
"SLC18A2",
"SLC6A3"
)

DA2<- c(
"TBB3",
"DDC",
"DCX",
"NR4A2",
"PBX1",
"PITX3",
"EN1",
"TH",
"BNC2",
"SLC18A2",
"SLC6A3",
"LM03",
"ALDH1A1",
"SOX6"
)


split_list <- list(IPSC =IPSC,Rgl=unique(c(Rgl3,Rgl1)),Prog=unique(c(ProgM,ProgFPs)),NProg=unique(c(NProg,NbM)),DA=unique(c(DA0,DA1,DA2)))
unlist_feat <- unlist(split_list)
colnames(Controls@meta.data)
Controls@meta.data <- Controls@meta.data[,c(1:20)]
Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D26","D35","D50"))
Controls <-AddModuleScore(Controls,features = split_list,name = "DA_s")
colnames(Controls@meta.data)[22:26]<-names(split_list)

pa <- DotPlot(Controls,features = c("IPSC","Rgl","Prog","NProg","DA"  ),group.by="Timepoints")

toplot <- pa$data[,c(2,3,4)]
#toplot$features.plot <- factor(toplot$features.plot,levels=c("Rgl3","Rgl1" ,   "ProgM", "NProg", "NbM",   "DA0",   "DA1"  ,"DA2"  ))
toplot$features.plot <- factor(toplot$features.plot,levels=c("IPSC","Rgl","Prog","NProg","DA"  ))



scores_meta <-  Controls@meta.data[22:26]
scores_meta_assign <- apply(scores_meta, 1,function(x) {
    if(x[which.max(x)]<0){
        col_n <- "Undefined"
    }else{
        col_n <- colnames(scores_meta)[which.max(x)]
    }
    return(col_n)
})
Controls$Phase_Dif <- factor(scores_meta_assign,
    levels=c("Undefined","IPSC","Rgl","Prog","NProg","DA"))


table_pl <- rbind(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D06"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D06"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D10"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D10"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D15"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D15"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D21"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D21"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D26"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D26"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D35"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D35"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D50"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D50"]))
)
rownames(table_pl) <-  c("IPSCs","D06","D10","D15","D21","D26","D35"  ,"D50")
toplot <- reshape2::melt(table_pl)
colnames(toplot) <- c("TimePoint","Group","value")






pdf("Cell_Percentage.pdf",width=12)


DotPlot(Controls,
    features = rev(unique(unlist_feat)),
    group.by = "Timepoints")+coord_flip()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+ 
    theme(text = element_text(size = 12))
dev.off()

TH<-c(1.92222E-05,	
    2.17257E-05,
    4.86911E-05,
    0.000110361,
    0.016417411,
    0.023048918)
TH_se <-c (1.05E-05,
3.08E-06,
1.77E-05,
1.52E-05,
7.11E-04,
9.72E-03)
TH_sd <-c(
1.82E-05,
6.16E-06,
3.54E-05,
2.64E-05,
1.23E-03,
1.37E-02)
TH_a<-c(
4.01E-05,    1.02E-05,    7.32E-06,    
1.98E-05,    3.01E-05,    1.54E-05,    2.16E-05,
1.89E-05,    9.66E-05,    2.55E-05,    5.37E-05,
1.37E-04,    1.11E-04,    8.39E-05,    
1.50E-02,    1.72E-02,    1.70E-02,    
1.33E-02,    3.28E-02)        
TH_t <- c(
0,0,0,
6,6,6,6,
10,10,10,10,
14,14,14,
21,21,21,
50,50)



ALDH1A1<-c(
0.000132344,
1.14092E-05,
3.75322E-05,
0.000430088,
0.001287628,
0.001430329)
ALDH1A1_a<-c(
1.43E-04,    1.25E-04,    1.29E-04,    
2.16E-06,    1.64E-05,    2.60E-06,    2.45E-05,
5.39E-05,    3.12E-05,    2.57E-05,    3.93E-05,
4.09E-04,    4.89E-04,    3.93E-04,    
1.69E-03,    1.10E-03,    1.07E-03,    
1.30E-03,    1.56E-03)        
ALDH1A1_t <- c(
0,0,0,
6,6,6,6,
10,10,10,10,
14,14,14,
21,21,21,
50,50)
ALDH1A1_sd <-c(
9.34E-06,
1.09E-05,
1.23E-05,
5.13E-05,
3.47E-04,
1.80E-04)
ALDH1A1_se <-c(
5.39E-06,
5.47E-06,
6.14E-06,
2.96E-05,
2.00E-04,
1.27E-04)



LMX1A_a<- c(
6.90E-06,    1.37033E-05, 3.56934E-06, 3.43233E-06, 
1.99E-03,    0.00175891,  0.001541975, 0.003139975, 0.001536432,
7.29E-03,    0.00659595,  0.008470973, 0.007002762, 0.007078652,
3.47E-02,    0.034986945, 0.046455696, 0.022744845, 
7.96E-02,    0.056238185, 0.083394383, 0.099128269, 
3.56E-02,    0.031324503, 0.039897134)     
LMX1A_t <- c(
0,0,0,0,
6,6,6,6,6,
10,10,10,10,10,
14,14,14,14,
21,21,21,21,
50,50,50)
length(LMX1A_a)
length(LMX1A_t)
LMX1A<- c(6.90166E-06,
0.001994323,
0.007287084,
0.034729162,
0.079586946,
0.035610819)
LMX1A_sd <-c(
5.8908E-06,
0.000770761,
0.000817219,
0.011857527,
0.021697056,
0.006061766)
LMX1A_se <-c(3.40105E-06,
0.000385381,
0.00040861,
0.006845947,
0.012526801,
0.004286316)




mat <- cbind(TH,ALDH1A1,LMX1A)
rownames(mat) <- c("0","6","10","14","21","50")
df <- as.data.frame(reshape2::melt(mat))
colnames(df) <- c("Timepoint","Feature","Expression")


TH_mat <- as.data.frame(cbind(TH_a,TH_t))
colnames(TH_mat)<-c("Expression", "TH_t")

LMX1A_mat <- as.data.frame(cbind(LMX1A_a,LMX1A_t))
colnames(LMX1A_mat)<-c("Expression", "LMX1A_t")

ALDH1A1_mat <- as.data.frame(cbind(ALDH1A1_a,ALDH1A1_t))
colnames(ALDH1A1_mat)<-c("Expression", "ALDH1A1_t")


fig3b1 <- df %>% filter(Feature=="TH") %>% 
  add_column(sd=TH_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    geom_point(size=1,aes(TH_t,Expression),data=TH_mat)+
    theme_minimal()+
    ggtitle("TH/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))


fig3b2 <- df %>% filter(Feature=="ALDH1A1") %>% 
  add_column(sd=ALDH1A1_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    geom_point(size=1,aes(ALDH1A1_t,Expression),data=ALDH1A1_mat)+
    theme_minimal()+
    ggtitle("ALDH1A1/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b3<- df %>% filter(Feature=="LMX1A") %>% 
  add_column(sd=LMX1A_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    geom_point(size=1,aes(LMX1A_t,Expression),data=LMX1A_mat)+
    theme_minimal()+
    ggtitle("LMX1A/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b <- ggarrange(plotlist = list(fig3b1,fig3b2,fig3b3),ncol=3)


fig3b1 <- df %>% filter(Feature=="TH") %>% 
  add_column(se=TH_se)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=3)+
    geom_point(size=1,aes(TH_t,Expression),data=TH_mat)+
    theme_minimal()+
    ggtitle("TH/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))


fig3b2 <- df %>% filter(Feature=="ALDH1A1") %>% 
  add_column(se=ALDH1A1_se)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=3)+
    geom_point(size=1,aes(ALDH1A1_t,Expression),data=ALDH1A1_mat)+
    theme_minimal()+
    ggtitle("ALDH1A1/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b3<- df %>% filter(Feature=="LMX1A") %>% 
  add_column(se=LMX1A_se)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=3)+
    geom_point(size=1,aes(LMX1A_t,Expression),data=LMX1A_mat)+
    theme_minimal()+
    ggtitle("LMX1A/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b2 <- ggarrange(plotlist = list(fig3b1,fig3b2,fig3b3),ncol=3)


fig3b
fig3b2
pdf("QPCR_SD_SE.pdf",height=4,width=12)
fig3b
fig3b2
dev.off()


library(tidyverse)
fig3b1 <- df %>% filter(Feature=="TH") %>% 
  add_column(sd=TH_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=1)+
    geom_point(size=1,aes(TH_t,TH),data=TH_mat)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    theme_minimal()+
    ggtitle("TH/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b2 <- df %>% filter(Feature=="ALDH1A1") %>% 
  add_column(sd=ALDH1A1_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=1)+
    geom_point(size=1,aes(ALDH1A1_t,ALDH1A1),data=ALDH1A1_mat)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    theme_minimal()+
    ggtitle("ALDH1A1/GAPDH")+ 
    theme(text = element_text(size = 15))+
  theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))

fig3b3 <- df %>% filter(Feature=="LMX1A")%>% 
  add_column(sd=LMX1A_sd)%>%
  ggplot(aes(Timepoint,Expression)) +
    geom_line() +
    geom_point(size=1)+
    geom_point(size=1,aes(LMX1A_t,LMX1A),data=LMX1A_mat)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=3)+
    theme_minimal()+
    ggtitle("LMX1A/GAPDH")+ 
    theme(text = element_text(size = 15))+
    theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))


fig3b <- ggarrange(plotlist = list(fig3b1,fig3b2,fig3b3),ncol=3)






fig3d <- DotPlot(Controls,
    features = rev(names(split_list)),
    group.by = "Timepoints")+coord_flip()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+ 
    theme(axis.text=element_text(size=11))+
    theme(legend.title=element_text(size=10))+
    theme(legend.position = "right")+xlab("")+ylab("")


fig3e <-ggplot(data=toplot, aes(x=TimePoint, y=value, fill=Group)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    # geom_text(aes(label=round(value*100,0)), vjust=1.6, color="white",
    #           position = position_dodge(0.9), size=2.5)+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")+ 
    theme(axis.text=element_text(size=11),axis.text.x = element_text(angle = 45))
    # +theme(legend.position ="top")



# ==== Fig 3 MDA DIFFERENTIATION

mDA2 <- c(
"PTCH1"	,
"FZD7"	,
"FABP7",
"CTNNB1",
"SHH"	,
"SOX2"	,
"FOXA2"	,
"OTX2"	,
"PBX1"	,
"HES1"	,
"NTN1",
"SOX6"	,
"DCX"	,
"DDC"	,
"TCF12"	,
"ALCAM",
"SLIT2"	,
"LMO3",
"LMX1A"	,
"WNT5A"	,
"RSPO2"	,
"MSX1"	,
"NEUROG2",
"MSX1"	,
"DMRTA2",	
"PITX2",
"ASCL1",
"PITX3"	,
"SLIT1",
"TH",
"NR4A2"
)
DefaultAssay(Controls) <- "RNA"
Controls<-ScaleData(Controls,mDA2)
fig3c <- DoHeatmap(Controls,features=mDA2,
    angle=0,hjust=0.5,
    group.by = "Timepoints",
    raster = FALSE)+
    scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 11))

fig3b <- ggarrange(plotlist = list(fig3b1,fig3b2,fig3b3),ncol=3,heights=c(5,5,5),widths=c(5,5,5))

pdf("Fig3B.pdf",width=10,height=4)
fig3b
dev.off()
pdf("Fig3C.pdf",width=10,height=6)
fig3c
dev.off()
pdf("Fig3D.pdf",width=10,height=6)
fig3d
dev.off()
pdf("Fig3E.pdf",width=10,height=6)
fig3e
dev.off()


fig3bc <- ggarrange(plotlist =list(NULL,NULL,NULL ,fig3c),
    labels=c("","B","","C"),
    ncol=4,heights=c(8,4,8,8),widths=c(0.1,8,1,8),
    vjust=1.6,hjust=0.3)

fig3de <- ggarrange(plotlist =list(NULL,fig3d,NULL ,fig3e),
    labels=c("","D","","E"),
    ncol=4,heights=c(6,6,6,6),widths=c(0.1,8,1,8),
    vjust=1.3,hjust=0.3)

pdf("Figure32.pdf",width=15)
ggarrange(plotlist =list(NULL ,fig3bc,NULL ,fig3de),ncol=1,
    vjust=1,hjust=0.3,heights=c(2,8,1,6))
dev.off()






# ======================= FIGURE 4 ==========================
mDA <- c("TCF12", "ALCAM", "PITX2", "DDC", "ASCL1")
mDA <- mDA[mDA %in% rownames(Combined@assays$RNA@counts)]
p0 <- DimPlot(Combined,group.by = "condition",cols=color_cond)
results <- ICSWrapper::scatter_gene(object = Combined,features = c("TH","KCNJ6"),ncol =1 ,nrow =2)+ theme_void()
fig4b <- ggarrange(plotlist=list(p0, results), widths=c(2/3, 1/2), ncol=2, nrow=1)
pdf("Fig4B.pdf",width=10,height=6)
fig4b
dev.off()


# ======================= FIGURE 5 ==========================
dirs_pairs <- list.dirs("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021",full.names = TRUE )[-1]
dirs_pairs <- grep('IPSC|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)
dirs_pairs <- grep('D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)


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
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.01 & abs(ctrl_l1$avg_logFC) >0.1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.01 & abs(pink_l1$avg_logFC) >0.1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]])
}
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)

genes_fig5a <- c(cntrl_intesect,pink_intesect)
Conserved_M <- subset(Combined,subset= Timepoints %in% c("Day06","Day15","Day21"))

Conserved_M$condition<- factor(Conserved_M$condition)
Conserved_M$condition1<- as.vector(Conserved_M$condition)
Conserved_M$condition1 <- gsub("Control", "C", Conserved_M$condition1)
Conserved_M$condition1 <- gsub("PINK1", "P", Conserved_M$condition1)
Conserved_M <- NormalizeData(Conserved_M)
Conserved_M <- ScaleData(Conserved_M,genes_fig5a)
color_h = CustomPalette(low = "#ff0905", 
                mid = "#000000", 
                high = "#ffff05", 
                k = 100)
             
fig5a <- DoHeatmap(Conserved_M,features=genes_fig5a,
    angle=0,hjust=0.5,
    group.by = "condition1",
    raster = FALSE)+
    scale_fill_gradientn(colours = color_h)+ 
    theme(text = element_text(size = 12))

pdf("Fig5A.pdf",width=10,height=12)
fig5a
dev.off()



# ===================== OPEN FILES TAKE THE P.Adj- FC Genes

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

day06 <- c(df_return_nt_cntrl[[1]])
day15 <- c(df_return_nt_cntrl[[2]])
day21 <- c(df_return_nt_cntrl[[3]])

p1<- ggVennDiagram(list(Day06=day06, Day15=day15,Day21=day21),label_txtWidth = 70, label_color = "firebrick", label_size = 15,set_size=15 )+  
    theme(text = element_text(size=30))


p1 <- ggVennDiagram::ggVennDiagram(list(Day06=day06, Day15=day15,Day21=day21), show_intersect = TRUE)

fig5b1 <- p1 +
geom_label(aes(x = 2, y = 8.5, label = "Day06"),size=6,label.size = NA, fill = "white")+
geom_label(aes(x = -3.2, y = -4.6, label = "Day15"),size=6,label.size = NA, fill = "white")+
geom_label(aes(x = 7.2, y = -4.6, label = "Day21"),size=6,label.size = NA, fill = "white")+ 
geom_label(aes(x = 2, y = 2.5, label = "SNHG5"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 2.1, label = "DLK1"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 1.7, label = "PGK1"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 1.3, label = "CCDC144NL.AS1"),size=5.5,label.size = NA)+
geom_label(aes(x = 2, y = 0.9, label = "LGI1"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 0.4, label = "(3%)"),label.size = NA,size=5.5, fill = "white")+NoLegend()+
ggtitle("Downregulated")+ 
theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))#,hjust = 0.5))


day06 <- c(df_return_nt_pink[[1]])
day15 <- c(df_return_nt_pink[[2]])
day21 <- c(df_return_nt_pink[[3]])
p1 <- ggVennDiagram::ggVennDiagram(list(Day06=day06, Day15=day15,Day21=day21), show_intersect = TRUE)



fig5b2 <-p1 +
geom_label(aes(x = 2, y = 8.5, label = "Day06"),size=6,label.size = NA, fill = "white")+
geom_label(aes(x = -3.2, y = -4.6, label = "Day15"),size=6,label.size = NA, fill = "white")+
geom_label(aes(x = 7.2, y = -4.6, label = "Day21"),size=6,label.size = NA, fill = "white")+
geom_label(aes(x = 2, y = 1.6, label = "MTRNR2L1\n"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 1.2, label = "S100A6"),label.size = NA,size=5.5, fill = "white")+
geom_label(aes(x = 2, y = 0.6, label = "(1%)"),label.size = NA,size=5.5, fill = "white")+NoLegend()+
ggtitle("Upregulated")+ 
theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))


Fig5b <- ggarrange(plotlist=list(fig5b1,fig5b2),ncol=2) 
pdf("Fig5B.pdf",width=15,height=10)
Fig5b
dev.off()


Fig5b+theme(element_text(size=10))



# ================== Figure 5
enrich <- read.csv("Enrichment.csv")
enrich$Gene_Ratio <- enrich$observed.gene.count /  enrich$background.gene.count 
ggplot(enrich,aes(x=Gene_Ratio,y=term.description,
    size=observed.gene.count,
    color=false.discovery.rate))+
    geom_point()+
    scale_color_viridis(option="plasma")

enrich$Database <-  as.vector(enrich$X.term.ID)
enrich$Database[grep("hsa0",as.vector(enrich$X.term.ID))] <- "KEGG pathways"
enrich$Database[5:8] <- "Molecular Function (Gene Ontology)"
enrich$Database[grep("HSA-",as.vector(enrich$X.term.ID))] <- "Reactome Pathways"
enrich$Database[grep("KW-",as.vector(enrich$X.term.ID))] <- "Annotated Keywords (UniProt)"
enrich$Database[grep("GO",as.vector(enrich$X.term.ID))] <- "Biological Process (Gene Ontology)"

enrich$fdr <-enrich$false.discovery.rate
toplot <- enrich[enrich$Database!="Biological Process (Gene Ontology)",]
enrich$term.description <- gsub("metabolic process", "MP", enrich$term.description)


p1<-ggplot(enrich,aes(x=observed.gene.count,y=term.description,
    size=Gene_Ratio,
    color=fdr))+
    geom_point()+
    scale_color_continuous(low="red", high="blue",
            guide=guide_colorbar(reverse=TRUE))+
    ggforce::facet_col(vars(Database), scales = "free_y", space = "free")
    ylab("") + 
    theme(strip.text.y = element_text(angle=0))


p2<-ggplot(toplot,aes(x=observed.gene.count,y=term.description,
    size=Gene_Ratio,
    color=-log10(fdr)))+
    geom_point()+
    scale_color_continuous(low="red", high="blue",
            guide=guide_colorbar(reverse=TRUE))+
    # scale_color_viridis(option="plasma",direction=1)+
    ggforce::facet_col(vars(Database), scales = "free_y", space = "free")
    ylab("") + 
    theme(strip.text.y = element_text(angle=0))

toplot <- toplot[toplot$Database!="Annotated Keywords (UniProt)",]

p3<-ggplot(toplot,aes(x=observed.gene.count,y=term.description,
    size=Gene_Ratio,
    color=-log10(fdr)))+
    geom_point()+
    scale_color_continuous(low="red", high="blue",
            guide=guide_colorbar(reverse=TRUE))+
    # scale_color_viridis(option="plasma",direction=1)+
    ggforce::facet_col(vars(Database), scales = "free_y", space = "free")
    ylab("") + 
    theme(strip.text.y = element_text(angle=0))


pdf("Fig5c1v1.pdf",width=10,height=12);p1;dev.off()
pdf("Fig5c2.pdf",width=10);p2;dev.off()
pdf("Fig5c3.pdf",width=10);p3;dev.off()


# ================== Figure 6


#=========================================================================================
#=========================================================================================
#=========================================================================================
#=========================================================================================







# ================ DIFF CONTROL TIMEPOINTS 
# ================ DIFF CONTROL TIMEPOINTS 
# ================ DIFF CONTROL TIMEPOINTS 
DefaultAssay(Combined) <-"RNA"
Combined <- NormalizeData(Combined)
Combined <- ScaleData(Combined)
Controls_DM <- subset(Combined,subset= Treatment !="PINK")
DefaultAssay(Controls_DM) <- "RNA"
Idents(Controls_DM)<-as.factor(Controls_DM$condition)
pbmc.markers <- FindAllMarkers(object = Controls_DM,
                            assay ="RNA",min.pct =0.1,
                            logfc.threshold=0.1,
                            only.pos = TRUE,
                            test.use = "MAST",latent.vars = c("nCount_RNA"))


pbmc.markers$gene <- rownames(pbmc.markers)
qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(Controls_DM@assays$RNA@counts)[1])
pbmc.markers$qvalue <- qvalue

top <- pbmc.markers[pbmc.markers$p_val_adj<0.01 & abs(pbmc.markers$avg_logFC) >0.1,]
to_heat <- top %>% group_by(cluster) %>% top_n(n = 15, wt = abs(avg_logFC)) 
to_heat <- to_heat[order(to_heat$cluster),]
Controls_DM$condition <- as.factor(Controls_DM$condition) 

color_h = CustomPalette(low = "#ff0905", 
                mid = "#000000", 
                high = "#ffff05", 
                k = 100)
             
pdf("Figure3a2.pdf",width= 12)
DefaultAssay(Controls_DM) <- "RNA"
Controls_DM<-ScaleData(Controls_DM,to_heat$gene)
fig3a <- DoHeatmap(Controls_DM,features=to_heat$gene,
    angle=0,hjust=0.5,
    group.by = "Timepoints",
    raster = FALSE)+
    scale_fill_gradientn(colours = color_h)+
    theme(text = element_text(size = 11))
dev.off()









Idents(Combined) <- Combined$condition
mclapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,object=Combined,cl_combinations=cl_combinations,mc.cores=1)



# ====================================== COMMONS ===================================
dirs_pairs <- list.dirs("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021",full.names = TRUE )[-1]
dirs_pairs <- grep('IPSC|D06.*D06|D15.*D15|D21.*D21',dirs_pairs,value = TRUE)

df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()




Timepoint <- c("D06","D15","D21","IPSCs")
for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    timp <- Timepoint[iter]
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
    print(file)
    l1 <- read.table(file,header=TRUE)
    
    l1$cluster <- l1$avg_logFC
    l1$cluster[ l1$avg_logFC2<0] <- "PINK"
    l1$cluster[ l1$avg_logFC2>0] <- "Control"
    l1$p.adj_BY <-p.adjust(l1$p_val,method="BY")
    l1$p.adj_bonferroni <-p.adjust(l1$p_val,method="bonferroni")
    l1$p.adj_BH <- p.adjust(l1$p_val,method="BH")
    ctrl_l1 <- l1[grep("Control",l1$cluster),]
    pink_l1 <- l1[grep("PINK",l1$cluster),]
	all_l1 <-  l1
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.05 & abs(ctrl_l1$avg_logFC2) >1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.05 & abs(pink_l1$avg_logFC2) >1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]]) 
    
    pdf(paste0(timp,"Volcanos.pdf"))
    l1$avg_logFC2 <- -l1$avg_logFC2
    l1$avg_logFC <-  -l1$avg_logFC   
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+ggtitle(timp))
    min_vlc <- min(l1$avg_logFC)
    max_vlc <- max(l1$avg_logFC)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC',
                        y = 'p_val_adj',subtitle = "FC (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change"))
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p.adj_BY',subtitle = "BY FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change")+ggtitle(timp))

    dev.off()

}




  
library(UpSetR )
glmer_res <- readRDS("../glmer_pval_mat_10percent.rds")
glmer_res$max_FC2 <- -glmer_res$max_FC/log(2)
glmer_res$FC2 <- -glmer_res$FC/log(2)
glmer_res$p_val_adj <- p.adjust(glmer_res$pvalues_l,method="bonferroni",n=dim(Combined@assays$RNA@counts)[1]) 

pdf("GLMER_Volcanos.pdf")
min_vlc <- min(glmer_res$max_FC2)
max_vlc <- max(glmer_res$max_FC2)
plot(EnhancedVolcano(glmer_res,
                        lab = rownames(glmer_res),
                        x = 'max_FC2',
                        y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Log2 fold change")+ggtitle("GLMM"))

dev.off()


glmm_genes  <- rownames(glmer_res[glmer_res$p_val_adj<0.05 & abs(glmer_res$max_FC2)> 1,])
df_return_nt_all[[5]] <- glmm_genes
names(df_return_nt_all) <- c("D06","D15","D21","IPSCs","GLMM")
pdf("Common.pdf",width=12)
UpSetR::upset(fromList(df_return_nt_all),text.scale =1.9)
dev.off()

unlist(df_return_nt_cntrl)

# # ============= Intersect Common Genes
cntrl_intesect <- Reduce(intersect, df_return_nt_cntrl)
print(cntrl_intesect)
pink_intesect <- Reduce(intersect, df_return_nt_pink)
print(pink_intesect)

intersect(df_return_nt_all$D21,df_return_nt_all$GLMM )
setdiff(b,a)
setdiff(unlist(df_return_nt_all[1:4]),df_return_nt_all$GLMM )

setdiff(df_return_nt_all$GLMM ,unlist(df_return_nt_all[1:4]))






df_return_nt_cntrl <- list()
df_return_nt_pink <- list()
df_return_nt_all <- list()




Timepoint <- c("D06","D15","D21","IPSCs")
for (iter in 1:length(dirs_pairs)){
    dirs_iter <- dirs_pairs[iter]
    timp <- Timepoint[iter]
    file <- paste0(dirs_iter ,"/", dir(dirs_iter, "*.tsv"))
    print(file)
    l1 <- read.table(file,header=TRUE)
    
    l1$cluster <- l1$avg_logFC
    l1$cluster[ l1$avg_logFC2<0] <- "PINK"
    l1$cluster[ l1$avg_logFC2>0] <- "Control"
    l1$p.adj_BY <-p.adjust(l1$p_val,method="BY")
    l1$p.adj_bonferroni <-p.adjust(l1$p_val,method="bonferroni")
    l1$p.adj_BH <- p.adjust(l1$p_val,method="BH")
    ctrl_l1 <- l1[grep("Control",l1$cluster),]
    pink_l1 <- l1[grep("PINK",l1$cluster),]
	all_l1 <-  l1
    df_return_nt_cntrl[[iter]] <- as.vector(ctrl_l1[ctrl_l1$p_val_adj<0.01 & abs(ctrl_l1$avg_logFC) >0.1,"gene"])
    df_return_nt_pink[[iter]] <- as.vector(pink_l1[pink_l1$p_val_adj<0.01 & abs(pink_l1$avg_logFC) >0.1,"gene"])
    print(length(df_return_nt_cntrl[[iter]]))
    print(length(df_return_nt_pink[[iter]]))
    df_return_nt_all[[iter]] <- c(df_return_nt_cntrl[[iter]] ,df_return_nt_pink[[iter]]) 
    
    pdf(paste0(timp,"Volcanos.pdf"))
    l1$avg_logFC2 <- -l1$avg_logFC2
    l1$avg_logFC <-  -l1$avg_logFC   
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p_val_adj',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+ggtitle(timp))
    min_vlc <- min(l1$avg_logFC)
    max_vlc <- max(l1$avg_logFC)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC',
                        y = 'p_val_adj',subtitle = "FC (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change"))
    min_vlc <- min(l1$avg_logFC2)
    max_vlc <- max(l1$avg_logFC2)
    plot(EnhancedVolcano(l1,
                        lab = l1$gene,
                        x = 'avg_logFC2',
                        y = 'p.adj_BY',subtitle = "BY FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Neutral Log fold change")+ggtitle(timp))

    dev.off()

}







r_d00 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021/Control_IPSCs_PINK1_IPSCs/2021-04-12_TO_EXP_each_condition_Control_IPSCs_PINK1_IPSCs.tsv",header=T)
r_d06 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021/Control_D06_PINK1_D06/2021-04-12_TO_EXP_each_condition_Control_D06_PINK1_D06.tsv",header=T)
r_d15 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021/Control_D15_PINK1_D15/2021-04-12_TO_EXP_each_condition_Control_D15_PINK1_D15.tsv",header=T)
r_d21 <- read.table("C:/Users/dimitrios.kyriakis/Desktop/PINK1/12.04.2021/Control_D21_PINK1_D21/2021-04-12_TO_EXP_each_condition_Control_D21_PINK1_D21.tsv",header=T)

r_d00_f <- r_d00[r_d00$p_val_adj<0.05 & abs(r_d00$avg_logFC2) > 0.5,]
r_d00_f_c <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2>0])
r_d00_f_p <- as.vector(r_d00_f$gene[r_d00_f$avg_logFC2<0])
r_d06_f <- r_d06[r_d06$p_val_adj<0.05 & abs(r_d06$avg_logFC2) > 0.5,]
r_d06_f_c <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2>0])
r_d06_f_p <- as.vector(r_d06_f$gene[r_d06_f$avg_logFC2<0])
r_d15_f <- r_d15[r_d15$p_val_adj<0.05 & abs(r_d15$avg_logFC2) > 0.5,]
r_d15_f_c <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2>0])
r_d15_f_p <- as.vector(r_d15_f$gene[r_d15_f$avg_logFC2<0])
r_d21_f <- r_d21[r_d21$p_val_adj<0.05 & abs(r_d21$avg_logFC2) > 0.5,]
r_d21_f_c <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2>0])
r_d21_f_p <- as.vector(r_d21_f$gene[r_d21_f$avg_logFC2<0])


control_m<- unique(c(r_d00_f_c,r_d06_f_c,r_d15_f_c,r_d21_f_c))
length(r_d00_f_c)
length(r_d06_f_c)
length(r_d15_f_c)
length(r_d21_f_c)
length(control_m)
pink_m <- unique(c(r_d00_f_p,r_d06_f_p,r_d15_f_p,r_d21_f_p))
length(r_d00_f_p)
length(r_d06_f_p)
length(r_d15_f_p)
length(r_d21_f_p)
length(pink_m)

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
GOclusterplot <- 0
Keggclusterplot<-0
Reactomeclusterplot<-0
GOclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
Keggclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichKEGG")
Reactomeclusterplot <- compareCluster(ENTREZID~Condition,data=mydf, fun = "enrichPathway")


pdf("GO_Enrichment_DF_Clusters_log2FC_0.5.pdf",width=22,height=11)
dotplot(GOclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dotplot(Reactomeclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")
dev.off()

Fig5c <-dotplot(Keggclusterplot)+ ggplot2::facet_grid(~Condition, scales = "free")+ 
    theme(text = element_text(size = 15))
Fig5b <- ggarrange(plotlist=list(fig5b1,fig5b2),ncol=2) 
Fig5bc <-ggarrange(plotlist=list(Fig5b,Fig5c),ncol=1,labels=c("B","C"),vjust=1,hjust=0.5) 
Fig5 <- ggarrange(plotlist=list(NULL,Fig5a,NULL,Fig5bc,NULL),
    widths=c(0.1,8,0.5,8,1),
    ncol=5,
    labels=c("","A","","",""),vjust=1,hjust=0.3)

pdf("Figure5.pdf",height=12,width=25)
ggarrange(plotlist=list(NULL,Fig5),ncol=1,heights=c(0.05,10))
dev.off()










Conserved_M <- subset(Combined,subset= condition !="Control_D10")
DefaultAssay(Conserved_M) <- "RNA"
Conserved_M <-NormalizeData(Conserved_M)
Conserved_M <-ScaleData(Conserved_M)
Conserved_M$Timepoints <- as.vector(Conserved_M$condition)
Conserved_M$Timepoints[grep("IPSC",Conserved_M$condition)] <-"IPSCs"
Conserved_M$Timepoints[grep("D06",Conserved_M$condition)] <-"Day06"
Conserved_M$Timepoints[grep("D15",Conserved_M$condition)] <-"Day15"
Conserved_M$Timepoints[grep("D21",Conserved_M$condition)] <-"Day21"
Idents(Conserved_M) <- "Treatment"
markers <- FindConservedMarkers(Conserved_M,ident.1 = "Control",ident.2 = "PINK",grouping.var = "Timepoints",test.use="MAST",latent.vars="nCount_RNA",logfc.threshold=0.0)

index_fc <- c( sign(markers$Day06_avg_logFC)==sign(markers$Day15_avg_logFC) & sign(markers$Day15_avg_logFC)==sign(markers$Day21_avg_logFC))
sub_markers <- markers[markers$max_pval < 0.01 & index_fc,]
sub_markers$avg_FC <- rowMeans(sub_markers[,c("IPSCs_avg_logFC","Day06_avg_logFC","Day15_avg_logFC","Day21_avg_logFC")])
sub_markers2 <- sub_markers[abs(sub_markers$avg_FC) >0.1,]
dim(sub_markers2)
#172
write.table(sub_markers2,"Conserved_all_alt.txt")




# pdf("Skypin_plot.pdf",width=12)
# ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
#     geom_bar(stat="identity", color="black", position=position_dodge() )+
#     geom_text(aes(label=round(pct.exp,0)), vjust=1.6, color="white",
#               position = position_dodge(0.9), size=2.5)+
#     theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")

# ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
#     geom_bar(stat="identity", color="black", position=position_dodge() )+
#     theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")

# dev.off()

# ggplot(data=toplot, aes(x=id, y=pct.exp, fill=features.plot)) +
#     geom_bar(stat="identity", color="black", position=position_dodge() )+
#     geom_text(aes(label=round(pct.exp,0)), vjust=1.6, color="white",
#               position = position_dodge(0.9), size=2.5)+
#     theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")
# DoHeatmap(Controls,unlist_feat,raster=FALSE,group.by = "Timepoints")
# pheatmap::pheatmap(ord_mat,cluster_rows = F,cluster_cols = F,scale = "none",
#  display_numbers = round(ord_mat,2))#annotation_row = as.data.frame(annot_map[,1]),
# dev.off()








# ====================================== PROTEOMICS =========================================
## ----Venn------------------------------------------------------------------------------------
library(VennDiagram)
library(EnhancedVolcano)

#  ================= Figure 7 
# FIG 7A
proteomics_volc_D25 <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/Proteomics_D25_Volcano.csv")
proteomics_volc_D40 <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/Proteomics_D40_Volcano.csv")

proteomics_volc_D25$p.value <- exp(-proteomics_volc_D25$neglog_pvalue )
pdf("GLMER_Volcanos.pdf")
min_vlc <- min(proteomics_volc_D25$log_fc)
max_vlc <- max(proteomics_volc_D25$log_fc)


plot(EnhancedVolcano(proteomics_volc_D25,
                        lab = proteomics_volc_D25$firstGene,
                        x = 'log_fc',
                        y = 'p.value',subtitle = "FC2 (FCcutoff=1)",
                        xlim = c(min_vlc-1,max_vlc +1),FCcutoff = 1)+xlab("Log2 fold change")+ggtitle("GLMM"))

dev.off()






# ================== FIG 7B
proteomics <- read.csv("C:/Users/dimitrios.kyriakis/Desktop/selected.proteins.short.csv")

# PINK1_D25_R1 	PINK1 	D25 	R1 	LFQ.intensity.D25.66.8
# PINK1_D25_R2 	PINK1 	D25 	R2 	LFQ.intensity.D25.66.7
# Ctrl_D25_R1 	Ctrl 	D25 	R1 	LFQ.intensity.D25.C8
# Ctrl_D25_R2 	Ctrl 	D25 	R2 	LFQ.intensity.D25.C7
# PINK1_D40_R1 	PINK1 	D40 	R1 	LFQ.intensity.D40.66.8
# PINK1_D40_R2 	PINK1 	D40 	R2 	LFQ.intensity.D40.66.7
# Ctrl_D40_R1 	Ctrl 	D40 	R1 	LFQ.intensity.D40.C.8
# Ctrl_D40_R2 	Ctrl 	D40 	R2 	LFQ.intensity.D40.C.7


col_keep <- c("Gene.names",
"LFQ.intensity.D25.66.7",
"LFQ.intensity.D25.66.8",       
"LFQ.intensity.D25.C7",      
"LFQ.intensity.D25.C8",        
"LFQ.intensity.D40.66.7",       
"LFQ.intensity.D40.66.8",       
"LFQ.intensity.D40.C.7",       
"LFQ.intensity.D40.C.8")

proteomics_filt <- proteomics[,col_keep]
melted_prot <- melt(proteomics_filt)

melted_prot$Day <- as.vector(melted_prot$variable)
melted_prot$Day[grep("D25",melted_prot$Day)] <- "Day25" 
melted_prot$Day[grep("D40",melted_prot$Day)] <- "Day40" 

# Replicate
melted_prot$Replicate <- as.vector(melted_prot$variable)
melted_prot$Replicate[grep("66.8",melted_prot$Replicate)] <- "R1" 
melted_prot$Replicate[grep("66.7",melted_prot$Replicate)] <- "R2"
melted_prot$Replicate[grep("C8",melted_prot$Replicate)] <- "R1" 
melted_prot$Replicate[grep("C7",melted_prot$Replicate)] <- "R2" 
melted_prot$Replicate[grep("C.8",melted_prot$Replicate)] <- "R1" 
melted_prot$Replicate[grep("C.7",melted_prot$Replicate)] <- "R2" 

# Condition
melted_prot$Condition <- as.vector(melted_prot$variable)
melted_prot$Condition[grep("66.8",melted_prot$Condition)] <- "PINK1" 
melted_prot$Condition[grep("66.7",melted_prot$Condition)] <- "PINK1"
melted_prot$Condition[grep("C8",melted_prot$Condition)] <- "Ctrl" 
melted_prot$Condition[grep("C7",melted_prot$Condition)] <- "Ctrl" 
melted_prot$Condition[grep("C.8",melted_prot$Condition)] <- "Ctrl" 
melted_prot$Condition[grep("C.7",melted_prot$Condition)] <- "Ctrl" 



library(ggpval)
p1 <- ggplot(melted_prot %>% filter(Gene.names=="CSRP2") ,aes(x=Day,y=value,fill=Condition))+
    geom_bar(stat='identity', position=position_dodge())+
    scale_fill_manual(values=c('lightgray','black'))+
    facet_wrap(~Replicate) + xlab("")+ggtitle("CSRP2") +
    theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
    legend.title  = element_blank())      
p2 <- ggplot(melted_prot %>% filter(Gene.names=="VWA5A") ,aes(x=Day,y=value,fill=Condition))+
    geom_bar(stat='identity', position=position_dodge())+
    scale_fill_manual(values=c('lightgray','black'))+
    facet_wrap(~Replicate) + xlab("")+ggtitle("VWA5A") +
    theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5,size =15, face = "bold"),
    legend.title  = element_blank())    

p3 <- ggplot(melted_prot %>% filter(Gene.names=="COMT") ,aes(x=Day,y=value,fill=Condition))+
    geom_bar(stat='identity', position=position_dodge())+
    scale_fill_manual(values=c('lightgray','black'))+
    facet_wrap(~Replicate) + xlab("")+ggtitle("COMT") +
    theme(text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
    legend.title  = element_blank())      
p4 <- ggplot(melted_prot %>% filter(Gene.names=="DLK1") ,aes(x=Day,y=value,fill=Condition))+
    geom_bar(stat='identity', position=position_dodge())+
    scale_fill_manual(values=c('lightgray','black'))+
    facet_wrap(~Replicate) + xlab("")+ggtitle("DLK1") +
    theme(text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),
    legend.title  = element_blank())   

pdf("Fig7B.pdf",width=8,height=12)
ggarrange(plotlist=list(p1,p2,p3,p4),ncol=1)
dev.off()



