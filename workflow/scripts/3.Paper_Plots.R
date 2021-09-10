args = commandArgs(trailingOnly=TRUE)
print(args)

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
library(dplyr)
library(gridExtra)
set.seed(123)
list.of.packages <- c("tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
library(tibble)

project ="IPSCs_pink1"
colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")[c(5,1,2,3,4,9,6,7,8)]
color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6])[c(5,1,2,3,4,9,6,7,8)]
color_clust <- c(brewer.pal(12,"Paired")[-11],"black","gray","magenta4","seagreen4",brewer.pal(9,"Set1")[-6],brewer.pal(8,"Dark2"))
color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

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
input_file <- args[1]
threads <- as.numeric(args[2])
output_group_A <- args[3]
# ---------------------------------------------

Combined <- readRDS(input_file)
Combined$Timepoints <- as.vector(Combined$date)
Combined$Timepoints[grep("IPSC",Combined$date)] <-"iPSCs"
Combined$Timepoints[grep("D06",Combined$date)] <-"Day06"
Combined$Timepoints[grep("D10",Combined$date)] <-"Day10"
Combined$Timepoints[grep("D15",Combined$date)] <-"Day15"
Combined$Timepoints[grep("D21",Combined$date)] <-"Day21"
Combined$Timepoints <- factor(Combined$Timepoints,c("iPSCs","Day06","Day10","Day15","Day21"))

# Controls <- readRDS("../Reviewes/All_controls.rds")
# Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21","D26","D35","D50"))

Controls <- subset(Combined,subset= Treatment !="PINK1")
Controls_AD10 <- subset(Combined,subset= Timepoints !="Day10")
Controls_CD10 <- subset(Controls,subset= Timepoints !="Day10")

DefaultAssay(Combined) <- "RNA"
DefaultAssay(Controls) <- "RNA"
DefaultAssay(Controls_AD10) <- "RNA"
DefaultAssay(Controls_CD10) <- "RNA"

Combined <- ScaleData(Combined,features=rownames(Combined))
Controls <- ScaleData(Controls,features=rownames(Controls))
Controls_AD10 <- ScaleData(Controls_AD10,features=rownames(Controls_AD10))
Controls_CD10 <- ScaleData(Controls_CD10,features=rownames(Controls_CD10))

# ================= Figure 1 B ============================
fig1_markers <- c("MYC","POU5F1",
"PTCH1","FZD7",
"HES1","OTX2",
"SLIT2","LMX1A",
"DCX","DDC")

pdf("result/Paper_Figures/Figure1B.pdf",width= 8)
DefaultAssay(Controls_AD10) <- "RNA"
Controls_AD10$Timepoints <- factor(Controls_AD10$Timepoints)
DoHeatmap(Controls_AD10,fig1_markers,hjust=0.5,angle=0,
    group.by = "Timepoints", raster = FALSE)+
    scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))

Controls_CD10$Timepoints <- factor(Controls_CD10$Timepoints)
DoHeatmap(Controls_CD10,fig1_markers,hjust=0.5,angle=0,
    group.by = "Timepoints", raster = FALSE)+
    scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))
dev.off()
# ---------------------------------------------------------

# ==================== Fig 2 B ============================
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
pdf("result/Paper_Figures/Fig2B.pdf",width= 8)
DefaultAssay(Combined) <- "RNA"
# DoHeatmap(Combined,ipsc_markers,hjust=0.5,angle=0,
#     group.by = "Timepoints", raster = FALSE)+
#     scale_fill_viridis(option="inferno")+ 
#     theme(text = element_text(size = 20))
DoHeatmap(Controls_CD10,ipsc_markers,hjust=0.5,angle=0,group.by = "Timepoints", raster = FALSE)+ scale_fill_viridis(option="inferno")+ 
    theme(text = element_text(size = 20))
dev.off()
# -------------------------------------------------------




# ==================== Fig 3B ============================

TH<-c(1.92222E-05,	
    2.17257E-05,
    4.86911E-05,
    0.000110361,
    0.016417411,
    0.023048918)

TH_se<-c(1.05E-05,
3.08E-06,
1.77E-05,
1.52E-05,
7.11E-04,
9.72E-03)

TH_sd <-c(1.82E-05,
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

print(df)

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


pdf("result/Paper_Figures/Figure_3B_SD_SE.pdf",height=4,width=12)
fig3b
fig3b2
dev.off()
# ---------------------------------------------------------



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
# "ROBO1"	,
# "ROBO2"	,
# "NTN1"	,
# "SLIT1"	,
# "OTX2"	,
# "KCNH6" ,	
"PITX2",
"ASCL1",
"PITX3"	,
"NEUROD1",
"CORIN",
"SLIT1",
"TH",
"NR4A2"
# "CORIN",
# "NEUROD1",
# "TUBB3",
# "NR4A2"	
)

# ===================== Figure 3D, 3E =====================
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

print("OK")
split_list <- list(IPSC =IPSC,Rgl=unique(c(Rgl3,Rgl1)),Prog=unique(c(ProgM,ProgFPs)),NProg=unique(c(NProg,NbM)),DA=unique(c(DA0,DA1,DA2)))
unlist_feat <- unlist(split_list)
colnames(Controls@meta.data)
Controls@meta.data <- Controls@meta.data[,c(1:20)]
Controls$Timepoints <- factor(unlist(lapply(as.vector(Controls$condition),function(x){strsplit(x,"_")[[1]][2]})),levels=c("IPSCs","D06","D10","D15","D21"))
Controls <-AddModuleScore(Controls,features = split_list,name = "DA_s")
colnames(Controls@meta.data)[22:26]<-names(split_list)
scores_meta <-  Controls@meta.data[22:26]
scores_meta_assign <- apply(scores_meta, 1,function(x) {
    if(x[which.max(x)]<0){
        col_n <- "Undefined"
    }else{
        col_n <- colnames(scores_meta)[which.max(x)]
    }
    return(col_n)
})
print("OK1")
Controls$Phase_Dif <- factor(scores_meta_assign,
    levels=c("Undefined","IPSC","Rgl","Prog","NProg","DA"))
table_pl <- rbind(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="IPSCs"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D06"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D06"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D10"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D10"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D15"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D15"])),
table(Controls$Phase_Dif[Controls$Timepoints=="D21"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D21"]))#,
# table(Controls$Phase_Dif[Controls$Timepoints=="D26"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D26"])),
# table(Controls$Phase_Dif[Controls$Timepoints=="D35"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D35"])),
# table(Controls$Phase_Dif[Controls$Timepoints=="D50"]) / sum(table(Controls$Phase_Dif[Controls$Timepoints=="D50"]))
)
rownames(table_pl) <-  c("IPSCs","D06","D10","D15","D21")#,"D26","D35"  ,"D50")
toplot <- reshape2::melt(table_pl)
colnames(toplot) <- c("TimePoint","Group","value")
print("OK3")
pdf("result/Paper_Figures/Figure_3D_3C_3E.pdf",width=12)
ggplot(data=toplot, aes(x=TimePoint, y=value, fill=Group)) +
    geom_bar(stat="identity", color="black", position=position_dodge() )+
    theme_minimal()+xlab("Timepoints")+ylab("% cell percentage")+ 
    theme(axis.text=element_text(size=11),axis.text.x = element_text(angle = 45))
DotPlot(Controls,
    features =c("Undefined","IPSC","Rgl","Prog","NProg","DA"),group.by="Timepoints")+coord_flip()
DoHeatmap(Controls,
    features = mDA2,
    group.by = "Timepoints")+ scale_fill_viridis(option="inferno")
dev.off()
# ----------------------------------------------------------




# ======================= FIGURE 4 ==========================
mDA <- c("TCF12", "ALCAM", "PITX2", "DDC", "ASCL1")
mDA <- mDA[mDA %in% rownames(Combined@assays$RNA@counts)]
Combined@reductions$umap@cell.embeddings[,2] <- Combined@reductions$umap@cell.embeddings[,2]*(-1)
p0 <- DimPlot(Combined,group.by = "condition",cols=color_cond)
results <- ICSWrapper::scatter_gene(object = Combined,features = c("TH","KCNJ6"),ncol =1 ,nrow =2)+ theme_void()
plot5 <- grid.arrange(p0, results, widths=c(2/3, 1/2), ncol=2, nrow=1)
plot5 <- ggarrange(plotlist=list(p0, results), widths=c(2/3, 1/2), ncol=2, nrow=1)
pdf("result/Paper_Figures/Figure4.pdf",width=12,height=10)
plot5
dev.off()
# ---------------------------------------------------------











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











# ================================== MODULE SCORES 
# ================================== MODULE SCORES 







EN1_GAPDH<- c(4.31E-07,4.73E-08,9.81E-08,
            3.12E-07,7.49E-07,2.21E-06,
            3.58E-06,2.90E-06,2.67E-06,
            5.22E-06,3.44E-06,2.76E-06)
EN1_GAPDH_t <- c(0,0,0,
                6,6,6,
                10,10,10,
                21,21,21)

OTX2_GAPDH<- c(0.000717182,	0.000217518,	0.000234211,
                0.007851153,	0.005740741,	0.007124842,
                0.006518804,	0.008503679,	0.006132597,
                0.014738903,	0.018468354,	0.012615979,
                0.009584121,	0.008388278,	0.010460772)
OTX2_GAPDH_t <- c(0,0,0,
                6,6,6,
                10,10,10,
                15,15,15,
                21,21,21)

EN1_mat <- as.data.frame(cbind(EN1_GAPDH,EN1_GAPDH_t))
colnames(EN1_mat)<-c("Expression", "EN1_t")
OTX2_mat <- as.data.frame(cbind(OTX2_GAPDH,OTX2_GAPDH_t))
colnames(OTX2_mat)<-c("Expression", "OTX2_t")

stderr <- function(x, na.rm=TRUE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= stderr(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

df1 <- data_summary(EN1_mat, varname="Expression", 
                    groupnames=c("EN1_GAPDH_t"))

df2 <- data_summary(OTX2_mat, varname="Expression", 
                    groupnames=c("OTX2_t"))
# Convert dose to a factor variable

p1 <- ggplot(df1, aes(EN1_GAPDH_t,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=2)+
    geom_point(size=1,aes(EN1_GAPDH_t,Expression),data=EN1_mat)+
    theme_minimal()+
    ggtitle("EN1/GAPDH")+ 
    theme(text = element_text(size = 15))+
    theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))+xlab("TimePoint")

p2 <- ggplot(df2, aes(OTX2_t,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-sd, ymax=Expression+sd), width=2)+
    geom_point(size=1,aes(OTX2_t,Expression),data=OTX2_mat)+
    theme_minimal()+
    ggtitle("OTX2/GAPDH")+ 
    theme(text = element_text(size = 15))+
    theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))+xlab("TimePoint")
p3 <- ggplot(df1, aes(EN1_GAPDH_t,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=2)+
    geom_point(size=1,aes(EN1_GAPDH_t,Expression),data=EN1_mat)+
    theme_minimal()+
    ggtitle("EN1/GAPDH")+ 
    theme(text = element_text(size = 15))+
    theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))+xlab("TimePoint")

p4 <- ggplot(df2, aes(OTX2_t,Expression)) +
    geom_line() +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Expression-se, ymax=Expression+se), width=2)+
    geom_point(size=1,aes(OTX2_t,Expression),data=OTX2_mat)+
    theme_minimal()+
    ggtitle("OTX2/GAPDH")+ 
    theme(text = element_text(size = 15))+
    theme(plot.title = element_text(size=11,face="bold",hjust = 0.5))+xlab("TimePoint")

EN1_OTX2_p12 <- ggarrange(plotlist = list(p1,p2),ncol=2)
EN1_OTX2_p34 <- ggarrange(plotlist = list(p3,p4),ncol=2)
pdf("result/Paper_Figures/QPCR_SD_SE_EN1_OTX2.pdf",height=4,width=8)
EN1_OTX2_p12
EN1_OTX2_p34
dev.off()


