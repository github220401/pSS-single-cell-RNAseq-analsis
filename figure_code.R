library(LIANLAB)
library(Seurat)
library(ggplot2)
library(scRepertoire)
library(GSEABase)
library(pheatmap)
library(limma)
library(biomaRt)
library("ggsci")
library("ggplot2")
library("gridExtra")
color_used <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
color <- c("#DFCEE6","#DED179","#CA9687","#CA9687","#A6C5D8","#D07F79","#77A2EB","#76B1F6","#CA81C1")
get_distinct_hues <- function(ncolor,s=0.5,v=0.95,seed=40) {
  golden_ratio_conjugate <- 0.618033988749895
  set.seed(seed)
  h <- runif(1)
  H <- vector("numeric",ncolor)
  for(i in seq_len(ncolor)) {
    h <- (h + golden_ratio_conjugate) %% 1
    H[i] <- h
  }
  hsv(H,s=s,v=v)
}

col=get_distinct_hues(2500)

#figure 1
my_all_cells <- readRDS(file = 'my_all_cells.rds')

Idents(my_all_cells) <- 'rename'
my_all_cells@active.ident <- factor(my_all_cells@active.ident,levels = c(
  'ANO6_Acinar_cell','AQP5_Acinar_cell','SOX10_Acinar_cell','Basal_Duct_1','Basal_Duct_2',
  'KIT_duct_cell','Naive_B','Memory_B','CD11c_Memory_B','PlasmaBlast','Plasma','NK',
  'MAIT','gd_T','CD8_gd_T','CD4T_naive','CD4TCM','CD4Treg','GZMK_CD4TCM','GZMB_CD4TEM',
  'CD4TRM','CD8T_naive','Proliferating_CD8T','GZMK_CD8TCM','GZMB_CD8TEM','CD69_CD8TRM','CD103_CD8TRM',
  'CD14_CD11b_Mo','CD14_CD163_Mo','CD14_Mo','CX3CR1_CD14_Mo','CX3CR1_CD16_Mo','Neutrophil',
  'Mast_cell','cDC_2'
))
DimPlot(my_all_cells,col = colorful[[10]])
DefaultAssay(my_all_cells) <- 'RNA'
myfindmarkers(my_all_cells,filename = 'fig1c',colors = colorful[[10]])

#figure 2
my_T_cells <- readRDS(file = 'my_T_cells.rds')

#figure 2A
combined <- expression2List(my_T_cells,split.by = 'group')
pdf(file = 'figure 2A.pdf',width = 5,height = 5)
clonalHomeostasis(combined, cloneCall = "aa")
dev.off()
write.csv(clonalHomeostasis(combined, cloneCall = "aa",exportTable = T),file = 'figure 2A.csv')

#figure 2B
a <- as.data.frame(table(ss_all_T$CTaa))
b <- a[a$Freq>1,]

hc <- subset(my_T_cells@meta.data,group%in%'hc')
hc.clonal <- subset(hc,CTaa %in% b$Var1)
nrow(hc.clonal) / nrow(hc)
nrow(hc.clonal[hc.clonal$samples%in%'HC5',]) / nrow(hc[hc$samples%in%'HC5',])
nrow(hc.clonal[hc.clonal$samples%in%'HC6',]) / nrow(hc[hc$samples%in%'HC6',])
nrow(hc.clonal[hc.clonal$samples%in%'HC7',]) / nrow(hc[hc$samples%in%'HC7',])
nrow(hc.clonal[hc.clonal$samples%in%'HC8',]) / nrow(hc[hc$samples%in%'HC8',])

ss_pbmc <- subset(my_T_cells@meta.data,group%in%'SS_PBMC')
ss_pbmc.clonal <- subset(ss_pbmc,CTaa %in% b$Var1)
table(ss_pbmc$samples)
nrow(ss_pbmc.clonal) / nrow(ss_pbmc)
nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_0',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_0',])
nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_1',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_1',])
nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_2',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_2',])

ss_sg <- subset(my_T_cells@meta.data,group%in%'SS_SG')
ss_sg.clonal <- subset(ss_sg,CTaa %in% b$Var1)
table(ss_sg.clonal$samples)
nrow(ss_sg.clonal) / nrow(ss_sg)
nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_0',]) / nrow(ss_sg[ss_sg$samples%in%'SG_0',])
nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_1',]) / nrow(ss_sg[ss_sg$samples%in%'SG_1',])
nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_2',]) / nrow(ss_sg[ss_sg$samples%in%'SG_2',])

data <- as.data.frame(matrix(nrow = 10,ncol = 5))
colnames(data) <- c('samples','group','clonal.percent','clonal.percent.average','clonal.percent.sd')
data$samples <- c('hc5','hc6','hc7','hc8','s0_pbmc','s1_pbmc','s2_pbmc','s0_sg','s1_sg','s2_sg')
data$group <- rep(c('hc','pbmc','sg'),time = c(4,3,3))
data$clonal.percent[data$samples%in%'hc5'] <- nrow(hc.clonal[hc.clonal$samples%in%'HC5',]) / nrow(hc[hc$samples%in%'HC5',])
data$clonal.percent[data$samples%in%'hc6'] <- nrow(hc.clonal[hc.clonal$samples%in%'HC6',]) / nrow(hc[hc$samples%in%'HC6',])
data$clonal.percent[data$samples%in%'hc7'] <- nrow(hc.clonal[hc.clonal$samples%in%'HC7',]) / nrow(hc[hc$samples%in%'HC7',])
data$clonal.percent[data$samples%in%'hc8'] <- nrow(hc.clonal[hc.clonal$samples%in%'HC8',]) / nrow(hc[hc$samples%in%'HC8',])
data$clonal.percent[data$samples%in%'s0_pbmc'] <- nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_0',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_0',])
data$clonal.percent[data$samples%in%'s1_pbmc'] <- nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_1',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_1',])
data$clonal.percent[data$samples%in%'s2_pbmc'] <- nrow(ss_pbmc.clonal[ss_pbmc.clonal$samples%in%'PBMC_2',]) / nrow(ss_pbmc[ss_pbmc$samples%in%'PBMC_2',])
data$clonal.percent[data$samples%in%'s0_sg'] <- nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_0',]) / nrow(ss_sg[ss_sg$samples%in%'SG_0',])
data$clonal.percent[data$samples%in%'s1_sg'] <- nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_1',]) / nrow(ss_sg[ss_sg$samples%in%'SG_1',])
data$clonal.percent[data$samples%in%'s2_sg'] <- nrow(ss_sg.clonal[ss_sg.clonal$samples%in%'SG_2',]) / nrow(ss_sg[ss_sg$samples%in%'SG_2',])
data$clonal.percent.average[data$group%in%'hc'] <- mean(data$clonal.percent[data$group%in%'hc'])
data$clonal.percent.average[data$group%in%'pbmc'] <- mean(data$clonal.percent[data$group%in%'pbmc'])
data$clonal.percent.average[data$group%in%'sg'] <- mean(data$clonal.percent[data$group%in%'sg'])
data$clonal.percent.sd[data$group%in%'hc'] <- sd(data$clonal.percent[data$group%in%'hc'])
data$clonal.percent.sd[data$group%in%'pbmc'] <- sd(data$clonal.percent[data$group%in%'pbmc'])
data$clonal.percent.sd[data$group%in%'sg'] <- sd(data$clonal.percent[data$group%in%'sg'])

pdf(file = 'figure 2B.pdf',width = 5.13,height = 3.48)
ggplot(data = data) + geom_boxplot(aes(x = group, y = clonal.percent,
                                              fill = factor(group)))+
  scale_fill_manual(values=colorful[[1]])+
  scale_y_continuous(labels = scales::percent)+
  theme_grey(base_size = 18)+
  labs(x='cluster',y='clonal cell percent in each cluster')+
  theme(axis.text= element_text(colour = "black",face='bold'),
        axis.title = element_text(face = 'bold'),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size = 10))+
  geom_point(data = data,aes(group,clonal.percent),size=2,pch=19)
dev.off()

write.csv(data,file = 'figure 2B.csv')

#figure 2C
combined <- expression2List(subset(my_T_cells,group%in%c('SS_PBMC','SS_SG')),split.by = 'cluster')
pdf(file = 'figure 2 C.pdf',width = 4.5,height = 5)
clonalHomeostasis(combined, cloneCall = "aa")
dev.off()
write.csv(clonalHomeostasis(combined, cloneCall = "aa",exportTable = T),file = 'figure 2C.csv')

#figure 2D
meta <- my_T_cells@meta.data
meta$patients <- 'other'
meta$patients[meta$samples%in%c('HC5','HC6','HC7','HC8')] <- 'hc'
meta$patients[meta$samples%in%c('PBMC_0','PBMC_1','PBMC_2','SG_0','SG_1','SG_2')] <- 'ss'
table(meta$patients)
meta$Cell <- 'other'
meta$Cell[meta$rename%in%c('CD4T_naive','CD4TCM','CD4TEM','CD4Treg','CD4TRM','GZMB_CD4TEM','GZMK_CD4TCM')] <- 'CD4'
meta$Cell[meta$rename%in%c('CD103_CD8TRM','CD69_CD8TRM','CD8_gd_T','CD8T_naive','GZMB_CD8TEM','GZMK_CD8TCM','Proliferating_CD8T')] <- 'CD8'
table(meta$Cell)
meta$patients_Cell <- paste0(meta$patients,'_',meta$Cell)
table(meta$patients_Cell)

hc.4 <- subset(meta,patients_Cell%in%'hc_CD4')
hc.4.clonal <- subset(hc.4,CTaa %in% b$Var1)
nrow(hc.4.clonal) / nrow(hc.4)
nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC5',]) / nrow(hc.4[hc.4$samples%in%'HC5',])
nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC6',]) / nrow(hc.4[hc.4$samples%in%'HC6',])
nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC7',]) / nrow(hc.4[hc.4$samples%in%'HC7',])
nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC8',]) / nrow(hc.4[hc.4$samples%in%'HC8',])

ss.4 <- subset(meta,patients_Cell%in%'ss_CD4')
ss.4.clonal <- subset(ss.4,CTaa %in% b$Var1)
nrow(ss.4.clonal) / nrow(ss.4)
nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_0','SG_0'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_0','SG_0'),])
nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_1','SG_1'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_1','SG_1'),])
nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_2','SG_2'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_2','SG_2'),])

hc.8 <- subset(meta,patients_Cell%in%'hc_CD8')
hc.8.clonal <- subset(hc.8,CTaa %in% b$Var1)
nrow(hc.8.clonal) / nrow(hc.8)
nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC5',]) / nrow(hc.8[hc.8$samples%in%'HC5',])
nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC6',]) / nrow(hc.8[hc.8$samples%in%'HC6',])
nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC7',]) / nrow(hc.8[hc.8$samples%in%'HC7',])
nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC8',]) / nrow(hc.8[hc.8$samples%in%'HC8',])

ss.8 <- subset(meta,patients_Cell%in%'ss_CD8')
ss.8.clonal <- subset(ss.8,CTaa %in% b$Var1)
nrow(ss.8.clonal) / nrow(ss.8)
nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_0','SG_0'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_0','SG_0'),])
nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_1','SG_1'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_1','SG_1'),])
nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_2','SG_2'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_2','SG_2'),])

data <- as.data.frame(matrix(nrow = 14,ncol = 6))
colnames(data) <- c('samples','cluster','group','clonal.percent','clonal.percent.average','clonal.percent.sd')
data$samples <- rep(c('hc5','hc6','hc7','hc8','s0','s1','s2'),time = 2)
data$cluster <- rep(c('CD4','CD8'),each = 7)
data$group <- rep(c('hc.4','ss.4','hc.8','ss.8'),time = c(4,3,4,3))
data$clonal.percent[data$samples%in%'hc5' & data$cluster%in%'CD4'] <- nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC5',]) / nrow(hc.4[hc.4$samples%in%'HC5',])
data$clonal.percent[data$samples%in%'hc6' & data$cluster%in%'CD4'] <- nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC6',]) / nrow(hc.4[hc.4$samples%in%'HC6',])
data$clonal.percent[data$samples%in%'hc7' & data$cluster%in%'CD4'] <- nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC7',]) / nrow(hc.4[hc.4$samples%in%'HC7',])
data$clonal.percent[data$samples%in%'hc8' & data$cluster%in%'CD4'] <- nrow(hc.4.clonal[hc.4.clonal$samples%in%'HC8',]) / nrow(hc.4[hc.4$samples%in%'HC8',])
data$clonal.percent[data$samples%in%'s0' & data$cluster%in%'CD4'] <- nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_0','SG_0'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_0','SG_0'),])
data$clonal.percent[data$samples%in%'s1' & data$cluster%in%'CD4'] <- nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_1','SG_1'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_1','SG_1'),])
data$clonal.percent[data$samples%in%'s2' & data$cluster%in%'CD4'] <- nrow(ss.4.clonal[ss.4.clonal$samples%in%c('PBMC_2','SG_2'),]) / nrow(ss.4[ss.4$samples%in%c('PBMC_2','SG_2'),])
data$clonal.percent[data$samples%in%'hc5' & data$cluster%in%'CD8'] <- nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC5',]) / nrow(hc.8[hc.8$samples%in%'HC5',])
data$clonal.percent[data$samples%in%'hc6' & data$cluster%in%'CD8'] <- nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC6',]) / nrow(hc.8[hc.8$samples%in%'HC6',])
data$clonal.percent[data$samples%in%'hc7' & data$cluster%in%'CD8'] <- nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC7',]) / nrow(hc.8[hc.8$samples%in%'HC7',])
data$clonal.percent[data$samples%in%'hc8' & data$cluster%in%'CD8'] <- nrow(hc.8.clonal[hc.8.clonal$samples%in%'HC8',]) / nrow(hc.8[hc.8$samples%in%'HC8',])
data$clonal.percent[data$samples%in%'s0' & data$cluster%in%'CD8'] <- nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_0','SG_0'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_0','SG_0'),])
data$clonal.percent[data$samples%in%'s1' & data$cluster%in%'CD8'] <- nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_1','SG_1'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_1','SG_1'),])
data$clonal.percent[data$samples%in%'s2' & data$cluster%in%'CD8'] <- nrow(ss.8.clonal[ss.8.clonal$samples%in%c('PBMC_2','SG_2'),]) / nrow(ss.8[ss.8$samples%in%c('PBMC_2','SG_2'),])

data$clonal.percent.average[data$group%in%'hc.4'] <- mean(data$clonal.percent[data$group%in%'hc.4'])
data$clonal.percent.average[data$group%in%'hc.8'] <- mean(data$clonal.percent[data$group%in%'hc.8'])
data$clonal.percent.average[data$group%in%'ss.4'] <- mean(data$clonal.percent[data$group%in%'ss.4'])
data$clonal.percent.average[data$group%in%'ss.8'] <- mean(data$clonal.percent[data$group%in%'ss.8'])
data$clonal.percent.sd[data$group%in%'hc.4'] <- mean(data$clonal.percent[data$group%in%'hc.4'])
data$clonal.percent.sd[data$group%in%'hc.8'] <- mean(data$clonal.percent[data$group%in%'hc.8'])
data$clonal.percent.sd[data$group%in%'ss.4'] <- mean(data$clonal.percent[data$group%in%'ss.4'])
data$clonal.percent.sd[data$group%in%'ss.8'] <- mean(data$clonal.percent[data$group%in%'ss.8'])

pdf(file = 'figure 2D.pdf',width = 6.5,height = 5)
ggplot(data = data) + geom_boxplot(aes(x = group, y = clonal.percent,
                                       fill = factor(group)))+
  scale_fill_manual(values=colorful[[1]][5:10])+
  scale_y_continuous(labels = scales::percent)+
  theme_grey(base_size = 18)+
  labs(x='cluster',y='clonal cell percent in each cluster')+
  theme(axis.text= element_text(colour = "black",face='bold'),
        axis.title = element_text(face = 'bold'),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size = 10))+
  geom_point(data = data,aes(group,clonal.percent),size=2,pch=19)
dev.off()

write.csv(data,file = 'figure 2D.csv')


#figure 2E
Idents(my_T_cells) <- 'cluster'
pdf(file = 'figure 2E.pdf',width = 6.5,height = 5)
DimPlot(my_T_cells,cols = color_used[4:5])
dev.off()

slot(my_T_cells, "meta.data")$cloneType <- factor(slot(my_T_cells, "meta.data")$cloneType,
                                                levels = c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)",
                                                           "Medium (0.001 < X <= 0.01)", "Small (1e-04 < X <= 0.001)",  NA))
pdf("figure 2E.pdf",width = 7.5,height = 5)
DimPlot(my_T_cells, group.by = "cloneType",label = F) +
  scale_color_manual(values = c("coral1", "mediumpurple1","deepskyblue","burlywood3"), na.value="grey")
dev.off()

#figure 2F
my_CD8_cells <- readRDS(file = 'my_CD8_cells.rds')
combined2 <- expression2List(my_CD8_cells,split.by = 'samples')
table3 = clonalDiversity(combined2, cloneCall = "aa", group = "group",exportTable = T)

pdf(file = 'figure 2F Chao.pdf',width = 3.5,height = 5)
ggplot(table3, aes(x = group,y = Chao,fill = group)) +
  geom_bar(width = 0.9,stat="identity") +
  scale_fill_manual(values=c("#86B67E","#CA81C1","#D36E46"))+
  theme_classic() +
  theme(panel.grid = element_blank())
dev.off()

pdf(file = 'figure 2F ACE.pdf',width = 3.5,height = 5)
ggplot(table3, aes(x = group,y = ACE,fill = group)) +
  geom_bar(width = 0.9,stat="identity") +
  scale_fill_manual(values=c("#86B67E","#CA81C1","#D36E46"))+
  theme_classic() +
  theme(panel.grid = element_blank())
dev.off()

write.csv(clonalDiversity(combined2, cloneCall = "aa", group = "group",exportTable = T),file = 'figure 2F.csv')

#figure 2G
pdf("figure 2G.pdf",width = 6.5,height = 5)
clonalOverlap(combined2, cloneCall = "aa", method = "morisita")
dev.off()

#figure 2H
combined2 <- expression2List(my_CD8_cells,split.by = 'group')

p1 <- compareClonotypes(combined2, numbers = 2000, samples = c("SS_PBMC",'SS_SG'),
                        cloneCall="aa", graph = "alluvial")
p1 + scale_fill_manual(values=col) + theme(legend.position = 'none')

#figure 3
#figure 3A
Idents(my_T_cells) <- "rename"
DefaultAssay(my_T_cells) <- 'RNA'
T.cell = RenameIdents(subset(my_T_cells,rename != 'CD8_gd_T'),'GZMK_CD8TCM'='PB_CD8TM','GZMB_CD8TEM'='PB_CD8TM',"GZMB_CD4TEM"="PB_CD4TM","GZMK_CD4TCM"="PB_CD4TM","CD4TCM"="PB_CD4TM",'CD4TEM' = 'PB_CD4TM')
T.cell@active.ident = factor(T.cell@active.ident, levels = c("CD4T_naive",'CD8T_naive',"CD4Treg",'Proliferating_CD8T',"PB_CD4TM",'PB_CD8TM',
                                                             "CD4TRM",'CD69_CD8TRM','CD103_CD8TRM'))
Trm_up_gene <- readxl::read_xlsx("trm_up_1.xlsx")
gene <- as.list(Trm_up_gene)
View(gene)

T.cell.1 <- AddModuleScore(
  object = T.cell,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'CD_Features'
)

color <- c("#DFCEE6","#DED179","#CA9687","#86B67E","#B05F5F","#A6C5D8","#D36E46","#76B1F6","#CA81C1")
VlnPlot(T.cell.1,features = 'CD_Features1',cols = color,pt.size = 0.1)

#figure 3B
my_ss_cells <- readRDS(file = 'my_ss_cells.rds')
sg <- subset(my_ss_cells,group%in%'SS_SG')
data.input = sg@assays$RNA@data
meta = sg@meta.data
cell.use = rownames(meta)[!meta$rename%in%c('CD11c_Memory_B','Neutrophil','Naive_B','CD4T_naive','CD4TCM','CD4Treg','CD8T_naive','GZMB_CD4TEM','GZMB_CD8TEM','GZMK_CD4TCM','GZMK_CD8TCM','MAIT','NK','CX3CR1_CD16_Mo','cDC_2','Proliferating_CD8T')]

data.input = data.input[, cell.use]
meta = meta[cell.use, ]
unique(meta$rename)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "rename")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "rename")
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
pathways.show <- c("IFN-II")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
pdf(file = 'figure 3B down.pdf',width = 16,height = 10)
ht1 + ht2
dev.off()

pdf(file = 'figure 3B up.pdf',width = 12,height = 10)
plotGeneExpression(cellchat,signaling = pathways.show,type = 'dot')
dev.off()

#figure 3C
GZMK.sg <- subset(subset(my_CD8_cells,rename%in%c('CD103_CD8TRM','CD69_CD8TRM')),group%in%'SS_SG')

create_supercell_matrix=function(seurat_object,ncells=50){
  dm=as.data.frame(seurat_object@assays$RNA@data)
  split.clusters=split(as.data.frame(t(dm)),seurat_object@active.ident)
  supercell_matrix=data.frame()
  for (i in 1:length(split.clusters)) {
    split.matrix=split.clusters[[i]]
    along=1:nrow(split.matrix)
    supercells=aggregate(split.matrix,list(ceiling(along/ncells)),mean)
    supercells$clusters=rep(names(split.clusters)[[i]],nrow(supercells))
    supercell_matrix=rbind(supercell_matrix,supercells)
  }
  cluster=supercell_matrix$clusters
  supercell_matrix=subset(supercell_matrix,select=-c(Group.1,clusters))
  supercell_matrix=as.data.frame(t(supercell_matrix))
  return(list(supercell_matrix=supercell_matrix,cluster=cluster))
}

Idents(GZMK.sg) <- 'rename'
total_matrix=create_supercell_matrix(seurat_object = GZMK.sg,ncells = 50)
supercell_cluster=total_matrix$cluster
total_matrix=total_matrix$supercell_matrix

group1='CD103_CD8TRM'
group_cell=supercell_cluster
group_cell=ifelse(supercell_cluster=='CD103_CD8TRM',"group1","group2")
group_cell=factor(group_cell,levels = c("group1","group2"))
design <- model.matrix(~ factor(group_cell))
vsname=paste0("group2","vs","group1")
colnames(design) <- c("group1",vsname)
fit <- lmFit(total_matrix, design)
fit <- eBayes(fit)
allGeneSets <- topTable(fit, coef=vsname,number=Inf)
gene=rownames(allGeneSets)
df=allGeneSets
genelist=df$logFC
names(genelist)=rownames(df)
genelist=sort(genelist,decreasing = T)

gmt=read.gmt('mygeneset.gmt')
gsea=GSEA(geneList = genelist,TERM2GENE = gmt,minGSSize = 6,eps = 0,pvalueCutoff = 1)
pdf(file = 'figure 3C-1.pdf',width = 6.5,height = 5)
gseaplot2(gsea,geneSetID = 'CD8_Cytotoxic')
dev.off()
pdf(file = 'figure 3C-2.pdf',width = 6.5,height = 5)
gseaplot2(gsea,geneSetID = 'REACTOME_INTERFERON_GAMMA_SIGNALING')
dev.off()

#figure 3F
my_ss_CD8 <- readRDS(file = 'my_ss_CD8.rds')
b=FindMarkers(SS.CD8.integrated,ident.1 = "CD69_CD8TRM",ident.2 = "CD103_CD8TRM",logfc.threshold = 0.25)
b=b[order(b$avg_log2FC,decreasing = F),]
myvolcano(b,gene.plot = 100,width = 6.5,height = 5,text.repel = F,pt.size = 2,max.overlaps=100,filename = 'figure 3F')
write.csv(b,file = 'figure 3F.csv')

#figure 5
#figure 5A
Idents(my_ss_CD8) <- 'rename'
pdf("figure 5A.pdf",width = 6.5,height = 5)
DimPlot(SS.CD8.integrated,cols = c("#E41A1C",'#377EB8',"#4DAF4A","#FF7F00","#984EA3", "#FFFF33"))
dev.off()

#figure 5B
Idents(my_ss_CD8) <- 'tissue'
pdf("figure 5B.pdf",width = 6.5,height = 5)
DimPlot(my_ss_CD8,cols = color_used[10:11])
dev.off()

#figure 5C
my_ss_CD8 <- RenameIdents(my_ss_CD8,'CD103_CD8TRM'='TRM','CD69_CD8TRM'='TRM')
a <- subset(my_ss_CD8,rename!='Proliferating_CD8T')
a@active.ident = factor(a@active.ident, levels =c("CD8T_naive","GZMK_CD8TCM","GZMB_CD8TEM","TRM"))

Trm_up_gene <- readxl::read_xlsx("trm_up_1.xlsx")
gene <- as.list(Trm_down_gene)
Trm_up_score <- AddModuleScore(
  object = a,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'CD_Features'
)
VlnPlot(Trm_up_score,features = 'CD_Features1', cols = color)

Trm_down_gene <- readxl::read_xlsx("trm_down_1.xlsx")
gene <- as.list(Trm_down_gene)
Trm_down_score <- AddModuleScore(
  object = a,
  features = gene,
  ctrl = 100, #默认值是100
  name = 'CD_Features'
)
VlnPlot(Trm_down_score,features = 'CD_Features1', cols = color)

#figure 5D
naive <- subset(my_ss_CD8@meta.data,rename%in%'CD8T_naive')
GZMK <- subset(my_ss_CD8@meta.data,rename%in%'GZMK_CD8TCM')
GZMB <- subset(my_ss_CD8@meta.data,rename%in%'GZMB_CD8TEM')
CD103 <- subset(my_ss_CD8@meta.data,rename%in%'CD103_CD8TRM')
CD69 <- subset(my_ss_CD8@meta.data,rename%in%'CD69_CD8TRM')

data.clonal.percent <- as.data.frame(matrix(nrow = 15,ncol = 4))
colnames(data.clonal.percent) <- c('cluster','sample','sample.cluster','clonalcell.percent')
data.clonal.percent$cluster <- rep(c('naive','GZMK','GZMB','CD69','CD103'),time = 3)
data.clonal.percent$sample <- rep(c('s0','s1','s2'),each = 5)
data.clonal.percent$sample.cluster <- paste0(data.clonal.percent$sample,'_',data.clonal.percent$cluster)

for (i in unique(data.clonal.percent$sample)) {
  naive.sample <- subset(naive,patient%in%i)
  GZMK.sample <- subset(GZMK,patient%in%i)
  GZMB.sample <- subset(GZMB,patient%in%i)
  CD69.sample <- subset(CD69,patient%in%i)
  CD103.sample <- subset(CD103,patient%in%i)
  #Proliferating.sample <- subset(Proliferating,patient%in%i)

  data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','naive')] <- length(rownames(naive.sample)[!is.na(naive.sample$CTaa)]) / nrow(naive.sample)
  data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMK')] <- length(rownames(GZMK.sample)[!is.na(GZMK.sample$CTaa)]) / nrow(GZMK.sample)
  data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMB')] <- length(rownames(GZMB.sample)[!is.na(GZMB.sample$CTaa)]) / nrow(GZMB.sample)
  data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','CD69')] <- length(rownames(CD69.sample)[!is.na(CD69.sample$CTaa)]) / nrow(CD69.sample)
  data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','CD103')] <- length(rownames(CD103.sample)[!is.na(CD103.sample$CTaa)]) / nrow(CD103.sample)
  #data.clonal.percent$clonalcell.percent[data.clonal.percent$sample.cluster%in%paste0(i,'_','Proliferating')] <- length(rownames(Proliferating.sample)[!is.na(Proliferating.sample$CTaa)]) / nrow(Proliferating.sample)

}

data.clonal.percent$average <- 0
data.clonal.percent$average[data.clonal.percent$cluster%in%'naive'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$average[data.clonal.percent$cluster%in%'GZMK'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$average[data.clonal.percent$cluster%in%'GZMB'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])
data.clonal.percent$average[data.clonal.percent$cluster%in%'CD69'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD69']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD69'])
data.clonal.percent$average[data.clonal.percent$cluster%in%'CD103'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD103']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD103'])
#data.clonal.percent$average[data.clonal.percent$cluster%in%'Proliferating'] <- sum(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'Proliferating']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'Proliferating'])
data.clonal.percent$sd <- 0
data.clonal.percent$sd[data.clonal.percent$cluster%in%'naive'] <- sd(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$sd[data.clonal.percent$cluster%in%'GZMK'] <- sd(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$sd[data.clonal.percent$cluster%in%'GZMB'] <- sd(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])
data.clonal.percent$sd[data.clonal.percent$cluster%in%'CD69'] <- sd(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD69'])
data.clonal.percent$sd[data.clonal.percent$cluster%in%'CD103'] <- sd(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'CD103'])

naive.na = subset(naive,CTaa!='NA')
GZMK.na = subset(GZMK,CTaa!='NA')
GZMB.na = subset(GZMB,CTaa!='NA')

data.clonal.percent$share.percent.103 <- 0
for (i in unique(data.clonal.percent$sample)) {
  naive.na.sample <- subset(naive.na,patient%in%i)
  GZMK.na.sample <- subset(GZMK.na,patient%in%i)
  GZMB.na.sample <- subset(GZMB.na,patient%in%i)
  #Proliferating.na.sample <- subset(Proliferating.na,patient%in%i)

  data.clonal.percent$share.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','naive')] <- nrow(naive.na.sample[naive.na.sample$CTaa%in%CD103$CTaa,]) / nrow(naive.na.sample)
  data.clonal.percent$share.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMK')] <- nrow(GZMK.na.sample[GZMK.na.sample$CTaa%in%CD103$CTaa,]) / nrow(GZMK.na.sample)
  data.clonal.percent$share.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMB')] <- nrow(GZMB.na.sample[GZMB.na.sample$CTaa%in%CD103$CTaa,]) / nrow(GZMB.na.sample)
  #data.clonal.percent$share.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','Proliferating')] <- nrow(Proliferating.na.sample[Proliferating.na.sample$CTaa%in%CD103$CTaa,]) / nrow(Proliferating.na.sample)
}

data.clonal.percent$share.percent.69 <- 0
for (i in unique(data.clonal.percent$sample)) {
  naive.na.sample <- subset(naive.na,patient%in%i)
  GZMK.na.sample <- subset(GZMK.na,patient%in%i)
  GZMB.na.sample <- subset(GZMB.na,patient%in%i)
  #Proliferating.na.sample <- subset(Proliferating.na,patient%in%i)

  data.clonal.percent$share.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','naive')] <- nrow(naive.na.sample[naive.na.sample$CTaa%in%CD69$CTaa,]) / nrow(naive.na.sample)
  data.clonal.percent$share.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMK')] <- nrow(GZMK.na.sample[GZMK.na.sample$CTaa%in%CD69$CTaa,]) / nrow(GZMK.na.sample)
  data.clonal.percent$share.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMB')] <- nrow(GZMB.na.sample[GZMB.na.sample$CTaa%in%CD69$CTaa,]) / nrow(GZMB.na.sample)
  #data.clonal.percent$share.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','Proliferating')] <- nrow(Proliferating.na.sample[Proliferating.na.sample$CTaa%in%CD69$CTaa,]) / nrow(Proliferating.na.sample)
}

data.clonal.percent$average.share.percent.103 <- 0
data.clonal.percent$average.share.percent.103[data.clonal.percent$cluster%in%'naive'] <- sum(data.clonal.percent$share.percent.103[data.clonal.percent$cluster%in%'naive']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$average.share.percent.103[data.clonal.percent$cluster%in%'GZMK'] <- sum(data.clonal.percent$share.percent.103[data.clonal.percent$cluster%in%'GZMK']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$average.share.percent.103[data.clonal.percent$cluster%in%'GZMB'] <- sum(data.clonal.percent$share.percent.103[data.clonal.percent$cluster%in%'GZMB']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])

data.clonal.percent$average.share.percent.69 <- 0
data.clonal.percent$average.share.percent.69[data.clonal.percent$cluster%in%'naive'] <- sum(data.clonal.percent$share.percent.69[data.clonal.percent$cluster%in%'naive']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$average.share.percent.69[data.clonal.percent$cluster%in%'GZMK'] <- sum(data.clonal.percent$share.percent.69[data.clonal.percent$cluster%in%'GZMK']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$average.share.percent.69[data.clonal.percent$cluster%in%'GZMB'] <- sum(data.clonal.percent$share.percent.69[data.clonal.percent$cluster%in%'GZMB']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])

CD103.na = subset(CD103,CTaa!='NA')
CD69.na = subset(CD69,CTaa!='NA')

data.clonal.percent$share.kind.percent.103 <- 0
for (i in unique(data.clonal.percent$sample)) {
  naive.na.sample <- subset(naive.na,patient%in%i)
  GZMK.na.sample <- subset(GZMK.na,patient%in%i)
  GZMB.na.sample <- subset(GZMB.na,patient%in%i)
  #Proliferating.na.sample <- subset(Proliferating.na,patient%in%i)
  CD103.na.sample <- subset(CD103.na,patient%in%i)

  data.clonal.percent$share.kind.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','naive')] <- length(intersect(unique(naive.na.sample$CTaa),unique(CD103.na.sample$CTaa))) / length(unique(CD103.na.sample$CTaa))
  data.clonal.percent$share.kind.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMK')] <- length(intersect(unique(GZMK.na.sample$CTaa),unique(CD103.na.sample$CTaa))) / length(unique(CD103.na.sample$CTaa))
  data.clonal.percent$share.kind.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMB')] <- length(intersect(unique(GZMB.na.sample$CTaa),unique(CD103.na.sample$CTaa))) / length(unique(CD103.na.sample$CTaa))
  #data.clonal.percent$share.kind.percent.103[data.clonal.percent$sample.cluster%in%paste0(i,'_','Proliferating')] <- length(intersect(unique(Proliferating.na.sample$CTaa),unique(CD103.na.sample$CTaa))) / length(unique(CD103.na.sample$CTaa))
}

data.clonal.percent$share.kind.percent.69 <- 0
for (i in unique(data.clonal.percent$sample)) {
  naive.na.sample <- subset(naive.na,patient%in%i)
  GZMK.na.sample <- subset(GZMK.na,patient%in%i)
  GZMB.na.sample <- subset(GZMB.na,patient%in%i)
  #Proliferating.na.sample <- subset(Proliferating.na,patient%in%i)
  CD69.na.sample <- subset(CD69.na,patient%in%i)

  data.clonal.percent$share.kind.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','naive')] <- length(intersect(unique(naive.na.sample$CTaa),unique(CD69.na.sample$CTaa))) / length(unique(CD69.na.sample$CTaa))
  data.clonal.percent$share.kind.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMK')] <- length(intersect(unique(GZMK.na.sample$CTaa),unique(CD69.na.sample$CTaa))) / length(unique(CD69.na.sample$CTaa))
  data.clonal.percent$share.kind.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','GZMB')] <- length(intersect(unique(GZMB.na.sample$CTaa),unique(CD69.na.sample$CTaa))) / length(unique(CD69.na.sample$CTaa))
  #data.clonal.percent$share.kind.percent.69[data.clonal.percent$sample.cluster%in%paste0(i,'_','Proliferating')] <- length(intersect(unique(Proliferating.na.sample$CTaa),unique(CD69.na.sample$CTaa))) / length(unique(CD69.na.sample$CTaa))
}

data.clonal.percent$average.share.kind.percent.103 <- 0
data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'naive'] <- sum(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'naive']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMK'] <- sum(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'GZMK']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMB'] <- sum(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'GZMB']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])

data.clonal.percent$average.share.kind.percent.69 <- 0
data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'naive'] <- sum(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'naive']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMK'] <- sum(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'GZMK']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMB'] <- sum(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'GZMB']) / length(data.clonal.percent$clonalcell.percent[data.clonal.percent$cluster%in%'GZMB'])

data.clonal.percent$sd.share.kind.percent.103 <- 0
data.clonal.percent$sd.share.kind.percent.103[data.clonal.percent$cluster%in%'naive'] <- sd(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$sd.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMK'] <- sd(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$sd.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMB'] <- sd(data.clonal.percent$share.kind.percent.103[data.clonal.percent$cluster%in%'GZMB'])

data.clonal.percent$sd.share.kind.percent.69 <- 0
data.clonal.percent$sd.share.kind.percent.69[data.clonal.percent$cluster%in%'naive'] <- sd(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'naive'])
data.clonal.percent$sd.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMK'] <- sd(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'GZMK'])
data.clonal.percent$sd.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMB'] <- sd(data.clonal.percent$share.kind.percent.69[data.clonal.percent$cluster%in%'GZMB'])

mydata <- as.data.frame(matrix(nrow = 6,ncol = 3))
colnames(mydata) <- c('share_cluster','cluster','share_kind_percent')
mydata$share_cluster <- rep(c('share_kind_gzmk','share_kind_gzmb','other'),time = 2)
mydata$cluster <- rep(c('CD103','CD69'),each=3)
mydata$share_kind_percent[mydata$cluster%in%'CD103'&mydata$share_cluster%in%'share_kind_gzmk'] <- data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMK'][1]
mydata$share_kind_percent[mydata$cluster%in%'CD103'&mydata$share_cluster%in%'share_kind_gzmb'] <- data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'GZMB'][1]
#mydata$share_kind_percent[mydata$cluster%in%'CD103'&mydata$share_cluster%in%'share_kind_pro'] <- data.clonal.percent$average.share.kind.percent.103[data.clonal.percent$cluster%in%'Proliferating'][1]
mydata$share_kind_percent[mydata$cluster%in%'CD103'&mydata$share_cluster%in%'other'] <- 1-sum(mydata[1:2,3])

mydata$share_kind_percent[mydata$cluster%in%'CD69'&mydata$share_cluster%in%'share_kind_gzmk'] <- data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMK'][1]
mydata$share_kind_percent[mydata$cluster%in%'CD69'&mydata$share_cluster%in%'share_kind_gzmb'] <- data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'GZMB'][1]
#mydata$share_kind_percent[mydata$cluster%in%'CD69'&mydata$share_cluster%in%'share_kind_pro'] <- data.clonal.percent$average.share.kind.percent.69[data.clonal.percent$cluster%in%'Proliferating'][1]
mydata$share_kind_percent[mydata$cluster%in%'CD69'&mydata$share_cluster%in%'other'] <- 1-sum(mydata[4:5,3])

pdf(file = 'figure 5D.pdf',width = 5,height = 5)
ggplot(mydata, aes(x = cluster,y = share_kind_percent,fill = share_cluster)) +
  geom_bar(width = 0.3,stat="identity",position="stack") +
  scale_fill_manual(values=c("grey","mediumpurple1","chartreuse3","cornflowerblue"))+
  #scale_fill_brewer(palette = "Blues") +
  scale_y_continuous(labels = scales::percent)+  theme_bw() +
  labs(x = '',y = ' variety percent of each cluster share clone in TRM') +
  guides(fill = guide_legend(ncol = 1, bycol = TRUE, override.aes = list(size = 5))) +#定义图例的布局，1列，排序，图例中色块的大小增大5倍
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 14,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 10),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 12,angle = 45,vjust = 0.5), #x轴标签偏转45°，并下降0.5
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10))
dev.off()

###########
my.so <- ProjectDim(my_ss_CD8, reduction = "pca")
expression_matrix <- my.so@assays$RNA@counts
cell_metadata <- my.so@meta.data
if (all.equal(colnames(expression_matrix), rownames(cell_metadata))) {
  print(sprintf("Cell identifiers match"))
} else {
  print(sprintf("Cell identifier mismatch - %i cells in expression matrix, %i cells in metadata",
                ncol(expression_matrix), nrow(cell_metadata)))
  print("If the counts are equal, sort differences will throw this error")
}
gene_annotation <- data.frame(gene_short_name = rownames(my.so@assays$RNA), row.names = rownames(my.so@assays$RNA))
if (all.equal(rownames(expression_matrix), rownames(gene_annotation))) {
  print(sprintf("Gene identifiers all match"))
} else {
  print(sprintf("Gene identifier mismatch - %i genes in expression matrix, %i gene in gene annotation",
                nrow(expression_matrix), nrow(gene_annotation)))
  print("If the counts are equal, sort differences will throw this error")
}
my.cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
reducedDim(my.cds, type = "PCA") <- my.so@reductions$pca@cell.embeddings
my.cds@preprocess_aux$prop_var_expl <- my.so@reductions$pca@stdev
plot_pc_variance_explained(my.cds)
my.cds@int_colData@listData$reducedDims$UMAP <- my.so@reductions$umap@cell.embeddings
my.cds@clusters$UMAP_so$clusters <- my.so@meta.data$gt_tp_cell_type_integrated_.0.9
my.cds <- cluster_cells(my.cds, reduction_method = "UMAP", resolution = 1e-4)
my.cds <- learn_graph(my.cds)
cds <- order_cells(my.cds)

#figure 5E
GZMK <- subset(my_ss_CD8@meta.data,rename%in%'GZMK_CD8TCM')
GZMK.na = subset(GZMK,CTaa!='NA')

#share clone
SG_clone <- subset(my_ss_CD8@meta.data,rename%in%'CD69_CD8TRM')
SG_clone.na <- subset(SG_clone,CTaa != 'NA')
GZMK <- subset(my_ss_CD8@meta.data,rename%in%'GZMK_CD8TCM')
GZMK.na = subset(GZMK,CTaa!='NA')
my.cds@colData$rename1 <- 'other'
my.cds@colData$rename1[rownames(my.cds@colData)%in%rownames(GZMK.na)[GZMK.na$CTaa%in%SG_clone.na$CTaa]] <- 'clonal_expanded'

pdf(file = 'figure 5E.pdf',width = 6.5,height = 5)
plot_cells(my.cds,color_cells_by = 'rename1',cell_size = 1) + scale_color_brewer(palette="Set1")
dev.off()

#figure 5F
pdf(file = 'figure 5F.pdf',width = 6.5,height = 5)
plot_cells(cds, color_cells_by = "pseudotime", group_label_size = 0)
dev.off()

#figure 5G
GZMK_share_cell <- rownames(my.cds@colData)[my.cds@colData$rename1%in%'clonal_expanded']
my_ss_CD8$rename1 <- 'other'
my_ss_CD8$rename1[rownames(my_ss_CD8@meta.data)%in%GZMK_share_cell] <- 'clonal_expanded'
Idents(my_ss_CD8) <- 'rename1'
DefaultAssay(my_ss_CD8) <- 'RNA'
b=FindMarkers(subset(SS.CD8.integrated,group%in%'SS_PBMC'),ident.1 = "clonal_expanded",ident.2 = "other",logfc.threshold = 0.25)
b=b[order(b$avg_log2FC,decreasing = F),]
myvolcano(b,gene.plot = 100,width = 6.5,height = 5,text.repel = T,pt.size = 2,max.overlaps=100,filename = 'b.label.T')
write.csv(b,file = 'figure 5G.csv')

#figure S1
my_B_cells <- readRDS('my_B_cells.rds')

#figure S1A
pdf(file = 'figure S1A.pdf',width = 6.5,height = 5)
DimPlot(my_B_cells,col = color_used)
dev.off()

#figure S1B
DefaultAssay(my_B_cells) <- 'RNA'
features <- c('IGHD', 'IGHM', 'TCL1A', 'CXCR4',
              'CD82', 'HLA-A',  'HLA-C', 'GPR183', 'MARCKS','XBP1',
              'MZB1', 'SDC1', 'TIMP1','JCHAIN', 'IL6R', 'CD38', 'TNFRSF17', 'AQP3')
pdf(file = 'figS1B.pdf',width =20, height = 5.18)
f3=DotPlot(object = my_B_cells, features = features, cols = c("blue","firebrick2"),
           dot.scale = 5) + RotatedAxis()+
  theme(axis.text= element_text(colour = "black",face='bold'),
        axis.title = element_text(face = 'bold'),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.5,size = 10))
print(f3)
dev.off()

#figure S1C
source('function.R')
percentplot_group(my_B_cells,filename = 'figure S1B', colors = color_used)

#figure S1D
pdf(file = 'figure S1D.pdf',width = 6.5,height = 5)
FeaturePlot(my_B_cells,features = c('IGHM','IGHA1','IGHG1','IGHE'),ncol = 2)
dev.off()

#the file were given by the immcantation pipeline used filtered_contig.fasta and filtered_contig_annotations.csv
data.hc <- read.table('10X_hc_clone-pass.tsv',sep='\t',header = 1)
data.s0 <- read.table('10X_s0_clone-pass.tsv',sep='\t',header = 1)
data.s1 <- read.table('10X_s1_clone-pass.tsv',sep='\t',header = 1)
data.s2 <- read.table('10X_s2_clone-pass.tsv',sep='\t',header = 1)
unique(data.hc$c_call)

data.hc$cell_id <- paste0('hc_',data.hc$cell_id)
data.s0$cell_id <- paste0('s0_',data.s0$cell_id)
data.s1$cell_id <- paste0('s1_',data.s1$cell_id)
data.s2$cell_id <- paste0('s2_',data.s2$cell_id)
data.hc$sample_id <- 'hc'
data.s0$sample_id <- 's0'
data.s1$sample_id <- 's1'
data.s2$sample_id <- 's2'

intersect(data.hc$clone_id,data.s0$clone_id)
intersect(data.hc$clone_id,data.s1$clone_id)#[1] "48_611"  "756_626"
intersect(data.hc$clone_id,data.s2$clone_id)
intersect(data.s0$clone_id,data.s1$clone_id)
intersect(data.s1$clone_id,data.s2$clone_id)#[1] "358_398"
data <- rbind(data.hc,rbind(data.s0,rbind(data.s1,data.s2)))

hc.cell <- subset(my_all_cells@meta.data,group%in%'hc')
rownames(hc.cell) <- paste0('hc_',rownames(hc.cell))
data$sample_id[data$cell_id%in%rownames(hc.cell)[hc.cell$samples%in%'HC5']] <- 'hc1'
data$sample_id[data$cell_id%in%rownames(hc.cell)[hc.cell$samples%in%'HC6']] <- 'hc2'
data$sample_id[data$cell_id%in%rownames(hc.cell)[hc.cell$samples%in%'HC7']] <- 'hc3'
data$sample_id[data$cell_id%in%rownames(hc.cell)[hc.cell$samples%in%'HC8']] <- 'hc4'
data <- data[!data$sample_id%in%'hc',]
data$group[data$sample_id%in%c('hc1','hc2','hc3','hc4')]<- 'hc'
data$group[data$sample_id%in%c('s0','s1','s2')]<- 'ss'

data$clone_id[data$sample_id%in%'s1' & data$clone_id == '48_611'] <- '48_612'
data$clone_id[data$sample_id%in%'s1' & data$clone_id == '756_626'] <- '756_627'
data$clone_id[data$sample_id%in%'s2' & data$clone_id == '358_398'] <- '358_399'

#figure S1F
sample_curve <- alphaDiversity(data, group="group", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0, nboot=100,min_n = 40)

main_title <- 'Sample diversity'
y_value <- 'd'
y_min <- "d_lower"
y_max <- "d_upper"
y_label <- expression(""^q * D)
p1 <- ggplot(sample_curve@diversity, aes_string(x = "q",
                                        y = y_value, group = sample_curve@group_by)) + ggtitle(main_title) +
  baseTheme() + xlab("q") + ylab(y_label) + geom_ribbon(aes_string(ymin = y_min,
                                                                   ymax = y_max, fill = sample_curve@group_by), alpha = 0.4) +
  geom_line(aes_string(color = sample_curve@group_by),size = 1)
p1 <- p1 + scale_color_manual(name = legend_title,
                                labels = group_labels, values = colors) + scale_fill_manual(name = legend_title,
                                                                                            labels = group_labels, values = colors)
pdf(file = 'figure S1F.pdf',width = 6.5,height = 5)
p1
dev.off()

#figure S1G
curve <- estimateAbundance(data, group="group", ci=0, nboot=100, clone="clone_id")
group_labels <- setNames(curve@groups, curve@groups)
legend_title <- 'Rank Abundance'
colors <- c('green','blue','red')
p1 <- ggplot(curve@abundance, aes_string(x = "rank",
                                        y = "p", group = curve@group_by)) + ggtitle(legend_title) +
  baseTheme() + xlab("Rank") + ylab("Abundance") +
  scale_x_log10(limits = NULL, breaks = scales::trans_breaks("log10",function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(labels = scales::percent) +
  geom_ribbon(aes_string(ymin = "lower", ymax = "upper",fill = curve@group_by), alpha = 0.4) +
  geom_line(aes_string(color = curve@group_by),size=1)
p1 <- p1 + scale_color_manual(name = legend_title,labels = group_labels, values = colors) +
  scale_fill_manual(name = legend_title,labels = group_labels, values = colors)
pdf(file = 'figure S1G.pdf',width = 6.5,height = 5)
p1
dev.off()
