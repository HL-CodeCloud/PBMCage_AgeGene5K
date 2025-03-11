library(Seurat)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(dittoSeq)
library(ggsignif)
library(tidyr)
library(tidyverse)
library(dplyr)
library(clusterProfiler)



### Further identify subsets of CD8 T and NK cells
#####################################
###

plotlist=list()

### Identify the subsets of CD8 Tcm
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD8_Tcm=subset(object, reassigned %in% c("CD8 Tcm"))

library(dplyr)
CD8_Tcm=CD8_Tcm %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD8_Tcm_try=FindClusters2(CD8_Tcm, cluster.range=c(3,4), by=0.1, res=0.3, verbose=T)
DimPlot(CD8_Tcm_try, label=T, label.size=5, repel=T) + NoLegend()

# plot T cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
library(patchwork)
CD8_Tcm_try=AddModuleScore_UCell(CD8_Tcm_try, features=T_cell_scoring)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD8 Tcm", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
plotlist[[1]]=
  p_title /
  SCpubr::do_ViolinPlot(CD8_Tcm_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD8_Tcm_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: early differentiation
VlnPlot(CD8_Tcm_try, features=c("GZMK","LYAR","NELL2"), pt.size=0)
# cluster1: cytotoxic Tcm (late differentiation maybe)
VlnPlot(CD8_Tcm_try, features=c("NKG7","GNLY","GZMH","GZMB"), pt.size=0)
# cluster2: adhesion
VlnPlot(CD8_Tcm_try, features=c("VIM","PLP2","LGALS1"), pt.size=0)
# cluster3: transcriptional regulation
VlnPlot(CD8_Tcm_try, features=c("MALAT1","NEAT1","DDX17"), pt.size=0)

library(scRNAtoolVis)
features_=list(early_diff=c("GZMK","LYAR","NELL2"), 
               cytotoxic=c("NKG7","GNLY","GZMB"), 
               adhesion=c("VIM","PLP2","LGALS1"), 
               transcrip_reg=c("MALAT1","NEAT1","DDX17"))
plotlist[[2]]=
  jjDotPlot(CD8_Tcm_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD8 Tcm")

# add the annotation
obj_annotation=data.frame(cell=colnames(object), reassigned=object$reassigned)
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tcm_try, idents=0), "CD8 Tcm.LYAR+", reassigned)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tcm_try, idents=1), "CD8 Tcm.GZMB+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tcm_try, idents=2), "CD8 Tcm.LGALS1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tcm_try, idents=3), "CD8 Tcm.MALAT1+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")


### Identify the subsets of GZMB+ CD8 Tem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD8_Tem.GZMB=subset(object, reassigned %in% c("GZMB+ CD8 Tem"))

library(dplyr)
CD8_Tem.GZMB=CD8_Tem.GZMB %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD8_Tem.GZMB_try=FindClusters2(CD8_Tem.GZMB, cluster.range=c(3,4), by=0.1, res=0.1, verbose=T)
DimPlot(CD8_Tem.GZMB_try, label=T, label.size=5, repel=T) + NoLegend()
table(CD8_Tem.GZMB_try$seurat_clusters)
CD8_Tem.GZMB_try$seurat_clusters=ifelse(CD8_Tem.GZMB_try$seurat_clusters %in% c(0,2,3), 0, 1)
Idents(CD8_Tem.GZMB_try)="seurat_clusters"
CD8_Tem.GZMB_try$seurat_clusters=as.factor(CD8_Tem.GZMB_try$seurat_clusters)

# plot T cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD8_Tem.GZMB_try=AddModuleScore_UCell(CD8_Tem.GZMB_try, features=T_cell_scoring)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="GZMB+ CD8 Tem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
plotlist[[3]]=
  p_title /
  SCpubr::do_ViolinPlot(CD8_Tem.GZMB_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD8_Tem.GZMB_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: mature cytotoxic
VlnPlot(CD8_Tem.GZMB_try, features=c("FCGR3A","PRF1","GZMB"), pt.size=0)
# cluster1: less mature, but more proliferative
VlnPlot(pbmc.seu_merged, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)

library(scRNAtoolVis)
features_=list(mature_cytotoxic=c("FCGR3A","PRF1","GZMB"),
               proliferation=c("C1orf56","HNRNPH1","CDC42SE1"))
plotlist[[4]]=
  jjDotPlot(CD8_Tem.GZMB_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("GZMB+ CD8 Tem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMB_try, idents=0), "CD8 Tem.GZMB+ FCGR3A+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMB_try, idents=1), "CD8 Tem.GZMB+ C1orf56+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets.pdf")
for (i in 1:length(plotlist)) plot(plotlist[[i]])
dev.off()


### Identify the subsets of GZMK+ CD8 Tem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD8_Tem.GZMK=subset(object, reassigned %in% c("GZMK+ CD8 Tem"))

library(dplyr)
CD8_Tem.GZMK=CD8_Tem.GZMK %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD8_Tem.GZMK_try=FindClusters2(CD8_Tem.GZMK, cluster.range=c(4,5), by=0.1, res=0.2, verbose=T)
DimPlot(CD8_Tem.GZMK_try, label=T, label.size=5, repel=T) + NoLegend()
table(CD8_Tem.GZMK_try$seurat_clusters)

# plot T cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD8_Tem.GZMK_try=AddModuleScore_UCell(CD8_Tem.GZMK_try, features=T_cell_scoring)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="GZMK+ CD8 Tem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
plotlist[[5]]=
  p_title /
  SCpubr::do_ViolinPlot(CD8_Tem.GZMK_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD8_Tem.GZMK_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: terminal effector T (Tte)-like cells (also a subset of conventional CD8 Tmem)
VlnPlot(CD8_Tem.GZMK_try, features=c("GZMB","GNLY","NKG7","ZEB2","GZMA"), pt.size=0)
# cluster1: CCR7- GZMKhi, the conventional Tem cells (cytolytic and generally lacked the memory markers LEF1, CD27, CD28, and CD127)
VlnPlot(CD8_subset_try, features=c("CCR7","GZMK","LEF1","CD27","CD28"), pt.size=0)
# cluster2: NKG2C+GZMB–XCL1+, a unique CD8 Tmem that decreases during aging
VlnPlot(CD8_Tem.GZMK_try, features=c("ZNF683","KLRC2","GZMB","LEF1","XCL1"), pt.size=0)
# cluster3: CD8+ TEMRA cells (terminally differentiated, highly cytotoxic, exhausted at some extent, expressing IFN-I genes)
VlnPlot(CD8_Tem.GZMK_try, features=c("GZMH","IFI44L","MX1","EPSTI1","OAS3"), pt.size=0)
# cluster4: lineage infidelity to B cells
VlnPlot(CD8_Tem.GZMK_try, features=c("CD3D","CD8A","MS4A1","IGHD","IGKC"), pt.size=0)

library(scRNAtoolVis)
features_=list(Tte_like=c("GZMB","GNLY","NKG7","ZEB2","GZMA"),
               conv_Tem=c("CCR7","GZMK","LEF1","CD27","CD28"),
               NKG2Cpos.Tem=c("ZNF683","KLRC2","GZMB","LEF1","XCL1"),
               CD8_Temra=c("GZMH","IFI44L","MX1","EPSTI1","OAS3"),
               lin_infidelity=c("CD3G","CD8A","MS4A1","IGHD","IGKC"))
plotlist[[6]]=
  jjDotPlot(CD8_Tem.GZMK_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("GZMK+ CD8 Tem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMK_try, idents=0), "CD8 Tem.GZMK+ NKG7+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMK_try, idents=1), "CD8 Tem.GZMKhi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMK_try, idents=2), "CD8 Tem.GZMK+ XCL1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMK_try, idents=3), "CD8 Tem.GZMK+ IFI44L+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8_Tem.GZMK_try, idents=4), "CD8 Tem.lin_infidel.MS4A1+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets.pdf")
for (i in 1:length(plotlist)) plot(plotlist[[i]])
dev.off()


### Identify the subsets of CD56dim NK
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD56dim.NK=subset(object, reassigned %in% c("CD56dim NK"))

library(dplyr)
CD56dim.NK=CD56dim.NK %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD56dim.NK_try=FindClusters2(CD56dim.NK, cluster.range=c(5,8), by=0.1, res=0.3, verbose=T)
DimPlot(CD56dim.NK_try, label=T, label.size=5, repel=T) + NoLegend()
table(CD56dim.NK_try$seurat_clusters)
CD56dim.NK_try$seurat_clusters[CD56dim.NK_try$seurat_clusters==6]=5
CD56dim.NK_try$seurat_clusters[CD56dim.NK_try$seurat_clusters==7]=0
Idents(CD56dim.NK_try)="seurat_clusters"
CD56dim.NK_try$seurat_clusters=as.factor(CD56dim.NK_try$seurat_clusters)
levels(CD56dim.NK_try$seurat_clusters)=c("0","1","2","3","4","5","5","0")

# plot NK cell subsets markers
NK_cell_scoring=list(immature=c("ICAM1+","NCAM1+","B3GAT1-","SELL+","IL2RA+","IL7R+","IL18R1+","KIT+","CCR5+","CCR7+","CXCR1-","CXCR3+","CXCR4+","CX3CR1-","CMKLR1-","CD27+","CD226+","KLRF1+","KIR2DS1-","KIR2DS2-","KIR2DS3-","KIR2DL4-","KIR2DS4-","KIR2DS5-","KIR3DS1-","KLRC2-","NCR3+","NCR2+","NCR1+","KLRC1+","KLRD1+","KIR2DL1-","KIR2DL2-","KIR2DL3-","KIR2DL4-","KIR2DL5-","KIR3DL1-","KIR3DL2-","LILRB1-","PDCD1-"),
                     maturing=c("ITGAL+","NCAM1+","IL2RA-","IL7R-","IL18R1+","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","KIR2DS1-","KIR2DS2-","KIR2DS3-","KIR2DL4-","KIR2DS4-","KIR2DS5-","KIR3DS1-","NCR3+","NCR1+","KLRC1+","KLRD1+"),
                     double_pos=c("ITGAL+","NCAM1+","IL2RA-","IL7R-","IL18R1+","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","NCR3+","NCR1+","KLRC1+","KLRD1+"),
                     mature=c("ITGAL+","NCAM1+","IL2RA-","IL7R-","IL18R1+","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","NCR3+","NCR1+","KLRC1-","KLRD1-"),
                     terminally_diff=c("ITGAL+","NCAM1+","B3GAT1+","SELL+","IL2RA-","IL7R-","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","NCR3+","NCR1+","KLRC1-","KLRD1-"),
                     memory_like=c("ITGAL+","NCAM1+","B3GAT1+","SELL+","IL2RA-","IL7R-","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","NCR3+","NCR1+","KLRC1-","KLRD1-"),
                     CD56neg_CD16dim=c("ITGAL+","NCAM1-","SELL+","IL2RA-","IL7R-","KIT-","CCR5-","CCR7-","CXCR1+","CXCR3+","CXCR4+","CX3CR1+","CMKLR1+","CD27-","CD226+","KLRF1+","FCGR3A+","NCR1+"),
                     CD56dim_CD16neglo=c("NCAM1+","B3GAT1-","SELL+","IL2RA+","IL7R+","IL18R1+","KIT-","CCR5+","CCR7-","CXCR1-","CXCR3+","CXCR4+","CX3CR1-","CMKLR1-","CD226+","KLRF1+","FCGR3A-","KIR2DS1-","KIR2DS2-","KIR2DS3-","KIR2DL4-","KIR2DS4-","KIR2DS5-","KIR3DS1-","NCR3+","NCR1+","KLRC1+","KLRD1+","KIR2DL1-","KIR2DL2-","KIR2DL3-","KIR2DL4-","KIR2DL5-","KIR3DL1-","KIR3DL2-"),
                     activation=c("NCR2","FAS","FASLG","CD40LG","TNFSF10","ITGAM","CD2","CD58"),
                     degranulation=c("LAMP1","LAMP2"),
                     circulating=c("ITGAM","ITGB2"))

NK_cell_scoring=list(adhesion=c("ITGAL","CD2","CD58","ITGAM","ICAM1","NCAM1","B3GAT1","SELL"),
                     cytokine=c("IL2RA","IL2RB","IL2RG","IL7R","IL18R1","KIT","CCR5","CCR7","CXCR1","CXCR3","CXCR4","CX3CR1","CMKLR1"),
                     death_R=c("FAS","FASLG","CD40LG","TNFSF10"),
                     degranulation=c("LAMP1","LAMP2"),
                     activating_coR=c("CD244","CD59","CD226","KLRF1","SLAMF6"),
                     activating_R=c("SLAMF7","FCGR3A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DL4","KIR2DS4","KIR2DS5","KIR3DS1","KLRC2","KLRK1","NCR3","NCR2","NCR1","CD96"),
                     inhibition=c("KLRC1","KLRD1","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5","KIR3DL1","KIR3DL2","LILRB1","PDCD1","SIGLEC7","CD300A","TIGIT"))
library(UCell)
CD56dim.NK_try=AddModuleScore_UCell(CD56dim.NK_try, features=NK_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD56dim NK", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD56dim.NK_try,
                        features=paste0(names(NK_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(NK_cell_scoring)),
                        ylab=names(NK_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD56dim.NK_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: cytotoxic 
VlnPlot(CD56dim.NK_try, features=c("MYOM2","FCER1G","PRF1","SPON2"), pt.size=0)
# cluster1: less activated but more proliferative
VlnPlot(CD56dim.NK_try, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster2: adaptive NK (different from adaptive-like because KLRC2+, which is a marker of human cytomegalovirus/HCMV infection)
VlnPlot(CD56dim.NK_try, features=c("KLRC2","CD3E","CCL5","VIM"), pt.size=0)
# cluster3: transcriptionally active (this is a homeostatic activation; apart from below also express "NFKBIA","CD69","CXCR4","ZFP36")
VlnPlot(CD56dim.NK_try, features=c("FOS","FOSB","JUN","JUNB"), pt.size=0)
# cluster4: adaptive-like NK (because KLRC2-; this is possibly a subset derived from infections)
VlnPlot(CD56dim.NK_try, features=c("CD3E","VIM","IL32","CD52","KLRC2","FCER1G"), pt.size=0)
# cluster5: IFN-I+ NK (reported to have been found in early severe COVID-19 NK cells, so it's related to infections)
VlnPlot(CD56dim.NK_try, features=c("MX1","ISG15","ISG20"), pt.size=0)


library(scRNAtoolVis)
features_=list(cytotoxic=c("MYOM2","FCER1G","PRF1","SPON2"),
               proliferative=c("C1orf56","HNRNPH1","CDC42SE1"),
               adaptive=c("KLRC2","CD3E","CCL5","VIM"),
               transcript_active=c("FOS","FOSB","JUN","JUNB"),
               adaptive_like=c("CD3E","VIM","IL32","CD52","KLRC2","FCER1G"),
               ISG_pos=c("MX1","ISG15","ISG20"))
p2=
  jjDotPlot(CD56dim.NK_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD56dim NK")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=0), "NK.CD56dim FCER1G+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=1), "NK.CD56dim C1orf56+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=2), "NK.CD56dim FCER1Glo KLRC2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=3), "NK.CD56dim FOS+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=4), "NK.CD56dim FCER1Glo KLRC2-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56dim.NK_try, idents=5), "NK.CD56dim MX1+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_2.pdf")
p1
p2
dev.off()


### Identify the subsets of CD56hi NK
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD56hi.NK=subset(object, reassigned %in% c("CD56hi NK"))

library(dplyr)
CD56hi.NK=CD56hi.NK %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD56hi.NK_try=FindClusters2(CD56hi.NK, cluster.range=c(2,4), by=0.1, res=0.2, verbose=T)
DimPlot(CD56hi.NK_try, label=T, label.size=5, repel=T) + NoLegend()

# plot NK cell subsets markers
NK_cell_scoring=list(adhesion=c("ITGAL","CD2","CD58","ITGAM","ICAM1","NCAM1","B3GAT1","SELL"),
                     cytokine=c("IL2RA","IL2RB","IL2RG","IL7R","IL18R1","KIT","CCR5","CCR7","CXCR1","CXCR3","CXCR4","CX3CR1","CMKLR1"),
                     death_R=c("FAS","FASLG","CD40LG","TNFSF10"),
                     degranulation=c("LAMP1","LAMP2"),
                     activating_coR=c("CD244","CD59","CD226","KLRF1","SLAMF6"),
                     activating_R=c("SLAMF7","FCGR3A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DL4","KIR2DS4","KIR2DS5","KIR3DS1","KLRC2","KLRK1","NCR3","NCR2","NCR1","CD96"),
                     inhibition=c("KLRC1","KLRD1","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5","KIR3DL1","KIR3DL2","LILRB1","PDCD1","SIGLEC7","CD300A","TIGIT"))
library(UCell)
CD56hi.NK_try=AddModuleScore_UCell(CD56hi.NK_try, features=NK_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD56hi NK", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD56hi.NK_try,
                        features=paste0(names(NK_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(NK_cell_scoring)),
                        ylab=names(NK_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD56hi.NK_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: CD56bright NK with responses to cytokines and regulation on cytokine production
VlnPlot(CD56hi.NK_try, features=c("COTL1","XCL1","LTB","GZMK"), pt.size=0)
# cluster1: transitional NK between CD56dim and CD56hi
VlnPlot(CD56hi.NK_try, features=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"), pt.size=0)
# cluster2: CD56bright NK with chemokine-secreting property (may play a critical role in recruitment of T cells and other immune cells)
VlnPlot(CD56hi.NK_try, features=c("CCL4L2","CCL3","CCL4"), pt.size=0)
# cluster3: CD56bright NK which is S100B+ cytotoxic (in the CD56hiCD16-CD49a-KIR-NK cluster, one of subsets of conventional CD56hi NK)
VlnPlot(NK_subset_try, features=c("S100B","MYOM2","BEX3","SOX4","ITGA1","KIR2DL4"), pt.size=0)


library(scRNAtoolVis)
features_=list(cytokine_resp=c("COTL1","XCL1","LTB","GZMK"),
               transitional=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"),
               chemokine_prod=c("CCL4L2","CCL3","CCL4"),
               cytotoxic=c("S100B","MYOM2","BEX3","SOX4","ITGA1","KIR2DL4"))
p2=
  jjDotPlot(CD56hi.NK_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD56hi NK")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56hi.NK_try, idents=0), "NK.CD56bright COTL1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56hi.NK_try, idents=1), "NK.CD56bright FCGR3A+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56hi.NK_try, idents=2), "NK.CD56bright CCL4L2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD56hi.NK_try, idents=3), "NK.CD56bright S100B+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_3.pdf")
p1
p2
dev.off()


### Identify the subsets of CD8 Tnaive
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD8.Tnaive=subset(object, reassigned %in% c("CD8 Tnaive"))

library(dplyr)
CD8.Tnaive=CD8.Tnaive %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD8.Tnaive_try=FindClusters2(CD8.Tnaive, cluster.range=c(3,5), by=0.1, res=0.2, verbose=T)
DimPlot(CD8.Tnaive_try, label=T, label.size=5, repel=T) + NoLegend()

VlnPlot(CD8.Tnaive_try, features=c("LEF1","CD27","SELL","IL6R","IL7R","S100B","CD38"), pt.size=0)
# plot CD8 Tnaive cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD8.Tnaive_try=AddModuleScore_UCell(CD8.Tnaive_try, features=T_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD8 Tnaive", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD8.Tnaive_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD8.Tnaive_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: preparing for TCR-mediated events (molecules important for TCR-mediated survival, activation, signaling)
VlnPlot(CD8.Tnaive_try, features=c("TNFAIP3","DUSP1","GADD45B","CITED2"), pt.size=0)
# cluster1: homeostatic proliferation (spontaneous proliferation and in partial activation)
VlnPlot(CD8.Tnaive_try, features=c("LEF1","SELL","C1orf56","HNRNPH1"), pt.size=0)
# cluster2: cell regulation (events eg. proliferation/cytokine production/migration/apoptosis...)
VlnPlot(CD8_subset_try, features=c("ITGB1","LGALS1","S100A11","ANXA1","IL10RA"), pt.size=0)
# cluster3: transcriptional repressors, essential for maintaining cellular homeostasis and regulating
VlnPlot(CD8.Tnaive_try, features=c("TSHZ2","FHIT","RNASET2"), pt.size=0)


library(scRNAtoolVis)
features_=list(init_activation=c("TNFAIP3","DUSP1","GADD45B","CITED2"),
               homeostatic_prolif=c("LEF1","SELL","C1orf56","HNRNPH1"),
               cell_reg=c("ITGB1","LGALS1","S100A11","ANXA1","IL10RA"),
               transcrip_repress=c("TSHZ2","FHIT","RNASET2"))
p2=
  jjDotPlot(CD8.Tnaive_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD8 Tnaive")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8.Tnaive_try, idents=0), "CD8 Tnaive.TNFAIP3+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8.Tnaive_try, idents=1), "CD8 Tnaive.C1orf56+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8.Tnaive_try, idents=2), "CD8 Tnaive.LGALS1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD8.Tnaive_try, idents=3), "CD8 Tnaive.TSHZ2+", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_3.pdf")
p1
p2
dev.off()


### Identify the subsets of Proliferating T
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Proliferating.T=subset(object, reassigned %in% c("Proliferating T"))

library(dplyr)
Proliferating.T=Proliferating.T %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Proliferating.T_try=FindClusters2(Proliferating.T, cluster.range=c(2,3), by=0.1, res=0.2, verbose=T)
DimPlot(Proliferating.T_try, label=T, label.size=5, repel=T) + NoLegend()

VlnPlot(Proliferating.T_try, features=c("CD3E","CD8A","CD4","LEF1","CD27","SELL","IL7R","CD38"), pt.size=0)
VlnPlot(Proliferating.T_try, features=c("IFNG","TNF","TBX21","PRDM1","ID2","IRF4","IL4","IL5","IL13","GATA3","IL9","IL10","IL17A","IL21","RORC","RORA","IL6","TGFB1","FOXP3"), pt.size=0)
# plot Proliferating T cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
Proliferating.T_try=AddModuleScore_UCell(Proliferating.T_try, features=T_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Proliferating T", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Proliferating.T_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Proliferating.T_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: highly proliferating CD4 T cells (with high infiltrating properties)
VlnPlot(Proliferating.T_try, features=c("MKI67","RRM2","CDCA2","UBE2C","KIF23","MXD3"), pt.size=0)
# cluster1: less proliferating CD4 T cells (proliferation under regulation/suppression)
VlnPlot(Proliferating.T_try, features=c("TSHZ2","PIK3IP1","ZFC3H1","PAG1"), pt.size=0)
# cluster2: highly proliferating CD8+ effector CTL, most are Tc1
VlnPlot(Proliferating.T_try, features=c("GZMH","NKG7","ID2","IFNG"), pt.size=0) 


library(scRNAtoolVis)
features_=list(prolif_CD4T=c("MKI67","RRM2","CDCA2","UBE2C","KIF23","MXD3"),
               prolif_underCtrl=c("TSHZ2","PIK3IP1","ZFC3H1","PAG1"),
               CD8_CTL=c("GZMH","NKG7","ID2","IFNG"))
p2=
  jjDotPlot(Proliferating.T_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Proliferating T")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.T_try, idents=0), "CD4 T.MKI67hi CDCA2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.T_try, idents=1), "CD4 T.MKI67int TSHZ2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.T_try, idents=2), "CD8 CTL", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_3.pdf")
p1
p2
dev.off()


### Identify the subsets of Proliferating NK
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Proliferating.NK=subset(object, reassigned %in% c("Proliferating NK"))

library(dplyr)
Proliferating.NK=Proliferating.NK %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Proliferating.NK_try=FindClusters2(Proliferating.NK, cluster.range=c(3,5), by=0.1, res=0.2, verbose=T)
DimPlot(Proliferating.NK_try, label=T, label.size=5, repel=T) + NoLegend()

VlnPlot(Proliferating.NK_try, features=c("CD3E","CD8A","CD4","LEF1","CD27","SELL","IL7R","CD38"), pt.size=0)
VlnPlot(Proliferating.NK_try, features=c("IFNG","TNF","TBX21","PRDM1","ID2","IRF4","IL4","IL5","IL13","GATA3","IL9","IL10","IL17A","IL21","RORC","RORA","IL6","TGFB1","FOXP3"), pt.size=0)
# plot Proliferating NK cell subsets markers
NK_cell_scoring=list(adhesion=c("ITGAL","CD2","CD58","ITGAM","ICAM1","NCAM1","B3GAT1","SELL"),
                     cytokine=c("IL2RA","IL2RB","IL2RG","IL7R","IL18R1","KIT","CCR5","CCR7","CXCR1","CXCR3","CXCR4","CX3CR1","CMKLR1"),
                     death_R=c("FAS","FASLG","CD40LG","TNFSF10"),
                     degranulation=c("LAMP1","LAMP2"),
                     activating_coR=c("CD244","CD59","CD226","KLRF1","SLAMF6"),
                     activating_R=c("SLAMF7","FCGR3A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DL4","KIR2DS4","KIR2DS5","KIR3DS1","KLRC2","KLRK1","NCR3","NCR2","NCR1","CD96"),
                     inhibition=c("KLRC1","KLRD1","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5","KIR3DL1","KIR3DL2","LILRB1","PDCD1","SIGLEC7","CD300A","TIGIT"))
library(UCell)
Proliferating.NK_try=AddModuleScore_UCell(Proliferating.NK_try, features=NK_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Proliferating NK", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Proliferating.NK_try,
                        features=paste0(names(NK_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(NK_cell_scoring)),
                        ylab=names(NK_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Proliferating.NK_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: less proliferating but highly active CD56dim CD16+ NK cells
VlnPlot(Proliferating.NK_try, features=c("FCGR3A","NCAM1","IL2RB","SPON2","FCER1G","GZMB"), pt.size=0)
# cluster1: proliferating NKT cells
VlnPlot(Proliferating.NK_try, features=c("MKI67","RRM2","CD3G","CD8B","KLRD1","KLRB1"), pt.size=0)
# cluster2: highly proliferating CD56dim CD16+ NK cells
VlnPlot(Proliferating.NK_try, features=c("CDCA3","AURKB","UBE2C","TOP2A"), pt.size=0)
# cluster3: CD56dim CD16+ NK cells under regulation of cell cycling
VlnPlot(Proliferating.NK_try, features=c("HMGB2","BIRC5","SMC4","CENPF"), pt.size=0)
# cluster4: proliferating NK cells with lineage infidelity to B cells
VlnPlot(Proliferating.NK_try, features=c("GNLY","GZMB","CD22","MS4A1","IGHD"), pt.size=0)
VlnPlot(Proliferating.NK_try, features=c("NKG7","KLRD1","TYROBP","GNLY","FCER1G","PRF1","CD247","KLRF1","CST7","GZMB"), pt.size=0)
VlnPlot(Proliferating.NK_try, features=c("CD79A","RALGPS2","CD79B","MS4A1","BANK1","CD74","TNFRSF13C","HLA-DQA1","IGHM","MEF2C"), pt.size=0)
highlighted=WhichCells(Proliferating.NK_try, idents=4)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()


library(scRNAtoolVis)
features_=list(active_NK=c("FCGR3A","NCAM1","IL2RB","SPON2","FCER1G","GZMB"),
               prolif_NKT=c("MKI67","RRM2","CD3G","CD8B","KLRD1","KLRB1"),
               prolif_NK=c("CDCA3","AURKB","UBE2C","TOP2A"),
               cycling_reg_NK=c("HMGB2","BIRC5","SMC4","CENPF"),
               lin_infidelity=c("GNLY","GZMB","CD22","MS4A1","IGHD"))
p2=
  jjDotPlot(Proliferating.NK_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Proliferating NK")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.NK_try, idents=0), "NK.MKI67- RRM2lo GZMBhi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.NK_try, idents=1), "NKT.MKI67+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.NK_try, idents=2), "NK.MKI67hi RRM2hi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.NK_try, idents=3), "NK.MKI67hi RRM2-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Proliferating.NK_try, idents=4), "NK.lin_infidel.MS4A1+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_3.pdf")
p1
p2
dev.off()


### Identify the subsets of Treg
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Treg=subset(object, reassigned %in% c("Treg"))

library(dplyr)
Treg=Treg %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Treg_try=FindClusters2(Treg, cluster.range=c(2,3), by=0.1, res=0.2, verbose=T)
DimPlot(Treg_try, label=T, label.size=5, repel=T) + NoLegend()

# plot Proliferating T cell subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
Treg_try2=AddModuleScore_UCell(Treg_try2, features=T_cell_scoring)


p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Treg", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Treg_try2,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Treg_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: subsets of activated/effector Tregs:
# (HLA-DRhi subset + FOXP3hi effector subset which is FOXP3hi PTPRChi DDX17hi)
# Path I is defined as the developmen path from naïve Tregs ending with the FOXP3hi subset
VlnPlot(Treg_try, features=c("FOXP3","IL2RA","HLA-DRA","PTPRC","DDX17"), pt.size=0)
# cluster1: naive Tregs (at naïve status with CCR7hi TCF7hi HLA-DRlow FOXP3low profile)
VlnPlot(Treg_try, features=c("CCR7","TCF7"), pt.size=0)
# cluster2: contaminated IgA+ PC
VlnPlot(Treg_try, features=c("JCHAIN","IGKC","IGHA1","IGHA2"), pt.size=0)
highlighted=WhichCells(Treg_try, idents=2)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

# remove the contaminated PCs
Treg_try2=subset(Treg_try, seurat_clusters!=2)
DimPlot(Treg_try2, label=T, label.size=5, repel=T) + NoLegend()


library(scRNAtoolVis)
features_=list(Treg_activ.eff=c("FOXP3","IL2RA","HLA-DRA","DDX17"),
               Treg_naive=c("CCR7","TCF7"))
p2=
  jjDotPlot(Treg_try2,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)",
            scale=F) +
  ggtitle("Treg")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Treg_try, idents=0), "Treg.FOXP3hi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Treg_try, idents=1), "Treg.TCF7+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Treg_try, idents=2), "IgA+ PC", CD8T_NK_subsets))
write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")

pdf("~/Project_PBMCage/Plots/CD8T_NK_subsets_3.pdf")
p1
p2
dev.off()

#####################################



### Plot the cell proportion of the CD8 and NK subsets across agecuts
#####################################
###

###
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(rownames(obj_annotation)==colnames(object)) # check
object=AddMetaData(object, metadata=obj_annotation$CD8T_NK_subsets, col.name="CD8T_NK_subsets")

ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)

Age_cuts=list(c(19:35), c(36:50), c(51:65), c(66:80), c(81:90), c(91:97))
AgecutAnnot=ifelse(data$Age %in% Age_cuts[[1]], "19~35", 
                   ifelse(data$Age %in% Age_cuts[[2]], "36~50",
                          ifelse(data$Age %in% Age_cuts[[3]], "51~65",
                                 ifelse(data$Age %in% Age_cuts[[4]], "66~80",
                                        ifelse(data$Age %in% Age_cuts[[5]], "81~90",
                                               ifelse(data$Age %in% Age_cuts[[6]], "91~97", NA))))))
table(AgecutAnnot, useNA="ifany") # check
data[["AgeCut"]]=AgecutAnnot

# make a function
library(ggpubr)
draw_proportion=function(df, celltype.drawn) {
  df.used=subset(df, CellType==celltype.drawn)
  ggplot(df.used, aes(x=Age, y=Freq, color=Sex)) +
    geom_point(size=1, shape=1) +
    scale_color_manual(values=col.sex) +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=12),
          legend.title=element_text(size=12),
          legend.text=element_text(size=10)) +
    labs(title=paste0("Percent of ", celltype.drawn))+
    ylab("Cell percent (%)") +
    geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
    stat_cor(show.legend=FALSE, size=3.5)
}

# plot cell freq vs. ages for each celltype
celltypes_ordered=names(table(data$CellType))
celltypes_ordered=celltypes_ordered[grepl("(CD8)|(NK)|(CD4 T\\.)|(Treg\\.)", celltypes_ordered)]
celltypes_ordered # check the types
plist=list()
for (i in 1:length(celltypes_ordered)){
  current_type=celltypes_ordered[i]
  plist[[i]]=draw_proportion(df=data, celltype.drawn=current_type)
}
pdf("~/Project_PBMCage/Plots/All_cells_eachTypeProportion_vs._age_CD8TNKsubsets.pdf", width=40, height=25)
cowplot::plot_grid(plotlist=plist, nrow=5, align="hv")
dev.off()


# make another function for beeswarm plot
library(ggbeeswarm)
library(ggpubr)
draw_subsets_forAgeCut=
  function(object, celltype.drawn) {
    ID_=object$donor_id
    Age_=object$age
    CellType_=object$CD8T_NK_subsets
    cellNum=table(CellType_, ID_)
    cellProp=t(t(cellNum)/rowSums(t(cellNum)))
    data=t(cellProp*100)
    data=as.data.frame(data)
    colnames(data)=c("ID","CellType","Freq")
    ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
    data=merge(data, ID_Age_match, by="ID")
    data$CellType=as.factor(data$CellType)
    
    Age_cuts=list(c(19:35), c(36:50), c(51:65), c(66:80), c(81:90), c(91:97))
    AgecutAnnot=ifelse(data$Age %in% Age_cuts[[1]], "19~35", 
                       ifelse(data$Age %in% Age_cuts[[2]], "36~50",
                              ifelse(data$Age %in% Age_cuts[[3]], "51~65",
                                     ifelse(data$Age %in% Age_cuts[[4]], "66~80",
                                            ifelse(data$Age %in% Age_cuts[[5]], "81~90",
                                                   ifelse(data$Age %in% Age_cuts[[6]], "91~97", NA))))))
    table(AgecutAnnot, useNA="ifany") # check
    data[["AgeCut"]]=AgecutAnnot
    
    df.used=subset(data, CellType==celltype.drawn)
    p=
      ggplot(df.used, aes(x=AgeCut, y=Freq)) +
      geom_violin(color="gray30") +
      theme_classic() +
      labs(x=NULL, y="Cell percent (%)") +
      theme(axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10),
            axis.title.y=element_text(size=11),
            legend.position="none") +
      labs(title=paste0("Percent of ", celltype.drawn)) +
      stat_compare_means(method="anova")
    
    median_1=unlist(lapply(split(df.used$Freq, df.used$AgeCut), median))
    p1_l=data.frame(AgeCut=names(median_1), Freq=median_1)
    p=p + geom_line(data=p1_l, aes(x=AgeCut, y=Freq, group=1), linetype="dashed", color="black")
    
    
    # plot cell counts instead of freq
    ID_=object$donor_id
    Age_=object$age
    CellType_=object$CD8T_NK_subsets
    cellNum=table(CellType_, ID_)
    data2=as.data.frame(t(cellNum))
    colnames(data2)=c("ID","CellType","Counts")
    ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
    data2=merge(data2, ID_Age_match, by="ID")
    data2$CellType=as.factor(data2$CellType)
    
    Age_cuts=list(c(19:35), c(36:50), c(51:65), c(66:80), c(81:90), c(91:97))
    AgecutAnnot=ifelse(data2$Age %in% Age_cuts[[1]], "19~35", 
                       ifelse(data2$Age %in% Age_cuts[[2]], "36~50",
                              ifelse(data2$Age %in% Age_cuts[[3]], "51~65",
                                     ifelse(data2$Age %in% Age_cuts[[4]], "66~80",
                                            ifelse(data2$Age %in% Age_cuts[[5]], "81~90",
                                                   ifelse(data2$Age %in% Age_cuts[[6]], "91~97", NA))))))
    table(AgecutAnnot, useNA="ifany") # check
    data2[["AgeCut"]]=AgecutAnnot
    
    df.used2=subset(data2, CellType==celltype.drawn)
    p2=
      ggplot(df.used2, aes(x=AgeCut, y=Counts)) +
      geom_violin(color="gray30") +
      theme_classic() +
      labs(x=NULL, y="Cell counts") +
      theme(axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10),
            axis.title.y=element_text(size=11),
            legend.position="none") +
      labs(title=paste0("Counts of ", celltype.drawn)) +
      stat_compare_means(method="anova")
    
    median_2=unlist(lapply(split(df.used2$Counts, df.used2$AgeCut), median))
    p2_l=data.frame(AgeCut=names(median_2), Counts=median_2)
    p2=p2 + geom_line(data=p2_l, aes(x=AgeCut, y=Counts, group=1), linetype="dashed", color="black")
    
    return(list(p, p2))
  }

# plot cell freq and counts vs. agecuts for each celltype
celltypes_ordered=names(table(data$CellType))
celltypes_ordered==celltypes_ordered[grepl("(CD8)|(NK)|(CD4 T\\.)|(Treg\\.)", celltypes_ordered)]
celltypes_ordered # check the types
plist=list()
for (i in 1:length(celltypes_ordered)){
  current_type=celltypes_ordered[i]
  ps=draw_subsets_forAgeCut(object=object, celltype.drawn=current_type)
  plist=c(plist, list(ps[[1]]), list(ps[[2]]))
}
pdf("~/Project_PBMCage/Plots/All_cells_eachTypeCountsProp_vs._age_CD8TNKsubsets.pdf", width=40, height=25)
cowplot::plot_grid(plotlist=plist, nrow=6, align="hv")
dev.off()

#####################################



### Take CD8 T and NK for cytotoxicity scoring
#####################################
###

### Take the celltypes
# add celltype annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(rownames(obj_annotation)==colnames(object)) # check
object=AddMetaData(object, metadata=obj_annotation$CD8T_NK_subsets, col.name="CD8T_NK_subsets")
# take the subsets
all_celltypes=names(table(object$CD8T_NK_subsets))[!grepl("lin_infidel", names(table(object$CD8T_NK_subsets)))]

### Make functions for scoring
# a function for getting the df to draw
score_df_for_drawing=function(rdspath, object=NULL, celltypeKEY_for_grepl=NULL, gensetlist=NULL) {
  if (length(rdspath)==1) {
    subset_=readRDS(rdspath)
  } else {
    typenames=names(table(object$CD8T_NK_subsets))[grepl(celltypeKEY_for_grepl, names(table(object$CD8T_NK_subsets)))]
    subset_=subset(object, (CD8T_NK_subsets %in% typenames) & (CD8T_NK_subsets %in% all_celltypes))
    Idents(subset_)="CD8T_NK_subsets"
    library(UCell)
    subset_=AddModuleScore_UCell(subset_, features=gensetlist)
    filename=gsub("[[:punct:]]","",celltypeKEY_for_grepl)
    saveRDS(subset_, paste0("~/Project_PBMCage/Tempt_RDS/Scoring_",filename,".rds"))
  }
  
  score_df=subset_[[]][, colnames(subset_[[]])[grepl("_UCell$", colnames(subset_[[]]))]]
  score_df=cbind(score_df, subset_[[]][, c("age","CD8T_NK_subsets")])
  Age_cuts=list(c(19:35), c(36:50), c(51:65), c(66:80), c(81:90), c(91:97))
  AgecutAnnot=ifelse(score_df$age %in% Age_cuts[[1]], "19~35", 
                     ifelse(score_df$age %in% Age_cuts[[2]], "36~50",
                            ifelse(score_df$age %in% Age_cuts[[3]], "51~65",
                                   ifelse(score_df$age %in% Age_cuts[[4]], "66~80",
                                          ifelse(score_df$age %in% Age_cuts[[5]], "81~90",
                                                 ifelse(score_df$age %in% Age_cuts[[6]], "91~97", NA))))))
  score_df$AgeCut=AgecutAnnot
  score_df$AgeCut=as.factor(score_df$AgeCut)
  
  return(score_df)
}
# a function for plotting on ages
geneset_scoring_forAge=function(dataframe, celltype_as_LegendName, scale) {
  score_df=dataframe
  library(ComplexHeatmap)
  score_df_agemean=aggregate(score_df[,(1:(ncol(score_df)-3))], list(Age=score_df$age, Celltype=score_df$CD8T_NK_subsets), FUN=mean)
  Celltypes_in_df=names(table(score_df$CD8T_NK_subsets))
  df_full_list=list()
  for (i in 1:length(Celltypes_in_df)) {
    df_=score_df_agemean[score_df_agemean$Celltype==Celltypes_in_df[i],]
    age_not_there=as.integer(levels(score_df$age)[sapply(levels(score_df$age), function(x) !(x %in% df_$Age))])
    df_empty=data.frame(matrix(nrow=length(age_not_there), ncol=ncol(df_)))
    colnames(df_empty)=colnames(df_)
    df_empty$Age=age_not_there
    tryCatch({df_empty$Celltype=Celltypes_in_df[i]}, error=function(msg) {print("Pass.")})
    df_full=data.table::rbindlist(list(df_, df_empty))
    df_full=df_full[order(df_full$Age),]
    df_full_list[[i]]=df_full
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  df_full_combined=df_full_combined[, 3:ncol(df_full_combined)]
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[2]*100)/100), abs(floor(quantiles[4]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/2, 0, ranges/2), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-colMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/4, 0, ranges/4), c("#00007E", "white", "#7E0000"))}
  
  p_age=list()
  codes=""
  for (i in 1:length(df_full_list)) {
    df_=as.data.frame(df_full_list[[i]])
    rownames(df_)=df_$Age
    df_to_draw=df_[,3:ncol(df_)]
    colnames(df_to_draw)=gsub("_UCell", "", colnames(df_to_draw))
    celltype_of_df=names(table(df_$Celltype))
    df_to_draw=t(df_to_draw)
    df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)
    
    Age_cuts=list(c(19:35), c(36:50), c(51:65), c(66:80), c(81:90), c(91:97))
    column_split_index_in_df=sapply(as.numeric(colnames(df_to_draw)), function(e) which(sapply(Age_cuts, function(agerange) is.element(e, agerange))))
    agecuts_in_df=c("19~35","36~50","51~65","66~80","81~90","91~97")[column_split_index_in_df]
    ha_=HeatmapAnnotation(
      age=anno_block(labels=agecuts_in_df[!duplicated(agecuts_in_df)],
                     labels_gp=gpar(col="black", fontsize=8),
                     gp=gpar(col="transparent"),
                     height=unit(2,"mm")))
    if (i==1) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                name=paste0("UCell of ", celltype_as_LegendName),
                heatmap_legend_param=list(title_position="leftcenter-rot",
                                          title_gp=gpar(fontsize=8, fontface="bold")),
                column_title=NULL,
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=F,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                column_split=column_split_index_in_df,
                top_annotation=ha_,
                row_names_gp=gpar(fontsize=8),
                width=ncol(df_to_draw)*unit(3,"mm"), 
                height=nrow(df_to_draw)*unit(3,"mm"))}
    if (i>=2) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                show_heatmap_legend=FALSE,
                column_title=NULL,
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=F,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                column_split=column_split_index_in_df,
                row_names_gp=gpar(fontsize=8),
                width=ncol(df_to_draw)*unit(3,"mm"), 
                height=nrow(df_to_draw)*unit(3,"mm"))
    }
    codes=paste0(codes, " %v% p_age[[", i, "]]")
  }
  codes=gsub("^ \\%v\\% ","",codes)
  eval(parse(text=paste0("ht_list=", codes)))
  p_age_total=draw(ht_list, ht_gap=unit(1, "mm"))
  
  return(p_age_total)
}
# a function for plotting on AgeCuts
geneset_scoring_forAgeCuts=function(dataframe, celltype_as_LegendName, scale) {
  score_df=dataframe
  library(ComplexHeatmap)
  score_df_agemean=aggregate(score_df[,(1:(ncol(score_df)-3))], list(Age_cut=score_df$AgeCut, Celltype=score_df$CD8T_NK_subsets), FUN=mean)
  Celltypes_in_df=names(table(score_df$CD8T_NK_subsets))
  df_full_list=list()
  for (i in 1:length(Celltypes_in_df)) {
    df_=score_df_agemean[score_df_agemean$Celltype==Celltypes_in_df[i],]
    AgeCut_not_there=levels(score_df$AgeCut)[sapply(levels(score_df$AgeCut), function(x) !(x %in% df_$Age_cut))]
    df_empty=data.frame(matrix(nrow=length(AgeCut_not_there), ncol=ncol(df_)))
    colnames(df_empty)=colnames(df_)
    df_empty$Age_cut=AgeCut_not_there
    tryCatch({df_empty$Celltype=Celltypes_in_df[i]}, error=function(msg) {print("Pass.")})
    df_full=data.table::rbindlist(list(df_, df_empty))
    df_full=df_full[match(levels(score_df$AgeCut), df_full$Age_cut), ]
    df_full_list[[i]]=df_full
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  df_full_combined=df_full_combined[, 3:ncol(df_full_combined)]
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[2]*100)/100), abs(floor(quantiles[4]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/2, 0, ranges/2), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-colMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/4, 0, ranges/4), c("#00007E", "white", "#7E0000"))}
  
  p_age=list()
  codes=""
  for (i in 1:length(df_full_list)) {
    df_=as.data.frame(df_full_list[[i]])
    rownames(df_)=df_$Age_cut
    df_to_draw=df_[,3:ncol(df_)]
    colnames(df_to_draw)=gsub("_UCell", "", colnames(df_to_draw))
    celltype_of_df=names(table(df_$Celltype))
    df_to_draw=t(df_to_draw)
    df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)
    
    if (i==1) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                name=paste0("UCell of ", celltype_as_LegendName),
                heatmap_legend_param=list(title_position="leftcenter-rot",
                                          title_gp=gpar(fontsize=8, fontface="bold")),
                column_title=NULL,
                column_names_gp=gpar(fontsize=8),
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=T,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
                width=ncol(df_to_draw)*unit(3,"mm"), 
                height=nrow(df_to_draw)*unit(3,"mm"))}
    if (i>=2) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                show_heatmap_legend=FALSE,
                column_title=NULL,
                column_names_gp=gpar(fontsize=8),
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=T,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
                width=ncol(df_to_draw)*unit(3,"mm"), 
                height=nrow(df_to_draw)*unit(3,"mm"))
    }
    codes=paste0(codes, " %v% p_age[[", i, "]]")
  }
  codes=gsub("^ \\%v\\% ","",codes)
  eval(parse(text=paste0("ht_list=", codes)))
  p_age_total=draw(ht_list, ht_gap=unit(1, "mm"))
  
  return(p_age_total)
}

### Score
# NK
NK_cell_scoring=list(maturation=c("ITGAL","CXCR1","CX3CR1","CMKLR1","FCGR3A"),
                     adhesion=c("ITGAL","CD2","CD58","ITGAM","ICAM1","NCAM1","B3GAT1","SELL"),
                     cytokine=c("IL2RA","IL2RB","IL2RG","IL7R","IL18R1","KIT","CCR5","CCR7","CXCR1","CXCR3","CXCR4","CX3CR1","CMKLR1"),
                     activation=c("CD244","CD59","CD226","KLRF1","SLAMF6","SLAMF7","FCGR3A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DL4","KIR2DS4","KIR2DS5","KIR3DS1","KLRC2","KLRK1","NCR3","NCR2","NCR1","CD96"),
                     cytotoxic=c("FAS","FASLG","CD40LG","TNFSF10","LAMP1","LAMP2","GNLY","NKG7","GZMB","PRF1"),
                     inhibition=c("KLRC1","KLRD1","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5","KIR3DL1","KIR3DL2","LILRB1","PDCD1","SIGLEC7","CD300A","TIGIT"))
NK_df=score_df_for_drawing(rdspath="~/Project_PBMCage/Tempt_RDS/Scoring_NK.rds", 
                           object=object, celltypeKEY_for_grepl="NK\\.", gensetlist=NK_cell_scoring)
NK_scoring_plots_forAge=geneset_scoring_forAge(dataframe=NK_df, celltype_as_LegendName="NK cells", scale=F)
NK_scoring_plots_forAgeCuts=geneset_scoring_forAgeCuts(dataframe=NK_df, celltype_as_LegendName="NK cells", scale=F)
# CD8T and NKT
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB","AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","TNF","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
CD8T.NKT_df=score_df_for_drawing(rdspath="~/Project_PBMCage/Tempt_RDS/Scoring_CD8NKT.rds",
                                 object=object, celltypeKEY_for_grepl="(CD8)|(NKT)", gensetlist=T_cell_scoring)
CD8T.NKT_scoring_plots_forAge=geneset_scoring_forAge(dataframe=CD8T.NKT_df, celltype_as_LegendName="CD8 T and NKT cells", scale=F)
CD8T.NKT_scoring_plots_forAgeCuts=geneset_scoring_forAgeCuts(dataframe=CD8T.NKT_df, celltype_as_LegendName="CD8 T and NKT cells", scale=F)

pdf("~/Project_PBMCage/Plots/Scoring_CD8_and_NK_rowmeans.pdf", width=15, height=25)
NK_scoring_plots_forAge
NK_scoring_plots_forAgeCuts
CD8T.NKT_scoring_plots_forAge
CD8T.NKT_scoring_plots_forAgeCuts
dev.off()

# pdf("~/Project_PBMCage/Plots/Scoring_CD8_and_NK_scale.pdf", width=15, height=25)
# NK_scoring_plots_forAge
# NK_scoring_plots_forAgeCuts
# CD8T.NKT_scoring_plots_forAge
# CD8T.NKT_scoring_plots_forAgeCuts
# dev.off()

#####################################





### Identify the subsets of other celltypes
#####################################
###

### Identify the subsets of Bnaive
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Bnaive=subset(object, reassigned %in% c("Bnaive"))

library(dplyr)
library(Seurat)
Bnaive=Bnaive %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Bnaive_try=FindClusters2(Bnaive, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(Bnaive_try, label=T, label.size=5, repel=T) + NoLegend()
table(Bnaive_try$seurat_clusters)
Bnaive_try$seurat_clusters[Bnaive_try$seurat_clusters==7]=0
levels(Bnaive_try$seurat_clusters)=c("0","1","2","3","4","5","6","0")
Idents(Bnaive_try)="seurat_clusters"

VlnPlot(Bnaive_try, features=c("BACH2","BANK1","BLK","BTLA","CD79A","CD79B","FCRL1","FCRL3","HVCN1","RALGPS2"), pt.size=0)
VlnPlot(Bnaive_try, features=c("IGHM","IGHD","MME","CD19","CD24","CD27","CD38"), pt.size=0)
# plot Bnaive cell subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
Bnaive_try=AddModuleScore_UCell(Bnaive_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Bnaive", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Bnaive_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Bnaive_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: common B naive cells
VlnPlot(Bnaive_try, features=c("IGHM","IGHD","MS4A1","TCL1A","IL4R"), pt.size=0)
# cluster1: IgMhi transitional B cells with MZB (marginal zone B) genes which will later develop into MZB cells
VlnPlot(Bnaive_try, features=c("IGHM","PLD4","MZB1"), pt.size=0)
# cluster2: B intermediate (expressing naive markers - IL4R+TCL1A+, memory markers - CD24+CD27+, unswithced - IgMhi, switched - IgG+IgA+)
VlnPlot(CD4T_obj_subset, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","IGHG3"), pt.size=0)
# cluster3: B cells (express CD27 than Bnaive) with lineage infidelity to T cells (express CD3 but not CD4, CD8, or any marker of other T cells)
VlnPlot(Bnaive_try, features=c("CD79A","CD79B","IGHM","IGHD","CD3G","TRBC1","CD2","LAT"), pt.size=0)
# cluster4: B cells (IgD lower than Bnaive) with lineage infidelity to NKT cells
VlnPlot(Bnaive_try, features=c("CD79A","CD79B","IGHM","IGHD","KLRD1","FGFBP2","KLRF1","TRBC1"), pt.size=0)
# cluster5: IgM+IgD+CD27+ Bmem
VlnPlot(Bnaive_try, features=c("CD19","IGHD","CD27","IGHM","CD24"), pt.size=0)
# cluster6: pDCs
VlnPlot(Bnaive_try, features=c("EPHB1","SCN9A","PTCRA","SMPD3"), pt.size=0)

highlighted=WhichCells(Bnaive_try, idents=0)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(B_naive=c("IGHM","IGHD","MS4A1","TCL1A","IL4R"),
               IgMhi_transB=c("IGHM","PLD4","MZB1"),
               B_interm=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","IGHG3"),
               B_infidel_T=c("CD79A","CD79B","IGHM","IGHD","CD3G","TRBC1","CD2","LAT"),
               B_infidel_NKT=c("CD79A","CD79B","IGHM","IGHD","KLRD1","FGFBP2","KLRF1","TRBC1"),
               IgDdim_IgG.Bmem=c("CD19","IGHD","CD27","IGHM","CD24"),
               pDCs=c("EPHB1","SCN9A","PTCRA","SMPD3"))
p2=
  jjDotPlot(Bnaive_try, 
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Bnaive")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=0), "Bnaive", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=1), "Btransitional.IGHMhi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=2), "Bintermediate.IL4R+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=3), "B.lin_infidel.TRBC1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=4), "B.lin_infidel.TRBC1+ KLRD1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=5), "Bmem.IGHM+ IGHD+ CD27+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive_try, idents=6), "DC.IL3RA+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of lambda+ Bnaive
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Bnaive.lambda=subset(object, reassigned %in% c("lambda+ Bnaive"))

library(dplyr)
library(Seurat)
Bnaive.lambda=Bnaive.lambda %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Bnaive.lambda_try=FindClusters2(Bnaive.lambda, cluster.range=3, by=0.1, res=0.2, verbose=T)
DimPlot(Bnaive.lambda_try, label=T, label.size=5, repel=T) + NoLegend()
Idents(Bnaive.lambda_try)="seurat_clusters"

VlnPlot(Bnaive.lambda_try, features=c("BACH2","BANK1","BLK","BTLA","CD79A","CD79B","FCRL1","FCRL3","HVCN1","RALGPS2"), pt.size=0)
VlnPlot(Bnaive.lambda_try, features=c("IGHM","IGHD","MME","CD19","CD24","CD27","CD38"), pt.size=0)
# plot Bnaive.lambda subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
Bnaive.lambda_try=AddModuleScore_UCell(Bnaive.lambda_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Bnaive.lambda", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Bnaive.lambda_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Bnaive.lambda_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: lambda B naive cells
VlnPlot(Bnaive.lambda_try, features=c("IGHM","IGHD","MS4A1","TCL1A","IL4R","IGLC6","IGLC7","IGLC3"), pt.size=0)
# cluster1: B cells at the intermediate status between naive (IGHD+IL4R+TCL1A+) and Bmem (CD27+IGHA1+IGHG2+), slightly closer to Bnaive compared to Cluster2
VlnPlot(Bnaive.lambda_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"), pt.size=0)
# cluster2: B cells at the intermediate status between naive (IGHD+TCL1A+) and Bmem (CD27+IGHA1+IGHG2+CD24+) but closer to Bmem
VlnPlot(Bnaive.lambda_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"), pt.size=0)

highlighted=WhichCells(Bnaive.lambda_try, idents=1)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(lambda.Bnaive=c("IGHM","IGHD","MS4A1","TCL1A","IL4R","IGLC6","IGLC7","IGLC3"),
               lambda.Bint_early=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"),
               lambda.Bint_late=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"))
p2=
  jjDotPlot(Bnaive.lambda_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Bnaive.lambda")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive.lambda_try, idents=0), "Bnaive.lambda+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive.lambda_try, idents=1), "Bintermediate.lambda+ IL4R+ CD24-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Bnaive.lambda_try, idents=2), "Bintermediate.lambda+ IL4R- CD24+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of Classical Bmem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Classical_Bmem=subset(object, reassigned %in% c("Classical Bmem"))

library(dplyr)
library(Seurat)
Classical_Bmem=Classical_Bmem %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Classical_Bmem_try=FindClusters2(Classical_Bmem, cluster.range=c(2,5), by=0.1, res=0.2, verbose=T)
DimPlot(Classical_Bmem_try, label=T, label.size=5, repel=T) + NoLegend()
table(Classical_Bmem_try$seurat_clusters)
Classical_Bmem_try$seurat_clusters[Classical_Bmem_try$seurat_clusters==3]=0
levels(Classical_Bmem_try$seurat_clusters)=c("0","1","2","0")
Idents(Classical_Bmem_try)="seurat_clusters"

VlnPlot(Classical_Bmem_try, features=c("BACH2","BANK1","BLK","BTLA","CD79A","CD79B","FCRL1","FCRL3","HVCN1","RALGPS2"), pt.size=0)
VlnPlot(Classical_Bmem_try, features=c("IGHM","IGHD","MME","CD19","CD24","CD27","CD38"), pt.size=0)
# plot Classical Bmem subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
Classical_Bmem_try=AddModuleScore_UCell(Classical_Bmem_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="Classical Bmem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(Classical_Bmem_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(Classical_Bmem_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
# cluster0: B intermediate cells transforming from B naive to classical Bmem (close to classical Bmem since IL4R-TCL1A-CD27dim IGHMdim IGHDdim)
VlnPlot(Classical_Bmem_try, features=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R"), pt.size=0)
# cluster1: classical Bmem (IgM-IgD-CD27+)
VlnPlot(B_subset_try, features=c("IGHM","IGHD","CD27","HOPX","S100A10","IGHG1","IGHG2","IGHG3","IGHA1"), pt.size=0)
# cluster2: B intermediate cells transforming from B naive to IgM+IgD+ Bmem (closer to IgM+IgD+ Bmem since CD27+ IGHMhi IGHD+ PLD4+ MZB1+)
VlnPlot(Classical_Bmem_try, features=c("IGHM","IGHD","CD27","PLD4","MZB1","CD1C"), pt.size=0)

highlighted=WhichCells(Classical_Bmem_try, idents=2)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Bint_to_cBmem=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R"),
               cBmem=c("IGHM","IGHD","CD27","HOPX","S100A10"),
               Bint_to_IgGBmem=c("IGHM","IGHD","CD27","PLD4","MZB1","CD1C"))
p2=
  jjDotPlot(Classical_Bmem_try, 
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Classical Bmem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Classical_Bmem_try, idents=0), "Bintermediate.IGHDlo", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Classical_Bmem_try, idents=1), "Bmem.IGHM- IGHD- CD27+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Classical_Bmem_try, idents=2), "Bintermediate.IGHMhi", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of DN Bmem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
DN_Bmem=subset(object, reassigned %in% c("DN Bmem"))

library(dplyr)
library(Seurat)
DN_Bmem=DN_Bmem %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
DN_Bmem_try=FindClusters2(DN_Bmem, cluster.range=c(2,5), by=0.1, res=0.2, verbose=T)
DimPlot(DN_Bmem_try, label=T, label.size=5, repel=T) + NoLegend()
table(DN_Bmem_try$seurat_clusters)
DN_Bmem_try$seurat_clusters[DN_Bmem_try$seurat_clusters==4]=0
levels(DN_Bmem_try$seurat_clusters)=c("0","1","2","3","0")
Idents(DN_Bmem_try)="seurat_clusters"

VlnPlot(DN_Bmem_try, features=c("BACH2","BANK1","BLK","BTLA","CD79A","CD79B","FCRL1","FCRL3","HVCN1","RALGPS2"), pt.size=0)
VlnPlot(DN_Bmem_try, features=c("IGHM","IGHD","MME","CD19","CD24","CD27","CD38"), pt.size=0)
# plot DN_Bmem subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
DN_Bmem_try=AddModuleScore_UCell(DN_Bmem_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="DN Bmem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(DN_Bmem_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(DN_Bmem_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

VlnPlot(DN_Bmem_try, features=c("CD19","CD27","CD69","FOS","FOSB","IGHM","IGHD","CD1C","CD24","PLD4","TNFRSF13B","GPR183","IGHG1","IGHA1","IGHA2",
                                "ANXA2","S100A10","S100A4"), pt.size=0)
VlnPlot(DN_Bmem_try, features=c("IGHM","IGHD","CD27","IGHG1","IGHG2","IGHG3"), pt.size=0)
# cluster0: B intermediate cells transforming from B naive to IgM+IgD+CD5-CD27- Bmem (close to Bmem since IL4R-TCL1A-CD27-IGHM+IGHD+)
VlnPlot(DN_Bmem_try, features=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R","CD5"), pt.size=0)
# cluster1: under-regulated IgMhi/dim IgD-CD27- Bmem, of which have been isotype-switched since IGHM downregulating, IGHA1+IGHG1+IGHG2+
VlnPlot(DN_Bmem_try, features=c("CD19","IGHM","IGHD","CD27"), pt.size=0)
# cluster2: activated IgM+IgD+CD5-CD27- Bmem (with activation markers and IgM and IgD synthesis)
VlnPlot(DN_Bmem_try, features=c("IGHM","IGHD","CD27","CD5","CD69","FOSB"), pt.size=0)
# cluster3: newly generated IgM-IgD+CD27+ Bmem cells
VlnPlot(DN_Bmem_try, features=c("IGHM","IGHD","IGHG1","IGHG2","IGHG3","CD27","RGS13","VPREB3"), pt.size=0)

highlighted=WhichCells(DN_Bmem_try, idents=3)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Bint_to_CD27negBmem=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R","CD5"),
               CD27negBmem_underCtr=c("CD19","IGHM","IGHD","CD27"),
               CD27negBmem_activat=c("IGHM","IGHD","CD27","CD5","CD69","FOSB"),
               IgMneg.IgDpos.CD27neg_Bmem=c("IGHM","IGHD","CD27"))
p2=
  jjDotPlot(DN_Bmem_try,
            scale=F,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("DN Bmem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(DN_Bmem_try, idents=0), "Bintermediate.CD27-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(DN_Bmem_try, idents=1), "Bmem.CD27-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(DN_Bmem_try, idents=2), "Bmem.CD27- CD69+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(DN_Bmem_try, idents=3), "Bmem.IGHM- IGHD+ CD27+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of IgDdim IgG Bmem and IgDhi IgG Bmem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
IgD_IgG_Bmem=subset(object, reassigned %in% c("IgDdim IgG Bmem", "IgDhi IgG Bmem"))

library(dplyr)
library(Seurat)
IgD_IgG_Bmem=IgD_IgG_Bmem %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
IgD_IgG_Bmem_try=FindClusters2(IgD_IgG_Bmem, cluster.range=2, by=0.1, res=0.2, verbose=T)
DimPlot(IgD_IgG_Bmem_try, label=T, label.size=5, repel=T) + NoLegend()
Idents(IgD_IgG_Bmem_try)="seurat_clusters"

VlnPlot(DN_Bmem_try, features=c("BACH2","BANK1","BLK","BTLA","CD79A","CD79B","FCRL1","FCRL3","HVCN1","RALGPS2"), pt.size=0)
VlnPlot(IgD_IgG_Bmem_try, features=c("IGHM","IGHD","MME","CD19","CD24","CD27","CD38","IGHG1","IGHG2","IGHG3"), pt.size=0)
# plot DN_Bmem subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
DN_Bmem_try=AddModuleScore_UCell(DN_Bmem_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="DN Bmem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(DN_Bmem_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(IgD_IgG_Bmem_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

VlnPlot(IgD_IgG_Bmem_try, features=c("CD19","CD27","CD69","FOS","FOSB","IGHM","IGHD","CD1C","CD24","PLD4","TNFRSF13B","GPR183","IGHG1","IGHA1","IGHA2",
                                "ANXA2","S100A10","S100A4"), pt.size=0)

# cluster0: (optional: IgM+) IgD+ CD27+ IgG+ Bmem
VlnPlot(IgD_IgG_Bmem_try, features=c("IGHM","IGHD","CD27","IGHG1","IGHG2","IGHG3"), pt.size=0)
# cluster1: (optional: IgM+) IgDdim CD27+ IgG+ Bmem
VlnPlot(IgD_IgG_Bmem_try, features=c("IGHM","IGHD","CD27","IGHG1","IGHG2","IGHG3"), pt.size=0)

highlighted=WhichCells(IgD_IgG_Bmem_try, idents=1)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(IgDpos_IgG_Bmem=c("IGHM","IGHD","CD27","IGHG1","IGHG2","IGHG3"),
               IgDlo_IgG_Bmem=c("IGHM","IGHD","CD27","IGHG1","IGHG2","IGHG3"))
p2=
  jjDotPlot(IgD_IgG_Bmem_try,
            scale=F,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("IgG Bmem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(IgD_IgG_Bmem_try, idents=0), "Bmem.IGHG1+ IGHD+ CD27+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(IgD_IgG_Bmem_try, idents=1), "Bmem.IGHG1+ IGHDlo CD27+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CD14 Mono
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD14_Mono=subset(object, reassigned %in% c("CD14 Mono"))

library(dplyr)
library(Seurat)
CD14_Mono=CD14_Mono %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD14_Mono_try=FindClusters2(CD14_Mono, cluster.range=c(3,6), by=0.1, res=0.2, verbose=T)
DimPlot(CD14_Mono_try, label=T, label.size=5, repel=T) + NoLegend()
Idents(CD14_Mono_try)="seurat_clusters"

VlnPlot(CD14_Mono_try, features=c("S100A9","VCAN","LYZ","RNASE2","CD163","CSF3R","CD14","NAIP","F13A1","S100A12","CD36","TREM1"), pt.size=0)
VlnPlot(CD14_Mono_try, features=c("FCGR3A","LST1","AIF1","CTSS","MTSS1","TCF7L2","IFITM3","MS4A7","LIRB2","CSF1R","IFITM2","HCK","CTSL1","MAFB","TNFRSF1B","SIGLEC10","FGR",'LIRA2',"NEAT1","RHOC","EMR2"), pt.size=0)
# plot CD14 Monocytes subsets markers
Monocyte_scoring=list(cytotoxic=c("PRF1","GNLY","CTSW","KLRD1","NKG7","IL2RB","GZMA","ZAP70","KLRF1","GZMH","IL32","IKZF3","LCK","CD96","TGFBR3","CD2","CD247"),
                      cycl_diff_traffic=c("G0S2","NAMPT","FCGR3B","SRGN","TNFRSF10C","MXD1","CXCR2","VNN2","CXCR1","FCGR2A","CLEC4E","LITAF"),
                      APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                      APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                      CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                      check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                      HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                      inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                      MHCI=c("B2M","HLA-A","TAP1"),
                      parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                      IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                      IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
CD14_Mono_try=AddModuleScore_UCell(CD14_Mono_try, features=Monocyte_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD14 Mono", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD14_Mono_try,
                        features=paste0(names(Monocyte_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(Monocyte_scoring)),
                        ylab=names(Monocyte_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD14_Mono_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(CD14_Mono_try, features=c("CD14","VIM","S100A8","S100A9","S100A12","SOD2","CYBA","NAMPT","EIF5A","C1orf56","ISG15","MX1","MX2",
                                  "HLA-DRA","HLA-DQA1","HLA-DPB1","CD74"), pt.size=0)

# cluster0: CD14 classical Mono
VlnPlot(CD14_Mono_try, features=c("S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14","FCGR3A"), pt.size=0)
# cluster1: CD14 classical Mono with some mature DC phenotypes
VlnPlot(Monocytes_subset_try, features=c("CD83","CD14","CCL17","CCL13","HLA-DPA1","HLA-DQA1","HLA-DRB1"), pt.size=0)
# cluster2: CD14 classical monocytes with lineage infidelity to T cells (maybe CD4 T)
VlnPlot(CD14_Mono_try, features=c("IL7R","CD3G","TRAC","CD2","TRBC1"), pt.size=0)
# cluster3: CD14++CD16+ Intermediate monocytes with cytotoxic properties
VlnPlot(Monocytes_subset_try, features=c("S100A9","VCAN","LYZ","FCGR3A","RHOC","KLRD1","KLRF1","GNLY"), pt.size=0)
# cluster4 CD14 classical monocytes with lineage infidelity to B cells
VlnPlot(CD14_Mono_try, features=c("S100A9","VCAN","LYZ","FCGR3A","JCHAIN","TCL1A","MS4A1","IGHD"), pt.size=0)

highlighted=WhichCells(CD14_Mono_try, idents=4)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(classical_Mono=c("S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14","FCGR3A"),
               cMono_with_DCpheno=c("CD83","CD14","CCL17","CCL13","HLA-DPA1","HLA-DQA1","HLA-DRB1"),
               cMono_infide_toT=c("IL7R","CD3G","TRAC","CD2","TRBC1"),
               IntMono_cytotoxic=c("S100A9","VCAN","LYZ","FCGR3A","RHOC","KLRD1","KLRF1","GNLY"),
               cMono_infide_toB=c("S100A9","VCAN","LYZ","FCGR3A","JCHAIN","TCL1A","MS4A1","IGHD"))
p2=
  jjDotPlot(CD14_Mono_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD14 Mono")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD14_Mono_try, idents=0), "Mono.CD14+ CD16-", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD14_Mono_try, idents=1), "Mono.CD14+ CD16- HLA+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD14_Mono_try, idents=2), "Mono.CD14.lin_infidel.TRBC1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD14_Mono_try, idents=3), "Mono.CD14+ CD16+ GNLY+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD14_Mono_try, idents=4), "Mono.CD14.lin_infidel.MS4A1+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CD16 Mono
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD16_Mono=subset(object, reassigned %in% c("CD16 Mono"))

library(dplyr)
library(Seurat)
CD16_Mono=CD16_Mono %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD16_Mono_try=FindClusters2(CD16_Mono, cluster.range=c(3,6), by=0.1, res=0.2, verbose=T)
DimPlot(CD16_Mono_try, label=T, label.size=5, repel=T) + NoLegend()
Idents(CD16_Mono_try)="seurat_clusters"

VlnPlot(CD16_Mono_try, features=c("S100A9","VCAN","LYZ","RNASE2","CD163","CSF3R","CD14","NAIP","F13A1","S100A12","CD36","TREM1"), pt.size=0)
VlnPlot(CD16_Mono_try, features=c("FCGR3A","LST1","AIF1","CTSS","MTSS1","TCF7L2","IFITM3","MS4A7","LIRB2","CSF1R","IFITM2","HCK","CTSL1","MAFB","TNFRSF1B","SIGLEC10","FGR",'LIRA2',"NEAT1","RHOC","EMR2"), pt.size=0)

# plot CD16 Monocytes subsets markers
Monocyte_scoring=list(cytotoxic=c("PRF1","GNLY","CTSW","KLRD1","NKG7","IL2RB","GZMA","ZAP70","KLRF1","GZMH","IL32","IKZF3","LCK","CD96","TGFBR3","CD2","CD247"),
                      cycl_diff_traffic=c("G0S2","NAMPT","FCGR3B","SRGN","TNFRSF10C","MXD1","CXCR2","VNN2","CXCR1","FCGR2A","CLEC4E","LITAF"),
                      APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                      APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                      CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                      check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                      HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                      inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                      MHCI=c("B2M","HLA-A","TAP1"),
                      parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                      IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                      IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
CD16_Mono_try=AddModuleScore_UCell(CD16_Mono_try, features=Monocyte_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD16 Mono", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD16_Mono_try,
                        features=paste0(names(Monocyte_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(Monocyte_scoring)),
                        ylab=names(Monocyte_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD16_Mono_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","MS4A7","CX3CR1","CLEC10A","FCER1A","CST3","CD74","C1QA","C1QB","C1QC","CD68"), pt.size=0)

# cluster0: CD14-CD16+ nonclassical Mono
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","LST1","AIF1","CTSS","CTSL"), pt.size=0)
# cluster1: CD14-CD16+ nonclassical Mono with cell proliferation property (EIF5A and C1orf56)
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","EIF5A","C1orf56"), pt.size=0)
# cluster2: CD14-CD16+ nonclassical monocytes with cytotoxic properties
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","KLRD1","KLRF1","GNLY"), pt.size=0)
# cluster3: CD14-CD16+ nonclassical monocytes with lineage infidelity to T cells (maybe CD4 T)
VlnPlot(CD16_Mono_try, features=c("IL7R","CD3G","TRAC","CD2","TRBC1"), pt.size=0)
# cluster4: CD14+CD16+ Intermediate monocytes with C1qhi as a proinflammatory cytokine-secreting and enhanced phagocytotic marker
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","C1QA","C1QB","C1QC"), pt.size=0)
# cluster5 :CD14-CD16+ nonclassical monocytes with lineage infidelity to B cells
VlnPlot(CD16_Mono_try, features=c("CD14","FCGR3A","JCHAIN","TCL1A","MS4A1","IGHD"), pt.size=0)

highlighted=WhichCells(CD16_Mono_try, idents=5)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(nonclassical_Mono=c("CD14","FCGR3A","LST1","AIF1","CTSS","CTSL"),
               nMono_prolif=c("CD14","FCGR3A","EIF5A","C1orf56"),
               nMono_cytotoxic=c("CD14","FCGR3A","KLRD1","KLRF1","GNLY"),
               nMono_infide_toT=c("IL7R","CD3G","TRAC","CD2","TRBC1"),
               IntMono_C1qhi=c("CD14","FCGR3A","C1QA","C1QB","C1QC"),
               nMono_infide_toB=c("CD14","FCGR3A","JCHAIN","TCL1A","MS4A1","IGHD"))
p2=
  jjDotPlot(CD16_Mono_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD16 Mono")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=0), "Mono.CD14- CD16+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=1), "Mono.CD14- CD16+ C1orf56+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=2), "Mono.CD14- CD16+ GNLY+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=3), "Mono.CD16.lin_infidel.TRBC1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=4), "Mono.CD14+ CD16+ C1Q+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD16_Mono_try, idents=5), "Mono.CD16.lin_infidel.MS4A1+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of IgG+ PC, IgA+ PC, PB
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
PB.PC=subset(object, reassigned %in% c("IgG+ PC", "IgA+ PC", "PB") | CD8T_NK_subsets=="IgA+ PC")

library(dplyr)
library(Seurat)
PB.PC=PB.PC %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
PB.PC_try=FindClusters2(PB.PC, cluster.range=c(3,4), by=0.1, res=0.1, verbose=T)
DimPlot(PB.PC_try, label=T, label.size=5, repel=T) + NoLegend()
table(PB.PC_try$seurat_clusters)
Idents(PB.PC_try)="seurat_clusters"

# plot PB and PC subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
PB.PC_try=AddModuleScore_UCell(PB.PC_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="PB and PC", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(PB.PC_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(PB.PC_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

# cluster0: plasma cells
VlnPlot(PB.PC_try, features=c("CD38","IRF4","PRDM1","JCHAIN","XBP1"), pt.size=0)
# cluster1: plasma cells with lineage infidelity to NKT cells
VlnPlot(PB.PC_try, features=c("CD3G","TRAC","CD2","TRBC1","KLRD1","FGFBP2"), pt.size=0)
# cluster2: plasmablasts
VlnPlot(PB.PC_try, features=c("CXCR3","HLA-DRA","MKI67","CXCR4"), pt.size=0)

highlighted=WhichCells(PB.PC_try, idents=2)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(PC=c("CD38","IRF4","PRDM1","JCHAIN","XBP1"),
               PC_infide_toNKT=c("CD3G","TRAC","CD2","TRBC1","KLRD1","FGFBP2"),
               PB=c("CXCR3","HLA-DRA","MKI67","CXCR4"))
p2=
  jjDotPlot(PB.PC_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("PB and PC")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(PB.PC_try, idents=0), "Bpc", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(PB.PC_try, idents=1), "Bpc.lin_infidel.TRBC1+ KLRD1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(PB.PC_try, idents=2), "Bpb", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of pDC and cDC2
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
pDC_cDC2=subset(object, reassigned %in% c("cDC2", "pDC"))

library(dplyr)
library(Seurat)
pDC_cDC2=pDC_cDC2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
pDC_cDC2_try=FindClusters2(pDC_cDC2, cluster.range=c(4,6), by=0.1, res=0.4, verbose=T)
DimPlot(pDC_cDC2_try, label=T, label.size=5, repel=T) + NoLegend()
table(pDC_cDC2_try$seurat_clusters)
Idents(pDC_cDC2_try)="seurat_clusters"

# plot pDC and cDC2 subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
pDC_cDC2_try=AddModuleScore_UCell(pDC_cDC2_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="pDC and cDC2", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(pDC_cDC2_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(pDC_cDC2_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(pDC_cDC2_try, features=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","CLEC9A","THBD"), pt.size=0)

# cluster0: cDC2
VlnPlot(pDC_cDC2_try, features=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO"), pt.size=0)
# cluster1: pDC
VlnPlot(pDC_cDC2_try, features=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","CLEC4C"), pt.size=0)
# cluster2: cDC2 with lineage infidelity to NKT cells
VlnPlot(pDC_cDC2_try, features=c("CD3G","TRAC","CD2","TRBC1","KLRD1","FGFBP2"), pt.size=0)
# cluster3: pDC with lineage infidelity to T cells
VlnPlot(pDC_cDC2_try, features=c("CD3G","TRBC1","LCK","TRAC"), pt.size=0)
# cluster4: cDC2 with lineage infidelity to B cells
VlnPlot(pDC_cDC2_try, features=c("PAX5","MS4A1","IGHD"), pt.size=0)
# cluster5: cDC2 with cell proliferative property
VlnPlot(pDC_cDC2_try, features=c("PCLAF","TYMS","MKI67"), pt.size=0)

highlighted=WhichCells(pDC_cDC2_try, idents=5)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(cDC2=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO"),
               pDC=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2"),
               cDC2_infide_toNKT=c("CD3G","TRAC","CD2","TRBC1","KLRD1","FGFBP2"),
               pDC_infide_toT=c("CD3G","TRBC1","LCK","TRAC"),
               cDC2_infide_toB=c("PAX5","MS4A1","IGHD"),
               cDC2_prolif=c("PCLAF","TYMS","MKI67"))
p2=
  jjDotPlot(pDC_cDC2_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("pDC and cDC2")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=0), "DC.CD1C+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=1), "DC.IL3RA+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=2), "DC.CD1C.lin_infidel.TRBC1+ KLRD1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=3), "DC.IL3RA.lin_infidel.TRBC1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=4), "DC.CD1C.lin_infidel.MS4A1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(pDC_cDC2_try, idents=5), "DC.CD1C+ MKI67+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of AXL+ DC and cDC1
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
aDC_cDC1=subset(object, reassigned %in% c("AXL+ DC", "cDC1"))

library(dplyr)
library(Seurat)
aDC_cDC1=aDC_cDC1 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
aDC_cDC1_try=FindClusters2(aDC_cDC1, cluster.range=3, by=0.1, res=0.4, verbose=T)
DimPlot(aDC_cDC1_try, label=T, label.size=5, repel=T) + NoLegend()
table(aDC_cDC1_try$seurat_clusters)
Idents(aDC_cDC1_try)="seurat_clusters"

# plot AXL+ DC and cDC1 subsets markers
General_cell_scoring=list(APC_coinhibition=c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"),
                          APC_costimu=c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"),
                          CCR=c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14","CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8","IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG","IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2","CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2","IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP","CCL3L1","TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2","CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A","IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"),
                          check_point=c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244","CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"),
                          cytolytic=c("PRF1","GZMA"),
                          HLA=c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"),
                          inflam_promoting=c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"),
                          MHCI=c("B2M","HLA-A","TAP1"),
                          parainflammation=c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1","CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3","RETNLB","IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1","IFITM3","IL1RN"),
                          T_coinhibition=c("BTLA","C10orf54","CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"),
                          T_costimu=c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14"),
                          IFN_I_response=c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"),
                          IFN_II_response=c("GPR146","SELP","AHR"))
library(UCell)
aDC_cDC1_try=AddModuleScore_UCell(aDC_cDC1_try, features=General_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="AXL+ DC and cDC1", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(aDC_cDC1_try,
                        features=paste0(names(General_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(General_cell_scoring)),
                        ylab=names(General_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(aDC_cDC1_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(aDC_cDC1_try, features=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","CLEC9A","THBD","HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT",
                                 "ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2"), pt.size=0)

# cluster0: AXL+ mDC
VlnPlot(aDC_cDC1_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","FCER1A","CLEC10A"), pt.size=0)
# cluster1: cDC1
VlnPlot(aDC_cDC1_try, features=c("HLA-DRB1","HLA-DRA","THBD","CLEC9A","CADM1","XCR1","BATF3"), pt.size=0)
# cluster2: AXL+ pDC
VlnPlot(aDC_cDC1_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","TPM2","LRRC26"), pt.size=0)


highlighted=WhichCells(aDC_cDC1_try, idents=1)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(AXL.mDC=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","FCER1A","CLEC10A"),
               cDC1=c("HLA-DRB1","HLA-DRA","THBD","CLEC9A","CADM1","XCR1","BATF3"),
               AXL.pDC=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","TPM2","LRRC26"))
p2=
  jjDotPlot(aDC_cDC1_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("AXL+ DC and cDC1")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(aDC_cDC1_try, idents=0), "DC.AXL+ CLEC10A+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(aDC_cDC1_try, idents=1), "DC.THBD+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(aDC_cDC1_try, idents=2), "DC.AXL+ TPM2+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CD4 Tcm
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD4_Tcm=subset(object, reassigned %in% c("CD4 Tcm"))

library(dplyr)
library(Seurat)
CD4_Tcm=CD4_Tcm %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD4_Tcm_try=FindClusters2(CD4_Tcm, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4_Tcm_try, label=T, label.size=5, repel=T) + NoLegend()
table(CD4_Tcm_try$seurat_clusters)
CD4_Tcm_try$seurat_clusters[CD4_Tcm_try$seurat_clusters==5]=0
levels(CD4_Tcm_try$seurat_clusters)=c("0","1","2","3","4","0")
Idents(CD4_Tcm_try)="seurat_clusters"

VlnPlot(CD4_Tcm_try, features=c("CCR4","CCR6","CXCR5","CXCR3","KLRB1","GZMK","RORC","AHR","CCR10","GATA3","HLA-DRA","HLA-DRB1","CTLA4","PDCD1","LAG3","HAVCR2"), pt.size=0)
VlnPlot(CD4_Tcm_try, features=c("CCR7","PRF1","GNLY","FOXP3","EOMES","TBX21","GZMB","GZMK","CCR7","IL2RA","IL7R"), pt.size=0)
VlnPlot(CD4_Tcm_try, features=c("CD44","SELL","IL7R"), pt.size=0)

# plot CD4 Tcm subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD4_Tcm_try=AddModuleScore_UCell(CD4_Tcm_try, features=T_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD4 Tcm", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD4_Tcm_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD4_Tcm_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(CD4_Tcm_try, features=c("IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL"), pt.size=0)

# cluster0: CD4 Tcm
VlnPlot(CD4_subset_try, features=c("CCR7","SELL","GPR183","PIK3IP1"), pt.size=0)
# cluster1: CD4 Tcm with cell proliferative property
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster2: CD4 Tcm TH17-like
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","KLRB1","CTSH","AQP3","TIMP1"), pt.size=0)
# cluster3: CD4 Tem
VlnPlot(CD4_Tcm_try, features=c("CD44","SELL","IL7R","CCL5","FYB1"), pt.size=0)
# cluster4: CD8 Tcm with Th2 transcriptional signature, i.e., type2/Th2 memory
VlnPlot(CD8_subset_try, features=c("CD8A","CD8B","CD44","SELL","IL7R","GATA3","ANXA1","XBP1","LGALS1","MALAT1","GZMB","LYAR"), pt.size=0)


highlighted=WhichCells(CD4_Tcm_try, idents=4)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Tcm=c("CCR7","SELL","GPR183","PIK3IP1"),
               Tcm_prolif=c("CD44","SELL","IL7R","C1orf56","HNRNPH1","CDC42SE1"),
               Tcm_Th17like=c("CD44","SELL","IL7R","KLRB1","CTSH","AQP3","TIMP1"),
               Tem=c("CD44","SELL","IL7R","CCL5","FYB1"),
               CD8_Tcm_Th2like=c("CD8A","CD8B","CD44","SELL","IL7R","GATA3","ANXA1","XBP1"))
p2=
  jjDotPlot(CD4_Tcm_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD4 Tcm")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tcm_try, idents=0), "CD4 Tcm.PIK3IP1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tcm_try, idents=1), "CD4 Tcm.C1orf56+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tcm_try, idents=2), "CD4 Tcm.KLRB1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tcm_try, idents=3), "CD4 Tem", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tcm_try, idents=4), "CD8 Tcm.GATA3+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CD4 Tem
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
object=AddMetaData(object, obj_annotation$CD8T_NK_subsets, col.name="CD8T_NK_subsets2")
CD4_Tem=subset(object, CD8T_NK_subsets2 %in% c("CD4 Tem"))

library(dplyr)
library(Seurat)
CD4_Tem=CD4_Tem %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CD4_Tem_try=FindClusters2(CD4_Tem, cluster.range=c(5,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4_Tem_try, label=T, label.size=5, repel=T) + NoLegend()
table(CD4_Tem_try$seurat_clusters)
Idents(CD4_Tem_try)="seurat_clusters"

VlnPlot(CD4_Tem_try, features=c("CCR4","CCR6","CXCR5","CXCR3","KLRB1","GZMK","RORC","AHR","CCR10","GATA3","HLA-DRA","HLA-DRB1","CTLA4","PDCD1","LAG3","HAVCR2"), pt.size=0)
VlnPlot(CD4_Tem_try, features=c("CCR7","PRF1","GNLY","FOXP3","EOMES","TBX21","GZMB","GZMK","CCR7","IL2RA","IL7R"), pt.size=0)
VlnPlot(CD4_Tem_try, features=c("CD44","SELL","IL7R"), pt.size=0)

# plot CD4 Tem subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD4_Tem_try=AddModuleScore_UCell(CD4_Tem_try, features=T_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD4 Tem", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD4_Tem_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD4_Tem_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(CD4_Tem_try, features=c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL"), pt.size=0)

# cluster0: CD4 Tem with Th17 transcriptional signature (KLRB1+ GZMK+)
VlnPlot(CD4_Tem_try, features=c("CD44","SELL","IL7R","GZMK","KLRB1"), pt.size=0)
# cluster1: CD4 Tcm under cell regulation (events eg. proliferation/cytokine production/apoptosis...)
VlnPlot(CD4_Tem_try, features=c("CD44","SELL","IL7R","COTL1","ITGB1","LGALS1","S100A11","ANXA1"), pt.size=0)
# cluster2: CD4 Tem with transcriptional regulation property
VlnPlot(CD4_Tem_try, features=c("JUND","MALAT1","DDX17"), pt.size=0)
# cluster3: CD4 Tcm with upregulated HLA-DR
VlnPlot(CD4_Tem_try, features=c("CD44","SELL","IL7R","HLA-DRA","CD38","IL2RA"), pt.size=0)
# cluster4: ILCreg
VlnPlot(CD4_Tem_try, features=c("IL7R","SOX4","ID3"), pt.size=0)


highlighted=WhichCells(CD4_Tem_try, idents=4)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Tem_Th17=c("CD44","SELL","IL7R","GZMK","KLRB1"),
               Tem_cellreg=c("CD44","SELL","IL7R","COTL1","ITGB1","LGALS1","S100A11","ANXA1"),
               Tem_transreg=c("JUND","MALAT1","DDX17"),
               Tcm_HLA=c("CD44","SELL","IL7R","HLA-DRA","CD38","IL2RA"),
               ILCreg=c("IL7R","SOX4","ID3"))
p2=
  jjDotPlot(CD4_Tem_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD4 Tem")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tem_try, idents=0), "CD4 Tem.KLRB1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tem_try, idents=1), "CD4 Tcm.LGALS1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tem_try, idents=2), "CD4 Tem.MALAT1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tem_try, idents=3), "CD4 Tcm.HLA+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tem_try, idents=4), "ILC.SOX4+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CD4 Tnaive
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CD4_Tnaive=subset(object, CD8T_NK_subsets %in% c("CD4 Tnaive"))

library(dplyr)
library(Seurat)
CD4_Tnaive=CD4_Tnaive %>% 
  NormalizeData() %>%
  FindVariableFeatures()

CD4_Tnaive=SketchData(CD4_Tnaive, ncells=50000, method="LeverageScore", sketched.assay="sketch")
DefaultAssay(CD4_Tnaive)="sketch"
CD4_Tnaive=CD4_Tnaive %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:20) %>%
  RunUMAP(dims=1:20, return.model=T)
DimPlot(CD4_Tnaive)
source("~/Rscripts/FindCluster2_Functions.R")
CD4_Tnaive_try=FindClusters2(CD4_Tnaive, cluster.range=c(3,5), by=0.1, res=0.01, verbose=T)
CD4_Tnaive_try=FindClusters(CD4_Tnaive, resolution=0.2, verbose=T)
DimPlot(CD4_Tnaive_try, label=T, label.size=5, repel=T) + NoLegend()

CD4_Tnaive_try=ProjectData(
  CD4_Tnaive_try,
  assay="RNA",
  full.reduction="pca.full",
  sketched.assay="sketch",
  sketched.reduction="pca",
  umap.model="umap",
  dims=1:50,
  refdata=list(seurat_clusters_full="seurat_clusters")
)

DefaultAssay(CD4_Tnaive_try)="RNA"
Idents(CD4_Tnaive_try)="seurat_clusters_full"
DimPlot(CD4_Tnaive_try, label=T, label.size=5, repel=T) + NoLegend()


VlnPlot(CD4_Tnaive_try, features=c("CCR4","CCR6","CXCR5","CXCR3","KLRB1","GZMK","RORC","AHR","CCR10","GATA3","HLA-DRA","HLA-DRB1","CTLA4","PDCD1","LAG3","HAVCR2"), pt.size=0)
VlnPlot(CD4_Tnaive_try, features=c("CCR7","PRF1","GNLY","FOXP3","EOMES","TBX21","GZMB","GZMK","CCR7","IL2RA","IL7R"), pt.size=0)
VlnPlot(CD4_Tnaive_try, features=c("CD44","SELL","IL7R"), pt.size=0)

# plot CD4 Tnaive subsets markers
T_cell_scoring=list(naive=c("TCF7","SELL","LEF1","CCR7"),
                    activation=c("COTL1","TYROBP","IL2RB",
                                 "AKT1","AKT2","AKT3","ARAF","BRAF","CD247","CD28","CD3D","CD3E","CD3G","CD80","CD86","CDC42","CHUK","CSK","FOS","GRAP2","HRAS","IKBKB","ITPR1","JUN","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAP3K1","MAPK1","MAPK3","MAPK8","MAPK9","NCK1","NCK2","NFKBIA","NRAS","PAK1","PAK2","PAK3","PIK3C3","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PLCG1","PPP3CA","PPP3CB","PPP3CC","PRKCQ","PTPRC","RAC1","RAF1","SOS1","SOS2","VAV1","VAV2","VAV3","WAS","ZAP70"),
                    migration=c("CCR5","XCL1","CXCR3","CXCR5","CCR7-","GPR183"),
                    cytotoxic=c("GZMK","IFNG","GZMH","GNLY","PRF1","NKG7","GZMA","GZMB"),
                    inhibitory=c("LAG3","PDCD1","CTLA4","TIGIT","BTLA","HAVCR2"),
                    effector=c("IL2","GZMA","GNLY","PRF1","GZMB","GZMK","IFNG","NKG7"),
                    cytokine=c("IL2","TNF"),
                    regulatory=c("IL2RA","FOXP3","IKZF2"),
                    proliferating=c("MKI67","CDK1","STMN1"),
                    transcription=c("EOMES","HOPX","TBX21","ZEB2","ZNF683","HIF1A","ID2","TOX"),
                    costimulatory=c("CD28","TNFRSF14","ICOS","TNFRSF9"))
library(UCell)
CD4_Tnaive_try=AddModuleScore_UCell(CD4_Tnaive_try, features=T_cell_scoring)

library(ggplot2)
library(patchwork)
p_title=
  ggplot() +
  theme_void() +
  geom_text(aes(0,0, label="CD4 Tnaive", fontface="bold"), size=6) +
  theme(plot.margin=unit(c(0,1,3,1), "mm")) +
  xlab(NULL)
p1=
  p_title /
  SCpubr::do_ViolinPlot(CD4_Tnaive_try,
                        features=paste0(names(T_cell_scoring),"_UCell"),
                        plot_boxplot=FALSE,
                        xlab=rep("", length(T_cell_scoring)),
                        ylab=names(T_cell_scoring),
                        axis.title.face="plain",
                        axis.text.x.angle=0) +
  plot_layout(heights=unit(c(1,10), "null"))

# find overexpressed genes
markers=FindAllMarkers(CD4_Tnaive_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(CD4_Tnaive_try, features=c("CD3D","CD3E","CD3G","CD8A","CD8B","CD4","IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL"), pt.size=0)

# cluster0: CD4 Tnaive with proliferative property (under transformation to Tcm as CD44+ SELL+)
VlnPlot(CD4_subset_try, features=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster1: CD4 Tnaive under cell regulation of events including activation and proliferation
VlnPlot(CD4_subset_try, features=c("CD44","SELL","IL7R","COTL1","SOCS3","ANXA1","TAGLN2"), pt.size=0)
# cluster2: myeloid-derived DCs (express both Mono markers and DC markers)
VlnPlot(CD4_Tnaive_try, features=c("S100A9","LYZ","S100A8","CD1C"), pt.size=0)
# cluster3: CD4 Tnaive under cell regulation of quiescence (repressing differentiation), reported to be downreg during aging
VlnPlot(CD4_Tnaive_try, features=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"), pt.size=0)
# cluster4: myeloid-derived DCs (express both Mono markers and DC markers)
VlnPlot(CD4_Tnaive_try, features=c("S100A9","LYZ","S100A8","CD1C"), pt.size=0)


highlighted=WhichCells(CD4_Tnaive_try, idents=4)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Tnaive_prolif=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"),
               Tnaive_activate=c("CD44","SELL","IL7R","COTL1","SOCS3","ANXA1","TAGLN2"),
               mDC=c("S100A9","LYZ","S100A8","CD1C"),
               Tnaive_quies=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"),
               mDC2=c("S100A9","LYZ","S100A8","CD1C"))
p2=
  jjDotPlot(CD4_Tnaive_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CD4 Tnaive")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tnaive_try, idents=0), "CD4 Tnaive.C1orf56+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tnaive_try, idents=1), "CD4 Tnaive.SOCS3+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tnaive_try, idents=2), "DC.LYZ+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tnaive_try, idents=3), "CD4 Tnaive.SOX4+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CD4_Tnaive_try, idents=4), "DC.LYZ+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of CLP
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
CLP=subset(object, reassigned %in% c("CLP"))

library(dplyr)
library(Seurat)
CLP=CLP %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
CLP_try=FindClusters2(CLP, cluster.range=c(2,6), by=0.1, res=0.3, verbose=T)
DimPlot(CLP_try, label=T, label.size=5, repel=T) + NoLegend()
table(CLP_try$seurat_clusters)
Idents(CLP_try)="seurat_clusters"

# find overexpressed genes
markers=FindAllMarkers(CLP_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

# cluster0: CLP
VlnPlot(CLP_try, features=c("DNTT","PRSS2","EGFL7","CYTL1","CD34"), pt.size=0)
# cluster1: HSC
VlnPlot(CLP_try, features=c("SAMHD1","LINC00299","TESC","ITGB7","COTL1"), pt.size=0)


highlighted=WhichCells(CLP_try, idents=1)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(CLP=c("DNTT","PRSS2","EGFL7","CYTL1","CD34"),
               HSC=c("SAMHD1","LINC00299","TESC","ITGB7","COTL1"))
p2=
  jjDotPlot(CLP_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("CLP")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CLP_try, idents=0), "Progenitor.CLP", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(CLP_try, idents=1), "Progenitor.HSC", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of EMP
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
EMP=subset(object, reassigned %in% c("EMP"))

library(dplyr)
library(Seurat)
EMP=EMP %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
EMP_try=FindClusters2(EMP, cluster.range=c(2,6), by=0.1, res=0.3, verbose=T)
DimPlot(EMP_try, label=T, label.size=5, repel=T) + NoLegend()
table(EMP_try$seurat_clusters)
Idents(CLP_try)="seurat_clusters"

# find overexpressed genes
markers=FindAllMarkers(EMP_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

# cluster0: previously reported as presumptive multi-lineage progenitor, may actually be MPP
VlnPlot(EMP_try, features=c("EGFL7","CYTL1","CD34","CSF3R","SMIM24","C1QTNF4"), pt.size=0)
# cluster1: MEP (with genes reported to be Ery-lineage progenitor gene signatures/the first row; and those to be megakaryocytic-erythroid progenitors/the 2nd row)
VlnPlot(EMP_try, features=c("EPOR","KLF1","CSF2RB","APOC1","CNRIP1",
                            "CD44","TFRC","ITGA2B"), pt.size=0)
# cluster2: NK cell progenitors (NKPs in short)
VlnPlot(EMP_try, features=c("IL7R","ITGA4","CD27","CD44","IL2RB","CD7"), pt.size=0)

highlighted=WhichCells(EMP_try, idents=1)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(MPP=c("EGFL7","CYTL1","CD34","CSF3R","SMIM24","C1QTNF4"),
               EMP=c("EPOR","KLF1","CSF2RB","APOC1","CNRIP1","TFRC","ITGA2B"),
               NKP=c("IL7R","ITGA4","CD27","CD44","IL2RB","CD7"))
p2=
  jjDotPlot(EMP_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("EMP")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(EMP_try, idents=0), "Progenitor.MPP", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(EMP_try, idents=1), "Progenitor.MEP", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(EMP_try, idents=2), "Progenitor.NKP", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of Platelet and Mast cell
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Mast_plat=subset(object, reassigned %in% c("Mast cell","Platelet"))

library(dplyr)
library(Seurat)
Mast_plat=Mast_plat %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)

source("~/Rscripts/FindCluster2_Functions.R")
Mast_plat_try=FindClusters2(Mast_plat, cluster.range=c(10,20), by=0.1, res=1, verbose=T)
DimPlot(Mast_plat_try, label=T, label.size=5, repel=T) + NoLegend()
table(Mast_plat_try$seurat_clusters)
Idents(Mast_plat_try)="seurat_clusters"

Mast_plat2=subset(Mast_plat_try, seurat_clusters==1)
library(dplyr)
library(Seurat)
Mast_plat2=Mast_plat2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)
source("~/Rscripts/FindCluster2_Functions.R")
Mast_plat2_try=FindClusters2(Mast_plat2, cluster.range=c(2,3), by=0.1, res=0.5, verbose=T)
DimPlot(Mast_plat2_try, label=T, label.size=5, repel=T) + NoLegend()
table(Mast_plat2_try$seurat_clusters)
Idents(Mast_plat2_try)="seurat_clusters"
Mast_plat_try$seurat_clusters=as.character(Mast_plat_try$seurat_clusters)
Mast_plat_try[[]][WhichCells(Mast_plat2_try, idents=0),]$seurat_clusters=16
Mast_plat_try[[]][WhichCells(Mast_plat2_try, idents=1),]$seurat_clusters=16

Mast_plat3=subset(Mast_plat_try, seurat_clusters==6)
library(dplyr)
library(Seurat)
Mast_plat3=Mast_plat3 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)
source("~/Rscripts/FindCluster2_Functions.R")
Mast_plat3_try=FindClusters2(Mast_plat3, cluster.range=2, by=0.1, res=0.5, verbose=T)
DimPlot(Mast_plat3_try, label=T, label.size=5, repel=T) + NoLegend()
table(Mast_plat3_try$seurat_clusters)
Idents(Mast_plat3_try)="seurat_clusters"
Mast_plat_try$seurat_clusters=as.character(Mast_plat_try$seurat_clusters)
Mast_plat_try[[]][WhichCells(Mast_plat3_try, idents=0),]$seurat_clusters=17

# find overexpressed genes
markers=FindAllMarkers(Mast_plat_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(Mast_plat_try, features=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9"), pt.size=0)

# cluster0: Platelet mixed with CD4 Tcm
VlnPlot(Mast_plat_try, features=c("IL7R","SELL","CD44","CD27","AQP3","IL7R","MAL","LTB"), pt.size=0)
# cluster1: Platelet mixed with CD14 Mono
VlnPlot(Mast_plat_try, features=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9",
                                  "S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14"), pt.size=0)
# cluster2: Platelet mixed with resting CD4 Tnaive
VlnPlot(Mast_plat_try, features=c("MAL","CCR7","TCF7","IL7R","FHIT","LEF1","NOSIP","LDHB","TRABD2A"), pt.size=0)
# cluster3: Platelet mixed with CD8 Tem
VlnPlot(Mast_plat_try, features=c("IL7R","SELL","CD44","CD8A","CD8B","CD3D","CD3G","GZMH","KLRD1","NKG7","GZMK"), pt.size=0)
# cluster4: Platelet
VlnPlot(Mast_plat_try, features=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9"), pt.size=0)
# cluster5: Platelet mixed with activated CD4 Tnaive
VlnPlot(Mast_plat_try, features=c("MAL","CCR7","TCF7","IL7R","FHIT","LEF1","NOSIP","LDHB","TIMP1","CCL5"), pt.size=0)
# cluster6: Platelet mixed with transitional CD56dim NK (previously reported as "active" and "transitional" hNK_Bm2)
VlnPlot(Mast_plat_try, features=c("GNLY","NKG7","GZMB","SELL","KLRC1","XCL2"), pt.size=0)
# cluster7: Platelet mixed with CD4 Tnaive which have some immunoregulatory properties
VlnPlot(Mast_plat_try, features=c("MAL","CCR7","TCF7","IL7R","FHIT","LEF1","CCL5","BCL11B","TSHZ2","PIM2"), pt.size=0)
# cluster8: Platelet mixed with Bnaive
VlnPlot(Mast_plat_try, features=c("IGHM","IGHD","MS4A1","TCL1A","CD79A","CD79B"), pt.size=0)
# cluster9: B intermediate (expressing naive markers - IL4R+TCL1A+, memory markers - CD24+CD27+, unswithced - IgMhi, switched - IgG+IgA+)
VlnPlot(Mast_plat_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","IGHG3"), pt.size=0)
# cluster10: B1-like cells (...adult peripheral blood that express CD20, CD27, and CD43 as human B1 cells,...CD20+ CD5+ B-cells...human B1-like cells subset partly overlaps with that of circulating human CD27+ IgM+ IgD+ B-cells/IgM MBCs)
VlnPlot(Mast_plat_try, features=c("IGHM","IGHD","CD27","IGHD","CD5","MS4A1"), pt.size=0)
# cluster11: Platelet mixed with nonclassical CD16+ Mono
VlnPlot(Mast_plat_try, features=c("CD14","FCGR3A","CDKN1C","HES4","TCF7L2","MS4A7","IFITM3"), pt.size=0)
# cluster12: Platelet mixed with CD8 Tcm
VlnPlot(Mast_plat_try, features=c("IL7R","SELL","CD44","CD27","CD8A","CD8B","CD3D","CD3G"), pt.size=0)
# cluster13: Mast cells
VlnPlot(Mast_plat_try, features=c("TPSAB1","TPSB2"), pt.size=0)
# cluster14: Platelet mixed with cDC2
VlnPlot(Mast_plat_try, features=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO"), pt.size=0)
# cluster15: Platelet mixed with pDC
VlnPlot(Mast_plat_try, features=c("PLD4","SERPINF1","LILRA4","IL3RA","SMPD3"), pt.size=0)
# cluster16: CD14 Mono (slightly contaminated by platelets)
VlnPlot(Mast_plat_try, features=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9",
                                  "S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14"), pt.size=0)
# cluster17: Platelet mixed with cytotoxic CD56dim NK (previously reported as "mature" and "terminal" hNK_Bm1)
VlnPlot(Mast_plat_try, features=c("GNLY","NKG7","TRDC","PRF1","FGFBP2","KLRF1","SELL","KLRC1","XCL2"), pt.size=0)


highlighted=WhichCells(Mast_plat_try, idents=10)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(Platelet.mix_CD4Tcm=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9"),
               Platelet.mix_CD14Mono=c("IL7R","SELL","CD44","CD27","AQP3","IL7R","MAL","LTB"),
               Platelet.mix_restingCD4Tn=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9","S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14"),
               Platelet.mix_CD8Tem=c("IL7R","SELL","CD44","CD8A","CD8B","CD3D","CD3G","GZMH","KLRD1","NKG7","GZMK"),
               Platelet=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9"),
               Platelet.mix_activCD4Tn=c("MAL","CCR7","TCF7","IL7R","FHIT","LEF1","NOSIP","LDHB","TIMP1","CCL5"),
               Platelet.mix_transNK=c("GNLY","NKG7","GZMB","SELL","KLRC1","XCL2"),
               Platelet.mix_immunoregCD4Tn=c("MAL","CCR7","TCF7","IL7R","FHIT","LEF1","CCL5","BCL11B","TSHZ2","PIM2"),
               Platelet.mix_Bnaive=c("IGHM","IGHD","MS4A1","TCL1A","CD79A","CD79B"),
               Bint=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","IGHG3"),
               B1_like=c("IGHM","IGHD","CD27","IGHD","CD5","MS4A1"),
               Platelet.mix_CD16Mono=c("CD14","FCGR3A","CDKN1C","HES4","TCF7L2","MS4A7","IFITM3"),
               Platelet.mix_CD8Tcm=c("IL7R","SELL","CD44","CD27","CD8A","CD8B","CD3D","CD3G"),
               Mast_cell=c("TPSAB1","TPSB2"),
               Platelet.mix_cDC2=c("FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO"),
               Platelet.mix_pDC=c("PLD4","SERPINF1","LILRA4","IL3RA","SMPD3"),
               CD14_Mono.mix_platelet=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9","S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14"),
               Platelet.mix_CD56dimNK=c("GNLY","NKG7","TRDC","PRF1","FGFBP2","KLRF1","SELL","KLRC1","XCL2"))
p2=
  jjDotPlot(Mast_plat_try,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Mast cells and platelets")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=0), "Platelet.mix_CD4Tcm", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=1), "Platelet.mix_CD14Mono", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=2), "Platelet.mix_CD4Tnaive.TRABD2A+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=3), "Platelet.mix_CD8Tem", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=4), "Platelet", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=5), "Platelet.mix_CD4Tnaive.CCL5hi", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=6), "Platelet.mix_transNK", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=7), "Platelet.mix_CD4Tnaive.PIM2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=8), "Platelet.mix_Bnaive", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=9), "Bintermediate.IL4R+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=10), "B1.CD5+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=11), "Platelet.mix_CD16Mono", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=12), "Platelet.mix_CD8Tcm", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=13), "Mast", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=14), "Platelet.mix_cDC2", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=15), "Platelet.mix_pDC", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=16), "Mono.CD14.lin_infidel.PF4+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Mast_plat_try, idents=17), "Platelet.mix_CD56dimNK", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of MAIT and gdT
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
MAIT_gdT=subset(object, reassigned %in% c("MAIT","gdT"))

library(dplyr)
library(Seurat)
MAIT_gdT=MAIT_gdT %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
MAIT_gdT_try=FindClusters2(MAIT_gdT, cluster.range=c(3,5), by=0.1, res=0.2, verbose=T)
DimPlot(MAIT_gdT_try, label=T, label.size=5, repel=T) + NoLegend()
table(MAIT_gdT_try$seurat_clusters)
Idents(MAIT_gdT_try)="seurat_clusters"

# find overexpressed genes
markers=FindAllMarkers(MAIT_gdT_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(MAIT_gdT_try, features=c("SLC4A10","KLRB1","TRAV1-2"), pt.size=0)

# cluster0: effector gdT as SELL-CD27-CD44+
VlnPlot(MAIT_gdT_try, features=c("SELL","CD44","CD3D","CD3E","CD3G",
                                 "TRDC","TRGC1","TRGC2"), pt.size=0)
# cluster1: gdT as TRDC+TRGC1+TRGC2+, most likely Tcm as SELL+CD27+CD44+
VlnPlot(MAIT_gdT_try, features=c("SELL","CD27","CD44","IL7R","CD3D","CD3E","CD3G",
                                 "TRDC","TRGC1","TRGC2","KLRC1"), pt.size=0)
# cluster2: MAIT cells
VlnPlot(MAIT_gdT_try, features=c("CEBPD","KLRB1","SLC4A10","GZMK","NCR3","IL7R"), pt.size=0)
# cluster3: Platelet mixed with gdT
VlnPlot(MAIT_gdT_try, features=c("PPBP","PF4","NRGN","GNG11",
                                 "TRDC","TRGC1","TRGC2"), pt.size=0)

highlighted=WhichCells(MAIT_gdT_try, idents=2)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(gdT_eff=c("SELL","CD44","CD3D","CD3E","CD3G","TRDC","TRGC1","TRGC2"),
               gdT_Tcm=c("SELL","CD27","CD44","IL7R","CD3D","CD3E","CD3G","TRDC","TRGC1","TRGC2","KLRC1"),
               MAIT=c("CEBPD","KLRB1","SLC4A10","GZMK","NCR3","IL7R"),
               gdT_infide_toPlatelet=c("PPBP","PF4","NRGN","GNG11","TRDC","TRGC1","TRGC2"))
p2=
  jjDotPlot(MAIT_gdT_try,
            scale=F,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("MAIT and gdT")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(MAIT_gdT_try, idents=0), "gdTeff", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(MAIT_gdT_try, idents=1), "gdTcm", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(MAIT_gdT_try, idents=2), "MAIT", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(MAIT_gdT_try, idents=3), "gdT.lin_infidel.PF4+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of ISAGhi T and NKT
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
ISAsigT_NKT=subset(object, reassigned %in% c("ISAGhi T","NKT"))

library(dplyr)
library(Seurat)
ISAsigT_NKT=ISAsigT_NKT %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
ISAsigT_NKT_try=FindClusters2(ISAsigT_NKT, cluster.range=c(3,5), by=0.1, res=0.3, verbose=T)
DimPlot(ISAsigT_NKT_try, label=T, label.size=5, repel=T) + NoLegend()
table(ISAsigT_NKT_try$seurat_clusters)
Idents(ISAsigT_NKT_try)="seurat_clusters"

# find overexpressed genes
markers=FindAllMarkers(ISAsigT_NKT_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters
VlnPlot(MAIT_gdT_try, features=c("SLC4A10","KLRB1","TRAV1-2"), pt.size=0)

# cluster0: NKT
VlnPlot(ISAsigT_NKT_try, features=c("CD3D","CD3G","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G"), pt.size=0)
# cluster1: IFN-I T at the status of Tcm
VlnPlot(ISAsigT_NKT_try, features=c("LIMS1","MAL","LTB","AQP3","MX1","ISG15","ISG20"), pt.size=0)
# cluster2: IFN-I T at the status of transcriptional repression
VlnPlot(ISAsigT_NKT_try, features=c("MX1","ISG15","ISG20","TSHZ2","FHIT","RNASET2"), pt.size=0)
# cluster3: CD8 Tnaive with cell regulation on cell migration
VlnPlot(ISAsigT_NKT_try, features=c("LINC02446","CD8A","CD8B","CCR7","LEF1","RGS10","AIF1","ACTN1"), pt.size=0)
# cluster4: CD56bright NK with responses to cytokines and regulation on cytokine production
VlnPlot(ISAsigT_NKT_try, features=c("XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A",
                                    "COTL1","XCL1","LTB","GZMK"), pt.size=0)

highlighted=WhichCells(ISAsigT_NKT_try, idents=0)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(NKT=c("CD3D","CD3G","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G"),
               IFNpos_Tcm=c("LIMS1","MAL","LTB","AQP3","MX1","ISG15","ISG20"),
               IFNpos_Ttranscriprepress=c("MX1","ISG15","ISG20","TSHZ2","FHIT","RNASET2"),
               CD8Tnaive_migration=c("LINC02446","CD8A","CD8B","CCR7","LEF1","RGS10","AIF1","ACTN1"),
               gd_NKT=c("CD3D","CD3G","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G","TRDC","TRGC1","TRGC2","KLRC1"))
p2=
  jjDotPlot(ISAsigT_NKT_try,
            # scale=F,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("ISAsigT and NKT")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(ISAsigT_NKT_try, idents=0), "NKT", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(ISAsigT_NKT_try, idents=1), "CD4 Tcm.ISG15+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(ISAsigT_NKT_try, idents=2), "CD4 Tcm.ISG15+ TSHZ2+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(ISAsigT_NKT_try, idents=3), "CD8 Tnaive.AIF1+", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(ISAsigT_NKT_try, idents=4), "NK.CD56bright COTL1+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")



### Identify the subsets of Erythroid cell, dnT, ILC
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_withCD8TNK.rds")
table(object$reassigned_rough)
table(object$reassigned_interm)
table(object$reassigned)
Ery_dnT_ILC=subset(object, reassigned %in% c("dnT","Erythroid cell","ILC"))

library(dplyr)
library(Seurat)
Ery_dnT_ILC=Ery_dnT_ILC %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)

source("~/Rscripts/FindCluster2_Functions.R")
Ery_dnT_ILC_try=FindClusters2(Ery_dnT_ILC, cluster.range=3, by=0.1, res=0.2, verbose=T)
DimPlot(Ery_dnT_ILC_try, label=T, label.size=5, repel=T) + NoLegend()
table(Ery_dnT_ILC_try$seurat_clusters)
Idents(Ery_dnT_ILC_try)="seurat_clusters"

# find overexpressed genes
markers=FindAllMarkers(Ery_dnT_ILC_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

# define the clusters

# cluster0: dnT
VlnPlot(Ery_dnT_ILC_try, features=c("CD3D","CD3G","CD8A","CD4",
                                    "GZMK","NUCB2","CD8B","GPR183","TCF7","LYAR","MALAT1","C12orf57","LEF1","LDHB"), pt.size=0)
# cluster1: Erythroid
VlnPlot(Ery_dnT_ILC_try, features=c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","TRIM58","SELENBP1"), pt.size=0)
# cluster2: ILC (activation markers: TYROBP, FCER1G, most likley ILC2 as GATA3+)
VlnPlot(Ery_dnT_ILC_try, features=c("IL7R","FXYD7","TRDC","IL2RA","KLRB1","TNFRSF18","TNFRSF4","TYROBP","FCER1G","GATA3"), pt.size=0)


highlighted=WhichCells(Ery_dnT_ILC_try, idents=2)
DimPlot(object, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()

library(scRNAtoolVis)
features_=list(dnT=c("CD3D","CD3G","CD8A","CD4","GZMK","NUCB2","CD8B","GPR183","TCF7","LYAR","MALAT1","C12orf57","LEF1","LDHB"),
               Eryth=c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","TRIM58","SELENBP1"),
               ILC2=c("IL7R","FXYD7","TRDC","IL2RA","KLRB1","TNFRSF18","TNFRSF4","TYROBP","FCER1G","GATA3"))
p2=
  jjDotPlot(Ery_dnT_ILC_try,
            # scale=F,
            gene=unlist(features_),
            ytree=F,
            base_size=10,
            bar.legendTitle="Mean expression",
            point.lengdTitle="Fraction (%)") +
  ggtitle("Eryth, dnT, and ILC")

# add the annotation
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
obj_annotation=obj_annotation %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Ery_dnT_ILC_try, idents=0), "gdT", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Ery_dnT_ILC_try, idents=1), "Eryth", CD8T_NK_subsets)) %>%
  mutate(CD8T_NK_subsets=ifelse(cell %in% WhichCells(Ery_dnT_ILC_try, idents=2), "ILC.GATA3+", CD8T_NK_subsets))

write.table(obj_annotation, "~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt", sep="\t")
dev.off()

#####################################



### Draw 122-subset DimPlot
#####################################
###

###
obj_annotation=read.delim("~/Project_PBMCage/Results/all_cells_annoted_CD8T_NK_subsets.txt")
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")
table(rownames(obj_annotation)==colnames(object)) # check
object=AddMetaData(object, metadata=obj_annotation$CD8T_NK_subsets, col.name="Celltype")
saveRDS(object, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

pdf("~/Project_PBMCage/Plots/122subsets.pdf", height=12, width=35)
DimPlot(object, group.by="Celltype", raster=FALSE)
dev.off()

#####################################



### Draw violin plot of freq vs ages
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets2
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)

# make a function
library(ggpubr)
draw_proportion_perAge=function(df, celltype.drawn) {
  df.used=subset(data, CellType==celltype.drawn)
  df.used=df.used %>%
    dplyr::add_count(Age, name="Age_n") %>%
    as.data.frame()
  df.used$Age=as.factor(df.used$Age)
  p=
    ggplot(df.used, aes(x=Age, y=Freq)) +
    geom_boxplot(color="gray80", outlier.size=0.5) +
    # geom_point(df.used %>% subset(Age_n==1), aes(x=Age, y=Freq)) +
    theme_classic() +
    labs(x=NULL, y="Cell percent (%)") +
    theme(axis.text.x=element_text(size=3),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0("Percent of ", celltype.drawn)) +
    stat_compare_means(method="anova", label.x=10, label.y=max(df.used$Freq))
  
  median_1=unlist(lapply(split(df.used$Freq, df.used$Age), median))
  p1_l=data.frame(Age=names(median_1), Freq=median_1)
  p=p + geom_line(data=p1_l, aes(x=Age, y=Freq, group=1), linetype="dashed", color="black")
  plot(p)
}

# plot cell freq vs. ages for each celltype
celltypes_ordered=names(table(data$CellType))
celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B", celltypes_ordered)],
                      CD4_T=celltypes_ordered[grepl("^CD4 ", celltypes_ordered)],
                      CD8_T=celltypes_ordered[grepl("^CD8 ", celltypes_ordered)],
                      CD4_T=celltypes_ordered[grepl("^DC", celltypes_ordered)],
                      Erythroid=celltypes_ordered[grepl("^Eryth", celltypes_ordered)],
                      gdT=celltypes_ordered[grepl("^gdT", celltypes_ordered)],
                      ILC=celltypes_ordered[grepl("^ILC", celltypes_ordered)],
                      MAIT=celltypes_ordered[grepl("^MAIT", celltypes_ordered)],
                      Mast_cell=celltypes_ordered[grepl("^Mast", celltypes_ordered)],
                      Mono=celltypes_ordered[grepl("^Mono", celltypes_ordered)],
                      NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
                      NKT=celltypes_ordered[grepl("^NKT", celltypes_ordered)],
                      Platelet=celltypes_ordered[grepl("^Platelet", celltypes_ordered)],
                      HSPC=celltypes_ordered[grepl("^Progenitor", celltypes_ordered)],
                      Treg=celltypes_ordered[grepl("^Treg", celltypes_ordered)])
length_max=max(sapply(celltypes_groups, length))
plistcowplot=list()
for (i in 1:length(celltypes_groups)) {
  current_types=celltypes_groups[[i]]
  plist=list()
  for (j in 1:length(current_types)) {
    current_type=current_types[j]
    plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
  }
  length_of_plist=length(plist)
  codes=""
  for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
  for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")
  codes=gsub("^, ","",codes)
  eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", nrow=5, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
}
pdf("~/Project_PBMCage/Plots/122subsets_Proportion_vs._age.pdf", width=25, height=16)
for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
dev.off()

#####################################



### Draw violin plot of counts vs ages
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets2
cellNum=table(CellType_, ID_)
data=as.data.frame(t(cellNum))
colnames(data)=c("ID","CellType","Counts")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)

# make a function
library(ggpubr)
draw_proportion_perAge=function(df, celltype.drawn) {
  df.used=subset(data, CellType==celltype.drawn)
  df.used=df.used %>%
    dplyr::add_count(Age, name="Age_n") %>%
    as.data.frame()
  df.used$Age=as.factor(df.used$Age)
  p=
    ggplot(df.used, aes(x=Age, y=Counts)) +
    geom_boxplot(color="gray80", outlier.size=0.5) +
    theme_classic() +
    labs(x=NULL, y="Cell counts") +
    theme(axis.text.x=element_text(size=3),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0("Counts of ", celltype.drawn)) +
    stat_compare_means(method="anova", label.x=10, label.y=max(df.used$Counts))
  
  median_1=unlist(lapply(split(df.used$Counts, df.used$Age), median))
  p1_l=data.frame(Age=names(median_1), Counts=median_1)
  p=p + geom_line(data=p1_l, aes(x=Age, y=Counts, group=1), linetype="dashed", color="black")
  plot(p)
}

# plot cell counts vs. ages for each celltype
celltypes_ordered=names(table(data$CellType))
celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B", celltypes_ordered)],
                      CD4_T=celltypes_ordered[grepl("^CD4 ", celltypes_ordered)],
                      CD8_T=celltypes_ordered[grepl("^CD8 ", celltypes_ordered)],
                      CD4_T=celltypes_ordered[grepl("^DC", celltypes_ordered)],
                      Erythroid=celltypes_ordered[grepl("^Eryth", celltypes_ordered)],
                      gdT=celltypes_ordered[grepl("^gdT", celltypes_ordered)],
                      ILC=celltypes_ordered[grepl("^ILC", celltypes_ordered)],
                      MAIT=celltypes_ordered[grepl("^MAIT", celltypes_ordered)],
                      Mast_cell=celltypes_ordered[grepl("^Mast", celltypes_ordered)],
                      Mono=celltypes_ordered[grepl("^Mono", celltypes_ordered)],
                      NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
                      NKT=celltypes_ordered[grepl("^NKT", celltypes_ordered)],
                      Platelet=celltypes_ordered[grepl("^Platelet", celltypes_ordered)],
                      HSPC=celltypes_ordered[grepl("^Progenitor", celltypes_ordered)],
                      Treg=celltypes_ordered[grepl("^Treg", celltypes_ordered)])
length_max=max(sapply(celltypes_groups, length))
plistcowplot=list()
for (i in 1:length(celltypes_groups)) {
  current_types=celltypes_groups[[i]]
  plist=list()
  for (j in 1:length(current_types)) {
    current_type=current_types[j]
    plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
  }
  length_of_plist=length(plist)
  codes=""
  for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
  for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")
  codes=gsub("^, ","",codes)
  eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", nrow=5, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
}
pdf("~/Project_PBMCage/Plots/122subsets_Counts_vs._age.pdf", width=25, height=16)
for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
dev.off()

#####################################



### Analyze the lineage infidelity
#####################################
###

### Plot the scores of lineage markers in the infidel cells by heatmap
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
infidel_types=names(table(object$CD8T_NK_subsets2))
infidel_types=infidel_types[grepl("lin_infidel|mix", infidel_types)]
infidel_subsets=subset(object, CD8T_NK_subsets2 %in% infidel_types)
infidel_subsets[["Possible lineage infidelity"]]=infidel_subsets$CD8T_NK_subsets2
# Score lineage markers
Lin_markers=list(HSPC=c("PRSS57","CYTL1","EGFL7","GATA2","CD34","SMIM24","AVP","LAPTM4B"),
                 B=c("MS4A1","IGHM","IGHD","CD79A","CD79B","BANK1","TNFRSF13C"),
                 NK=c("TRGC1-","TRGV9-","TRDV2-","CD3G-","TRAC-","KLRC1","KLRC2","KLRD1","KLRF1","NCAM1","NCR1"),
                 "T"=c("GPR183-","TRGC1-","TRGV9-","TRDV2-","CD4","CD8A","CD28","CD38","CTLA4","LEF1","TCF7","TRAC"),
                 Mono=c("FCGR1A","S100A9","S100A8","LYZ","VCAN","S100A12","IL1B","FCN1","G0S2"),
                 DC=c("AXL","CCDC88A","CLEC10A","CLEC4C","CLEC9A"),
                 Platelet=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","RGS18","GP9"),
                 Eryth=c("HBD","HBM","AHSP","ALAS2","CA1","SLC4A1","IFIT1B","TRIM58","SELENBP1","TMCC2"),
                 Mast=c("CMA1","MS4A2","TPSAB1","TPSB2"))
library(UCell)
library(ggplot2)
library(patchwork)
p=SCpubr::do_EnrichmentHeatmap(infidel_subsets,
                               input_gene_list=Lin_markers,
                               group.by="Possible lineage infidelity",
                               flavor="UCell",
                               scale_scores=FALSE)
p

### Plot the freq of infidel cells in the all that type of cells
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes_r=names(table(object$CD8T_NK_subsets2_rough))
celltypes_inter=names(table(object$CD8T_NK_subsets2_inter))
celltypes_d=names(table(object$CD8T_NK_subsets2))
infidel_types=celltypes_d[grepl("lin_infidel|mix", celltypes_d)]
# lin_infidel counts
B_cell_counts=ncol(subset(object, CD8T_NK_subsets2_inter %in% c("B.intermediate","B.mem","B.naive","B.other","B.pbpc")))
T_cell_counts=ncol(subset(object, CD8T_NK_subsets2_inter %in% c("CD4T.naive","CD4T.prolif","CD4T.reg","CD4T.Tcm","CD4T.Tem","CD8T.naive","CD8T.Tcm","CD8T.Tem")))
DC_counts=ncol(subset(object, CD8T_NK_subsets2_inter %in% c("DC.allsubtypes")))
Mono_counts=ncol(subset(object, CD8T_NK_subsets2_inter %in% c("Mono.classical","Mono.intermediate","Mono.nonclassical")))
NK_counts=ncol(subset(object, CD8T_NK_subsets2_inter %in% c("NK.immature","NK.mature","NK.prolif")))
NKT_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("NKT","NKT.MKI67+")))
gdT_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("gdT","gdTcm","gdTeff")))
Plt_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Platelet")))
B.infidel.T_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("B.lin_infidel.TRBC1+","CD8 Tem.lin_infidel.MS4A1+")))
B.infidel.NK_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("NK.lin_infidel.MS4A1+")))
B.infidel.NKT_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("B.lin_infidel.TRBC1+ KLRD1+","Bpc.lin_infidel.TRBC1+ KLRD1+")))
B.infidel.DC_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("DC.CD1C.lin_infidel.MS4A1+")))
B.infidel.Mono_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Mono.CD14.lin_infidel.MS4A1+","Mono.CD16.lin_infidel.MS4A1+")))
B.infidel.Plt_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Platelet.mix_Bnaive")))
DC.infidel.NKT_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("DC.CD1C.lin_infidel.TRBC1+ KLRD1+")))
DC.infidel.T_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("DC.IL3RA.lin_infidel.TRBC1+")))
DC.infidel.Plt_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Platelet.mix_cDC2","Platelet.mix_pDC")))
Mono.infidel.T_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Mono.CD14.lin_infidel.TRBC1+","Mono.CD16.lin_infidel.TRBC1+")))
Plt.infidel.T_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Platelet.mix_CD4Tcm","Platelet.mix_CD4Tnaive.CCL5hi","Platelet.mix_CD4Tnaive.PIM2+","Platelet.mix_CD4Tnaive.TRABD2A+","Platelet.mix_CD8Tcm","Platelet.mix_CD8Tem")))
Plt.infidel.gdT_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("gdT.lin_infidel.PF4+")))
Plt.infidel.Mono_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Mono.CD14.lin_infidel.PF4+","Platelet.mix_CD14Mono","Platelet.mix_CD16Mono")))
Plt.infidel.NK_counts=ncol(subset(object, CD8T_NK_subsets2 %in% c("Platelet.mix_CD56dimNK","Platelet.mix_transNK")))
# make a matrix with the counts
count_matrix=data.frame("B"=c(B_cell_counts,B.infidel.T_counts,B.infidel.DC_counts,B.infidel.Mono_counts,B.infidel.NK_counts,B.infidel.NKT_counts,0,B.infidel.Plt_counts),
                        "T"=c(B.infidel.T_counts,T_cell_counts,DC.infidel.T_counts,Mono.infidel.T_counts,0,0,0,Plt.infidel.T_counts),
                        "DC"=c(B.infidel.DC_counts,DC.infidel.T_counts,DC_counts,0,0,DC.infidel.NKT_counts,0,DC.infidel.Plt_counts),
                        "Mono"=c(B.infidel.Mono_counts,Mono.infidel.T_counts,0,Mono_counts,0,0,0,Plt.infidel.Mono_counts),
                        "NK"=c(B.infidel.NK_counts,0,0,0,NK_counts,0,0,Plt.infidel.NK_counts),
                        "NKT"=c(B.infidel.NKT_counts,0,DC.infidel.NKT_counts,0,0,NKT_counts,0,0),
                        "gdT"=c(0,0,0,0,0,0,gdT_counts,Plt.infidel.gdT_counts),
                        "Plt"=c(B.infidel.Plt_counts,Plt.infidel.T_counts,DC.infidel.Plt_counts,Plt.infidel.Mono_counts,Plt.infidel.NK_counts,0,Plt.infidel.gdT_counts,Plt_counts))
rownames(count_matrix)=c("B","T","DC","Mono","NK","NKT","gdT","Plt")
write.table(count_matrix, "~/Project_PBMCage/Results/LinInfidel_counts.txt", sep="\t")
count_matrix_freq=t(t(count_matrix)/colSums(count_matrix))
write.table(count_matrix_freq, "~/Project_PBMCage/Results/LinInfidel_freq.txt", sep="\t")
# plot the result with heatmap
library(ComplexHeatmap)
count_matrix_freq=read.delim("~/Project_PBMCage/Results/LinInfidel_freq.txt")
col_fun=circlize::colorRamp2(c(0, 0.01, 0.5), c("white", "#B66DFF", "#490092"))
Heatmap(count_matrix_freq,
        name="Freq",
        show_heatmap_legend=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        row_names_gp=gpar(fontsize=8),
        column_names_gp=gpar(fontsize=8),
        column_title="",
        col=col_fun,
        height=nrow(count_matrix_freq)*unit(5,"mm"),
        width=ncol(count_matrix_freq)*unit(5,"mm"))
count_matrix=read.delim("~/Project_PBMCage/Results/LinInfidel_counts.txt")
col_fun=circlize::colorRamp2(c(0, 1000, 8000), c("white", "#B66DFF", "#490092"))
Heatmap(count_matrix,
        name="Counts",
        show_heatmap_legend=TRUE,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        row_names_gp=gpar(fontsize=8),
        column_names_gp=gpar(fontsize=8),
        column_title="",
        col=col_fun,
        height=nrow(count_matrix_freq)*unit(5,"mm"),
        width=ncol(count_matrix_freq)*unit(5,"mm"))

### Plot the expression of top genes in each lin_infidel celltypes
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes_d=names(table(object$CD8T_NK_subsets2))
infidel_types=celltypes_d[grepl("lin_infidel|mix", celltypes_d)]
infidel_subsets=subset(object, CD8T_NK_subsets2 %in% infidel_types)
celltypes_infidel=names(table(infidel_subsets$CD8T_NK_subsets2))
# create the top gene list
data=edgeR::Seurat2PB(infidel_subsets, sample="orig.ident", cluster="CD8T_NK_subsets2")
counts=data$counts
counts=counts[!(grepl("^RP[SL]",rownames(counts)) | grepl("^RP[0-9]|^LINC|^MT-",rownames(counts))),]
colnames(counts)=gsub(".*_cluster","",colnames(counts))
CellType_=infidel_subsets$CD8T_NK_subsets2
cellNum=as.data.frame(table(CellType_))
table(colnames(counts)==cellNum$CellType_) # check
avg_expr=t(t(counts)/cellNum$Freq)
avg_expr[avg_expr==0]=NA
avg_expr=as.data.frame(log2(t(t(avg_expr)/colMeans(avg_expr, na.rm=TRUE))))
avg_expr$gene=rownames(avg_expr)
colnames(avg_expr)=c("B_T","B_NKT","Bpc_NKT","CD8Tem_B",
                     "cDC1_B","cDC_NKT","pDC_T","gdT_plt",
                     "CD14Mono_B","CD14Mono_plt","CD14Mono_T","CD16Mono_B",
                     "CD16Mono_T","NK_B","Plt_B","Plt_CD14Mono",
                     "Plt_CD16Mono","Plt_CD4Tcm","Plt_CD4Tn.1","Plt_CD4Tn.2","Plt_CD4Tn.3","Plt_iNK","Plt_CD8Tcm","Plt_CD8Tem",
                     "Plt_cDC2","Plt_pDC","Plt_trNK","gene")
avg_expr_long=avg_expr %>%
  tidyr::pivot_longer(cols=colnames(avg_expr)[1:(ncol(avg_expr)-1)], names_to="cluster", values_to="avg_log2FC")
avg_expr_long=na.omit(avg_expr_long)
avg_expr_long=avg_expr_long %>% group_by(cluster) %>% slice_max(avg_log2FC, n=2500) %>% as.data.frame(.)
log2FC.cutoff=0.25
avg_expr_long=avg_expr_long %>%
  mutate(p_val=ifelse(abs(avg_log2FC)<log2FC.cutoff, 1, 
                      ifelse(T, 1e-3, 0.5)),
         p_val_adj=ifelse(abs(avg_log2FC)<log2FC.cutoff, 1, 
                          ifelse(T, 1e-4, 0.1)))
# Plot the vulcano plot and label the marker genes
mygene=list("B"=c("CD79A","RALGPS2","CD79B","MS4A1","BANK1","CD74","TNFRSF13C","HLA-DQA1","IGHM","MEF2C"),
            "T"=c("IL7R","MAL","LTB","CD4","LDHB","TPT1","TRAC","TMSB10","CD3D","CD3G","CD8B","CD8A","CD3D","TMSB10","HCST","CD3G","LINC02446","CTSW","CD3E","TRAC",
                  "CCL5","GZMH","CD8A","TRAC","KLRD1","NKG7","GZMK","CST7","CD8B","TRGC2",
                  "CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB",
                  "IL7R","TMSB10","CD4","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL",
                  "TCF7","CD4","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1"),
            "NK"=c("NKG7","KLRD1","TYROBP","GNLY","FCER1G","PRF1","CD247","KLRF1","CST7","GZMB",
                   "GNLY","TYROBP","NKG7","FCER1G","GZMB","TRDC","PRF1","FGFBP2","SPON2","KLRF1",
                   "XCL2","FCER1G","SPINK2","TRDC","KLRC1","XCL1","SPTSSB","PPP1R9A","NCAM1","TNFRSF11A"),
            "Bpc"=c("MZB1","JCHAIN","TNFRSF17","ITM2C","DERL3","TXNDC5","POU2AF1","IGHA1","TXNDC11","CD79A"),
            "DC"=c("CLEC9A","DNASE1L3","C1orf54","IDO1","CLNK","CADM1","FLT3","ENPP1","XCR1","NDRG2",
                   "FCER1A","HLA-DQA1","CLEC10A","CD1C","ENHO","PLD4","GSN","SLC38A1","NDRG2","AFF3",
                   "ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3"),
            "gdT"=c("TRDC","TRGC1","TRGC2","KLRC1","NKG7","TRDV2","CD7","TRGV9","KLRD1","KLRG1"),
            "plt"=c("PPBP","PF4","NRGN","GNG11","CAVIN2","TUBB1","CLU","HIST1H2AC","RGS18","GP9"),
            "Mono"=c("S100A9","CTSS","S100A8","LYZ","VCAN","S100A12","IL1B","CD14","G0S2","FCN1",
                     "CDKN1C","FCGR3A","PTPRC","LST1","IER5","MS4A7","RHOC","IFITM3","AIF1","HES4")
            )
mygene=unlist(mygene)
mygene=mygene[!duplicated(mygene)]
library(scRNAtoolVis)
pdf("~/Project_PBMCage/Plots/LinInfidel_MarkerExpr.pdf", width=50, height=10)
jjVolcano(diffData=avg_expr_long,
          myMarkers=mygene,
          log2FC.cutoff=log2FC.cutoff,
          aesCol=c('#E5B17E','#7EC3E5'), # Sig.down and Sig.up color
          back.col='grey99',
          tile.col=paletteer::paletteer_c("grDevices::Cyan-Magenta", 30),
          pSize=0.25, # point size
          size=2.5, # label font size
          base_size=10, # x or y title/axis font size
          fontface='italic',
          legend.position="none") +
  labs(y="log2[Expr/Mean(Expr)]", x=NULL)
dev.off()

#####################################



### Determine individual outlier(s) from each age
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

library(Seurat)
library(edgeR)
library(DESeq2)

# PCA of donor_id grouped by age
pb_obj=Seurat2PB(object, sample="donor_id", cluster="age")
colData=list(donor_id=as.factor(pb_obj$samples$sample), age=as.factor(pb_obj$samples$cluster))

se_obj=SummarizedExperiment(assays=list(counts=pb_obj$counts), colData=colData)
dds=DESeqDataSet(se_obj, design=~donor_id)
dds=estimateSizeFactors(dds)
vsd=varianceStabilizingTransformation(dds, blind=TRUE)
sampleDists=dist(t(assay(vsd)))
plotPCA(vsd, intgroup=c("age"))

# PCA of celltypes grouped by donor_id
library(ggplot2)
library("pheatmap")
library("RColorBrewer") 
col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255) 

Age_=levels(object$age)
plot_pca=plot_heatmap=list()
processbar=txtProgressBar(min=0, max=length(Age_), style=3, width=60, char="=")
for (i in 1:length(Age_)) {
  age_=Age_[i]
  age_subset=subset(object, age==age_)
  pb_obj=Seurat2PB(age_subset, sample="donor_id", cluster="CD8T_NK_subsets2")
  colData=list(donor_id=as.factor(pb_obj$samples$sample), celltype=as.factor(gsub("\\.|\\+|\\-| ","_",pb_obj$samples$cluster)))
  se_obj=SummarizedExperiment(assays=list(counts=pb_obj$counts), colData=colData)
  dds=DESeqDataSet(se_obj, design=~donor_id+celltype)
  vsd=varianceStabilizingTransformation(dds, blind=TRUE)
  sampleDists=dist(t(assay(vsd)))
  PCA_obj=plotPCA(vsd, intgroup=c("donor_id"))
  PCA_result=PCA_obj$data
  labels=PCA_obj$labels
  percentVar_y=labels$y
  percentVar_x=labels$x
  plot_pca[[i]]=
    ggplot(PCA_result, aes(x=PC1, y=PC2, color=donor_id))+
    geom_point(size=3) +
    xlab(percentVar_x) +
    ylab(percentVar_y) +
    coord_fixed() +
    ggtitle(paste0("Age ", age_)) +
    theme_light() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  
  pb_obj_celltypeCondensed=Seurat2PB(age_subset, sample="donor_id", cluster="orig.ident")
  colData_celltypeCondensed=list(donor_id=as.factor(pb_obj_celltypeCondensed$samples$sample))
  counts=pb_obj_celltypeCondensed$counts
  colnames(counts)=gsub("_clusterSeuratProject","",colnames(counts))
  se_obj_celltypeCondensed=SummarizedExperiment(assays=list(counts=counts), colData=colData_celltypeCondensed)
  dds_celltypeCondensed=DESeqDataSet(se_obj_celltypeCondensed, design=~donor_id)
  vsd_celltypeCondensed=varianceStabilizingTransformation(dds_celltypeCondensed, blind=TRUE)
  sampleDists_celltypeCondensed=dist(t(assay(vsd_celltypeCondensed)))
  # plotPCA(vsd_celltypeCondensed, intgroup=c("donor_id"))
  sampleDistMatrix=as.matrix(sampleDists_celltypeCondensed)
  colnames(sampleDistMatrix)=NULL
  plot_heatmap_=
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists_celltypeCondensed,
             clustering_distance_cols=sampleDists_celltypeCondensed,
             col=col,
             main=paste0("Age ", age_),
             cellwidth=14,
             cellheight=14,
             treeheight_row=2,
             treeheight_col=2)
  plot_heatmap[[i]]=plot_heatmap_[[4]]
  setTxtProgressBar(processbar, i)
}
close(processbar)

pdf("~/Project_PBMCage/Plots/Sample_Distance_Heatmap.pdf", height=100, width=100)
cowplot::plot_grid(plotlist=plot_heatmap, nrow=9, align="hv")
dev.off()

pdf("~/Project_PBMCage/Plots/Sample_Distance_PCA.pdf", height=100, width=100)
cowplot::plot_grid(plotlist=plot_pca, nrow=9, align="hv")
dev.off()
# age with only 1 sample will return an error, check:
for (age_ in 94:97) {
  check_subset=subset(object[[]], age==age_)
  print(table(check_subset$donor_id)[table(check_subset$donor_id)!=0])
}
save(plot_heatmap, plot_pca, file="~/Project_PBMCage/Tempt_RDS/sample_filter_pca.RData")

# PCA of donor_id per age, all cells pooled
Age_=levels(object$age)
Age_=Age_[1:74]
plot_pca_allcelltypetogether=list()
processbar=txtProgressBar(min=0, max=length(Age_), style=3, width=60, char="=")
for (i in 1:length(Age_)) {
  age_=Age_[i]
  age_subset=subset(object, age==age_)
  pb_obj_celltypeCondensed=Seurat2PB(age_subset, sample="donor_id", cluster="orig.ident")
  colData_celltypeCondensed=list(donor_id=as.factor(pb_obj_celltypeCondensed$samples$sample))
  counts=pb_obj_celltypeCondensed$counts
  colnames(counts)=gsub("_clusterSeuratProject","",colnames(counts))
  se_obj_celltypeCondensed=SummarizedExperiment(assays=list(counts=counts), colData=colData_celltypeCondensed)
  dds_celltypeCondensed=DESeqDataSet(se_obj_celltypeCondensed, design=~donor_id)
  vsd_celltypeCondensed=varianceStabilizingTransformation(dds_celltypeCondensed, blind=TRUE)
  PCA_obj=plotPCA(vsd_celltypeCondensed, intgroup=c("donor_id"))
  PCA_result=PCA_obj$data
  labels=PCA_obj$labels
  percentVar_y=labels$y
  percentVar_x=labels$x
  plot_pca_allcelltypetogether[[i]]=
    ggplot(PCA_result, aes(x=PC1, y=PC2, color=donor_id))+
    geom_point(size=3) +
    xlab(percentVar_x) +
    ylab(percentVar_y) +
    coord_fixed() +
    ggtitle(paste0("Age ", age_)) +
    theme_light() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  setTxtProgressBar(processbar, i)
}
close(processbar)

pdf("~/Project_PBMCage/Plots/Sample_Distance_PCA_allcelltypetogether.pdf", height=100, width=100)
cowplot::plot_grid(plotlist=plot_pca_allcelltypetogether, nrow=9, align="hv")
dev.off()

#####################################



### Determine individual outlier(s) from each age by WGCNA dendrogram - at the detailed level of celltype
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes=names(table(object$CD8T_NK_subsets2))
sample_keep=list()

### Filter per celltype
processbar=txtProgressBar(min=0, max=length(celltypes), style=3, width=60, char="=")
for (i in 1:length(celltypes)) {
  # Prepare the expression matrix
  celltype_subset=subset(object, CD8T_NK_subsets2==celltypes[i])
  pb_obj=edgeR::Seurat2PB(celltype_subset, sample="age", cluster="donor_id")
  mat=as.data.frame(pb_obj$counts)
  data=log2(mat+1)
  datExpr0=as.data.frame(t(data))
  # Cluster the samples to check if there is any discrete sample
  ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
  AllkeepSamples=c()
  for (j in 1:length(ages)) {
    sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
    datExpr_=datExpr0[sample_of_this_age,]
    if (nrow(datExpr_)>=3) {
      sampleTree=hclust(dist(datExpr_), method="average")
      heights=sampleTree$height
      upper_bound=median(heights)+3*mad(heights, constant=1)
      # Remove discrete samples
      clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
      table(clust)
      keepSamples=clust==1
      print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
    } else {
      keepSamples=rep(TRUE, nrow(datExpr_))
      print("No filter.")
    }
    AllkeepSamples=c(AllkeepSamples, keepSamples)
  }
  # Filter out discrete samples
  datExpr0=datExpr0[AllkeepSamples,]
  sample_keep[[i]]=rownames(datExpr0)
  
  setTxtProgressBar(processbar, i)
}
close(processbar)

### Filter by taking all celltypes as a whole
pb_obj=edgeR::Seurat2PB(object, sample="age", cluster="donor_id")
mat=as.data.frame(pb_obj$counts)
data=log2(mat+1)
datExpr0=as.data.frame(t(data))
ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
AllkeepSamples=c()
for (j in 1:length(ages)) {
  sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
  datExpr_=datExpr0[sample_of_this_age,]
  if (nrow(datExpr_)>=3) {
    sampleTree=hclust(dist(datExpr_), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
  } else {
    keepSamples=rep(TRUE, nrow(datExpr_))
  }
  AllkeepSamples=c(AllkeepSamples, keepSamples)
}
datExpr0=datExpr0[AllkeepSamples,]
sample_keep[[length(celltypes)+1]]=rownames(datExpr0)

### Output the sample names that will keep for further analysis
names(sample_keep)=c(celltypes, "All")
saveRDS(sample_keep, "~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA.rds")

#####################################



### Combine the 122 subtypes into intermediate levels
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes_meta=object[[]]
names(table(celltypes_meta$CD8T_NK_subsets2))
celltypes_meta$CD8T_NK_subsets2_inter=celltypes_meta$CD8T_NK_subsets2
# check B
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^B.*.lin_",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bintermediate",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bmem",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bnaive",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bpb$|^Bpc$",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Btrans|B1\\.",celltypes_meta$CD8T_NK_subsets2)]))

celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^B.*.lin_",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.lin_infidel",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bintermediate",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.intermediate",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bmem",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.mem",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bnaive",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.naive",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Bpb$|^Bpc$",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.pbpc",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Btrans|B1\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "B.other",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check CD4 T
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 T\\.MKI67",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tcm\\.",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tem\\.",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tnaive\\.",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Treg\\.",celltypes_meta$CD8T_NK_subsets2)]))

celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 T\\.MKI67",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD4T.prolif",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tcm\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD4T.Tcm",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tem\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD4T.Tem",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD4 Tnaive\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD4T.naive",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Treg\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD4T.reg",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check CD8 T
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 CTL|CD8 Tem\\.",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 Tcm\\.",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 Tnaive\\.",celltypes_meta$CD8T_NK_subsets2)]))

celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 CTL|CD8 Tem\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD8T.Tem",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 Tcm\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD8T.Tcm",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^CD8 Tnaive\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "CD8T.naive",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check DCs
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^DC.*.lin_",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^DC\\.",celltypes_meta$CD8T_NK_subsets2)]))
celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^DC.*.lin_",celltypes_meta$CD8T_NK_subsets2)])),
                                       "DC.lin_infidel",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^DC\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "DC.allsubtypes",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check Monocytes
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono.*.lin_",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\- CD16\\+",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\+ CD16\\-",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\+ CD16\\+",celltypes_meta$CD8T_NK_subsets2)]))
celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono.*.lin_",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Mono.lin_infidel",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\- CD16\\+",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Mono.nonclassical",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\+ CD16\\-",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Mono.classical",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mono\\.CD14\\+ CD16\\+",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Mono.intermediate",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check other cells
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Platelet\\.mix_",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Eryth$|^Platelet$",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mast$|^ILC\\.|^NK\\.lin",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Progenitor\\.",celltypes_meta$CD8T_NK_subsets2)]))
celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Platelet\\.mix_",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Other.pltcontamin",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Eryth$|^Platelet$",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Other.RBC",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Mast$|^ILC\\.|^NK\\.lin",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Other.innateimmune",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^Progenitor\\.",celltypes_meta$CD8T_NK_subsets2)])),
                                       "Other.progenitor",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# check other T cells
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^gdT|^MAIT|^NKT",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.CD56bright",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.CD56dim",celltypes_meta$CD8T_NK_subsets2)]))
names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.MKI67",celltypes_meta$CD8T_NK_subsets2)]))

celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^gdT|^MAIT|^NKT",celltypes_meta$CD8T_NK_subsets2)])),
                                       "OtherT",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.CD56bright",celltypes_meta$CD8T_NK_subsets2)])),
                                       "NK.immature",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.CD56dim",celltypes_meta$CD8T_NK_subsets2)])),
                                       "NK.mature",
                                       CD8T_NK_subsets2_inter)) %>%
  mutate(CD8T_NK_subsets2_inter=ifelse(CD8T_NK_subsets2_inter %in% names(table(celltypes_meta$CD8T_NK_subsets2[grepl("^NK\\.MKI67",celltypes_meta$CD8T_NK_subsets2)])),
                                       "NK.prolif",
                                       CD8T_NK_subsets2_inter))
table(celltypes_meta$CD8T_NK_subsets2_inter)

# save the intermediate level of annotation
write.table(celltypes_meta, "~/Project_PBMCage/Results/all_cells_annoted_metadata_w_multiplelevels.txt", sep="\t")

###
celltypes_meta=read.delim("~/Project_PBMCage/Results/all_cells_annoted_metadata_w_multiplelevels.txt")
celltypes_meta$CD8T_NK_subsets2_rough=celltypes_meta$CD8T_NK_subsets2_inter
celltypes_meta=celltypes_meta %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("lin_infidel|contamin",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "Unassigned",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^B\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "B cells",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^CD4T\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "CD4T cells",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^CD8T\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "CD8T cells",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^DC\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "DCs",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^Mono\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "Monocytes",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^NK\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "NK cells",
                                       CD8T_NK_subsets2_rough)) %>%
  mutate(CD8T_NK_subsets2_rough=ifelse(CD8T_NK_subsets2_rough %in% names(table(celltypes_meta$CD8T_NK_subsets2_inter[grepl("^Other\\.",celltypes_meta$CD8T_NK_subsets2_inter)])),
                                       "Other cells",
                                       CD8T_NK_subsets2_rough))
table(celltypes_meta$CD8T_NK_subsets2_rough)

# save the rough level of annotation
write.table(celltypes_meta, "~/Project_PBMCage/Results/all_cells_annoted_metadata_w_multiplelevels.txt", sep="\t")

#####################################



### Transcriptional factors
#####################################
###

library(decoupleR)
library(tidyr)
net=get_collectri(organism="human", split_complexes=FALSE)

### Load
library(Seurat)
THEObj=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
assigned_types=unique(THEObj$CD8T_NK_subsets2)[!grepl("lin_infidel|mix", unique(THEObj$CD8T_NK_subsets2))]
THEObj=subset(THEObj, CD8T_NK_subsets2 %in% assigned_types)
THEObj=NormalizeData(THEObj)
# add the agecut metadata
AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}
agecut_df=THEObj[[]]
agecut_df$age=as.factor(agecut_df$age)
agecut_df$agecut=agecut_df$age
for (i in 1:length(AGECUT)) {
  agecut_df=agecut_df %>%
    mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
}
table(rownames(agecut_df)==colnames(THEObj)) # check
THEObj=AddMetaData(THEObj, metadata=agecut_df$agecut, col.name="agecut")

# Take subsets
CD8T=subset(THEObj, CD8T_NK_subsets2_rough=="CD8T cells")
CD8T_celltypes=unique(CD8T$CD8T_NK_subsets2)
CD8T_celltypes=sort(CD8T_celltypes)

DF_list=list()
for (c in 1:length(CD8T_celltypes)) {
  celltype_=CD8T_celltypes[c]
  
  message(paste0("----------", celltype_, "----------"))
  
  obj_celltype=subset(CD8T, CD8T_NK_subsets2==celltype_)
  ages=sort(unique(obj_celltype$age))
  
  for (i in 1:length(ages)) {
    message(paste0(celltype_, " - ", ages[i], " Start!"))
    
    obj_celltype_age=subset(obj_celltype, age==ages[i])
    
    mat=as.matrix(obj_celltype_age@assays$RNA@data)
    acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
    
    score=acts %>% 
      pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
      textshape::column_to_rownames("source")
    
    pvalue=acts %>% 
      pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
      textshape::column_to_rownames("source")
    
    df_score=t(score) %>%
      as.data.frame() %>%
      mutate(cluster=obj_celltype_age$CD8T_NK_subsets2) %>%
      pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
      group_by(cluster, source) %>%
      summarize(mean=mean(score, na.rm=T))
    
    df_pvalue=t(pvalue) %>%
      as.data.frame() %>%
      mutate(cluster=obj_celltype_age$CD8T_NK_subsets2) %>%
      pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
      group_by(cluster, source) %>%
      summarize(mean=mean(pvalue, na.rm=T))
    
    df_score_pvalue=df_score
    df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
    colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
    df_score_pvalue$age=ages[i]
    
    message(paste0(celltype_, " - ", ages[i], " Done!"))
  }
  
  DF_list=c(DF_list, list(df_score_pvalue))
  
  message(paste0("==========", celltype_, "=========="))
}

DF_CD8T=data.table::rbindlist(DF_list)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")














# create sub-objs for each agecut
obj_slim_list=list()
for (i in 1:length(AGECUT_name)) {
  message(paste0("AgeCut No.", i, " is loading..."))
  subset_=subset(THEObj, agecut==AGECUT_name[i])
  subset_=subset_ %>%
    NormalizeData() %>%
    FindVariableFeatures()
  subset_=SketchData(
    subset_,
    ncells=5000L,
    sketched.assay="sketch",
    method="LeverageScore",
    var.name="leverage.score",
    seed=2024L,
    cast="dgCMatrix",
    verbose=TRUE)
  assay.v5=GetAssay(subset_, assay="sketch")
  obj_slim=CreateSeuratObject(assay.v5)
  obj_slim[[]]=THEObj[[]][match(colnames(obj_slim), colnames(THEObj)),]
  obj_slim=obj_slim %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    FindNeighbors(dims=1:30) %>%
    RunUMAP(dims=1:30)
  obj_slim_list[[i]]=obj_slim
  message(paste0("AgeCut No.", i, " is done with the sketch building."))
}
names(obj_slim_list)=AGECUT_name
saveRDS(obj_slim_list, "~/Project_PBMCage/Tempt_RDS/Sketch_SeuratObj.rds")

### TF analysis
library(decoupleR)
library(tidyr)
net=get_collectri(organism="human", split_complexes=FALSE)
for (i in 1:length(obj_slim_list)) {
  message(paste0("TF analysis on AgeCut No.", i, "..."))
  obj_slim=obj_slim_list[[i]]
  mat=as.matrix(obj_slim@assays$RNA@data)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  message(paste0("TF analysis on AgeCut No.", i, " is done."))
  
  obj_slim[["tfsulm"]]=acts %>% 
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source") %>%
    CreateAssayObject(.)
  DefaultAssay(obj_slim)="tfsulm"
  obj_slim=ScaleData(obj_slim)
  obj_slim@assays$tfsulm@data=obj_slim@assays$tfsulm@scale.data
  
  obj_slim_list[[i]]=obj_slim
  message(paste0("AgeCut No.", i, " is done with the tfsulm results saving."))
}
names(obj_slim_list)=AGECUT_name
saveRDS(obj_slim_list, "~/Project_PBMCage/Tempt_RDS/Sketch_SeuratObj.rds")

### Get the top_n tfs in each agecut
library(dplyr)
n_tfs=25
TFs_list=list()
for (i in 1:length(obj_slim_list)) {
  obj_slim=obj_slim_list[[i]]
  
  # extract activities from the obj
  df=t(as.matrix(obj_slim@assays$tfsulm@data)) %>%
    as.data.frame() %>%
    mutate(cluster=obj_slim$CD8T_NK_subsets2_rough) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score))
  # get the top_n tfs with the most variable means across clusters
  tfs=df %>%
    group_by(source) %>%
    summarize(std=sd(mean)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  
  TFs_list[[i]]=tfs
}
names(TFs_list)=names(obj_slim_list)

### Get the activities based on the top_n tfs across all the agecuts
TFs=unlist(TFs_list)
TFs=TFs[!duplicated(TFs)]
for (i in 1:length(obj_slim_list)) {
  obj_slim=obj_slim_list[[i]]
  
  # extract activities from the obj
  df=t(as.matrix(obj_slim@assays$tfsulm@data)) %>%
    as.data.frame() %>%
    mutate(cluster=obj_slim$CD8T_NK_subsets2_rough) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score))
  
  # subset the extracted dataframe to the TFs
  top_acts_mat=df %>%
    filter(source %in% TFs) %>%
    pivot_wider(id_cols="cluster", names_from="source", values_from="mean") %>%
    textshape::column_to_rownames("cluster") %>%
    as.matrix()
  
  top_acts[[i]]=top_acts_mat
}
names(top_acts)=names(obj_slim_list)
saveRDS(top_acts, "~/Project_PBMCage/Tempt_RDS/Sketch_SeuratObj_TFs_dataframe.rds")

top_acts=readRDS("~/Project_PBMCage/Tempt_RDS/Sketch_SeuratObj_TFs_dataframe.rds")

library(ComplexHeatmap)
# plot
Heatmap(top_acts_mat)



#####################################

### Determine individual outlier(s) from each age by WGCNA dendrogram - at the intermediate level of celltype
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes_meta=read.delim("~/Project_PBMCage/Results/all_cells_annoted_metadata_w_multiplelevels.txt")
table(rownames(celltypes_meta)==colnames(object))
object=AddMetaData(object, metadata=celltypes_meta$CD8T_NK_subsets2_inter, col.name="CD8T_NK_subsets2_inter")
object=AddMetaData(object, metadata=celltypes_meta$CD8T_NK_subsets2_rough, col.name="CD8T_NK_subsets2_rough")
saveRDS(object, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

celltypes=names(table(object$CD8T_NK_subsets2_inter))
sample_keep=list()
### Filter per celltype at intermediate level
processbar=txtProgressBar(min=0, max=length(celltypes), style=3, width=60, char="=")
for (i in 1:length(celltypes)) {
  # Prepare the expression matrix
  celltype_subset=subset(object, CD8T_NK_subsets2_inter==celltypes[i])
  pb_obj=edgeR::Seurat2PB(celltype_subset, sample="age", cluster="donor_id")
  mat=as.data.frame(pb_obj$counts)
  data=log2(mat+1)
  datExpr0=as.data.frame(t(data))
  # Cluster the samples to check if there is any discrete sample
  ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
  AllkeepSamples=c()
  for (j in 1:length(ages)) {
    sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
    datExpr_=datExpr0[sample_of_this_age,]
    if (nrow(datExpr_)>=3) {
      sampleTree=hclust(dist(datExpr_), method="average")
      heights=sampleTree$height
      upper_bound=median(heights)+3*mad(heights, constant=1)
      # Remove discrete samples
      clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
      table(clust)
      keepSamples=clust==1
      print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
    } else {
      keepSamples=rep(TRUE, nrow(datExpr_))
      print("No filter.")
    }
    AllkeepSamples=c(AllkeepSamples, keepSamples)
  }
  # Filter out discrete samples
  datExpr0=datExpr0[AllkeepSamples,]
  sample_keep[[i]]=rownames(datExpr0)
  
  setTxtProgressBar(processbar, i)
}
close(processbar)

### Filter by taking all celltypes as a whole but remove platelet-contaminated=Other.pltcontamin cells
pb_obj=edgeR::Seurat2PB(subset(object, CD8T_NK_subsets2_inter!="Other.pltcontamin"), sample="age", cluster="donor_id")
mat=as.data.frame(pb_obj$counts)
data=log2(mat+1)
datExpr0=as.data.frame(t(data))
ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
AllkeepSamples=c()
for (j in 1:length(ages)) {
  sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
  datExpr_=datExpr0[sample_of_this_age,]
  if (nrow(datExpr_)>=3) {
    sampleTree=hclust(dist(datExpr_), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
  } else {
    keepSamples=rep(TRUE, nrow(datExpr_))
  }
  AllkeepSamples=c(AllkeepSamples, keepSamples)
}
datExpr0=datExpr0[AllkeepSamples,]
sample_keep[[length(celltypes)+1]]=rownames(datExpr0)

### Output the sample names that will keep for further analysis
names(sample_keep)=c(celltypes, "All_except_pltcontamin")
saveRDS(sample_keep, "~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA_intermediate.rds")

#####################################



### Determine individual outlier(s) from each age by WGCNA dendrogram - at the rough level of celltype
#####################################
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
celltypes=names(table(object$CD8T_NK_subsets2_rough))
sample_keep=list()
### Filter per celltype at intermediate level
processbar=txtProgressBar(min=0, max=length(celltypes), style=3, width=60, char="=")
for (i in 1:length(celltypes)) {
  # Prepare the expression matrix
  celltype_subset=subset(object, CD8T_NK_subsets2_rough==celltypes[i])
  pb_obj=edgeR::Seurat2PB(celltype_subset, sample="age", cluster="donor_id")
  mat=as.data.frame(pb_obj$counts)
  data=log2(mat+1)
  datExpr0=as.data.frame(t(data))
  # Cluster the samples to check if there is any discrete sample
  ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
  AllkeepSamples=c()
  for (j in 1:length(ages)) {
    sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
    datExpr_=datExpr0[sample_of_this_age,]
    if (nrow(datExpr_)>=3) {
      sampleTree=hclust(dist(datExpr_), method="average")
      heights=sampleTree$height
      upper_bound=median(heights)+3*mad(heights, constant=1)
      # Remove discrete samples
      clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
      table(clust)
      keepSamples=clust==1
      print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
    } else {
      keepSamples=rep(TRUE, nrow(datExpr_))
      print("No filter.")
    }
    AllkeepSamples=c(AllkeepSamples, keepSamples)
  }
  # Filter out discrete samples
  datExpr0=datExpr0[AllkeepSamples,]
  sample_keep[[i]]=rownames(datExpr0)
  
  setTxtProgressBar(processbar, i)
}
close(processbar)

### Filter by taking all celltypes as a whole but remove mixed/lin_infidel=unassigned cells
pb_obj=edgeR::Seurat2PB(subset(object, CD8T_NK_subsets2_rough!="Unassigned"), sample="age", cluster="donor_id")
mat=as.data.frame(pb_obj$counts)
data=log2(mat+1)
datExpr0=as.data.frame(t(data))
ages=levels(as.factor(gsub("_cluster.*","",colnames(mat))))
AllkeepSamples=c()
for (j in 1:length(ages)) {
  sample_of_this_age=grepl(paste0("^",ages[j],"_"), rownames(datExpr0))
  datExpr_=datExpr0[sample_of_this_age,]
  if (nrow(datExpr_)>=3) {
    sampleTree=hclust(dist(datExpr_), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
  } else {
    keepSamples=rep(TRUE, nrow(datExpr_))
  }
  AllkeepSamples=c(AllkeepSamples, keepSamples)
}
datExpr0=datExpr0[AllkeepSamples,]
sample_keep[[length(celltypes)+1]]=rownames(datExpr0)

### Output the sample names that will keep for further analysis
names(sample_keep)=c(celltypes, "All_except_unassigned")
saveRDS(sample_keep, "~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA_rough.rds")

#####################################



### Visualize the outliers of sample to be removed by WGCNA analysis at 3 different annotation levels
#####################################
###

### Detailed level
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
donorid_per_age=table(object[[]]$age, object[[]]$donor_id)
donorid_per_age=as.data.frame(donorid_per_age)
donorid_per_age=subset(donorid_per_age, Freq!=0)
donorid_per_age=donorid_per_age[,c(1,2)]
colnames(donorid_per_age)=c("age","donor_id")
sample_per_age=lapply(split(donorid_per_age$donor_id, donorid_per_age$age), length)
names(sample_per_age)
length(sample_per_age)

sample_keep_detailed=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA.rds")
celltypes=names(sample_keep_detailed)
No_to_remove_allcelltypes=list()
for (i in 1:length(celltypes)) {
  SampleKept_FORcurrentcelltype=sample_keep_detailed[[i]]
  df_=data.frame(age=as.factor(gsub("_cluster.*","",SampleKept_FORcurrentcelltype)), donor_id=gsub(".*_cluster","",SampleKept_FORcurrentcelltype))
  SampleKept_FORcurrentcelltype=lapply(split(df_$donor_id, df_$age), length)
  no_to_remove=list()
  for (age_ in names(sample_per_age)) {
    no_to_remove[[age_]]=sample_per_age[[age_]]-SampleKept_FORcurrentcelltype[[age_]]
  }
  no_to_remove=unlist(no_to_remove)
  No_to_remove_allcelltypes[[i]]=data.frame(age=names(no_to_remove),
                                            no_removed=no_to_remove,
                                            celltype=celltypes[i])
}
names(No_to_remove_allcelltypes)=celltypes
df_no_to_remove_detailed=data.table::rbindlist(No_to_remove_allcelltypes)
df_no_to_remove_detailed$level="Detailed"

### Intermediate level
sample_keep_intermediate=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA_intermediate.rds")
celltypes=names(sample_keep_intermediate)
No_to_remove_allcelltypes=list()
for (i in 1:length(celltypes)) {
  SampleKept_FORcurrentcelltype=sample_keep_intermediate[[i]]
  df_=data.frame(age=as.factor(gsub("_cluster.*","",SampleKept_FORcurrentcelltype)), donor_id=gsub(".*_cluster","",SampleKept_FORcurrentcelltype))
  SampleKept_FORcurrentcelltype=lapply(split(df_$donor_id, df_$age), length)
  no_to_remove=list()
  for (age_ in names(sample_per_age)) {
    no_to_remove[[age_]]=sample_per_age[[age_]]-SampleKept_FORcurrentcelltype[[age_]]
  }
  no_to_remove=unlist(no_to_remove)
  No_to_remove_allcelltypes[[i]]=data.frame(age=names(no_to_remove),
                                            no_removed=no_to_remove,
                                            celltype=celltypes[i])
}
names(No_to_remove_allcelltypes)=celltypes
df_no_to_remove_intermediate=data.table::rbindlist(No_to_remove_allcelltypes)
df_no_to_remove_intermediate$level="Intermediate"

### Rough level
sample_keep_rough=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA_rough.rds")
celltypes=names(sample_keep_rough)
No_to_remove_allcelltypes=list()
for (i in 1:length(celltypes)) {
  SampleKept_FORcurrentcelltype=sample_keep_rough[[i]]
  df_=data.frame(age=as.factor(gsub("_cluster.*","",SampleKept_FORcurrentcelltype)), donor_id=gsub(".*_cluster","",SampleKept_FORcurrentcelltype))
  SampleKept_FORcurrentcelltype=lapply(split(df_$donor_id, df_$age), length)
  no_to_remove=list()
  for (age_ in names(sample_per_age)) {
    no_to_remove[[age_]]=sample_per_age[[age_]]-SampleKept_FORcurrentcelltype[[age_]]
  }
  no_to_remove=unlist(no_to_remove)
  No_to_remove_allcelltypes[[i]]=data.frame(age=names(no_to_remove),
                                            no_removed=no_to_remove,
                                            celltype=celltypes[i])
}
names(No_to_remove_allcelltypes)=celltypes
df_no_to_remove_rough=data.table::rbindlist(No_to_remove_allcelltypes)
df_no_to_remove_rough$level="Rough"

### Save all the no_to_remove
df_no_to_remove=data.table::rbindlist(list(df_no_to_remove_detailed, df_no_to_remove_intermediate, df_no_to_remove_rough))
write.table(df_no_to_remove, "~/Project_PBMCage/Results/sample_filter_No_PerAgePerCelltypePerLevel.txt", sep="\t")

### Visualize the removal
donorid_per_age=table(object[[]]$age, object[[]]$donor_id)
donorid_per_age=as.data.frame(donorid_per_age)
donorid_per_age=subset(donorid_per_age, Freq!=0)
donorid_per_age=donorid_per_age[,c(1,2)]
colnames(donorid_per_age)=c("age","donor_id")
sample_per_age=lapply(split(donorid_per_age$donor_id, donorid_per_age$age), length)
df_no_total=data.frame(age=names(unlist(sample_per_age)), no_all=unlist(sample_per_age))
df_no_to_remove=read.delim("~/Project_PBMCage/Results/sample_filter_No_PerAgePerCelltypePerLevel.txt")
# detailed
df_to_draw=subset(df_no_to_remove, level=="Detailed")
df_to_draw=merge(df_no_total, df_to_draw, by="age")
df_to_draw[df_to_draw==0]=NA
df_to_draw$percent_removed=df_to_draw$no_removed/df_to_draw$no_all
library(ggplot2)
plot_d_notoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=no_removed)) +
  geom_point() +
  scale_size_continuous(name="N removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))
plot_d_percenttoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=percent_removed)) +
  geom_point() +
  scale_size_continuous(name="% removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

# intermediate
df_to_draw=subset(df_no_to_remove, level=="Intermediate")
df_to_draw=merge(df_no_total, df_to_draw, by="age")
df_to_draw[df_to_draw==0]=NA
df_to_draw$percent_removed=df_to_draw$no_removed/df_to_draw$no_all
library(ggplot2)
plot_int_notoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=no_removed)) +
  geom_point() +
  scale_size_continuous(name="N removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))
plot_int_percenttoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=percent_removed)) +
  geom_point() +
  scale_size_continuous(name="% removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

# rough
df_to_draw=subset(df_no_to_remove, level=="Rough")
df_to_draw=merge(df_no_total, df_to_draw, by="age")
df_to_draw[df_to_draw==0]=NA
df_to_draw$percent_removed=df_to_draw$no_removed/df_to_draw$no_all
library(ggplot2)
plot_r_notoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=no_removed)) +
  geom_point() +
  scale_size_continuous(name="N removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))
plot_r_percenttoremove=
  ggplot(df_to_draw, aes(x=age, y=celltype, size=percent_removed*100)) +
  geom_point() +
  scale_size_continuous(name="% removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

pdf("~/Project_PBMCage/Plots/Sample_outlier_removed.pdf", height=20, width=30)
cowplot::plot_grid(plotlist=list(plot_d_notoremove, plot_d_percenttoremove), ncol=2, align="hv")
cowplot::plot_grid(plotlist=list(plot_int_notoremove, plot_int_percenttoremove), ncol=2, align="hv")
cowplot::plot_grid(plotlist=list(plot_r_notoremove, plot_r_percenttoremove), ncol=2, align="hv")
dev.off()


#####################################



### Determine individual outlier(s) from each age by WGCNA-like celltype freq
#####################################
###

### Construct the Celltype Freq (at the detailed level) - Donor_id dataframe
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets2
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
data_short=tidyr::pivot_wider(data[,c(1,2,3)], names_from="CellType", values_from="Freq")
data_short[1:5,1:5]
data_short_wAge=merge(data_short, ID_Age_match[,c(1,2)], by="ID")
dim(data_short_wAge) # check

library(WGCNA)
ages=levels(as.factor(data_short_wAge$Age))
Samples_removed_EachAge=Samples_keep=c()
for (j in 1:length(ages)) {
  data_freq=subset(data_short_wAge, Age==ages[j])
  data_freq=data_freq[,1:(ncol(data_freq)-1)]
  if (nrow(data_freq)>=3) {
    sampleTree=hclust(dist(data_freq), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    # Remove discrete samples
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
    print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
  } else {
    keepSamples=rep(TRUE, nrow(data_freq))
    print("No filter.")
  }
  Samples_removed_EachAge[j]=sum(as.integer(keepSamples==FALSE))
  Samples_keep=c(Samples_keep, as.character(data_freq$ID[keepSamples]))
}
Sample_removed_df_detailed=data.frame(age=ages, no_of_removed=Samples_removed_EachAge, level="Detailed")
saveRDS(Samples_keep,"~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ.rds")

### Construct the Celltype Freq (at the intermediate level) - Donor_id dataframe
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets2_inter
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
data_short=tidyr::pivot_wider(data[,c(1,2,3)], names_from="CellType", values_from="Freq")
data_short[1:5,1:5]
data_short_wAge=merge(data_short, ID_Age_match[,c(1,2)], by="ID")
dim(data_short_wAge) # check

library(WGCNA)
ages=levels(as.factor(data_short_wAge$Age))
Samples_removed_EachAge=Samples_keep=c()
for (j in 1:length(ages)) {
  data_freq=subset(data_short_wAge, Age==ages[j])
  data_freq=data_freq[,1:(ncol(data_freq)-1)]
  if (nrow(data_freq)>=3) {
    sampleTree=hclust(dist(data_freq), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    # Remove discrete samples
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
    print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
  } else {
    keepSamples=rep(TRUE, nrow(data_freq))
    print("No filter.")
  }
  Samples_removed_EachAge[j]=sum(as.integer(keepSamples==FALSE))
  Samples_keep=c(Samples_keep, as.character(data_freq$ID[keepSamples]))
}
Sample_removed_df_inter=data.frame(age=ages, no_of_removed=Samples_removed_EachAge, level="Intermediate")
saveRDS(Samples_keep,"~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ_intermediate.rds")

### Construct the Celltype Freq (at the rough level) - Donor_id dataframe
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
ID_=object$donor_id
Age_=object$age
CellType_=object$CD8T_NK_subsets2_rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt")
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
data_short=tidyr::pivot_wider(data[,c(1,2,3)], names_from="CellType", values_from="Freq")
data_short[1:5,1:5]
data_short_wAge=merge(data_short, ID_Age_match[,c(1,2)], by="ID")
dim(data_short_wAge) # check

library(WGCNA)
ages=levels(as.factor(data_short_wAge$Age))
Samples_removed_EachAge=Samples_keep=c()
for (j in 1:length(ages)) {
  data_freq=subset(data_short_wAge, Age==ages[j])
  data_freq=data_freq[,1:(ncol(data_freq)-1)]
  if (nrow(data_freq)>=3) {
    sampleTree=hclust(dist(data_freq), method="average")
    heights=sampleTree$height
    upper_bound=median(heights)+3*mad(heights, constant=1)
    # Remove discrete samples
    clust=cutreeStatic(sampleTree, cutHeight=upper_bound, minSize=1)
    table(clust)
    keepSamples=clust==1
    print(paste0("cutHeight=",round(upper_bound,2), ";  ",sum(as.integer(keepSamples==FALSE))," sample(s) got removed."))
  } else {
    keepSamples=rep(TRUE, nrow(data_freq))
    print("No filter.")
  }
  Samples_removed_EachAge[j]=sum(as.integer(keepSamples==FALSE))
  Samples_keep=c(Samples_keep, as.character(data_freq$ID[keepSamples]))
}
Sample_removed_df_rough=data.frame(age=ages, no_of_removed=Samples_removed_EachAge, level="Rough")
saveRDS(Samples_keep,"~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ_rough.rds")

Sample_removed_df=data.table::rbindlist(list(Sample_removed_df_detailed, Sample_removed_df_inter, Sample_removed_df_rough))
saveRDS(Sample_removed_df, "~/Project_PBMCage/Tempt_RDS/Sample_filter_Freq.rds")

#####################################



### Visualize the outliers of sample to be removed by Freq_of_celltype per donor analysis at 3 different annotation levels
#####################################
###

### Visualize the removal
# Plot counts
Sample_removed_df=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_Freq.rds")
Sample_removed_df[Sample_removed_df==0]=NA
library(ggplot2)
plot_notoremove=
  ggplot(Sample_removed_df, aes(x=age, y=level, size=no_of_removed)) +
  geom_point() +
  scale_size_continuous(name="N removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

# Plot freq
Sample_removed_df=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_Freq.rds")
Sample_removed_df=as.data.frame(Sample_removed_df)
# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
donorid_per_age=table(object[[]]$age, object[[]]$donor_id)
donorid_per_age=as.data.frame(donorid_per_age)
donorid_per_age=subset(donorid_per_age, Freq!=0)
donorid_per_age=donorid_per_age[,c(1,2)]
colnames(donorid_per_age)=c("age","donor_id")
sample_per_age=unlist(lapply(split(donorid_per_age$donor_id, donorid_per_age$age), length))
sample_per_age_df=data.frame(age=names(sample_per_age), no_of_all=sample_per_age)
df_combined=merge(Sample_removed_df, sample_per_age_df, by="age")
df_combined[df_combined==0]=NA
plot_percenttoremove=
  ggplot(df_combined, aes(x=age, y=level, size=no_of_removed/no_of_all*100)) +
  geom_point() +
  scale_size_continuous(name="% removed") +
  scale_x_discrete(breaks=c(19:97)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

pdf("~/Project_PBMCage/Plots/Sample_outlier_removed_BasedOnFreq.pdf", height=6, width=30)
cowplot::plot_grid(plotlist=list(plot_notoremove, plot_percenttoremove), ncol=2, align="hv")
dev.off()

#####################################



### Compare the outlier determination by Gene and Freq
#####################################
###

###
Sample_determined_byGene=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_WGCNA.rds")
Sample_determined_byGene=Sample_determined_byGene[["All"]]
Sample_determined_byFreq_d=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ.rds")
Sample_determined_byFreq_inter=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ_intermediate.rds")
Sample_determined_byFreq_r=readRDS("~/Project_PBMCage/Tempt_RDS/Sample_filter_FREQ_rough.rds")

# object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
Sample_all=table(object$age, object$donor_id)
Sample_all=as.data.frame(Sample_all)
Sample_all[Sample_all==0]=NA
Sample_all=na.omit(Sample_all)
colnames(Sample_all)=c("age","donor_id","counts")

byGene=match(Sample_all$donor_id, gsub(".*cluster","",Sample_determined_byGene), nomatch=NA)
byGene[is.na(byGene)]="Removed"; byGene[byGene!="Removed"]=NA
byFreq_d=match(Sample_all$donor_id, Sample_determined_byFreq_d, nomatch=NA)
byFreq_d[is.na(byFreq_d)]="Removed"; byFreq_d[byFreq_d!="Removed"]=NA
byFreq_inter=match(Sample_all$donor_id, Sample_determined_byFreq_inter, nomatch=NA)
byFreq_inter[is.na(byFreq_inter)]="Removed"; byFreq_inter[byFreq_inter!="Removed"]=NA
byFreq_r=match(Sample_all$donor_id, Sample_determined_byFreq_r, nomatch=NA)
byFreq_r[is.na(byFreq_r)]="Removed"; byFreq_r[byFreq_r!="Removed"]=NA

Sample_all$byGene=byGene
Sample_all$byFreq_d=byFreq_d
Sample_all$byFreq_inter=byFreq_inter
Sample_all$byFreq_r=byFreq_r

Sample_all=Sample_all[order(Sample_all$age),]
Sample_all$sampleNo=data.table::rowid(Sample_all$age)
Sample_all_longdf=reshape2::melt(Sample_all, id=colnames(Sample_all)[c(1,2,3,length(colnames(Sample_all)))])
library(ggplot2)
ggplot(subset(Sample_all_longdf, value=="Removed"), aes(x=age, y=sampleNo, shape=variable, size=variable, color=variable)) +
  geom_point() +
  scale_shape_manual(values=c(19,1,1,1)) +
  scale_size_manual(values=c(1,2,3,4)) +
  scale_color_manual(values=c("blue",rep("gray30",3))) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  ylab("# of sample")



#####################################
