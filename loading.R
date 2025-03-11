library(SeuratDisk)
library(SeuratData)
library(Seurat)
library(edgeR)
library(Mfuzz)
library(maSigPro)
library(vcd)
library(ggplot2)
library(ggbeeswarm)
library(scRNAtoolVis)
library(ClusterGVis)
library(GseaVis)
library(pheatmap)
library(dittoSeq)
library(ggsignif)
library(tidyr)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(harmony)
library(ggpubr)
library(patchwork)
library(cowplot)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(CellChat)
library(data.table)
library(paletteer)
library(ComplexHeatmap)
library(SCpubr)
library(RColorBrewer)

col.p=scater:::.get_palette("tableau20")
col.p30=ggpubr::get_palette("locuszoom", 30)
col.p25=ggpubr::get_palette("locuszoom", 25)
col.p10=ggpubr::get_palette("locuszoom", 10)
col.pMedium=scater:::.get_palette("tableau10medium")
col.pDark=scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight=scater:::.get_palette("tableau20")[2*(1:10)]
col.sex=scater:::.get_palette("tableau20")[c(1,11)]
col.age=paletteer::paletteer_c("ggthemes::Gray", 78)

### DotPlot
DotPlot(sub_obj_try, 
        features=c("HOPX","PDE4D","IGHE","SELL","EMP3","CIB1","PASP","CD72","DAPP1","LTB","HCK",
                   "ZEB2","RHOB","TNFRSF1B","FCRL3","FCRL5","FCR","MPP6","TAGLN2","IGHA2","AHNAK",
                   "S100A4","CRIP2","ITGB1","JCHAIN","VIM","PLPP5","FCER2","IL4R","CRIP1","LGALS1",
                   "IGHA1","CTSH","S100A10","IGHG2","VPREB3","PPP1R14A","PCDH9","PLD4","IGHM","MT-ATP8",
                   "IGHD","SOX4","AL139020.1","IGLL5","TCL1A"),
        group.by="myannotation",
        cols=c("lightblue", "red")
) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        axis.title=element_blank(),
        
        legend.margin=margin(4,4,4,4),
        legend.position="top",
        legend.direction="horizontal", 
        legend.box="horizontal",
        legend.background=element_rect(fill="transparent",
                                       linewidth=0.5,
                                       linetype="solid",
                                       colour="black"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10)) +
  guides(size=guide_legend(title.position="left", 
                           title.hjust=0.5,
                           label.position="bottom",
                           title="Percent Expressed"))

### Scatter Plot
ggplot(data, aes(x=Age, y=Freq, color=CellType)) +
  geom_point(size=1, shape=1, alpha=0.3) +
  scale_color_manual(values=col.p10) +
  guides(color=guide_legend(nrow=1, override.aes=list(alpha=1, size=2))) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=12),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  labs(title="Percent of B memory cell subsets")+
  ylab("Cell percent (%)") +
  geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
  stat_cor(show.legend=FALSE, size=3.5)

### Convert to Replicates
ss.edesign$Replicates=match(as.factor(ss.edesign$Time), levels(as.factor(ss.edesign$Time)))
ss.edesign$Replicates=as.factor(ss.edesign$Replicates)