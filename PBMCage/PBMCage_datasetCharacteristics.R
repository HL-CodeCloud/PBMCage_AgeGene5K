

### Rb and mt gene control
#####################################
###

###
library(Seurat)
object_PBMCage=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
object_PBMCage_df=data.frame(dataset="Yazar et al. (2022) dataset", percent.rb=object_PBMCage$percent.rb, percent.mt=object_PBMCage$percent.mt)

object_17PBMC=readRDS("~/Project_PBMCage/GSE213516_17PBMC/17_pbmc_annotated.rds")
object_17PBMC_df=data.frame(dataset="GSE213516", percent.rb=object_17PBMC$percent.rb, percent.mt=object_17PBMC$percent.mt)

object_SC=readRDS("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated.rds")
object_SC_df=data.frame(dataset="SC cohort", percent.rb=object_SC$percent.rb, percent.mt=object_SC$percent.mt)

object_CovidCtrl=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
object_CovidCtrl_df=data.frame(dataset="GSE158055", percent.rb=object_CovidCtrl$percent.rb, percent.mt=object_CovidCtrl$percent.mt)

object_immunity=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/all_pbmcs_metadata.csv")
object_immunity_df=data.frame(dataset="Terekhova et al. (2023) dataset", percent.rb=object_immunity$percent.ribo, percent.mt=object_immunity$percent.mt)

### 5 datasets combined
object_df_5all=data.table::rbindlist(list(object_PBMCage_df, object_17PBMC_df, object_SC_df, object_CovidCtrl_df, object_immunity_df))
data.table::fwrite(object_df_5all, "~/Project_PBMCage/Tempt_RDS/Confirm_5_datasets_mtANDrb_control.csv.gz", sep="\t")

### Plot
DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Confirm_5_datasets_mtANDrb_control.csv.gz", sep="\t")
colnames(DF)=c("dataset","ribosomal genes","mitochondrial genes")
DF=DF %>% tidyr::pivot_longer(cols=c("ribosomal genes","mitochondrial genes"), names_to="parameter", values_to="percent")
DF$dataset=forcats::fct_relevel(DF$dataset, c("Yazar et al. (2022) dataset","Terekhova et al. (2023) dataset","GSE213516","GSE158055","SC cohort"))
plot_=
  ggplot(DF, aes(x=dataset, y=percent)) +
  facet_wrap(~parameter, scales="free") +
  ggrastr::geom_violin_rast(linewidth=0.3) +
  theme_classic() +
  labs(x=NULL, y="Percent (%)") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        strip.text=element_text(size=10),
        title=element_text(size=11)) +
  scale_x_discrete(labels=label_wrap_gen(10))
  
pdf("~/Project_PBMCage/Plots/Confirm_5_datasets_mtANDrb.pdf", height=2, width=8.5)
plot(plot_)
dev.off()

#####################################



### Draw DimPlot
#####################################
###

library(Seurat)
library(paletteer)
library(ggplot2)

###

object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

full_plot=
  DimPlot(object, group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Blue", object[[]] %>% subset(Annot.rough=="B cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))),
                 paletteer_c("ggthemes::Orange", object[[]] %>% subset(Annot.rough=="CD4T cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))),
                 paletteer_c("ggthemes::Green", object[[]] %>% subset(Annot.rough=="CD8T cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))),
                 paletteer_c("ggthemes::Red", object[[]] %>% subset(Annot.rough=="DCs") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))),
                 rev(paletteer_c("grDevices::Purp", object[[]] %>% subset(Annot.rough=="Monocytes") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x))))),
                 paletteer_c("ggthemes::Brown", object[[]] %>% subset(Annot.rough=="NK cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))),
                 rev(paletteer_c("grDevices::Magenta", object[[]] %>% subset(Annot.rough=="Other cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x))))),
                 paletteer_c("ggthemes::Gray", object[[]] %>% subset(Annot.rough=="OtherT") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
  theme_bw() +
  theme(axis.title=element_text(size=15),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) +
  labs(title=NULL) +
  guides(color=guide_legend(ncol=3, bycol=TRUE))
pdf("~/Project_PBMCage/Plots/DimPlot.pdf", height=7, width=7)
full_plot + NoLegend()
dev.off()

# legend for B cells
lg_Bcells=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="B cells"), group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Blue", object[[]] %>% subset(Annot.rough=="B cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="B cells", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_Bcells.pdf", height=3, width=2.25)
ggpubr::as_ggplot(lg_Bcells)
dev.off()

# legend for CD4T
lg_CD4T=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="CD4T cells"), group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Orange", object[[]] %>% subset(Annot.rough=="CD4T cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="CD4T cells", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_CD4T.pdf", height=3.6, width=2.25)
ggpubr::as_ggplot(lg_CD4T)
dev.off()

# legend for CD8T
lg_CD8T=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="CD8T cells"), group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Green", object[[]] %>% subset(Annot.rough=="CD8T cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="CD8T cells", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_CD8T.pdf", height=3.6, width=2.25)
ggpubr::as_ggplot(lg_CD8T)
dev.off()

# legend for DCs
lg_DCs=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="DCs"), group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Red", object[[]] %>% subset(Annot.rough=="DCs") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="DCs", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_DCs.pdf", height=1.8, width=2.25)
ggpubr::as_ggplot(lg_DCs)
dev.off()

# legend for Monocytes
lg_Mono=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="Monocytes"), group.by="Annot.detailed", raster=TRUE, 
          cols=rev(paletteer_c("grDevices::Purp", object[[]] %>% subset(Annot.rough=="Monocytes") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="Monocytes", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_Mono.pdf", height=2, width=2.25)
ggpubr::as_ggplot(lg_Mono)
dev.off()

# legend for NK
lg_NK=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="NK cells"), group.by="Annot.detailed", raster=TRUE, 
          cols=c(paletteer_c("ggthemes::Brown", object[[]] %>% subset(Annot.rough=="NK cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="NK cells", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_NKcells.pdf", height=3.2, width=2.4)
ggpubr::as_ggplot(lg_NK)
dev.off()

# legend for Other cells
lg_Othercells=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="Other cells"), group.by="Annot.detailed", raster=TRUE, 
          cols=rev(paletteer_c("grDevices::Magenta", object[[]] %>% subset(Annot.rough=="Other cells") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x)))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="Other cells", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_Othercells.pdf", height=2.4, width=2.25)
ggpubr::as_ggplot(lg_Othercells)
dev.off()

# legend for OtherT
lg_OtherT=ggpubr::get_legend(
  DimPlot(object %>% subset(Annot.rough=="OtherT"), group.by="Annot.detailed", raster=TRUE, 
          cols=paletteer_c("ggthemes::Gray", object[[]] %>% subset(Annot.rough=="OtherT") %>% dplyr::select(Annot.detailed) %>% sapply(., function(x) length(unique(x))))) +
    theme_bw() +
    theme(legend.text=element_text(size=9),
          axis.title=element_text(size=15),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(title=NULL) +
    guides(color=guide_legend(ncol=1, title="OtherT", title.theme=element_text(size=10, face="bold.italic")))
)
pdf("~/Project_PBMCage/Plots/DimPlot_legend_OtherT.pdf", height=2, width=2.25)
ggpubr::as_ggplot(lg_OtherT)
dev.off()

#####################################




### Marker genes with all the cell types
#####################################
###

###
library(Seurat)
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Find markers for Annot.rough
Idents(object)="Annot.rough"
markers_rough=FindAllMarkers(object)
saveRDS(markers_rough, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds")

### Find markers for each Annot.inter
markers_inter=list()
celltype.rough=unique(object$Annot.rough)
for (i in 1:length(celltype.rough)) {
  subset_=subset(object, Annot.rough==celltype.rough[i])
  Idents(subset_)="Annot.inter"
  markers_inter[[i]]=FindAllMarkers(subset_)
}
names(markers_inter)=celltype.rough
saveRDS(markers_inter, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_inter.rds")

### Find markers for each Annot.detailed
markers_detailed=list()
celltype.inter=unique(object$Annot.inter)
for (i in 1:length(celltype.inter)) {
  subset_=subset(object, Annot.inter==celltype.inter[i])
  Idents(subset_)="Annot.detailed"
  markers_detailed[[i]]=FindAllMarkers(subset_)
}
names(markers_detailed)=celltype.inter
saveRDS(markers_detailed, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed.rds")

### Make celltype-specific gene lists
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds")
markers_inter=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_inter.rds")
markers_detailed=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed.rds")
# at the rough level
markers_rough_perCluster=markers_rough %>% 
  split(.$cluster) %>%
  lapply(., function(df) df %>% 
           subset(p_val_adj<0.05 & avg_log2FC>0) %>%
           select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
           mutate(celltype.level="Annot.rough")) %>%
  data.table::rbindlist(.) %>%
  mutate(Annot.detailed=NA, Annot.inter=NA) %>%
  dplyr::rename(Annot.rough=cluster) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))
markers_rough_perCluster %>% split(.$Annot.rough) %>% lapply(., nrow) # check if all the rough celltypes have their marker genes

# at the inter level
markers_inter_perCluster=markers_inter %>% 
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
               mutate(celltype.level="Annot.inter")) %>%
      data.table::rbindlist(.)
  })
lapply(lapply(markers_inter_perCluster, function(DF) DF %>% split(.$cluster)), function(DF) lapply(DF, nrow)) # check if all the inter celltypes have their marker genes
markers_inter_perCluster=lapply(1:length(markers_inter_perCluster), 
                                   function(idx) markers_inter_perCluster[[idx]] %>% 
                                     mutate(Annot.rough=names(markers_inter_perCluster)[idx])) %>%
  data.table::rbindlist(.) %>%
  dplyr::rename(Annot.inter=cluster) %>%
  mutate(Annot.detailed=NA) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))

# at the detailed level
inter_to_remove=names(unlist(lapply(markers_detailed, nrow))[unlist(lapply(markers_detailed, nrow))==0]) # determine the skipped celltypes at the inter level
markers_detailed_perCluster=markers_detailed %>% 
  .[!(names(.) %in% inter_to_remove)] %>% # remove the skipped celltypes
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
               mutate(celltype.level="Annot.detailed")) %>%
      data.table::rbindlist(.)
  })
lapply(lapply(markers_detailed_perCluster, function(DF) DF %>% split(.$cluster)), function(DF) lapply(DF, nrow)) # check if all the detailed celltypes have their marker genes
markers_detailed_perCluster=lapply(1:length(markers_detailed_perCluster), 
                                   function(idx) markers_detailed_perCluster[[idx]] %>% 
                                     mutate(Annot.inter=names(markers_detailed_perCluster)[idx])) %>%
  data.table::rbindlist(.) %>%
  mutate(Annot.rough=gsub("\\..*","",Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="B","B cells",
                            ifelse(Annot.rough=="CD4T","CD4T cells",
                                   ifelse(Annot.rough=="CD8T","CD8T cells",
                                          ifelse(Annot.rough=="DC","DCs",
                                                 ifelse(Annot.rough=="Mono","Monocytes",
                                                        ifelse(Annot.rough=="Other","Other cells",
                                                               ifelse(Annot.rough=="OtherT","OtherT",
                                                                      ifelse(Annot.rough=="NK","NK cells",NA))))))))) %>%
  dplyr::rename(Annot.detailed=cluster) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))

# merge at all 3 levels
merged_results=data.table::rbindlist(list(markers_rough_perCluster, markers_inter_perCluster, markers_detailed_perCluster))
data.table::fwrite(merged_results, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")

# ### Plot the marker genes for Annot.rough
# markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds")
# allgenes=markers_rough$gene
# # get the marker genes
# marker_genes=markers_rough %>% subset(! gene %in% allgenes[grepl("^RP",allgenes)]) %>%
#   subset(p_val_adj<0.05) %>%
#   group_by(cluster) %>% slice_max(avg_log2FC, n=10) %>% ungroup() %>%
#   select(gene) %>% unlist() %>% unname()
# # plot
# source("~/Rscripts/AverageHeatmap2.R")
# Idents(object)="Annot.rough"
# plot_=
#   AverageHeatmap2(object=object,
#                   markerGene=marker_genes,
#                   row_title=NULL,
#                   clusterAnnoName=FALSE,
#                   fontsize=5,
#                   annot_fontsize=8,
#                   htCol=c("lightblue", "white", "coral3"),
#                   width=4,
#                   height=16)
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_averageExpr_rough.pdf", height=7, width=3)
# draw(plot_)
# dev.off()
# 
#####################################




### Marker genes with updated annotation (excluding eryth/mast and merging some populations)
#####################################
###

###
library(Seurat)
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q","Mono.inter.GNLY"), "Mono.inter", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

### Find markers for Annot.rough
Idents(object)="Annot.rough"
markers_rough=FindAllMarkers(object)
saveRDS(markers_rough, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough_updated.rds")

### Find markers for each Annot.inter
markers_inter=list()
celltype.rough=unique(object$Annot.rough)
for (i in 1:length(celltype.rough)) {
  subset_=subset(object, Annot.rough==celltype.rough[i])
  Idents(subset_)="Annot.inter"
  markers_inter[[i]]=FindAllMarkers(subset_)
}
names(markers_inter)=celltype.rough
saveRDS(markers_inter, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_inter_updated.rds")

### Find markers for each Annot.detailed
markers_detailed=list()
celltype.inter=unique(object$Annot.inter)
for (i in 1:length(celltype.inter)) {
  subset_=subset(object, Annot.inter==celltype.inter[i])
  Idents(subset_)="Annot.detailed"
  markers_detailed[[i]]=FindAllMarkers(subset_)
}
names(markers_detailed)=celltype.inter
saveRDS(markers_detailed, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed_updated.rds")

# ### Find markers for each Annot.detailed per Annot.rough
# markers_detailed=list()
# celltype.rough=unique(object$Annot.rough)
# for (i in 1:length(celltype.rough)) {
#   subset_=subset(object, Annot.rough==celltype.rough[i])
#   Idents(subset_)="Annot.detailed"
#   markers_detailed[[i]]=FindAllMarkers(subset_)
#   print("Done.")
# }
# names(markers_detailed)=celltype.rough
# saveRDS(markers_detailed, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed_perRough_updated.rds")

### Make celltype-specific gene lists
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough_updated.rds")
markers_inter=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_inter_updated.rds")
markers_detailed=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed_updated.rds")
# at the rough level
markers_rough_perCluster=markers_rough %>% 
  split(.$cluster) %>%
  lapply(., function(df) df %>% 
           subset(p_val_adj<0.05 & avg_log2FC>0) %>%
           select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
           mutate(celltype.level="Annot.rough")) %>%
  data.table::rbindlist(.) %>%
  mutate(Annot.detailed=NA, Annot.inter=NA) %>%
  dplyr::rename(Annot.rough=cluster) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))
markers_rough_perCluster %>% split(.$Annot.rough) %>% lapply(., nrow) # check if all the rough celltypes have their marker genes

# at the inter level
markers_inter_perCluster=markers_inter %>% 
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
               mutate(celltype.level="Annot.inter")) %>%
      data.table::rbindlist(.)
  })
lapply(lapply(markers_inter_perCluster, function(DF) DF %>% split(.$cluster)), function(DF) lapply(DF, nrow)) # check if all the inter celltypes have their marker genes
markers_inter_perCluster=lapply(1:length(markers_inter_perCluster), 
                                function(idx) markers_inter_perCluster[[idx]] %>% 
                                  mutate(Annot.rough=names(markers_inter_perCluster)[idx])) %>%
  data.table::rbindlist(.) %>%
  dplyr::rename(Annot.inter=cluster) %>%
  mutate(Annot.detailed=NA) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))

# at the detailed level
inter_to_remove=names(unlist(lapply(markers_detailed, nrow))[unlist(lapply(markers_detailed, nrow))==0]) # determine the skipped celltypes at the inter level
markers_detailed_perCluster=markers_detailed %>% 
  .[!(names(.) %in% inter_to_remove)] %>% # remove the skipped celltypes
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
               mutate(celltype.level="Annot.detailed")) %>%
      data.table::rbindlist(.)
  })
lapply(lapply(markers_detailed_perCluster, function(DF) DF %>% split(.$cluster)), function(DF) lapply(DF, nrow)) # check if all the detailed celltypes have their marker genes
markers_detailed_perCluster=lapply(1:length(markers_detailed_perCluster), 
                                   function(idx) markers_detailed_perCluster[[idx]] %>% 
                                     mutate(Annot.inter=names(markers_detailed_perCluster)[idx])) %>%
  data.table::rbindlist(.) %>%
  mutate(Annot.rough=gsub("\\..*","",Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="B","B cells",
                            ifelse(Annot.rough=="CD4T","CD4T cells",
                                   ifelse(Annot.rough=="CD8T","CD8T cells",
                                          ifelse(Annot.rough=="DC","DCs",
                                                 ifelse(Annot.rough=="Mono","Monocytes",
                                                        ifelse(Annot.rough=="Other","Other cells",
                                                               ifelse(Annot.rough=="OtherT","OtherT",
                                                                      ifelse(Annot.rough=="NK","NK cells",NA))))))))) %>%
  dplyr::rename(Annot.detailed=cluster) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough","pct.1","pct.2")))

# merge at all 3 levels
merged_results=data.table::rbindlist(list(markers_rough_perCluster, markers_inter_perCluster, markers_detailed_perCluster))
data.table::fwrite(merged_results, "~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t")

#####################################




### Plot marker genes from experimental/review data
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(B_cell=c("CD19","MS4A1","CD79A","CD79B","CD27","CD38","PTPRC","CD74","CD22"),
                        CD4T_cell=c("CD4","IL7R","CD3D","CCR7","LTB","SELL","CD3E","FOXP3","CD27"),
                        CD8T_cell=c("CD8A","CD3D","CD8B","CD3E","CCR7","CD28","CD4"),
                        DCs=c("FCER1A","CD1C","ITGAX","CD83","CLEC10A","IL3RA","CD86","CST3","LILRA4"),
                        Mono=c("CD14","FCGR3A","LYZ","HLA-DRA","HLA-DRB1","HLA-DRB5","S100A8","CD163","FCGR1A","CST3"),
                        NK=c("NKG7","GNLY","KLRF1","NCAM1","FCGR3A","B3GAT1","KLRD1"),
                        Othercell=c("IL7R","KLRB1","CD200R1","PTGDR2", # ILC
                                    "PPBP","ITGA2B","SELP","CD40LG","GP1BA","ITGB3","CD63","SRC","CCL5","ICAM2", # platelet
                                    "CD34","IFNAR1","GYPA","CD38","ITGA2B"), # progenitors
                        OtherT=c("NCAM1","CD3D","CD3E","CD8A","FCGR3A","CD2","KLRB1","KLRD1", # NKT
                                 "TRGV9","TRDC","TRDV2","CD8B","TRGC1","CCL5","CD2","CD3D","CD3E", # gdT
                                 "SLC4A10","KLRB1","KLRB1","TRAV1-2","CCR6","DPP4","GZMK","IL17RA") # MAIT
                        )
# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")
marker_genes=split(merged_results$gene, merged_results$Annot.rough)
marker_genes=lapply(1:length(marker_genes), function(idx) intersect(marker_genes[[idx]], experiment_markers[[idx]])) %>%
  Reduce(c, .) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
object_donor=Seurat::AggregateExpression(object, group.by="Annot.rough", return.seurat=T)
object_donor$Annot.rough=colnames(object_donor)
expr_matrix=Seurat::FetchData(object_donor, vars=marker_genes, layer="data")
# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=T, cluster_columns=F, show_row_dend=F,
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2, 0, 2), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9), column_names_rot=45,
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(5,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_averageExpr_rough.pdf", height=2, width=15)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in B cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(B.naive=c(
                                  # "IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3",
                                  "CD38","CD27","IGHD","IGHM","CD24"), # CD38+ CD27− IGHD+ IGHMlo CD24lo
                        B.naive.lambda=c("IGHM","IGHD","MS4A1","TCL1A","IL4R","IGLC6","IGLC7","IGLC3"),
                        B.mem=c(
                                # "MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781",
                                "CD38","CD27","CD24","CR2"), # CD38− CD27− CD24+ CR2+; CD38lo CD27+ CD24+ CR2−
                        B.mem.CD27neg=c("CD19","MS4A1","CD27","IGHD","IGHM",
                                        "CD1C","CD38","CD27","MS4A1","CR2","FCRL5","FCRL4","IL10RB","CXCR3","TBX21","ZEB2"), # CD1C+ CD38− CD27− CD20hi CD21lo FCRL5+ FCRL4+ IL10RB+ CXCR3+ TBX21+ ZEB2+
                        B.trans=c(
                                  # "CD19","IGHD","CD27","IGHM","IGHD","MME",
                                  "CD38","CD27","IGHD","IGHM","CD24","MME"), # CD38hi CD27- IGHD+ IGHM+ CD24+ MME+
                        B.inter=c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B"), # this's not the activated B cells in the Immunity paper which are CXCR5+ JUN+ FOS+ NR4A1+
                        B.pb=c("CXCR3","HLA-DRA","KLF4","MKI67","MS4A1"),
                        B.pc=c("PTPRC","CXCR4","IGHA1","IGHG1","IGHA2","IGHG2","IGHG3","IGHM","TNFRSF17","SDC1",
                               "CD38","CD27","CD24","PRDM1","XBP1") # CD38hi CD27hi CD24- PRDM1+ XBP1+)
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="B cells")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="B cells")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data") %>%
  .[c("B.naive.kappa","B.naive.lambda","B.trans","B.inter.IL4R","B.inter.IGHM","B.inter.CD27neg","B.mem.IGHM","B.mem.IGHA","B.mem.CD27neg","B.pb","B.pc"),]
# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-3, 0, 3), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_Bcells.pdf", height=2.5, width=7)
draw(plot_rough)
dev.off()

#####################################




### Plot marker genes from experimental/review data at the detailed level in CD4T cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(CD4T.naive=c(
                                     # "TCF7","CCR7","IL7R","FHIT","LEF1","MAL","NOSIP","LDHB","PIK3IP1",
                                     "CCR7"),
                        CD4T.naive.HNRNPH1=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"),
                        CD4T.naive.COTL1=c("CD44","SELL","IL7R","COTL1","ANXA1","TAGLN2"),
                        CD4T.naive.SOX4=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"),
                        CD4T.Tcm=c(
                                   # "IL7R","ITGB1","LTB","TRAC","AQP3","LDHB","IL32","MAL","CCR7","SELL","GPR183","PIK3IP1",
                                   "CXCR5", # Tfh
                                   "CXCR3","CCR6","KLRB1","GZMK", # Th1
                                   "CCR4","CCR6","GATA3"), # Th2
                        CD4T.Tcm.COTL1=c(
                                         # "CD44","SELL","IL7R","COTL1","ITGB1","LGALS1","S100A11","ANXA1",
                                         "CTLA4","PDCD1","LAG3","HAVCR2","GZMK"),
                        CD4T.Tcm.HLA=c(
                                       # "CD44","SELL","IL7R","HLA-DRA","CD38","IL2RA",
                                       "HLA-DRB1"),
                        CD4T.Tcm.ISG=c("LIMS1","MAL","LTB","AQP3","MX1","ISG15","ISG20"),
                        CD4T.Tcm.HNRNPH1=c("CD44","SELL","IL7R","C1orf56","HNRNPH1","CDC42SE1"),
                        CD4T.Tcm.KLRB1=c(
                                         # "CD44","SELL","IL7R","KLRB1","CTSH","AQP3","TIMP1",
                                         "CXCR3","CCR6","KLRB1","GZMK", # Th1/Th17
                                         "CXCR3","CCR6","RORC","AHR"), # Th17
                        CD4T.Tem.KLRB1=c(
                                         # "IL7R","CCL5","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3","CD44","SELL","IL7R","GZMK","KLRB1",
                                         "PRF1","GNLY","EOMES","TBX21","GZMB","GZMK","CCR7"),
                        CD4T.Tem.MALAT1=c(
                                          # "IL7R","CCL5","GZMK","IL32","GZMA","KLRB1","TRAC","LTB","AQP3","JUND","MALAT1","DDX17",
                                          "CCR4","CCR6","CCR10", # Th22
                                          "PRF1","GNLY","EOMES","TBX21","GZMB","GZMK"), # Temra
                        CD4T.Tfh=c("CXCR5"),
                        CD4T.Th1=c("CXCR3","CCR6","KLRB1","GZMK"), # CXCR3+ CCR6- KLRB1- GZMK+
                        CD4T.Th2=c("CCR4","CCR6","GATA3"), # CCR4+ CCR6− GATA3+
                        CD4T.Th1Th17=c("CXCR3","CCR6","KLRB1","GZMK"), # CXCR3+ CCR6+ KLRB1+ GZMK+
                        CD4T.Th17=c("CXCR3","CCR6","RORC","AHR"), # CXCR3- CCR6+ RORC+ AHR+
                        CD4T.Th22=c("CCR4","CCR6","CCR10"), # CCR4+ CCR6+ CCR10+
                        CD4T.HLADR.Tmem=c("HLA-DRB1"),
                        CD4T.Tex=c("CTLA4","PDCD1","LAG3","HAVCR2","GZMK"),
                        CD4T.cytotoxic=c("PRF1","GNLY","EOMES","TBX21","GZMB","GZMK"),
                        CD4T.Tte=c("PRF1","GNLY","EOMES","TBX21","GZMB","GZMK","CCR7"), # PRF1+ GNLY+ EOMES+ TBX21+ GZMB+ GZMK+ CCR7-
                        CD4T.Ttemra=c("PRF1","GNLY","EOMES","TBX21","GZMB","GZMK"), # PRF1+ GNLY+ EOMES+ TBX21+ GZMB+ GZMK-
                        CD4T.Treg.cytotoxic=c("PRF1","GNLY","EOMES","TBX21","GZMB","FOXP3","IL2RA","IL7R"), # PRF1+ GNLY+ EOMES+ TBX21+ GZMB+ FOXP3+ CD25+ IL7Rlow
                        CD4T.Treg.naive=c("FOXP3","IL2RA","IL7R","CCR7","SELL","LEF1","TCF7"), # FOXP3+ CD25+ IL7Rlow CCR7+ SELL+ LEF1+ TCF7+
                        CD4T.Treg.mem=c("FOXP3","IL2RA","IL7R","IL2RA","CCR4","CTLA4","HLA-DRB1"), # FOXP3+ CD25+ IL7Rlow IL2RA+ CCR4+ CTLA4+ HLA-DRB1+
                        CD4T.Treg.RORC=c("FOXP3","IL2RA","IL7R","KLRB1","RORC","CCR4","CTLA4","HLA-DRB1","IKZF2") # FOXP3+ CD25+ IL7Rlow KLRB1+ RORC+ CCR4+ CTLA4+ HLA-DRB1- IKZF2-
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="CD4T cells")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="CD4T cells")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-SOX4","_SOX4",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-3, 0, 3), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_CD4T.pdf", height=3, width=8.5)
draw(plot_rough)
dev.off()

#####################################




### Plot marker genes from experimental/review data at the detailed level in CD8T cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(CD8T.naive=c("SELL","CCR7"), # CD62L+CCR7+
                        CD8T.naive.LINC02446=c("CD8B","LDHB","LEF1","LINC02446","S100B","ID2","TCF7","VIM","CCR7"),
                        CD8T.naive.HNRNPH1=c("LEF1","SELL","C1orf56","HNRNPH1"),
                        CD8T.naive.TSHZ2=c("TSHZ2","FHIT","RNASET2"),
                        CD8T.Tcm=c("LEF1","TCF7"),
                        CD8T.Tcm.LYAR=c("GZMK","LYAR","NELL2"),
                        CD8T.Tcm.GZMB=c("NKG7","GNLY","GZMH","GZMB"),
                        CD8T.Tcm.GATA3=c(
                                         "CD8B","C1orf162","IL7R","GATA3","YBX3","KRT1","CD8A","CTSW","INPP4B","LTB",
                                         "CCR4","IL2RA","IL4","IL2"), # annotated CD8T.Tcm.CCR4pos in Immunity
                        CD8T.Tcm.MALAT1=c("MALAT1","NEAT1","DDX17"),
                        CD8T.CTL=c(
                                   "GZMH","NKG7","ID2","IFNG","MKI67",
                                   "CD38","HLA-DRB1","MKI67"), # annotated proliferative CD8T in Immunity
                        CD8T.Tem=c("CCL5","GZMB","GZMK"),
                        CD8T.NKT_like=c("GZMB","CD8A","CD8B","NCAM1","TYROBP","ZNF683","IKZF2"),
                        CD8T.Tem.GZMB_FCGR3A=c("FCGR3A","PRF1","GZMB",
                                               "GZMB","ZNF683", # CD8T.Tem.GZMB
                                               "GZMB","IKZF2", # CD8T.Temra
                                               "GZMB","CD8A","CD8B","NCAM1","TYROBP","ZNF683","IKZF2"), # CD8T.NKT_like
                        CD8T.Tem.GZMB_HNRNPH1=c("GZMB","C1orf56","HNRNPH1","CDC42SE1"),
                        CD8T.Tem.GZMK_NKG7=c("GZMB","GNLY","NKG7","ZEB2","GZMA","GZMK"),
                        CD8T.Tem.GZMKhi=c("CCR7","GZMK","LEF1","CD27","CD28",
                                          "GZMK"), # annotated CD8T.Tem.GZMK in Immunity
                        CD8T.Temra=c("GZMB","IKZF2"),
                        CD8T.Tem.GZMB=c("GZMB","ZNF683"),
                        CD8T.Tem.GZMK=c("GZMK"),
                        CD8T.prolif=c("CD38","HLA-DRB1","MKI67"),
                        # CD8T.Trm=c("IL7R","CD69","ITGAE","CCR9","KLRB1","SELL","CCR7"), # SELL−CCR7−
                        CD8T.KLRC2=c("ZNF683","KLRC2","GZMB","LEF1","XCL1", # GZMB-
                                     "KLRC2","GZMB"), # GZMB-
                        CD8T.Tem.ISG=c("MX1","ISG15","ISG20")
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="CD8T cells")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="CD8T cells")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-3, 0, 3), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_CD8T.pdf", height=3, width=10)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in CD8T.Tem.GZMKhi vs. CD8T.Tem.GZMK_NKG7 cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(C0_Tem.like=c("GZMK"),
                        C3_Tte.like=c("GZMB","GNLY","ZEB2","GZMA")
)
experiment_markers_orig=list(`Co-stimulatory M.`=c("CD27","CD28"),
                        `Cytotoxicity M.`=c("PRF1","GZMB"),
                        `Exhaustion/Apoptosis M.`=c("PDCD1","FAS"),
                        `Cytokine`=c("IFNG","TNF"),
                        `Transcription Factors`=c("TBX21","EOMES"),
                        `Surface Inhibitory M.`=c("HAVCR2","LAG3"))
experiment_markers=experiment_markers_orig %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols=colnames(.), names_to="cluster", values_to="gene")
experiment_markers$cluster=rep(names(experiment_markers_orig), times=2)

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object_subset=subset(object, Annot.detailed %in% c("CD8T.Tem.GZMK_NKG7","CD8T.Tem.GZMKhi","CD8T.Tem.GZMB_FCGR3A","CD8T.Tem.GZMB_HNRNPH1"))
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=Reduce(union, experiment_markers), layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
source("~/Rscripts/jjDotPlot2.R")
plot_jjplot=
  jjDotPlot2(object=object_subset,
            markerGene=experiment_markers,
            id="Annot.detailed",
            plot.margin=c(5.2,1,1,1),
            anno=T,
            xtree=F,
            rescale=T,
            rescale.min=0,
            rescale.max=1,
            yPosition=5,
            tile.geom=T,
            textSize=12,
            bar.legendTitle="Mean expression", bar.width=5,
            point.lengdTitle="Fraction of cells (%)") +
  theme(
    axis.text.y=element_text(size=10),
    axis.text.x=element_text(size=10),
    legend.title=element_text(size=10),
    legend.text=element_text(size=9),
    legend.direction="horizontal",
    legend.box="horizontal",
    legend.position="bottom"
  )

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_CD8T.Tem.GZMKhi_vs_CD8T.Tem.GZMK_NKG7.pdf", height=5.2, width=8)
plot_jjplot
dev.off()

# plot
library(ComplexHeatmap)
plot_CD8T.Tem=
  Heatmap(expr_matrix,
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Mean expr.",
                 legend_height=unit(2, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(0, 2), c("white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.5,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_CD8T.Tem=draw(plot_CD8T.Tem, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_CD8T.Tem.GZMKhi_vs_CD8T.Tem.GZMK_NKG7.pdf", height=3, width=10.5)
draw(plot_CD8T.Tem)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in DCs
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(ASDC.mDC=c("PPP1R14A","LILRA4","AXL","IL3RA","SCT","SCN9A","LGMN","CLEC4C","DNASE1L3","GAS6"),
                        ASDC.pDC=c("LILRA4", "CLEC4C", "SCT", "EPHB1", "AXL", "PROC", "LRRC26", "SCN9A", "LTK", "DNASE1L3"),
                        cDC2_1=c("FCER1A", "CD14", "CLEC10A", "CTSS", "ENHO", "CD1C", "MRC1", "FCGR2B", "PID1", "IL13RA1"),
                        cDC2_2=c("FCER1A", "BASP1", "CD1C", "CD74", "CLEC10A", "HLA-DPA1", "ENHO", "HLA-DPB1", "PLD4", "HLA-DQA1","MKI67"),
                        cDC1=c("BLTA", "CADM1", "CD8A", "CLEC9A", "ITGAE", "ITGAX", "LY75", "THBD", "XCR1", "BATF3","ID2","IRF8","ZBTB46"),
                        pDC=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","MZB1","SPIB","IRF4","SMPD3"))

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="DCs")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="DCs")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-3, 0, 3), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_DCs.pdf", height=2, width=8)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in Monocytes
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(Mono.CD14=c("S110A4","S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14","FCGR3A"), # CD14+ FCGR3A-
                        Mono.CD14.HLA=c("HLA-DPA1","HLA-DQA1","HLA-DRB1"),
                        Mono.CD16=c("MS4A7","S100A4","FCGR3A","LST1","AIF1","CTSS","CTSL"), # CD14- FCGR3A+
                        Mono.CD16.GNLY=c("FCGR3A","KLRD1","KLRF1","GNLY","LYZ","NKG7"),
                        Macrophage_like_CD64pos_Mono=c("FCGR3A","FCGR1A","FCGR1B"))

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="Monocytes")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="Monocytes")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2, 0, 2), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_Monocytes.pdf", height=2, width=5)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in NK cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(NK.prolif.MCM=c("FCGR3A","NCAM1","IL2RB","SPON2","FCER1G","GZMB","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MKI67"),
                        NK.prolif.MKI67=c("MKI67"), # MKI67+
                        NK.CD56dim.FCER1G=c("MYOM2","FCER1G","PRF1","SPON2",
                                            "NCAM1"), # annotated NK.CD56dim.CD57int in Immunity
                        NK.CD56dim.HNRNPH1=c("C1orf56","HNRNPH1","CDC42SE1"),
                        NK.CD56dim.FCER1Glo_KLRC2pos=c("KLRC2","CD3E","CCL5","VIM",
                                                       "PRF1","CXCR1","CX3CR1","KLRG1","GZMB","TBX21","JUND","FOSB","ZBTB38"), # annotated NK.CD56dim.CD57+ in Immunity
                        NK.CD56dim.FCER1Glo_KLRC2neg=c("CD3E","VIM","IL32","CD52","KLRC2","FCER1G"),
                        NK.CD56dim.FOS=c("FOS","FOSB","JUN","JUNB",
                                         "PRF1","CXCR1","CX3CR1","KLRG1","GZMB","TBX21","JUND","FOSB"), # annotated NK.CD56dim.CD57low in Immunity
                        NK.CD56dim.ISG=c("MX1","ISG15","ISG20"),
                        NK.CD56hi.FCGR3A=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"),
                        NK.CD56hi.S100B=c("S100B","MYOM2","BEX3","SOX4","ITGA1","KIR2DL4"),
                        NK.CD56hi.CCL3=c("CCL4L2","CCL3","CCL4"),
                        NK.CD56hi.COLT1=c("COTL1","XCL1","LTB","GZMK"),
                        NK.CD56dim.CD57neg=c("PRF1","CXCR1","CX3CR1","KLRG1","JUND","FOSB"), # PRF1- CXCR1- CX3CR1- KLRG1hi JUNDint FOSBint
                        NK.CD56dim.CD57low=c("PRF1","CXCR1","CX3CR1","KLRG1","GZMB","TBX21","JUND","FOSB"), # PRF1low CXCR1low CX3CR1low KLRG1int GZMB+ TBX21+ JUNDhi FOSBhi
                        NK.CD56dim.CD57int=c("PRF1","CXCR1","CX3CR1","KLRG1","GZMB","TBX21","JUND","FOSB","ZBTB38"), # PRF1int CXCR1int CX3CR1int KLRG1low GZMB+ TBX21+ JUND- FOSB- ZBTB38lo
                        NK.CD56dim.CD57hi=c("PRF1","CXCR1","CX3CR1","KLRG1","GZMB","TBX21","JUND","FOSB","ZBTB38") # PRF1hi CXCR1hi CX3CR1hi KLRG1- GZMB+ TBX21+ JUNDlo FOSBlo ZBTB38hi
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="NK cells")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="NK cells")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-4, 0, 4), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_NKcells.pdf", height=2.5, width=9.5)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in OtherT cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(gdT.CMC1=c(
                                   # "TRDC","TIGIT","KLRC2","TRGC2","IKZF2","GCSAM","FCRL6","TRDV1","CST7","CMC1",
                                   "TRDV1","GZMK","GZMB","CXCR4","CX3CR1","KLRC2","KLRF1","TIGIT","TBX21","JUN","ZNF683","EOMES"), # annotated Vd1 GZMK+, Vd1 GZMB+ in Immunity
                        gdT.KLRC1=c(
                                    # "TRDC", "TRGC1", "TRGV9", "TRDV2", "KLRD1", "IL7R", "KLRC1", "DUSP2", "GNLY", "KLRG1",
                                    "TRDV2","TRDV9","IL7R","GZMK","CXCR4","KLRB1","KLRC1","TBX21","JUN"), # annotated Vd2 GZMK+, Vd2 GZMB+ in Immunity
                        gdT.LEF1=c(
                                   # "RTKN2", "TRDC", "TRGC2", "LEF1", "IKZF2", "SOX4", "ZNF331", "ARID5B", "NUCB2", "CRTAM",
                                   "CCR7","IL7R","CXCR4","TCF7","JUN","ZNF683"), # annotated gd naive in Immunity
                        MAIT=c("KLRB1", "NKG7", "GZMK", "SLC4A10", "NCR3", "CTSW", "IL7R", "KLRG1", "CEBPD", "DUSP2",
                               "KLRB1", "NKG7", "GZMK", "IL7R", "SLC4A10", "GZMA", "CXCR6", "PRSS35", "RBM24", "NCR3",
                               "TRAV1-2"),
                        NKT=c("CD3D","CD3E","CD3G","NCAM1","CD8A","CD8B","FCGR3A","KLRB1","KLRD1","NKG7",
                              "GNLY","CD7","CTSW","CCL5","CD2","ZFP36L2","CXCR4","LGALS1","CEBPB","S100A6","MKI67")
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="OtherT")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="OtherT")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=c(marker_genes,"TRAV1-2"), layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2, 0, 2), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_OtherT.pdf", height=2, width=8)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in Other cells
#####################################
###

library(ggplot2)
library(dplyr)

### Prepare the markers
experiment_markers=list(ILC.reg=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS",
                                  "IL7R","SOX4","ID3"),
                        ILC2=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS",
                               "TBX21","GATA3","RORC","AHR"),
                        plt=c("PF4","CAVIN2","GP9"),
                        HSC=c(
                              "CRHBP","AVP","MYCT1","BEX1","NPR3","CRYGD","MSRB3","CD34","NPDC1","MLLT3",
                              "SAMHD1"),
                        CLP=c("ACY3","PRSS2","C1QTNF4","SPINK2","SMIM24","NREP","CD34","DNTT","FLT3","SPNS3",
                              "DNTT","PRSS2"),
                        MEP=c("EPOR","KLF1","CSF2RB","APOC1","CNRIP1",
                              "CD44","TFRC","ITGA2B"),
                        MPP=c("EGFL7","CD34","CSF3R","SMIM24","C1QTNF4"),
                        NKP=c("IL7R","CD27","CD44","IL2RB","CD7") # Id2+IL-7Rα+CD25–α4β7–NKG2A/C/E+Bcl11b–
                        )

# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="Other cells")
marker_genes=intersect(merged_results$gene, unlist(experiment_markers)) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")

object_subset=subset(object, Annot.rough=="Other cells")
object_subset_pseudo=Seurat::AggregateExpression(object_subset, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

# plot
library(ComplexHeatmap)
plot_rough=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_height=unit(3.5, "cm"), grid_width=unit(0.2, "cm"),
                 direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2.5, 0, 2.5), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.2,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_rough=draw(plot_rough, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_Othercells.pdf", height=2, width=8.5)
draw(plot_rough)
dev.off()

#####################################



### Plot marker genes from experimental/review data at the detailed level in Progenitors and Other cells, seperately
#####################################
###

library(ggplot2)
library(dplyr)
library(ComplexHeatmap)

### Prepare the markers
experiment_markers=list(ILC.reg=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS",
                                  "IL7R","SOX4","ID3"),
                        ILC2=c("KIT","TRDC","TTLL10","LINC01229","SOX4","KLRB1","TNFRSF18","TNFRSF4","IL1R1","HPGDS",
                               "TBX21","GATA3","RORC","AHR"),
                        plt=c("PF4","CAVIN2","GP9"),
                        HSC=c(
                          "CRHBP","AVP","MYCT1","BEX1","NPR3","CRYGD","MSRB3","CD34","NPDC1","MLLT3",
                          "SAMHD1"),
                        CLP=c("ACY3","PRSS2","C1QTNF4","SPINK2","SMIM24","NREP","CD34","DNTT","FLT3","SPNS3",
                              "DNTT","PRSS2"),
                        MEP=c("EPOR","KLF1","CSF2RB","APOC1","CNRIP1",
                              "CD44","TFRC","ITGA2B"),
                        MPP=c("EGFL7","CD34","CSF3R","SMIM24","C1QTNF4"),
                        NKP=c("IL7R","CD27","CD44","IL2RB","CD7") # Id2+IL-7Rα+CD25–α4β7–NKG2A/C/E+Bcl11b–
)
# load the cell-specific genes
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t") %>%
  subset(Annot.rough=="Other cells")
marker_genes_prog=intersect(merged_results$gene, unlist(experiment_markers[4:8])) %>% .[!duplicated(.)]
marker_genes_others=intersect(merged_results$gene, unlist(experiment_markers[1:3])) %>% .[!duplicated(.)]

### Plot the marker genes
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
object=subset(object, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")
object_subset=subset(object, Annot.rough=="Other cells")
object_progenitor=subset(object_subset, Annot.inter=="Other.progenitor")
object_other.than.progenitor=subset(object_subset, Annot.inter!="Other.progenitor")

# plot the progenitor subpopulations
object_subset_pseudo=Seurat::AggregateExpression(object_progenitor, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes_prog, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))
plot_prog=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, 
          heatmap_legend_param=
            list(title="Z-score",
                 legend_width=unit(3.5, "cm"), grid_height=unit(0.2, "cm"),
                 direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2.5, 0, 2.5), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.5,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_prog=draw(plot_prog, heatmap_legend_side="top")

# plot the other subpopulations
object_subset_pseudo=Seurat::AggregateExpression(object_other.than.progenitor, group.by="Annot.detailed", return.seurat=T)
colnames(object_subset_pseudo)
expr_matrix=Seurat::FetchData(object_subset_pseudo, vars=marker_genes_others, layer="data")
rownames(expr_matrix)=gsub("-","_",rownames(expr_matrix))

library(ComplexHeatmap)
plot_other=
  Heatmap(scale(expr_matrix),
          name=" ",
          cluster_rows=F, cluster_columns=T, show_row_dend=F, show_column_dend=F, show_heatmap_legend=F,
          heatmap_legend_param=
            list(title="Z-score",
                 legend_width=unit(3.5, "cm"), grid_height=unit(0.2, "cm"),
                 direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
                 title_gp=gpar(fontface="plain", fontsize=10)),
          col=circlize::colorRamp2(c(-2.5, 0, 2.5), c("skyblue3", "white", "darkred")),
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          show_column_names=T,
          row_title_gp=gpar(fontsize=10), row_title_rot=0,
          height=unit(3.5,"mm")*nrow(expr_matrix),
          width=unit(3.5,"mm")*ncol(expr_matrix)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
plot_other=draw(plot_other, heatmap_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_Othercells.pdf", height=2, width=6.5)
draw(plot_prog)
draw(plot_other)
dev.off()

#####################################



### Find highly expressed genes for each celltype by analyzing the mean expression and the percentage of expressing cells
#####################################
###

library(Seurat)
# library(scCustomize)

pseudo_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudo_inter=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudo_detailed=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")

# for rough
rough_celltypes=unique(pseudo_rough$Annot.rough)
genelist_rough=list()
for (i in 1:length(rough_celltypes)) {
  subset_=subset(pseudo_rough, Annot.rough==rough_celltypes[i])
  rough_data=LayerData(subset_, assay="RNA", layer="counts")
  
  # filter by mean expression of genes
  Mean.expr_detailed=rowMeans(rough_data)
  filter_mean=Mean.expr_detailed>=5
  
  # filter by percentage of expressing donors
  filter_percent=rowSums(rough_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
  
  # get the genes
  genelist_rough[[i]]=rownames(rough_data)[filter_mean & filter_percent]
  print(paste0("Annot.rough celltype ",i, "/",length(rough_celltypes)," is done."))
}
names(genelist_rough)=rough_celltypes

# for inter
inter_celltypes=unique(pseudo_inter$Annot.inter)
genelist_inter=list()
for (i in 1:length(inter_celltypes)) {
  subset_=subset(pseudo_inter, Annot.inter==inter_celltypes[i])
  inter_data=LayerData(subset_, assay="RNA", layer="counts")
  
  # filter by mean expression of genes
  Mean.expr_detailed=rowMeans(inter_data)
  filter_mean=Mean.expr_detailed>=5
  
  # filter by percentage of expressing donors
  filter_percent=rowSums(inter_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
  
  # get the genes
  genelist_inter[[i]]=rownames(inter_data)[filter_mean & filter_percent]
  print(paste0("Annot.inter celltype ",i, "/",length(inter_celltypes)," is done."))
}
names(genelist_inter)=inter_celltypes

# for detailed
detailed_celltypes=unique(pseudo_detailed$Annot.detailed)
genelist_detailed=list()
for (i in 1:length(detailed_celltypes)) {
  subset_=subset(pseudo_detailed, Annot.detailed==detailed_celltypes[i])
  detailed_data=LayerData(subset_, assay="RNA", layer="counts")
  
  # filter by mean expression of genes
  Mean.expr_detailed=rowMeans(detailed_data)
  filter_mean=Mean.expr_detailed>=5
  
  # filter by percentage of expressing donors
  filter_percent=rowSums(detailed_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
  
  # get the genes
  genelist_detailed[[i]]=rownames(detailed_data)[filter_mean & filter_percent]
  print(paste0("Annot.detailed celltype ",i, "/",length(detailed_celltypes)," is done."))
}
names(genelist_detailed)=detailed_celltypes

# combine
genelist_all=list(Annot.rough=genelist_rough, Annot.inter=genelist_inter, Annot.detailed=genelist_detailed)
saveRDS(genelist_all, "~/Project_PBMCage/Tempt_RDS/HighlyExpressedGenes_All.3.CelltypeLevels.rds")

#####################################



### Jaccard similarity in highly expressed genes across celltypes
#####################################
###

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

### Load the data
genelist_all=readRDS("~/Project_PBMCage/Tempt_RDS/HighlyExpressedGenes_All.3.CelltypeLevels.rds")
highlyexpr_genes=genelist_all$Annot.rough

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Calculate the jaccard similarity for pos genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(highlyexpr_genes)) {
  the_first_celltype=highlyexpr_genes[[i]]
  
  for (j in 1:length(highlyexpr_genes)) {
    the_second_celltype=highlyexpr_genes[[j]]
    
    jaccard_similarity=jaccard(the_first_celltype, the_second_celltype)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(highlyexpr_genes)[i],"_vs_",names(highlyexpr_genes)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name
Jaccard_SIM_pos=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)

### Merge the pos and neg results
Jaccard_SIM_pos=Jaccard_SIM_pos %>% 
  # subset(!grepl("Eryth|Mast")) %>%
  mutate(celltype_1=gsub("_vs_.*","",pair), celltype_2=gsub(".*_vs_","",pair))

### Plot the results
Jaccard_SIM_df=Jaccard_SIM_pos %>%
  select(-pair) %>%
  tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype_1")
range(Jaccard_SIM_df[Jaccard_SIM_df!=1]) # check

col_fun=circlize::colorRamp2(c(0, 0.5, 1), c("white", "brown3", "#561515"))
ht2=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=T, cluster_rows=T, show_row_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=F,
          column_names_side="top",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(Jaccard_SIM_df), 
          width=unit(5, "mm")*ncol(Jaccard_SIM_df))

lgd2=Legend(col_fun=col_fun, 
            title="Jaccard similarity",
            legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_final_merged=draw(ht2, annotation_legend_list=lgd2, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/HighlyExprGenes_JaccardSimilarty_amongCelltypes_rough.pdf", height=3, width=2.5)
draw(ht_final_merged)
dev.off()

#####################################



### Analyze the similarity between celltypes based on gene expr PCA
#####################################
###

library(ComplexHeatmap)
library(dplyr)

PCA_Result=DistResult=list()
### dist analysis on detailed celltypes
# pseudo_detailed=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# pseudo_detailed=pseudo_detailed %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA()
# pca_data=Embeddings(pseudo_detailed, reduction="pca")
# PCA_Result[["detailed"]]=pca_data

pca_data=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")$detailed
pca_data=pca_data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>%
  mutate(Annot.inter=gsub("_.*","",gsub("-.*","",cell_id))) %>%
  select(-cell_id) %>%
  group_by(across(Annot.inter)) %>%
  summarize_at(vars(colnames(.)[grepl("PC_",colnames(.))]), mean) %>%
  tibble::column_to_rownames("Annot.inter") %>%
  .[,1:50] %>%
  .[!grepl("Other\\.Eryth|Other\\.Mast",rownames(.)),]
cell_dist=dist(pca_data) %>% as.matrix(.)
DistResult[["detailed"]]=cell_dist

### dist analysis on inter celltypes
# pseudo_inter=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
# pseudo_inter=pseudo_inter %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA()
# pca_data=Embeddings(pseudo_inter, reduction="pca")
# PCA_Result[["inter"]]=pca_data
pca_data=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")$inter
pca_data=pca_data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>%
  mutate(Annot.inter=gsub("_.*","",gsub("-.*","",cell_id))) %>%
  select(-cell_id) %>%
  group_by(across(Annot.inter)) %>%
  summarize_at(vars(colnames(.)[grepl("PC_",colnames(.))]), mean) %>%
  tibble::column_to_rownames("Annot.inter") %>%
  .[,1:50] %>%
  .[!grepl("Other\\.Eryth|Other\\.Mast",rownames(.)),]
cell_dist=dist(pca_data) %>% as.matrix(.)
DistResult[["inter"]]=cell_dist

### dist analysis on detailed celltypes
# pseudo_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
# pseudo_rough=pseudo_rough %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA()
# pca_data=Embeddings(pseudo_rough, reduction="pca")
# PCA_Result[["rough"]]=pca_data
pca_data=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")$rough
pca_data=pca_data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>%
  mutate(Annot.inter=gsub("_.*","",gsub("-.*","",cell_id))) %>%
  select(-cell_id) %>%
  group_by(across(Annot.inter)) %>%
  summarize_at(vars(colnames(.)[grepl("PC_",colnames(.))]), mean) %>%
  tibble::column_to_rownames("Annot.inter") %>%
  .[,1:50] %>%
  .[!grepl("Other cells",rownames(.)),]
cell_dist=dist(pca_data) %>% as.matrix(.)
DistResult[["rough"]]=cell_dist

### Save
# saveRDS(PCA_Result, "~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")

### Plot the distance between cell types
mat_=DistResult[["rough"]]
col_fun=circlize::colorRamp2(c(0,60), c("brown3","white"))
ht_=
  Heatmap(mat_,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=T, cluster_rows=T,
          row_names_gp=gpar(fontsize=10),
          column_names_side="top", show_column_names=T, show_row_names=F,
          column_dend_height=unit(1,"cm"), row_dend_width=unit(1,"cm"),
          row_dend_side="right",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=nrow(mat_)*unit(4.5, "mm"), 
          width=ncol(mat_)*unit(4.5, "mm"))
lgd_=Legend(col_fun=col_fun, 
            title="Distance",
            legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_final_merged=draw(ht_, annotation_legend_list=lgd_, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_PseudobulkObj_PCA_distance_Inter.pdf", height=3, width=2.5)
ht_final_merged
dev.off()

#####################################



### Analyze the PCA including ages for each celltype at the inter level
#####################################
###

## dist analysis on inter celltypes
pseudo_inter=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudo_inter_convert=Seurat::FindVariableFeatures(pseudo_inter)
var_gene=Seurat::VariableFeatures(pseudo_inter_convert)[1:100]

inter_celltypes=unique(pseudo_inter_convert$Annot.inter)
df_all=list()
for (i in 1:length(inter_celltypes)) {
  pseudo_inter_subset=pseudo_inter_convert %>% subset(Annot.inter==inter_celltypes[i])
  data_normalized=Seurat::FetchData(pseudo_inter_subset, vars=c("age",var_gene), layer="data")
  data_normalized=data_normalized[,colSums(data_normalized)!=0]
  data.pca=stats::princomp(scale(data_normalized))
  plot_=factoextra::fviz_cos2(data.pca, choice="var", axes=1:2)
  df_all[[i]]=plot_$data %>% subset(name=="age") %>% mutate(Annot.inter=inter_celltypes[i])
  print(inter_celltypes[i])
}
df_all=data.table::rbindlist(df_all)

ggplot(df_all, aes(x=Annot.inter, y=cos2)) +
  geom_point()


pseudo_inter_convert=pseudo_inter %>% subset(Annot.inter=="CD8T.naive")
pseudo_inter_convert=Seurat::FindVariableFeatures(pseudo_inter_convert)
var_gene=Seurat::VariableFeatures(pseudo_inter_convert)[1:100]
data_normalized=Seurat::FetchData(pseudo_inter_convert, vars=c("age",var_gene), layer="data")
data.pca=stats::princomp(scale(data_normalized))
summary(data.pca)
factoextra::fviz_cos2(data.pca, choice = "var", axes = 1:5)
pseudo_inter_convert=pseudo_inter %>% subset(Annot.inter=="B.naive")
pseudo_inter_convert=Seurat::FindVariableFeatures(pseudo_inter_convert)
var_gene=Seurat::VariableFeatures(pseudo_inter_convert)[1:500]
data_normalized=Seurat::FetchData(pseudo_inter_convert, vars=c("age",var_gene), layer="data")
data.pca=stats::princomp(scale(data_normalized))
summary(data.pca)
factoextra::fviz_cos2(data.pca, choice = "var", axes = 1:5) +
  scale_x_discrete(labels=function(x) ifelse(x=="age","age",""))




pseudo_inter=pseudo_inter %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
pca_data=Embeddings(pseudo_inter, reduction="pca")

pca_data=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")$inter
pca_data=pca_data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>%
  mutate(Annot.inter=gsub("_.*","",gsub("-.*","",cell_id))) %>%
  select(-cell_id) %>%
  group_by(across(Annot.inter)) %>%
  summarize_at(vars(colnames(.)[grepl("PC_",colnames(.))]), mean) %>%
  tibble::column_to_rownames("Annot.inter") %>%
  .[,1:50] %>%
  .[!grepl("Other\\.Eryth|Other\\.Mast",rownames(.)),]
cell_dist=dist(pca_data) %>% as.matrix(.)
DistResult[["inter"]]=cell_dist

### Plot the distance between cell types
pca_data=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_PseudobulkObj_PCAextraction_ALL3CelltypeLevels.rds")$inter
pca_data_arranged=pca_data %>% as.data.frame() %>% tibble::rownames_to_column("cell_id") %>%
  mutate(Annot.inter=gsub("_.*","",gsub("-.*","",cell_id))) %>%
  mutate(age=as.integer(gsub(".*-[0-9]{1,4}_","",cell_id))) %>%
  select(-cell_id) %>%
  split(.$Annot.inter) %>%
  lapply(., function(df) df %>% group_by(age) %>%
           summarize_at(colnames(.)[grepl("PC_",colnames(.))], function(x) mean(x, na.rm=T)) %>%
           tibble::column_to_rownames("age") %>%
           apply(., 2, function(x) sd(x, na.rm=T))
         )
df_=pca_data_arranged %>% as.data.frame() %>%
  tibble::rownames_to_column("PC") %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("PC", colnames(.))], names_to="Annot.inter", values_to="SD_PC") %>%
  mutate(marker=ifelse(Annot.inter %in% c("B.naive","CD4T.naive","CD8T.naive"), "marked", NA))

df_$PC=forcats::fct_relevel(df_$PC, unique(df_$PC))
ggplot(df_, aes(x=PC, y=SD_PC)) +
  geom_line(data=. %>% subset(!is.na(marker)), aes(group=Annot.inter), linewidth=0.5, color="brown3", linetype="solid") +
  geom_line(data=. %>% subset(is.na(marker)), aes(group=Annot.inter), linewidth=0.25, color="grey50", linetype="dashed") +
  ggrepel::geom_text_repel(data=. %>% subset(!is.na(marker) & PC=="PC_6"), aes(label=Annot.inter),
                           min.segment.length=0, box.padding=0.3) +
  theme_light() +
  labs(y="Standard deviation of PC", x="PC") +
  theme(plot.background=element_rect(fill="white", color="transparent"),
        title=element_text(size=9),
        strip.text=element_text(size=10),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_x_discrete(labels=function(x) gsub("PC_","",x))

ggplot(pca_data_arranged %>% subset(grepl("B\\.",Annot.inter)), 
       aes(x=PC_1, y=PC_2, color=age, shape=Annot.inter)) +
  geom_point()

range(df_, na.rm=T)
col_fun=circlize::colorRamp2(c(-1,1), c("dodgerblue3","brown3"))
mat_=pca_data_arranged
# ht_=
  Heatmap(t(scale(t(df_))),
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10),
          column_names_side="top", show_column_names=T, show_row_names=T,
          column_dend_height=unit(1,"cm"), row_dend_width=unit(1,"cm"),
          row_dend_side="right",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=nrow(df_)*unit(4.5, "mm"), 
          width=ncol(df_)*unit(4.5, "mm"))
lgd_=Legend(col_fun=col_fun, 
            title="Distance",
            legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_final_merged=draw(ht_, annotation_legend_list=lgd_, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_PseudobulkObj_PCA_distance_Inter.pdf", height=3, width=2.5)
ht_final_merged
dev.off()


#####################################










