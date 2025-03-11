# 
# ### Draw DimPlot
# #####################################
# ###
# 
# ###
# library(Seurat)
# object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# 
# pdf("~/Project_PBMCage/Plots/DimPlot.pdf", height=7, width=14)
# DimPlot(object, group.by="Annot.detailed", raster=TRUE) +
#   theme_void() +
#   theme(axis.text.x=element_blank(), 
#         axis.ticks.x=element_blank(), 
#         axis.text.y=element_blank(), 
#         axis.ticks.y=element_blank()) +
#   labs(title=NULL) +
#   guides(color=guide_legend(ncol=3, bycol=TRUE))
# dev.off()
# 
# #####################################



### Marker genes
#####################################
###

###
library(Seurat)
library(dplyr)
object=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")

### Find markers for Annot.rough
Idents(object)="Annot.rough"
markers_rough=FindAllMarkers(object)
saveRDS(markers_rough, "~/Project_PBMCage/Results/Immunity_results/Cell_Markers_rough.rds")

### Find markers for each Annot.inter
markers_inter=list()
celltype.rough=unique(object$Annot.rough)
for (i in 1:length(celltype.rough)) {
  subset_=subset(object, Annot.rough==celltype.rough[i])
  Idents(subset_)="Annot.inter"
  markers_inter[[i]]=FindAllMarkers(subset_)
}
saveRDS(markers_inter, "~/Project_PBMCage/Results/Immunity_results/Cell_Markers_inter.rds")

### Find markers for each Annot.detailed
markers_detailed=list()
celltype.inter=unique(object$Annot.inter)
for (i in 1:length(celltype.inter)) {
  subset_=subset(object, Annot.inter==celltype.inter[i])
  Idents(subset_)="Annot.detailed"
  markers_detailed[[i]]=FindAllMarkers(subset_)
}
names(markers_detailed)=unique(object$Annot.inter)
saveRDS(markers_detailed, "~/Project_PBMCage/Results/Immunity_results/Cell_Markers_detailed.rds")

### Make celltype-specific gene lists
markers_rough=readRDS("~/Project_PBMCage/Results/Immunity_results/Cell_Markers_rough.rds")
markers_inter=readRDS("~/Project_PBMCage/Results/Immunity_results/Cell_Markers_inter.rds")
markers_detailed=readRDS("~/Project_PBMCage/Results/Immunity_results/Cell_Markers_detailed.rds")
# at the rough level
markers_rough_perCluster=markers_rough %>% 
  split(.$cluster) %>%
  lapply(., function(df) df %>% 
           subset(p_val_adj<0.05 & avg_log2FC>0) %>%
           select(cluster, gene, avg_log2FC, p_val_adj) %>%
           mutate(celltype.level="Annot.rough")) %>%
  data.table::rbindlist(.) %>%
  mutate(Annot.detailed=NA, Annot.inter=NA) %>%
  dplyr::rename(Annot.rough=cluster) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough")))
markers_rough_perCluster %>% split(.$Annot.rough) %>% lapply(., nrow) # check if all the rough celltypes have their marker genes

# at the inter level
markers_inter_perCluster=markers_inter %>% 
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj) %>%
               mutate(celltype.level="Annot.inter")) %>%
      data.table::rbindlist(.)
  }) %>%
  .[!unlist(lapply(., function(df) nrow(df)==0))] # remove the empty inter
lapply(lapply(markers_inter_perCluster, function(DF) DF %>% split(.$cluster)), function(DF) lapply(DF, nrow)) # check if all the inter celltypes have their marker genes
markers_inter_perCluster=lapply(1:length(markers_inter_perCluster), 
                                function(idx) markers_inter_perCluster[[idx]] %>% 
                                  mutate(Annot.rough=names(markers_inter_perCluster)[idx])) %>%
  data.table::rbindlist(.) %>%
  dplyr::rename(Annot.inter=cluster) %>%
  mutate(Annot.rough=gsub("\\..*","",Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="B","B cells",
                            ifelse(Annot.rough=="CD4T","CD4T cells",
                                   ifelse(Annot.rough=="CD8T","CD8T cells",
                                          ifelse(Annot.rough=="DC","DCs",
                                                 ifelse(Annot.rough=="Mono","Monocytes",
                                                        ifelse(Annot.rough=="Other","Other cells",
                                                               ifelse(Annot.rough=="OtherT","OtherT",
                                                                      ifelse(Annot.rough=="NK","NK cells",NA))))))))) %>%
  mutate(Annot.detailed=NA) %>%
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough")))

# at the detailed level
inter_to_remove=names(unlist(lapply(markers_detailed, nrow))[unlist(lapply(markers_detailed, nrow))==0]) # determine the skipped celltypes at the inter level
markers_detailed_perCluster=markers_detailed %>% 
  .[!(names(.) %in% inter_to_remove)] %>% # remove the skipped celltypes
  lapply(., function(DF) {
    DF %>% split(.$cluster) %>%
      lapply(., function(df) df %>% 
               subset(p_val_adj<0.05 & avg_log2FC>0) %>%
               select(cluster, gene, avg_log2FC, p_val_adj) %>%
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
  relocate(any_of(c("gene", "avg_log2FC", "p_val_adj", "celltype.level", "Annot.detailed", "Annot.inter", "Annot.rough")))

# merge at all 3 levels
merged_results=data.table::rbindlist(list(markers_rough_perCluster, markers_inter_perCluster, markers_detailed_perCluster))
data.table::fwrite(merged_results, "~/Project_PBMCage/Results/Immunity_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")


# ### Plot the marker genes for Annot.rough
# markers_rough=readRDS("~/Project_PBMCage/Results/Immunity_results/Cell_Markers_rough.rds")
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
# ### Plot the marker genes for Annot.detailed
# Markers_detailed=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_detailed.rds")
# plot_=list()
# for (i in 1:length(Markers_detailed)) {
#   markers_detailed=Markers_detailed[[i]]
#   plot_[[i]]=
#     markerVocalno(markers=markers_detailed,
#                   topn=5,
#                   labelCol=ggsci::pal_d3("category20")(20))
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/MakerGenes_averageExpr_detailed_",i,".pdf"), height=6, width=32)
#   plot(plot_[[i]])
#   dev.off()
# }
# 
# #####################################
# 
# 
# 
# ### Find highly expressed genes for each celltype by analyzing the mean expression and the percentage of expressing cells
# #####################################
# ###
# 
# library(Seurat)
# # library(scCustomize)
# 
# pseudo_rough_immunity=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
# pseudo_inter_immunity=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
# pseudo_detailed_immunity=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
# 
# # for rough
# rough_celltypes=unique(pseudo_rough_immunity$Annot.rough)
# genelist_rough=list()
# for (i in 1:length(rough_celltypes)) {
#   subset_=subset(pseudo_rough_immunity, Annot.rough==rough_celltypes[i])
#   rough_data=LayerData(subset_, assay="RNA", layer="counts")
#   
#   # filter by mean expression of genes
#   Mean.expr_detailed=rowMeans(rough_data)
#   filter_mean=Mean.expr_detailed>=5
#   
#   # filter by percentage of expressing donors
#   filter_percent=rowSums(rough_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
#   
#   # get the genes
#   genelist_rough[[i]]=rownames(rough_data)[filter_mean & filter_percent]
#   print(paste0("Annot.rough celltype ",i, "/",length(rough_celltypes)," is done."))
# }
# names(genelist_rough)=rough_celltypes
# 
# # for inter
# inter_celltypes=unique(pseudo_inter_immunity$Annot.inter)
# genelist_inter=list()
# for (i in 1:length(inter_celltypes)) {
#   subset_=subset(pseudo_inter_immunity, Annot.inter==inter_celltypes[i])
#   inter_data=LayerData(subset_, assay="RNA", layer="counts")
#   
#   # filter by mean expression of genes
#   Mean.expr_detailed=rowMeans(inter_data)
#   filter_mean=Mean.expr_detailed>=5
#   
#   # filter by percentage of expressing donors
#   filter_percent=rowSums(inter_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
#   
#   # get the genes
#   genelist_inter[[i]]=rownames(inter_data)[filter_mean & filter_percent]
#   print(paste0("Annot.inter celltype ",i, "/",length(inter_celltypes)," is done."))
# }
# names(genelist_inter)=inter_celltypes
# 
# # for detailed
# detailed_celltypes=unique(pseudo_detailed_immunity$Annot.detailed)
# genelist_detailed=list()
# for (i in 1:length(detailed_celltypes)) {
#   subset_=subset(pseudo_detailed_immunity, Annot.detailed==detailed_celltypes[i])
#   detailed_data=LayerData(subset_, assay="RNA", layer="counts")
#   
#   # filter by mean expression of genes
#   Mean.expr_detailed=rowMeans(detailed_data)
#   filter_mean=Mean.expr_detailed>=5
#   
#   # filter by percentage of expressing donors
#   filter_percent=rowSums(detailed_data!=0)/ncol(subset_)>=0.90 # more than 90% samples expressing
#   
#   # get the genes
#   genelist_detailed[[i]]=rownames(detailed_data)[filter_mean & filter_percent]
#   print(paste0("Annot.detailed celltype ",i, "/",length(detailed_celltypes)," is done."))
# }
# names(genelist_detailed)=detailed_celltypes
# 
# # combine
# genelist_all=list(Annot.rough=genelist_rough, Annot.inter=genelist_inter, Annot.detailed=genelist_detailed)
# saveRDS(genelist_all, "~/Project_PBMCage/Immunity/Tempt_RDS/HighlyExpressedGenes_All.3.CelltypeLevels.rds")
# 
# #####################################
# 
# 
# 
# 
# 
