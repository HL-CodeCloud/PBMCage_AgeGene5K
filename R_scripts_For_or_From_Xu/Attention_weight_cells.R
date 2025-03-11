
### 1-way Attention Weights
#####################################
###

object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Map the weight to the rough celltypes
# add the layer5 info
layer5=read.table("~/Project_PBMCage/For_or_From_Xu/layer5_weights.txt")
layer5_matchTo_object=layer5[match(rownames(object[[]]), layer5$V1), ]
object$layer5=layer5_matchTo_object$V2
object_CleanByLayer5=na.omit(object[[]])
pdf("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_rough.pdf", height=5, width=8)
ggplot(object_CleanByLayer5, aes(x=Annot.rough, y=layer5)) +
  geom_point(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention Weight") +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        legend.position="none")
dev.off()
# add the layer6 info
layer6=read.table("~/Project_PBMCage/For_or_From_Xu/layer6_weights.txt")
layer6_matchTo_object=layer6[match(rownames(object[[]]), layer6$V1), ]
object$layer5=layer6_matchTo_object$V2
object_CleanByLayer6=na.omit(object[[]])
ggplot(object_CleanByLayer6, aes(x=Annot.rough, y=layer5)) +
  geom_point(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention Weight") +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        legend.position="none")

### Map the weight to the inter and detailed celltypes
# layer5
ggplot(object_CleanByLayer5, aes(x=Annot.inter, y=layer5)) +
  geom_point(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention Weight") +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        legend.position="none")
# layer6
layer6_inter=
  ggplot(object_CleanByLayer6 %>% subset(!grepl("Other\\.Eryth|Other\\.Mast",Annot.inter)), aes(x=Annot.inter, y=layer5)) +
  ggrastr::geom_point_rast(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention weight") +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        legend.position="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_celltype_1way_inter.pdf", height=3.5, width=6.5)
plot(layer6_inter)
dev.off()

pdf("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_detailed.pdf", height=5, width=20)
# layer5
ggplot(object_CleanByLayer5, aes(x=Annot.detailed, y=layer5)) +
  geom_point(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention Weight") +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        legend.position="none")
# layer6
ggplot(object_CleanByLayer6, aes(x=Annot.detailed, y=layer5)) +
  geom_point(size=1, shape=1, color="gray50") +
  theme_classic() +
  labs(x=NULL, y="Attention Weight") +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        legend.position="none")
dev.off()

#####################################



### 2-way Attention Weights
#####################################
###

### Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Prepare the correlation matrix from each donor
# function to merge rows and columns with the same name
merge_cor_matrix=function(cor_matrix, celltype_level, rownames_to_order) {
  
  original_colnames=gsub("\\.","-",colnames(cor_matrix))
  celltype_colnames=THEObj[[]][match(original_colnames, colnames(THEObj)), celltype_level]
  
  merged_matrix=matrix(0, nrow=length(rownames_to_order), ncol=dim(cor_matrix)[2])
  for (i in 1:length(rownames_to_order)) {
    if (rownames_to_order[i] %in% celltype_colnames) {
      idx_to_merge=which(celltype_colnames==rownames_to_order[i])
      if (!(length(idx_to_merge) %in% c(0,1))) {
        merged_matrix[i, ]=colMeans(cor_matrix[idx_to_merge, ], na.rm=T)
      } else if (length(idx_to_merge)==1) {
        merged_matrix[i, ]=t(cor_matrix[idx_to_merge, ])
      }
    }
  }
  merged_matrix_full=matrix(0, nrow=length(rownames_to_order), ncol=length(rownames_to_order))
  for (j in 1:length(rownames_to_order)) {
    if (rownames_to_order[j] %in% celltype_colnames) {
      idx_to_merge=which(celltype_colnames==rownames_to_order[j])
      if (!(length(idx_to_merge) %in% c(0,1))) {
        merged_matrix_full[, j]=rowMeans(merged_matrix[, idx_to_merge], na.rm=T)
      } else if (length(idx_to_merge)==1) {
        merged_matrix_full[, j]=merged_matrix[, idx_to_merge]
      }
    }
  }
  
  rownames(merged_matrix_full)=colnames(merged_matrix_full)=rownames_to_order
  
  return(merged_matrix_full)
}

# get the files
mydir="~/Project_PBMCage/For_or_From_Xu/AttentionWeightMtx/"
myfiles=list.files(mydir)

# merge according to the detailed/inter/rough celltypes
all_celltype_detailed=unique(THEObj$Annot.detailed)
all_celltype_detailed=sort(all_celltype_detailed)
all_celltype_inter=unique(THEObj$Annot.inter)
all_celltype_inter=sort(all_celltype_inter)
all_celltype_rough=unique(THEObj$Annot.rough)
all_celltype_rough=sort(all_celltype_rough)

myx_age=Detailed_cor_matrix=Inter_cor_matrix=Rough_cor_matrix=list()
for (i in 1:length(myfiles)) {
  
  message(paste0("------ The No.",i, " file starts... ------"))
  
  myx=read.csv(file=paste0(mydir,myfiles[i]), header=T, sep=",")
  myx_age[[i]]=THEObj$age[match(gsub("\\.","-",colnames(myx)), colnames(THEObj))]; myx_age[[i]]=unique(na.omit(myx_age[[i]]))
  
  merged_cor_matrix_d=merge_cor_matrix(cor_matrix=myx, celltype_level="Annot.detailed", rownames_to_order=all_celltype_detailed)
  Detailed_cor_matrix[[i]]=merged_cor_matrix_d
  merged_cor_matrix_i=merge_cor_matrix(cor_matrix=myx, celltype_level="Annot.inter", rownames_to_order=all_celltype_inter)
  Inter_cor_matrix[[i]]=merged_cor_matrix_i
  merged_cor_matrix_r=merge_cor_matrix(cor_matrix=myx, celltype_level="Annot.rough", rownames_to_order=all_celltype_rough)
  Rough_cor_matrix[[i]]=merged_cor_matrix_r
  
  message(paste0("====== The No.",i, " file Complete! ======"))
}

names(Detailed_cor_matrix)=names(Inter_cor_matrix)=names(Rough_cor_matrix)=unlist(myx_age)
Detailed_cor_matrix=Detailed_cor_matrix[order(names(Detailed_cor_matrix))]
Inter_cor_matrix=Inter_cor_matrix[order(names(Inter_cor_matrix))]
Rough_cor_matrix=Rough_cor_matrix[order(names(Rough_cor_matrix))]

saveRDS(Detailed_cor_matrix, "~/Project_PBMCage/For_or_From_Xu/CorMatrix_detailed.rds")
saveRDS(Inter_cor_matrix, "~/Project_PBMCage/For_or_From_Xu/CorMatrix_intermediate.rds")
saveRDS(Rough_cor_matrix, "~/Project_PBMCage/For_or_From_Xu/CorMatrix_rough.rds")

### Plot the heatmap
library(ComplexHeatmap)
# load
Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_detailed.rds")
Inter_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_intermediate.rds")
Rough_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_rough.rds")

# plot detailed_cor_matrix
MATRIX=Reduce(`+`, Detailed_cor_matrix)/length(Detailed_cor_matrix)
plot_detailed=
  Heatmap(
    MATRIX, 
    col=rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 30)),
    name="Attention\nWeight",
    heatmap_legend_param=list(title_gp=gpar(fontsize=10, fontface="bold")),
    cluster_columns=F, 
    cluster_rows=F, 
    row_names_side="left", 
    column_names_side="top",
    row_names_gp=gpar(fontsize=5.5),
    column_names_gp=gpar(fontsize=5.5),
    height=nrow(MATRIX)*unit(2.5,"mm"),
    width=ncol(MATRIX)*unit(2.5,"mm")
  )

# plot inter_cor_matrix
MATRIX=Reduce(`+`, Inter_cor_matrix)/length(Inter_cor_matrix)
col_fun=circlize::colorRamp2(c(0,0.004), c("white", "#A50F15"))
col_celltype=c(paletteer::paletteer_d("ggsci::category20_d3"), paletteer::paletteer_d("ggsci::category20b_d3"))
names(col_celltype)=colnames(MATRIX)
plot_inter=
  Heatmap(
    MATRIX, 
    col=col_fun,
    name=" ", show_heatmap_legend=FALSE,
    cluster_columns=F, 
    cluster_rows=F, 
    show_row_names=F,
    show_column_names=F,
    row_names_gp=gpar(fontsize=10),
    column_names_gp=gpar(fontsize=10),
    height=nrow(MATRIX)*unit(2,"mm"),
    width=ncol(MATRIX)*unit(2.5,"mm"),
    top_annotation=HeatmapAnnotation(
      `cell type`=colnames(MATRIX),
      simple_anno_size=unit(0.2, "cm"), show_annotation_name=FALSE,
      col=list(`cell type`=col_celltype[!is.na(names(col_celltype))]),
      show_legend=c(`cell type`=FALSE)),
    left_annotation=rowAnnotation(
      `cell type`=colnames(MATRIX),
      simple_anno_size=unit(0.2, "cm"), show_annotation_name=FALSE,
      col=list(`cell type`=col_celltype[!is.na(names(col_celltype))]),
      show_legend=c(`cell type`=FALSE))
  )
lgd=Legend(col_fun=col_fun, 
           title="Attention weight",
           legend_width=unit(5, "cm"), grid_height=unit(0.3, "cm"),
           direction="horizontal", title_position="lefttop", title_gp=gpar(fontsize=11, fontface="bold"), labels_gp=gpar(fontsize=9),
)
topanno_lgd=Legend(labels=colnames(MATRIX), 
                   legend_gp = gpar(fill=col_celltype[!is.na(names(col_celltype))]),
                   title = "cell type", ncol = 3,
                   title_gp=gpar(fontsize=10, fontface="bold")
                   )
plot_inter_final=draw(plot_inter, annotation_legend_list=list(lgd, topanno_lgd), annotation_legend_side="right")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_2way_eachLevels.pdf", height=3.5, width=7)
plot(plot_inter_final)
dev.off()

# plot rough_cor_matrix
MATRIX=Reduce(`+`, Rough_cor_matrix)/length(Rough_cor_matrix)
plot_rough=
  Heatmap(
    MATRIX, 
    col=rev(paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 30)),
    name="Attention\nWeight",
    heatmap_legend_param=list(title_gp=gpar(fontsize=10, fontface="bold")),
    cluster_columns=F, 
    cluster_rows=F, 
    row_names_side="left", 
    column_names_side="top",
    row_names_gp=gpar(fontsize=5.5),
    column_names_gp=gpar(fontsize=5.5),
    height=nrow(MATRIX)*unit(2.5,"mm"),
    width=ncol(MATRIX)*unit(2.5,"mm")
  )

pdf("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_2way_eachLevels.pdf", height=12, width=12)
draw(plot_detailed)
draw(plot_inter)
draw(plot_rough)
dev.off()

#####################################



### Cluster attention matrix of celltype-by-celltype
#####################################
###

library("factoextra")
library("FactoMineR")

### For rough level
# arrange the matrix
Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_rough.rds")
MATRIX=Reduce(`+`, Detailed_cor_matrix)/length(Detailed_cor_matrix)
MATRIX_wo_zerosd=MATRIX[which(apply(MATRIX, 2, var)!=0), which(apply(MATRIX, 2, var)!=0)]
# run pca and clustering
pc_cell=prcomp(MATRIX_wo_zerosd,
               center=TRUE,
               scale.=TRUE)
results=pc_cell$x
fviz_nbclust(results, FUNcluster=kmeans, k.max=7) # check the opitmal number of clusters
km=eclust(results, "kmeans", hc_metric="eucliden", k=3)
# plot
km_data=km$clust_plot$data
plot_rough=
  ggplot(km_data, aes(x=x, y=y, color=cluster, shape=cluster, fill=cluster)) +
  geom_point() +
  geom_polygon(data=. %>% group_by(cluster) %>% do(.[chull(.[2:3]), ]), alpha=0.25, linewidth=0.1) +
  ggrepel::geom_text_repel(aes(label=name), box.padding=0.5, min.segment.length=2) +
  labs(x=km$clust_plot$labels$x, y=km$clust_plot$labels$y, title="Cluster plot - level1") +
  scale_shape_manual(values=c(1,0,2)) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::few_Light")) +
  scale_fill_manual(values=paletteer::paletteer_d("ggthemes::few_Light")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_celltype_2way_cluster_rough.pdf", height=3.5, width=3.5)
plot_rough
dev.off()

### For inter level
# arrange the matrix
Inter_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_intermediate.rds")
MATRIX=Reduce(`+`, Inter_cor_matrix)/length(Inter_cor_matrix)
rownames(MATRIX)=colnames(MATRIX)=paste0(rownames(MATRIX),rep(" ",5)) # to avoid squeezing the labels when plotting
MATRIX_wo_zerosd=MATRIX[which(apply(MATRIX, 2, var)!=0), which(apply(MATRIX, 2, var)!=0)]
# remove DCs, Other cells, Monocytes
MATRIX_wo_zerosd=MATRIX_wo_zerosd[!grepl("^DC\\.|^Other\\.|^Mono\\.",rownames(MATRIX_wo_zerosd)),]
# run pca and clustering
res.pca=PCA(MATRIX_wo_zerosd, ncp=10, graph=FALSE)
res.hcpc=HCPC(res.pca, graph=FALSE, nb.clust=6)
plot_dend=
  fviz_dend(res.hcpc, 
            lwd=0.25,
            cex=0.8,
            palette="jco",
            rect=TRUE, rect_fill=TRUE,
            rect_border="jco",
            labels_track_height=15,
            main="Cluster dendrogram - level2", 
            xlab="", ylab="Height", 
            ggtheme=theme_classic()
  ) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_celltype_2way_cluster_inter.pdf", height=3.5, width=7.25)
plot_dend
dev.off()
    









### Compare the 2-way celltype-celltype attention matrix from each donor
#######################################
###

library(EFAtools)
library(factoextra)

### Analyze the detailed celltype-celltype matrices
Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_detailed.rds")
DF=list()
for (i in 1:length(Detailed_cor_matrix)) {
  mtx1=Detailed_cor_matrix[[i]]
  for (j in 1:length(Detailed_cor_matrix)) {
    mtx2=Detailed_cor_matrix[[j]]
    tempt=COMPARE(mtx1, mtx2, reorder="none")
    df_=data.frame(mean.abs.diff=tempt$mean_abs_diff,
                   median.abs.diff=tempt$median_abs_diff,
                   min.abs.diff=tempt$min_abs_diff,
                   max.abs.diff=tempt$max_abs_diff,
                   group=paste0(names(Detailed_cor_matrix)[i],"=",i, "_vs_", names(Detailed_cor_matrix)[j],"=",j)
    )
    DF=c(DF, list(df_))
  }
  print(paste0("Done with : ",i))
}
DF_merged=data.table::rbindlist(DF) %>% mutate(compare.element="celltype.detailed")
data.table::fwrite(DF_merged, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_celltype.detailed_2way_CompareMatrixResults.txt.gz", sep="\t")

### Analyze the intermediate celltype-celltype matrices
Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_intermediate.rds")
DF=list()
for (i in 1:length(Detailed_cor_matrix)) {
  mtx1=Detailed_cor_matrix[[i]]
  for (j in 1:length(Detailed_cor_matrix)) {
    mtx2=Detailed_cor_matrix[[j]]
    tempt=COMPARE(mtx1, mtx2, reorder="none")
    df_=data.frame(mean.abs.diff=tempt$mean_abs_diff,
                   median.abs.diff=tempt$median_abs_diff,
                   min.abs.diff=tempt$min_abs_diff,
                   max.abs.diff=tempt$max_abs_diff,
                   group=paste0(names(Detailed_cor_matrix)[i],"=",i, "_vs_", names(Detailed_cor_matrix)[j],"=",j)
    )
    DF=c(DF, list(df_))
  }
  print(paste0("Done with : ",i))
}
DF_merged=data.table::rbindlist(DF) %>% mutate(compare.element="celltype.inter")
data.table::fwrite(DF_merged, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_celltype.inter_2way_CompareMatrixResults.txt.gz", sep="\t")

### Analyze the rough celltype-celltype matrices
Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_rough.rds")
DF=list()
for (i in 1:length(Detailed_cor_matrix)) {
  mtx1=Detailed_cor_matrix[[i]]
  for (j in 1:length(Detailed_cor_matrix)) {
    mtx2=Detailed_cor_matrix[[j]]
    tempt=COMPARE(mtx1, mtx2, reorder="none")
    df_=data.frame(mean.abs.diff=tempt$mean_abs_diff,
                   median.abs.diff=tempt$median_abs_diff,
                   min.abs.diff=tempt$min_abs_diff,
                   max.abs.diff=tempt$max_abs_diff,
                   group=paste0(names(Detailed_cor_matrix)[i],"=",i, "_vs_", names(Detailed_cor_matrix)[j],"=",j)
    )
    DF=c(DF, list(df_))
  }
  print(paste0("Done with : ",i))
}
DF_merged=data.table::rbindlist(DF) %>% mutate(compare.element="celltype.rough")
data.table::fwrite(DF_merged, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_celltype.rough_2way_CompareMatrixResults.txt.gz", sep="\t")


### Rearrange the results
DF_merged_detailed=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_celltype.inter_2way_CompareMatrixResults.txt.gz", sep="\t")
DF_merged_detailed=DF_merged_detailed %>%
  dplyr::select(mean.abs.diff, group) %>%
  mutate(group1=gsub("_vs_.*","",group), group2=gsub(".*_vs_","",group)) %>%
  dplyr::select(-group) %>%
  tidyr::pivot_wider(names_from="group2", values_from="mean.abs.diff") %>%
  tibble::column_to_rownames("group1")
# cluster based on the dist
fviz_nbclust(DF_merged_detailed, FUNcluster=kmeans, k.max=8) # check the optimal number of clusters, be 2 means 1 community
DF_merged_detailed_dist=as.dist(DF_merged_detailed)
hc.cut=hcut(DF_merged_detailed_dist, k=2, hc_method="ward.D2")
fviz_cluster(hc.cut, data=DF_merged_detailed_dist)
fviz_dend(hc.cut, FUNcluster=kmeans, k.max=8)


