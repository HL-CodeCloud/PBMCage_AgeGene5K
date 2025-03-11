
###############################################################################################################
##################################### Cell-by-Cell Attention Weight #####################################
###############################################################################################################

### Arrange the attention weights (2000 cells * 2000 cells * 99 donors)
#####################################
###

library(dplyr)

### Load the matrices
mydir="~/Project_PBMCage/For_or_From_Xu/AttentionMatrix_full/layer6_cell2cell.zip"
myfiles=paste0("layer6/",c(0:98),".csv")
Matrices=lapply(myfiles, function(file) readr::read_csv(unzip(mydir, file)))

### Match the donor ids
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage_donor_id.sex.age=metadata_PBMCage %>% dplyr::select(c(donor_id,sex,age,cellid))
cellid_1st=lapply(Matrices, function(mtx) colnames(mtx)[1]) %>% unlist() %>% unname()
donor_id_match=data.frame(cellid=cellid_1st) %>% left_join(., metadata_PBMCage_donor_id.sex.age, by="cellid")
names(Matrices)=donor_id_match$donor_id

### Save the matrices
saveRDS(list(MATRIX=Matrices, METADATA=list(donor_id_match)), 
        "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices.rds")

# ### Reorder the matrices by Annot.detailed
# CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices.rds")
# CELL_matrices_metadata=CELL_matrices_metadata[["MATRIX"]]
# # remove the empty columns (in the case of <2000 cells)
# CELL_matrices_metadata=lapply(CELL_matrices_metadata, function(mtx) mtx[!grepl("^0",colnames(mtx)),!grepl("^0",colnames(mtx))])
# metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
# metadata_PBMCage=metadata_PBMCage %>% dplyr::select(cellid, Annot.detailed)
# # order
# mtx_list_=list()
# for (i in 1:length(CELL_matrices_metadata)) {
#   df_=data.frame(cellid=colnames(CELL_matrices_metadata[[i]])) %>% left_join(., metadata_PBMCage, by="cellid")
#   order_=order(df_$Annot.detailed)
#   mtx_reordered=CELL_matrices_metadata[[i]][order_,order_] %>% as.matrix()
#   colnames(mtx_reordered)=rownames(mtx_reordered)=make.names(df_$Annot.detailed[order_], unique=T)
#   mtx_list_[[i]]=mtx_reordered
#   print(paste0("Done with ",i))
# }
# # expand
# r_all=Reduce(f=intersect, lapply(mtx_list_, rownames))
# c_all=Reduce(f=intersect, lapply(mtx_list_, colnames))
# element_names=outer(r_all, c_all, FUN=paste)
# mtx_expand=list()
# for (i in 1:length(mtx_list_)) {
#   indx=element_names %in% outer(rownames(mtx_list_[[i]]), colnames(mtx_list_[[i]]), FUN=paste)
#   C=matrix(NA, nrow=length(r_all), ncol=length(c_all), dimnames=list(r_all, c_all))
#   C[indx]=mtx_list_[[i]]
#   mtx_expand[[i]]=C
#   print(paste0("Done with ",i))
# }
# names(mtx_expand)=names(CELL_matrices_metadata)
# saveRDS(mtx_expand, "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices_reorderBydetailed_expanded.rds")

# ### Reorder the matrices by Annot.rough
# CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices.rds")
# CELL_matrices_metadata=CELL_matrices_metadata[["MATRIX"]]
# # remove the empty columns (in the case of <2000 cells)
# CELL_matrices_metadata=lapply(CELL_matrices_metadata, function(mtx) mtx[!grepl("^0",colnames(mtx)),!grepl("^0",colnames(mtx))])
# metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
# metadata_PBMCage=metadata_PBMCage %>% dplyr::select(cellid, Annot.rough)
# # order
# mtx_list_=list()
# for (i in 1:length(CELL_matrices_metadata)) {
#   df_=data.frame(cellid=colnames(CELL_matrices_metadata[[i]])) %>% left_join(., metadata_PBMCage, by="cellid")
#   order_=order(df_$Annot.rough)
#   mtx_reordered=CELL_matrices_metadata[[i]][order_,order_] %>% as.matrix()
#   colnames(mtx_reordered)=rownames(mtx_reordered)=make.names(df_$Annot.rough[order_], unique=T)
#   mtx_list_[[i]]=mtx_reordered
#   print(paste0("Done with ",i))
# }
# # expand
# r_all=Reduce(f=union, lapply(mtx_list_, rownames))
# c_all=Reduce(f=union, lapply(mtx_list_, colnames))
# element_names=outer(r_all, c_all, FUN=paste)
# mtx_expand=list()
# for (i in 1:length(mtx_list_)) {
#   indx=element_names %in% outer(rownames(mtx_list_[[i]]), colnames(mtx_list_[[i]]), FUN=paste)
#   C=matrix(NA, nrow=length(r_all), ncol=length(c_all), dimnames=list(r_all, c_all))
#   C[indx]=mtx_list_[[i]]
#   mtx_expand[[i]]=C
#   print(paste0("Done with ",i))
# }
# names(mtx_expand)=names(CELL_matrices_metadata)
# saveRDS(mtx_expand, "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices_reorderByrough_expanded.rds")

#####################################



### Calculate the distance of the matrices by distribution between donors
#####################################
###

library(dplyr)

### Load the cell-by-cell matrices
CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices.rds")
CELL_matrices=CELL_matrices_metadata[["MATRIX"]]
DF=list()
for (i in 1:length(CELL_matrices)) {
  mtx1=CELL_matrices[[i]] %>% as.matrix()
  dens1=density(as.vector(mtx1))
  vec1=dens1$y*diff(dens1$x)[1]/sum(dens1$y*diff(dens1$x)[1])
  for (j in 1:length(CELL_matrices)) {
    mtx2=CELL_matrices[[j]] %>% as.matrix()
    dens2=density(as.vector(mtx2))
    vec2=dens2$y*diff(dens2$x)[1]/sum(dens2$y*diff(dens2$x)[1])
    res=philentropy::KL(rbind(vec1, vec2), unit="log")
    df_=data.frame(kl_divergence=unname(res),
                   group=paste0(names(CELL_matrices)[i],"_vs_", names(CELL_matrices)[j])
    )
    DF=c(DF, list(df_))
    print(paste0("Done with j=",j))
  }
  print(paste0("====== Done with i=",i,". ======"))
}
DF_merged=data.table::rbindlist(DF) %>% mutate(compare.element="cellid")
data.table::fwrite(DF_merged, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_CellByCell_CompareMatrix_acrossDonor_KLdivergence.txt.gz", sep="\t")

#####################################



### Arrange the matrix distance between donors and do hierarchical clustering
#####################################
###

### Rearrange the distance
DF_merged=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_CellByCell_CompareMatrix_acrossDonor_KLdivergence.txt.gz", sep="\t")
DF_merged=DF_merged %>%
  dplyr::select(kl_divergence, group) %>%
  mutate(group1=gsub("_vs_.*","",group), group2=gsub(".*_vs_","",group)) %>%
  dplyr::select(-group) %>%
  tidyr::pivot_wider(names_from="group2", values_from="kl_divergence") %>%
  tibble::column_to_rownames("group1")

### Plot distance
dist_m=as.dist(DF_merged)
dist_mi=1/dist_m # one over, as qgraph takes similarity matrices as input
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_CompareMatrixAcrossDonor_distance.pdf", height=3.5, width=3.5)
qgraph::qgraph(dist_mi, layout='spring', vsize=3)
dev.off()

### Determine the number of clusters
# ... turns out to be n=2
NbClust::NbClust(data=NULL, diss=dist_m, distance=NULL, method="ward.D2", index="silhouette")

### Perform clustering
hc.cut=factoextra::hcut(dist_m, k=2, hc_method="ward.D2")

### Plot the cluster dendrogram
# arrange the results
dend=as.dendrogram(hc.cut)
dend_data=ggdendro::dendro_data(dend, type="rectangle")
# prepare the meta info
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage_donor_id.sex.age=metadata_PBMCage %>% dplyr::select(c(donor_id,sex,age)) %>% subset(!duplicated(.))

plot_dend=
  ggplot(dend_data$segments) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.25)+
  geom_text(data=dend_data$labels %>% 
              left_join(., metadata_PBMCage_donor_id.sex.age %>% dplyr::rename(label=donor_id)) %>%
              mutate(age_sex=ifelse(sex=="female",-age,age)) %>%
              mutate(age_sex_name=paste0(ifelse(sex=="female","F","M"),"  ",age,paste0(rep(" ",2), collapse=" "))), 
            aes(x, y, label=age_sex_name, color=age_sex),
            hjust=1, angle=90, size=3) +
  coord_cartesian(ylim=c(-1.5,1.3e1)) +
  labs(x=NULL, y="Height", title="Cluster dendrogram") +
  scale_y_continuous(labels=function(x) format(x, scientific=FALSE)) +
  scale_color_gradient2(low="brown3", mid="white", high="dodgerblue3") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none",
        legend.spacing.x=unit(0.1,"cm"),
        legend.direction="vertical") +
  guides(color=guide_legend(labels=function(x) abs(x)))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_CompareMatrixAcrossDonor_dend.pdf", height=3.5, width=12)
plot_dend
dev.off()

### Check the distribution of sex & age in the clusters
# ... turns out to be no sig. diff. between clusters
df_=data.frame(donor_id=names(hc.cut$cluster), cluster=hc.cut$cluster) %>%
  left_join(., metadata_PBMCage_donor_id.sex.age, by="donor_id")
# plot the ages in each cluster
wilcox_res=
  wilcox.test((df_ %>% subset(cluster=="1"))$age, (df_ %>% subset(cluster=="2"))$age)
plot_age_dens=
  ggplot(df_, aes(x=age)) +
  geom_density(aes(group=cluster, fill=as.character(cluster)), alpha=0.5, linewidth=0.25) +
  theme_classic() +
  geom_vline(xintercept=df_ %>% group_by(cluster) %>% summarize_at("age", function(x) exp(mean(log(x)))) %>% .[["age"]],
             color=c("lightgreen","darkgreen"),
             linetype="dashed") +
  ggrepel::geom_text_repel(data=df_ %>% group_by(cluster) %>% summarize_at("age", function(x) exp(mean(log(x)))),
                           aes(x=age, y=0.08, label=floor(age)),
                           color="black", box.padding=0.5, point.padding=0.5, min.segment.length=0, size=3) +
  scale_fill_manual(values=c("lightgreen","darkgreen")) +
  theme_classic() +
  labs(x="age",y="density") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none") +
  annotate(geom="text", x=38, y=0.045, label=paste0("p=",signif(wilcox_res$p.value, 2)), size=3)

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_CompareMatrixAcrossDonor_clustering_AgeDistribution.pdf", height=2.5, width=1.5)
plot_age_dens
dev.off()

# plot the ages in each cluster
z_test_res=
  prop.test(c(nrow(df_ %>% subset(cluster=="1" & sex=="female")), nrow(df_ %>% subset(cluster=="2" & sex=="female"))), 
            n=c(nrow(df_ %>% subset(cluster=="1")), nrow(df_ %>% subset(cluster=="2"))),
            correct=F)
plot_sex_dens=
  ggplot(df_ %>% mutate(sex=ifelse(sex=="female","F","M")), aes(x=sex)) +
  geom_bar(aes(group=cluster, fill=as.character(cluster)), alpha=0.5, linewidth=0.25) +
  theme_classic() +
  scale_fill_manual(values=c("lightgreen","darkgreen")) +
  theme_classic() +
  labs(x="sex",y="count") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="right",
        legend.key.size=unit(3,"mm"),
        # legend.spacing.y=unit(2,"mm"),
        legend.margin=margin(l=3),
        legend.box.spacing=margin(3)
        ) +
  annotate(geom="text", x="F", y=70, hjust=0, label=paste0("p=",round(z_test_res$p.value,2)), size=3) +
  guides(fill=guide_legend(title="cluster", override.aes=list(size=3)))
    
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_CompareMatrixAcrossDonor_clustering_SexDistribution.pdf", height=2.5, width=1.55)
plot_sex_dens
dev.off()

#####################################



# ### DON'T DO THIS!!!: Take the mean of matrices from all donors, cell-by-cell instead of celltype-by-celltype
# # ... BECAUSE EACH MATRIX HAS DIFFERENT FEATURES (CELLS); THEY CANNOT BE POOLED
# #####################################
# ###
# 
# ### For matrices sorted and expanded by Annot.rough (ie., colnames=c("B cells","B cells.1","CD4T cells","CD4T cells.1","CD4T cells.2",...))
# CELL_matrices=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices_reorderByrough_expanded.rds")
# # check
# unique(lapply(CELL_matrices, nrow) %>% unlist() %>% unname()==lapply(CELL_matrices, ncol) %>% unlist() %>% unname())
# CELL_matrices_mean=matrix(NA, nrow=nrow(CELL_matrices[[1]]), ncol=ncol(CELL_matrices[[1]]))
# for (i in 1:nrow(CELL_matrices[[1]])) {
#   for (j in 1:ncol(CELL_matrices[[1]])) {
#     CELL_matrices_mean[i, j]=lapply(CELL_matrices, function(mtx) mtx[i, j]) %>% unlist() %>% mean(., na.rm=T)
#   }
#   print(paste0("Avg. Done with row : ", i, " in ", nrow(CELL_matrices[[1]])))
# }
# colnames(CELL_matrices_mean)=rownames(CELL_matrices_mean)=colnames(CELL_matrices[[1]])
# saveRDS(CELL_matrices_mean, "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_mean_ByRough.rds")
# 
# #####################################



# ### DON'T DO THIS!!!: Cluster the avg. cell-by-cell attention matrix
# # ... BECAUSE EACH MATRIX HAS DIFFERENT FEATURES (CELLS); THEY CANNOT BE POOLED
# #####################################
# ### 
# 
# CELL_matrices_mean=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_mean_ByRough.rds")
# sum(is.na(CELL_matrices_mean))/sum(!is.na(CELL_matrices_mean))
# # impute
# CELL_matrices_mean[is.na(CELL_matrices_mean)]=0
# 
# library("factoextra")
# library("FactoMineR")
# 
# ### For rough level
# # arrange the matrix
# Detailed_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_rough.rds")
# MATRIX=Reduce(`+`, Detailed_cor_matrix)/length(Detailed_cor_matrix)
# MATRIX_wo_zerosd=CELL_matrices_mean[which(apply(CELL_matrices_mean, 2, var)!=0), which(apply(CELL_matrices_mean, 2, var)!=0)]
# dim(MATRIX_wo_zerosd)
# # run pca and clustering
# pc_cell=prcomp(MATRIX_wo_zerosd,
#                center=TRUE,
#                scale.=TRUE)
# results=pc_cell$x
# factoextra::fviz_nbclust(results, FUNcluster=kmeans, k.max=100) # check the optimal number of clusters, n=32
# km=factoextra::eclust(results, "kmeans", hc_metric="eucliden", k=32)
# # plot
# km_data=km$clust_plot$data
# # plot_rough=
# ggplot(km_data, aes(x=x, y=y, color=cluster, shape=cluster, fill=cluster)) +
#   geom_point() +
#   geom_polygon(data=. %>% group_by(cluster) %>% do(.[chull(.[2:3]), ]), alpha=0.25, linewidth=0.1) +
#   ggrepel::geom_text_repel(aes(label=name), box.padding=0.5, min.segment.length=2) +
#   labs(x=km$clust_plot$labels$x, y=km$clust_plot$labels$y, title="Cluster plot - level1") +
#   # scale_shape_manual(values=c(1,0,2)) +
#   # scale_color_manual(values=paletteer::paletteer_d("ggthemes::few_Light")) +
#   # scale_fill_manual(values=paletteer::paletteer_d("ggthemes::few_Light")) +
#   theme_classic() +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=10),
#         legend.position="none")
# pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_celltype_2way_cluster_rough.pdf", height=3.5, width=3.5)
# plot_rough
# dev.off()
# 
# ### For inter level
# # arrange the matrix
# Inter_cor_matrix=readRDS("~/Project_PBMCage/For_or_From_Xu/CorMatrix_intermediate.rds")
# MATRIX=Reduce(`+`, Inter_cor_matrix)/length(Inter_cor_matrix)
# rownames(MATRIX)=colnames(MATRIX)=paste0(rownames(MATRIX),rep(" ",5)) # to avoid squeezing the labels when plotting
# MATRIX_wo_zerosd=MATRIX[which(apply(MATRIX, 2, var)!=0), which(apply(MATRIX, 2, var)!=0)]
# # remove DCs, Other cells, Monocytes
# MATRIX_wo_zerosd=MATRIX_wo_zerosd[!grepl("^DC\\.|^Other\\.|^Mono\\.",rownames(MATRIX_wo_zerosd)),]
# # run pca and clustering
# res.pca=PCA(MATRIX_wo_zerosd, ncp=10, graph=FALSE)
# res.hcpc=HCPC(res.pca, graph=FALSE, nb.clust=6)
# plot_dend=
#   fviz_dend(res.hcpc, 
#             lwd=0.25,
#             cex=0.8,
#             palette="jco",
#             rect=TRUE, rect_fill=TRUE,
#             rect_border="jco",
#             labels_track_height=15,
#             main="Cluster dendrogram - level2", 
#             xlab="", ylab="Height", 
#             ggtheme=theme_classic()
#   ) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         title=element_text(size=10),
#         legend.position="none")
# pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_celltype_2way_cluster_inter.pdf", height=3.5, width=7.25)
# plot_dend
# dev.off()
# 
# #####################################



### Cluster the cells by their association with other cells in each matrix
#####################################
###

### Load the matrices
CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_cellbycell_100donor_matrices.rds")
CELL_matrices_metadata=CELL_matrices_metadata[["MATRIX"]]

### Create function
celltype_clustering=function(matrix, cell.annot.df) {
  colnames(cell.annot.df)=c("cellid","Annot")
  # remove the empty columns (in the case of <2000 cells)
  matrix=matrix[!grepl("^0",colnames(matrix)),!grepl("^0",colnames(matrix))]
  # order and annotate
  df_=data.frame(cellid=colnames(matrix)) %>% left_join(., cell.annot.df, by="cellid")
  order_=order(df_$Annot)
  mtx_reordered=matrix[order_,order_] %>% as.matrix()
  colnames(mtx_reordered)=rownames(mtx_reordered)=make.names(df_$Annot[order_], unique=T)
  # remove the columns with sd=0
  MATRIX_wo_zerosd=mtx_reordered[which(apply(mtx_reordered, 2, var)!=0), which(apply(mtx_reordered, 2, var)!=0)]
  # run pca and clustering
  pc_cell=prcomp(MATRIX_wo_zerosd,
                 center=TRUE,
                 scale.=TRUE)
  results=pc_cell$x
  try_catch_df=function(results) {
    tryCatch({
      k_best=factoextra::fviz_nbclust(results, FUNcluster=kmeans, k.max=20) # check the optimal number of clusters
      km=factoextra::eclust(results, "kmeans", hc_metric="eucliden",
                            k=which(k_best$data$y==max(k_best$data$y)))
      df_cluster=data.frame(cell=names(km$cluster), cluster=km$cluster)
      return(df_cluster)
    }, error=function(msg) return(data.frame()))
  }
  df_cluster=try_catch_df(results)
  
  return(df_cluster)
}


### Run PCA clustering for each matrices with Annot.rough annotation
# load Annot.rough cell types
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% dplyr::select(cellid, Annot.rough)
# get the clustering DFs
Cluster_DFs=list()
for (i in 1:length(CELL_matrices_metadata)) {
  MATRIX_=CELL_matrices_metadata[[i]]
  Cluster_DFs[[i]]=celltype_clustering(MATRIX_, metadata_PBMCage)
  print(paste0("Done with ",i," in ",length(CELL_matrices_metadata)))
}
names(Cluster_DFs)=names(CELL_matrices_metadata)
saveRDS(Cluster_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_rough_eachDonor.rds")

### Analyze the PCA clustering at rough and plot
Cluster_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_rough_eachDonor.rds")
Cluster_DFs=lapply(1:length(Cluster_DFs), function(idx) Cluster_DFs[[idx]] %>% mutate(donor_id=names(Cluster_DFs)[idx]))
Cluster_DFs=data.table::rbindlist(Cluster_DFs, fill=TRUE)
Cluster_DFs_arranged=Cluster_DFs %>% mutate(cell=gsub("\\."," ",gsub("\\.[0-9]+$","",cell))) %>%
  group_by(donor_id, cluster, cell) %>%
  mutate(cellN_per_type=n()) %>%
  ungroup() %>% group_by(donor_id, cluster) %>%
  mutate(cellN_per_cluster=n()) %>%
  mutate(freq_per_type=cellN_per_type/cellN_per_cluster) %>%
  ungroup() %>% subset(!duplicated(.)) %>%
  group_by(cluster, donor_id) %>%
  mutate(dom_celltype=ifelse(freq_per_type==max(freq_per_type),cell,NA))
match_domcelltype_w_cluster=Cluster_DFs_arranged %>%
  subset(!is.na(dom_celltype)) %>%
  dplyr::select(cluster, donor_id, dom_celltype) %>%
  rename(dom_celltype_fill=dom_celltype)
fill_domcelltype=Cluster_DFs_arranged %>% 
  left_join(., match_domcelltype_w_cluster, by=c("cluster","donor_id"))

pca_result=
  ggplot(fill_domcelltype %>% subset(cell!="Unassigned"), 
         aes(x=dom_celltype_fill, y=freq_per_type, color=cell)) +
  geom_point(size=0.5, alpha=0.5) +
  geom_point(data=. %>%
               ungroup() %>% group_by(cell, dom_celltype_fill) %>%
               mutate(occur_in_donor=n()/99) %>%
               summarize(freq_per_type=mean(freq_per_type),
                         occur_in_donor=max(occur_in_donor)),
             aes(size=occur_in_donor), 
             alpha=1) +
  labs(x="Cluster by dominant cell type", y="Proportion in cluster", title=NULL) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  scale_size_continuous(range=c(0.5,4)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="right",
        legend.key.size=unit(5,"mm"),
        legend.spacing.y=unit(2.5,"mm"),
        legend.direction="vertical") +
  guides(color=guide_legend(override.aes=list(size=2), title="Cell type"),
         size=guide_legend(override.aes=list(color="black"), title="Dominating percent\nacross donors"))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_clustering_AnnotRough.pdf", height=3.5, width=6.5)
pca_result
dev.off()


### Run PCA clustering for each matrices with Annot.inter annotation
# load Annot.inter cell types
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% dplyr::select(cellid, Annot.inter)
# get the clustering DFs
Cluster_DFs=list()
for (i in 1:length(CELL_matrices_metadata)) {
  MATRIX_=CELL_matrices_metadata[[i]]
  Cluster_DFs[[i]]=celltype_clustering(MATRIX_, metadata_PBMCage)
  print(paste0("Done with ",i," in ",length(CELL_matrices_metadata)))
}
names(Cluster_DFs)=names(CELL_matrices_metadata)
saveRDS(Cluster_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_inter_eachDonor.rds")

### Analyze the PCA clustering at inter and plot
Cluster_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_inter_eachDonor.rds")
Cluster_DFs=lapply(1:length(Cluster_DFs), function(idx) Cluster_DFs[[idx]] %>% mutate(donor_id=names(Cluster_DFs)[idx]))
Cluster_DFs=data.table::rbindlist(Cluster_DFs, fill=TRUE)
Cluster_DFs_arranged=Cluster_DFs %>% mutate(cell=gsub("\\.[0-9]+$","",cell)) %>%
  group_by(donor_id, cluster, cell) %>%
  mutate(cellN_per_type=n()) %>%
  ungroup() %>% group_by(donor_id, cluster) %>%
  mutate(cellN_per_cluster=n()) %>%
  mutate(freq_per_type=cellN_per_type/cellN_per_cluster) %>%
  ungroup() %>% subset(!duplicated(.)) %>%
  group_by(cluster, donor_id) %>%
  mutate(dom_celltype=ifelse(freq_per_type==max(freq_per_type),cell,NA))
match_domcelltype_w_cluster=Cluster_DFs_arranged %>%
  subset(!is.na(dom_celltype)) %>%
  dplyr::select(cluster, donor_id, dom_celltype) %>%
  rename(dom_celltype_fill=dom_celltype)
fill_domcelltype=Cluster_DFs_arranged %>% 
  left_join(., match_domcelltype_w_cluster, by=c("cluster","donor_id"))

pca_result_inter=
  ggplot(fill_domcelltype %>% subset(cell!="Unassigned"), 
         aes(x=dom_celltype_fill, y=freq_per_type, color=cell)) +
  geom_point(size=0.5, alpha=0.5) +
  geom_point(data=. %>%
               ungroup() %>% group_by(cell, dom_celltype_fill) %>%
               mutate(occur_in_donor=n()/99) %>%
               summarize(freq_per_type=mean(freq_per_type),
                         occur_in_donor=max(occur_in_donor)),
             aes(size=occur_in_donor), 
             alpha=1) +
  ggrepel::geom_text_repel(data=. %>%
                             ungroup() %>% group_by(cell, dom_celltype_fill) %>%
                             mutate(occur_in_donor=n()/99) %>%
                             summarize(freq_per_type=mean(freq_per_type),
                                       occur_in_donor=max(occur_in_donor)) %>%
                             group_by(dom_celltype_fill) %>% 
                             subset(occur_in_donor>min(occur_in_donor) | freq_per_type>1/3*mean(freq_per_type)) %>%
                             subset(freq_per_type>0.05),
                           aes(label=cell), color="black", size=3, max.overlaps=18) +
  labs(x="Cluster by dominant cell type", y="Proportion in cluster", title=NULL) +
  scale_color_manual(values=c(paletteer::paletteer_d("ggsci::category20_d3"), 
                              paletteer::paletteer_d("ggsci::category20b_d3"))) +
  scale_size_continuous(range=c(0.5,4)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="right",
        legend.key.size=unit(5,"mm"),
        legend.spacing.y=unit(2.5,"mm"),
        legend.direction="horizontal") +
  guides(color=guide_legend(override.aes=list(size=2), title="Cell type", nrow=12),
         size=guide_legend(override.aes=list(color="black"), title="Dominating percent across donors"))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_clustering_AnnotInter.pdf", height=3.5, width=15)
pca_result_inter
dev.off()


### Run PCA clustering for each matrices with Annot.detailed annotation
# load Annot.detailed cell types
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% dplyr::select(cellid, Annot.detailed)
# get the clustering DFs
Cluster_DFs=list()
for (i in 1:length(CELL_matrices_metadata)) {
  MATRIX_=CELL_matrices_metadata[[i]]
  Cluster_DFs[[i]]=celltype_clustering(MATRIX_, metadata_PBMCage)
  print(paste0("Done with ",i," in ",length(CELL_matrices_metadata)))
}
names(Cluster_DFs)=names(CELL_matrices_metadata)
saveRDS(Cluster_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_detailed_eachDonor.rds")

### Analyze the clustering result at detailed
Cluster_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_cellbycell_clustering_detailed_eachDonor.rds")
Cluster_DFs=lapply(1:length(Cluster_DFs), function(idx) Cluster_DFs[[idx]] %>% mutate(donor_id=names(Cluster_DFs)[idx]))
Cluster_DFs=data.table::rbindlist(Cluster_DFs, fill=TRUE)
Cluster_DFs_arranged=Cluster_DFs %>% mutate(cell=gsub("\\.[0-9]+$","",cell)) %>%
  group_by(donor_id, cluster, cell) %>%
  mutate(cellN_per_type=n()) %>%
  ungroup() %>% group_by(donor_id, cluster) %>%
  mutate(cellN_per_cluster=n()) %>%
  mutate(freq_per_type=cellN_per_type/cellN_per_cluster) %>%
  ungroup() %>% subset(!duplicated(.)) %>%
  group_by(cluster, donor_id) %>%
  mutate(dom_celltype=ifelse(freq_per_type==max(freq_per_type),cell,NA))
match_domcelltype_w_cluster=Cluster_DFs_arranged %>%
  subset(!is.na(dom_celltype)) %>%
  dplyr::select(cluster, donor_id, dom_celltype) %>%
  rename(dom_celltype_fill=dom_celltype)
fill_domcelltype=Cluster_DFs_arranged %>% 
  left_join(., match_domcelltype_w_cluster, by=c("cluster","donor_id"))
fill_domcelltype_occur=fill_domcelltype %>%
  ungroup() %>% group_by(cell, dom_celltype_fill) %>%
  mutate(occur_in_donor=n()/99) %>%
  summarize(freq_per_type=mean(freq_per_type),
            occur_in_donor=max(occur_in_donor))
# use weighted cell proportion in each cluster as the association
fill_domcelltype_sorted=fill_domcelltype_occur %>% 
  ungroup() %>%
  mutate(rank=freq_per_type*occur_in_donor) %>%
  subset(!grepl("lin_infidel|pltcontamin",dom_celltype_fill) & !grepl("lin_infidel|pltcontamin",cell)) %>%
  tidyr::pivot_wider(names_from="cell", values_from="rank", id_cols="dom_celltype_fill") %>%
  tibble::column_to_rownames("dom_celltype_fill")
fill_domcelltype_sorted=fill_domcelltype_sorted[,colnames(fill_domcelltype_sorted) %in% rownames(fill_domcelltype_sorted)]
fill_domcelltype_sorted=fill_domcelltype_sorted[order(rownames(fill_domcelltype_sorted)),order(colnames(fill_domcelltype_sorted))]
unique(colnames(fill_domcelltype_sorted)==rownames(fill_domcelltype_sorted)) # check for dimnames for distance
fill_domcelltype_sorted_triangle=matrix(NA, 
                                        nrow=nrow(fill_domcelltype_sorted), ncol=ncol(fill_domcelltype_sorted), 
                                        dimnames=list(rownames(fill_domcelltype_sorted), colnames(fill_domcelltype_sorted))
                                        )
for (i in 1:nrow(fill_domcelltype_sorted)) { # take average by diagonal
  for (j in 1:ncol(fill_domcelltype_sorted)) {
    fill_domcelltype_sorted_triangle[i,j]=mean(c(fill_domcelltype_sorted[i,j], fill_domcelltype_sorted[j,i]), na.rm=T)
  }
}
fill_domcelltype_sorted_triangle[is.na(fill_domcelltype_sorted_triangle)]=0
# the association between cluster and cell type
association_=as.dist(fill_domcelltype_sorted_triangle)

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_CellByCell_clustering_AnnotDetailed.pdf", height=4.5, width=4.5)
qgraph::qgraph(association_, layout='spring', shape="rectangle", vsize=9, vsize2=1.5, label.cex=0.5, 
               edge.color="lightblue3", color="transparent", borders=FALSE, repulsion=1,
               labels=attr(association_,"Labels"), label.scale=FALSE)
dev.off()

#####################################






###############################################################################################################
##################################### Gene-by-Gene Attention Weight #####################################
###############################################################################################################

### Arrange the attention weights (5000 genes * 5000 cells * 99 donors)
#####################################
###

library(dplyr)

### Load the matrices
mydir="~/Project_PBMCage/For_or_From_Xu/AttentionMatrix_full/layer6_newCT_gene2gene.zip"
myfiles=paste0("layer6_newCT_2k/",c(0:98),".csv")
Matrices=lapply(myfiles, function(file) readr::read_delim(unzip(mydir, file), delim=" ", col_names=FALSE))

### Get the donor ids
donor_ids=lapply(Matrices, function(mtx) unique(mtx$X1)) %>% unlist() %>% unname()
names(Matrices)=donor_ids

### Match the gene symbols
metadata_PBMCage=read.delim("~/Project_PBMCage/For_or_From_Xu/GeneFeatures_ENStoSYMBOL.txt", sep="\t") %>%
  tibble::rownames_to_column("X2") %>%
  dplyr::select(X2, feature_name)
Matrices=lapply(Matrices, 
                function(mtx) {
                  mtx=mtx %>%
                    left_join(., metadata_PBMCage, by="X2") %>%
                    dplyr::select(-c(X1,X2)) %>% 
                    tibble::remove_rownames() %>% tibble::column_to_rownames("feature_name")
                  colnames(mtx)=rownames(mtx)
                  mtx
                })

### Save the matrices
saveRDS(Matrices, "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_genebygene_100donor_matrices.rds")

#####################################



### Calculate the distance of the matrices by distribution between donors
#####################################
###

library(dplyr)

### Load the cell-by-cell matrices
CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_genebygene_100donor_matrices.rds")
CELL_matrices=CELL_matrices_metadata
DF=list()
for (i in 1:length(CELL_matrices)) {
  mtx1=CELL_matrices[[i]] %>% as.matrix()
  dens1=density(as.vector(mtx1))
  vec1=dens1$y*diff(dens1$x)[1]/sum(dens1$y*diff(dens1$x)[1])
  for (j in 1:length(CELL_matrices)) {
    mtx2=CELL_matrices[[j]] %>% as.matrix()
    dens2=density(as.vector(mtx2))
    vec2=dens2$y*diff(dens2$x)[1]/sum(dens2$y*diff(dens2$x)[1])
    res=philentropy::KL(rbind(vec1, vec2), unit="log")
    df_=data.frame(kl_divergence=unname(res),
                   group=paste0(names(CELL_matrices)[i],"_vs_", names(CELL_matrices)[j])
    )
    DF=c(DF, list(df_))
    print(paste0("Done with j=",j))
  }
  print(paste0("====== Done with i=",i,". ======"))
}
DF_merged=data.table::rbindlist(DF) %>% mutate(compare.element="gene")
data.table::fwrite(DF_merged, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_GeneByGene_CompareMatrix_acrossDonor_KLdivergence.txt.gz", sep="\t")

#####################################



### Arrange the matrix distance between donors and do hierarchical clustering
#####################################
###

### Rearrange the distance
DF_merged=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_GeneByGene_CompareMatrix_acrossDonor_KLdivergence.txt.gz", sep="\t")
DF_merged=DF_merged %>%
  dplyr::select(kl_divergence, group) %>%
  mutate(group1=gsub("_vs_.*","",group), group2=gsub(".*_vs_","",group)) %>%
  dplyr::select(-group) %>%
  tidyr::pivot_wider(names_from="group2", values_from="kl_divergence") %>%
  tibble::column_to_rownames("group1")

### Plot distance
dist_m=as.dist(DF_merged)
dist_mi=1/dist_m # one over, as qgraph takes similarity matrices as input
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_CompareMatrixAcrossDonor_distance.pdf", height=3.5, width=3.5)
qgraph::qgraph(dist_mi, layout='spring', vsize=3)
dev.off()

### Determine the number of clusters
# ... turns out to be n=2
NbClust::NbClust(data=NULL, diss=dist_m, distance=NULL, method="ward.D2", index="silhouette")

### Perform clustering
hc.cut=factoextra::hcut(dist_m, k=2, hc_method="ward.D2")

### Plot the cluster dendrogram
# arrange the results
dend=as.dendrogram(hc.cut)
dend_data=ggdendro::dendro_data(dend, type="rectangle")
# prepare the meta info
metadata_PBMCage=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
metadata_PBMCage_donor_id.sex.age=metadata_PBMCage %>% dplyr::select(c(donor_id,sex,age)) %>% subset(!duplicated(.))

plot_dend=
  ggplot(dend_data$segments) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend), linewidth=0.25)+
  geom_text(data=dend_data$labels %>% 
              left_join(., metadata_PBMCage_donor_id.sex.age %>% dplyr::rename(label=donor_id)) %>%
              mutate(age_sex=ifelse(sex=="female",-age,age)) %>%
              mutate(age_sex_name=paste0(ifelse(sex=="female","F","M"),"  ",age,paste0(rep(" ",2), collapse=" "))), 
            aes(x, y, label=age_sex_name, color=age_sex),
            hjust=1, angle=90, size=3) +
  coord_cartesian(ylim=c(-0.02, 0.3)) +
  labs(x=NULL, y="Height", title="Cluster dendrogram") +
  scale_y_continuous(labels=function(x) format(x, scientific=FALSE)) +
  scale_color_gradient2(low="brown3", mid="white", high="dodgerblue3") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none",
        legend.spacing.x=unit(0.1,"cm"),
        legend.direction="vertical") +
  guides(color=guide_legend(labels=function(x) abs(x)))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_CompareMatrixAcrossDonor_dend.pdf", height=3.5, width=12)
plot_dend
dev.off()

### Check the distribution of sex & age in the clusters
# ... turns out to be no sig. diff. between clusters
df_=data.frame(donor_id=names(hc.cut$cluster), cluster=hc.cut$cluster) %>%
  left_join(., metadata_PBMCage_donor_id.sex.age, by="donor_id")
# plot the ages in each cluster
wilcox_res=
  wilcox.test((df_ %>% subset(cluster=="1"))$age, (df_ %>% subset(cluster=="2"))$age)
plot_age_dens=
  ggplot(df_, aes(x=age)) +
  geom_density(aes(group=cluster, fill=as.character(cluster)), alpha=0.5, linewidth=0.25) +
  theme_classic() +
  geom_vline(xintercept=df_ %>% group_by(cluster) %>% summarize_at("age", function(x) exp(mean(log(x)))) %>% .[["age"]],
             color=c("lightgreen","darkgreen"),
             linetype="dashed") +
  ggrepel::geom_text_repel(data=df_ %>% group_by(cluster) %>% summarize_at("age", function(x) exp(mean(log(x)))),
                           aes(x=age, y=0.06, label=floor(age)),
                           color="black", box.padding=0.5, point.padding=0.5, nudge_x=3, min.segment.length=0, size=3.5) +
  scale_fill_manual(values=c("lightgreen","darkgreen")) +
  theme_classic() +
  labs(x="age",y="density") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="none") +
  annotate(geom="text", x=30, y=0.045, label=paste0("p = ",signif(wilcox_res$p.value, 2)), size=3.5)

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_CompareMatrixAcrossDonor_clustering_AgeDistribution.pdf", height=3.5, width=2.5)
plot_age_dens
dev.off()

# plot the ages in each cluster
z_test_res=
  prop.test(c(nrow(df_ %>% subset(cluster=="1" & sex=="female")), nrow(df_ %>% subset(cluster=="2" & sex=="female"))), 
            n=c(nrow(df_ %>% subset(cluster=="1")), nrow(df_ %>% subset(cluster=="2"))),
            correct=F)
plot_sex_dens=
  ggplot(df_, aes(x=sex)) +
  geom_bar(aes(group=cluster, fill=as.character(cluster)), alpha=0.5, linewidth=0.25) +
  theme_classic() +
  scale_fill_manual(values=c("lightgreen","darkgreen")) +
  theme_classic() +
  labs(x="sex",y="count") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="right") +
  annotate(geom="text", x="male", y=50, label=paste0("p = ",round(z_test_res$p.value,2)), size=3.5) +
  guides(fill=guide_legend(title="cluster"))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_CompareMatrixAcrossDonor_clustering_SexDistribution.pdf", height=3.5, width=2.5)
plot_sex_dens
dev.off()

####################################



### Take the mean of matrices from all donors
#####################################
###

CELL_matrices=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_genebygene_100donor_matrices.rds")
# check
length(lapply(CELL_matrices, function(mtx) colnames(mtx)) %>% .[!duplicated(.)])

CELL_matrices_mean=matrix(NA, nrow=nrow(CELL_matrices[[1]]), ncol=ncol(CELL_matrices[[1]]))
for (i in 1:nrow(CELL_matrices[[1]])) {
  for (j in 1:ncol(CELL_matrices[[1]])) {
    CELL_matrices_mean[i, j]=lapply(CELL_matrices, function(mtx) mtx[i, j]) %>% unlist() %>% mean(., na.rm=T)
  }
  print(paste0("Avg. Done with row : ", i, " in ", nrow(CELL_matrices[[1]])))
}
colnames(CELL_matrices_mean)=rownames(CELL_matrices_mean)=colnames(CELL_matrices[[1]])
saveRDS(CELL_matrices_mean, "~/Project_PBMCage/For_or_From_Xu/AttentionWeight_genebygene_100donor_matrices_mean.rds")

#####################################



### Cluster the genes by their association with other genes in the mean matrix
#####################################
###

### Load the matrices
CELL_matrices_metadata=readRDS("~/Project_PBMCage/For_or_From_Xu/AttentionWeight_genebygene_100donor_matrices_mean.rds")

### Cluster the genes by the gene-by-gene attention matrix
unique(rownames(CELL_matrices_metadata)==colnames(CELL_matrices_metadata)) # check
# run pca and clustering
pc_cell=prcomp(CELL_matrices_metadata,
               center=TRUE,
               scale.=TRUE)
saveRDS(pc_cell, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x
try_catch_df=function(results) {
  tryCatch({
    k_best=factoextra::fviz_nbclust(results, FUNcluster=kmeans, k.max=20) # check the optimal number of clusters
    km=factoextra::eclust(results, "kmeans", hc_metric="eucliden",
                          k=which(k_best$data$y==max(k_best$data$y)))
    df_cluster=data.frame(gene=names(km$cluster), cluster=km$cluster)
    return(df_cluster)
  }, error=function(msg) return(data.frame()))
}
df_cluster=try_catch_df(results)
saveRDS(df_cluster, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor.rds")

### Plot the PCA result
Cluster_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor.rds")
pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
# draw fviz_pca_ind to get the percent of Dim1 and Dim2
# factoextra::fviz_pca_ind(pc_cell, repel=TRUE)
results=pc_cell$x
df_plot=results[,c(1,2,3)] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% left_join(., Cluster_DFs)
# check PC1, PC2, PC3 and found that the superplane to separate the cluster1&2 is unnatural, so consider only 1 cluster
plotly::plot_ly(x=df_plot$PC1, y=df_plot$PC2, z=df_plot$PC3, type="scatter3d", mode="markers", color=df_plot$cluster)
# plot the PC1-PC2 as only 1 cluster
pca_plot=
  ggplot(df_plot, aes(x=PC1, y=PC2)) +
  geom_point(size=0.1, alpha=0.5) +
  labs(x="PC1 (25.2%)",y="PC2 (5.7%)") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_clustering_PCA.pdf", height=2.5, width=2.5)
pca_plot
dev.off()

### Analyze the PCA clustering (2 clusters) and do enrichment
Cluster_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor.rds")
Cluster_genes=Cluster_DFs %>% split(.$cluster) %>% lapply(., function(df) df[["gene"]])
# only two clusters
cluster_enrich=lapply(Cluster_genes, function(gene) clusterProfiler::enrichGO(gene, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP"))
saveRDS(cluster_enrich, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_enrichGO.rds")

### Plot
plot_data=lapply(cluster_enrich, function(item) clusterProfiler::dotplot(item, showCategory=10)$data)
names(plot_data)=paste0("Cluster",names(plot_data))
# plot with ggplot2
MODULE_ENRICH_Plots=list()
for (i in 1:length(plot_data)) {
  module_=plot_data[[i]] %>% dplyr::select(Description, GeneRatio, p.adjust, Count)
  module_descrip_levels=levels(module_$Description)
  module_$Description=as.character(module_$Description)
  module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
  module_$Description=forcats::fct_relevel(module_$Description, formatters::wrap_string(module_descrip_levels, width=50, collapse="\n"))
  x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10
  plot_=
    ggplot(module_, aes(x=GeneRatio, y=Description, color=p.adjust, size=Count)) +
    geom_point() +
    scale_colour_gradient(limits=quantile(module_$p.adjust, na.rm=T)[c(1,5)],
                          low="brown4", high="pink",
                          labels=sprintf(fmt="%0.01e", quantile(module_$p.adjust, na.rm=T)[c(1,5)]),
                          breaks=quantile(module_$p.adjust, na.rm=T)[c(1,5)]) +
    scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
                                   round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                   max(module_$Count, na.rm=T)),
                          labels=c(min(module_$Count, na.rm=T),
                                   round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                   max(module_$Count, na.rm=T)),
                          range=c(2,6)) +
    guides(size=guide_legend(order=1)) +
    labs(title=names(plot_data)[i], y=NULL, x="GeneRatio") +
    theme_minimal() +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          title=element_text(size=11),
          plot.title.position="plot",
          plot.title=element_text(hjust=0.5),
          legend.position="right",
          legend.box="vertical",
          legend.direction="vertical",
          legend.key.width=unit(0.5,"cm"),
          legend.margin=margin(t=20, b=20)
    ) +
    guides(color=guide_colorbar(label.theme=element_text(angle=90, size=9, vjust=0.5, hjust=0.5),
                                title.theme=element_text(angle=90, size=10, vjust=1, hjust=0), order=1),
           size=guide_legend(label.theme=element_text(size=9),
                             title.theme=element_text(size=10), order=2))
  
  MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
}

plot_both=cowplot::plot_grid(plotlist=MODULE_ENRICH_Plots, align="hv", nrow=2)

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_GeneByGene_clustering_enrichments.pdf", width=7, height=7.5)
plot_both
dev.off()

### Rank the genes with PC1
pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x
PC1=results[,1] %>% as.data.frame() %>% rename(PC1=".") %>% tibble::rownames_to_column("gene") %>%
  arrange(desc(PC1)) %>%
  tibble::deframe()
res=clusterProfiler::gseGO(PC1, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
saveRDS(res, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
# plot the gseGO result
res=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
plot_data=clusterProfiler::dotplot(res, showCategory=10)$data
module_=plot_data %>% dplyr::select(Description, GeneRatio, NES, Count)
module_descrip_levels=levels(module_$Description)
module_$Description=as.character(module_$Description)
module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
module_$Description=forcats::fct_relevel(module_$Description, formatters::wrap_string(module_descrip_levels, width=50, collapse="\n"))
x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10
plot_=
  ggplot(module_, aes(x=GeneRatio, y=Description, color=NES, size=Count)) +
  geom_point() +
  scale_colour_gradient(limits=quantile(module_$NES, na.rm=T)[c(1,5)],
                        low="brown4", high="pink",
                        labels=sprintf(fmt="%0.2f", quantile(module_$NES, na.rm=T)[c(1,5)]),
                        breaks=quantile(module_$NES, na.rm=T)[c(1,5)]) +
  scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
                                 round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                 max(module_$Count, na.rm=T)),
                        labels=c(min(module_$Count, na.rm=T),
                                 round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                 max(module_$Count, na.rm=T)),
                        range=c(2,6)) +
  guides(size=guide_legend(order=1)) +
  labs(title=NULL, y=NULL, x="GeneRatio") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"),
        legend.key.height=unit(0.5,"cm"),
        legend.margin=margin(t=10, b=10)
  ) +
  guides(color=guide_colorbar(label.theme=element_text(angle=90, size=9, vjust=0.5, hjust=0.5),
                              title.theme=element_text(angle=90, size=10, vjust=1, hjust=0), order=1),
         size=guide_legend(label.theme=element_text(size=9),
                           title.theme=element_text(size=10), order=2))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_dotplot.pdf", height=3, width=6)
plot(plot_)
dev.off()

# map the gseGO terms to the 12 superclusters (by numbers)
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
number_of_mapped_terms=lapply(GOs_both, function(cluster) length(res@result$ID[res@result$ID %in% cluster])) %>%
  data.frame(supercluster=names(.), number=unlist(.)) %>%
  dplyr::select(supercluster, number)

# map the gseGO terms to the 12 superclusters (by p.adjust)
NES_abstract_of_mapped_terms=full_join(res@result,
                                       lapply(1:length(GOs_both), function(idx) data.frame(supercluster=rep(names(GOs_both)[idx],length(GOs_both[[idx]])), ID=GOs_both[[idx]])) %>% data.table::rbindlist(), 
                                       by=c("ID")) %>% 
  subset(!is.na(ID) & !is.na(supercluster) & !is.na(Description)) %>%
  mutate(orientation=ifelse(NES<0,"downreg","upreg")) %>%
  group_by(supercluster, orientation) %>%
  slice_min(p.adjust, n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(supercluster) %>%
  slice_max(abs(NES), n=1, with_ties=FALSE) %>%
  summarize(NES_abstract=mean(NES))

# merge the number and p.adjust
mtx_=full_join(NES_abstract_of_mapped_terms, number_of_mapped_terms, by="supercluster") %>% tibble::column_to_rownames("supercluster")
mtx_col_reorder=names(table(rownames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_col_reorder,]

# plot the number and p.adjust with heatmap
range(mtx_[,1]) # check
col_fun_cor=circlize::colorRamp2(c(-3, 0, 3), c("dodgerblue3", "white", "brown3"))
ht_cor=
  Heatmap(mtx_ %>% dplyr::select(NES_abstract),
          name="NES", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="NES", title_gp=gpar(fontface="plain"),
            legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=1, lwd=0.25),
          col=col_fun_cor,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          show_row_names=F,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title=NULL, column_title_gp=gpar(fontsize=11),
          width=1*unit(6, "mm"),
          height=nrow(mtx_)*unit(6, "mm")
  )

range(mtx_[,2]) # check
col_fun_pval=circlize::colorRamp2(c(0, 110), c("white","black"))
ht_pval=
  Heatmap(mtx_ %>% dplyr::select(number),
          name="termN", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="termN", title_gp=gpar(fontface="plain"),
            legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=1, lwd=0.25),
          col=col_fun_pval,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          show_row_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title=NULL, column_title_gp=gpar(fontsize=11),
          width=1*unit(6, "mm"),
          height=nrow(mtx_)*unit(6, "mm")
  )

ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_cor + ht_pval, heatmap_legend_side="left", merge_legend=TRUE)

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_MappingToSuperclusters.pdf", height=3, width=2.5)
ht_final
dev.off()

#####################################



### Map the genes with PC1 (excluding ribosomal-related) to 12 superclusters
#####################################
###

### Rank the genes with PC1
pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x
PC1=results[,1] %>% as.data.frame() %>% rename(PC1=".") %>% tibble::rownames_to_column("gene") %>%
  arrange(desc(PC1)) %>%
  tibble::deframe()
res=clusterProfiler::gseGO(PC1, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
saveRDS(res, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
# plot the gseGO result
res=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
plot_data=clusterProfiler::dotplot(res, showCategory=10)$data
module_=plot_data %>% dplyr::select(Description, GeneRatio, NES, Count)
module_descrip_levels=levels(module_$Description)
module_$Description=as.character(module_$Description)
module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
module_$Description=forcats::fct_relevel(module_$Description, formatters::wrap_string(module_descrip_levels, width=50, collapse="\n"))
x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10
plot_=
  ggplot(module_, aes(x=GeneRatio, y=Description, color=NES, size=Count)) +
  geom_point() +
  scale_colour_gradient(limits=quantile(module_$NES, na.rm=T)[c(1,5)],
                        low="brown4", high="pink",
                        labels=sprintf(fmt="%0.2f", quantile(module_$NES, na.rm=T)[c(1,5)]),
                        breaks=quantile(module_$NES, na.rm=T)[c(1,5)]) +
  scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
                                 round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                 max(module_$Count, na.rm=T)),
                        labels=c(min(module_$Count, na.rm=T),
                                 round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                 max(module_$Count, na.rm=T)),
                        range=c(2,6)) +
  guides(size=guide_legend(order=1)) +
  labs(title=NULL, y=NULL, x="GeneRatio") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"),
        legend.key.height=unit(0.5,"cm"),
        legend.margin=margin(t=10, b=10)
  ) +
  guides(color=guide_colorbar(label.theme=element_text(angle=90, size=9, vjust=0.5, hjust=0.5),
                              title.theme=element_text(angle=90, size=10, vjust=1, hjust=0), order=1),
         size=guide_legend(label.theme=element_text(size=9),
                           title.theme=element_text(size=10), order=2))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_dotplot.pdf", height=3, width=6)
plot(plot_)
dev.off()

# map the gseGO terms to the 12 superclusters (by numbers)
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
number_of_mapped_terms=lapply(GOs_both, function(cluster) length(res@result$ID[res@result$ID %in% cluster])) %>%
  data.frame(supercluster=names(.), number=unlist(.)) %>%
  dplyr::select(supercluster, number)

# map the gseGO terms to the 12 superclusters (by p.adjust)
NES_abstract_of_mapped_terms=full_join(res@result,
                                       lapply(1:length(GOs_both), function(idx) data.frame(supercluster=rep(names(GOs_both)[idx],length(GOs_both[[idx]])), ID=GOs_both[[idx]])) %>% data.table::rbindlist(), 
                                       by=c("ID")) %>% 
  subset(!is.na(ID) & !is.na(supercluster) & !is.na(Description)) %>%
  mutate(orientation=ifelse(NES<0,"downreg","upreg")) %>%
  group_by(supercluster, orientation) %>%
  slice_max(abs(NES), n=1, with_ties=FALSE) %>%
  # slice_max(abs(NES)*(-log10(p.adjust)), n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(supercluster) %>%
  slice_min(p.adjust, n=1, with_ties=FALSE) %>%
  summarize(NES_abstract=sum(NES))





### Map the gene-by-gene PCA clustering results to different cell types
#####################################
###

pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x[,c(1,2)]
markergenes=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels_updated.csv.gz", sep="\t")
markergenes_rough=markergenes %>% subset(celltype.level=="Annot.rough") %>% 
  dplyr::select(gene, Annot.rough, pct.1, pct.2) %>%
  right_join(., results %>% as.data.frame() %>% tibble::rownames_to_column("gene"), by="gene") %>%
  subset(!is.na(Annot.rough))
plot_pct_cut=
  ggplot(markergenes_rough, aes(x=PC1, y=PC2, color=pct.1)) +
  facet_wrap(~Annot.rough, nrow=2) +
  geom_point(size=0.1, shape=19) +
  scale_color_gradient(low="ivory3",high="brown3") +
  ggrepel::geom_text_repel(data=. %>% slice_max(pct.1, n=4, with_ties=FALSE, by="Annot.rough"), 
                           aes(label=gene), size=3, color="black", 
                           box.padding=0.4, point.padding=0.4) +
  theme_classic() +
  theme(title=element_text(size=9),
        strip.text=element_text(size=10),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.2,"cm"),
        legend.text=element_text(size=9)) +
  guides(color=guide_colorbar(title="% of expr. cells", title.position="top", 
                              title.theme=element_text(angle=90, vjust=1)))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes.pdf", height=3, width=6)
plot_pct_cut
dev.off()

#####################################



### Analyze similarity between the celltype-specific PC1-top genes (derived from marker genes and gene-by-gene matrix PCA clustering) across cell types
#####################################
###

### Load the celltype marker genes and the PC1 results
merged_results=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")
pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x[,c(1,2)]
# take marker genes
marker_genes_PC_res=merged_results %>% 
  subset(!grepl("^RP[SL]|^MT-",gene)) %>%
  subset(avg_log2FC>=0.5 & celltype.level=="Annot.rough") %>%
  dplyr::select(gene, Annot.rough) %>%
  right_join(., results %>% as.data.frame() %>% tibble::rownames_to_column("gene"), by="gene") %>%
  subset(!is.na(Annot.rough)) %>%
  slice_max(abs(PC1), n=300, by="Annot.rough", with_ties=FALSE) %>%
  split(.$Annot.rough) %>%
  .[names(table(merged_results$Annot.rough))] %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))

### Jaccard similarity analysis on PC1-top genes with cell-specific expr across cell types
# the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
SIM=SIM_Name=c()
celltypes_=names(marker_genes_PC_res)
for (i in 1:length(celltypes_)) {
  for (j in 1:length(celltypes_)) {
    similarity_result=jaccard(marker_genes_PC_res[[i]], marker_genes_PC_res[[j]])
    SIM=c(SIM, similarity_result)
    SIM_Name=c(SIM_Name, paste0(celltypes_[i],"_vs_",celltypes_[j]))
  }
}
names(SIM)=SIM_Name
# plot
Jaccard_SIM_pos=data.frame(pair=names(SIM), jaccard=SIM) %>%
  mutate(celltype_1=gsub("_vs_.*","",pair), celltype_2=gsub(".*_vs_","",pair)) %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype_1")
range(Jaccard_SIM_pos) # check
col_fun_down=circlize::colorRamp2(c(0, 1), c("white","brown3"))
ht_genesim=
  Heatmap(Jaccard_SIM_pos,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title="Similarity",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"), 
            direction="vertical", 
            title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun_down,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10), show_row_names=F,
          column_title="Individual genes",
          column_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(Jaccard_SIM_pos),
          width=unit(5, "mm")*ncol(Jaccard_SIM_pos))
ht_genesim=draw(ht_genesim, heatmap_legend_side="left")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes_individualGeneSimilarity.pdf", height=3, width=2.5)
draw(ht_genesim)
dev.off()

### Take the overlap between marker genes and top-PC1 genes for functional similarity analysis (mclusterSim)
marker_genes_PC_res_entrez=
  lapply(marker_genes_PC_res, 
         function(x) clusterProfiler::bitr(x, fromType="SYMBOL", OrgDb="org.Hs.eg.db", toType="ENTREZID") %>% .[["ENTREZID"]])
names(marker_genes_PC_res_entrez)=names(marker_genes_PC_res)
# analyze functional similarity with mclusterSim
d=GOSemSim::godata('org.Hs.eg.db', ont="BP")
res=GOSemSim::mclusterSim(marker_genes_PC_res_entrez, semData=d, measure="Wang", combine="BMA")
saveRDS(res, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes_mclusterSim_results.rds")

### Plot
SIM=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes_mclusterSim_results.rds")
range(SIM) # check
col_fun_down=circlize::colorRamp2(c(0.6, 1), c("white","brown3"))
ht_genesim=
  Heatmap(SIM,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title="Similarity",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"), 
            direction="vertical", 
            title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun_down,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10), show_row_names=F,
          column_title="Functional clusters",
          column_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(SIM),
          width=unit(5, "mm")*ncol(SIM))
ht_genesim=draw(ht_genesim, heatmap_legend_side="left")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes_mClusterSimilarity.pdf", height=3, width=2.5)
draw(ht_genesim)
dev.off()

#####################################



### Confirm the key role of MYC with gene-by-gene matrix
#####################################
###

pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x[,c(1,2)]
# mark the cell-cell communication key molecules
key_genes=c("HLA-A","HLA-B","HLA-C","HLA-E","CD8B","CD8A",
            "KLRK1","KLRD1","KLRC1","CD74","CXCR4","CD44",
            "MIF","CLEC1A","CLEC1B","LCK","CD99")
results_marked=results %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
  mutate(mark=ifelse(grepl("^RP[SL]",gene) & !(gene %in% key_genes), "ribosome", 
                     ifelse(gene %in% key_genes, "communication", "others")))
results_marked$mark=forcats::fct_relevel(results_marked$mark, c("communication","ribosome","others"))
plot_=
  ggplot(results_marked, aes(x=PC1, y=PC2, color=mark, size=mark)) +
  ggrastr::geom_point_rast() +
  ggrepel::geom_text_repel(data=. %>% subset(gene %in% key_genes), aes(label=gene), color="black", size=3.5, force=10) +
  theme_classic() +
  scale_color_manual(values=c("brown1","dodgerblue1","snow3")) +
  scale_size_manual(values=c(1,1,0.05)) +
  theme(title=element_text(size=11),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        legend.position="top",
        legend.text=element_text(size=11)) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2), ncol=1),
         size="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToCellCellCommunicationGenes.pdf", height=5.5, width=2.5)
plot(plot_)
dev.off()

#####################################



### Map the gene-by-gene PCA clustering results to TFs
#####################################
###

library(dplyr)
library(ggplot2)

pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
results=pc_cell$x[,c(1,2)] %>% as.data.frame() %>% tibble::rownames_to_column("gene")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB=c("NFKB1","NFKB2"),KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
markergenes=results %>% subset(gene %in% unlist(chosen_TFs))

# plot_pct_cut=
  ggplot(results, aes(x=PC1, y=PC2)) +
  geom_point(data=. %>% subset(!(gene %in% unlist(chosen_TFs))), size=0.1, shape=1, alpha=0.1) +
  geom_point(data=markergenes, size=2, color="brown3") +
  ggrepel::geom_text_repel(data=markergenes, 
                           aes(label=gene), size=3, color="black", 
                           box.padding=0.4, point.padding=0.4) +
  theme_classic() +
  theme(title=element_text(size=9),
        strip.text=element_text(size=10),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.key.width=unit(0.2,"cm"),
        legend.text=element_text(size=9)) +
  guides(color=guide_colorbar(title="% of expr. cells", title.position="top", 
                              title.theme=element_text(angle=90, vjust=1)))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_MappingToMarkerGenes.pdf", height=3, width=6)
plot_pct_cut
dev.off()

#####################################
