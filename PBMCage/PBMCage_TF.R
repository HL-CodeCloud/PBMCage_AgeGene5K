
###############################################################################################################
##################################### Bcell TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for B cells
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="B cells")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])
  
  print(paste0("MetaObject ", i, " Creation: Starts......"))
  
  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)
  
  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)
  
  # run the metacell combination
  age_subset=MetacellsByGroups(
    age_subset,
    group.by=c("Annot.detailed","age"),
    reduction="umap",
    k=3,
    min_cells=10,
    max_shared=1,
    target_metacells=25,
    verbose=TRUE,
    slot="counts",
    ident.group="Annot.detailed")
  
  # get the metacell obj
  obj_age[[i]]=GetMetacellObject(age_subset)
  obj_age[[i]]=JoinLayers(obj_age[[i]])
  obj_age[[i]]=NormalizeData(obj_age[[i]])
  # Check the metacells
  print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
  print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))
  
  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_Bcells.rds")

#####################################



### Transcriptional factors of B cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_Bcells.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcells.rds")

#####################################



### WordCloud of the TF scores regardless of ages in B cells
#####################################
###

### Get the top 10 TF genes in each cluster * age
library(dplyr)
library(tidyr)
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcells.rds")
ages=unique(DF_CD8T$age)
ages=sort(ages)
subset_=DF_CD8T %>%
  subset(pvalue<=0.05) %>%
  slice_max(order_by=score, n=10, by=c(cluster, age))
Gene_to_mark=unique(subset_$source)
### Process the dataframe for plotting
# !: no observable difference in score along ages
# use the sig. and top. genes and,
# ... color by the cluster which has the highest score for each source
unique(subset_$cluster) # check
subset_=subset(subset_, pvalue<=0.05 & source %in% Gene_to_mark) %>%
  subset(!grepl("B\\.trans",cluster)) %>% # remove B.trans as it is a minor population
  group_by(source, cluster) %>% summarize_at("score", mean) %>% 
  mutate(Annot.inter=ifelse(grepl("B\\.inter", cluster), "B.inter", 
                            ifelse(grepl("B\\.mem", cluster), "B.mem",
                                   ifelse(grepl("B\\.naive", cluster), "B.naive", "B.pbpc")))) %>%
  group_by(source, Annot.inter) %>% summarize_at("score", mean)

subset_$Annot.inter=forcats::fct_relevel(subset_$Annot.inter, c("B.naive","B.inter","B.mem","B.pbpc"))

### Plot
library(ggwordcloud)
plot_=
  ggplot(subset_, aes(label=source, size=score)) +
  facet_wrap(~Annot.inter) +
  geom_text_wordcloud(area_corr=TRUE) +
  scale_size_area(max_size=15) +
  theme_minimal() +
  theme(strip.text=element_text(size=10),
        strip.background=element_rect(fill=paletteer::paletteer_d("ggsci::category20c_d3")[11]))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_Bcells_wordcloud.pdf", height=4, width=5)
plot(plot_)
dev.off()

### Calculate the dtw distance in ordered TFs between neighbouring ages in each celltype
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcells.rds")
DF_CD8T=DF_CD8T %>% 
  mutate(cluster=ifelse(grepl("B\\.inter", cluster), "B.inter", 
                        ifelse(grepl("B\\.mem", cluster), "B.mem",
                               ifelse(grepl("B\\.trans", cluster), "B.trans", 
                                      ifelse(grepl("B\\.naive", cluster), "B.naive", "B.pbpc")))))
DF_CD8T$cluster=forcats::fct_relevel(DF_CD8T$cluster, c("B.trans","B.naive","B.inter","B.mem","B.pbpc"))

# make a function to go through the ages for dtw distance
windowsliding_dtw=function(vector_list) {
  dtw_dis_list=list()
  for (j in 2:length(vector_list)) {
    dtw_dis_list[[j-1]]=dtw::dtw(vector_list[[j]], vector_list[[j-1]], distance.only=TRUE)$distance
  }
  return(dtw_dis_list)
}
# order the TFs
source_idx=DF_CD8T %>% arrange(source); source_idx=unique(source_idx$source)
# remove the celltype with only 1 age
celltpye_chosen=DF_CD8T %>% select(cluster, age) %>% subset(!duplicated(.)) %>% ungroup() %>% 
  select(-age) %>% group_by(cluster) %>% dplyr::count() %>%
  subset(n>=2)
celltpye_chosen=celltpye_chosen$cluster
# arrange the df
DF_CD8T_clean=DF_CD8T %>% subset(cluster %in% celltpye_chosen) %>%
  mutate(order=match(source, source_idx)) %>% split(.$cluster) %>%
  lapply(., function(df) df %>% arrange(order) %>% select(age, source, score)) %>%
  lapply(., function(df) df %>% split(.$age))
# calculate
DTW_df_list=list()
for(i in 1:length(DF_CD8T_clean)) {
  celltype_=names(DF_CD8T_clean)[i]
  
  score_everyAge=DF_CD8T_clean[[i]] %>% lapply(., function(df) df$score)
  dtwlist=windowsliding_dtw(score_everyAge)
  dtwlist=unlist(dtwlist)
  names(dtwlist)=names(score_everyAge)[2:length(score_everyAge)]
  
  DTW_df_list[[i]]=data.frame(age=names(dtwlist), distance=dtwlist, celltypes=celltype_)
}
DTW_df=data.table::rbindlist(DTW_df_list)
DTW_df$age=as.integer(DTW_df$age)
DTW_df$celltypes=forcats::fct_relevel(DTW_df$celltypes, c("B.trans","B.naive","B.inter","B.mem","B.pbpc"))

plot_=
  ggplot(subset(DTW_df, celltypes!="B.trans" & celltypes!="B.inter"), # remove B.trans and B.inter to match celltypes in CD8T
         aes(x=age, y=distance, linetype=celltypes)) +
  geom_point(shape="+", size=1, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dotdash","dotted","dashed")) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[1]) +
  labs(x="age", y="DTW distance") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position=c(0.5,0.9),
        legend.key.width=unit(1.5,"cm"),
        legend.direction="horizontal") +
  scale_x_continuous(n.breaks=5) +
  ylim(c(200,1600)) +
  geom_smooth(aes(group=celltypes), method="loess", se=F, linewidth=0.8,
              color=paletteer::paletteer_d("ggsci::category20_d3")[1]) +
  guides(linetype=guide_legend(override.aes=list(linewidth=0.5, alpha=1), ncol=2, byrow=TRUE, title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_Bcells_DTWdistance.betweenAge.pdf", width=3.6, height=4)
plot(plot_)
dev.off()

#####################################



### Reanalyze the B cell TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcells.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("B\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("B\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("B\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("B\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcell_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Bcell_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("ZEB2","XBP1","PRDM1","AICDA","IKZF3","PAX5","TBX21","TCF3","EBF1","STAT5A","BACH2","TCF4","IRF4","NOTCH2","HOXC4")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cluster!="inter" & cluster!="pb") %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("trans","naive","mem","pc"))
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free") +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-2.5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="B cells") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_Bcell_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### CD4T TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for CD4 T cells
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="CD4T cells")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])
  
  print(paste0("MetaObject ", i, " Creation: Starts......"))
  
  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)
  
  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)
  
  # run the metacell combination
  age_subset=MetacellsByGroups(
    age_subset,
    group.by=c("Annot.detailed","age"),
    reduction="umap",
    k=3,
    min_cells=10,
    max_shared=1,
    target_metacells=25,
    verbose=TRUE,
    slot="counts",
    ident.group="Annot.detailed")
  
  # get the metacell obj
  obj_age[[i]]=GetMetacellObject(age_subset)
  obj_age[[i]]=JoinLayers(obj_age[[i]])
  obj_age[[i]]=NormalizeData(obj_age[[i]])
  # Check the metacells
  print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
  print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))
  
  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_CD4T.rds")

#####################################



### Transcriptional factors of CD4T cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_CD4T.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T.rds")

#####################################



### WordCloud of the TF scores regardless of ages in CD4T cells
#####################################
###

### Get the top 10 TF genes in each cluster * age
library(dplyr)
library(tidyr)
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T.rds")
ages=unique(DF_CD8T$age)
ages=sort(ages)
subset_=DF_CD8T %>%
  subset(pvalue<=0.05) %>%
  slice_max(order_by=score, n=10, by=c(cluster, age))
Gene_to_mark=unique(subset_$source)
### Process the dataframe for plotting
# !: no observable difference in score along ages
# use the sig. and top. genes and,
# ... color by the cluster which has the highest score for each source
unique(subset_$cluster) # check
subset_=subset(subset_, pvalue<=0.05 & source %in% Gene_to_mark) %>%
  subset(!grepl("CD4T\\.prolif",cluster)) %>% # remove B.trans as it is a minor population
  group_by(source, cluster) %>% summarize_at("score", mean) %>% 
  mutate(Annot.inter=ifelse(grepl("CD4T\\.naive", cluster), "CD4T.naive", 
                            ifelse(grepl("CD4T\\.Tcm", cluster), "CD4T.Tcm",
                                   ifelse(grepl("CD4T\\.Tem", cluster), "CD4T.Tem", "CD4T.Treg")))) %>%
  group_by(source, Annot.inter) %>% summarize_at("score", mean)

subset_$Annot.inter=forcats::fct_relevel(subset_$Annot.inter, c("CD4T.naive","CD4T.Tcm","CD4T.Tem","CD4T.Treg"))

### Plot
library(ggwordcloud)
plot_=
  ggplot(subset_, aes(label=source, size=score)) +
  facet_wrap(~Annot.inter) +
  geom_text_wordcloud(area_corr=TRUE) +
  scale_size_area(max_size=15) +
  theme_minimal() +
  theme(strip.text=element_text(size=10),
        strip.background=element_rect(fill=paletteer::paletteer_d("ggsci::category20c_d3")[7]))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD4T_wordcloud.pdf", height=4, width=4)
plot(plot_)
dev.off()

### Calculate the dtw distance in ordered TFs between neighbouring ages in each celltype
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T.rds")
DF_CD8T=DF_CD8T %>% 
  mutate(cluster=ifelse(grepl("CD4T\\.naive", cluster), "CD4T.naive", 
                        ifelse(grepl("CD4T\\.Tcm", cluster), "CD4T.Tcm",
                               ifelse(grepl("CD4T\\.Tem", cluster), "CD4T.Tem", 
                                      ifelse(grepl("CD4T\\.Treg", cluster), "CD4T.Treg", "CD4T.prolif")))))
DF_CD8T$cluster=forcats::fct_relevel(DF_CD8T$cluster, c("CD4T.naive","CD4T.Tcm","CD4T.Tem","CD4T.Treg"))

# make a function to go through the ages for dtw distance
windowsliding_dtw=function(vector_list) {
  dtw_dis_list=list()
  for (j in 2:length(vector_list)) {
    dtw_dis_list[[j-1]]=dtw::dtw(vector_list[[j]], vector_list[[j-1]], distance.only=TRUE)$distance
  }
  return(dtw_dis_list)
}
# order the TFs
source_idx=DF_CD8T %>% arrange(source); source_idx=unique(source_idx$source)
# remove the celltype with only 1 age
celltpye_chosen=DF_CD8T %>% select(cluster, age) %>% subset(!duplicated(.)) %>% ungroup() %>% 
  select(-age) %>% group_by(cluster) %>% dplyr::count() %>%
  subset(n>=2)
celltpye_chosen=celltpye_chosen$cluster
# arrange the df
DF_CD8T_clean=DF_CD8T %>% subset(cluster %in% celltpye_chosen) %>%
  mutate(order=match(source, source_idx)) %>% split(.$cluster) %>%
  lapply(., function(df) df %>% arrange(order) %>% select(age, source, score)) %>%
  lapply(., function(df) df %>% split(.$age))
# calculate
DTW_df_list=list()
for(i in 1:length(DF_CD8T_clean)) {
  celltype_=names(DF_CD8T_clean)[i]
  
  score_everyAge=DF_CD8T_clean[[i]] %>% lapply(., function(df) df$score)
  dtwlist=windowsliding_dtw(score_everyAge)
  dtwlist=unlist(dtwlist)
  names(dtwlist)=names(score_everyAge)[2:length(score_everyAge)]
  
  DTW_df_list[[i]]=data.frame(age=names(dtwlist), distance=dtwlist, celltypes=celltype_)
}
DTW_df=data.table::rbindlist(DTW_df_list)
DTW_df$age=as.integer(DTW_df$age)

plot_=
  ggplot(subset(DTW_df, celltypes!="CD4T.prolif"), # remove CD4T.prolif as it has few observations
         aes(x=age, y=distance, linetype=celltypes)) +
  geom_point(shape="+", size=1, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dotdash","dotted","dashed")) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
  labs(x="age", y="DTW distance") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position=c(0.5,0.9),
        legend.key.width=unit(1.5,"cm"),
        legend.direction="horizontal") +
  scale_x_continuous(n.breaks=5) +
  geom_smooth(aes(group=celltypes), method="loess", se=F, linewidth=0.8,
              color=paletteer::paletteer_d("ggsci::category20_d3")[2]) +
  guides(linetype=guide_legend(override.aes=list(linewidth=0.5, alpha=1), ncol=2, byrow=TRUE, title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD4T_DTWdistance.betweenAge.pdf", width=3.5, height=4)
plot(plot_)
dev.off()

#####################################



### Reanalyze the CD4T TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD4T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD4T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD4T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD4T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD4T_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("BATF","RFX5","KLF14","TGIF2","ETV1","RUNX2","NKX6-1","EOMES","SPI1","FLI1","KLF5","TBX5","TGIF1","SMAD3","RUNX1",
             "GRHL2","ZBTB12","STAT4","SMAD2","SOX4","ETS1","GFY","PTF1A","NRF1","YY1","ELK1","IRF2","TFAP2A","TCF4","TAL1",
             "FOSL2","JUN","PRDM1","GATA3","BACH2","GABPA","ETS1","ELK4","ELF1","ELK1","NFYA","SP1","NRF1")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cluster!="Treg") %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("naive","prolif","Tcm","Tem"))
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free") +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-2.5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="CD4T cells") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD4T_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### CD8T TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for CD8 T cells
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="CD8T cells")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])

  print(paste0("MetaObject ", i, " Creation: Starts......"))

  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)

  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)

  # run the metacell combination
  age_subset=MetacellsByGroups(
    age_subset,
    group.by=c("Annot.detailed","age"),
    reduction="umap",
    k=3,
    min_cells=10,
    max_shared=1,
    target_metacells=25,
    verbose=TRUE,
    slot="counts",
    ident.group="Annot.detailed")

  # get the metacell obj
  obj_age[[i]]=GetMetacellObject(age_subset)
  obj_age[[i]]=JoinLayers(obj_age[[i]])
  obj_age[[i]]=NormalizeData(obj_age[[i]])
  # Check the metacells
  print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
  print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))

  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_CD8T.rds")

#####################################



### Transcriptional factors of CD8T cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_CD8T.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")

#####################################



### Analyze the TF results of CD8T
#####################################
###

### Get the top 10 TF genes in each cluster * age
library(dplyr)
library(tidyr)
library(rstatix)
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")
ages=unique(DF_CD8T$age)
ages=sort(ages)
subset_=DF_CD8T %>%
  subset(pvalue<=0.05) %>%
  slice_max(order_by=score, n=10, by=c(cluster, age))
Gene_to_mark=unique(subset_$source)
# get the dataframe with the sig. and top genes
DF_subset_forPlot=subset(DF_CD8T, pvalue<=0.05 & source %in% Gene_to_mark)

### Make a function for plotting on ages
geneset_scoring_forAge=function(dataframe, celltype_as_LegendName, scale) {
  library(ComplexHeatmap)
  df_full_list=list()
  CD8T_NK_subsets_celltypes=unique(dataframe$cluster)
  
  for (i in 1:length(CD8T_NK_subsets_celltypes)) {
    cluster_subset=dataframe %>% subset(cluster==CD8T_NK_subsets_celltypes[i]) %>% select(-pvalue) %>% arrange(age)
    
    matrix_orig=pivot_wider(cluster_subset, names_from=age, values_from=score)
    df_to_draw=matrix_orig %>% select(-c(cluster, source)) %>% as.data.frame()
    rownames(df_to_draw)=matrix_orig$source
    
    age_not_there=unique(dataframe$age)[!(unique(dataframe$age) %in% unique(cluster_subset$age))]
    if (length(age_not_there)!=0) {
      for (j in age_not_there) {
        matrix_orig[[j]]=NA
      }
    }
    matrix_=matrix_orig[, sort(unique(dataframe$age))]
    matrix_=as.data.frame(matrix_)
    rownames(matrix_)=matrix_orig$source
    df_full_list[[i]]=matrix_
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges, 0, ranges), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-rowMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/2, 0, ranges/2), c("#00007E", "white", "#7E0000"))}
  
  p_age=list()
  codes=""
  for (i in 1:length(df_full_list)) {
    df_to_draw=as.data.frame(df_full_list[[i]])
    if (scale==T) {df_to_draw=t(scale(t(df_to_draw)))} else {df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)}
    
    if (i==1) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                name=paste0("TF score of ", celltype_as_LegendName),
                heatmap_legend_param=list(title_position="leftcenter-rot",
                                          title_gp=gpar(fontsize=8, fontface="bold")),
                column_title=NULL,
                row_title=CD8T_NK_subsets_celltypes[i],
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=T,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8),
                column_names_side="top",
                width=ncol(df_to_draw)*unit(3,"mm"), 
                height=nrow(df_to_draw)*unit(3,"mm"))}
    if (i>=2) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                show_heatmap_legend=FALSE,
                column_title=NULL,
                row_title=CD8T_NK_subsets_celltypes[i],
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=F,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
                column_names_gp=gpar(fontsize=8),
                column_names_side="top",
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

### Plot the TF genes
CD8T_TF_plots=geneset_scoring_forAge(dataframe=DF_subset_forPlot, celltype_as_LegendName="CD8 T cells", scale=T)

# pdf("~/Project_PBMCage/Plots/Scoring_CD8T_TFgenes.pdf", width=12.5, height=40)
# draw(CD8T_TF_plots)
# dev.off()

pdf("~/Project_PBMCage/Plots/Scoring_CD8T_TFgenes_scale.pdf", width=12.5, height=40)
draw(CD8T_TF_plots)
dev.off()

#####################################



### WordCloud of the TF scores regardless of ages in CD8T
#####################################
###

### Get the top 10 TF genes in each cluster * age
library(dplyr)
library(tidyr)
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")
ages=unique(DF_CD8T$age)
ages=sort(ages)
subset_=DF_CD8T %>%
  subset(pvalue<=0.05) %>%
  slice_max(order_by=score, n=10, by=c(cluster, age))
Gene_to_mark=unique(subset_$source)
### Process the dataframe for plotting
# !: no observable difference in score along ages
# use the sig. and top. genes and,
# ... color by the cluster which has the highest score for each source
subset_=subset(subset_, pvalue<=0.05 & source %in% Gene_to_mark) %>%
  group_by(source, cluster) %>% summarize_at("score", mean) %>% 
  mutate(Annot.inter=ifelse(grepl("CD8T\\.CTL", cluster), "CD8T.CTL", 
                            ifelse(grepl("CD8T\\.Tem", cluster), "CD8T.Tem",
                                   ifelse(grepl("CD8T\\.naive", cluster), "CD8T.naive", "CD8T.Tcm")))) %>%
  group_by(source, Annot.inter) %>% summarize_at("score", mean)
subset_$Annot.inter=forcats::fct_relevel(subset_$Annot.inter, c("CD8T.naive","CD8T.Tcm","CD8T.Tem","CD8T.CTL"))

### Plot
library(ggwordcloud)
plot_=
  ggplot(subset_, aes(label=source, size=score)) +
  facet_wrap(~Annot.inter) +
  geom_text_wordcloud(area_corr=TRUE) +
  scale_size_area(max_size=15) +
  theme_minimal() +
  theme(strip.text=element_text(size=10),
        strip.background=element_rect(fill=paletteer::paletteer_d("ggsci::category20c_d3")[13]))
  
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD8T_wordcloud.pdf", height=4, width=4)
plot(plot_)
dev.off()

### Calculate the dtw distance in ordered TFs between neighbouring ages in each celltype
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")
DF_CD8T=DF_CD8T %>% mutate(cluster=ifelse(grepl("CD8T\\.CTL", cluster), "CD8T.CTL", 
                                          ifelse(grepl("CD8T\\.Tem", cluster), "CD8T.Tem",
                                                 ifelse(grepl("CD8T\\.naive", cluster), "CD8T.naive", "CD8T.Tcm"))))
DF_CD8T$cluster=forcats::fct_relevel(DF_CD8T$cluster, c("CD8T.naive","CD8T.Tcm","CD8T.Tem","CD8T.CTL"))
# make a function to go through the ages for dtw distance
windowsliding_dtw=function(vector_list) {
  dtw_dis_list=list()
  for (j in 2:length(vector_list)) {
    dtw_dis_list[[j-1]]=dtw::dtw(vector_list[[j]], vector_list[[j-1]], distance.only=TRUE)$distance
  }
  return(dtw_dis_list)
}
# order the TFs
source_idx=DF_CD8T %>% arrange(source); source_idx=unique(source_idx$source)
# remove the celltype with only 1 age
celltpye_chosen=DF_CD8T %>% select(cluster, age) %>% subset(!duplicated(.)) %>% ungroup() %>% 
  select(-age) %>% group_by(cluster) %>% dplyr::count() %>%
  subset(n>=2)
celltpye_chosen=celltpye_chosen$cluster
# arrange the df
DF_CD8T_clean=DF_CD8T %>% subset(cluster %in% celltpye_chosen) %>%
  mutate(order=match(source, source_idx)) %>% split(.$cluster) %>%
  lapply(., function(df) df %>% arrange(order) %>% select(age, source, score)) %>%
  lapply(., function(df) df %>% split(.$age))
# calculate
DTW_df_list=list()
for(i in 1:length(DF_CD8T_clean)) {
  celltype_=names(DF_CD8T_clean)[i]
  
  score_everyAge=DF_CD8T_clean[[i]] %>% lapply(., function(df) df$score)
  dtwlist=windowsliding_dtw(score_everyAge)
  dtwlist=unlist(dtwlist)
  names(dtwlist)=names(score_everyAge)[2:length(score_everyAge)]
  
  DTW_df_list[[i]]=data.frame(age=names(dtwlist), distance=dtwlist, celltypes=celltype_)
}
DTW_df=data.table::rbindlist(DTW_df_list)
DTW_df$age=as.integer(DTW_df$age)

plot_=
  ggplot(subset(DTW_df, celltypes!="CD8T.CTL"), # checked CD8T.CTL has too few age points, so removed
         aes(x=age, y=distance, linetype=celltypes)) +
  geom_point(shape="+", size=1, alpha=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dotdash","dotted","dashed")) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[3]) +
  labs(x="age", y="DTW distance") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position=c(0.5,0.9),
        legend.key.width=unit(1.5,"cm"),
        legend.direction="horizontal") +
  scale_x_continuous(n.breaks=5) +
  ylim(c(200,1600)) +
  geom_smooth(aes(group=celltypes), method="loess", se=F, linewidth=0.8,
              color=paletteer::paletteer_d("ggsci::category20_d3")[3]) +
  guides(linetype=guide_legend(override.aes=list(linewidth=0.5, alpha=1), ncol=2, byrow=TRUE, title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD8T_DTWdistance.betweenAge.pdf", width=3.6, height=4)
plot(plot_)
dev.off()

#####################################



### Enrich the target genes of the sig.&Top TFs
#####################################
###

### Get the top 10 TF genes in each cluster * age
library(decoupleR)
library(dplyr)
library(clusterProfiler)
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")
ages=unique(DF_CD8T$age)
ages=sort(ages)
subset_=DF_CD8T %>%
  subset(pvalue<=0.05) %>%
  slice_max(order_by=score, n=10, by=c(cluster, age))
Gene_to_mark=unique(subset_$source)

### Get the targets
net=get_collectri(organism="human", split_complexes=FALSE)
net_subset=subset(net, source %in% Gene_to_mark)
saveRDS(net_subset, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T_TargetGenes.rds")
Genes_list=list()
for (i in 1:length(Gene_to_mark)) {
  Genes_list_actORreg=list()
  net_subset_source_act=subset(net_subset, source==Gene_to_mark[i] & mor==1)
  net_subset_source_reg=subset(net_subset, source==Gene_to_mark[i] & mor==(-1))
  Genes_list_actORreg=c(list(net_subset_source_act$target), list(net_subset_source_reg$target))
  names(Genes_list_actORreg)=paste0(Gene_to_mark[i],"_", c("act.","reg."))
  Genes_list=c(Genes_list, Genes_list_actORreg)
}

### Erich
EGO_compare=compareCluster(geneCluster=Genes_list, 
                           fun=enrichGO, 
                           OrgDb="org.Hs.eg.db", 
                           keyType="SYMBOL", 
                           pvalueCutoff=0.05,
                           minGSSize=1,
                           maxGSSize=500,
                           ont="BP")
saveRDS(EGO_compare, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T_TargetGenes_enrichments.rds")
EGO_compare_filter=gofilter(EGO_compare, level=4)

### Plot
plot_EGO=
  dotplot(EGO_compare_filter, showCategory=5, font.size=10) +
  labs(x=NULL, title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        legend.title=element_text(size=9))

#####################################



### ClusterGVis plotting of the target genes of the sig.&Top TFs
#####################################
###

library(Seurat)
library(dplyr)

### Extract the expr and make the Expression DF
# obj
net_subset=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T_TargetGenes.rds")
all_genes_to_extract=unique(c(net_subset$source, net_subset$target))
all_genes_to_extract=all_genes_to_extract[!grepl("^RP[SL]",all_genes_to_extract)]
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T=subset(THEObj, Annot.rough=="CD8T cells")

# expr
pseudo_CD8T=AggregateExpression(CD8T, assays="RNA", return.seurat=T, group.by=c("Annot.detailed","age"))
pseudo_CD8T$Annot.detailed=pseudo_CD8T$orig.ident
pseudo_CD8T$age=as.numeric(gsub(".*_","",rownames(pseudo_CD8T[[]])))
Genes_Expr=FetchData(pseudo_CD8T,
                     vars=c("Annot.detailed","age",
                            all_genes_to_extract), 
                     layer="data")

# expr-AUC (time-series)
library(pracma)
celltypes=unique(Genes_Expr$Annot.detailed)
genes=colnames(Genes_Expr)[3:ncol(Genes_Expr)]
THE_AUC_per_Gene=list()
for (i in 1:length(genes)) {
  subset_gene_=Genes_Expr[,c(1,2,2+i)]
  the_AUC_=c()
  for (j in 1:length(celltypes)) {
    subset_=subset(subset_gene_, Annot.detailed==celltypes[j])
    the_AUC_[j]=trapz(subset_$age, subset_[[genes[i]]])
  }
  THE_AUC_per_Gene[[i]]=the_AUC_
}
names(THE_AUC_per_Gene)=genes
AUC_df=as.data.frame(THE_AUC_per_Gene)
rownames(AUC_df)=celltypes

# gene list
net_subset_exist=subset(net_subset, (source %in% colnames(AUC_df)) & (target %in% colnames(AUC_df)))
source_genes=unique(net_subset_exist$source)
gene_list=list()
for (i in 1:length(source_genes)) {
  act=net_subset_exist %>% subset(source==source_genes[i] & mor==1)
  reg=net_subset_exist %>% subset(source==source_genes[i] & mor==(-1))
  act_genes=c(source_genes[i], act$target); act_genes=act_genes[!duplicated(act_genes)]
  if (table(act_genes %in% genes))
  reg_genes=c(source_genes[i], reg$target); reg_genes=reg_genes[!duplicated(reg_genes)]
  genes_under_this_TF=list(act_genes, reg_genes)
  names(genes_under_this_TF)=paste0(source_genes[i], c("_act.","_reg."))
  gene_list=c(gene_list, genes_under_this_TF)
}

# cluster
cluster=unlist(lapply(1:length(gene_list), function(i) rep(names(gene_list)[i], length(gene_list[[i]]))))

# gene
gene=unlist(gene_list)
names(gene)=NULL

# expression df
Expression_DF=list()
for (i in 1:nrow(AUC_df)) {
  Expression_DF[[i]]=as.numeric(AUC_df[i,])[match(gene, colnames(AUC_df))]
}
names(Expression_DF)=rownames(AUC_df)
Expression_DF=as.data.frame(Expression_DF)
Expression_DF$gene=gene
Expression_DF$cluster=cluster


### Cluster, enrich, and plot for each TF
library(ClusterGVis)
CK=ENRICH=CK_PLOT=list()
for (i in 1:length(unique(Expression_DF$cluster))) {
  
  # analyze on each TF
  df_for_plot=subset(Expression_DF, cluster==unique(Expression_DF$cluster)[i])
  rownames(df_for_plot)=df_for_plot$gene
  df_for_plot=df_for_plot %>%
    select(-gene, -cluster) %>%
    subset(rowSums(.)!=0)
  
  if (nrow(df_for_plot)>=100) {
    
    # cluster
    ck=clusterData(exp=df_for_plot,
                   cluster.method="kmeans",
                   cluster.num=5)
    CK[[i]]=ck
    
    # enrich
    library(org.Hs.eg.db)
    library(clusterProfiler)
    enrich=compareCluster(split(ck$long.res$gene, ck$long.res$cluster_name), 
                          fun=enrichGO,
                          keyType="SYMBOL",
                          OrgDb="org.Hs.eg.db",
                          ont="BP",
                          pvalueCutoff=0.05)
    ENRICH[[i]]=enrich
    
    # process the enrichment results
    enrich_filter=gofilter(enrich, level=4)
    enrich_result_df=enrich_filter@compareClusterResult
    enrich_result_df=enrich_result_df[!duplicated(enrich_result_df$ID), ]
    rownames(enrich_result_df)=enrich_result_df$ID
    enrich_result_df=enrich_result_df[,c(1,3,6)]
    colnames(enrich_result_df)=c("group","Description","pvalue")
    enrich_result_df$group=paste0("C",gsub("cluster ","",gsub(" \\(.*","",enrich_result_df$group)))
    enrich_result_df=enrich_result_df %>% slice_min(order_by=pvalue, n=5, by=group) # take the top5 terms
    
    # plot
    CK_PLOT[[i]]=
      visCluster(object=ck,
                 plot.type="both",
                 column_names_rot=90,
                 annoTerm.data=enrich_result_df,
                 line.side="left",
                 show_row_dend=F,
                 add.sampleanno=F,
                 add.box=T,
                 add.line=T)
  } else {
    
    # plot
    col_fun=circlize::colorRamp2(c(-2,0,2), c("#08519C", "white", "#A50F15"))
    
    CK_PLOT[[i]]=
      Heatmap(as.matrix(df_for_plot),
              name="Z-score",
              cluster_columns=T,
              show_column_dend=F,
              show_column_names=T,
              column_names_side="top",
              cluster_rows=F,
              show_row_dend=F,
              show_row_names=F,
              col=col_fun)
  }
  
  message(paste0("------ Cluster: ", unique(Expression_DF$cluster)[i], " (", i, ") Done. ------"))
}


### Enrich all the genes in each TF
library(org.Hs.eg.db)
library(clusterProfiler)

# enrich
TF_ENRICH=list()
gene_split=split(Expression_DF$gene, Expression_DF$cluster)
for (i in 1:length(gene_split)) {
  gene_=gene_split[[i]]
  TF_enrich=enrichGO(gene_,
                     keyType="SYMBOL",
                     OrgDb="org.Hs.eg.db",
                     ont="BP",
                     pvalueCutoff=0.05,
                     minGSSize=1,
                     maxGSSize=500)
  TF_enrich_filter=gofilter(TF_enrich, level=4)
  TF_ENRICH[[i]]=TF_enrich_filter
}

# plot
TF_ENRICH_Plots=list()
for (i in 1:length(gene_split)) {
  TF_ENRICH_Plots[[i]]=
    dotplot(TF_ENRICH[[i]], showCategory=20) + 
    ggtitle(paste0("\n\n\n\n\n\n\n\n\n\n\n", names(gene_split)[i], " genes")) +
    theme(axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11),
          title=element_text(size=12))
}

### Plot the heatmap and the general enrichment for each TF
library(patchwork)
library(ComplexHeatmap)
pdf("~/Project_PBMCage/Plots/PBMCage_TF_CD8T.pdf", height=10, width=15)
for (i in 1:length(CK_PLOT)) {
  plot(CK_PLOT[[i]])
  plot( TF_ENRICH_Plots[[i]] + (ggplot()+theme_void()) + (ggplot()+theme_void()) )
}
dev.off()

#####################################



### Reanalyze the CD8T TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD8T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD8T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD8T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("CD8T\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_CD8T_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("BATF","RFX5","KLF14","TGIF2","ETV1","RUNX2","NKX6-1","EOMES","SPI1","FLI1","KLF5","TBX5","TGIF1","SMAD3","RUNX1",
             "GRHL2","ZBTB12","STAT4","SMAD2","SOX4","ETS1","GFY","PTF1A","NRF1","YY1","ELK1","IRF2","TFAP2A","TCF4","TAL1",
             "FOSL2","JUN","PRDM1","GATA3","BACH2","GABPA","ETS1","ELK4","ELF1","ELK1","NFYA","SP1","NRF1")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("naive","CTL","Tcm","Tem"))
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free") +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-10, segment.color="grey60", segment.size=0.2, min.segment.length=0.5, max.overlaps=5) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="CD8T cells") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_CD8T_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### DC TF ANALYSIS #####################################
###############################################################################################################

### Transcriptional factors of DCs
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)

THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T=subset(THEObj, Annot.rough=="DCs")
CD8T_ages=names(table(CD8T$age)[table(CD8T$age)>0])
obj_age=list()
for (i in 1:length(CD8T_ages)) {
  obj_age[[i]]=subset(CD8T, age==CD8T_ages[i])
  obj_age[[i]]=Seurat::NormalizeData(obj_age[[i]])
}
names(obj_age)=CD8T_ages

net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_DCs.rds")

#####################################



### Analyze the DC TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_DCs.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("DC\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("DC\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("DC\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("DC\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_DC_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_DC_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")

# mark the classical TFs
CD4_CD8_TF=c("ID2","IRF8","IRF4","ZEB2","KLF4","BATF3","IKZF1","MEF2C","SPI1","GFI1","STAT3","CEBPA","CEBPB","RELB","RELA","REL","IRF4","RBPJ","BCL6","TCF3","TCF4","TCF12","BCL11B","SPIB","STAT5A", # in general
             "TCF4","ZEB2","IRF8","IRF4", # pDC
             "ID2","IRF8","BATF3", # cDC1
             "ID2","ZEB2","IRF4","NOTCH2","KLF4", # cDC2
             "ZEB2","IRF4","KLF4", # AXL DC
             "MAFB","KLF4" # mDC
)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster=ifelse(grepl("AXL",cluster),"AXL_DC",cluster)) %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  group_by(cluster.source) %>%
  filter(p.adj==min(p.adj)) %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("cDC1","cDC2_2","pDC","AXL_DC"))

# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free") +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-2.5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="DCs") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_DC_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### Monocyte TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for Monocytes
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="Monocytes")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])
  
  print(paste0("MetaObject ", i, " Creation: Starts......"))
  
  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)
  
  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)
  
  if (max(table(age_subset$Annot.detailed))<10) { # to fulfill the requirment of min_cells>=10 for metacell combination
    obj_age[[i]]=age_subset
  } else {
    # run the metacell combination
    age_subset=MetacellsByGroups(
      age_subset,
      group.by=c("Annot.detailed","age"),
      reduction="umap",
      k=3,
      min_cells=10,
      max_shared=1,
      target_metacells=25,
      verbose=TRUE,
      slot="counts",
      ident.group="Annot.detailed")
    
    # get the metacell obj
    obj_age[[i]]=GetMetacellObject(age_subset)
    obj_age[[i]]=JoinLayers(obj_age[[i]])
    obj_age[[i]]=NormalizeData(obj_age[[i]])
    # Check the metacells
    print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
    print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))
  }
  
  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_Monocytes.rds")

#####################################



### Transcriptional factors of Monocytes
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_Monocytes.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  if ("layers" %in% slotNames(the_obj_age@assays$RNA)) {
    mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  } else {
    mat=as.matrix(the_obj_age@assays$RNA$data)
  }
  
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Monocytes.rds")

#####################################



### Reanalyze the Monocyte TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Monocytes.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Mono\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Mono\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Mono\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Mono\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Monocytes_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Monocytes_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("NR4A1","SPI1","CEBPA","CEBPB","FOS","JUN","MNDA","KLF2","FOXQ1")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("nonclassical","inter","classical"))
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free", ncol=2) +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="Monocytes") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_Monocytes_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### NK cell TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for NK cells
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="NK cells")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])
  
  print(paste0("MetaObject ", i, " Creation: Starts......"))
  
  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)
  
  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)
  
  # run the metacell combination
  age_subset=MetacellsByGroups(
    age_subset,
    group.by=c("Annot.detailed","age"),
    reduction="umap",
    k=3,
    min_cells=10,
    max_shared=1,
    target_metacells=25,
    verbose=TRUE,
    slot="counts",
    ident.group="Annot.detailed")
  
  # get the metacell obj
  obj_age[[i]]=GetMetacellObject(age_subset)
  obj_age[[i]]=JoinLayers(obj_age[[i]])
  obj_age[[i]]=NormalizeData(obj_age[[i]])
  # Check the metacells
  print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
  print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))
  
  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_NK.rds")

#####################################



### Transcriptional factors of NK cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_NK.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_NK.rds")

#####################################



### Reanalyze the NK TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_NK.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("NK\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("NK\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("NK\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("NK\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_NK_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_NK_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("NFIL3","TCF7","ETS1","STAT5A","ID2","TBX21","EOMES","TOX","PRDM1","ZEB2","GATA3","SMAD4","SMAD3","FOXO1",
             "STAT4","STAT1","STAT2","IRF9","ZBTB32","IRF8","RUNX1","RUNX2","RUNX3","KLF12","FLI1","ZHX2","BCL11B","KLF2","RFX7","KLF2",
             "TCF7","LEF1","BACH2", # NK CD56hi
             "PRDM1") # NK CD56dim
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("CD56hi","CD56dim","prolif"))
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free", ncol=2) +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-2.5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="NK cells") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_NK_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### OtherT TF ANALYSIS #####################################
###############################################################################################################

### Make metacells for OtherT cells
#####################################
###

###
library(Seurat)
library(dplyr)
library(hdWGCNA)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# Take subsets
CD8T=subset(THEObj, Annot.rough=="OtherT")
CD8T_age=unique(CD8T$age)
CD8T_age=sort(CD8T_age)

# Make metacells for each age
obj_age=list()
for (i in 1:length(CD8T_age)) {
  age_subset=subset(CD8T, age==CD8T_age[i])
  
  print(paste0("MetaObject ", i, " Creation: Starts......"))
  
  # remove the celltype with only 0 or 1 cell
  More_than_oneCell_Celltype=names(table((age_subset[[]]$Annot.detailed)))[!table(age_subset[[]]$Annot.detailed) %in% c(0,1)]
  age_subset=subset(age_subset, Annot.detailed %in% More_than_oneCell_Celltype)
  
  # take the features
  age_subset=SetupForWGCNA(
    age_subset,
    group.by="Annot.detailed",
    gene_select="fraction",
    fraction=0.2,
    wgcna_name="WGCNA")
  length(age_subset@misc$WGCNA$wgcna_genes)
  
  if (max(table(age_subset$Annot.detailed))<10) { # to fulfill the requirment of min_cells>=10 for metacell combination
    obj_age[[i]]=age_subset
  } else {
    # run the metacell combination
    age_subset=MetacellsByGroups(
      age_subset,
      group.by=c("Annot.detailed","age"),
      reduction="umap",
      k=3,
      min_cells=10,
      max_shared=1,
      target_metacells=25,
      verbose=TRUE,
      slot="counts",
      ident.group="Annot.detailed")
    
    # get the metacell obj
    obj_age[[i]]=GetMetacellObject(age_subset)
    obj_age[[i]]=JoinLayers(obj_age[[i]])
    obj_age[[i]]=NormalizeData(obj_age[[i]])
    # Check the metacells
    print(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age))
    print(sum(table(obj_age[[i]]$Annot.detailed, obj_age[[i]]$age)))
  }
  
  print(paste0("====== MetaObject ", i, " Creation: Done! ======"))
}
names(obj_age)=CD8T_age
saveRDS(obj_age, "~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_OtherT.rds")

#####################################



### Transcriptional factors of OtherT cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/CD8T_NK_detailed/MetaCellObj_OtherT.rds")
CD8T_ages=names(obj_age)
net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  if ("layers" %in% slotNames(the_obj_age@assays$RNA)) {
    mat=as.matrix(the_obj_age@assays$RNA@layers$data)
  } else {
    mat=as.matrix(the_obj_age@assays$RNA$data)
  }
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherT.rds")

#####################################



### Reanalyze the OtherT TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherT.rds")

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("OtherT\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("OtherT\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("OtherT\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("OtherT\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherT_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherT_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("TCF12","TCF3","NOTCH1","MYB","RUNX1","STAT5A","GATA3","BCL11B","EGR1","EGR2","EGR3","ID3","SOX13","RUNX3","NR4A1","NR4A2","NR4A3","ETV5","KLF2","RELB","HES1","ZBTB16", # gdT
             "IKZF2","RORC","EOMES","STAT1","RUNX3","PRDM1","GATA3","ZBTB16","TBX21", # MAIT
             "ZBTB16","NFKB1","NFKB2","RELA","RELB","REL","EGR1","EGR2","EGR3","RORC","RUNX1","GATA3","ZBTB7B","MYC","TCF12","TCF3","TCF4","TBX21","ID2") # NKT
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free", ncol=2) +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-2.5, segment.color="grey60", segment.size=0.2) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="OtherT") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_OtherT_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################




###############################################################################################################
##################################### Other TF ANALYSIS #####################################
###############################################################################################################

### Transcriptional factors of Other cells
#####################################
###

###
library(decoupleR)
library(tidyr)
library(dplyr)

THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T=subset(THEObj, Annot.rough=="Other cells")
CD8T_ages=names(table(CD8T$age)[table(CD8T$age)>0])
obj_age=list()
for (i in 1:length(CD8T_ages)) {
  obj_age[[i]]=subset(CD8T, age==CD8T_ages[i])
  obj_age[[i]]=Seurat::NormalizeData(obj_age[[i]])
}
names(obj_age)=CD8T_ages

net=get_collectri(organism="human", split_complexes=FALSE)

### Determine the TF for each age
DF_SCORE_PVALUE=list()
for (i in 1:length(CD8T_ages)) {
  age_=CD8T_ages[i]
  
  message(paste0("Analysis on ", age_, ": Starts......"))
  
  the_obj_age=obj_age[[i]]
  mat=as.matrix(the_obj_age@assays$RNA@data)
  rownames(mat)=rownames(the_obj_age)
  colnames(mat)=colnames(the_obj_age)
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  score=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="score") %>%
    textshape::column_to_rownames("source")
  
  pvalue=acts %>%
    pivot_wider(id_cols="source", names_from="condition", values_from="p_value") %>%
    textshape::column_to_rownames("source")
  
  df_score=t(score) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="score") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(score, na.rm=T))
  
  df_pvalue=t(pvalue) %>%
    as.data.frame() %>%
    mutate(cluster=the_obj_age$Annot.detailed) %>%
    pivot_longer(cols=-cluster, names_to="source", values_to="pvalue") %>%
    group_by(cluster, source) %>%
    summarize(mean=mean(pvalue, na.rm=T))
  
  df_score_pvalue=df_score
  df_score_pvalue$pvalue=df_pvalue$mean[match(df_score$source, df_pvalue$source)]
  colnames(df_score_pvalue)=c("cluster","source","score","pvalue")
  df_score_pvalue$age=age_
  
  DF_SCORE_PVALUE=c(DF_SCORE_PVALUE, list(df_score_pvalue))
  
  message(paste0("--- Analysis on ", age_, ": Done! ---"))
}

DF_CD8T=data.table::rbindlist(DF_SCORE_PVALUE)
saveRDS(DF_CD8T, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Othercells.rds")

#####################################



### Reanalyze the Other cell - except for progenitors TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Othercells.rds")
DF_CD8T=subset(DF_CD8T, !grepl("progenitor",cluster))

### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Other\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Other\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Other\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("\\..*","",gsub("Other\\.","",cluster))) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherCell.ExceptForProgenitors_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherCell.ExceptForProgenitors_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("TBX21","EOMES","GATA3","RORC","RORA","STAT4","STAT3","STAT5A","STAT6","BATF","IKZF1","RUNX3","MAF","BCL11B","ZBTB46","BACH2","IRF4","GFI1", # ILC
             "GATA1","GATA2","FLI1","GABP","KLF1","ZBTB7A","TAL1","ZFPM1","LDB1","BCL11A","SPI1","MYB", # Eryth
             "RUNX1","GATA1","STAT3","NFKB1","PPARG","FLI1","GFI1B","MECOM","ETV6","NFE2","IKZF5", # platelet
             "SPI1","CEBPA","STAT5A","GATA1","GATA2","GATA3","MITF","NOTCH2","JUN","JDP2","JUND","JUNB","FOS","FOSB","FOSL1","FOSL2","NFE2","NRF1","NRF2","USF2","ATF2","NFATC1","ATF3") # mast
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("Eryth","Plt","Mast","ILC"))

# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free") +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=2.5, segment.color="grey60", segment.size=0.2, max.overlaps=8) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-10, segment.color="grey60", segment.size=0.2, max.overlaps=5) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="Mast cells, ILCs, erythrocytes, and platelets") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_OtherCell.ExceptForProgenitors_interSpecific_TFs.pdf", height=4, width=5.5)
plot(plot_)
dev.off()

#####################################



### Reanalyze the progenitor TF, focusing on celltype-specific ones and their changes across age
#####################################
###

library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)

### Load
DF_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_Othercells.rds")
DF_CD8T=subset(DF_CD8T, grepl("progenitor",cluster))


### Statistical test to get the celltype-specific TFs
# test data
tempt=DF_CD8T %>%
  mutate(cluster=gsub("Other\\.","",cluster)) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  nest()
test_function=function(i) {
  tryCatch({
    tempt2=tempt %>%
      purrr::pluck("data", i) %>%
      rstatix::wilcox_test(., score~cluster, ref.group="all") %>%
      adjust_pvalue(method="bonferroni")
    tempt2$source=tempt$source[i]
    return(tempt2)
  }, error=function(msg) return(NA))
}
DF_CD8T.wtest=list()
for (i in 1:nrow(tempt)) {
  DF_CD8T.wtest[[i]]=test_function(i)
}
DF_CD8T.wtest=DF_CD8T.wtest[which(unlist(lapply(DF_CD8T.wtest, function(x) class(x)[1]))!="logical")] # remove NA
DF_CD8T.wtest=DF_CD8T.wtest %>% data.table::rbindlist(.) %>% rename(cluster=group2)
# log2FC data
DF_CD8T.logFC=DF_CD8T %>%
  mutate(cluster=gsub("Other\\.","",cluster)) %>% 
  subset(pvalue<0.05) %>%
  group_by(source) %>%
  mutate(log2FC=log2(score/mean(score, na.rm=TRUE))) %>%
  ungroup() %>% group_by(cluster, source) %>%
  summarize_at("log2FC", mean)
DF_CD8T.wtest_logFC=left_join(DF_CD8T.wtest, DF_CD8T.logFC, by=c("cluster","source"))
# mark the cell-specific TFs
DF_CD8T.cellspecific=DF_CD8T.wtest_logFC %>%
  subset(p.adj<0.05) %>%
  group_by(source) %>%
  slice_min(p.adj, n=1) # use p.adj instead of log2FC which is coarsely calculated
celltypes.group.list=paste0(DF_CD8T.cellspecific$cluster, ".", DF_CD8T.cellspecific$source)
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  mutate(cell_specific=ifelse(cluster.source %in% celltypes.group.list, "cell-specific", NA))
# mark the top50 TFs in each celltype (can be common ones)
top50_TFs=DF_CD8T %>%
  mutate(cluster=gsub("Other\\.","",cluster)) %>% 
  subset(pvalue<0.05) %>% 
  group_by(cluster, source) %>% summarize_at("score",mean) %>%
  ungroup() %>% group_by(cluster) %>% slice_max(score, n=50) %>% 
  ungroup() %>%
  mutate(cluster.source.top=paste0(cluster,".",source)) %>%
  select(cluster.source.top) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50=ifelse(cluster.source %in% top50_TFs, "top50", NA))
# mark the top50 pseudo-cell-specific TFs (higher than the average of all celltypes, but can be shared by several celltypes)
top50_shared=DF_CD8T.wtest_logFC %>%
  group_by(cluster) %>%
  slice_min(p.adj, n=50) %>%
  ungroup() %>%
  select(cluster.source) %>%
  tibble::deframe()
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(top50_HigherThanAvg=ifelse(cluster.source %in% top50_shared, "top50_HigherThanAvg", NA))
# select the columns needed
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  select(cluster, statistic, p.adj, source, log2FC, cluster.source, cell_specific, top50, top50_HigherThanAvg)

### Correlation analysis to get the age-related changes in each TF in each celltype
cluster_TF_list=DF_CD8T %>%
  mutate(cluster=gsub("Other\\.","",cluster)) %>% 
  subset(pvalue<0.05) %>%
  subset(source %in% DF_CD8T.wtest$source) %>%
  split(list(.$cluster, .$source))
cluster_TF_list=cluster_TF_list[lapply(cluster_TF_list, nrow)>=3]
cluster_TF_list=cluster_TF_list %>% 
  lapply(., function(df) cor.test(df[["score"]], as.numeric(df[["age"]]), method="spearman", exact=F))
cluster_TF.cor_result=data.frame(cluster.source=names(cluster_TF_list),
                                 rho=unlist(lapply(cluster_TF_list, function(x) x$estimate)) %>% unname(),
                                 rho_pval=unlist(lapply(cluster_TF_list, function(x) x$p.value)) %>% unname()) %>%
  mutate(source=gsub(".*\\.","",cluster.source))
cluster_TF.tfs=paste0(paste0("\\.",unique(cluster_TF.cor_result$source)), collapse="|")
cluster_TF.cor_result=cluster_TF.cor_result %>%
  mutate(cluster=gsub(cluster_TF.tfs,"",cluster.source))

### Merge statistical analysis and correlation analysis
DF_CD8T.wtest_logFC_cor=DF_CD8T.wtest_logFC %>% 
  left_join(cluster_TF.cor_result, by=c("cluster.source","cluster","source")) %>%
  relocate(c("cluster.source","cluster","source","statistic","p.adj","log2FC","rho","rho_pval","top50","top50_HigherThanAvg","cell_specific"))

### Save the statistical analysis results and correlation analysis results
write.table(DF_CD8T.wtest_logFC_cor, "~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherCell.Progenitors_cellspecificTF_corAnalysis.txt", sep="\t")

### Plot cell-specific TFs
# get the celltype-specific TFs
DF_CD8T.wtest_logFC=read.delim("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_OtherCell.Progenitors_cellspecificTF_corAnalysis.txt", sep="\t")
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  subset(cell_specific=="cell-specific"|top50_HigherThanAvg=="top50_HigherThanAvg")
# mark the classical TFs
CD4_CD8_TF=c("GFI1","IKZF1","TCF3","MYC", # HSC
             "IKZF1","TCF3","SPI1","EBF1","PAX5","BCL11A","FOXO1","NOTCH1","NOTCH2","TCF12","GATA3","TCF7","TCF4", # CLP
             "TCF3","TAL1","GATA1", # MEP
             "TCF3","TAL1","GATA1","SPI1","IKZF1", # MPP
             "TBX21","NFIL3","ID2","TOX","EOMES","GATA3","RORA","BCL11B","RORC","ZNF683") # NKP
DF_CD8T.wtest_logFC=DF_CD8T.wtest_logFC %>%
  mutate(mark=ifelse(source %in% CD4_CD8_TF, source, NA)) %>% 
  mutate(mark_remain=ifelse(source %in% CD4_CD8_TF, NA, source)) %>%
  group_by(cluster) %>% arrange(desc(p.adj)) %>% tibble::rowid_to_column("id")
DF_CD8T.wtest_logFC$cluster=forcats::fct_relevel(DF_CD8T.wtest_logFC$cluster, c("progenitor.HSC","progenitor.MPP","progenitor.CLP","progenitor.NKP","progenitor.MEP"))

# plot
plot_=
  ggplot(DF_CD8T.wtest_logFC, aes(x=id, y=-log10(p.adj))) +
  facet_wrap(~cluster, scales="free", ncol=2) +
  ggrastr::geom_point_rast(shape=20, alpha=0.3, size=0.5, color="gray50") +
  ggrepel::geom_text_repel(aes(label=mark), size=3, nudge_x=-1, nudge_y=1.5, segment.color="grey60", segment.size=0.2) +
  ggrepel::geom_text_repel(aes(label=mark_remain), size=3, nudge_x=1, nudge_y=-5, segment.color="grey60", segment.size=0.2, max.overlaps=5) +
  theme_classic() +
  labs(x=expression(avg.log[2]~FC), y=expression(-log[10]~p.adj), title="Progenitors") +
  theme(axis.text.x=element_blank(),
        # axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        # axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        title=element_text(size=11, face="bold"),
        strip.background=element_rect(fill="gray90"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_OtherCell.Progenitors_interSpecific_TFs.pdf", height=6, width=5.5)
plot(plot_)
dev.off()

#####################################







### Age-correlation analysis of TFs in all celltypes
#####################################
###

### Load
celltypes_=c("Bcell","CD4T","CD8T","NK","OtherT","DC","Monocytes","OtherCell.Progenitors","OtherCell.ExceptForProgenitors")
cluster_DF=list()
for (i in 1:length(celltypes_)) {
  cluster_DF[[i]]=read.delim(paste0("~/Project_PBMCage/Tempt_RDS/TFanalysis_results_",celltypes_[i],"_cellspecificTF_corAnalysis.txt"), sep="\t")
}
names(cluster_DF)=celltypes_
# remove B.inter and B.pb
cluster_DF[["Bcell"]]=cluster_DF[["Bcell"]] %>% subset(cluster!="inter" & cluster!="pb")
# remove CD4T.Treg
cluster_DF[["CD4T"]]=cluster_DF[["CD4T"]] %>% subset(cluster!="Treg")
# merge AXL_mDC and AXL_pDC
cluster_DF[["DC"]]=cluster_DF[["DC"]] %>% 
  mutate(cluster=ifelse(grepl("AXL",cluster),"AXL_DC",cluster)) %>%
  mutate(cluster.source=paste0(cluster,".",source)) %>%
  group_by(cluster.source) %>%
  filter(p.adj==min(p.adj)) %>%
  filter((rho_pval<=0.05 & rho==max(rho)) | (rho_pval>0.05 & rho_pval==min(rho_pval))) %>%
  as.data.frame()

### Get top50 or top50_HigherThanAvg or cell_specific TFs in all Annot.inter celltype
union_TFs=cluster_DF %>% 
  lapply(., function(df) df %>% subset(!(is.na(top50) & is.na(top50_HigherThanAvg) & is.na(cell_specific)))) %>%
  Reduce(rbind, .) %>%
  select(source) %>% subset(!duplicated(.)) %>% tibble::deframe()

exist_TFs_inUnion=cluster_DF %>%
  lapply(., function(df) df %>% subset(source %in% union_TFs) %>% select(source) %>% subset(!duplicated(.)) %>% tibble::deframe()) %>%
  Reduce(intersect, .)

# if take only the key TFs published in each celltype
published_TFs=list(Bcell=c("ZEB2","XBP1","PRDM1","AICDA","IKZF3","PAX5","TBX21","TCF3","EBF1","STAT5A","BACH2","TCF4","IRF4","NOTCH2","HOXC4"),
                   CD4T=c("BATF","RFX5","KLF14","TGIF2","ETV1","RUNX2","NKX6-1","EOMES","SPI1","FLI1","KLF5","TBX5","TGIF1","SMAD3","RUNX1","GRHL2","ZBTB12","STAT4","SMAD2","SOX4","ETS1","GFY","PTF1A","NRF1","YY1","ELK1","IRF2","TFAP2A","TCF4","TAL1","FOSL2","JUN","PRDM1","GATA3","BACH2","GABPA","ETS1","ELK4","ELF1","ELK1","NFYA","SP1","NRF1"),
                   CD8T=c("BATF","RFX5","KLF14","TGIF2","ETV1","RUNX2","NKX6-1","EOMES","SPI1","FLI1","KLF5","TBX5","TGIF1","SMAD3","RUNX1","GRHL2","ZBTB12","STAT4","SMAD2","SOX4","ETS1","GFY","PTF1A","NRF1","YY1","ELK1","IRF2","TFAP2A","TCF4","TAL1","FOSL2","JUN","PRDM1","GATA3","BACH2","GABPA","ETS1","ELK4","ELF1","ELK1","NFYA","SP1","NRF1"),
                   DC=c("ID2","IRF8","IRF4","ZEB2","KLF4","BATF3","IKZF1","MEF2C","SPI1","GFI1","STAT3","CEBPA","CEBPB","RELB","RELA","REL","IRF4","RBPJ","BCL6","TCF3","TCF4","TCF12","BCL11B","SPIB","STAT5A","TCF4","ZEB2","IRF8","IRF4","ID2","IRF8","BATF3","ID2","ZEB2","IRF4","NOTCH2","KLF4","ZEB2","IRF4","KLF4","MAFB","KLF4"),
                   Mono=c("NR4A1","SPI1","CEBPA","CEBPB","FOS","JUN","MNDA","KLF2","FOXQ1"),
                   NK=c("NFIL3","TCF7","ETS1","STAT5A","ID2","TBX21","EOMES","TOX","PRDM1","ZEB2","GATA3","SMAD4","SMAD3","FOXO1","STAT4","STAT1","STAT2","IRF9","ZBTB32","IRF8","RUNX1","RUNX2","RUNX3","KLF12","FLI1","ZHX2","BCL11B","KLF2","RFX7","KLF2","TCF7","LEF1","BACH2","PRDM1"),
                   OtherT=c("TCF12","TCF3","NOTCH1","MYB","RUNX1","STAT5A","GATA3","BCL11B","EGR1","EGR2","EGR3","ID3","SOX13","RUNX3","NR4A1","NR4A2","NR4A3","ETV5","KLF2","RELB","HES1","ZBTB16","IKZF2","RORC","EOMES","STAT1","RUNX3","PRDM1","GATA3","ZBTB16","TBX21","ZBTB16","NFKB1","NFKB2","RELA","RELB","REL","EGR1","EGR2","EGR3","RORC","RUNX1","GATA3","ZBTB7B","MYC","TCF12","TCF3","TCF4","TBX21","ID2"),
                   OtherCell=c("TBX21","EOMES","GATA3","RORC","RORA","STAT4","STAT3","STAT5A","STAT6","BATF","IKZF1","RUNX3","MAF","BCL11B","ZBTB46","BACH2","IRF4","GFI1","GATA1","GATA2","FLI1","GABP","KLF1","ZBTB7A","TAL1","ZFPM1","LDB1","BCL11A","SPI1","MYB","RUNX1","GATA1","STAT3","NFKB1","PPARG","FLI1","GFI1B","MECOM","ETV6","NFE2","IKZF5","SPI1","CEBPA","STAT5A","GATA1","GATA2","GATA3","MITF","NOTCH2","JUN","JDP2","JUND","JUNB","FOS","FOSB","FOSL1","FOSL2","NFE2","NRF1","NRF2","USF2","ATF2","NFATC1","ATF3"),
                   HSPC=c("GFI1","IKZF1","TCF3","MYC","IKZF1","TCF3","SPI1","EBF1","PAX5","BCL11A","FOXO1","NOTCH1","NOTCH2","TCF12","GATA3","TCF7","TCF4","TCF3","TAL1","GATA1","TCF3","TAL1","GATA1","SPI1","IKZF1","TBX21","NFIL3","ID2","TOX","EOMES","GATA3","RORA","BCL11B","RORC","ZNF683")
                   )
published_TFs=unlist(published_TFs) %>% .[!duplicated(.)]
exist_TFs_inUnion=intersect(exist_TFs_inUnion, published_TFs)

### Take the rho data with the selected TFs
cluster_DF_cor=cluster_DF %>%
  lapply(., function(df) df %>% subset(source %in% exist_TFs_inUnion) %>% 
           mutate(rho=ifelse(rho_pval>0.05, NA, rho)) %>%
           select(rho, source, cluster) %>%
           tidyr::pivot_wider(names_from="cluster",values_from="rho") %>%
           tibble::column_to_rownames("source")) %>%
  lapply(., function(df) {df[is.na(df)]=0; t(df)})

# ### Get the corresponding rho_pval
# cluster_DF_cor_pval=cluster_DF %>%
#   lapply(., function(df) df %>% subset(source %in% exist_TFs_inUnion) %>% select(rho_pval, source, cluster) %>%
#            mutate(rho_pval=ifelse(rho_pval<0.01,"**",
#                                   ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", NA))) %>%
#            tidyr::pivot_wider(names_from="cluster",values_from="rho_pval") %>%
#            tibble::column_to_rownames("source")) %>%
#   lapply(., function(df) t(df))

### Plot individually
library(ComplexHeatmap)
col_fun=circlize::colorRamp2(c(-0.5, 0, 0.5), c("dodgerblue3", "white", "brown3"))

HT_list=list()
row_Titles=c("B cells","CD4T cells","CD8T cells","NK","Innate-like T","DCs","Mono","HSPC","Other cells")
for (i in 1:length(cluster_DF_cor)) {
  ht=
    Heatmap(cluster_DF_cor[[i]],
            name="mat", show_heatmap_legend=F,
            col=col_fun,
            rect_gp=gpar(col="grey50", lwd=0.1),
            cluster_rows=T, cluster_columns=T,
            show_row_dend=F, show_column_dend=F,
            column_names_rot=90, column_names_centered=F,
            column_names_gp=gpar(fontsize=8),
            row_names_gp=gpar(fontsize=9),
            row_title=stringr::str_wrap(row_Titles[i], width=9, whitespace_only=F),
            row_title_gp=gpar(fontsize=9, fontface="bold"),
            width=unit(4, "mm")*ncol(cluster_DF_cor[[i]]),
            height=unit(4, "mm")*nrow(cluster_DF_cor[[i]])
    )
  HT_list[[i]]=ht
}

lgd=Legend(col_fun=col_fun,
           title=NULL, title_gp=gpar(fontsize=11, fontface="bold"),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="lefttop", labels_gp=gpar(fontsize=8),
)

# codes=paste0(paste0("HT_list[[",1:length(HT_list),"]]"), collapse=" %v% ")
# eval(parse(text=paste0("ht_final=draw(",codes,", annotation_legend_list=lgd, annotation_legend_side='top')")))

ht_opt$ANNOTATION_LEGEND_PADDING=unit(4, "mm")
codes=paste0(paste0("HT_list[[",1:5,"]]"), collapse=" %v% ")
eval(parse(text=paste0("ht_final1=draw(",codes,")")))
codes=paste0(paste0("HT_list[[",6:9,"]]"), collapse=" %v% ")
eval(parse(text=paste0("ht_final2=draw(",codes,", annotation_legend_list=lgd, annotation_legend_side='top')")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_TF_AgeCor.Heatmap.pdf", height=4, width=7.5)
ht_final1
ht_final2
dev.off()

#####################################



###############################################################################################################
##################################### TF REANALYSIS BASED ON PSEUDOBULK_DETIALED #####################################
###############################################################################################################



# ### TF analysis on age-related and celltype-specific genes
# #####################################
# ###
# 
# ###
# library(decoupleR)
# library(tidyr)
# library(dplyr)
# 
# ### Load
# net=get_collectri(organism="human", split_complexes=FALSE)
# obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# marker_genes=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")
# corr_genes=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# 
# ### Run for each Annot.rough
# rough_celltypes=names(table(obj_age$Annot.rough))
# short_names_of_rough_celltypes=c("B","CD4T","CD8T","DC","Mono","NK","Other","OtherT")
# 
# for (i in 2:length(rough_celltypes)) {
#   obj_age_subset=obj_age %>% subset(Annot.rough==rough_celltypes[i])
#   # select genes are up/down-regulated with age and specifically expressed in this celltype
#   the_genes=marker_genes %>% subset(Annot.rough==rough_celltypes[i] & p_val_adj<0.05 & avg_log2FC>0) %>% select(gene) %>% tibble::deframe()
#   corr_genes_subset=corr_genes %>% subset(grepl(paste0(short_names_of_rough_celltypes[i],"\\."),celltypes) & analysis=="all" & rho_pval<0.05 & gene %in% the_genes)
#   
#   # run
#   data_df=Seurat::FetchData(obj_age_subset, vars=corr_genes_subset$gene, layer="data")
#   mat=as.matrix(t(data_df))
#   rownames(mat)=colnames(data_df)
#   colnames(mat)=rownames(data_df)
#   acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
#   
#   print(paste0("Analysis on ",rough_celltypes[i], " is done."))
#   
#   data.table::fwrite(acts, paste0("~/Project_PBMCage/Results/PBMC_results/TF_basedon_PseudobulkAnnotDetailed_",short_names_of_rough_celltypes[i],".csv.gz"), sep="\t")
# }
# 
# ### Combine the results
# short_names_of_rough_celltypes=c("B","CD4T","CD8T","DC","Mono","NK","Other","OtherT")
# rough_celltypes=c("B cells","CD4T cells","CD8T cells","DCs","Monocytes","NK cells","Other cells","OtherT")
# ACTs=list()
# for (i in 1:length(short_names_of_rough_celltypes)) {
#   ACTs[[i]]=data.table::fread(paste0("~/Project_PBMCage/Results/PBMC_results/TF_basedon_PseudobulkAnnotDetailed_",short_names_of_rough_celltypes[i],".csv.gz"), sep="\t")
# }
# ACTs=lapply(1:length(ACTs), function(idx) ACTs[[idx]] %>% mutate(celltypes=rough_celltypes[idx]))
# ACTs=data.table::rbindlist(ACTs)
# # save full results
# data.table::fwrite(ACTs, "~/Project_PBMCage/Results/PBMC_results/TF_basedon_PseudobulkAnnotDetailed_ALL.csv.gz", sep="\t")
# 
# #####################################



# ### TF analysis on enriched shared negatively age-correlated genes, i.e., rb-related genes
# #####################################
# ###
# 
# library(decoupleR)
# library(tidyr)
# library(dplyr)
# 
# ### Load
# net=get_collectri(organism="human", split_complexes=FALSE)
# obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# 
# # extract the enriched ribosome-related age-cor. genes
# ego_cluster=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
# ego_results=ego_cluster@compareClusterResult %>% 
#   # subset(p.adjust<0.05) %>%
#   subset(grepl("ribo|rRNA",Description))
# ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
# names(ego_genes)=ego_results$Description
# rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]
# 
# ### Run
# print("Loading...")
# data_df=Seurat::FetchData(obj_age, vars=rb_related_genes, layer="data")
# mat=as.matrix(t(data_df))
# rownames(mat)=colnames(data_df)
# colnames(mat)=rownames(data_df)
# 
# print("Start running...")
# acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
# 
# print("Done.")
# saveRDS(acts, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_ie.RbRelatedGenes_TF.rds")
# 
# ### Load the result
# acts=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_ie.RbRelatedGenes_TF.rds")
# # check the TFs with high scores
# acts_subset=acts %>%
#   subset(p_value<0.05) %>%
#   group_by(source) %>%
#   summarize_at("score", mean) %>%
#   arrange(desc(score)) %>%
#   subset(score>0)
# acts_subset
# 
### Process the pseudobulk obj for plotting
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
obj_age=obj_age %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()
# umap
obj_age_subset=obj_age %>%
  FindNeighbors(dims=1:10) %>%
  RunUMAP(dims=1:10)
DimPlot(obj_age_subset, reduction="umap", label=TRUE, pt.size=0.1, group.by="Annot.inter") + NoLegend() # check
saveRDS(obj_age_subset, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
# plot the umap
p1=
  DimPlot(obj_age_subset, reduction="umap", label=TRUE, pt.size=3, group.by="Annot.inter", label.size=2.8, repel=T, raster=T) +
  NoLegend() +
  scale_color_manual(values=paletteer::paletteer_c("grDevices::Pastel 1", 34)) +
  labs(title=NULL) +
  theme(axis.line=element_line(linewidth=0.25),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        title=element_blank())
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudobulk_detailed_umap.pdf", height=4.5, width=4.5)
p1
dev.off()
# 
# ### Add TF data to the obj
# obj_age_subset[['tfsulm']]=acts %>%
#   tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
#   tibble::column_to_rownames('source') %>%
#   Seurat::CreateAssayObject(.)
# # check the top TFs
# DefaultAssay(obj_age_subset)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age_subset, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3') +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists)
# # select the TFs for plotting
# acts_subset$source
# idx_storngAndCoverall=c(1,3)
# idx_mildAndCoverall=c(24,21,16)
# idx_mildAndCoverless=c(20,8)
# 
# ## Plot
# plot_idx1=
#   FeaturePlot(obj_age_subset, features=acts_subset$source[1], raster=T) +
#   scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                          limits=c(0,8),
#                          breaks=c(0,2,4,6)) +
#   labs(title="ribosome", subtitle=acts_subset$source[1]) +
#   theme(axis.line=element_blank(),
#         axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         plot.title=element_text(size=12.5, hjust=0),
#         plot.subtitle=element_text(size=12, hjust=0.5),
#         legend.key.width=unit(0.4,"cm"),
#         legend.text=element_text(size=9))
# # plot_idx3=
#   FeaturePlot(obj_age_subset, features=acts_subset$source[3], raster=T) +
#   scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                          limits=c(0,8),
#                          breaks=c(0,2,4,6)) +
#   labs(title=acts_subset$source[3]) +
#   labs(title=" ", subtitle=acts_subset$source[3]) +
#   theme(axis.line=element_blank(),
#         axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         plot.title=element_text(size=12.5, hjust=0),
#         plot.subtitle=element_text(size=12, hjust=0.5),
#         legend.key.width=unit(0.4,"cm"),
#         legend.text=element_text(size=9))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_rb_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=3.5, width=4)
# cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
# dev.off()
# 
# #####################################



# ### TF analysis on enriched genes related to "mt.","protein proc.","tRNA","cytoplasmic tr.","RNA splicing","cell activ.","prolif./immunity"
# #####################################
# ###
# 
# library(decoupleR)
# library(tidyr)
# library(dplyr)
# 
# ### Load
# net=get_collectri(organism="human", split_complexes=FALSE)
# obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# 
# # extract the enriched genes
# Super_cluster=c("mitochon","proteo|ubiquit","tRNA","cytoplasmic translation","RNA splicing","activa|adhesion","prolif|immun")
# Fullnames=c("mt.","protein proc.","tRNA","cytoplasmic tr.","RNA splicing","cell activ.","prolif./immunity")
# ego_objs=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")
# ego_results=ego_objs %>% lapply(., function(obj) obj@compareClusterResult) %>%
#   data.table::rbindlist()
# ego_subset=lapply(1:length(Super_cluster), function(idx) ego_results %>% 
#                     subset(grepl(Super_cluster[idx],Description)) %>% 
#                     select(Cluster,Description,geneID) %>%
#                     mutate(Super_cluster=Fullnames[idx])) %>%
#   data.table::rbindlist(.) %>%
#   split(.$Super_cluster) %>%
#   lapply(., function(df) df %>% split(.$Cluster) %>%
#            lapply(., function(dff) dff %>% select(geneID) %>% tibble::deframe(.) %>% paste0(.,collapse="/")) %>%
#            lapply(., function(genelist) strsplit(genelist, split="/")[[1]] %>% .[!duplicated(.)])
#   ) %>%
#   lapply(., function(genelist) genelist %>% unlist() %>% unname() %>% .[!duplicated(.)]) %>%
#   .[match(Fullnames, names(.))]
# 
# ### Run
# short_names_for_saving=names(ego_subset); short_names_for_saving=gsub(" |\\.|/","",short_names_for_saving)
# for (i in 1:length(ego_subset)) {
#   print(paste0("======",i,"======"))
#   print("Loading...")
#   data_df=Seurat::FetchData(obj_age, vars=ego_subset[[i]], layer="data")
#   mat=as.matrix(t(data_df))
#   rownames(mat)=colnames(data_df)
#   colnames(mat)=rownames(data_df)
#   
#   print("Start running...")
#   acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
#   
#   print("Done.")
#   saveRDS(acts, paste0("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_OtherTermsRelated_TF_",short_names_for_saving[i],".rds"))
# }
# 
# #####################################



### TF analysis on 12 superclusters
#####################################
###

library(decoupleR)
library(tidyr)
library(dplyr)

### Load
net=get_collectri(organism="human", split_complexes=FALSE)
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")

# extract the enriched genes
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters.rds")
exist_genes_merged=lapply(exist_genes_per_orientation, 
                          function(x) lapply(x, function(y) 
                            lapply(1:length(y), 
                                   function(idx) y[[idx]] %>%
                                     subset(abs(rho)>0.2) %>%
                                     mutate(superclusters=gsub("\\|.*","",names(y)[idx]))
                                   ) %>% 
                              data.table::rbindlist()) %>% 
                            data.table::rbindlist(.)) %>%
  data.table::rbindlist(.) %>%
  split(.$superclusters) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe() %>% .[!duplicated(.)])

### Run
short_names_for_saving=gsub(" |\\.|/|-","",names(exist_genes_merged))
for (i in 1:length(exist_genes_merged)) {
  print(paste0("======",i,"======"))
  print("Loading...")
  data_df=Seurat::FetchData(obj_age, vars=exist_genes_merged[[i]], layer="data")
  mat=as.matrix(t(data_df))
  rownames(mat)=colnames(data_df)
  colnames(mat)=rownames(data_df)
  
  print("Start running...")
  acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)
  
  print("Done.")
  saveRDS(acts, paste0("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_",short_names_for_saving[i],"_updated.rds"))
}

#####################################



### Plot the TF results on enriched genes related to "nucl. metab."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_nuclmetab_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>1), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be NR2F1 which is the highest apart from MYC

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene %in% c("NR2F1"))
cor_data_subset %>% subset(gene %in% c("VDR"))

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="NR2F1", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         breaks=c(0,2,4),
                         labels=c(0,2,4)
  ) +
  labs(
    title=NULL, # title="nucl. metab.",
    subtitle="NR2F1"
    ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_nuclmetab_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "anti-virus"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)
library(ggplot2)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_antivirus_updated.rds")

### Check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough, aes(x=source, y=score, color=Annot.rough)) +
  geom_point()
# ... turns out to be TP53 which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="TP53")
# ... turns out that TP53 is decreasing with age

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3') +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx2=
  FeaturePlot(obj_age, features="TP53", reduction="umap", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3'
                         # limits=c(-2,4),
                         # breaks=c(-2,0,2,4)
                         ) +
  labs(
    title=NULL, # title="nucl. metab.",
    subtitle="TP53"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

# plot_idx3=
  FeaturePlot(obj_age_subset, features=acts_subset$source[3], raster=T) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         limits=c(-2,6),
                         breaks=c(-2,0,2,4)) +
  labs(title=" ", subtitle=acts_subset$source[3]) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_antivirus_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx2), align="hv", nrow=1)
dev.off()

# check the remaining sources with score<0
acts_subset_low=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score<0)
acts_subset_low
plotlists_low=list()
for (i in 1:length(acts_subset_low$source)) {
  gene_=acts_subset_low$source[i]
  plotlists_low[[i]]=
    FeaturePlot(obj_age_subset, features=gene_, raster=T) +
    scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3') +
    labs(title=paste0("i=",i,"; ",gene_)) +
    theme(axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          title=element_blank(),
          legend.key.width=unit(0.4,"cm"),
          legend.text=element_text(size=9))
}
cowplot::plot_grid(plotlist=plotlists_low[21:40])

#####################################



### Plot the TF results on enriched genes related to "autophagy"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_autophagy_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough, aes(x=source, y=score, color=Annot.rough)) +
  geom_point()
# ... turns out to be ESR1 which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="ESR1")
# ... turns out that ESR1 is stable across age

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx2=
  FeaturePlot(obj_age, features="ESR1", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,8),
                         breaks=c(0,4,8),
                         labels=c(0,4,8)
                         ) +
  labs(
    title=NULL, # title="nucl. metab.",
    subtitle="ESR1"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))
# plot_idx3=
#   FeaturePlot(obj_age_subset, features="XBP1", raster=T) +
#   scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                          limits=c(-2,6),
#                          breaks=c(-2,0,2,4)) +
#   labs(title=" ", subtitle="XBP1") +
#   theme(axis.line=element_blank(),
#         axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         plot.title=element_text(size=12.5, hjust=0),
#         plot.subtitle=element_text(size=12, hjust=0.5),
#         legend.key.width=unit(0.4,"cm"),
#         legend.text=element_text(size=9))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_autophagy_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx2), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "cell division"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_celldivision_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough, aes(x=source, y=score, color=Annot.rough)) +
  geom_point(data=. %>% subset(score>0)) +
  theme(axis.text.x=element_text(size=5, angle=90, hjust=0.5))
# ... turns out that CREM is the key apart from MYC, as compared to STAT6 (insig.), REL (low rho), and EWSR1 (insig.)
# ... based on the kinase results, choose E2F rather than ETV4
# ... although CREM is much higher, it is more related to immune responses and cytokine production through cAMP signaling,
#..instead of cell cycling or proliferation

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene %in% c("E2F2","E2F3","E2F4"))
cor_data_subset %>% subset(gene %in% c("CREM","ETV4"))

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
# prepare acts with the group of E2F
acts_group=subset(acts, source %in% c("E2F2","E2F3","E2F4")) %>%
  group_by(statistic, condition) %>%
  summarize_at(c("score","p_value"), mean) %>%
  mutate(source="E2F_all") %>%
  relocate(colnames(acts))
acts_join=data.table::rbindlist(list(acts, acts_group))
obj_age[['tfsulm']]=acts_join %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="E2F-all", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,3),
                         breaks=c(0,1,2),
                         labels=c(0,1,2)
                         ) +
  labs(
    title=NULL, # title="nucl. metab.",
    subtitle="E2F"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_celldivision_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "cell response"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_cellresponse_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough %>% subset(score>2), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be HMGA2 which is the highest in all celltypes apart from MYC, ZBTB4

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="HMGA2")
# ... turns out that HMGA2 is stable across age

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="HMGA2", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,2),
                         # breaks=c(-2,0,2,4)
  ) +
  labs(
    title=NULL, # title="cell response",
    subtitle="HMGA2"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_cellresponse_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "cellular org."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_cellularorg_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>1), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be AR which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene %in% c("RFXAP","RFXANK","RFX5","CIITA"))
# ... turns out that CIITA/RFX is the key

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
# prepare acts with the group of CIITA/RFX
acts_group=subset(acts, source %in% c("RFXAP","RFXANK","RFX5","CIITA")) %>%
  group_by(statistic, condition) %>%
  summarize_at(c("score","p_value"), mean) %>%
  mutate(source="CIITA_PFX") %>%
  relocate(colnames(acts))
acts_join=data.table::rbindlist(list(acts, acts_group))
obj_age[['tfsulm']]=acts_join %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="CIITA-PFX", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         breaks=c(0,4,8)
  ) +
  labs(
    title=NULL, # title="cellular org.",
    subtitle="CIITA/RFX"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_cellularorg_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "energymetab"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_energymetab_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough, aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be AR which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="AR")
# ... turns out that AR is stable across age

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="AR", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         # breaks=c(-2,0,2,4)
                         ) +
  labs(
    title=NULL, # title="energy metab.",
    subtitle="AR"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_energymetab_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "immune proc."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_immuneproc_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>1), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be NFKB which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene %in% c("NFKB1","NFKB2"))

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_AfterPCA.rds")
# prepare acts with the group of CIITA/RFX
acts_group=subset(acts, source %in% c("NFKB1","NFKB2")) %>%
  group_by(statistic, condition) %>%
  summarize_at(c("score","p_value"), mean) %>%
  mutate(source="NFKB_both") %>%
  relocate(colnames(acts))
acts_join=data.table::rbindlist(list(acts, acts_group))
obj_age[['tfsulm']]=acts_join %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="NFKB-both", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         # breaks=c(-2,0,2,4)
  ) +
  labs(
    title=NULL, # title="immune proc.",
    subtitle="NFKB"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_immuneproc_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "leuk. devel."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_leukdevel_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>1), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out that KMT2A is the key apart from MYC and HBP1 (an inhibitor of MYC)

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="KMT2A")

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
FeaturePlot(obj_age, features="KMT2A", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,2),
                         # breaks=c(-2,0,2,4)
  ) +
  labs(
    title=NULL, # title="leuk. devel.",
    subtitle="KMT2A"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_leukdevel_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "prog. death"
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)
library(ggplot2)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_progdeath_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset
# ATF6 is the highest

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>0), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5))
# ... turns out to be GABPA which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="GABPA")
# ... turns out that GABPA is the key

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="GABPA", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         breaks=c(0,2,4),
                         labels=c(0,2,4)
  ) +
  labs(
    title=NULL, # title="prog. death",
    subtitle="GABPA"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_progdeath_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "prot. metab."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_protmetab_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(subset(acts_arranged_rough, score>1), aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out that ZBTB4 is the highest TF in all the celltypes apart from MYC

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="ZBTB4")

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="ZBTB4", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         # breaks=c(-2,0,2,4)
  ) +
  labs(
    title=NULL, # title="prot. metab.",
    subtitle="ZBTB4"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_protmetab_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Plot the TF results on enriched genes related to "ribosome syn."
#####################################
###

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

library(Seurat)
library(dplyr)

### Load the result
acts=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_ribosomesyn_updated.rds")
# check the TFs with high scores
acts_subset=acts %>%
  subset(p_value<0.05) %>%
  group_by(source) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
acts_subset

### Determine the cell-specific TFs for THE supercluster
# add Annot.detailed, inter, and rough info
acts_arranged=acts %>%
  mutate(Annot.detailed=gsub("_[0-9][0-9]{0,3}-[0-9].*","",condition))
celltypes_detailed=unique(acts_arranged$Annot.detailed)
acts_arranged=acts_arranged %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# summarize scores at the 3 levels
acts_arranged_detailed=acts_arranged %>%
  group_by(source, Annot.detailed) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_inter=acts_arranged %>%
  group_by(source, Annot.inter) %>%
  summarize_at(c("score","p_value"), mean)
acts_arranged_rough=acts_arranged %>%
  group_by(source, Annot.rough) %>%
  summarize_at(c("score","p_value"), mean)
ggplot(acts_arranged_rough, aes(x=source, y=score, color=Annot.rough)) +
  geom_point() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0.5))
# ... turns out to be MYC which is the highest in all celltypes

### Check the age-association of the TFs
cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
cor_data_subset=cor_data %>%
  subset(analysis=="all") %>%
  subset(gene %in% acts_arranged_rough$source)
cor_data_subset %>% subset(gene=="MYC")
# ... turns out that MYC is the key

### Add TF data to the obj
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated_AfterPCA.rds")
obj_age[['tfsulm']]=acts %>%
  tidyr::pivot_wider(id_cols='source', names_from='condition', values_from='score') %>%
  tibble::column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
# check the top TFs
DefaultAssay(obj_age)="tfsulm"
# plotlists=list()
# for (i in 1:length(acts_subset$source)) {
#   gene_=acts_subset$source[i]
#   plotlists[[i]]=
#     FeaturePlot(obj_age, features=gene_, raster=T) +
#     scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
#                            limits=c(-2,6),
#                            breaks=c(-2,0,2,4)) +
#     labs(title=paste0("i=",i,"; ",gene_)) +
#     theme(axis.title.x=element_text(size=10),
#           axis.title.y=element_text(size=10),
#           axis.text.x=element_text(size=9),
#           axis.text.y=element_text(size=9),
#           title=element_blank(),
#           legend.key.width=unit(0.4,"cm"),
#           legend.text=element_text(size=9))
# }
# cowplot::plot_grid(plotlist=plotlists[1:20])

### Plot
plot_idx1=
  FeaturePlot(obj_age, features="MYC", raster=T, pt.size=3) +
  scale_colour_gradient2(low='dodgerblue3', mid='white', high='brown3',
                         # limits=c(-2,6),
                         # breaks=c(-2,0,2,4)
  ) +
  labs(
    title=NULL, # title="ribosome syn.",
    subtitle="MYC"
  ) +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        plot.title=element_text(size=12, hjust=0),
        plot.subtitle=element_text(size=12, hjust=0.5),
        legend.key.width=unit(0.2,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.text=element_text(size=11),
        legend.box.spacing=unit(0,"mm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_ribosomesyn_related_MappingTo.Pseudobulk_detailed_umap.pdf", height=1.5, width=1.5)
cowplot::plot_grid(plotlist=list(plot_idx1), align="hv", nrow=1)
dev.off()

#####################################



### Manhattan Plot of TF results in each subprocesses
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(c("TP53"),c("ESR1","MYC"),c("E2F2","E2F3","E2F4"),c("HMGA2"),
                c("RFXAP","RFXANK","RFX5","CIITA"),c("AR"),c("NFKB1","NFKB2"),c("KMT2A"),
                c("NR2F1"),c("GABPA"),c("ZBTB4"),c("MYC","TP53"))

for (i in 1:length(short_names)) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_",short_names[i],"_updated.rds"))
  # add Annot.detailed, inte, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  # summarize the scores at detailed
  acts_perdetailed=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source, Annot.detailed) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0)
  # get the TFs with top average scores
  acts_top=acts_perdetailed %>%
    group_by(source) %>%
    summarize_at("score", max) %>%
    subset(score>0) %>%
    ungroup() %>%
    filter(score>=1*mean(score)) %>%
    arrange(desc(score)) %>%
    slice_max(score, n=5)
  # plot
  acts_perdetailed_celltypes=names(table(acts_perdetailed$Annot.detailed))
  acts_perdetailed_arrange=acts_perdetailed
  acts_perdetailed_arrange$Annot.detailed=forcats::fct_relevel(acts_perdetailed_arrange$Annot.detailed, acts_perdetailed_celltypes)
  
  plot_=
    ggplot(acts_perdetailed_arrange, aes(x=source, y=score, fill=Annot.detailed, label=source)) +
    geom_bar(stat="identity", position=position_dodge(width=0.8, preserve="single"), width=0.8, linewidth=0) +
    ggrepel::geom_text_repel(data=. %>%
                               group_by(source) %>%
                               filter(score==max(score) & (source %in% c(acts_top$source,chosen_TFs[[i]]))) %>%
                               tidyr::complete(Annot.detailed, fill=list(score=0, source="")),
                             position=position_dodge(width=0.8, preserve="single"), 
                             size=3, force=20,
                             box.padding=0.1, point.padding=0.1,
                             min.segment.length=0) +
    scale_fill_manual(values=paletteer::paletteer_c("ggthemes::Orange-Gold", length(acts_perdetailed_celltypes))) +
    labs(x="Transcription factor", y="Score", title=Fullnames[i]) +
    theme_classic() +
    theme(legend.position="none",
          plot.title=element_text(size=12),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10))
  
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_",short_names[i],"_related_MappingTo.Pseudobulk_detailed_Manhattan.pdf"), 
      height=2, width=4)
  plot(plot_)
  dev.off()
  
  Sys.sleep(3)
}

#####################################



### Plot the score ranking of TF for the 12 superclusters by score*freq
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

ACTs=list()
for (i in 1:length(short_names)) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_",short_names[i],"_updated.rds"))
  # add Annot.detailed, inter, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  # summarize the scores per TF
  acts_perTF=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0) %>%
    mutate(supercluster=Fullnames[i])
  ACTs[[i]]=acts_perTF
}
ACTs_all=data.table::rbindlist(ACTs)

### Plot the rank of source by weighted.score=sum(score*freq)
# add the times of occurance of each source in 12 supercluster as another indicator of importance
freq_source=ACTs_all %>% group_by(source) %>% count()
ACTs_all=ACTs_all %>% left_join(., freq_source, by="source") %>%
  mutate(score.freq=score*n/12)
# order by sum(score*freq) and take only the top15 TFs
ACTs_ordered=ACTs_all %>% group_by(source) %>% 
  summarize_at("score.freq", sum) %>%
  slice_max(score.freq, n=15)
ACTs_ordered$source=forcats::fct_relevel(ACTs_ordered$source, ACTs_ordered$source)
# plot
plot_=
  ggplot(ACTs_ordered, aes(x=source, y=score.freq)) +
    geom_point() +
    labs(x=NULL, y="Weighted score") +
    theme_classic() +
    theme(legend.position="none",
          plot.title=element_text(size=12),
          axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10)) +
  scale_x_discrete(guide=guide_axis(n.dodge=1))

pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_12superclusters_MappingTo.Pseudobulk_detailed_Ranking.pdf"), 
    height=1.5, width=3)
plot(plot_)
dev.off()

#####################################



### Plot the score ranking of TF for the immune proc. supercluster by score (filtered with occurrence)
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")

ACTs=list()
for (i in 7) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_",short_names[i],"_updated.rds"))
  # add Annot.detailed, inter, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  
  # summarize the scores per TF
  acts_perTF=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source, Annot.detailed) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0)
  
  # add occurrence for each source on cell types
  celltypecount_all=length(unique(acts_perTF$Annot.detailed)) # count the types of detailed celltype
  acts_occurrence_type=acts_perTF %>%
    ungroup() %>%
    group_by(source) %>%
    count() %>%
    rename(celltypen=n)
  
  # add occurrence for each source on cells
  cellcount_all=length(unique(acts_rearrange$condition)) # count the total cells
  acts_occurrence_cell=acts_rearrange %>%
    ungroup() %>%
    group_by(source) %>%
    count() %>%
    rename(celln=n)
  
  acts_perTF.occur=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source) %>%
    summarize_at("score", mean) %>%
    subset(score>0) %>% 
    left_join(., acts_occurrence_type, by="source") %>%
    left_join(., acts_occurrence_cell, by="source") %>%
    arrange(desc(score)) %>%
    subset(celln>cellcount_all*0.5) # filter the sources with occurrence in most of the detailed celltypes
}
### Plot the rank of TFs for immune proc.
acts_top_TF=acts_perTF.occur %>% 
  slice_max(score, n=15)
acts_top_TF$source=forcats::fct_relevel(
  acts_top_TF$source, acts_top_TF %>%
    dplyr::select(source) %>%
    tibble::deframe())
plot_=
  ggplot(acts_top_TF %>% arrange(desc(score)), aes(x=source, y=score)) +
  geom_point() +
  # geom_segment(aes(x="NFKB2", xend="NFKB2", y=3.3, yend=3.1), linewidth=0.1, arrow=arrow(length=unit(0.04, "npc")), color="coral4") +
  labs(x=NULL, y="Average score", title="immune proc.") +
  theme_classic() +
  theme(legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))
pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_immuneproc_related_MappingTo.Pseudobulk_detailed_Ranking.pdf"), 
    height=1.7, width=3)
plot(plot_)
dev.off()

#####################################



### Plot the alterations of the key TFs with age at Annot.rough
#####################################
### 

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Extract the key TFs
acts_list=list()

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB=c("NFKB1","NFKB2"),KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
for (i in 1:length(short_names)) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms_TF_",short_names[i],".rds"))
  # add Annot.detailed, inte, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  # summarize the scores at detailed
  acts_perdetailed=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source, Annot.detailed) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0)
  # get the TFs with top average scores
  acts_top=acts_perdetailed %>%
    group_by(source) %>%
    summarize_at("score", max) %>%
    subset(score>0) %>%
    ungroup() %>%
    filter(score>=1*mean(score)) %>%
    arrange(desc(score)) %>%
    slice_max(score, n=5)
  
  acts_list[[i]]=c(chosen_TFs[[i]], acts_top$source)
}
acts_list_total=Reduce(c, acts_list) %>% .[!duplicated(.)]

### Make the df
TF_expr=FetchData(THEObj, vars=c("donor_id","age","sex","Annot.detailed","Annot.inter","Annot.rough",
                                 acts_list_total), layer="data")
TF_expr=TF_expr %>% tibble::rownames_to_column("cell_id")
data.table::fwrite(TF_expr, "~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")



### Plot all the 8 Annot.rough celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean)

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot only CD4T and CD8T
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)

SCORE_DFs=SCORE_DFs[c("CD4T cells","CD8T cells","NK cells","OtherT","B cells","DCs","Monocytes","Other cells")]
ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.rough) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)

  ht_list[[i]]=
    Heatmap(score_df_t, na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Z-score",
              legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            left_annotation=
              rowAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                                col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                                annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                            show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=T, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=F,
            row_names_gp=gpar(fontsize=9),
            row_title=names(SCORE_DFs)[i],
            row_title_gp=gpar(fontsize=10), row_title_side="left", row_title_rot=90,
            column_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(1,"mm")*nrow(score_df_t),
            width=unit(5,"mm")*ncol(score_df_t))
}
ht_list_=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("grey50","white")),
               at=1:length(age_shared), 
               labels=rev(rownames(score_df_t)[1:length(age_shared)]),
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, heatmap_legend_side="top", merge_legend=FALSE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.WholeObj_Ucell_heatmap.pdf", height=7.2, width=3.8)
draw(whole_ht)
dev.off()



### Plot all the 8 Annot.rough celltypes for NFKB and MYC only
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(agecut, Annot.rough, MYC, NFKB2)

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot only CD4T and CD8T
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)

SCORE_DFs=SCORE_DFs[c("CD4T cells","CD8T cells","NK cells","OtherT","B cells","DCs","Monocytes","Other cells")]
ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.rough) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(t(score_df_t), na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Yazar et al. (2022)\nZ-score",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            top_annotation=
              HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                            col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                            annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                            show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=c(rep(F,(length(SCORE_DFs)-1)),T)[i], 
            show_column_names=F,
            column_names_gp=gpar(fontsize=9),
            column_title=names(SCORE_DFs)[i],
            column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=90,
            row_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(5,"mm")*nrow(t(score_df_t)),
            width=unit(0.75,"mm")*ncol(t(score_df_t)))
}
ht_list_=Reduce(`+`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared), 
               labels=rownames(score_df_t)[1:length(age_shared)], labels_rot=90,
               legend_width=unit(7, "cm"), grid_height=unit(0.2, "cm"),
               title="age", title_position="topcenter", title_gp=gpar(fontface="plain", fontsize=10),
               direction="horizontal", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$legend_gap=unit(2,"cm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, annotation_legend_side="top", 
              heatmap_legend_side="top", merge_legend=TRUE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.WholeObj_Ucell_heatmap_NFKBandMYC.pdf", height=2.1, width=5.5)
draw(whole_ht)
dev.off()



### Plot all the 8 Annot.rough celltypes for other TFs than NFKB or MYC
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),
                # NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(-all_of(c("MYC", "NFKB2")))

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot only CD4T and CD8T
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)

SCORE_DFs=SCORE_DFs[c("CD4T cells","CD8T cells","NK cells","OtherT","B cells","DCs","Monocytes","Other cells")]
ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.rough) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(t(score_df_t), na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Z-score",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            top_annotation=
              HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                                col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                                annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                                show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=c(rep(F,(length(SCORE_DFs)-1)),T)[i], 
            show_column_names=F,
            column_names_gp=gpar(fontsize=9),
            column_title=names(SCORE_DFs)[i],
            column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=90,
            row_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(4,"mm")*nrow(t(score_df_t)),
            width=unit(1,"mm")*ncol(t(score_df_t)))
}
ht_list_=Reduce(`+`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared), 
               labels=rownames(score_df_t)[1:length(age_shared)], labels_rot=90,
               legend_width=unit(8, "cm"), grid_height=unit(0.2, "cm"),
               title="age", title_position="topcenter", title_gp=gpar(fontface="plain", fontsize=10),
               direction="horizontal", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$legend_gap=unit(2,"cm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, annotation_legend_side="top", 
              heatmap_legend_side="top", merge_legend=TRUE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.WholeObj_Ucell_heatmap_Other.than.NFKBorMYC.pdf", height=4, width=7)
draw(whole_ht)
dev.off()

#####################################



### Plot the alterations of the key TFs with age at Annot.inter within CD4T and CD8T
#####################################
### 

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Extract the key TFs
acts_list=list()

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB2=c("NFKB2"),KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))

### Plot all the CD4T and CD8T Annot.inter celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
TF_expr=subset(TF_expr, Annot.rough %in% c("CD4T cells","CD8T cells"))
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.inter*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.inter) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean)
age_shared=TF_expr_sumByagecut %>% group_by(Annot.inter, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.inter)

ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.inter) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(score_df_t, na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Z-score",
              legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            left_annotation=
              rowAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                            col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                            annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                            show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=T, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=F,
            row_names_gp=gpar(fontsize=9),
            row_title=names(SCORE_DFs)[i],
            row_title_gp=gpar(fontsize=10), row_title_side="left", row_title_rot=0,
            column_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(1,"mm")*nrow(score_df_t),
            width=unit(10,"mm")*ncol(score_df_t))
}
ht_list_=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("grey50","white")),
               at=1:length(age_shared), 
               labels=rev(rownames(score_df_t)[1:length(age_shared)]),
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, heatmap_legend_side="top", merge_legend=FALSE, ht_gap=unit(0.2, "cm"))

pdf("~/tempt.pdf", height=30, width=10)
draw(whole_ht)
dev.off()

#####################################



### Plot the alterations of the key TFs with age at Annot.inter within NK and OtherT
#####################################
### 

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Extract the key TFs
acts_list=list()

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB2=c("NFKB2"),KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))

### Plot all the NK and OtherT Annot.inter celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
TF_expr=subset(TF_expr, Annot.rough %in% c("NK cells","OtherT"))
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.inter*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.inter) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean)
age_shared=TF_expr_sumByagecut %>% group_by(Annot.inter, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.inter)

ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.inter) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(score_df_t, na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Z-score",
              legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            left_annotation=
              rowAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                            col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                            annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                            show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=T, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=F,
            row_names_gp=gpar(fontsize=9),
            row_title=names(SCORE_DFs)[i],
            row_title_gp=gpar(fontsize=10), row_title_side="left", row_title_rot=0,
            column_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(1,"mm")*nrow(score_df_t),
            width=unit(10,"mm")*ncol(score_df_t))
}
ht_list_=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("grey50","white")),
               at=1:length(age_shared), 
               labels=rev(rownames(score_df_t)[1:length(age_shared)]),
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, heatmap_legend_side="top", merge_legend=FALSE, ht_gap=unit(0.2, "cm"))

pdf("~/tempt.pdf", height=30, width=10)
draw(whole_ht)
dev.off()

#####################################



### Plot the alterations of MYC only with age at all Annot.detailed within CD4T, CD8T, NK, and OtherT
#####################################
### 

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Extract the key TFs
acts_list=list()

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
# chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
#                 `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB2=c("NFKB2"),KMT2A=c("KMT2A"),
#                 NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
chosen_TFs=list(MYC=c("MYC"))

### Plot all the CD4T and CD8T Annot.inter celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
TF_expr=subset(TF_expr, Annot.rough %in% c("CD4T cells","CD8T cells","NK cells","OtherT"))
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.inter*agecut)
TF_expr_sumByagecut=TF_expr %>% group_by(Donor_id, Annot.detailed) %>% 
  summarize_at(colnames(.)[!grepl("agecut|Annot\\.|sex|cell_id", colnames(.))], mean)
age_shared=TF_expr_sumByagecut %>% group_by(Annot.detailed, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% 
  tidyr::pivot_wider(id_cols="agecut", names_from="Annot.detailed", values_from="MYC") %>%
  tibble::column_to_rownames("agecut")
SCORE_DFs=scale(SCORE_DFs)
# ht_=
  Heatmap(SCORE_DFs, na_col="transparent",
          name=" ", 
          # rect_gp=gpar(col="white", lwd=0.5),
          heatmap_legend_param=list(
            title="Z-score",
            legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ), 
          left_annotation=
            rowAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                          col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                          annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                          show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
          col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
          cluster_columns=F, cluster_rows=F, show_column_dend=F, show_row_dend=F,
          show_row_names=F,
          row_names_gp=gpar(fontsize=9),
          # row_title=names(SCORE_DFs)[i],
          row_title_gp=gpar(fontsize=10), row_title_side="left", row_title_rot=0,
          column_names_gp=gpar(fontsize=9),
          # column_title=names(SCORE_DFs)[i],
          # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
          height=unit(1,"mm")*nrow(SCORE_DFs),
          width=unit(5,"mm")*ncol(SCORE_DFs))


ht_list_=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("grey50","white")),
               at=1:length(age_shared), 
               labels=rev(rownames(score_df_t)[1:length(age_shared)]),
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, heatmap_legend_side="top", merge_legend=FALSE, ht_gap=unit(0.2, "cm"))

pdf("~/tempt.pdf", height=30, width=10)
draw(whole_ht)
dev.off()

#####################################



### Plot TFs across ages/ TFs cor with ages at Annot.rough/detailed
#####################################
###

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Extract the key TFs
acts_list=list()

Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
# chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
#                 `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),NFKB2=c("NFKB2"),KMT2A=c("KMT2A"),
#                 NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
chosen_TFs=list(NFKB2=c("NFKB2"), MYC=c("MYC"))

### Plot all the Annot.rough celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
TF_expr=TF_expr %>% 
  mutate(cell_group=ifelse(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells","OtherT"), 
                           "Effector lymphocytes", 
                           "Supportive cells"))
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.|cell_group",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

### Plot MYC/NFKB2 expr across ages at Annot.rough
TF_expr_sumByagecut=TF_expr %>% group_by(donor_id, age, sex, Annot.rough, cell_group) %>% 
  summarize_at(colnames(.)[!grepl("age|Annot\\.|sex|_id|cell_group", colnames(.))], mean)
MYC_Yazar=
  ggplot(TF_expr_sumByagecut, aes(x=age, y=MYC, color=Annot.rough)) +
  facet_wrap(~cell_group) +
  coord_cartesian(ylim=c(0,0.5)) +
  ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=Annot.rough), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  ggpubr::stat_cor(aes(group=Annot.rough), method="spearman", 
                   label.y=c(0.5,0.5,0.45,0.45,0.4,0.4,0.35,0.35)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="Yazar et al. (2022) dataset", y="Expr. of MYC", x="age") +
  theme_bw() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1)))
NFKB2_Yazar=
  ggplot(TF_expr_sumByagecut, aes(x=age, y=NFKB2, color=Annot.rough)) +
    facet_wrap(~cell_group) +
    coord_cartesian(ylim=c(0,0.3)) +
    ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
    geom_smooth(aes(group=Annot.rough), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
    ggpubr::stat_cor(aes(group=Annot.rough), method="spearman", 
                     label.y=c(0.3,0.3,0.275,0.275,0.25,0.25,0.225,0.225)) +
    scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
    labs(title="Yazar et al. (2022) dataset", y="Expr. of NFKB2", x="age") +
    theme_bw() +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.title.y=element_text(size=10),
          axis.title.x=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9)) +
    guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1)))

### Volcano Plot of rho vs. pval (correlation between MYC expr and age) at Annot.detailed
TF_expr_sum.detailed=TF_expr %>% 
  group_by(donor_id, age, sex, Annot.rough, Annot.inter, Annot.detailed, cell_group) %>% 
  summarize_at(colnames(.)[!grepl("age|Annot\\.|sex|_id|cell_group", colnames(.))], mean)
celltypes=unique(TF_expr_sum.detailed$Annot.detailed)
cor=pval=c()
for (i in 1:length(celltypes)) {
  celltype_=celltypes[[i]]
  tempt_=
    TF_expr_sum.detailed %>% subset(Annot.detailed==celltype_) %>% 
    dplyr::ungroup() %>%
    dplyr::select(age, MYC)
  res_=cor.test(tempt_$age, tempt_$MYC, method="spearman")
  cor=c(cor, res_$estimate %>% unname())
  pval=c(pval, res_$p.value %>% unname())
}
df_=data.frame(cor=cor, pval=pval, Annot.detailed=celltypes) %>% 
  mutate(Annot.rough=paste0(gsub("\\..*","",Annot.detailed), " cells")) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="DC cells","DCs",
                            ifelse(Annot.rough=="Mono cells","Monocytes",
                                   ifelse(Annot.rough=="OtherT cells","OtherT",Annot.rough))))
plot_Annot.detailed_MYC=
  ggplot(df_, aes(x=cor, y=-log10(pval), color=Annot.rough)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="MYC expression across age", y=expression(-log[10]~pval), x=expression(Spearman~rho)) +
  theme_classic() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1))) +
  ggrepel::geom_text_repel(data=. %>% group_by(Annot.rough) %>% subset(cor<0 & pval<0.05) %>% slice_min(pval, n=3),
                           aes(label=Annot.detailed), min.segment.length=0, box.padding=0.1, point.padding=0.1, 
                           max.overlaps=8, force=5, force_pull=0.1,
                           show.legend=F, size=3.5)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.corWage_pvalXcor_MYC_Annot.detailed.pdf", height=2.5, width=5)
plot(plot_Annot.detailed_MYC)
dev.off()

### Volcano Plot of rho vs. pval (correlation between NFKB2 expr and age) at Annot.detailed
celltypes=unique(TF_expr_sum.detailed$Annot.detailed)
cor=pval=c()
for (i in 1:length(celltypes)) {
  celltype_=celltypes[[i]]
  tempt_=
    TF_expr_sum.detailed %>% subset(Annot.detailed==celltype_) %>% 
    dplyr::ungroup() %>%
    dplyr::select(age, NFKB2)
  res_=cor.test(tempt_$age, tempt_$NFKB2, method="spearman")
  cor=c(cor, res_$estimate %>% unname())
  pval=c(pval, res_$p.value %>% unname())
}
df_=data.frame(cor=cor, pval=pval, Annot.detailed=celltypes) %>% 
  mutate(Annot.rough=paste0(gsub("\\..*","",Annot.detailed), " cells")) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="DC cells","DCs",
                            ifelse(Annot.rough=="Mono cells","Monocytes",
                                   ifelse(Annot.rough=="OtherT cells","OtherT",Annot.rough))))
plot_Annot.detailed_NFKB2=
  ggplot(df_, aes(x=cor, y=-log10(pval), color=Annot.rough)) +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="NFKB2 expression across age", y=expression(-log[10]~pval), x=expression(Spearman~rho)) +
  theme_classic() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text=element_text(size=10)) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1))) +
  ggrepel::geom_text_repel(data=. %>% group_by(Annot.rough) %>% subset(cor>0 & pval<0.05 & Annot.rough!="Other cells") %>% 
                             slice_min(pval, n=2),
                           aes(label=Annot.detailed), min.segment.length=0, box.padding=0.1, point.padding=0.1, 
                           # max.overlaps=8, force=3, force_pull=0.1,
                           show.legend=F, size=3.5)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.corWage_pvalXcor_NFKB2_Annot.detailed.pdf", height=2.5, width=5)
plot(plot_Annot.detailed_NFKB2)
dev.off()

### Scatter Plot of MYC/NFKB2 expr across ages at Annot.detailed
# MYC_Yazar_inter=
  ggplot(TF_expr_sum.detailed, aes(x=age, y=MYC, color=Annot.detailed)) +
  facet_wrap(~cell_group) +
  coord_cartesian(ylim=c(0,1)) +
  # ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=cell_group), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  geom_smooth(aes(group=Annot.detailed), method="lm", se=F, show.legend=F, linewidth=0.1) +
  ggpubr::stat_cor(aes(group=cell_group), method="spearman") +
  # scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="Yazar et al. (2022) dataset", y="Expr. of MYC", x="age") +
  theme_bw() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  # guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1))) +
  guides(color="none")
# NFKB2_Yazar=
ggplot(TF_expr_sum.inter %>% subset(cell_group=="Effector lymphocytes"), aes(x=age, y=NFKB2, color=Annot.inter)) +
  ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=Annot.inter), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  ggpubr::stat_cor(aes(group=Annot.inter), method="spearman") +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title="Yazar et al. (2022) dataset", y="Expr. of NFKB2", x="age") +
  theme_bw() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1)))

### Plot all the 8 Annot.rough celltypes for NFKB and MYC only
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(agecut, Annot.rough, MYC, NFKB2)

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot only CD4T and CD8T
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)

SCORE_DFs=SCORE_DFs[c("CD4T cells","CD8T cells","NK cells","OtherT","B cells","DCs","Monocytes","Other cells")]
ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.rough) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(t(score_df_t), na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Yazar et al. (2022)\nZ-score",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            top_annotation=
              HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                                col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                                annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                                show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=c(rep(F,(length(SCORE_DFs)-1)),T)[i], 
            show_column_names=F,
            column_names_gp=gpar(fontsize=9),
            column_title=names(SCORE_DFs)[i],
            column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=90,
            row_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(5,"mm")*nrow(t(score_df_t)),
            width=unit(0.75,"mm")*ncol(t(score_df_t)))
}
ht_list_=Reduce(`+`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared), 
               labels=rownames(score_df_t)[1:length(age_shared)], labels_rot=90,
               legend_width=unit(7, "cm"), grid_height=unit(0.2, "cm"),
               title="age", title_position="topcenter", title_gp=gpar(fontface="plain", fontsize=10),
               direction="horizontal", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$legend_gap=unit(2,"cm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, annotation_legend_side="top", 
              heatmap_legend_side="top", merge_legend=TRUE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.WholeObj_Ucell_heatmap_NFKBandMYC.pdf", height=2.1, width=5.5)
draw(whole_ht)
dev.off()



### Plot all the 8 Annot.rough celltypes for other TFs than NFKB or MYC
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),
                # NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
TF_expr=data.table::fread("~/Project_PBMCage/Results/PBMC_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|age|sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(-all_of(c("MYC", "NFKB2")))

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot only CD4T and CD8T
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)

SCORE_DFs=SCORE_DFs[c("CD4T cells","CD8T cells","NK cells","OtherT","B cells","DCs","Monocytes","Other cells")]
ht_list=list()
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% 
    select(-Annot.rough) %>%
    tibble::column_to_rownames("agecut") %>%
    .[age_shared, ]
  # score_df_t=data.table::transpose(score_df_up); colnames(score_df_t)=rownames(score_df_up); rownames(score_df_t)=colnames(score_df_up)
  score_df_t=score_df_up
  score_df_t=scale(score_df_t)
  # score_df_t=score_df_t-rowMeans(score_df_t)
  
  ht_list[[i]]=
    Heatmap(t(score_df_t), na_col="transparent",
            name=" ", 
            # rect_gp=gpar(col="white", lwd=0.5),
            heatmap_legend_param=list(
              title="Z-score",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            top_annotation=
              HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                                col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                                annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                                show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=c(rep(F,(length(SCORE_DFs)-1)),T)[i], 
            show_column_names=F,
            column_names_gp=gpar(fontsize=9),
            column_title=names(SCORE_DFs)[i],
            column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=90,
            row_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(4,"mm")*nrow(t(score_df_t)),
            width=unit(1,"mm")*ncol(t(score_df_t)))
}
ht_list_=Reduce(`+`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared), 
               labels=rownames(score_df_t)[1:length(age_shared)], labels_rot=90,
               legend_width=unit(8, "cm"), grid_height=unit(0.2, "cm"),
               title="age", title_position="topcenter", title_gp=gpar(fontface="plain", fontsize=10),
               direction="horizontal", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$legend_gap=unit(2,"cm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, annotation_legend_side="top", 
              heatmap_legend_side="top", merge_legend=TRUE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/TF_KeyTFs.WholeObj_Ucell_heatmap_Other.than.NFKBorMYC.pdf", height=4, width=7)
draw(whole_ht)
dev.off()

#####################################



###############################################################################################################
##################################### TF Analysis on WGCNA hemostasis module #####################################
###############################################################################################################


### TF analysis on WGCNA hemostasis module based on the WGCNA obj which was derived from Pseudobulk_rough
#####################################
###

library(decoupleR)
library(tidyr)
library(dplyr)

### Load
obj_age=readRDS("~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")
net=get_collectri(organism="human", split_complexes=FALSE)

# extract the module genes
Module_genes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>%
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()

### Run
print("Loading...")
data_df=Seurat::FetchData(obj_age, vars=Module_genes, layer="data")
mat=as.matrix(t(data_df))
rownames(mat)=colnames(data_df)
colnames(mat)=rownames(data_df)

print("Start running...")
acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)

print("Done.")
saveRDS(acts, paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.rough_HubGenes_TF_hemostasis.rds"))

#####################################



### TF analysis on WGCNA hemostasis module based on the Pseudobulk_detailed which is used for all other TF analysis
### ... (is repetitively running in the next section)
#####################################
###

library(decoupleR)
library(tidyr)
library(dplyr)

### Load
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
net=get_collectri(organism="human", split_complexes=FALSE)

# extract the module genes
Module_genes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>%
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()

### Run
print("Loading...")
data_df=Seurat::FetchData(obj_age, vars=Module_genes, layer="data")
mat=as.matrix(t(data_df))
rownames(mat)=colnames(data_df)
colnames(mat)=rownames(data_df)

print("Start running...")
acts=run_ulm(mat=mat, net=net, .source="source", .target="target", .mor="mor", minsize=5)

print("Done.")
saveRDS(acts, paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_hemostasis.rds"))

#####################################



### Manhattan Plot
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

acts=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_hemostasis.rds")
acts=acts %>% subset(grepl("CD8T\\.",condition))
# add Annot.detailed, inte, and rough info
acts_rearrange=acts %>%
  subset(p_value<0.05) %>%
  mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
acts_rearrange$Annot.rough=Annot.rough_
acts_rearrange$Annot.inter=Annot.inter_
# summarize the scores at detailed
acts_perdetailed=acts_rearrange %>%
  subset(p_value<0.05) %>%
  group_by(source, Annot.detailed) %>%
  summarize_at("score", mean) %>%
  arrange(desc(score)) %>%
  subset(score>0)
# get the TFs with top average scores
acts_top=acts_perdetailed %>%
  group_by(source) %>%
  summarize_at("score", max) %>%
  subset(score>0) %>%
  ungroup() %>%
  filter(score>=1*mean(score)) %>%
  arrange(desc(score)) %>%
  slice_max(score, n=5)
# plot
acts_perdetailed_celltypes=names(table(acts_perdetailed$Annot.detailed))
acts_perdetailed_arrange=acts_perdetailed
acts_perdetailed_arrange$Annot.detailed=forcats::fct_relevel(acts_perdetailed_arrange$Annot.detailed, acts_perdetailed_celltypes)
plot_=
  ggplot(acts_perdetailed_arrange, aes(x=source, y=score, fill=Annot.detailed, label=source)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8, preserve="single"), width=0.8, linewidth=0) +
  ggrepel::geom_text_repel(data=. %>%
                             group_by(source) %>%
                             filter((score==max(score) & (source %in% acts_top$source))) %>%
                             tidyr::complete(Annot.detailed, fill=list(score=0, source="")) %>%
                             filter(score==max(score)),
                           position=position_dodge(width=0.8, preserve="single"), 
                           size=3, force=20,
                           box.padding=0,
                           min.segment.length=0) +
  scale_fill_manual(values=paletteer::paletteer_c("ggthemes::Blue-Green Sequential", length(acts_perdetailed_celltypes))) +
  labs(x="Transcription factor", y="Score", title="hemostasis") +
  theme_classic() +
  theme(legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))

pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_hemostasis_related_MappingTo.Pseudobulk_detailed_Manhattan.pdf"), 
    height=2, width=4.5)
plot(plot_)
dev.off()

#####################################






###############################################################################################################
##################################### TF ANALYSIS ON WGCNA MODULES, BASED ON PSEUDOBULK_DETIALED #####################################
###############################################################################################################

### TF analysis on WGCNA modules based on the Pseudobulk_detailed which is used for all other TF analysis
#####################################
###

library(decoupleR)
library(tidyr)
library(dplyr)

### Load
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
net=get_collectri(organism="human", split_complexes=FALSE)

# extract the module genes
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe()
  )

# name with terms
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)
names(Module_genes)=ordered_terms

### Run
for (i in 1:length(Module_genes)) {
  module_name=gsub("\\/| |-","_",names(Module_genes)[i])
  
  print("Loading...")
  data_df=Seurat::FetchData(obj_age, vars=Module_genes[[i]], layer="data")
  mat=as.matrix(t(data_df))
  rownames(mat)=colnames(data_df)
  colnames(mat)=rownames(data_df)
  
  print("Start running...")
  run_ulm_test=function(mat_) {
    tryCatch({
      acts=run_ulm(mat=mat_, net=net, .source="source", .target="target", .mor="mor", minsize=5)
      return(acts)
    },
    error=function(msg) return(NULL))
  }
  acts=run_ulm_test(mat)
  print("Done.")
  
  if (!is.null(acts)) {
    saveRDS(acts, paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",module_name,".rds"))
  }
}

#####################################



### Manhattan Plot of TF results in each WGCNA module
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

### Get the names of the WGCNA-TF files
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe())
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)

# take only those with enriched terms
ordered_terms=ordered_terms[1:11]
filenames=gsub("\\/| |-","_",ordered_terms)

### Arrange and plot
for (i in 1:length(ordered_terms)) {
  file_name_=filenames[i]
  file_path_=paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",file_name_,".rds")
  if (file.exists(file_path_)) {
    acts=readRDS(file_path_)
    # add Annot.detailed, inte, and rough info
    acts_rearrange=acts %>%
      subset(p_value<0.05) %>%
      mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
    Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
    Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
    Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
    acts_rearrange$Annot.rough=Annot.rough_
    acts_rearrange$Annot.inter=Annot.inter_
    # summarize the scores at detailed
    acts_perdetailed=acts_rearrange %>%
      subset(p_value<0.05) %>%
      group_by(source, Annot.detailed) %>%
      summarize_at("score", mean) %>%
      arrange(desc(score)) %>%
      subset(score>0)
    # get the TFs with top average scores
    acts_top=acts_perdetailed %>%
      group_by(source) %>%
      summarize_at("score", max) %>%
      subset(score>0) %>%
      ungroup() %>%
      filter(score>=1*mean(score)) %>%
      arrange(desc(score)) %>%
      slice_max(score, n=5)
    # plot
    acts_perdetailed_celltypes=names(table(acts_perdetailed$Annot.detailed))
    acts_perdetailed_arrange=acts_perdetailed
    acts_perdetailed_arrange$Annot.detailed=forcats::fct_relevel(acts_perdetailed_arrange$Annot.detailed, acts_perdetailed_celltypes)
    
    plot_=
      ggplot(acts_perdetailed_arrange, aes(x=source, y=score, fill=Annot.detailed, label=source)) +
      geom_bar(stat="identity", position=position_dodge(width=0.8, preserve="single"), width=0.8, linewidth=0) +
      ggrepel::geom_text_repel(data=. %>%
                                 group_by(source) %>%
                                 filter(score==max(score) & (source %in% acts_top$source)) %>% 
                                 tidyr::complete(Annot.detailed, fill=list(score=0, source="")),
                               position=position_dodge(width=0.8, preserve="single"),
                               size=3, force=20,
                               box.padding=0,
                               min.segment.length=0) +
      scale_fill_manual(values=paletteer::paletteer_c("ggthemes::Blue-Green Sequential", length(acts_perdetailed_celltypes))) +
      labs(x="Transcription factor", y="Score", title=ordered_terms[i]) +
      theme_classic() +
      theme(legend.position="none",
            plot.title=element_text(size=12),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.text.y=element_text(size=9),
            axis.title.x=element_text(size=10),
            axis.title.y=element_text(size=10))
    
    pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_WGCNA_",file_name_,"_related_MappingTo.Pseudobulk_detailed_Manhattan.pdf"), 
        height=2, width=4)
    plot(plot_)
    dev.off()
    
    Sys.sleep(3)
  }
}

#####################################



### Plot the score ranking of TF for the WGCNA modules by score*freq
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

### Get the names of the WGCNA-TF files
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe())
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)

# # take only those with enriched terms
# ordered_terms=ordered_terms[1:11]
# filenames=gsub("\\/| |-","_",ordered_terms)

# remove the terms without successful TF inference
remove_or_not=c()
for (i in 1:length(ordered_terms)) {
  file_path=paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",gsub("\\/| |-","_",ordered_terms)[i],".rds")
  remove_or_not=c(remove_or_not, file.exists(file_path))
}
ordered_terms=ordered_terms[remove_or_not]
filenames=gsub("\\/| |-","_",ordered_terms)

ACTs=list()
for (i in 1:length(ordered_terms)) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",filenames[i],".rds"))
  # add Annot.detailed, inter, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  # summarize the scores per TF
  acts_perTF=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0) %>%
    mutate(supercluster=ordered_terms[i])
  ACTs[[i]]=acts_perTF
}
ACTs_all=data.table::rbindlist(ACTs)

### Plot the rank of source by weighted.score=sum(score*freq)
ACTs_all=data.table::rbindlist(ACTs)
# add the times of occurance of each source in WGCNA modules as another indicator of importance
freq_source=ACTs_all %>% group_by(source) %>% count()
ACTs_all=ACTs_all %>% left_join(., freq_source, by="source") %>%
  mutate(score.freq=score*n/length(filenames))
# order by sum(score*freq) and take only the top15 TFs
ACTs_ordered=ACTs_all %>% group_by(source) %>% 
  summarize_at("score.freq", sum) %>%
  slice_max(score.freq, n=15)
ACTs_ordered$source=forcats::fct_relevel(ACTs_ordered$source, ACTs_ordered$source)
# plot
plot_=
  ggplot(ACTs_ordered, aes(x=source, y=score.freq)) +
  geom_point() +
  labs(x=NULL, y="Weighted score") +
  theme_classic() +
  theme(panel.background=element_rect(fill="transparent", color=NA),
        plot.background=element_rect(fill="transparent", color=NA),
        legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))
pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_allWGCNAmodules_MappingTo.Pseudobulk_detailed_Ranking.pdf"), 
    height=2, width=3)
plot(plot_)
dev.off()

#####################################



### Plot the score ranking of TF for the immune proc.-related WGCNA modules by score (filtered with occurrence)
#####################################
###

library(decoupleR)
library(Seurat)
library(dplyr)
library(ggplot2)

terms_list=c("greenyellow"="tRNA modification",
             "salmon"="chromosome localization/mitotic checkpoint",
             "purple"="cytoplasmic translation",
             "yellow"="oxidative phosphorylation",
             "midnightblue"="regulation of hemopoiesis",
             "green"="T cell differentiation",
             "turquoise"="defense response",
             "brown"="leukocyte-mediated cytotoxicity",
             "black"="antigen processing and presentation",
             "red"="hemostasis",
             "blue"="humoral immunity")
terms_list=terms_list[c(6:11)] # choose the WGCNA-TFs that are related to immune proc.
short_names=gsub("\\/| |-","_",terms_list)

# remove the terms without successful TF inference
remove_or_not=c()
for (i in 1:length(short_names)) {
  file_path=paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",gsub("\\/| |-","_",short_names)[i],".rds")
  remove_or_not=c(remove_or_not, file.exists(file_path))
}
terms_list=terms_list[remove_or_not]
short_names=gsub("\\/| |-","_",terms_list)

ACTs=list()
for (i in 1:length(short_names)) {
  acts=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/WGCNA_BasedOn.Pseudobulk.detailed_HubGenes_TF_",gsub("\\/| |-","_",short_names)[i],".rds"))
  # add Annot.detailed, inter, and rough info
  acts_rearrange=acts %>%
    subset(p_value<0.05) %>%
    mutate(Annot.detailed=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",condition))
  Annot.rough_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[1]) %>% unlist()
  Annot.inter_=lapply(strsplit(acts_rearrange$Annot.detailed, split="\\."), function(x) x[2]) %>% unlist()
  Annot.inter_=paste0(Annot.rough_, ".", Annot.inter_)
  acts_rearrange$Annot.rough=Annot.rough_
  acts_rearrange$Annot.inter=Annot.inter_
  
  # summarize the scores per TF
  acts_perTF=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source, Annot.detailed) %>%
    summarize_at("score", mean) %>%
    arrange(desc(score)) %>%
    subset(score>0)
  
  # add occurrence for each source on cell types
  celltypecount_all=length(unique(acts_perTF$Annot.detailed)) # count the types of detailed celltype
  acts_occurrence_type=acts_perTF %>%
    ungroup() %>%
    group_by(source) %>%
    count() %>%
    rename(celltypen=n)
  
  # add occurrence for each source on cells
  cellcount_all=length(unique(acts_rearrange$condition)) # count the total cells
  acts_occurrence_cell=acts_rearrange %>%
    ungroup() %>%
    group_by(source) %>%
    count() %>%
    rename(celln=n)
  
  acts_perTF.occur=acts_rearrange %>%
    subset(p_value<0.05) %>%
    group_by(source) %>%
    summarize_at("score", mean) %>%
    subset(score>0) %>% 
    left_join(., acts_occurrence_type, by="source") %>%
    left_join(., acts_occurrence_cell, by="source") %>%
    arrange(desc(score)) %>%
    subset(celln>cellcount_all*0.5) # filter the sources with occurrence in most of the detailed celltypes
  
  ACTs[[i]]=acts_perTF.occur
}
ACTs_all=data.table::rbindlist(ACTs)

### Plot the rank of TFs for immune proc.-related WGCNA modules by average
acts_top_TF=ACTs_all %>% 
  group_by(source) %>%
  summarize_at("score", mean) %>%
  slice_max(score, n=15)
acts_top_TF$source=forcats::fct_relevel(
  acts_top_TF$source, acts_top_TF %>%
    dplyr::select(source) %>%
    tibble::deframe())

plot_=
  ggplot(acts_top_TF %>% arrange(desc(score)), aes(x=source, y=score)) +
  geom_point() +
  # geom_segment(aes(x="NFKB2", xend="NFKB2", y=3.3, yend=3.1), linewidth=0.1, arrow=arrow(length=unit(0.04, "npc")), color="coral4") +
  labs(x=NULL, y="Average score", title="Immune-related") +
  theme_classic() +
  theme(legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))
pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/TF_WGCNA_immuneproc_related_MappingTo.Pseudobulk_detailed_Ranking.pdf"), 
    height=2, width=3)
plot(plot_)
dev.off()

#####################################






###############################################################################################################
##################################### TF ANALYSIS ON AgeGene5K, BASED ON PSEUDOBULK_DONOR #####################################
###############################################################################################################

### TF analysis on AgeGene5K based on the Pseudobulk_donor because AgeGene5K is not cell-specific
#####################################
###

library(decoupleR)
library(tidyr)
library(dplyr)

### Load
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_donor.rds")
net=get_collectri(organism="human", split_complexes=FALSE)
# extract the PC1 top genes from gene-by-gene matrix
pc_cell=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
PC1_top_genes=pc_cell$x[,c(1,2)] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
  slice_max(abs(PC1), n=300, with_ties=FALSE) %>%
  dplyr::select(gene) %>% tibble::deframe()

### Run
print("Loading...")
data_df=Seurat::FetchData(obj_age, vars=PC1_top_genes, layer="data")
mat=as.matrix(t(data_df))
rownames(mat)=colnames(data_df)
colnames(mat)=rownames(data_df)

print("Start running...")
run_ulm_test=function(mat_) {
  tryCatch({
    acts=run_ulm(mat=mat_, net=net, .source="source", .target="target", .mor="mor", minsize=5)
    return(acts)
  },
  error=function(msg) return(NULL))
}
acts=run_ulm_test(mat)
print("Done.")

if (!is.null(acts)) {
  saveRDS(acts, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_topPC1.300_TF.rds")
}

### Analyze the TF results on PC1-top300 genes from AgeGene5K
acts=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_topPC1.300_TF.rds")
acts_arrange=acts %>% group_by(source) %>% summarize(score_avg=mean(score), p_val_avg=mean(p_value)) %>%
  subset(score_avg>0)
acts_arrange=acts_arrange %>% arrange(desc(score_avg)) %>% slice_max(score_avg, n=15)
acts_arrange$source=forcats::fct_relevel(acts_arrange$source, acts_arrange$source)
plot_=
  ggplot(acts_arrange, aes(x=source, y=score_avg)) +
  geom_point() +
  theme_classic() +
  labs(x=NULL, y="Average score", title=" ") +
  theme(legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))

pdf(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_topPC1.300_TF_all.pdf"), 
    height=2.5, width=3)
plot(plot_)
dev.off()

#####################################



### TF analysis on immune-process-relevant genes within AgeGene5K based on the Pseudobulk_donor
#####################################
###

### Get the immune-related genes from AgeGene5K
res=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
GOs_terms_selected=GOs_both[[names(GOs_both)[grepl("immune proc\\.",names(GOs_both))]]]
res_selected=res@result %>% subset(ID %in% GOs_terms_selected)
genes_immune=lapply(res_selected$core_enrichment, function(x) strsplit(x, split="/")[[1]]) %>%
  Reduce(union,.)

### Load the data
library(decoupleR)
library(tidyr)
library(dplyr)
obj_age=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
obj_age=subset(obj_age, Annot.rough %in% c("CD4T cells","CD8T cells"))
obj_rough=Seurat::AggregateExpression(obj_age, group.by="donor_id", return.seurat=T)
net=get_collectri(organism="human", split_complexes=FALSE)

### Run
print("Loading...")
data_df=Seurat::FetchData(obj_rough, vars=genes_immune, layer="data")
mat=as.matrix(t(data_df))
rownames(mat)=colnames(data_df)
colnames(mat)=rownames(data_df)

print("Start running...")
run_ulm_test=function(mat_) {
  tryCatch({
    acts=run_ulm(mat=mat_, net=net, .source="source", .target="target", .mor="mor", minsize=5)
    return(acts)
  },
  error=function(msg) return(NULL))
}
acts=run_ulm_test(mat)
print("Done.")

if (!is.null(acts)) {
  saveRDS(acts, "~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_immuneprocess_TF.rds")
}

### Analyze the TF results on PC1-top300 genes from AgeGene5K
acts=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_immuneprocess_TF.rds")
acts_arrange=acts %>% group_by(source) %>% summarize(score_avg=mean(score), p_val_avg=mean(p_value)) %>%
  subset(score_avg>0)
acts_arrange=acts_arrange %>% arrange(desc(score_avg)) %>% slice_max(score_avg, n=15)
acts_arrange$source=forcats::fct_relevel(acts_arrange$source, acts_arrange$source)
plot_=
  ggplot(acts_arrange, aes(x=source, y=score_avg)) +
  geom_point() +
  theme_classic() +
  labs(x=NULL, y="Average score", title="immune proc.") +
  theme(legend.position="none",
        plot.title=element_text(size=12),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10))

pdf(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_TF_immuneprocess.pdf"), 
    height=2.5, width=3)
plot(plot_)
dev.off()

#####################################






###############################################################################################################
##################################### KINASE ANALYSIS #####################################
###############################################################################################################

### Match the sig. genes across celltypes to the human kinase list
#####################################
###

library(Seurat)
library(scRNAtoolVis)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

# calculate the sig. genes
Idents(THEObj)="Annot.rough"
all.markers=FindAllMarkers(THEObj)
write.table(all.markers, file="~/Project_PBMCage/Results/PBMC_results/DEGs_Across_rough_celltypes.txt", sep="\t")

# extract the human kinases from http://kinhub.org/kinases.html
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
mygene=kinase_df$HGNC.Name

# plot only the kinases and mark the top ones
all.markers=read.delim("~/Project_PBMCage/Results/PBMC_results/DEGs_Across_rough_celltypes.txt", sep="\t")
all.markers_subset=all.markers %>%
  subset(gene %in% mygene)
plot_=
  jjVolcano(diffData=all.markers_subset,
            topGeneN=5,
            tile.col=paletteer::paletteer_d("ggsci::category20_d3")[11:20],
            size=3,
            pSize=0.5,
            fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(legend.text=element_text(size=12),
        legend.position="top",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossCelltypes.pdf", width=8, height=10)
plot(plot_)
dev.off()

#####################################



### Get the sig. genes across ages in each celltype at the rough
#####################################
###

library(Seurat)
library(scRNAtoolVis)

Cor_data=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz")
sig.genes_df=Cor_data %>% subset(p_val_adj<0.05 & analysis=="all") %>% group_by(age, celltypes) %>% slice_max(abs.avg_log2F, n=10)
write.table(sig.genes_df, "~/Project_PBMCage/tempt.txt", sep="\t")

#####################################



### Plot the volcano plots of kinases based on the correlation analysis results
#####################################
###

library(dplyr)

### Extract the correlation data
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# filter the kinases
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
kinase_genes=kinase_df$HGNC.Name
RNA_Cor_kinases=subset(RNA_Cor, gene %in% kinase_genes & analysis=="all")

### Plot
volcano_plot_list=
  ggplot(RNA_Cor_kinases, aes(x=rho, y=-log10(rho_pval), color=celltypes)) +
  ggrastr::geom_point_rast(size=0.5, alpha=0.3) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title=NULL, y=expression(-log[10]~pval), x=expression(Spearman~rho)) +
  facet_wrap(~celltypes, nrow=3, scales="free_y") +
  theme_classic() +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text=element_text(size=10),
        legend.position=c(0.85,0.15)) +
  guides(color=guide_legend(title="Cell Type", override.aes=list(size=2, alpha=1))) +
  ggrepel::geom_text_repel(aes(label=gene), size=3, show.legend=FALSE, max.overlaps=8)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_correlation_volcano.pdf", width=7.5, height=10)
plot(volcano_plot_list)
dev.off()

#####################################



### Enrich the age-correlated kinases - ReactomePA enrichment
#####################################
###

library(dplyr)
library(clusterProfiler)

### Extract the correlation data
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# filter the kinases
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
kinase_genes=kinase_df$HGNC.Name
RNA_Cor_kinases=subset(RNA_Cor, gene %in% kinase_genes & analysis=="all" & rho_pval<0.05)
# map to ENSEMBL ID
ensemble_map=bitr(geneID=RNA_Cor_kinases$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop=TRUE)
colnames(ensemble_map)=c("gene","ENTREZID")
RNA_Cor_kinases=left_join(RNA_Cor_kinases, ensemble_map, by="gene")

### Enrich
Enrich_up=RNA_Cor_kinases %>% subset(rho>0) %>% split(.$celltypes) %>% 
  lapply(., function(df) df %>% dplyr::select(ENTREZID) %>% unlist() %>% unname() %>%
           ReactomePA::enrichPathway(pvalueCutoff=0.05, readable=TRUE))
Enrich_down=RNA_Cor_kinases %>% subset(rho<0) %>% split(.$celltypes) %>% 
  lapply(., function(df) df %>% dplyr::select(ENTREZID) %>% unlist() %>% unname() %>%
           ReactomePA::enrichPathway(pvalueCutoff=0.05, readable=TRUE))
detach("package:org.Hs.eg.db", unload=TRUE)
detach("package:AnnotationDbi", unload=TRUE) # to release the select function to dplyr

### Plot
MODULE_ENRICH_Plots=list()
# get the top 10 terms for plotting
MODULE_ENRICH_up=lapply(Enrich_up, function(x) x@result[1:3,] %>% 
                          mutate(GeneRatio_div=as.numeric(strsplit(GeneRatio, split="/")[[1]][2])))
MODULE_ENRICH_up=lapply(MODULE_ENRICH_up, function(x) x %>% mutate(GeneRatio=Count/GeneRatio_div))

MODULE_ENRICH_down=lapply(Enrich_down, function(x) x@result[1:3,] %>% 
                          mutate(GeneRatio_div=as.numeric(strsplit(GeneRatio, split="/")[[1]][2])))
MODULE_ENRICH_down=lapply(MODULE_ENRICH_down, function(x) x %>% mutate(GeneRatio=Count/GeneRatio_div))
# plot with ggplot2
for (i in 1:length(MODULE_ENRICH_up)) {
  # upreg
  module_up=MODULE_ENRICH_up[[i]] %>% arrange(GeneRatio) %>% select(Description, GeneRatio, p.adjust, Count) %>% 
    mutate(orientation="upreg.")
  module_up$Description=sapply(module_up$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
  module_up$Description=forcats::fct_relevel(module_up$Description, module_up$Description)

  # downreg
  module_down=MODULE_ENRICH_down[[i]] %>% arrange(GeneRatio) %>% select(Description, GeneRatio, p.adjust, Count) %>% 
    mutate(orientation="downreg.")
  module_down$Description=sapply(module_down$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
  module_down$Description=forcats::fct_relevel(module_down$Description, module_down$Description)
  module_down=rbind(module_down, data.frame(Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA, orientation="downreg."))
  
  # merge upreg and downreg df
  module_=rbind(module_up, module_down)
  x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/5
  
  plot_=
    ggplot(module_, aes(x=GeneRatio, y=Description, color=p.adjust, size=Count)) +
      facet_wrap(~orientation, ncol=2, scales="free_x") +
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
    labs(title=names(MODULE_ENRICH_up)[i], y=NULL) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          title=element_text(size=11),
          strip.text=element_text(size=10),
          legend.key.width=unit(dev.size()[1] / 50, "inches"),
          legend.direction="vertical",
          legend.box="vertical",
          legend.position="right") +
    guides(color=guide_colorbar(title.theme=element_text(size=9), order=2),
           size=guide_legend(title.theme=element_text(size=9), order=1)) +
    scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(x_axis_expansion, x_axis_expansion)),
                        n.breaks=3)
  
  MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
}

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_ReactomePathEnrich.pdf", height=3, width=7)
for (p in MODULE_ENRICH_Plots) plot(p)
dev.off()

#####################################



### Enrich the age-correlated kinases - KEGG enrichment
#####################################
###

library(dplyr)
library(clusterProfiler)

### Extract the correlation data
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# filter the kinases
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
kinase_genes=kinase_df$HGNC.Name
RNA_Cor_kinases=subset(RNA_Cor, gene %in% kinase_genes & analysis=="all" & rho_pval<0.05)

### Enrich
Enrich_both=RNA_Cor_kinases %>% split(.$celltypes) %>% 
  lapply(., function(df) df %>% dplyr::select(gene, rho) %>% arrange(desc(rho)) %>% tibble::deframe() %>%
           gseGO(OrgDb="org.Hs.eg.db", keyType="SYMBOL", pvalueCutoff=1))

### Plot
MODULE_ENRICH_Plots=list()
# plot with ggplot2
for (i in 1:length(Enrich_both)) {
  module_=Enrich_both[[i]]@result %>% select(Description, NES, pvalue) %>% 
    mutate(orientation=ifelse(NES>0,"upreg.","downreg.")) %>%
    slice_max(abs(NES), n=5, by="orientation") %>%
    split(~orientation) %>%
    lapply(., function(df) df[1:5,]) %>%
    Reduce(rbind, .) %>%
    arrange(NES)
  module_$Description=forcats::fct_relevel(module_$Description, module_$Description)
  
  plot_=
    ggplot(module_, aes(x=NES, y=Description, fill=pvalue)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low=paletteer::paletteer_d("ggsci::category20_d3")[i],
                        high=paletteer::paletteer_d("ggsci::category20_d3")[i+10],
                        breaks=c(
                          min(module_$pvalue), 
                          mean(min(module_$pvalue), max(module_$pvalue)), 
                          max(module_$pvalue)
                        ),
                        labels=c(
                          sprintf("%0.01e", min(module_$pvalue)), 
                          sprintf("%0.01e", mean(min(module_$pvalue), max(module_$pvalue))), 
                          sprintf("%0.01e", max(module_$pvalue))
                        )
                        )+
    geom_text(data=module_ %>% subset(orientation=="upreg."), inherit.aes=F, 
              aes(x=0, y=Description, label=stringr::str_wrap(Description, width=45)), lineheight=0.75,
              nudge_x=-max(module_$NES)/40, hjust=1) +
    geom_text(data=module_ %>% subset(orientation=="downreg."), inherit.aes=F, 
              aes(x=0, y=Description, label=stringr::str_wrap(Description, width=45)), lineheight=0.75,
              nudge_x=max(module_$NES)/40, hjust=0) +
    theme_classic() +
    labs(title=names(Enrich_up)[i]) +
    xlim(c(-max(abs(module_$NES)), max(abs(module_$NES)))) +
    theme(title=element_text(size=11),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.title.x=element_text(size=10),
          axis.text.x=element_text(size=9),
          legend.title=element_text(size=10, angle=90, vjust=1),
          legend.text=element_text(size=9, angle=90, hjust=0.5),
          legend.key.width=unit(0.3,"cm"))
  
  MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
}

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_GO.GSEAEnrich.pdf", height=5, width=7)
for (p in MODULE_ENRICH_Plots) plot(p)
dev.off()

#####################################



### Plot the Ucell scoring heatmap of the age-correlated kinases in each celltype at the rough level
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(ComplexHeatmap)

### Extract the correlation data
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# filter the kinases
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
kinase_genes=kinase_df$HGNC.Name
RNA_Cor_kinases_All=subset(RNA_Cor, gene %in% kinase_genes & analysis=="all" & rho_pval<0.05) %>% mutate(orientation=ifelse(rho<0,"downreg.","upreg."))
RNA_Cor_kinases_F=subset(RNA_Cor, gene %in% kinase_genes & analysis=="females" & rho_pval<0.05) %>% mutate(orientation=ifelse(rho<0,"downreg.","upreg."))
RNA_Cor_kinases_M=subset(RNA_Cor, gene %in% kinase_genes & analysis=="males" & rho_pval<0.05) %>% mutate(orientation=ifelse(rho<0,"downreg.","upreg."))

### Get the pos. or neg. age-correlated genes in each celltype
Pos_and_Neg_Genes=lapply(RNA_Cor_kinases_All %>% split(.$celltypes), function(df) split(df$gene, df$orientation))
names(Pos_and_Neg_Genes)

### Load the RNA expression
SlimmedObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed.rds")

### Score
SCORE_DFs=list()
for (i in 1:length(Pos_and_Neg_Genes)) {
  subset_=subset(SlimmedObj, Annot.rough==names(Pos_and_Neg_Genes)[i])
  subset_Ucell=AddModuleScore_UCell(subset_, features=Pos_and_Neg_Genes[[i]])
  # extract the ucell scores
  score_df=subset_Ucell[[]]
  colnames(score_df)=gsub("_UCell$","",colnames(score_df))
  score_df=score_df %>% 
    select(-c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, pool_number, percent.rb)) %>% 
    mutate(sex=ifelse(sex=="female","F","M"),
           cell_id=rownames(.))
  SCORE_DFs[[i]]=score_df
}
SCORE_DFs_all=data.table::rbindlist(SCORE_DFs)
data.table::fwrite(SCORE_DFs_all, "~/Project_PBMCage/Results/PBMC_results/Kinase_UCellScore_Rough.csv.gz", sep="\t")

### Make the df
SCORE_DFs_all=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Kinase_UCellScore_Rough.csv.gz", sep="\t")
# get the ages that occur in all the groups (orientation*celltype)
age_shared=SCORE_DFs_all %>% group_by(Annot.rough, sex, agecut) %>% count() %>% 
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
# plot
ht_list=list()
SCORE_DFs=SCORE_DFs_all %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>% subset(agecut %in% age_shared) %>% group_by(agecut, sex) %>% summarize_at(c("upreg.","downreg."), mean) %>% 
    select(-downreg.) %>%
    tidyr::pivot_wider(names_from="agecut", values_from="upreg.") %>%
    tibble::column_to_rownames("sex")
  score_df_down=score_df %>% subset(agecut %in% age_shared) %>% group_by(agecut, sex) %>% summarize_at(c("upreg.","downreg."), mean) %>% 
    select(-upreg.) %>%
    tidyr::pivot_wider(names_from="agecut", values_from="downreg.") %>%
    tibble::column_to_rownames("sex")
  score_df_both.orientation=rbind(score_df_up, score_df_down)
  rownames(score_df_both.orientation)=c("F","M","F ","M ")
  # score_df_both.orientation=t(scale(t(score_df_both.orientation)))
  score_df_both.orientation=score_df_both.orientation-rowMeans(score_df_both.orientation)

  # plot
  ht_list[[i]]=
    Heatmap(score_df_both.orientation,
            name=" ",
            heatmap_legend_param=list(
              title="UCell Score",
              legend_height=unit(8, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)),
            left_annotation=rowAnnotation(group=anno_block(width=unit(0.2, "cm"),
                                                           gp=gpar(fill=c("darkred","darkblue"), col="transparent"))),
            row_split=c(1,1,2,2),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F,
            row_names_gp=gpar(fontsize=9),
            show_column_names=F,
            row_title=unique(score_df$Annot.rough),
            row_title_gp=gpar(fontsize=9),
            height=0.1*nrow(score_df_up), 
            width=0.1*ncol(score_df_up))
}
ht_list=Reduce(`%v%`, ht_list)
lgd_up_down=Legend(labels=c("upreg.","downreg."), 
                   legend_gp=gpar(fill=c("darkred","darkblue"), col=c("transparent")), 
                   title=" ",
                   direction="vertical", labels_gp=gpar(fontsize=9))
ha=HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), 
                                     col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))), 
                     annotation_name_gp=gpar(fontsize=11))
whole_ht=draw(ha %v% ht_list, merge_legend=TRUE, annotation_legend_list=lgd_up_down)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_ExprHeatmap_FvsM.pdf", height=6, width=4)
draw(whole_ht)
dev.off()

#####################################



### KEA3 analysis based on the sig. genes across ages in CD8 T cells
#####################################
###

# extract the sig. genes across the ages in CD8 T cells for KEA3 analysis: https://maayanlab.cloud/kea3/
CD8T.markers=read.delim("~/Project_PBMCage/Results/PBMC_results/DEGs_Across_ages.txt", sep="\t")
CD8T.markers_top20=CD8T.markers %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n=10)
sig.genes=CD8T.markers_top20$gene
sig.genes=sig.genes[!duplicated(sig.genes)]
write.table(data.frame(sig.genes), "~/Project_PBMCage/tempt.txt", sep="\t")

# ///https://maayanlab.cloud/kea3/#results///
# Subnetworks: MeanRank Score (TopRank Score is exactly the same), showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_subnetworks_CD8T.png"]
# Bar charts: MeanRank Score, showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_BarCharts_MeanRankScore_CD8T.png"]
# Bar charts: TopRank score, showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_BarCharts_TopRankScore_CD8T.png"]
# Tables: MeanRank Score ["~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_MeanRank.tsv"]
# Tables: TopRank score ["~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_IntegratedScaledRank.tsv"]

# select the top kinases shared between MeanRank Score table and TopRank score table
MeanRank=read_tsv("~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_MeanRank.tsv")
IntegratedScaledRank=read_tsv("~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_IntegratedScaledRank.tsv")
MeanRank_top=MeanRank %>%
  slice_min(`Mean rank`, n=50)
IntegratedScaledRank_top=IntegratedScaledRank %>%
  slice_min(`Integrated scaled rank`, n=50)
merge_=inner_join(x=MeanRank_top, y=IntegratedScaledRank_top, by="Protein")

# take the sig. genes that are in the human kinase list from http://kinhub.org/kinases.html
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
mygene=kinase_df$HGNC.Name
CD8T.markers=read.delim("~/Project_PBMCage/Results/PBMC_results/DEGs_Across_ages.txt", sep="\t")
CD8T.markers_subset=CD8T.markers %>%
  subset(gene %in% mygene)

# compare the sig. genes that are in the human kinase list and the top kinases in the KEA3 analysis
shared_kinases=intersect(CD8T.markers_subset$gene, merge_$Protein)

# mark the shared kinases in the Kinase_DEG_acrossAges
CD8T.markers_subset$cluster=as.factor(CD8T.markers_subset$cluster)
plot_=
  jjVolcano(diffData=CD8T.markers_subset, 
            myMarkers=shared_kinases,
            tile.col=paletteer::paletteer_c("ggthemes::Orange-Blue Light Diverging", 78),
            size=3,
            pSize=0.5,
            fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(legend.text=element_text(size=12),
        legend.position="top",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_KEA3Marking.pdf", width=60, height=10)
plot(plot_)
dev.off()


#####################################



### KEA3 analysis based on the correlation results of all celltypes at the rough level
#####################################
###

# extract the age-correlated genes in each rough celltype for KEA3 analysis: https://maayanlab.cloud/kea3/
donor_cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")

Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Results/")
CD8T.markers_top20=CD8T.markers %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n=10)
sig.genes=CD8T.markers_top20$gene
sig.genes=sig.genes[!duplicated(sig.genes)]
write.table(data.frame(sig.genes), "~/Project_PBMCage/tempt.txt", sep="\t")

# ///https://maayanlab.cloud/kea3/#results///
# Subnetworks: MeanRank Score (TopRank Score is exactly the same), showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_subnetworks_CD8T.png"]
# Bar charts: MeanRank Score, showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_BarCharts_MeanRankScore_CD8T.png"]
# Bar charts: TopRank score, showing top 15 results.["~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_KEA3_BarCharts_TopRankScore_CD8T.png"]
# Tables: MeanRank Score ["~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_MeanRank.tsv"]
# Tables: TopRank score ["~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_IntegratedScaledRank.tsv"]

# select the top kinases shared between MeanRank Score table and TopRank score table
MeanRank=read_tsv("~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_MeanRank.tsv")
IntegratedScaledRank=read_tsv("~/Project_PBMCage/Results/PBMC_results/Kinase_KEA3_IntegratedScaledRank.tsv")
MeanRank_top=MeanRank %>%
  slice_min(`Mean rank`, n=50)
IntegratedScaledRank_top=IntegratedScaledRank %>%
  slice_min(`Integrated scaled rank`, n=50)
merge_=inner_join(x=MeanRank_top, y=IntegratedScaledRank_top, by="Protein")

# take the sig. genes that are in the human kinase list from http://kinhub.org/kinases.html
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
mygene=kinase_df$HGNC.Name
CD8T.markers=read.delim("~/Project_PBMCage/Results/PBMC_results/DEGs_Across_ages.txt", sep="\t")
CD8T.markers_subset=CD8T.markers %>%
  subset(gene %in% mygene)

# compare the sig. genes that are in the human kinase list and the top kinases in the KEA3 analysis
shared_kinases=intersect(CD8T.markers_subset$gene, merge_$Protein)

# mark the shared kinases in the Kinase_DEG_acrossAges
CD8T.markers_subset$cluster=as.factor(CD8T.markers_subset$cluster)
plot_=
  jjVolcano(diffData=CD8T.markers_subset, 
            myMarkers=shared_kinases,
            tile.col=paletteer::paletteer_c("ggthemes::Orange-Blue Light Diverging", 78),
            size=3,
            pSize=0.5,
            fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(legend.text=element_text(size=12),
        legend.position="top",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_DEG_acrossAges_KEA3Marking.pdf", width=60, height=10)
plot(plot_)
dev.off()

#####################################



### Score the human kinase list in CD8 T cells
#####################################
###

library(Seurat)
library(ComplexHeatmap)
library(UCell)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T.NKT.NK_detailedcelltypes=names(table(THEObj$Annot.detailed))[grepl("^CD8T\\.|^OtherT\\.NKT|^NK\\.", names(table(THEObj$Annot.detailed)))]
CD8T.NKT.NK_subset=subset(THEObj, Annot.detailed %in% CD8T.NKT.NK_detailedcelltypes)

# score
kinase_df=read.delim("~/Project_PBMCage/kinhub.org_human_kinases.txt", sep="\t")
mygene=kinase_df$HGNC.Name
Kinase_scoring=list(Kinase=mygene)
CD8T_Ucell=AddModuleScore_UCell(CD8T.NKT.NK_subset, features=Kinase_scoring)

# extract the ucell scores
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
score_df=cbind(score_df, CD8T_Ucell[[]][, c("age","agecut","Annot.detailed","Annot.inter","Annot.rough")])
score_df$age=as.factor(score_df$age)
levels(score_df$age)=names(table(score_df$age))
score_df$agecut=as.factor(score_df$agecut)
levels(score_df$agecut)=names(table(score_df$agecut))
colnames(score_df)=c("Kinase", colnames(score_df)[2:ncol(score_df)])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/KinaseScore_CD8T.NKT.NK.rds")

# make the matrix for plotting
score_df=readRDS("~/Project_PBMCage/Tempt_RDS/KinaseScore_CD8T.NKT.NK.rds")
score_df=score_df %>%
  group_by(age, Annot.detailed) %>%
  summarise_at("Kinase", mean)
score_df_wide=tidyr::pivot_wider(score_df, names_from="age", values_from="Kinase", id_cols="Annot.detailed") %>% as.data.frame()
rownames(score_df_wide)=score_df_wide$Annot.detailed
score_df_wide=select(score_df_wide, -Annot.detailed)
row.names=rownames(score_df_wide)
row_reorder=c(row.names[grepl("^CD8T\\.",row.names)],
              row.names[grepl("^OtherT\\.NKT",row.names)],
              row.names[grepl("^NK\\.",row.names)])
score_df_wide=score_df_wide[row_reorder, c(order(as.integer(colnames(score_df_wide))))]

# scale the matrix
matrix_=score_df_wide

# plot
heatmap_=
  Heatmap(matrix_,
          na_col="gray99",
          # col=col_fun,
          name="UCell of Kinases",
          heatmap_legend_param=list(title_position="leftcenter-rot",
                                    title_gp=gpar(fontsize=8, fontface="bold")),
          column_title=NULL,
          # row_title="CD8 T cells",
          # row_title_gp=gpar(fontsize=8),
          row_title_rot=0,
          show_column_names=T,
          cluster_columns=F,
          cluster_rows=F,
          show_row_dend=F,
          row_names_gp=gpar(fontsize=8),
          # row_names_side="left",
          column_names_gp=gpar(fontsize=8),
          column_names_side="top",
          width=ncol(matrix_)*unit(3,"mm"), 
          height=nrow(matrix_)*unit(3,"mm"))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_scoring_CD8T.NKT.NK.pdf", width=12.5, height=4)
draw(heatmap_, heatmap_legend_side="left")
dev.off()

#####################################



### Predict kinases based on the genes related to naive/activation/migration/cytotoxic/inhibitory/effector/regulatory/proliferating/transciption/costimulatory based on the GO gaf file
#####################################
###

library(dplyr)

### Determine the subprocesses related to these scored genes
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
enrichment_list=list()
for (i in 1:length(T_cell_scoring)) {
  enrichment_list[[i]]=clusterProfiler::enrichGO(gene=T_cell_scoring[[i]], ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
}
names(enrichment_list)=names(T_cell_scoring)

### Setup subprocess terms
sub_cluster=list(`naive|'T cell differentiation'`=c("GO:0030217"),
                 `activation|'activation of immune response'`=c("GO:0002253"),
                 `migration|'chemokine-mediated signaling pathway_cellular response to chemokine_leukocyte migration'`=c("GO:0070098","GO:1990869","GO:0050900"),
                 `cytotoxic|'leukocyte mediated cytotoxicity'`=c("GO:0001909"),
                 `inhibitory|'negative regulation of immune system process'`=c("GO:0002683"),
                 `effector|'lymphocyte mediated immunity'`=c("GO:0002449"),
                 `regulatory|'T cell tolerance induction_lymphocyte anergy_regulation of tolerance induction_tolerance induction dependent upon immune response'`=c("GO:0002517","GO:0002249","GO:0002643","GO:0002461"),
                 `proliferating|'positive regulation of mitotic cell cycle_positive regulation of cell cycle phase transition_positive regulation of T cell proliferation'`=c("GO:0045931","GO:1901989","GO:0042102"),
                 `transcription|'alpha-beta T cell differentiation_alpha-beta T cell activation'`=c("GO:0046632","GO:0046631"),
                 `costimulatory|'T cell costimulation'`=c("GO:0031295"))

### Extract all the GO terms related to each supercluster
GOSim::setOntology(ont="BP", loadIC=FALSE, DIR=NULL)
allGOs=list()
for (i in 1:length(sub_cluster)) {
  other_offsprings=GOSim::getOffsprings()[sub_cluster[[i]]] %>% Reduce(c,.)
  tb=BiocGenerics::toTable(GO.db::GOBPCHILDREN); colnames(tb)[c(1,2)]=c("child","parent")
  other_offsprings=c(sub_cluster[[i]], other_offsprings, subset(tb, parent %in% sub_cluster[i]) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]
  
  j=1
  for (j in 1:10) {
    other_offsprings=c(other_offsprings, subset(tb, parent %in% other_offsprings) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]
    j=j+1
  }
  
  allGOs=c(allGOs, list(other_offsprings))
}
names(allGOs)=names(sub_cluster)
saveRDS(allGOs, "~/Project_PBMCage/Tempt_RDS/GOterms_in_ScoringPhenotypeRelated_Subclusters.rds")

### Collect the genes from each supercluster
# filter with the obj features to make sure that the genes exist
seuratobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
features_=rownames(seuratobj)
gaf_file=mgsa::readGAF("~/enrichment/goa_human.gaf.gz", evidence=NULL)
Genes_in_subclusters=list()
for (i in 1:length(allGOs)) {
  go_list=allGOs[[i]]
  geneIndices=Reduce(c, gaf_file@sets[go_list]) %>% .[!duplicated(.)] %>% .[!is.na(.)]
  genes_=gaf_file@itemAnnotations[geneIndices,]$symbol %>% .[. %in% features_]
  Genes_in_subclusters[[i]]=genes_
}
names(Genes_in_subclusters)=names(allGOs)
# save
saveRDS(Genes_in_subclusters,
        "~/Project_PBMCage/Tempt_RDS/GOterms_in_ScoringPhenotypeRelated_Subclusters_genes.rds")

### Map age-corr. genes in each celltype to the subcluster genes
# get the downreg./upreg. age-correlated genes in each celltype
Cor_analysis_DF_rough=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF_rough$celltype.level="Annot.rough"
Cor_analysis_DF_inter=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
Cor_analysis_DF_inter$celltype.level="Annot.inter"
Cor_analysis_DF_detailed=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
Cor_analysis_DF_detailed$celltype.level="Annot.detailed"
Cor_analysis_DF_subset=data.table::rbindlist(list(Cor_analysis_DF_rough, Cor_analysis_DF_inter, Cor_analysis_DF_detailed)) %>%
  subset(analysis=="all" & rho_pval<0.05) %>%
  arrange(desc(abs(rho)))
# map to the genes in subclusters
Genes_in_subclusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_in_ScoringPhenotypeRelated_Subclusters_genes.rds")
exist_genes_per_subcluster=lapply(1:length(Genes_in_subclusters),
                                  function(idx) Cor_analysis_DF_subset %>% subset(gene %in% Genes_in_subclusters[[idx]]) %>%
                                    mutate(subclusters=gsub("\\|.*","",names(Genes_in_subclusters)[idx]))) %>%
  data.table::rbindlist()

saveRDS(exist_genes_per_subcluster, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_ScoringPhenotypeRelatedSubclusters.rds")

### Analyze CD4T and CD8T
exist_genes_per_subcluster=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_ScoringPhenotypeRelatedSubclusters.rds")
exist_genes_per_subcluster=exist_genes_per_subcluster %>%
  subset(grepl("CD4T|CD8T",celltypes)) %>%
  split(.$subclusters) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.) %>% .[!duplicated(.)])
for (i in 1:length(exist_genes_per_subcluster)) {
  tempt=as.data.frame(exist_genes_per_subcluster[[i]])
  colnames(tempt)=paste0(i, ": ", names(exist_genes_per_subcluster)[i])
  write.table(tempt, paste0("~/Project_PBMCage/Results/tempt_",i,".txt"), sep="\t")
}

### Analyze at https://maayanlab.cloud/kea3/

### Ucell scoring of the kinases predicted at KEA3 in each subprocesses
library(UCell)
library(Seurat)
files=list.files("~/Project_PBMCage/Results/PBMC_results/Kinase/")
tsv_files=list()
for (i in 1:length(files)) {
  tsv_files[[i]]=read.delim(paste0("~/Project_PBMCage/Results/PBMC_results/Kinase/",files[i]), sep="\t")
  tsv_files[[i]]=tsv_files[[i]] %>% dplyr::select(Rank, Protein) %>% slice_min(Rank, n=15) %>% 
    dplyr::select(Protein) %>% tibble::deframe() %>%
    gsub("\\*","",.)
}
names(tsv_files)=gsub("\\..*","",gsub(".*_","",files))
# score
PseudobulkdObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
subset_Ucell=AddModuleScore_UCell(PseudobulkdObj, features=tsv_files)
SCORE_DF=subset_Ucell[[]]
saveRDS(SCORE_DF, "~/Project_PBMCage/Tempt_RDS/Subprocesses_kinases_Ucellw.Top15kinases_pseudobulkDetailed.rds")
# process the score
score_df=readRDS("~/Project_PBMCage/Tempt_RDS/Subprocesses_kinases_Ucellw.Top15kinases_pseudobulkDetailed.rds")
score_df=score_df %>%
  dplyr::select(-orig.ident) %>%
  mutate(sex=ifelse(sex=="female","F","M"))
colnames(score_df)=gsub("_UCell$","",colnames(score_df))
data.table::fwrite(score_df, "~/Project_PBMCage/Results/PBMC_results/Subprocesses_kinases_top15kinases_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")

### Plot the age correlation of genes predicted at KEA3
files=list.files("~/Project_PBMCage/Results/PBMC_results/Kinase/")
tsv_files_list=list()
for (i in 1:length(files)) {
  tsv_files_list[[i]]=read.delim(paste0("~/Project_PBMCage/Results/PBMC_results/Kinase/",files[i]), sep="\t")
  tsv_files_list[[i]]=tsv_files_list[[i]] %>% dplyr::select(Rank, Protein) %>% slice_min(Rank, n=10) %>% 
    dplyr::select(Protein) %>% tibble::deframe() %>%
    gsub("\\*","",.)
}
names(tsv_files_list)=gsub("\\..*","",gsub(".*_","",files))
names(tsv_files_list)=sapply(names(tsv_files_list), function(x) ifelse(x=="costimu","costimulatory",
                                                                       ifelse(x=="prolif","proliferating",x)))
# select and reorder
# tsv_files_list=tsv_files_list[c("naive","proliferating","transcription","inhibitory","regulatory","cytotoxic","effector","migration")]
tsv_files=tsv_files_list %>% unlist() %>% unname() %>% .[!duplicated(.)]
# extract the correlation
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
cor_df_subset=cor_df %>%
  subset(grepl("CD4T|CD8T",celltypes) & analysis=="all") %>%
  subset(!grepl("CD4T\\.prolif|CD4T\\.Tcm\\.HNRNPH1|CD4T\\.Tcm\\.ISG|CD4T\\.Treg", celltypes)) %>%
  subset(!grepl("CD8T\\.naive\\.HNRNPH1|CD8T\\.naive\\.TSHZ2|CD8T\\.Tem.\\ISG", celltypes)) %>%
  subset(celltypes %in% c("CD4T.naive.COTL1pos_SOX4neg","CD4T.naive.HNRNPH1","CD4T.Tcm","CD4T.Tem.",
                          "CD8T.naive.LINC02446","CD8T.Tcm.GZMB","CD8T.Tem.GZMB_FCGR3A","CD8T.Tem.GZMK_XCL1","CD8T.Tem.GZMKhi"))
cor_df_subset=lapply(1:length(tsv_files_list), function(idx) cor_df_subset %>% subset(gene %in% tsv_files_list[[idx]]))
names(cor_df_subset)=names(tsv_files_list)
row_reorder=unique(cor_df_subset[[1]]$celltypes)[c(1,2,3,4,5,6,8,7)]
cor_df_subset_reorder=cor_df_subset[c("naive","migration","cytotoxic","effector","inhibitory","regulatory","proliferating","transcription")]
# plot
ht_list=list()
lapply(cor_df_subset, function(df) range(df$rho)) # check
col_fun_sub=circlize::colorRamp2(c(-0.3, 0, 0.15), c("dodgerblue3", "white", "brown3"))
for (i in 1:length(cor_df_subset_reorder)) {
  all_cell_cor=cor_df_subset_reorder[[i]] %>% dplyr::select(gene, celltypes, rho, rho_pval) %>%
    dplyr::select(-rho_pval) %>%
    tidyr::pivot_wider(names_from="gene", values_from="rho") %>%
    tibble::column_to_rownames("celltypes") %>%
    .[row_reorder,]
  all_cell_cor[is.na(all_cell_cor)]=0
  # all_cell_cor=all_cell_cor-rowMeans(all_cell_cor)
  
  all_cell_fdr=cor_df_subset_reorder[[i]] %>% dplyr::select(gene, celltypes, rho, rho_pval) %>%
    dplyr::select(-rho) %>%
    # mutate(rho_pval=ifelse(rho_pval<0.0001,"****",
    #                        ifelse(rho_pval>=0.0001 & rho_pval<0.001, "***",
    #                               ifelse(rho_pval>=0.001 & rho_pval<0.01, "**",
    #                                      ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", ""))))) %>%
    mutate(rho_pval=ifelse(rho_pval<0.05, "*", "")) %>%
    tidyr::pivot_wider(names_from="gene", values_from="rho_pval") %>%
    tibble::column_to_rownames("celltypes") %>%
    apply(., c(1,2), function(x) ifelse(is.na(x),"",x)) %>%
    as.matrix() %>%
    .[row_reorder,]
  
  # plot
  ht_list[[i]]=
    Heatmap(all_cell_cor,
            name="mat", show_heatmap_legend=T,
            heatmap_legend_param=list(
              title=expression(Spearman~rho),
              legend_height=unit(4.6, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
            ),
            rect_gp=gpar(col="white", lwd=0.5),
            col=col_fun_sub,
            cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
            column_names_rot=90,
            show_column_names=T,
            column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
            column_title=names(cor_df_subset_reorder)[i], 
            column_title_gp=gpar(fontsize=11),
            layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
              v=pindex(all_cell_fdr, i, j)
              grid.text(v, x, y, rot=0, just="center",
                        gp=gpar(fontsize=10, fontface="bold"))
            },
            width=ncol(all_cell_cor)*unit(4.5, "mm"),
            height=nrow(all_cell_cor)*unit(5.5, "mm"))
}

ht_list=Reduce(`+`, ht_list)
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Kinase_in_subprocesses_ExprHeatmap.pdf", height=2.8, width=17.5)
draw(whole_ht)
dev.off()

#####################################






###############################################################################################################
##################################### KINASE ANALYSIS ON AgeGene5K and WGCNA modules #####################################
###############################################################################################################

### Predict kinases with AgeGene5K
#####################################
###

### Function for KEA3 analysis
KEA3_analysis=function(genes) {
  library(httr)
  library(jsonlite)
  payload=list(query_name="myQuery", gene_set=genes)
  response=POST(url="https://maayanlab.cloud/kea3/api/enrich/", body=payload, encode="json")
  json=content(response, "text")
  results=fromJSON(json)[[1]] %>% dplyr::select(Rank, TF, Score, Overlapping_Genes)
  # merge JAKs and MAPKs respectively
  results_merged=results %>% 
    mutate(TF_merged=ifelse(grepl("^MAPK[0-9]",TF),"MAPK",ifelse(grepl("^JAK[0-9]",TF),"JAK",TF))) %>%
    mutate(Score=1/as.numeric(Score)) %>%
    subset(Overlapping_Genes!="")
  Overlapping_Genes_merged=split(results_merged$Overlapping_Genes,results_merged$TF_merged) %>% 
    lapply(., function(vec) vec[!duplicated(vec)] %>% paste0(., collapse=",")) %>% unlist()
  Overlapping_Genes_merged_df=data.frame(`TF_merged`=names(Overlapping_Genes_merged), `Overlapping_Genes_merged`=Overlapping_Genes_merged) %>%
    subset(!duplicated(.))
  Overlapping_Genes_N=split(Overlapping_Genes_merged_df$Overlapping_Genes_merged, Overlapping_Genes_merged_df$TF_merged) %>% 
    lapply(., function(vec) vec %>% strsplit(., split=",") %>% .[[1]] %>% .[!duplicated(vec)] %>% length(.)) %>% unlist()
  Overlapping_Genes_N_df=data.frame(`TF_merged`=names(Overlapping_Genes_N), `Overlapping_Genes_N`=Overlapping_Genes_N) %>%
    subset(!duplicated(.))
  results_merged_arranged=results_merged %>% group_by(TF_merged) %>%
    summarize_at("Score", mean, na.rm=T) %>%
    left_join(., Overlapping_Genes_merged_df) %>%
    left_join(., Overlapping_Genes_N_df) %>%
    subset(!duplicated(TF_merged))
  
  return(results_merged_arranged)
}

### Run the KEA3 analysis on the PC1 top300 genes from AgeGene5K
# load the PC1 top genes
PCA_res=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA.rds")
topgenes=PCA_res$x[,c(1,2)] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>%
  slice_max(abs(PC1), n=300) %>%
  dplyr::select(gene) %>% tibble::deframe()

RES_allmodules=KEA3_analysis(genes=topgenes) %>%
  slice_max(Overlapping_Genes_N, n=10, with_ties=F) %>%
  arrange(desc(Score)) %>%
  mutate(Rank=as.integer(row_number())) %>%
  arrange(Rank) %>%
  rbind(., data.frame(TF_merged=paste0(rep(" ",20),collapse=""), Score=-1, Overlapping_Genes_merged="", Overlapping_Genes_N=100, Rank=10))
RES_allmodules$TF_merged=forcats::fct_relevel(RES_allmodules$TF_merged, RES_allmodules$TF_merged)
plot_RES_allmodules=
  ggplot(RES_allmodules, aes(x=TF_merged, y=Score, color=Overlapping_Genes_N)) +
  geom_point(size=2) +
  ylim(c(0,max(RES_allmodules$Score))) +
  # geom_segment(inherit.aes=F, aes(x=TF_merged, y=Score, xend=TF_merged, yend=0), 
  #              linewidth=0.1, color="grey50", linetype="dashed") +
  scale_color_gradient(low="lightpink", high="brown3", n.breaks=4) +
  labs(title=" ", x=NULL, y="Score") +
  theme_classic() +
  theme(title=element_text(size=11, color="black"),
        axis.title.x=element_text(size=10),
        axis.ticks.x=element_line(colour=c(rep("black",(nrow(RES_allmodules)-1)),"transparent")),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=9),
        legend.position="right",
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.5,"cm")) +
  guides(color=guide_colorbar(title="overlapN",
                              label.theme=element_text(size=9), 
                              title.theme=element_text(size=10)))

pdf(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_topPC1.300_kinase_all.pdf"), 
    height=2.8, width=4)
plot(plot_RES_allmodules)
dev.off()

### Run the KEA3 analysis on the immune process related genes from AgeGene5K
# load the relevant genes
res=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO.rds")
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
GOs_terms_selected=GOs_both[[names(GOs_both)[grepl("immune proc\\.",names(GOs_both))]]]
res_selected=res@result %>% subset(ID %in% GOs_terms_selected)
genes_immune=lapply(res_selected$core_enrichment, function(x) strsplit(x, split="/")[[1]]) %>%
  Reduce(union,.)

RES_immune=KEA3_analysis(genes=genes_immune) %>%
  slice_max(Overlapping_Genes_N, n=10, with_ties=F) %>%
  arrange(desc(Score)) %>%
  mutate(Rank=as.integer(row_number())) %>%
  arrange(Rank) %>%
  rbind(., data.frame(TF_merged=paste0(rep(" ",20),collapse=""), Score=-1, Overlapping_Genes_merged="", Overlapping_Genes_N=100, Rank=10))
RES_immune$TF_merged=forcats::fct_relevel(RES_immune$TF_merged, RES_immune$TF_merged)
plot_RES_immune=
  ggplot(RES_immune, aes(x=TF_merged, y=Score, color=Overlapping_Genes_N)) +
  geom_point(size=2) +
  ylim(c(0,max(RES_immune$Score))) +
  # geom_segment(inherit.aes=F, aes(x=TF_merged, y=Score, xend=TF_merged, yend=0), 
  #              linewidth=0.1, color="grey50", linetype="dashed") +
  scale_color_gradient(low="lightpink", high="brown3", n.breaks=4) +
  labs(title="immune proc.", x=NULL, y="Score") +
  theme_classic() +
  theme(title=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.ticks.x=element_line(colour=c(rep("black",(nrow(RES_allmodules)-1)),"transparent")),
        axis.title.y=element_text(size=10),
        axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=9),
        legend.position="right",
        legend.key.height=unit(0.8,"cm"),
        legend.key.width=unit(0.5,"cm")) +
  guides(color=guide_colorbar(title="overlapN",
                              label.theme=element_text(size=9), 
                              title.theme=element_text(size=10)))

pdf(paste0("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/AttentionWeight_genebygene_clustering_meanDonor_PCA_gseGO_kinase_immuneprocess.pdf"), 
    height=2.8, width=4)
plot(plot_RES_immune)
dev.off()

#####################################



### Predict kinases in WGCNA modules
# * Calling KEA3 API may sometimes be extremely slow, in this case, try do this locally.
#####################################
###

library(dplyr)

### Load the genes from WGCNA modules
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe()
  )
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="regulation of T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)
names(Module_genes)=ordered_terms
# take only those with GO terms
Module_genes=Module_genes[terms_list]

### Function for KEA3 analysis
KEA3_analysis=function(genes) {
  library(httr)
  library(jsonlite)
  payload=list(query_name="myQuery", gene_set=genes)
  response=POST(url="https://maayanlab.cloud/kea3/api/enrich/", body=payload, encode="json")
  json=content(response, "text")
  results=fromJSON(json)[[1]] %>% dplyr::select(Rank, TF, Score, Overlapping_Genes)
  # merge JAKs and MAPKs respectively
  results_merged=results %>% 
    mutate(TF_merged=ifelse(grepl("^MAPK[0-9]",TF),"MAPK",ifelse(grepl("^JAK[0-9]",TF),"JAK",TF))) %>%
    mutate(Score=1/as.numeric(Score)) %>%
    subset(Overlapping_Genes!="")
  Overlapping_Genes_merged=split(results_merged$Overlapping_Genes,results_merged$TF_merged) %>% 
    lapply(., function(vec) vec[!duplicated(vec)] %>% paste0(., collapse=",")) %>% unlist()
  Overlapping_Genes_merged_df=data.frame(`TF_merged`=names(Overlapping_Genes_merged), `Overlapping_Genes_merged`=Overlapping_Genes_merged) %>%
    subset(!duplicated(.))
  Overlapping_Genes_N=split(Overlapping_Genes_merged_df$Overlapping_Genes_merged, Overlapping_Genes_merged_df$TF_merged) %>% 
    lapply(., function(vec) vec %>% strsplit(., split=",") %>% .[[1]] %>% .[!duplicated(vec)] %>% length(.)) %>% unlist()
  Overlapping_Genes_N_df=data.frame(`TF_merged`=names(Overlapping_Genes_N), `Overlapping_Genes_N`=Overlapping_Genes_N) %>%
    subset(!duplicated(.))
  results_merged_arranged=results_merged %>% group_by(TF_merged) %>%
    summarize_at("Score", mean, na.rm=T) %>%
    left_join(., Overlapping_Genes_merged_df) %>%
    left_join(., Overlapping_Genes_N_df) %>%
    subset(!duplicated(TF_merged))
  
  return(results_merged_arranged)
}

### Run the KEA3 analysis on each module and plot the results
RES_allmodules=list()
for (i in 1:length(Module_genes)) {
  RES_allmodules[[i]]=KEA3_analysis(genes=Module_genes[[i]]) %>%
    slice_max(Overlapping_Genes_N, n=10, with_ties=F) %>%
    arrange(desc(Score)) %>%
    mutate(Rank=as.integer(row_number())) %>%
    arrange(Rank) %>%
    mutate(module=names(Module_genes)[i])
  print(paste0(i," out of ",length(Module_genes), " is done."))
}
names(RES_allmodules)=names(Module_genes)
saveRDS(RES_allmodules, "~/Project_PBMCage/Tempt_RDS/Kinase_KEA3.analysis.results.rds")

### Plot per module
RES_allmodules=readRDS("~/Project_PBMCage/Tempt_RDS/Kinase_KEA3.analysis.results.rds")
plot_=list()
for (i in 1:length(RES_allmodules)) {
  tempt_df=RES_allmodules[[i]]
  tempt_df$TF_merged=forcats::fct_relevel(tempt_df$TF_merged, tempt_df$TF_merged)
  plot_[[i]]=
    ggplot(tempt_df, aes(x=TF_merged, y=Score, color=Overlapping_Genes_N)) +
    geom_point(size=2) +
    # geom_segment(inherit.aes=F, aes(x=TF_merged, y=Score, xend=TF_merged, yend=0), 
    #              linewidth=0.1, color="grey50", linetype="dashed") +
    scale_color_gradient(low="lightpink", high="brown3") +
    labs(title=names(RES_allmodules)[i], x=NULL, y="Score") +
    theme_classic() +
    theme(title=element_text(size=11),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
          axis.text.y=element_text(size=9),
          legend.position="top",
          legend.key.height=unit(0.3,"cm")) +
    guides(color=guide_colorbar(title="overlapN",
                                label.theme=element_text(size=9), 
                                title.theme=element_text(size=10)))
}
plot_all=cowplot::plot_grid(plotlist=plot_)

### Combine the KEA3 analysis on metabolic and immunological modules, respectively
RES_allmodules=readRDS("~/Project_PBMCage/Tempt_RDS/Kinase_KEA3.analysis.results.rds")

Merged_modules=
  list(
    `Metabolism-related`=c("cytoplasmic translation","tRNA modification","chromosome localization/mitotic checkpoint",
                 "oxidative phosphorylation","regulation of hemopoiesis"),
    # metabolism=c("cytoplasmic translation","tRNA modification","chromosome localization/mitotic checkpoint"),
    `Immune-related`=c("defense response","humoral immunity","regulation of T cell differentiation",
               "leukocyte-mediated cytotoxicity","antigen processing and presentation",
               "hemostasis")
       )
RES_modulecluster=list()
for (i in 1:length(Merged_modules)) {
  module_cluster=Merged_modules[[i]]
  tempt_=data.table::rbindlist(RES_allmodules[module_cluster]) %>%
    group_by(TF_merged) %>% mutate(count=n()) %>%
    summarize(Weighted.Score=mean(Score)*count,
              Weighted.OverlapN=round(mean(Overlapping_Genes_N)*count)) %>%
    subset(!duplicated(.)) %>%
    arrange(desc(Weighted.Score))
  tempt_$TF_merged=forcats::fct_relevel(tempt_$TF_merged, tempt_$TF_merged)
  RES_modulecluster[[i]]=tempt_[1:15,]
}
names(RES_modulecluster)=names(Merged_modules)
plot_=list()
for (i in 1:length(RES_modulecluster)) {
  tempt_df=RES_modulecluster[[i]]
  plot_[[i]]=
    ggplot(tempt_df, aes(x=TF_merged, y=Weighted.Score, color=Weighted.OverlapN)) +
    geom_point(size=2) +
    # geom_segment(inherit.aes=F, aes(x=TF_merged, y=Score, xend=TF_merged, yend=0), 
    #              linewidth=0.1, color="grey50", linetype="dashed") +
    scale_color_gradient(low="lightpink", high="brown3") +
    labs(title=names(RES_modulecluster)[i], x=NULL, y="Score") +
    theme_classic() +
    theme(title=element_text(size=10),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
          axis.text.y=element_text(size=9),
          legend.position="right",
          legend.key.height=unit(0.8,"cm"),
          legend.key.width=unit(0.5,"cm")) +
    guides(color=guide_colorbar(title="overlapN",
                                label.theme=element_text(size=9), 
                                title.theme=element_text(size=10)))
}
plot_all=cowplot::plot_grid(plotlist=plot_)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Kinase_WGCNAmodules_relatedTo.MetabolismAndImmuneproc_Ranking.pdf", 
    height=2.5, width=7.5)
plot(plot_all)
dev.off()

#####################################






###############################################################################################################
##################################### BioGRID PPI ANALYSIS ON key kinases + TFs #####################################
###############################################################################################################

### PPI analysis with BioGRID
#####################################
###

### Download the BioGRID database
download.file("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.238/BIOGRID-ALL-4.4.238.tab3.zip", 
              "~/Project_PBMCage/BIOGRID_ALL_LATEST.tab3.zip", method="curl")
BioGRID_data=readr::read_delim(unzip("~/Project_PBMCage/BIOGRID_ALL_LATEST.tab3.zip"), delim="\t")

### Select
genes_of_interest=c("NFKB1"=4790,
                    "NFKB2"=4791,
                    "REL"=5966,
                    "RELA"=5970,
                    "RELB"=5971,
                    "MYC"=4609,
                    "MAPK1"=5594,
                    "LCK"=3932,
                    "AKT1"=207,
                    "JAK1"=3716,
                    "JAK3"=3718,
                    "MAPK8"=5599,
                    "MAP3K7"=6885)
colnames(BioGRID_data)
data_selected=BioGRID_data %>% subset(`Entrez Gene Interactor A` %in% genes_of_interest &
                                        `Entrez Gene Interactor B` %in% genes_of_interest) %>%
  subset(`Organism ID Interactor A`==9606 & `Organism ID Interactor B`==9606)

### Count the evidences for each pair of interactions
df_highthroughput=data_selected %>% 
  dplyr::select(-c(`Entrez Gene Interactor A`, `Entrez Gene Interactor B`)) %>%
  subset(Throughput=="High Throughput") %>%
  group_by(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  mutate(Score=ifelse(Score=="-", NA, as.numeric(Score))) %>%
  summarize(count=n(), Mean.Score=mean(Score, na.rm=T))
df_lowthroughput=data_selected %>% 
  dplyr::select(-c(`Entrez Gene Interactor A`, `Entrez Gene Interactor B`)) %>%
  subset(Throughput=="Low Throughput") %>%
  group_by(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  mutate(Score=ifelse(Score=="-", NA, as.numeric(Score))) %>%
  summarize(count=n(), Mean.Score=mean(Score, na.rm=T))
df_both=data_selected %>% 
  dplyr::select(-c(`Entrez Gene Interactor A`, `Entrez Gene Interactor B`)) %>%
  group_by(`Official Symbol Interactor A`, `Official Symbol Interactor B`) %>%
  mutate(Score=ifelse(Score=="-", NA, as.numeric(Score))) %>%
  summarize(count=n(), Mean.Score=mean(Score, na.rm=T))

### Rearrange
df_both_arranged=df_both %>% dplyr::select(-Mean.Score) %>%
  ungroup() %>% group_by(`Official Symbol Interactor A`,`Official Symbol Interactor B`) %>%
  mutate(group=c(`Official Symbol Interactor A`,`Official Symbol Interactor B`) %>% .[order(.)] %>% paste0(., collapse="+")) %>%
  ungroup() %>% group_by(group) %>% summarize(count=sum(count)) %>%
  mutate(A=gsub("\\+.*","",group), B=gsub(".*\\+","",group)) %>%
  arrange(desc(count))
df_both_arranged$group=forcats::fct_relevel(df_both_arranged$group, df_both_arranged$group)
plot=
  ggplot(df_both_arranged, aes(x=group, y=count)) +
  geom_point(size=1) +
  labs(y="# of evidence") +
  theme_light() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_blank(),
        axis.title.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        title=element_text(size=10)) +
  geom_text(aes(label=count), nudge_y=5, size=3)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PPI_Kinase.TFs.interactions_evidenceNo.pdf", height=2.5, width=8)
plot
dev.off()

#####################################






###############################################################################################################
##################################### KEY TFs AND KINASES IN HSCs #####################################
###############################################################################################################

### AKT1, JAKs, MAPKs, MYC, LCK, and NFKB in HSCs
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)


### Load the data
pseudobulk_inter=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
# from Other cells take only progenitor
pseudobulk_inter=subset(pseudobulk_inter, Annot.rough!="Other cells" | Annot.inter=="Other.progenitor")
# merge
pseudobulk_reaggregate=AggregateExpression(pseudobulk_inter, return.seurat=T, group.by=c("donor_id","Annot.rough"))
pseudobulk_reaggregate[[]]=pseudobulk_reaggregate[[]] %>%
  mutate(donor_id=gsub("_.*","",rownames(.)),
         Annot.rough=gsub(".*_","",rownames(.))) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="Other cells","Progenitor",Annot.rough)) %>%
  left_join(., pseudobulk_inter[[]] %>% dplyr::select(sex, age, donor_id) %>% subset(!duplicated(.)), by="donor_id")

### Get expr data
gene_OI_MYC=c("AKT1","JAK1","JAK3","MAPK1","MAPK8","MAP3K7","MYC")
pseudo_data=FetchData(pseudobulk_reaggregate, vars=c(gene_OI_MYC, "donor_id","age","sex","Annot.rough"), layer="data")
pseudo_data_mean=pseudo_data %>%
  dplyr::select(-c(donor_id, sex)) %>%
  group_by(age, Annot.rough) %>%
  summarize_at(colnames(.)[!grepl("age|Annot",colnames(.))], mean, na.rm=T) %>%
  ungroup() %>%
  mutate(mean.expr=rowMeans(pick(where(is.numeric), -age)))
ggplot(pseudo_data_mean, aes(x=age, y=mean.expr, color=Annot.rough)) +
  geom_point() +
  geom_smooth(aes(group=Annot.rough)) +
  labs(title="AKT/MYC")


gene_OI_NFKB=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","LCK","NFKB1","NFKB2","REL","RELA","RELB")
pseudo_data=FetchData(pseudobulk_reaggregate, vars=c(gene_OI_NFKB, "donor_id","age","sex","Annot.rough"), layer="data")
pseudo_data_mean=pseudo_data %>%
  dplyr::select(-c(donor_id, sex)) %>%
  group_by(age, Annot.rough) %>%
  summarize_at(colnames(.)[!grepl("age|Annot",colnames(.))], mean, na.rm=T) %>%
  ungroup() %>%
  mutate(mean.expr=rowMeans(pick(where(is.numeric), -age)))
ggplot(pseudo_data_mean, aes(x=age, y=mean.expr, color=Annot.rough)) +
  geom_point() +
  geom_smooth(aes(group=Annot.rough)) +
  labs(title="LCK/NFKB")



pseudo_data_arranged=pseudo_data %>%
  dplyr::select(-c(donor_id, sex)) %>%
  group_by(age, Annot.inter) %>%
  summarize_at(colnames(.)[!grepl("age|Annot",colnames(.))], mean) %>%
  ungroup() %>%
  split(.$Annot.inter) %>%
  lapply(., function(df) df %>% dplyr::select(-c(age,Annot.inter))) %>%
  lapply(., function(df) apply(df, 2, sd, na.rm=T))
pseudo_data_arranged_sd=lapply(1:length(pseudo_data_arranged), 
                          function(idx) data.frame(Annot.inter=names(pseudo_data_arranged)[idx], gene=names(pseudo_data_arranged[[idx]]), sd=pseudo_data_arranged[[idx]])) %>%
  data.table::rbindlist(.)

ggplot(pseudo_data_arranged_sd, aes(x=Annot.inter, y=sd, color=gene)) +
  geom_point() +
  geom_line(aes(group=gene)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

###



