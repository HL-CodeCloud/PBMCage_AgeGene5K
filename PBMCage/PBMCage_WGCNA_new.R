gc()
library(dplyr)
library(ggplot2)
library(WGCNA)
library(hdWGCNA)
library(Seurat)


### ANALYZE all celltypes at Annot.rough level in one wgcna project
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj, group.by=c("age","sex","Annot.rough"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
group_names=lapply(colnames(pseudobulk_obj_subset_byage), function(x) strsplit(x, split="_")[[1]])
pseudobulk_obj_subset_byage$age=lapply(group_names, function(x) x[1]) %>% unlist()
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=lapply(group_names, function(x) x[2]) %>% unlist()
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
pseudobulk_obj_subset_byage$Annot.rough=lapply(group_names, function(x) x[3]) %>% unlist()
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
# cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05) %>% dplyr::select(gene) %>% tibble::deframe() %>% .[!duplicated(.)]
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="Annot.rough",
  gene_select="fraction",
  fraction=0.2, 
  # features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name=unique(WGCNAobj$Annot.rough),
  group.by="Annot.rough",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  setwd("~/Project_PBMCage/Tempt_RDS/wgcna_results/")
  enableWGCNAThreads(nThreads=8)
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough", 
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD4T.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD4T.naive") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tnaive.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tnaive.rds")
} else {print("No selected_power found.")}

### Plot all gene expr in this subpopulation
library(ComplexHeatmap)
WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tnaive.rds")
expr_df=WGCNAobj[["RNA"]]$data
expr_df_mean=abs(rowMeans(expr_df)) %>% .[.>=quantile(abs(.), 0.75)]
expr_df_selected=expr_df[names(expr_df_mean), ]
expr_df_selected=t(scale(t(expr_df_selected)))

age_shared=colnames(expr_df_selected)
ht_list=
  Heatmap(expr_df_selected,
          name=" ",
          show_heatmap_legend=F, 
          top_annotation=
            HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"),
                                              col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                              annotation_name_gp=gpar(fontsize=10, fontfamily="mono"), show_annotation_name=T),
          col=circlize::colorRamp2(c(-2,0,2), c("#2E5A87FF","white","#A90C38FF")),
          cluster_columns=F, cluster_rows=T, show_column_names=F, show_row_names=F, show_row_dend=F,
          column_title="CD4T.naive", column_title_gp=gpar(fontsize=10), column_title_rot=0, column_title_side="bottom",
          height=unit(0.05,"mm")*nrow(expr_df_selected),
          width=unit(0.25,"mm")*ncol(expr_df_selected))
draw(ht_list)

lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=c(1:length(age_shared))[seq(1,78,by=10)],
               labels=colnames(expr_df_selected)[seq(1,78,by=10)],
               legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
lgd_score=Legend(
  title="Z-score",
  col_fun=circlize::colorRamp2(c(-2,0,2), c("#2E5A87FF","white","#A90C38FF")),
  legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
  direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
  title_gp=gpar(fontface="plain", fontsize=10)
)
legend_total=packLegend(lgd_score, lgd_age, direction="horizontal", column_gap=unit(0.3, "cm"))
draw(legend_total)

library(ClusterGVis)
getClusters(expr_df_selected)
cm <- clusterData(exp = expr_df_selected,
                  cluster.method = "mfuzz",
                  cluster.num = 4)
visCluster(object = cm,
           plot.type = "line")

#####################################



### ANALYZE CD4T.Tem
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Tem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD4T.Tem") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tem.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD4T.Tcm
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Tcm")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD4T.Tcm") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tcm.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tcm.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD4T.Treg
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Treg")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD4T.Treg") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Treg.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Treg.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD8T.naive") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.naive.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.naive.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.Tcm
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.Tcm")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD8T.Tcm") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tcm.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tcm.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.Tem
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.Tem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="CD8T.Tem") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tem.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.CTL
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.CTL")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.1 & celltypes=="CD8T.CTL") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only 1 module apart from the grey one
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.CTL.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.CTL.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="B.naive") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.naive.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.naive.rds")
} else {print("No selected_power found.")}
# only grey module found.

#####################################



### ANALYZE B.mem
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.mem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="B.mem") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.mem.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.mem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.inter
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.inter")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="B.inter") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.inter.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.inter.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.trans
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.trans")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.1 & celltypes=="B.trans") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only the grey module
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.trans.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.trans.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.pbpc
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.pbpc")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.1 & celltypes=="B.pbpc") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only the grey module
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.pbpc.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_B.pbpc.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.CD56dim
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.CD56dim")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="NK.CD56dim") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56dim.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56dim.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.CD56hi
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.CD56hi")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="NK.CD56hi") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only the grey module
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56hi.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56hi.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.prolif
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.prolif")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.2 & celltypes=="NK.prolif") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05 or <0.1, will be only the grey module
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.prolif.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_NK.prolif.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE Mono.classical
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="Mono.classical")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="Mono.classical") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.classical.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.classical.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE Mono.nonclassical
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="Mono.nonclassical")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="Mono.nonclassical") %>% dplyr::select(gene) %>% tibble::deframe()
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.nonclassical.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.nonclassical.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE Mono.inter
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="Mono.inter")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
cor_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.2 & celltypes=="Mono.inter") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05 or <0.1, will be only the grey module
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="orig.ident",
  gene_select="fraction",
  fraction=0.05, 
  features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=SetDatExpr(
  WGCNAobj,
  group_name="Aggregate",
  group.by="orig.ident",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.inter.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
possible_power=c()
for (j in 1:length(plot_list)) {
  tryCatch({
    possible_power[j]=(plot_list[[j]][[1]] %>% subset(text_color=="white"))$Power
  }, error=function(msg) {print("No Power found.")})
}
if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
print(selected_power)
# construct the network
if (selected_power!=0) {
  WGCNAobj=ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_Mono.inter.rds")
} else {print("No selected_power found.")}

#####################################



### Extract the results from the WGCNA analysis at the rough level
#####################################
###

###
WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")

# MEs
WGCNAobj=xpectr::suppress_mw(hdWGCNA::ModuleEigengenes(WGCNAobj))
MEs=hdWGCNA::GetMEs(WGCNAobj) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("cell_id",colnames(.))], names_to="color", values_to="MEs") %>%
  mutate(wgcna_slot="Annot.rough") %>%
  subset(color!="grey")

# correlation
WGCNAobj=hdWGCNA::ModuleTraitCorrelation(
  WGCNAobj,
  traits=c("age","sex","age.sex"),
  features="MEs",
  cor_method="pearson",
  group.by="Annot.rough"
)
cor_values=
  WGCNAobj@misc[["wgcna"]]$mt_cor[["cor"]] %>%
  lapply(., function(df) df %>%
           as.data.frame() %>%
           tibble::rownames_to_column("traits") %>%
           tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="cor"))
cor_values=lapply(1:length(cor_values), function(idx) cor_values[[idx]] %>% mutate(celltypes=names(cor_values)[idx])) %>%
  data.table::rbindlist()
  
pval_values=
  WGCNAobj@misc[["wgcna"]]$mt_cor[["pval"]] %>%
  lapply(., function(df) df %>%
           as.data.frame() %>%
           tibble::rownames_to_column("traits") %>%
           tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="pval"))
pval_values=lapply(1:length(pval_values), function(idx) pval_values[[idx]] %>% mutate(celltypes=names(pval_values)[idx])) %>%
  data.table::rbindlist()
  
fdr_values=
  WGCNAobj@misc[["wgcna"]]$mt_cor[["fdr"]] %>%
  lapply(., function(df) df %>%
           as.data.frame() %>%
           tibble::rownames_to_column("traits") %>%
           tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="fdr"))
fdr_values=lapply(1:length(fdr_values), function(idx) fdr_values[[idx]] %>% mutate(celltypes=names(fdr_values)[idx])) %>%
  data.table::rbindlist()

cor_pval_fdr=full_join(full_join(cor_values, pval_values, by=c("traits","color","celltypes")), fdr_values, by=c("traits","color","celltypes"))
cor_pval_fdr$wgcna_slot="Annot.rough"

# intra-module kMEs
WGCNAobj=hdWGCNA::ModuleConnectivity(
  WGCNAobj,
  group.by='Annot.rough',
  corFnc="pearson",
  corOptions="use='p'",
  harmonized=FALSE,
  assay=NULL,
  slot="data",
  group_name=unique(WGCNAobj$Annot.rough))
Module_genes=hdWGCNA::GetModules(WGCNAobj) %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols=colnames(.)[grepl("kME_",colnames(.))], names_to="kME_group", values_to="kMEs") %>%
  mutate(wgcna_slot="Annot.rough") %>%
  subset(color!="grey" & kME_group!="kME_grey")

# total kME
hb=hdWGCNA::GetHubGenes(WGCNAobj, n_hubs=nrow(WGCNAobj@misc$wgcna$wgcna_modules)) %>%
  as.data.frame() %>%
  mutate(wgcna_slot="Annot.rough") %>%
  subset(module!="grey")

saveRDS(list(MEs, cor_pval_fdr, Module_genes, hb), "~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")

#####################################




### Extract the results from the WGCNA analysis at the inter level
#####################################
###

###
files_=list.files("~/Project_PBMCage/Tempt_RDS/wgcna_results/") %>% .[grepl("^WGCNAobj_[^a]",.)]

ME_list=COR_PVAL_FDR_list=Module_kME_list=Module_totkME_list=list()
for (i in 1:length(files_)) {
  WGCNAobj=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/wgcna_results/",files_[i]))

  # MEs
  WGCNAobj=xpectr::suppress_mw(hdWGCNA::ModuleEigengenes(WGCNAobj))
  print(paste0(">>>>>> ModuleEigengenes calculation of ",files_[i]," is done. <<<<<<"))
  
  MEs=hdWGCNA::GetMEs(WGCNAobj) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("cell_id",colnames(.))], names_to="color", values_to="MEs") %>%
    mutate(wgcna_slot=gsub("WGCNAobj_|\\.rds","",files_)[i]) %>%
    subset(color!="grey")
  ME_list=c(ME_list, list(MEs))
  
  print(paste0(">>>>>> MEs of ",files_[i]," is done. <<<<<<"))
  
  # correlation
  WGCNAobj=hdWGCNA::ModuleTraitCorrelation(
    WGCNAobj,
    traits=c("age","sex","age.sex"),
    features="MEs",
    cor_method="pearson",
    group.by="orig.ident"
  )
  cor_values=
    WGCNAobj@misc[["wgcna"]]$mt_cor[["cor"]][["all_cells"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("traits") %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="cor")
  pval_values=
    WGCNAobj@misc[["wgcna"]]$mt_cor[["pval"]][["all_cells"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("traits") %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="pval")
  fdr_values=
    WGCNAobj@misc[["wgcna"]]$mt_cor[["fdr"]][["all_cells"]] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("traits") %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="fdr")
  cor_pval_fdr=full_join(full_join(cor_values, pval_values, by=c("traits","color")), fdr_values, by=c("traits","color"))
  cor_pval_fdr$wgcna_slot=gsub("WGCNAobj_|\\.rds","",files_)[i]
  COR_PVAL_FDR_list=c(COR_PVAL_FDR_list, list(cor_pval_fdr))
  print(paste0(">>>>>> Correlation of ",files_[i]," is done. <<<<<<"))
  
  # intra-module kMEs
  WGCNAobj=hdWGCNA::ModuleConnectivity(
    WGCNAobj,
    group.by='orig.ident',
    corFnc="pearson",
    corOptions="use='p'",
    harmonized=FALSE,
    assay=NULL,
    slot="data",
    group_name="Aggregate")
  Module_genes=hdWGCNA::GetModules(WGCNAobj) %>%
    as.data.frame() %>%
    tidyr::pivot_longer(cols=colnames(.)[grepl("kME_",colnames(.))], names_to="kME_group", values_to="kMEs") %>%
    mutate(wgcna_slot=gsub("WGCNAobj_|\\.rds","",files_)[i]) %>%
    subset(color!="grey" & kME_group!="kME_grey")
  Module_kME_list=c(Module_kME_list, list(Module_genes))
  print(paste0(">>>>>> kMEs of ",files_[i]," is done. <<<<<<"))
  
  # total kME
  hb=hdWGCNA::GetHubGenes(WGCNAobj, n_hubs=nrow(WGCNAobj@misc$wgcna$wgcna_modules)) %>%
    as.data.frame() %>%
    mutate(wgcna_slot=gsub("WGCNAobj_|\\.rds","",files_)[i]) %>%
    subset(module!="grey")
  Module_totkME_list=c(Module_totkME_list, list(hb))
  print(paste0(">>>>>> tot-kMEs of ",files_[i]," is done. <<<<<<"))
  
  gc()
}
ME_list_total=data.table::rbindlist(ME_list)
COR_PVAL_FDR_list_total=data.table::rbindlist(COR_PVAL_FDR_list)
Module_kME_list_total=data.table::rbindlist(Module_kME_list)
Module_totkME_list_total=data.table::rbindlist(Module_totkME_list)

saveRDS(list(ME_list_total, COR_PVAL_FDR_list_total, Module_kME_list_total, Module_totkME_list_total), "~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")
rm(list=ls());gc()

#####################################



### Analyze the correlation results at the rough level
#####################################
###

### Enrich the modules
Module_genes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[3]]
gene_lists=split(Module_genes$gene_name, Module_genes$module) %>% lapply(., function(x) x[!duplicated(x)])
gene_lists=gene_lists[names(gene_lists)!="grey"]
go_objs=lapply(gene_lists, function(x) clusterProfiler::enrichGO(x, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL"))
names(go_objs)=names(gene_lists)
saveRDS(go_objs, "~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")

### Plot the correlation only for the modules with enriched biological meanings
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
meaningful_go_objs=go_objs[lapply(go_objs, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the biological meanings
terms_list=c("greenyellow"="tRNA modification",
             "salmon"="chromosome localization/mitotic checkpoint",
             "purple"="cytoplasmic translation",
             "yellow"="oxidative phosphorylation",
             "midnightblue"="regulation of hemopoiesis",
             "green"="regulation of T cell differentiation",
             "turquoise"="defense response",
             "brown"="leukocyte-mediated cytotoxicity",
             "black"="antigen processing and presentation",
             "red"="hemostasis",
             "blue"="humoral immunity")
# arrange the trait-correlation data
cor_pval_fdr=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr %>% subset(celltypes!="all_cells") %>% subset(color %in% names(terms_list))
cor_pval_fdr=cor_pval_fdr %>% left_join(data.frame(color=names(terms_list), term=terms_list), by="color")
  
cor_pval_fdr_arranged=cor_pval_fdr %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(traits, term, cor) %>% 
           tidyr::pivot_wider(names_from="term", values_from="cor") %>%
           tibble::column_to_rownames("traits") %>%
           .[, terms_list])
all_cell_fdr=cor_pval_fdr %>%
  mutate(pval=ifelse(abs(cor)>=0.4, pval, 1)) %>%
  subset(celltypes!="all_cells") %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(traits, term, pval) %>% 
           mutate(pval=ifelse(pval<0.0001,"****",
                              ifelse(pval>=0.0001 & pval<0.001, "***",
                                     ifelse(pval>=0.001 & pval<0.01, "**",
                                            ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
           # mutate(pval=ifelse(pval<0.05,"*","")) %>%
           tidyr::pivot_wider(names_from="term", values_from="pval") %>%
           tibble::column_to_rownames("traits") %>%
           .[, terms_list])

library(ComplexHeatmap)
ht_list=lapply(1:length(cor_pval_fdr_arranged), 
               function(idx) cor_pval_fdr_arranged[[idx]] %>% .[,terms_list] %>% 
                 Heatmap(.,
                         name=" ",
                         heatmap_legend_param=list(
                           title="Correlation",
                           legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
                           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
                           title_gp=gpar(fontface="plain", fontsize=10)
                         ),
                         col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
                         cluster_columns=F, cluster_rows=F, show_column_dend=F,
                         row_names_gp=gpar(fontsize=9),
                         column_names_gp=gpar(fontsize=9),
                         row_title=names(cor_pval_fdr_arranged)[idx],
                         layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
                           v=pindex(as.matrix(all_cell_fdr[[idx]]) %>% .[,terms_list], i, j)
                           grid.text(v, x, y, rot=0, just="center",
                                     gp=gpar(fontsize=10, fontface="bold"))
                         },
                         row_title_gp=gpar(fontsize=10), row_title_rot=0, row_title_side="left",
                         height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged[[idx]]),
                         width=unit(7,"mm")*ncol(cor_pval_fdr_arranged[[idx]])))

ht_list=Reduce(`%v%`, ht_list)
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_rough.pdf", width=5.2, height=6.5)
whole_ht
dev.off()

#####################################



### Get the hub genes in each module at the rough level
#####################################
###

library(dplyr)

MEs_hub=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
# previously determined terms
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
MEs_hub_subset=MEs_hub %>% subset(module %in% names(terms_list)) %>%
  group_by(module) %>% slice_max(kME, n=15) %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module")
top15genes_per_module=split(MEs_hub_subset$gene_name, MEs_hub_subset$term)

saveRDS(top15genes_per_module, "~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_hubGenes.rds")

### Compare with the results from the Immunity dataset
top15genes_per_module=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_hubGenes.rds")
top15_immunity=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_hubGenes.rds")
PBMC_shared_modules=top15genes_per_module[names(top15genes_per_module) %in% names(top15_immunity)] %>% .[intersect(names(top15genes_per_module), names(top15_immunity))]
Immunity_shared_modules=top15_immunity[names(top15_immunity) %in% names(top15genes_per_module)] %>% .[intersect(names(top15genes_per_module), names(top15_immunity))]
shared_hb_genes=lapply(1:length(PBMC_shared_modules), function(idx) intersect(PBMC_shared_modules[[idx]], Immunity_shared_modules[[idx]]))
names(shared_hb_genes)=names(PBMC_shared_modules)

### Get the hub genes of the rest of the modules that are not shared between the two datasets
# defense response / defense response
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="turquoise") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="turquoise") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1[1:30], MEs_hub2[1:30])

# leukocyte-mediated cytotoxicity / leukocyte-mediated cytotoxicity
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="brown") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="pink") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1[1:30], MEs_hub2[1:30])

# regulation of hemopoiesis / regulation of hemopoiesis
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="midnightblue") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="yellow") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1[1:30], MEs_hub2[1:30])

# antigen processing & presentation / MHC complex assembly
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="black") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="plum1") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, MEs_hub2)

# hemostasis
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()

# cytoplasmic translation+tRNA modification / cytoplasmic translation
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="purple"|module=="greenyellow") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hub2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="salmon") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1, MEs_hub2)

# oxidative phosphorylation
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["yellow"]]@result
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="yellow") %>% 
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub1[1:100],genes_selected)

# chromosome localization/mitotic checkpoint / chromosome segregation
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["salmon"]]@result
genes_selected=paste0(go_objs_subset[1:20,"geneID"], collapse="/")
genes_selected=strsplit(genes_selected, split="/")[[1]] %>% .[!duplicated(.)]
go_objs2=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset2=go_objs2[["skyblue"]]@result
genes_selected2=paste0(go_objs_subset2[1:20,"geneID"], collapse="/")
genes_selected2=strsplit(genes_selected2, split="/")[[1]] %>% .[!duplicated(.)]
MEs_hub_both=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="salmon") %>%
  rbind(., readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
          subset(module=="skyblue")) %>%
  arrange(desc(kME)) %>% dplyr::select(gene_name) %>% tibble::deframe()
MEs_hubs=intersect(MEs_hub_both[1:50], c(genes_selected, genes_selected2))

#####################################



### Plot the expr of gene modules in NK cells
#####################################
###

### Load the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")

### Load the gene modules
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
genemodules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module") %>%
  filter(!is.na(term))
genemodules=split(genemodules$gene_name, genemodules$term)

### Extract the expr of each module
EXPR_All=list()
for (i in 1:length(genemodules)) {
  expr=Seurat::FetchData(pseudobulk_data %>% subset(Annot.rough=="NK cells"), 
                         vars=c("sex","age", genemodules[[i]]), 
                         layer="data")
  expr.merged=rowMeans(expr[,3:ncol(expr)], na.rm=T)
  expr$merged=expr.merged
  expr=expr %>% dplyr::select(age, sex, merged) %>% mutate(term=names(genemodules)[i])
  EXPR_All[[i]]=expr
}
EXPR_All_Combo=data.table::rbindlist(EXPR_All) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  mutate(term=ifelse(term=="T cell differentiation","regulation of T cell differentiation",term)) %>%
  mutate(term=ifelse(grepl("chromosome localization",term), "chromosome localization", term)) %>%
  mutate(term=gsub("regulation of ","regulation of\n",term)) %>%
  mutate(term=gsub("leukocyte-mediated ","leukocyte-mediated\n",term)) %>%
  mutate(term=gsub("antigen processing ","antigen processing\n",term))

### Plot
terms_list=gsub("T cell differentiation","regulation of T cell differentiation",terms_list)
EXPR_All_Combo$term=forcats::fct_relevel(EXPR_All_Combo$term, gsub("regulation of ","regulation of\n",
                                                                   gsub("leukocyte-mediated ","leukocyte-mediated\n",
                                                                        gsub("antigen processing ","antigen processing\n",
                                                                             gsub("/mitotic checkpoint","",terms_list)))))
plot_=
  ggplot(EXPR_All_Combo, aes(x=age, y=merged, color=sex)) +
  facet_wrap(~term, scales="free", ncol=6) +
  ggrastr::geom_point_rast(size=0.05, shape=1, alpha=0.05) +
  geom_smooth(aes(group=sex), method="loess", show.legend=T, fill="grey93", linewidth=0.25) +
  ggpubr::stat_cor(aes(group=sex), method="spearman", show.legend=F, size=3.5) +
  theme_minimal() +
  labs(y="Avg. expression", subtitle="Yazar et al. (2022) dataset") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.925,0.2),
        legend.direction="horizontal") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WCCNA_genemoduleExpr_NK.pdf", height=4, width=14.5)
plot_
dev.off()

#####################################



### Plot the expr of the module "hemostasis" in CD8T and Other cells
#####################################
###

### Load the expr data
library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
celltypes_in_othercells_toremove=names(table(Seurat_obj$Annot.inter)) %>% .[!grepl("Other\\.[^P]",.)]
Seurat_obj=subset(Seurat_obj, Annot.inter %in% celltypes_in_othercells_toremove)
sex_df=data.frame(donor_id=Seurat_obj$donor_id,
                  sex=Seurat_obj$sex)
sex_df=sex_df[!duplicated(sex_df$donor_id),]
# pseudobulk at the rough level
celltypes=names(table(Seurat_obj$Annot.rough))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.rough==celltypes[i])
  pseudo_rough=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("donor_id","age"))
  pseudo_rough$donor_id=gsub("_.*","",colnames(pseudo_rough))
  pseudo_rough$age=as.numeric(gsub(".*_","",colnames(pseudo_rough)))
  pseudo_rough$Annot.rough=celltypes[i]
  PseudoList[[i]]=pseudo_rough
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough_Only.Plt.In.OtherCells.rds")
# -- add sex metadata
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough_Only.Plt.In.OtherCells.rds")
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("male","female")
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough_Only.Plt.In.OtherCells.rds")

### Load the gene modules
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough_Only.Plt.In.OtherCells.rds")
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
genemodules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module") %>%
  filter(!is.na(term))
genemodules=split(genemodules$gene_name, genemodules$term)
genemodules_=genemodules[["hemostasis"]]

### Extract the expr of each module
EXPR_All=list()
for (i in 1:2) {
  expr=Seurat::FetchData(pseudobulk_data %>% subset(Annot.rough==c("CD8T cells","Other cells")[i]), 
                         vars=c("sex","age", genemodules_), 
                         layer="data")
  expr.merged=rowMeans(expr[,3:ncol(expr)], na.rm=T)
  expr$merged=expr.merged
  expr=expr %>% dplyr::select(age, sex, merged) %>% tibble::rownames_to_column("Annot.rough")
  EXPR_All[[i]]=expr
}

# ### Plot with facets
# EXPR_All_Combo=data.table::rbindlist(EXPR_All) %>%
#   mutate(sex=ifelse(sex=="female","F","M")) %>%
#   mutate(Annot.rough=gsub("_[0-9]{1,4}-[0-9]{1,4}_.*","",Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough=="Other cells","Platelets",Annot.rough))
# EXPR_All_Combo$Annot.rough=forcats::fct_relevel(EXPR_All_Combo$Annot.rough, c("Platelets","CD8T cells"))
# plot_=
#   ggplot(EXPR_All_Combo, aes(x=age, y=merged, color=sex)) +
#   facet_wrap(~Annot.rough, scales="free", ncol=3) +
#   ggrastr::geom_point_rast(size=0.05, shape=1, alpha=0.05) +
#   geom_smooth(aes(group=sex), method="loess", show.legend=T, fill="grey93", linewidth=0.25) +
#   ggpubr::stat_cor(aes(group=sex), method="spearman", show.legend=F, size=3.5) +
#   theme_minimal() +
#   labs(y="Avg. expression") +
#   theme(plot.background=element_rect(color="transparent", fill="white"),
#         plot.subtitle=element_text(size=11),
#         panel.spacing=unit(5,"lines"),
#         axis.title=element_text(size=10),
#         axis.text=element_text(size=9),
#         legend.text=element_text(size=9),
#         legend.title=element_text(size=10),
#         strip.text=element_text(size=10),
#         legend.position="right") +
#   guides(color=guide_legend(override.aes=list(fill="transparent")))

### Plot seperately
EXPR_All_split=lapply(EXPR_All, function(df) df %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  mutate(Annot.rough=gsub("_[0-9]{1,4}-[0-9]{1,4}_.*","",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="Other cells","Platelets",Annot.rough)))
plot_1=
  ggplot(EXPR_All_split[[1]], aes(x=age, y=merged, color=sex)) +
  ggrastr::geom_point_rast(size=0.05, shape=1, alpha=0.05) +
  geom_smooth(aes(group=sex), method="loess", show.legend=T, fill="grey93", linewidth=0.25) +
  ggpubr::stat_cor(aes(group=sex), method="spearman", show.legend=F, size=3.5) +
  theme_minimal() +
  labs(y="Avg. expression", x="age", title="CD8T cells") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        plot.title=element_text(size=11),
        panel.spacing=unit(5,"lines"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="right") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))
plot_2=
  ggplot(EXPR_All_split[[2]], aes(x=age, y=merged, color=sex)) +
  ggrastr::geom_point_rast(size=0.05, shape=1, alpha=0.05) +
  geom_smooth(aes(group=sex), method="loess", show.legend=T, fill="grey93", linewidth=0.25) +
  ggpubr::stat_cor(aes(group=sex), method="spearman", show.legend=F, size=3.5) +
  theme_minimal() +
  labs(y="Avg. expression", x="age", title="Platelets") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        plot.title=element_text(size=11),
        panel.spacing=unit(5,"lines"),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="right") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WCCNA_genemoduleExpr_hemostasis.CD8T.pdf", height=2, width=4)
plot_1
dev.off()
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WCCNA_genemoduleExpr_hemostasis.Plt.pdf", height=2, width=4)
plot_2
dev.off()

# # check the expression of the genes within the hemostasis module in all the cell types
# EXPR_All=plots_all=list()
# for (i in 1:8) {
#   expr=Seurat::FetchData(pseudobulk_data %>% subset(Annot.rough==names(table(pseudobulk_data$Annot.rough))[i]), 
#                          vars=c("sex","age", genemodules_), 
#                          layer="data") %>% 
#     tibble::rownames_to_column("Annot.rough") %>%
#     mutate(Annot.rough=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",Annot.rough)) %>% 
#     tidyr::pivot_longer(cols=colnames(.)[!grepl("^sex$|^age|Annot\\.",colnames(.))], names_to="gene", values_to="expr") %>%
#     group_by(Annot.rough, gene) %>%
#     summarize_at("expr", mean) %>%
#     ungroup() %>%
#     group_by(Annot.rough) %>%
#     arrange(desc(expr), .by_group=T)
#   
#   EXPR_All[[i]]=expr
#   
#   data_df_1=EXPR_All[[i]] %>% ungroup() %>% slice_max(expr, n=20)
#   data_df_1$gene=forcats::fct_relevel(data_df_1$gene, data_df_1$gene)
#   plots_all[[i]]=
#     ggplot(data_df_1,
#            aes(x=gene, y=expr)) +
#     geom_bar(stat="identity", width=0.4, fill="grey70") +
#     theme_minimal() +
#     labs(y="Avg. expression", x=NULL, subtitle=ifelse(i==7, "Platelets", names(table(pseudobulk_data$Annot.rough))[i])) +
#     theme(plot.background=element_rect(color="transparent", fill="white"),
#           plot.subtitle=element_text(size=11),
#           axis.title=element_text(size=10),
#           axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
#           axis.text.y=element_text(size=9),
#           legend.text=element_text(size=9),
#           legend.title=element_text(size=10),
#           strip.text=element_text(size=10),
#           legend.position=c(0.85,0.1),
#           legend.direction="horizontal") +
#     guides(color=guide_legend(override.aes=list(fill="transparent")))
# }
# plots_all=cowplot::plot_grid(plotlist=plots_all, align="hv", ncol=1)

### Plot the top genes in CD8T vs. Other cells
EXPR_All=list()
for (i in 1:2) {
  expr=Seurat::FetchData(pseudobulk_data %>% subset(Annot.rough==c("CD8T cells","Other cells")[i]), 
                         vars=c("sex","age", genemodules_), 
                         layer="data") %>% 
    tibble::rownames_to_column("Annot.rough") %>%
    mutate(Annot.rough=gsub("_[0-9]{1,4}-[0-9]{1,4}.*","",Annot.rough)) %>% 
    tidyr::pivot_longer(cols=colnames(.)[!grepl("^sex$|^age|Annot\\.",colnames(.))], names_to="gene", values_to="expr") %>%
    group_by(Annot.rough, gene) %>%
    summarize_at("expr", mean) %>%
    ungroup() %>%
    group_by(Annot.rough) %>%
    arrange(desc(expr), .by_group=T)
    
  EXPR_All[[i]]=expr
}
EXPR_All[[1]] %>% ungroup() %>% slice_max(expr, n=20) %>% dplyr::select(gene) %>% tibble::deframe()
EXPR_All[[2]]

# plot the top genes within the "hemostasis" module in CD8T vs. other cells
data_df_1=EXPR_All[[1]] %>% ungroup() %>% slice_max(expr, n=20)
data_df_1$gene=forcats::fct_relevel(data_df_1$gene, data_df_1$gene)
plot1=
  ggplot(data_df_1,
         aes(x=gene, y=expr)) +
  geom_bar(stat="identity", width=0.4, fill="grey70") +
  theme_minimal() +
  labs(y="Expression", x=NULL, subtitle="CD8T cells") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.85,0.1),
        legend.direction="horizontal") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

data_df_2=EXPR_All[[2]] %>% ungroup() %>% slice_max(expr, n=20)
data_df_2$gene=forcats::fct_relevel(data_df_2$gene, data_df_2$gene)
plot2=
  ggplot(data_df_2,
         aes(x=gene, y=expr)) +
  geom_bar(stat="identity", width=0.4, fill="grey70") +
  theme_minimal() +
  labs(y="Expression", x=NULL, subtitle="Platelets") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.85,0.1),
        legend.direction="horizontal") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

plots_=cowplot::plot_grid(plotlist=list(plot2, plot1), align="hv", ncol=1)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_HubgenesExpr_Module.hemostasis_CD8T.Plt.pdf", height=6, width=5)
plots_
dev.off()

#####################################



### Plot the expr of the module "hemostasis" in CD8T and Other cells, newer
#####################################
###

# find the OtherCell-expressing genes in the module
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>% # hemostasis
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds") %>% 
  subset(cluster=="Other cells") %>%
  dplyr::select(gene) %>% tibble::deframe()
MEs_hubs_PLT=intersect(MEs_hub1, markers_rough) # it turns out to be related to blood coagulation
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["red"]]@result
genes_selected=strsplit(paste0(go_objs_subset[grepl("coagulation",go_objs_subset$Description),"geneID"], collapse="/"), split="/")[[1]]
MEs_hubs_plt=intersect(MEs_hubs_PLT, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

# find the CD8T-expressing genes in the module
MEs_hub1=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  subset(module=="red") %>% # hemostasis
  arrange(desc(kME)) %>%
  dplyr::select(gene_name) %>% tibble::deframe()
markers_rough=readRDS("~/Project_PBMCage/Results/PBMC_results/Cell_Markers_rough.rds") %>% 
  subset(cluster=="CD8T cells") %>%
  dplyr::select(gene) %>% tibble::deframe()
MEs_hubs_CD8T=intersect(MEs_hub1, markers_rough) # it turns out to be related to cell-cell adhesion
go_objs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
go_objs_subset=go_objs[["red"]]@result
genes_selected=strsplit(paste0(go_objs_subset[grepl("adhesion|junction|organization",go_objs_subset$Description),"geneID"], collapse="/"), split="/")[[1]]
MEs_hubs_cd8t=intersect(MEs_hubs_CD8T, genes_selected)[1:30] # take the top 30 hub genes from the enriched module

### Load the expr data
pseudobulk_rough=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough_Only.Plt.In.OtherCells.rds")
pseudobulk_subset=subset(pseudobulk_rough, Annot.rough %in% c("CD8T cells","Other cells"))
expr_data_plt=FetchData(pseudobulk_subset, vars=c("donor_id", "age", "Annot.rough", MEs_hubs_plt), layer="data") %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|age|Annot",colnames(.))], values_to="expr", names_to="gene") %>%
  mutate(Submodule="blood coagulation")
expr_data_cd8t=FetchData(pseudobulk_subset, vars=c("donor_id", "age", "Annot.rough", MEs_hubs_cd8t), layer="data") %>%
  as.data.frame() %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|age|Annot",colnames(.))], values_to="expr", names_to="gene") %>%
  mutate(Submodule="cell-cell junction")
expr_data_total=rbind(expr_data_plt, expr_data_cd8t)

### Plot gene expression
plot_=
ggplot(expr_data_total, aes(x=gene, y=expr, color=Submodule)) +
  ggrastr::geom_point_rast(size=0.5) +
  facet_wrap(~Annot.rough, ncol=1, scales="free_y", 
             labeller=as_labeller(c("CD8T cells"="CD8T cells","Other cells"="Platelets"))) +
  scale_color_manual(values=c("deepskyblue","grey50")) +
  theme_minimal() +
  labs(y="Expression", x=NULL) +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
        axis.text.y=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  guides(color=guide_legend(override.aes=list(size=2)))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_HubgenesExpr_Module.hemostasis_CD8T.Plt_new.pdf", height=3.5, width=14.5)
plot_
dev.off()

#####################################



### Calculate the correlation between expr of modules and ages at Annot.detailed
#####################################
###

### Load the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
celltypes=names(table(pseudobulk_data$Annot.detailed))

### Load the gene modules
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
genemodules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module") %>%
  filter(!is.na(term))
genemodules=split(genemodules$gene_name, genemodules$term)

### Extract the expr of each module
cor_All=pval_All=c()
for (i in 1:length(genemodules)) {
  for (j in 1:length(celltypes)) {
    for (s in c("male","female")) {
      expr=
        tryCatch({Seurat::FetchData(pseudobulk_data %>% subset(Annot.detailed==celltypes[j] & sex==s), 
                                  vars=c("age", genemodules[[i]]), 
                                  layer="data")},
                 error=function(msg) return(NULL))
      # to avoid error-raising small samples
      if (!is.null(expr)) {
        expr.merged=rowMeans(expr[,2:ncol(expr)], na.rm=T)
        expr$merged=expr.merged
        expr=expr %>% dplyr::select(age, merged) %>%
          tibble::rownames_to_column("Annot.detailed") %>%
          mutate(Annot.detailed=gsub("_[0-9].*","",Annot.detailed))
        # calculate the correlation
        cor_test=cor.test(expr$merged, expr$age, method="spearman")
        rho=cor_test$estimate; names(rho)=paste0(names(genemodules)[i],"___",celltypes[j],"===",s)
        pval=unname(cor_test$p.value)
        # save the correlation results
        cor_All=c(cor_All, rho)
        pval_All=c(pval_All, pval)
      }
    }
    print(paste0(celltypes[j], " is done."))
  }
  print(paste0("----- ",names(genemodules)[i], " is done. -----"))
}
# arrange df
Cor_All_df=data.frame(group=names(cor_All), cor=cor_All, pval=pval_All) %>%
  mutate(term=gsub("___.*","",group), Annot.detailed=gsub(".*___|===.*","",group), sex=gsub(".*===","",group)) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  dplyr::select(-group) %>%
  tibble::remove_rownames() %>%
  mutate(term=ifelse(term=="T cell differentiation","regulation of T cell differentiation",term))
data.table::fwrite(Cor_All_df, "~/Project_PBMCage/Results/PBMC_results/WGCNA.module.genes_meanExpr.cor.Wage_Annot.detailed.txt.gz", sep="\t")

#####################################



### Calculate the correlation between expr of superclusters and ages at Annot.detailed
#####################################
###

### Load the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
celltypes=names(table(pseudobulk_data$Annot.detailed))

### Load the gene modules and group into superclusters
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
genemodules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module") %>%
  filter(!is.na(term))
# group into superclusters
cell_activities_and_metabolism=c("tRNA modification","chromosome localization/mitotic checkpoint","cytoplasmic translation")
autophagy=c("oxidative phosphorylation")
leukocyte_development=c("regulation of hemopoiesis")
immune_system_processes=c("T cell differentiation","defense response","leukocyte-mediated cytotoxicity",
                          "antigen processing and presentation","hemostasis","humoral immunity")
genemodules=genemodules %>%
  mutate(term_grouped=ifelse(term %in% cell_activities_and_metabolism, "cell activity and metabolism",
                             ifelse(term %in% autophagy, "autophagy",
                                    ifelse(term %in% leukocyte_development, "leuk. devel.",
                                           ifelse(term %in% immune_system_processes, "immune proc.", term)))))
genemodules=split(genemodules$gene_name, genemodules$term_grouped)

### Extract the expr of each module
cor_All=pval_All=c()
for (i in 1:length(genemodules)) {
  for (j in 1:length(celltypes)) {
    for (s in c("male","female")) {
      expr=
        tryCatch({Seurat::FetchData(pseudobulk_data %>% subset(Annot.detailed==celltypes[j] & sex==s), 
                                    vars=c("age", genemodules[[i]]), 
                                    layer="data")},
                 error=function(msg) return(NULL))
      # to avoid error-raising small samples
      if (!is.null(expr)) {
        expr.merged=rowMeans(expr[,2:ncol(expr)], na.rm=T)
        expr$merged=expr.merged
        expr=expr %>% dplyr::select(age, merged) %>%
          tibble::rownames_to_column("Annot.detailed") %>%
          mutate(Annot.detailed=gsub("_[0-9].*","",Annot.detailed))
        # calculate the correlation
        cor_test=cor.test(expr$merged, expr$age, method="spearman")
        rho=cor_test$estimate; names(rho)=paste0(names(genemodules)[i],"___",celltypes[j],"===",s)
        pval=unname(cor_test$p.value)
        # save the correlation results
        cor_All=c(cor_All, rho)
        pval_All=c(pval_All, pval)
      }
    }
    print(paste0(celltypes[j], " is done."))
  }
  print(paste0("----- ",names(genemodules)[i], " is done. -----"))
}
# arrange df
Cor_All_df=data.frame(group=names(cor_All), cor=cor_All, pval=pval_All) %>%
  mutate(term=gsub("___.*","",group), Annot.detailed=gsub(".*___|===.*","",group), sex=gsub(".*===","",group)) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  dplyr::select(-group) %>%
  tibble::remove_rownames()
data.table::fwrite(Cor_All_df, "~/Project_PBMCage/Results/PBMC_results/WGCNA.moduleGrouped.genes_meanExpr.cor.Wage_Annot.detailed.txt.gz", sep="\t")

### Plot
library(ComplexHeatmap)
cor_res=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA.moduleGrouped.genes_meanExpr.cor.Wage_Annot.detailed.txt.gz")
# rename the terms according to the manuscript's final version
cor_res=cor_res %>%
  mutate(term=ifelse(term=="cell activity and metabolism","cellular metab. and homeostasis",term))
mat_=cor_res %>% 
  subset(sex=="M") %>%
  dplyr::select(-pval) %>% 
  tidyr::pivot_wider(id_cols="Annot.detailed", names_from="term", values_from="cor") %>%
  tibble::column_to_rownames("Annot.detailed")
mat_p=cor_res %>% 
  subset(sex=="M") %>%
  dplyr::select(-cor) %>% 
  tidyr::pivot_wider(id_cols="Annot.detailed", names_from="term", values_from="pval") %>%
  tibble::column_to_rownames("Annot.detailed") %>%
  apply(., c(1,2), function(x) -log10(x))
mat_annotated=rownames(mat_p[apply(mat_p, 1, max)>-log10(0.05),])
range(mat_)
col_fun=circlize::colorRamp2(c(-0.4, 0, 0.4), c("dodgerblue3", "white", "brown3"))
ha=rowAnnotation(foo=anno_mark(at=match(mat_annotated, rownames(mat_p)), 
                               labels=mat_annotated, labels_gp=gpar(fontsize=8)))
ht_=
  Heatmap(mat_, 
          name=" ",
          heatmap_legend_param=list(
            title=expression(Spearman~rho),
            legend_height=unit(2.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            grid.circle(x=x, y=y, r=ifelse(mat_p[i, j]>-log10(0.05), mat_p[i, j], 0)/2 * min(unit.c(width, height)), 
                        gp=gpar(fill=col_fun(mat_[i, j]), col=NA))
          }, 
          cluster_rows=FALSE, cluster_columns=TRUE,
          show_column_dend=FALSE,
          show_row_names=FALSE, show_column_names=TRUE,
          right_annotation=ha,
          row_names_side="left",
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9), column_names_rot=-20,
          height=unit(2,"mm")*nrow(mat_),
          width=unit(8,"mm")*ncol(mat_))
ht_opt$HEATMAP_LEGEND_PADDING=unit(-3,"mm")
ht_=draw(ht_, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WCCNA_genemoduleExpr_corWage_Annot.detailed.pdf", height=6.75, width=4)
ht_
dev.off()

#####################################



### Map the gene clusters at the rough level to the 12 superclusters by mgoSim on the enriched GO terms
#####################################
###

### Load the terms belonging to the 12 superclusters
supercluster_terms=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")
supercluster_terms_up=supercluster_terms[[1]]
supercluster_terms_down=supercluster_terms[[2]]
supercluster_goterms=lapply(1:length(supercluster_terms_up), function(idx) c(supercluster_terms_up[[idx]], supercluster_terms_down[[idx]]) %>% .[!duplicated(.)])
names(supercluster_goterms)=gsub("\\|.*","",names(supercluster_terms_up))

### Load the enriched terms from the wgcna results
wgcna_terms=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
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
wgcna_terms=wgcna_terms[names(terms_list)]
wgcna_goterms=lapply(wgcna_terms, function(x) x@result$ID)
names(wgcna_goterms)=terms_list

### go similarity analysis
# create the function to analyze the enrichment similarity
d=GOSemSim::godata('org.Hs.eg.db', ont="BP", computeIC=FALSE)
similarity=function(description_vector1, description_vector2) {
  sim_=GOSemSim::mgoSim(
    description_vector1,
    description_vector2,
    semData=d,
    measure="Wang", combine="BMA"
  )
  return (sim_)
}
all_sim_results=all_sim_results_name=c()
for (i in 1:length(supercluster_goterms)) {
  superterm=supercluster_goterms[[i]]
  for (j in 1:length(wgcna_goterms)) {
    moduleterm=wgcna_goterms[[j]]
    similarity_result=similarity(superterm, moduleterm)
    all_sim_results=c(all_sim_results, similarity_result)
    all_sim_results_name=c(all_sim_results_name, paste0(gsub("\\|.*","",names(supercluster_goterms)[i]),"_vs_",names(wgcna_goterms)[j]))
  }
}
names(all_sim_results)=all_sim_results_name
# process the df
Jaccard_SIM_pos=data.frame(pair=names(all_sim_results), jaccard=all_sim_results) %>%
  mutate(superterm=gsub("_vs_.*","",pair), moduleterm=gsub(".*_vs_","",pair))
saveRDS(Jaccard_SIM_pos, "~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms.vs.WGCNAmoduleTerms_GOtermSimilarity.rds")

# plot the similarity
Jaccard_SIM_pos=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms.vs.WGCNAmoduleTerms_GOtermSimilarity.rds")
Jaccard_SIM_df=Jaccard_SIM_pos %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="superterm", values_from="jaccard") %>%
  tibble::column_to_rownames("moduleterm")
range(Jaccard_SIM_df) # check
col_fun_up=circlize::colorRamp2(c(0,0.5), c("white", "brown3"))
ht_gene_sim=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=TRUE, 
          heatmap_legend_param=list(
            title="Similarity",
            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontsize=10, fontface="plain")
          ),
          col=col_fun_up,
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          cluster_columns=T, cluster_rows=T, 
          show_column_dend=F, show_row_dend=F, column_dend_side="top", row_dend_side="right",
          row_names_gp=gpar(fontsize=10), show_row_names=F, row_names_side="left", 
          column_title="", column_title_gp=gpar(fontsize=11),
          show_column_names=T, 
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(Jaccard_SIM_df), 
          width=unit(5, "mm")*ncol(Jaccard_SIM_df))

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_gene_sim=draw(ht_gene_sim, heatmap_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/KeyTerms.vs.WGCNAmoduleTerms_GOSimilarity.pdf", height=3.8, width=3)
ht_gene_sim
dev.off()

# plot the scaled similarity
Jaccard_SIM_pos=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_KeyTerms.vs.WGCNAmoduleTerms_GOtermSimilarity.rds")
Jaccard_SIM_df=Jaccard_SIM_pos %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="superterm", values_from="jaccard") %>%
  tibble::column_to_rownames("moduleterm")
Jaccard_SIM_df=scale(Jaccard_SIM_df)
range(Jaccard_SIM_df) # check
col_fun_up=circlize::colorRamp2(c(-2,0,2), c("dodgerblue3", "white", "brown3"))
ht_gene_sim=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=TRUE, 
          heatmap_legend_param=list(
            title="Scaled similarity",
            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontsize=10, fontface="plain")
          ),
          col=col_fun_up,
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          cluster_columns=T, cluster_rows=T, 
          show_column_dend=F, show_row_dend=F, column_dend_side="top", row_dend_side="right",
          row_names_gp=gpar(fontsize=10), show_row_names=T, row_names_side="right", 
          column_title="", column_title_gp=gpar(fontsize=11),
          show_column_names=T, 
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(Jaccard_SIM_df), 
          width=unit(5, "mm")*ncol(Jaccard_SIM_df))

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_gene_sim=draw(ht_gene_sim, heatmap_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/KeyTerms.vs.WGCNAmoduleTerms_GOSimilarity_scaled.pdf", height=3.8, width=5.5)
ht_gene_sim
dev.off()

#####################################



### Map the gene clusters at the rough level to the 12 superclusters by mclusterSim
#####################################
###

### Extract the WGCNA module genes
Module_genes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[3]]
gene_lists=split(Module_genes$gene_name, Module_genes$module) %>% lapply(., function(x) x[!duplicated(x)])
gene_lists=gene_lists[names(gene_lists)!="grey"]
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
gene_lists=gene_lists[names(terms_list)]
names(gene_lists)=terms_list
length(gene_lists)

### Extract the genes in each supercluster
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
names(Genes_in_SuperClusters)=gsub("\\|.*","",names(Genes_in_SuperClusters))
length(Genes_in_SuperClusters)

### Analyze functional similarity
# prepare the gene clusters
cluster_list1=
  lapply(gene_lists, 
         function(x) clusterProfiler::bitr(x, fromType="SYMBOL", OrgDb="org.Hs.eg.db", toType="ENTREZID") %>% .[["ENTREZID"]])
names(cluster_list1)=names(gene_lists)
cluster_list2=
  lapply(Genes_in_SuperClusters, 
         function(x) clusterProfiler::bitr(x, fromType="SYMBOL", OrgDb="org.Hs.eg.db", toType="ENTREZID") %>% .[["ENTREZID"]])
names(cluster_list2)=names(Genes_in_SuperClusters)
# analyze functional similarity with mclusterSim
d=GOSemSim::godata('org.Hs.eg.db', ont="BP")
RESULTs=RESULT_NAMES=c()
for (i in 1:length(cluster_list1)) {
  cluster_list1_=cluster_list1[[i]]
  for (j in 1:length(cluster_list2)) {
    cluster_list2_=cluster_list2[[j]]
    res=GOSemSim::mclusterSim(list(a=cluster_list1_, b=cluster_list2_), semData=d, measure="Wang", combine="BMA") %>% 
      unlist() %>% Reduce(c,.) %>% .[.!=1] %>% unique(.)
    name_=paste0(names(cluster_list1)[i], "___", names(cluster_list2)[j])
    RESULTs=c(RESULTs, res)
    RESULT_NAMES=c(RESULT_NAMES, name_)
  }
}
names(RESULTs)=RESULT_NAMES
saveRDS(RESULTs, "~/Project_PBMCage/Tempt_RDS/GOSuperClustersGenes.vs.WGCNAmoduleGenes_mClusterSimilarity.rds")

# plot the similarity
mClusterSim=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClustersGenes.vs.WGCNAmoduleGenes_mClusterSimilarity.rds")
mClusterSim_df=mClusterSim %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pair") %>%
  mutate(supercluster=gsub("___.*","",pair),
         module=gsub(".*___","",pair)) %>%
  rename(similarity=".") %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="supercluster", values_from="similarity") %>%
  tibble::column_to_rownames("module")

# reorder
terms_list=c("greenyellow"="tRNA modification",
             "salmon"="chromosome localization/mitotic checkpoint",
             "purple"="cytoplasmic translation",
             "yellow"="oxidative phosphorylation",
             "midnightblue"="regulation of hemopoiesis",
             "green"="regulation of T cell differentiation",
             "turquoise"="defense response",
             "brown"="leukocyte-mediated cytotoxicity",
             "black"="antigen processing and presentation",
             "red"="hemostasis",
             "blue"="humoral immunity")
mClusterSim_df=mClusterSim_df[names(table(rownames(mClusterSim_df)))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)],
                              unname(terms_list)]
range(mClusterSim_df) # check
col_fun_up=circlize::colorRamp2(c(0,0.5,1), c("dodgerblue3", "white", "brown3"))
ht_gene_sim=
  Heatmap(mClusterSim_df,
          name=" ",
          show_heatmap_legend=TRUE, 
          heatmap_legend_param=list(
            title="Similarity",
            grid_width=unit(0.2, "cm"), legend_height=unit(4, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontsize=10, fontface="plain")
          ),
          col=col_fun_up,
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          cluster_columns=F, cluster_rows=F, 
          show_column_dend=F, show_row_dend=F, column_dend_side="top", row_dend_side="right",
          row_names_gp=gpar(fontsize=10), show_row_names=T, row_names_side="right", 
          column_title="", column_title_gp=gpar(fontsize=11),
          show_column_names=T, 
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(4, "mm")*nrow(mClusterSim_df), 
          width=unit(4, "mm")*ncol(mClusterSim_df))

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_gene_sim=draw(ht_gene_sim, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/GOSuperClustersGenes.vs.WGCNAmoduleGenes_mClusterSimilarity.pdf", height=5, width=3.5)
ht_gene_sim
dev.off()

# plot the scaled similarity
mClusterSim=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClustersGenes.vs.WGCNAmoduleGenes_mClusterSimilarity.rds")
mClusterSim_df=mClusterSim %>%
  as.data.frame() %>%
  tibble::rownames_to_column("pair") %>%
  mutate(supercluster=gsub("___.*","",pair),
         module=gsub(".*___","",pair)) %>%
  rename(similarity=".") %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="supercluster", values_from="similarity") %>%
  tibble::column_to_rownames("module")
mClusterSim_df=mClusterSim_df[names(table(rownames(mClusterSim_df)))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)],
                              unname(terms_list)]
mClusterSim_df=scale(mClusterSim_df, scale=T, center=T)
range(mClusterSim_df) # check
col_fun_up=circlize::colorRamp2(c(-3,0,3), c("dodgerblue3", "white", "brown3"))
ht_gene_sim=
  Heatmap(mClusterSim_df,
          name=" ",
          show_heatmap_legend=TRUE, 
          heatmap_legend_param=list(
            title="Scaled similarity",
            grid_width=unit(0.2, "cm"), legend_height=unit(4, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontsize=10, fontface="plain")
          ),
          col=col_fun_up,
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          cluster_columns=F, cluster_rows=F, 
          show_column_dend=F, show_row_dend=F, column_dend_side="top", row_dend_side="right",
          row_names_gp=gpar(fontsize=10), show_row_names=T, row_names_side="right", 
          column_title="", column_title_gp=gpar(fontsize=11),
          show_column_names=T, 
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(4, "mm")*nrow(mClusterSim_df), 
          width=unit(4, "mm")*ncol(mClusterSim_df))

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_gene_sim=draw(ht_gene_sim, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/GOSuperClustersGenes.vs.WGCNAmoduleGenes_mClusterSimilarity_scaled.pdf", height=5, width=3.5)
ht_gene_sim
dev.off()

#####################################



### Enrich the WGCNA gene modules in each celltype at the inter level
#####################################
###

### Enrich the modules
Module_genes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[3]]
celltypes=unique(Module_genes$wgcna_slot)
go_objs_allcelltypes=list()
for (i in 1:length(celltypes)) {
  Module_genes_subset=Module_genes %>% subset(wgcna_slot==celltypes[i])
  gene_lists=split(Module_genes_subset$gene_name, Module_genes_subset$module) %>% lapply(., function(x) x[!duplicated(x)])
  gene_lists_full=lapply(gene_lists, length) %>% unlist() %>% .[.!=0]
  gene_lists=gene_lists[names(gene_lists) %in% names(gene_lists_full)]
  go_objs=lapply(gene_lists, function(x) clusterProfiler::enrichGO(x, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL"))
  names(go_objs)=names(gene_lists)
  go_objs_allcelltypes[[i]]=go_objs
}
names(go_objs_allcelltypes)=celltypes
saveRDS(go_objs_allcelltypes, "~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")

#####################################



### Plot the trait-correlation and enriched terms of CD4T.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Tnaive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="RNA/protein metabolism",
                    "red"="small molecule catabolic process",
                    "green"="homeostasis of number of cells",
                    "brown"="cytoplasmic translation",
                    "yellow"="cellular respiration")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD4T.Tnaive")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD4T.naive",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD4T.naive.pdf", width=2.5, height=3)
whole_ht
dev.off()

### Plot the modules of CD4T.naive
kMEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[4]]
kMEs_subset=kMEs %>% subset(wgcna_slot=="CD4T.Tnaive" & module %in% unique(meaningful_go_results$module)) %>% group_by(module) %>%
  slice_max(kME, n=10)

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Tcm
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Tcm"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="RNA splicing",
                    "turquoise"="RNA catabolic process",
                    "red"="MHCII complex assembly",
                    "brown"="cytoplasmic translation",
                    "yellow"="T cell activation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD4T.Tcm")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD4T.Tcm",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD4T.Tcm.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Tem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Tem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("brown"="homeostasis of number of cells",
                    "yellow"="cytoplasmic translation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD4T.Tem")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD4T.Tem",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD4T.Tem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Treg
#####################################
###

go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Treg"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################



### Plot the trait-correlation and enriched terms of CD8T.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.naive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="cytoplasmic translation",
                    "turquoise"="RNA/protein metabolism",
                    "red"="DNA repair",
                    "brown"="autophagy",
                    "yellow"="protein modification")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD8T.naive")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD8T.naive",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD8T.naive.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.Tcm
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.Tcm"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="RNA catabolic process",
                    "black"="cytoplasmic translation",
                    "red"="RNA splicing",
                    "green"="protein modification",
                    "magenta"="TCR signaling",
                    "purple"="lymphocyte mediated immunity")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD8T.Tcm")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD8T.Tcm",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD8T.Tcm.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.Tem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.Tem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="lymphocyte mediated immunity",
                    "turquoise"="RNA catabolic process",
                    "black"="lymphocyte activation",
                    "green"="TCR signaling",
                    "brown"="T cell activation",
                    "greenyellow"="MHCII complex assembly",
                    "purple"="cytoplasmic translation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD8T.Tem")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="CD8T.Tem",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_CD8T.Tem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.CTL
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.CTL"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################



### Plot the trait-correlation and enriched terms of B.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.naive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="RNA/protein metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.naive")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)),
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="B.naive",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_B.naive.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.inter
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.inter"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("red"="cytoplasmic translation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.inter")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)),
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="B.inter",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_B.inter.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.trans
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.trans"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="cellular localization")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.trans")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)),
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="B.inter",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_B.trans.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.mem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.mem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("green"="RNA/protein metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.mem")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)),
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="B.mem",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_B.mem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of B.pbpc
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.pbpc"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="energy metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.pbpc")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="B.pbpc",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_B.pbpc.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of NK.CD56dim
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.CD56dim"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="RNA/protein metabolism",
                    "turquoise"="protein modification",
                    "green"="cytoplasmic translation",
                    "yellow"="leukocyte mediated immunity",
                    "pink"="cellular response",
                    "lightcyan"="protein catabolic process")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="NK.CD56dim")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)),
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="NK.CD56dim",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_NK.CD56dim.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of NK.CD56hi
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.CD56hi"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="leukocyte activation",
                    "turquoise"="ATP synthesis",
                    "red"="nucleoside metabolic process",
                    "green"="MHCII complex assembly",
                    "brown"="chromatin organization",
                    "pink"="RNA/protein metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="NK.CD56hi")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="NK.CD56hi",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_NK.CD56hi.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of NK.prolif
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.prolif"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="inteferon production")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="NK.prolif")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="NK.prolif",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*length(cor_pval_fdr_arranged),
          width=unit(7,"mm")*1)

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_NK.prolif.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of Mono.classical
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["Mono.classical"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################




### Plot the trait-correlation and enriched terms of Mono.inter
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["Mono.inter"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################




### Plot the trait-correlation and enriched terms of Mono.nonclassical
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["Mono.nonclassical"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="innate immune activation",
                    "red"="RNA/protein metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="Mono.nonclassical")
cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))

cor_pval_fdr_arranged=cor_pval_fdr %>%
  dplyr::select(traits, color, cor) %>% 
  tidyr::pivot_wider(names_from="color", values_from="cor") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]
all_cell_fdr=cor_pval_fdr %>%
  dplyr::select(traits, color, pval) %>% 
  mutate(pval=ifelse(pval<0.0001,"****",
                     ifelse(pval>=0.0001 & pval<0.001, "***",
                            ifelse(pval>=0.001 & pval<0.01, "**",
                                   ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
  tidyr::pivot_wider(names_from="color", values_from="pval") %>%
  tibble::column_to_rownames("traits") %>%
  .[, unique(meaningful_go_results$module)]

library(ComplexHeatmap)
ht_list=
  Heatmap(cor_pval_fdr_arranged,
          name=" ",
          show_heatmap_legend=F,
          heatmap_legend_param=list(
            title="Correlation",
            legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
                                           annotation_name_rot=90),
          col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
          cluster_columns=F, cluster_rows=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=9),
          column_names_gp=gpar(fontsize=9),
          row_title="Mono.nonclassical",
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(as.matrix(all_cell_fdr), i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation.heatmap_Mono.nonclassical.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################