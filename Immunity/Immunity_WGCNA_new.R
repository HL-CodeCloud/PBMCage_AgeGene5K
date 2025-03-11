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
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
pseudobulk_obj_subset_byage=Seurat::AggregateExpression(pseudobulk_obj, group.by=c("age","sex","Annot.rough"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
group_names=lapply(colnames(pseudobulk_obj_subset_byage), function(x) strsplit(x, split="_")[[1]])
pseudobulk_obj_subset_byage$age=lapply(group_names, function(x) x[1]) %>% unlist()
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=lapply(group_names, function(x) x[2]) %>% unlist()
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
pseudobulk_obj_subset_byage$sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$sex, c("M","F"))
age_list=names(table(pseudobulk_obj_subset_byage$age))
sex_list=c("M","F")
age.sex_list=paste0(rep(age_list, each=2), ", ", rep(sex_list, times=length(age_list)))
pseudobulk_obj_subset_byage$age.sex=paste0(pseudobulk_obj_subset_byage$age, ", ", pseudobulk_obj_subset_byage$sex)
pseudobulk_obj_subset_byage$age.sex=forcats::fct_relevel(pseudobulk_obj_subset_byage$age.sex, age.sex_list)
pseudobulk_obj_subset_byage$Annot.rough=lapply(group_names, function(x) x[3]) %>% unlist()
# normalize expr
pseudobulk_obj_subset_byage=pseudobulk_obj_subset_byage %>%
  Seurat::NormalizeData() %>%
  Seurat::FindVariableFeatures() %>%
  Seurat::ScaleData()
pseudobulk_obj_subset_byage[["RNA"]]=as(pseudobulk_obj_subset_byage[["RNA"]], Class="Assay")

### Set up with selected genes
WGCNAobj=pseudobulk_obj_subset_byage
# cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05) %>% dplyr::select(gene) %>% tibble::deframe() %>% .[!duplicated(.)]
WGCNAobj=hdWGCNA::SetupForWGCNA(
  WGCNAobj,
  group.by="Annot.rough",
  gene_select="fraction",
  fraction=0.2, 
  # features=gene_selected,
  wgcna_name="wgcna")

### Set the expression data for each wgcna slots
WGCNAobj=hdWGCNA::SetDatExpr(
  WGCNAobj,
  group_name=unique(WGCNAobj$Annot.rough),
  group.by="Annot.rough",
  assay="RNA",
  use_metacells=FALSE,
  slot="data")
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")
print("Set expression for MetaCell analysis is done!")

### Choose the power
selected_power=0
WGCNAobj=hdWGCNA::TestSoftPowers(
  WGCNAobj,
  powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
plot_list=hdWGCNA::PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
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
  setwd("~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/")
  enableWGCNAThreads(nThreads=8)
  WGCNAobj=hdWGCNA::ConstructNetwork(
    WGCNAobj,
    soft_power=selected_power,
    corType="pearson",
    networkType="signed",
    TOMType="signed",
    overwrite_tom=TRUE,
    tom_outdir="TOM",
    tom_name="rough"
  )
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")
} else {print("No selected_power found.")}

#####################################




### ANALYZE CD4T.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.naive.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.naive.rds")
} else {print("No selected_power found.")}

### Plot all gene expr in this subpopulation
library(ComplexHeatmap)
WGCNAobj=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.naive.rds")
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
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Tem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tem.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD4T.Tcm
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Tcm")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tcm.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Tcm.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD4T.Treg
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD4T.Treg")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Treg.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD4T.Treg.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.naive.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.naive.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.Tcm
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.Tcm")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tcm.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tcm.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.Tem
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.Tem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tem.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.Tem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE CD8T.CTL
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="CD8T.CTL")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.CTL.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_CD8T.CTL.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.naive
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.naive")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.naive.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.naive.rds")
} else {print("No selected_power found.")}
# only grey module found.

#####################################



### ANALYZE B.mem
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.mem")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.mem.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.mem.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.inter
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.inter")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.1 & celltypes=="B.inter") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only the grey module
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.inter.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.inter.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.trans
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.trans")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="B.trans") %>% dplyr::select(gene) %>% tibble::deframe()
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.trans.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.trans.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE B.pbpc
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="B.pbpc")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.05 & celltypes=="B.pbpc") %>% dplyr::select(gene) %>% tibble::deframe()
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.pbpc.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_B.pbpc.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.CD56dim
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.CD56dim")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56dim.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56dim.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.CD56hi
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.CD56hi")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56hi.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.CD56hi.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE NK.prolif
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="NK.prolif")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
gene_selected=cor_df %>% subset(analysis=="all" & rho_pval<0.1 & celltypes=="NK.prolif") %>% dplyr::select(gene) %>% tibble::deframe() # if rho_pval<0.05, will be only the grey module
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.prolif.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_NK.prolif.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE Mono.classical
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="Mono.classical")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_Mono.classical.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_Mono.classical.rds")
} else {print("No selected_power found.")}

#####################################



### ANALYZE Mono.nonclassical
#####################################
###

### Create pseudobulk obj
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
pseudobulk_obj_subset=subset(pseudobulk_obj, Annot.inter=="Mono.nonclassical")
pseudobulk_obj_subset_byage=AggregateExpression(pseudobulk_obj_subset, group.by=c("age","sex"), assays="RNA", return.seurat=T)
# arrange metadata
pseudobulk_obj_subset_byage$orig.ident="Aggregate"
pseudobulk_obj_subset_byage$age=gsub("_.*","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$age=as.integer(as.character(pseudobulk_obj_subset_byage$age))
pseudobulk_obj_subset_byage$sex=gsub(".*_","",colnames(pseudobulk_obj_subset_byage))
pseudobulk_obj_subset_byage$sex=ifelse(pseudobulk_obj_subset_byage$sex=="Male","M","F")
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
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
saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_Mono.nonclassical.rds")
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
  saveRDS(WGCNAobj, "~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_Mono.nonclassical.rds")
} else {print("No selected_power found.")}

#####################################




### Extract the results from the WGCNA analysis at the rough level
#####################################
###

###
WGCNAobj=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/WGCNAobj_allcells_Rough.rds")

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

saveRDS(list(MEs, cor_pval_fdr, Module_genes, hb), "~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")

#####################################




### Extract the results from the WGCNA analysis at the inter level
#####################################
###

###
files_=list.files("~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/") %>% .[grepl("^WGCNAobj_[^a]",.)]

ME_list=COR_PVAL_FDR_list=Module_kME_list=Module_totkME_list=list()
for (i in 1:length(files_)) {
  WGCNAobj=readRDS(paste0("~/Project_PBMCage/Immunity/Tempt_RDS/wgcna_results/",files_[i]))

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

saveRDS(list(ME_list_total, COR_PVAL_FDR_list_total, Module_kME_list_total, Module_totkME_list_total), "~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")
rm(list=ls());gc()

#####################################




### Analyze the correlation results at the rough level
#####################################
###

### Enrich the modules
Module_genes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[3]]
gene_lists=split(Module_genes$gene_name, Module_genes$module) %>% lapply(., function(x) x[!duplicated(x)])
gene_lists=gene_lists[names(gene_lists)!="grey"]
go_objs=lapply(gene_lists, function(x) clusterProfiler::enrichGO(x, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL"))
names(go_objs)=names(gene_lists)
saveRDS(go_objs, "~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")

### Plot the correlation only for the modules with enriched biological meanings
go_objs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results_GOenrich.rds")
meaningful_go_objs=go_objs[lapply(go_objs, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the biological meanings
terms_list=c("salmon"="cytoplasmic translation",
             "black"="ncRNA processing",
             "green"="secondary metabolic process",
             "magenta"="mitochondrial gene expression",
             "grey60"="mitochondrial RNA metabolic process",
             "tan"="actin filament organization",
             "maroon"="cell-substrate junction assembly",
             "yellowgreen"="axoneme assembly",
             "skyblue"="chromosome segregation",
             "saddlebrown"="DNA replication",
             "greenyellow"="cytoskeleton organization",
             "lightcyan1"="meiosis",
             "lavenderblush3"="meiotic nuclear division",
             "lightcyan"="Golgi vesicle transport",
             "ivory"="transmembrane transport",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "navajowhite2"="immune response signaling",
             "brown"="regulation of T cell differentiation",
             "pink"="leukocyte-mediated cytotoxicity",
             "red"="T/NK-mediated immunity",
             "yellow"="regulation of hemopoiesis",
             "plum1"="MHC complex assembly")
# arrange the trait-correlation data
cor_pval_fdr=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[2]]
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
               function(idx) cor_pval_fdr_arranged[[idx]] %>% 
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
                           v=pindex(as.matrix(all_cell_fdr[[idx]]), i, j)
                           grid.text(v, x, y, rot=0, just="center",
                                     gp=gpar(fontsize=10, fontface="bold"))
                         },
                         row_title_gp=gpar(fontsize=10), row_title_rot=0, row_title_side="left",
                         height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged[[idx]]),
                         width=unit(8,"mm")*ncol(cor_pval_fdr_arranged[[idx]])))

ht_list=Reduce(`%v%`, ht_list)
ht_opt$HEATMAP_LEGEND_PADDING=unit(5,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_rough.pdf", width=10, height=6.1)
whole_ht
dev.off()

#####################################



### Get the hub genes in each module at the rough level
#####################################
###

library(dplyr)

MEs_hub=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
# previously determined terms
terms_list=c("salmon"="cytoplasmic translation",
             "black"="ncRNA processing",
             "green"="secondary metabolic process",
             "magenta"="mitochondrial gene expression",
             "grey60"="mitochondrial RNA metabolic process",
             "tan"="actin filament organization",
             "maroon"="cell-substrate junction assembly",
             "yellowgreen"="axoneme assembly",
             "skyblue"="chromosome segregation",
             "saddlebrown"="DNA replication",
             "greenyellow"="cytoskeleton organization",
             "lightcyan1"="meiosis",
             "lavenderblush3"="meiotic nuclear division",
             "lightcyan"="Golgi vesicle transport",
             "ivory"="transmembrane transport",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "navajowhite2"="immune response signaling",
             "brown"="T cell differentiation",
             "pink"="leukocyte-mediated cytotoxicity",
             "red"="T/NK-mediated immunity",
             "yellow"="regulation of hemopoiesis",
             "plum1"="MHC complex assembly")
MEs_hub_subset=MEs_hub %>% subset(module %in% names(terms_list)) %>%
  group_by(module) %>% slice_max(kME, n=15) %>%
  left_join(data.frame(module=names(terms_list), term=terms_list), by="module")
top15genes_per_module=split(MEs_hub_subset$gene_name, MEs_hub_subset$term)

saveRDS(top15genes_per_module, "~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_hubGenes.rds")

#####################################



### Plot the expr of gene modules in NK cells
#####################################
###

### Load the expr data
pseudobulk_data=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")

### Load the gene modules
terms_list=c("salmon"="cytoplasmic translation",
             "black"="ncRNA processing",
             "green"="secondary metabolic process",
             "magenta"="mitochondrial gene expression",
             "grey60"="mitochondrial RNA metabolic process",
             "tan"="actin filament organization",
             "maroon"="cell-substrate junction assembly",
             "yellowgreen"="axoneme assembly",
             "skyblue"="chromosome segregation",
             "saddlebrown"="DNA replication",
             "greenyellow"="cytoskeleton organization",
             "lightcyan1"="meiosis",
             "lavenderblush3"="meiotic nuclear division",
             "lightcyan"="Golgi vesicle transport",
             "ivory"="transmembrane transport",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "navajowhite2"="immune response signaling",
             "brown"="T cell differentiation",
             "pink"="leukocyte-mediated cytotoxicity",
             "red"="T/NK-mediated immunity",
             "yellow"="regulation of hemopoiesis",
             "plum1"="MHC complex assembly")
genemodules=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]] %>%
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
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
  mutate(term=ifelse(term=="T cell differentiation","regulation of T cell differentiation",term)) %>%
  mutate(term=stringr::str_wrap(term, 20))

### Plot
terms_list=gsub("T cell differentiation","regulation of T cell differentiation",terms_list)
EXPR_All_Combo$term=forcats::fct_relevel(EXPR_All_Combo$term, stringr::str_wrap(terms_list, 20))
plot_=
  ggplot(EXPR_All_Combo, aes(x=age, y=merged, color=sex)) +
  facet_wrap(~term, scales="free", ncol=6, labeller=as_labeller(c())) +
  ggrastr::geom_point_rast(size=0.1, shape=1, alpha=0.1) +
  geom_smooth(aes(group=sex), method="loess", show.legend=T, fill="grey93", linewidth=0.25) +
  ggpubr::stat_cor(aes(group=sex), method="spearman", show.legend=F, size=3.5) +
  theme_minimal() +
  scale_x_continuous(limits=c(19,100), breaks=c(20,40,60,80,100)) +
  scale_y_continuous(n.breaks=3) +
  labs(y="Avg. expression", subtitle="Terekhova et al. (2023) dataset") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        plot.subtitle=element_text(size=11),
        axis.title=element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.925,0.1),
        legend.direction="horizontal") +
  guides(color=guide_legend(override.aes=list(fill="transparent")))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WCCNA_genemoduleExpr_NK.pdf", height=6, width=14.5)
plot_
dev.off()

#####################################



### Enrich the WGCNA gene modules in each celltype at the inter level
#####################################
###

### Enrich the modules
Module_genes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[3]]
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
saveRDS(go_objs_allcelltypes, "~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")

#####################################



### Plot the trait-correlation and enriched terms of CD4T.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.naive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="lymphocyte mediated immunity",
                    "pink"="signal release")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD4T.naive")
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD4T.naive.pdf", width=2.5, height=3)
whole_ht
dev.off()

### Plot the modules of CD4T.naive
kMEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[4]]
kMEs_subset=kMEs %>% subset(wgcna_slot=="CD4T.naive" & module %in% unique(meaningful_go_results$module)) %>% group_by(module) %>%
  slice_max(kME, n=10)

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Tcm
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Tcm"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="cellular detoxification",
                    "yellow"="T cell activation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity/WGCNA_correlation.heatmap_CD4T.Tcm.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Tem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Tem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="MHCII complex assembly",
                    "brown"="T cell differentiation",
                    "purple"="cell junction organization")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD4T.Tem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD4T.Treg
#####################################
###

go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD4T.Treg"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="energy metabolism",
                    "blue"="cytoplasmic translation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD4T.Treg")
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
          row_title="CD4T.Treg",
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD4T.Treg.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.naive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("tan"="chemokine-mediated signaling pathway",
                    "orangered4"="cell migration",
                    "cyan"="cell junction assembly",
                    "white"="Wnt signaling pathway",
                    "darkred"="amino acid transport",
                    "darkgreen"="integrin-mediated signaling pathway",
                    "violet"="maintenance of location")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD8T.naive.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.Tcm
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.Tcm"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("green"="transporter activity",
                    "yellow"="electron transport chain",
                    "brown"="T cell activation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD8T.Tcm.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.Tem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.Tem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="actin filament organization",
                    "red"="MHCII complex assembly",
                    "brown"="T cell activation")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD8T.Tem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of CD8T.CTL
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["CD8T.CTL"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="kinase activity")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="CD8T.CTL")
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
          row_title="CD8T.CTL",
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_CD8T.CTL.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.naive
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.naive"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c() # no good term found
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

# ## Get the correlation of modules with only the meaning ones
# cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
# cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="B.naive")
# cor_pval_fdr=cor_pval_fdr %>% subset(color %in% unique(meaningful_go_results$module))
# 
# cor_pval_fdr_arranged=cor_pval_fdr %>%
#   dplyr::select(traits, color, cor) %>% 
#   tidyr::pivot_wider(names_from="color", values_from="cor") %>%
#   tibble::column_to_rownames("traits") %>%
#   .[, unique(meaningful_go_results$module)]
# all_cell_fdr=cor_pval_fdr %>%
#   dplyr::select(traits, color, pval) %>% 
#   mutate(pval=ifelse(pval<0.0001,"****",
#                      ifelse(pval>=0.0001 & pval<0.001, "***",
#                             ifelse(pval>=0.001 & pval<0.01, "**",
#                                    ifelse(pval>=0.01 & pval<0.05, "*", ""))))) %>%
#   tidyr::pivot_wider(names_from="color", values_from="pval") %>%
#   tibble::column_to_rownames("traits") %>%
#   .[, unique(meaningful_go_results$module)]
# 
# library(ComplexHeatmap)
# ht_list=
#   Heatmap(cor_pval_fdr_arranged,
#           name=" ",
#           show_heatmap_legend=F,
#           heatmap_legend_param=list(
#             title="Correlation",
#             legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
#             direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
#             title_gp=gpar(fontface="plain", fontsize=10)
#           ),
#           top_annotation=HeatmapAnnotation(terms=anno_text(enrichments_terms, just="left", location=unit(0, "npc"), show_name=FALSE, gp=gpar(fontsize=9)), 
#                                            annotation_name_rot=90),
#           col=circlize::colorRamp2(c(-1,0,1), c("dodgerblue3","white","brown3")),
#           cluster_columns=F, cluster_rows=F, show_column_dend=F,
#           row_names_gp=gpar(fontsize=9),
#           column_names_gp=gpar(fontsize=9),
#           row_title="B.naive",
#           layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
#             v=pindex(as.matrix(all_cell_fdr), i, j)
#             grid.text(v, x, y, rot=0, just="center",
#                       gp=gpar(fontsize=10, fontface="bold"))
#           },
#           row_title_gp=gpar(fontsize=10), row_title_rot=90, row_title_side="left",
#           height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
#           width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))
# 
# ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
# whole_ht=draw(ht_list, heatmap_legend_side="right")
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_B.naive.pdf", width=2.5, height=3)
# whole_ht
# dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.inter
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.inter"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################



### Plot the trait-correlation and enriched terms of B.trans
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.trans"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="RNA metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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
          height=unit(3.5,"mm")*nrow(cor_pval_fdr_arranged),
          width=unit(7,"mm")*ncol(cor_pval_fdr_arranged))

ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_B.trans.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################



### Plot the trait-correlation and enriched terms of B.mem
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.mem"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("red"="synapse structure or activity",
                    "pink"="lipid transport",
                    "brown"="RNA/protein metabolism")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_B.mem.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of B.pbpc
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["B.pbpc"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################




### Plot the trait-correlation and enriched terms of NK.CD56dim
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.CD56dim"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("blue"="protein metabolism",
                    "red"="telomere maintenance")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_NK.CD56dim.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################




### Plot the trait-correlation and enriched terms of NK.CD56hi
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.CD56hi"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################




### Plot the trait-correlation and enriched terms of NK.prolif
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["NK.prolif"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
# no module with enriched terms at p.adjust<0.05

#####################################




### Plot the trait-correlation and enriched terms of Mono.classical
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["Mono.classical"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="antigen processing",
                    "blue"="negative regulation of apoptosis")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
cor_pval_fdr=cor_pval_fdr_total %>% subset(wgcna_slot=="Mono.classical")
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
          row_title="Mono.classical",
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_Mono.classical.pdf", width=2.5, height=3)
whole_ht
dev.off()


#####################################




### Plot the trait-correlation and enriched terms of Mono.nonclassical
#####################################
###

### Manually identified the summarized enriched terms
MEs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[1]]
go_objs_allcelltypes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results_GOenrich.rds")
go_objs_subset=go_objs_allcelltypes[["Mono.nonclassical"]]
meaningful_go_objs=go_objs_subset[lapply(go_objs_subset, function(x) nrow(x@result %>% subset(p.adjust<0.05))!=0) %>% unlist()]
meaningful_go_results=meaningful_go_objs %>% lapply(., function(x) x@result %>% dplyr::select(Description, GeneRatio, pvalue, p.adjust) %>% slice_min(pvalue, n=10))
# summarize the enrichments according to the top10 terms
enrichments_terms=c("turquoise"="negative regulation of chromosome segregation",
                    "yellow"="cellular response",
                    "black"="regulation of meiotic cell cycle")
meaningful_go_results=lapply(1:length(meaningful_go_results), function(idx) meaningful_go_results[[idx]] %>% 
                               mutate(module=names(meaningful_go_results)[idx])) %>%
  data.table::rbindlist() %>%
  left_join(data.frame(module=names(enrichments_terms), marker=enrichments_terms), by="module") %>%
  subset(!is.na(marker))

## Get the correlation of modules with only the meaning ones
cor_pval_fdr_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WGCNA_Annot.inter_results.rds")[[2]]
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WGCNA_correlation.heatmap_Mono.nonclassical.pdf", width=2.5, height=3)
whole_ht
dev.off()

#####################################