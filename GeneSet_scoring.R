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

###
#####################################
### Score the geneset of ROS
###

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES=read.gmt("~/Project_PBMCage/REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES.v2023.2.Hs.gmt")
ETHANOL_METABOLISM_PRODUCTION_OF_ROS_BY_CYP2E1=read.gmt("~/Project_PBMCage/WP_ETHANOL_METABOLISM_PRODUCTION_OF_ROS_BY_CYP2E1.v2023.2.Hs.gmt")
CTRL_VS_ROS_INHIBITOR_TREATED_DC_DN=read.gmt("~/Project_PBMCage/GSE20727_CTRL_VS_ROS_INHIBITOR_TREATED_DC_DN.v2023.2.Hs.gmt")
CTRL_VS_ROS_INHIBITOR_TREATED_DC_UP=read.gmt("~/Project_PBMCage/GSE20727_CTRL_VS_ROS_INHIBITOR_TREATED_DC_UP.v2023.2.Hs.gmt")
H2O2_VS_ROS_INHIBITOR_TREATED_DC_DN=read.gmt("~/Project_PBMCage/GSE20727_H2O2_VS_ROS_INHIBITOR_TREATED_DC_DN.v2023.2.Hs.gmt")
H2O2_VS_ROS_INHIBITOR_TREATED_DC_UP=read.gmt("~/Project_PBMCage/GSE20727_H2O2_VS_ROS_INHIBITOR_TREATED_DC_UP.v2023.2.Hs.gmt")

# score the geneset of ROS
ROS_scoring=list(NOX2_ROS=REACTOME_ROS_AND_RNS_PRODUCTION_IN_PHAGOCYTES$gene,
                 Metabolism_ROS=ETHANOL_METABOLISM_PRODUCTION_OF_ROS_BY_CYP2E1$gene,
                 DN_postDPI=CTRL_VS_ROS_INHIBITOR_TREATED_DC_DN$gene,
                 UP_postDPI=CTRL_VS_ROS_INHIBITOR_TREATED_DC_UP$gene,
                 DN_DPI_comparedto_H2O2=H2O2_VS_ROS_INHIBITOR_TREATED_DC_DN$gene,
                 UP_DPI_comparedto_H2O2=H2O2_VS_ROS_INHIBITOR_TREATED_DC_UP$gene)
library(UCell)
library(patchwork)
object=AddModuleScore_UCell(object, features=ROS_scoring)
mydatadf=object[[]]

# add the agecut metadata
AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}
agecut_df=object[[]]
agecut_df$age=as.factor(agecut_df$age)
agecut_df$agecut=agecut_df$age
for (i in 1:length(AGECUT)) {
  agecut_df=agecut_df %>%
    mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
}
table(rownames(agecut_df)==colnames(object)) # check
object=AddMetaData(object, metadata=agecut_df$agecut, col.name="agecut")
mydatadf=object[[]]
saveRDS(mydatadf, "~/Project_PBMCage/Tempt_RDS/Scoring_Ucell_ROS.rds")


### Plot the scores related to ROS for each rough celltype
library(SCpubr)
score_list=colnames(object[[]])[49:54]
score_list_name=c("Score Assco. w/ NOX2 ROS",
                  "Score Assoc. w/ mtROS",
                  "Score Pos.Assoc. w/ ROS",
                  "Score Neg.Assoc. w/ ROS",
                  "Score Pos.Assoc. w/ H2O2",
                  "Score Neg.Assoc. w/ H2O2")
plot_list_=list(plot_list_1=list(),
                plot_list_2=list(),
                plot_list_3=list(),
                plot_list_4=list(),
                plot_list_5=list(),
                plot_list_6=list())

for (s in 1:length(score_list)) {
  plot_list_[[s]]=list()
  for (i in 1:length(unique(object$CD8T_NK_subsets2_rough))) {
    plot_list_[[s]][[i]]=
      do_BoxPlot(subset(object, CD8T_NK_subsets2_rough==unique(object$CD8T_NK_subsets2_rough)[i]), 
                 feature=score_list[s], 
                 use_silhouette=TRUE, 
                 group.by="age", 
                 font.size=8, 
                 plot.title=unique(object$CD8T_NK_subsets2_rough)[i], 
                 legend.position="none",
                 ylab=score_list_name[s], 
                 xlab="", 
                 axis.text.x.angle=90)
  }
}

pdf("~/Project_PBMCage/Plots/Scoring_RelatedToROS.pdf", height=5, width=6)
for (s in 1:length(score_list)) {
  for (p in plot_list_[[s]]) plot(p)
}
dev.off()


### Plot the scores related to ROS for all the cells together
score_list=colnames(object[[]])[49:54]
score_list_name=c("Score Assco. w/ NOX2 ROS",
                  "Score Assoc. w/ mtROS",
                  "Score Pos.Assoc. w/ ROS",
                  "Score Neg.Assoc. w/ ROS",
                  "Score Pos.Assoc. w/ H2O2",
                  "Score Neg.Assoc. w/ H2O2")
plot_list_=list()
for (i in 1:length(score_list)) {
  plot_list_[[i]]=
    do_BoxPlot(object, 
               feature=score_list[i], 
               use_silhouette=TRUE, 
               group.by="age", 
               font.size=8, 
               plot.title=score_list_name[i], 
               legend.position="none", 
               ylab=score_list_name[i], 
               xlab="", 
               axis.text.x.angle=90)
}

pdf("~/Project_PBMCage/Plots/Scoring_RelatedToROS_AllCells.pdf", height=5, width=6)
for (p in plot_list_) plot(p)
dev.off()


### Plot the scores related to ROS for each rough celltype split.by agecut
library(SCpubr)
score_list=colnames(object[[]])[49:54]
score_list_name=c("Score Assco. w/ NOX2 ROS",
                  "Score Assoc. w/ mtROS",
                  "Score Pos.Assoc. w/ ROS",
                  "Score Neg.Assoc. w/ ROS",
                  "Score Pos.Assoc. w/ H2O2",
                  "Score Neg.Assoc. w/ H2O2")
plot_list_=list(plot_list_1=list(),
                plot_list_2=list(),
                plot_list_3=list(),
                plot_list_4=list(),
                plot_list_5=list(),
                plot_list_6=list())

for (s in 1:length(score_list)) {
  plot_list_[[s]]=list()
  for (i in 1:length(unique(object$CD8T_NK_subsets2_rough))) {
    plot_list_[[s]][[i]]=
      do_BoxPlot(subset(object, CD8T_NK_subsets2_rough==unique(object$CD8T_NK_subsets2_rough)[i]), 
                 feature=score_list[s], 
                 use_silhouette=TRUE, 
                 group.by="agecut", 
                 font.size=8, 
                 plot.title=unique(object$CD8T_NK_subsets2_rough)[i], 
                 legend.position="none",
                 ylab=score_list_name[s], 
                 xlab="", 
                 axis.text.x.angle=90)
  }
}

pdf("~/Project_PBMCage/Plots/Scoring_Agecut_RelatedToROS.pdf", height=5, width=6)
for (s in 1:length(score_list)) {
  for (p in plot_list_[[s]]) plot(p)
}
dev.off()


### Plot the scores related to ROS for all the cells together split.by agecut
score_list=colnames(object[[]])[49:54]
score_list_name=c("Score Assco. w/ NOX2 ROS",
                  "Score Assoc. w/ mtROS",
                  "Score Pos.Assoc. w/ ROS",
                  "Score Neg.Assoc. w/ ROS",
                  "Score Pos.Assoc. w/ H2O2",
                  "Score Neg.Assoc. w/ H2O2")
plot_list_=list()
for (i in 1:length(score_list)) {
  plot_list_[[i]]=
    do_BoxPlot(object, 
               feature=score_list[i], 
               use_silhouette=TRUE, 
               group.by="agecut", 
               font.size=8, 
               plot.title=score_list_name[i], 
               legend.position="none", 
               ylab=score_list_name[i], 
               xlab="", 
               axis.text.x.angle=90)
}

pdf("~/Project_PBMCage/Plots/Scoring_Agecut_RelatedToROS_AllCells.pdf", height=5, width=6)
for (p in plot_list_) plot(p)
dev.off()

### Plot the scores related to ROS for each rough celltype, add correlation and equation
mydf=object[[]][,c(8,10,46,47,48,49,50,51,52,53,54,55)]
score_list=colnames(object[[]])[49:54]
score_list_name=c("Score Assco. w/ NOX2 ROS",
                  "Score Assoc. w/ mtROS",
                  "Score Pos.Assoc. w/ ROS",
                  "Score Neg.Assoc. w/ ROS",
                  "Score Pos.Assoc. w/ H2O2",
                  "Score Neg.Assoc. w/ H2O2")
plot_list_=list(plot_list_1=list(),
                plot_list_2=list(),
                plot_list_3=list(),
                plot_list_4=list(),
                plot_list_5=list(),
                plot_list_6=list())

for (s in 1:length(score_list)) {
  plot_list_[[s]]=list()
  for (i in 1:length(unique(object$CD8T_NK_subsets2_rough))) {
    tempt_df=subset(mydf, CD8T_NK_subsets2_rough==unique(object$CD8T_NK_subsets2_rough)[i])
    tempt_df$age=as.numeric(tempt_df$age)
    plot_list_[[s]][[i]]=
      ggplot(tempt_df, aes_string(x="age", y=score_list[s])) +
      geom_point(size=1, shape=1, color="gray90") +
      theme_classic() +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            axis.title.y=element_text(size=9),
            legend.title=element_text(size=9),
            legend.text=element_text(size=8)) +
      labs(title=unique(object$CD8T_NK_subsets2_rough)[i])+
      ylab(score_list_name[s]) +
      geom_smooth(method='lm', show.legend=FALSE, linewidth=2, fill="darkblue") +
      stat_cor(show.legend=FALSE, size=3.5) +
      stat_regline_equation(size=3.5, label.x=55)
  }
}

pdf("~/Project_PBMCage/Plots/Scoring_RelatedToROS_statcor.pdf", height=5, width=6)
for (s in 1:length(score_list)) {
  for (p in plot_list_[[s]]) plot(p)
}
dev.off()

for (s in 1:length(score_list)) {
  pdf(paste0("~/Project_PBMCage/Plots/Scoring_RelatedToROS_statcor_",s,".pdf"), height=5, width=6)
  for (p in plot_list_[[s]]) plot(p)
  dev.off()
}

#####################################



###
#####################################
### Score the T cell recombination genes
###

### Load
THEObj=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
assigned_types=unique(THEObj$CD8T_NK_subsets2)[!grepl("lin_infidel|mix", unique(THEObj$CD8T_NK_subsets2))]
THEObj=subset(THEObj, CD8T_NK_subsets2 %in% assigned_types)
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

### Make a function for plotting on ages
geneset_scoring_forAge=function(dataframe, celltype_as_LegendName, scale) {
  df_full_list=list()
  CD8T_NK_subsets_celltypes=unique(dataframe$CD8T_NK_subsets2)
  for (i in 1:length(CD8T_NK_subsets_celltypes)) {
    matrix_=subset(dataframe, CD8T_NK_subsets2==CD8T_NK_subsets_celltypes[i])
    age_not_there=unique(dataframe$age)[!(unique(dataframe$age) %in% unique(matrix_$age))]
    if (length(age_not_there)!=0) {
      matrix_foragenotthere=data.frame(CD8T_NK_subsets2=unique(matrix_$CD8T_NK_subsets2),
                                       CD8T_NK_subsets2_inter=unique(matrix_$CD8T_NK_subsets2_inter),
                                       CD8T_NK_subsets2_rough=unique(matrix_$CD8T_NK_subsets2_rough),
                                       age=age_not_there,
                                       BCL11B=rep(NA, length(age_not_there)),
                                       LEF1=rep(NA, length(age_not_there)),
                                       LIG4=rep(NA, length(age_not_there)),
                                       PRKDC=rep(NA, length(age_not_there)),
                                       TCF7=rep(NA, length(age_not_there)))
    } else {
      matrix_foragenotthere=data.frame()
    }
    matrix_=data.table::rbindlist(list(matrix_, matrix_foragenotthere))
    matrix_=matrix_ %>% arrange(age)
    age_=matrix_$age
    matrix_=t(matrix_[,c("BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7")])
    matrix_=as.data.frame(matrix_)
    colnames(matrix_)=age_
    rownames(matrix_)=c("BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7")
    df_full_list[[i]]=matrix_
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    col_fun=circlize::colorRamp2(c(-1,0,1), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-rowMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    col_fun=circlize::colorRamp2(c(-1,0,1), c("#00007E", "white", "#7E0000"))}
  
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
                name=paste0("Mean Exp. of ", celltype_as_LegendName),
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

### Score the recombination genes in CD8 T cell subtypes
CD8T=subset(THEObj, CD8T_NK_subsets2_rough=="CD8T cells")
Idents(CD8T)="age"
Recombination_genes=FetchData(CD8T, 
                              vars=c("CD8T_NK_subsets2","CD8T_NK_subsets2_inter","CD8T_NK_subsets2_rough",
                                     "age","agecut","sex",
                                     "BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7"), 
                              layer="data")
Recombination_genes_group=Recombination_genes %>%
  group_by(CD8T_NK_subsets2, CD8T_NK_subsets2_inter, CD8T_NK_subsets2_rough, age) %>%
  summarize(BCL11B=mean(BCL11B),
            LEF1=mean(LEF1),
            LIG4=mean(LIG4),
            PRKDC=mean(PRKDC),
            TCF7=mean(TCF7)) %>%
  as.data.frame()
CD8T_recombination_plots=geneset_scoring_forAge(dataframe=Recombination_genes_group, celltype_as_LegendName="CD8 T cells", scale=F)

### Score the recombination genes in Progenitors
progenitor=subset(THEObj, CD8T_NK_subsets2_inter=="Other.progenitor")
Idents(progenitor)="age"
Recombination_genes_progenitor=FetchData(progenitor, 
                              vars=c("CD8T_NK_subsets2","CD8T_NK_subsets2_inter","CD8T_NK_subsets2_rough",
                                     "age","agecut","sex",
                                     "BCL11B", "LEF1", "LIG4", "PRKDC", "TCF7"), 
                              layer="data")
Recombination_genes_group_progenitor=Recombination_genes_progenitor %>%
  group_by(CD8T_NK_subsets2, CD8T_NK_subsets2_inter, CD8T_NK_subsets2_rough, age) %>%
  summarize(BCL11B=mean(BCL11B),
            LEF1=mean(LEF1),
            LIG4=mean(LIG4),
            PRKDC=mean(PRKDC),
            TCF7=mean(TCF7)) %>%
  as.data.frame()
progenitor_recombination_plots=geneset_scoring_forAge(dataframe=Recombination_genes_group_progenitor, celltype_as_LegendName="HSPC", scale=F)

### Plot
pdf("~/Project_PBMCage/Plots/Scoring_VDJrecombination_CD8T.pdf", width=12, height=12)
draw(CD8T_recombination_plots)
draw(progenitor_recombination_plots)
dev.off()


