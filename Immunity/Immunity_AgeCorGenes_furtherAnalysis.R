
### Enrich age-correlated genes across all celltypes (common ones)
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")

upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Rough.rds")

### Plot
library(aPEAR)
source("~/Rscripts/aPEAR2.R")
theme_=aPEAR.theme
theme_$colorBy="pvalue"; theme_$nodeSize="count"; theme_$fontSize=4; theme_$repelLabels=TRUE
# load the enrichment result
geo_cluster=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
# plot_down=
  enrichmentNetwork2(ego_cluster@compareClusterResult %>% subset(p.adjust<0.05) %>%
                       split(.$Cluster) %>% lapply(., function(df) df %>% slice_min(pvalue, n=20)) %>% data.table::rbindlist(.) %>%
                       subset(!duplicated(.$ID)),
                     minClusterSize=4,
                     title=NULL,
                     theme=theme_,
                     plotOnly=FALSE,
                     palette="Reds")
# turns out to be exactly the same as in PBMCage dataset

#####################################



### Plot the Ucell scoring heatmap of the rb-related genes in each celltype at the rough level
#... based on the enrichment results above
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(ComplexHeatmap)

### Extract the correlated genes
cor_genes=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$genes
ego_cluster=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
ego_results=ego_cluster@compareClusterResult %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]}) %>%
  Reduce(c, .) %>%
  .[!duplicated(.)]

### Score the obj
PseudobulkdObj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Slimmed_5wCellsForWGCNA.rds")
subset_Ucell=AddModuleScore_UCell(PseudobulkdObj, features=list(ribo=ego_genes))
score_df=subset_Ucell[[]]
saveRDS(score_df, "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.rds")
# process the score
score_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.rds")
colnames(score_df)=gsub("_UCell$","",colnames(score_df))
score_df=score_df %>%
  dplyr::select(-orig.ident) %>%
  mutate(Sex=ifelse(Sex=="Female","F","M"),
         cell_id=rownames(.))
data.table::fwrite(score_df, "~/Project_PBMCage/Results/Immunity_results/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")
unlink("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.rds")

### Make the df
CommonAlteredGenes_score=data.table::fread("~/Project_PBMCage/Results/Immunity_results/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
  dplyr::select(-any_of("agecut")) %>%
  dplyr::rename(age=Age, sex=Sex, donor_id=Donor_id) %>%
  # mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:50), "31~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:70), "51~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:97), "71~81", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# take only the terms about ribosomes
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
  dplyr::select(any_of(c("age","donor_id","Annot.detailed","Annot.inter","Annot.rough","sex","cell_id","agecut",
                  colnames(CommonAlteredGenes_score)[grepl("^r",colnames(CommonAlteredGenes_score))])))

# get the most age-correlated term(s) to plot
cor_terms=CommonAlteredGenes_score %>% dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","Annot.rough","sex","cell_id","agecut")))
term_cor_results=lapply(2:ncol(cor_terms), function(idx) cor.test(cor_terms[["age"]], cor_terms[[idx]], method="spearman", exact=F))
rho=lapply(term_cor_results, function(x) x$estimate) %>% unlist(); names(rho)=colnames(cor_terms)[2:ncol(cor_terms)]
rho_pval=lapply(term_cor_results, function(x) x$p.value) %>% unlist(); names(rho_pval)=colnames(cor_terms)[2:ncol(cor_terms)]

# get the ages that occur in all the groups (orientation*celltype)
age_shared=CommonAlteredGenes_score %>% group_by(Annot.rough, sex, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot
ht_list=list()
SCORE_DFs=CommonAlteredGenes_score %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|^age|Annot\\.|^sex$|cell_id",colnames(.))], names_to="term", values_to="scores") %>%
    group_by(agecut, sex) %>%
    summarize_at("scores", mean) %>%
    tidyr::pivot_wider(names_from="agecut", values_from="scores") %>%
    tibble::column_to_rownames("sex") %>%
    .[,age_shared]
  score_df_up=score_df_up-rowMeans(score_df_up)

  # plot
  ht_list[[i]]=
    Heatmap(score_df_up,
            name=" ",
            heatmap_legend_param=list(
              title="Terekhova et al. (2023) dataset\nUCell Score change",
              legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F,
            row_names_gp=gpar(fontsize=9),
            show_column_names=F,
            row_title=unique(score_df$Annot.rough),
            row_title_gp=gpar(fontsize=10), row_title_rot=0,
            height=unit(3.5,"mm")*nrow(score_df_up),
            width=unit(3.5,"mm")*ncol(score_df_up))
}
ht_list=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared),
               labels=colnames(score_df_up)[1:length(age_shared)],
               legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ha=HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"),
                                     col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                     annotation_name_gp=gpar(fontsize=10, fontfamily="mono"), show_annotation_name=T)

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ha %v% ht_list, merge_legend=FALSE, annotation_legend_list=lgd_age, heatmap_legend_side="right")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_acrossAges_ExprHeatmap_FvsM.pdf", height=3.5, width=3.5)
draw(whole_ht)
dev.off()

#####################################



### Plot the Ucell scoring heatmap of the rb-related genes in each celltype at the rough level
#... based on the enrichment results above
#... with the same agecuts as in the PBMCage dataset
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(ComplexHeatmap)

### Make the df
CommonAlteredGenes_score=data.table::fread("~/Project_PBMCage/Results/Immunity_results/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
  dplyr::select(-any_of("agecut")) %>%
  dplyr::rename(age=Age, sex=Sex, donor_id=Donor_id) %>%
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
  dplyr::rename(agecut=agecut.half)

# take only the terms about ribosomes
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
  dplyr::select(any_of(c("age","donor_id","Annot.detailed","Annot.inter","Annot.rough","sex","cell_id","agecut",
                         colnames(CommonAlteredGenes_score)[grepl("^r",colnames(CommonAlteredGenes_score))])))

# get the most age-correlated term(s) to plot
cor_terms=CommonAlteredGenes_score %>% dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","Annot.rough","sex","cell_id","agecut")))
term_cor_results=lapply(2:ncol(cor_terms), function(idx) cor.test(cor_terms[["age"]], cor_terms[[idx]], method="spearman", exact=F))
rho=lapply(term_cor_results, function(x) x$estimate) %>% unlist(); names(rho)=colnames(cor_terms)[2:ncol(cor_terms)]
rho_pval=lapply(term_cor_results, function(x) x$p.value) %>% unlist(); names(rho_pval)=colnames(cor_terms)[2:ncol(cor_terms)]

# get the ages that occur in all the groups (orientation*celltype)
age_shared=CommonAlteredGenes_score %>% group_by(Annot.rough, sex, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
# plot
SCORE.HEATMAP=ht_list=list()
SCORE_DFs=CommonAlteredGenes_score %>% subset(agecut %in% age_shared) %>% split(.$Annot.rough)
for (i in 1:length(SCORE_DFs)) {
  score_df=SCORE_DFs[[i]]
  score_df_up=score_df %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|^age|Annot\\.|^sex$|cell_id",colnames(.))], names_to="term", values_to="scores") %>%
    group_by(agecut, sex) %>%
    summarize_at("scores", mean) %>%
    tidyr::pivot_wider(names_from="agecut", values_from="scores") %>%
    tibble::column_to_rownames("sex") %>%
    .[,age_shared]
  score_df_up=score_df_up-rowMeans(score_df_up)
  
  # save for turning point analysis
  SCORE.HEATMAP[[i]]=score_df_up["F",] %>% as.data.frame() %>%
    tidyr::pivot_longer(cols=colnames(.), names_to="agecut", values_to="score") %>%
    mutate(sex="F") %>%
    rbind(., score_df_up["M",] %>% as.data.frame() %>%
            tidyr::pivot_longer(cols=colnames(.), names_to="agecut", values_to="score") %>%
            mutate(sex="M")) %>%
    mutate(Annot.rough=names(SCORE_DFs)[i])
  
  # plot
  ht_list[[i]]=
    Heatmap(score_df_up,
            name=" ",
            heatmap_legend_param=list(
              title="UCell Score change",
              legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=F, cluster_rows=F,
            row_names_gp=gpar(fontsize=9),
            show_column_names=F,
            row_title=unique(score_df$Annot.rough),
            row_title_gp=gpar(fontsize=10), row_title_rot=0,
            height=unit(3.5,"mm")*nrow(score_df_up),
            width=unit(3.5,"mm")*ncol(score_df_up))
}
ht_list=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared),
               labels=colnames(score_df_up)[1:length(age_shared)],
               legend_height=unit(7.5, "cm"), grid_width=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ha=HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"),
                                     col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                     annotation_name_gp=gpar(fontsize=10, fontfamily="mono"), show_annotation_name=T)

ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ha %v% ht_list, merge_legend=FALSE, annotation_legend_list=lgd_age, heatmap_legend_side="right")

# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_acrossAges_ExprHeatmap_FvsM.pdf", height=3.5, width=3)
# draw(whole_ht)
# dev.off()

# plot the turning points of the ribosome scores (after scaled by heatmap) across Annot.rough 
SCORE.HEATMAP_arranged_B.T=SCORE.HEATMAP %>% data.table::rbindlist() %>% 
  mutate(Annot.rough_sex=paste0(Annot.rough, "_", sex)) %>%
  mutate(group=ifelse(Annot.rough %in% c("B cells","CD4T cells","CD8T cells"), Annot.rough, "others"))
turningpoint_plot=
  ggplot(SCORE.HEATMAP_arranged_B.T, aes(x=agecut, y=score, alpha=sex, color=Annot.rough)) +
  theme_classic() +
  facet_wrap(~group, ncol=1) +
  geom_point(size=0) +
  geom_line(data=SCORE.HEATMAP_arranged_B.T %>% 
              group_by(group, Annot.rough, sex, Annot.rough_sex, agecut) %>% summarize_at("score", mean), 
            aes(group=Annot.rough_sex)) +
  scale_x_discrete(labels=
                     c(rep("",3),"45~47",rep("",6))) +
  scale_alpha_manual(values=c(0.5,1)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(x=NULL, y="UCell Score Change", subtitle="Terekhova et al. (2023) dataset") +
  theme(axis.title=element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        legend.direction="vertical",
        legend.box="vertical",
        legend.position="right",
        plot.subtitle=element_text(size=11),
        plot.title.position="plot") +
  guides(color=guide_legend(ncol=1,byrow=TRUE))
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes_acrossAges_ExprHeatmap_FvsM_turningpoints.pdf", height=4, width=2.75)
plot(turningpoint_plot)
dev.off()

#####################################






### Enrich age-correlated genes after exclusion of ribosome-related ones
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
ego_results=ego_cluster@compareClusterResult %>%
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]

### Exclude them and do EnrichGO
# downregulated ones
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check
# ... turns out to be a part of metabolic processes shared by some celltypes while another part shared by other cell types...

# downregulated ones
upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check
# ... turns out that cytoplasmic translation/RNA splicing/etc. are common, while leukocyte-mediated immunity is shared by some cell types...

# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")

### Plot the enrichment results
ego_objs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")
# get the sorted data by clusterprofiler plotting
plot_down=dotplot(ego_objs$ego_down, showCategory=10)
plot_up=dotplot(ego_objs$ego_up, showCategory=4)
plot_data=list(plot_down$data, plot_up$data)
# plot with ggplot2
MODULE_ENRICH_Plots=list()
for (i in 1:length(plot_data)) {
  module_=plot_data[[i]] %>% dplyr::select(Cluster, Description, GeneRatio, p.adjust, Count)
  module_descrip_levels=levels(module_$Description)
  module_$Description=as.character(module_$Description)
  module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
  module_$Description=forcats::fct_relevel(module_$Description, formatters::wrap_string(module_descrip_levels, width=50, collapse="\n"))

  # since the downreg clusters do not contain CD4T cells, DCs, Other cells, OtherT, we add them manually
  # ... and the upreg clusters do not contain NK cells, we add them manually
  if (i==1) {
    module_=rbind(module_, data.frame(Cluster=paste0("CD4T cells\n(",length(ego_objs[["ego_up"]]@geneClusters[["CD4T cells"]]),")"), Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
    module_=rbind(module_, data.frame(Cluster=paste0("DCs\n(",length(ego_objs[["ego_up"]]@geneClusters[["DCs"]]),")"), Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
    module_=rbind(module_, data.frame(Cluster=paste0("Other cells\n(",length(ego_objs[["ego_up"]]@geneClusters[["Other cells"]]),")"), Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
    module_=rbind(module_, data.frame(Cluster=paste0("OtherT\n(",length(ego_objs[["ego_up"]]@geneClusters[["OtherT"]]),")"), Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
  } else {
    module_=rbind(module_, data.frame(Cluster=paste0("NK cells\n(",length(ego_objs[["ego_up"]]@geneClusters[["NK cells"]]),")"), Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
  }

  # order clusters to make them consistent
  cluster_levels=levels(module_$Cluster)
  reorder_index=gsub(" .*|\\n.*","",cluster_levels)
  reorder_index=match(reorder_index[order(reorder_index)], reorder_index)
  cluster_levels_reorder=cluster_levels[reorder_index]
  module_$Cluster=forcats::fct_relevel(module_$Cluster, cluster_levels_reorder)

  x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10

  plot_=
    ggplot(module_, aes(x=Cluster, y=Description, color=p.adjust, size=Count)) +
    geom_point() +
    scale_colour_gradient(limits=quantile(module_$p.adjust, na.rm=T)[c(1,5)],
                          low="brown4", high="pink",
                          labels=sprintf(fmt="%0.01e", quantile(module_$p.adjust, na.rm=T)[c(1,5)]),
                          breaks=quantile(module_$p.adjust, na.rm=T)[c(1,5)]) +
    # scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
    #                                round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
    #                                max(module_$Count, na.rm=T)),
    #                       labels=c(min(module_$Count, na.rm=T),
    #                                round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
    #                                max(module_$Count, na.rm=T)),
    #                       range=c(2,6)) +
    guides(size=guide_legend(order=1)) +
    labs(title=c("downreg.","upreg.")[i], y=NULL, x=NULL) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10, lineheight=0.75),
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

plot_both=cowplot::plot_grid(plotlist=MODULE_ENRICH_Plots, align="hv", ncol=1, rel_heights=c(1,1))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded_EnrichPlots_Rough.pdf", width=7, height=12)
plot_both
dev.off()

#####################################



### Correlate cytoplasmic translation UCell scores and cell proportion in Immunity
#####################################
###

library(Seurat)
library(UCell)

### UCell scoring of cytoplasmic translation
# extract the related genes
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
ego_objs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")
downreg_genes_enriched=ego_objs$ego_down@compareClusterResult %>% 
  subset(Description=="cytoplasmic translation" & p.adjust<0.05) %>%
  dplyr::select(geneID) %>%
  as.list() %>%
  lapply(., function(x) strsplit(x, split="/")) %>% Reduce(union,.) %>% Reduce(union,.)
upreg_genes_enriched=ego_objs$ego_up@compareClusterResult %>% 
  subset(Description=="cytoplasmic translation" & p.adjust<0.05) %>%
  dplyr::select(geneID) %>%
  as.list() %>%
  lapply(., function(x) strsplit(x, split="/")) %>% Reduce(union,.) %>% Reduce(union,.)
# score the obj
PseudobulkdObj_immunity=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
genes_tobescored=list(cytoplasmic.translation=c(downreg_genes_enriched, upreg_genes_enriched) %>% .[!duplicated(.)])
subset_Ucell_immunity=AddModuleScore_UCell(PseudobulkdObj_immunity, features=genes_tobescored)
score_df_immunity=subset_Ucell_immunity[[]]

### Load the cell proportion (use the monocyte-modified version to be consistent with the annotation to be published)
meta_info_combined_immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t") %>%
  dplyr::select(donor_id, sex, age, Annot.rough, Annot.inter, percent.inter) %>%
  subset(!duplicated(.))

### Merge the cell proportion and the Ucell scores
df_merged_immunity=score_df_immunity %>% 
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
  dplyr::left_join(., meta_info_combined_immunity, by=c("donor_id","age","Annot.inter","Annot.rough")) %>%
  subset(percent.inter>=0.05) %>%
  # subset(Annot.rough!="Monocytes" & Annot.rough!="DCs" & Annot.rough!="Other cells") %>%
  subset(Annot.rough %in% c("CD4T cells","CD8T cells"))
data.table::fwrite(df_merged_immunity, "~/Project_PBMCage/Results/Immunity_results/Ucell_cytoplasmic.tr_With_Cellproportion_Cor_Tcells.txt.gz", sep="\t")

# Immunity_plot=
#   ggplot(df_merged_immunity, aes(x=cytoplasmic.translation_UCell, y=percent.inter, color=Annot.rough)) +
#   geom_point(size=0.25, shape=20, alpha=0.25) +
#   theme_classic() +
#   labs(y="Cell percent (%)", x="UCell Score of Cytoplasmic tr.", title="Terekhova et al. (2023) dataset") +
#   theme(axis.text.x=element_text(size=8),
#         axis.text.y=element_text(size=8),
#         axis.title.x=element_text(size=9),
#         axis.title.y=element_text(size=9),
#         legend.position="right",
#         # legend.title=element_blank(),
#         title=element_text(size=10),
#         plot.title.position="plot",
#         plot.margin=unit(c(t=1, b=1, l=1, r=4), units="mm")) +
#   scale_color_manual(values=scales::hue_pal()(8)[c(2,3)]) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#   geom_smooth(aes(group=Annot.rough), method="lm", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.25) +
#   ggpubr::stat_cor(aes(group=Annot.rough), show.legend=FALSE, size=3.5, method="pearson")
# 
# pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded_Enrich_cytoplasmictr.vs.cellfreq.pdf", width=3.2, height=2)
# Immunity_plot
# dev.off()

#####################################



### Correlate cytoplasmic translation UCell scores and cell proportion in Immunity at Annot.inter
#####################################
###

library(Seurat)
library(UCell)

### UCell scoring of cytoplasmic translation
# extract the related genes
ego_objs=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")
downreg_genes_enriched=ego_objs$ego_down@compareClusterResult %>% 
  subset(Description=="cytoplasmic translation" & p.adjust<0.05) %>%
  dplyr::select(geneID) %>%
  as.list() %>%
  lapply(., function(x) strsplit(x, split="/")) %>% Reduce(union,.) %>% Reduce(union,.)
upreg_genes_enriched=ego_objs$ego_up@compareClusterResult %>% 
  subset(Description=="cytoplasmic translation" & p.adjust<0.05) %>%
  dplyr::select(geneID) %>%
  as.list() %>%
  lapply(., function(x) strsplit(x, split="/")) %>% Reduce(union,.) %>% Reduce(union,.)
# score the obj
PseudobulkdObj_immunity=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
genes_tobescored=list(cytoplasmic.translation=c(downreg_genes_enriched, upreg_genes_enriched) %>% .[!duplicated(.)])
subset_Ucell_immunity=AddModuleScore_UCell(PseudobulkdObj_immunity, features=genes_tobescored)
score_df_immunity=subset_Ucell_immunity[[]]

### Load the cell proportion (use the monocyte-modified version to be consistent with the annotation to be published)
meta_info_combined_immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t") %>%
  dplyr::select(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed, percent.detailed) %>%
  subset(!duplicated(.))

### Merge the cell proportion and the Ucell scores
df_merged_immunity=score_df_immunity %>% 
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
  dplyr::left_join(., meta_info_combined_immunity, by=c("donor_id","age","Annot.detailed","Annot.inter","Annot.rough")) %>%
  # subset(percent.inter>=0.05) %>%
  subset(Annot.rough %in% c("CD4T cells","CD8T cells"))
data.table::fwrite(df_merged_immunity, "~/Project_PBMCage/Results/Immunity_results/Ucell_cytoplasmic.tr_With_Cellproportion_Cor_Tcells_Annot.inter.txt.gz", sep="\t")

#####################################



### Get all the genes annotated to the key terms enriched based on the GO gaf file
#####################################
###

library(dplyr)

### Setup superclusters
Super_cluster=list(`energy metab.|'nucleoside triphosphate metabolic process_cellular respiration_electron transport chain_proton transmembrane transport_reactive nitrogen species metabolic process'`=
                     c("GO:0009141","GO:0045333","GO:0022900","GO:1902600","GO:2001057"),
                   `nucl. metab.|'nucleic acid metabolic process_nucleoside phosphate metabolic process'`=c("GO:0090304","GO:0006753"),
                   `prot. metab.|'cytoplasmic translation_protein maturation_vesicle-mediated transport_intracellular transport_protein modification process_regulation of protein modification process_proteolysis_protein catabolic process_regulation of protein stability'`=
                     c("GO:0002181","GO:0051604","GO:0016192","GO:0046907","GO:0036211","GO:0031399","GO:0006508","GO:0030163","GO:0031647"),
                   `ribosome syn.|'ribosome biogenesis'`=c("GO:0042254"),
                   `cellular org.|'organelle localization_cellular component organization_cellular macromolecule localization_maintenance of location'`=c("GO:0051640","GO:0016043","GO:0070727","GO:0051235"),
                   `cell division|'cell division_cell cycle'`=c("GO:0051301","GO:0007049"),
                   `autophagy|'autophagy'`=c("GO:0006914"),
                   `prog. death|'programmed cell death'`=c("GO:0012501"),
                   `anti-virus|'viral process_biological process involved in symbiotic interaction'`=c("GO:0016032","GO:0044403"),
                   `cell response|'response to stimulus'`=c("GO:0050896"),
                   `leuk. devel.|'hemopoiesis_regulation of cell differentiation'`=c("GO:0030097","GO:0045595"),
                   `immune proc.|'cell-cell adhesion_negative regulation of cell activation_platelet activation_follicular dendritic cell activation_leukocyte activation_cell activation involved in immune response_thrombocyte activation_positive regulation of cell activation_regulation of cell activation_leukocyte proliferation_MHC protein complex assembly_cytokine production_cell killing'`=
                     c("GO:0098609","GO:0050866","GO:0030168","GO:0045321","GO:0002263","GO:0071892","GO:0050867","GO:0050865","GO:0070661","GO:0002376","GO:0002396","GO:0001816","GO:0001906"))

### Extract all the GO terms related to each supercluster
GOSim::setOntology(ont="BP", loadIC=FALSE, DIR=NULL)
other_allGOs=list()
for (i in 1:length(Super_cluster)) {
  other_offsprings=GOSim::getOffsprings()[Super_cluster[[i]]] %>% Reduce(c,.)
  tb=toTable(GO.db::GOBPCHILDREN); colnames(tb)[c(1,2)]=c("child","parent")
  other_offsprings=c(other_offsprings, subset(tb, parent %in% Super_cluster[i]) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]

  j=1
  for (j in 1:10) {
    other_offsprings=c(other_offsprings, subset(tb, parent %in% other_offsprings) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]
    j=j+1
  }

  other_allGOs=c(other_allGOs, list(c(unname(Super_cluster[[i]]), other_offsprings) %>% .[!duplicated(.)]))
}
# combine the mt and the others
names(other_allGOs)=names(Super_cluster)
saveRDS(other_allGOs, "~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster.rds")

### Collect the genes from each supercluster
# filter with the obj features to make sure that the genes exist
seuratobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
features_=rownames(seuratobj)
gaf_file=mgsa::readGAF("~/enrichment/goa_human.gaf.gz", evidence=NULL)
Genes_in_SuperClusters=list()
for (i in 1:length(other_allGOs)) {
  go_list=other_allGOs[[i]]
  geneIndices=Reduce(c, gaf_file@sets[go_list]) %>% .[!duplicated(.)] %>% .[!is.na(.)]
  genes_=gaf_file@itemAnnotations[geneIndices,]$symbol %>% .[. %in% features_]
  Genes_in_SuperClusters[[i]]=genes_
}
names(Genes_in_SuperClusters)=names(Super_cluster)
# save
saveRDS(Genes_in_SuperClusters, "~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")

### Map age-corr. genes in each celltypes to the SuperCluster genes
# get the downreg./upreg. age-correlated genes in each celltype
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
Cor_analysis_DF_subset=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  mutate(orientation=ifelse(rho<0,"downreg.","upreg.")) %>%
  arrange(desc(abs(rho))) %>%
  split(.$orientation) %>%
  .[c("downreg.","upreg.")] %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% .[celltypes_] %>%
           lapply(., function(dff) dff %>% dplyr::select(gene, rho, rho_pval, celltypes, orientation)))
# map to the genes in superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=list()
for (o in 1:length(Cor_analysis_DF_subset)) {
  orientation_subset=Cor_analysis_DF_subset[[o]]
  exist_genes_per_celltype=list()
  for (i in 1:length(orientation_subset)) {
    celltype_subset=orientation_subset[[i]]
    exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters), function(idx) celltype_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]))
    names(exist_genes_per_supercluster)=names(Genes_in_SuperClusters)
    exist_genes_per_celltype[[i]]=exist_genes_per_supercluster
  }
  names(exist_genes_per_celltype)=names(orientation_subset)
  exist_genes_per_orientation[[o]]=exist_genes_per_celltype
}
names(exist_genes_per_orientation)=c("downreg.","upreg.")

saveRDS(exist_genes_per_orientation, "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters.rds")

#####################################



### Jaccard similarity in shared terms apart from rb-related across rough celltypes
#####################################
###

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Load the data
exist_genes_per_orientation=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters.rds")
# arrange the genesets
down_genes=exist_genes_per_orientation$downreg.
down_genes_supercluster=lapply(1:length(down_genes[[1]]), function(idx) lapply(down_genes, function(celltype_subset) celltype_subset[[idx]]))
names(down_genes_supercluster)=names(down_genes[[1]])

up_genes=exist_genes_per_orientation$upreg.
up_genes_supercluster=lapply(1:length(up_genes[[1]]), function(idx) lapply(up_genes, function(celltype_subset) celltype_subset[[idx]]))
names(up_genes_supercluster)=names(up_genes[[1]])

# Plot downreg.
ht_down=list()
for (supercluster in 1:length(down_genes_supercluster)) {
  genelists=down_genes_supercluster[[supercluster]]
  # plot only lymphocytes
  genelists=genelists[c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")]

  # calculate the jaccard similarity
  Jaccard_SIM=Jaccard_SIM_Name=c()
  for (i in 1:length(genelists)) {
    the_first_celltype=genelists[[i]] %>% dplyr::select(gene) %>% tibble::deframe()

    for (j in 1:length(genelists)) {
      the_second_celltype=genelists[[j]] %>% dplyr::select(gene) %>% tibble::deframe()

      jaccard_similarity=ifelse(length(the_first_celltype)==0|length(the_second_celltype)==0, NA,
                                jaccard(the_first_celltype, the_second_celltype))
      Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
      Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(genelists)[i],"_vs_",names(genelists)[j]))
    }
  }
  names(Jaccard_SIM)=Jaccard_SIM_Name
  Jaccard_SIM_pos=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)

  # merge the pos and neg results
  Jaccard_SIM_pos=Jaccard_SIM_pos %>%
    mutate(celltype_1=gsub("_vs_.*","",pair), celltype_2=gsub(".*_vs_","",pair))

  # plot the results
  Jaccard_SIM_df=Jaccard_SIM_pos %>%
    dplyr::select(-pair) %>%
    tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
    tibble::column_to_rownames("celltype_1")

  range(Jaccard_SIM_df) # check
  col_fun_down=circlize::colorRamp2(c(0, 0.5, 1), c("white","dodgerblue3","darkblue"))
  ht_down[[supercluster]]=
    Heatmap(Jaccard_SIM_df,
            name=" ",
            show_heatmap_legend=FALSE,
            border_gp=gpar(col="black", lty=2),
            col=col_fun_down,
            rect_gp=gpar(col="white", lwd=0.5),
            cluster_columns=F, cluster_rows=F,
            row_names_gp=gpar(fontsize=10), show_row_names=F,
            column_title=gsub("\\|.*","",names(down_genes_supercluster))[supercluster],
            column_title_gp=gpar(fontsize=11),
            column_names_side="bottom",
            column_names_rot=90,
            column_names_gp=gpar(fontsize=10),
            height=unit(5, "mm")*nrow(Jaccard_SIM_df),
            width=unit(5, "mm")*ncol(Jaccard_SIM_df))
}
ht_list_down=Reduce(`+`,ht_down)
draw(ht_list_down)

# Plot upreg.
ht_up=list()
for (supercluster in 1:length(up_genes_supercluster)) {
  genelists=up_genes_supercluster[[supercluster]]
  # plot only lymphocytes
  genelists=genelists[c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")]

  # calculate the jaccard similarity
  Jaccard_SIM=Jaccard_SIM_Name=c()
  for (i in 1:length(genelists)) {
    the_first_celltype=genelists[[i]] %>% dplyr::select(gene) %>% tibble::deframe()

    for (j in 1:length(genelists)) {
      the_second_celltype=genelists[[j]] %>% dplyr::select(gene) %>% tibble::deframe()

      jaccard_similarity=ifelse(length(the_first_celltype)==0|length(the_second_celltype)==0, NA,
                                jaccard(the_first_celltype, the_second_celltype))
      Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
      Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(genelists)[i],"_vs_",names(genelists)[j]))
    }
  }
  names(Jaccard_SIM)=Jaccard_SIM_Name
  Jaccard_SIM_pos=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)

  # merge the pos and neg results
  Jaccard_SIM_pos=Jaccard_SIM_pos %>%
    mutate(celltype_1=gsub("_vs_.*","",pair), celltype_2=gsub(".*_vs_","",pair))

  # plot the results
  Jaccard_SIM_df=Jaccard_SIM_pos %>%
    dplyr::select(-pair) %>%
    tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
    tibble::column_to_rownames("celltype_1")

  range(Jaccard_SIM_df) # check
  col_fun_up=circlize::colorRamp2(c(0, 0.5, 1), c("white", "brown3", "darkred"))
  ht_up[[supercluster]]=
    Heatmap(Jaccard_SIM_df,
            name=" ",
            show_heatmap_legend=FALSE,
            border_gp=gpar(col="black", lty=2),
            col=col_fun_up,
            rect_gp=gpar(col="white", lwd=0.5),
            cluster_columns=F, cluster_rows=F,
            row_names_gp=gpar(fontsize=10), show_row_names=F,
            column_title=gsub("\\|.*","",names(down_genes_supercluster))[supercluster],
            column_title_gp=gpar(fontsize=11),
            column_names_side="bottom",
            column_names_rot=90,
            column_names_gp=gpar(fontsize=10),
            height=unit(5, "mm")*nrow(Jaccard_SIM_df),
            width=unit(5, "mm")*ncol(Jaccard_SIM_df))
}
ht_list_up=Reduce(`+`,ht_up)
draw(ht_list_up)

### Add legend
lgd_up=Legend(col_fun=col_fun_up,
              title="Jaccard similarity\nupreg.",
              legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
)

lgd_down=Legend(col_fun=col_fun_down,
                title="Jaccard similarity\ndownreg.",
                legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
                direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
)

ht_final_up=draw(ht_list_up, annotation_legend_list=lgd_up, annotation_legend_side="left")
pdf("~/Project_PBMCage/Plots/Immunity/Bulkanalysis_RNAexpr.Cor.wAges_RbExcludedUpregGenes_JaccardSimilarity_Rough.pdf", height=2.5, width=10.5)
draw(ht_final_up)
dev.off()

ht_final_down=draw(ht_list_down, annotation_legend_list=lgd_down, annotation_legend_side="left")
pdf("~/Project_PBMCage/Plots/Immunity/Bulkanalysis_RNAexpr.Cor.wAges_RbExcludedDownGenes_RbExcluded_JaccardSimilarity_Rough.pdf", height=2.5, width=10.5)
draw(ht_final_down)
dev.off()

#####################################



### Pie plot to show the proportions of the enriched terms in all the superclusters per celltype
#####################################
###

library(dplyr)

### Load the GO terms in each supercluster
All_relatedGOterms_total=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster.rds")

### Load the enrichments including rb
ego_down=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
ego_up=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Rough.rds")$ego
enriched_terms=list(downreg.=ego_down@compareClusterResult %>% subset(p.adjust<1e-5) %>% split(.$Cluster) %>% lapply(., function(df) df %>% dplyr::select(ID) %>% tibble::deframe(.)),
                    upreg.=ego_up@compareClusterResult %>% subset(p.adjust<1e-5) %>% split(.$Cluster) %>% lapply(., function(df) df %>% dplyr::select(ID) %>% tibble::deframe(.)))

### Map the enriched terms to the superclusters
GOs_in_supercluster_list_bothorientation=list()
for (i in 1:length(enriched_terms)) {
  per_orientation=enriched_terms[[i]]
  GOs_in_supercluster_list=list()
  for (j in 1:length(per_orientation)) {
    per_celltype=per_orientation[[j]]
    GOs_in_supercluster=per_celltype[per_celltype %in% Reduce(c,All_relatedGOterms_total)]
    GOs_in_supercluster_list[[j]]=GOs_in_supercluster
  }
  names(GOs_in_supercluster_list)=paste0(names(per_orientation),",ct,")
  GOs_in_supercluster_list_bothorientation[[i]]=GOs_in_supercluster_list
}
names(GOs_in_supercluster_list_bothorientation)=paste0(names(enriched_terms),",o,")

### Make the df
length_terms_inSuperCl=sapply(All_relatedGOterms_total, length)
length_terms_exist=sapply(unlist(GOs_in_supercluster_list_bothorientation, recursive=FALSE), length)
length_terms_exist_df=data.frame(group=names(length_terms_exist), termN=length_terms_exist) %>%
  mutate(orientation=gsub(",o,.*","",group)) %>%
  mutate(celltypes=gsub(".*,o,\\.|,ct,","",group)) %>%
  dplyr::select(-group)

length_terms_enriched=sapply(unlist(enriched_terms, recursive=FALSE) %>% .[!duplicated(.)], length)
length_terms_enriched_df=data.frame(group=names(length_terms_enriched), termN=length_terms_enriched) %>%
  mutate(orientation=paste0(gsub("\\..*","",group),".")) %>%
  mutate(celltypes=gsub(".*\\.\\.","",group)) %>%
  dplyr::select(-group)

combine_df=length_terms_exist_df %>% left_join(., length_terms_enriched_df, by=c("orientation","celltypes"))
colnames(combine_df)[c(1,4)]=c("N_in_SuperClusters","N_of_all_enrichments")
combine_df=combine_df %>% group_by(orientation) %>%
  summarize_at(c("N_in_SuperClusters","N_of_all_enrichments"), function(x) sum(x, na.rm=T)) %>%
  dplyr::rename(`Key terms`=`N_in_SuperClusters`,
                All=`N_of_all_enrichments`) %>%
  mutate(Others=All-`Key terms`) %>%
  tidyr::pivot_longer(cols=c("Key terms","Others"), names_to="in_or_not", values_to="termN") %>%
  mutate(percent=termN/All)

plot_=
  ggplot(combine_df, aes(x="", y=termN, fill=in_or_not)) +
    geom_bar(stat="identity", width=0.5, color="white") +
    facet_wrap(~orientation) +
    geom_text(data=. %>% subset(in_or_not=="Others"), aes(y=termN/2, label=paste0(sprintf("%.2f",percent*100), "%")), color="black") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=c("grey40","grey80"), labels=c("Key terms","Others")) +
    theme_light() +
    theme(axis.text.x=element_text(size=9),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=9, angle=0),
          axis.title.y=element_blank(),
          strip.text=element_text(size=10, color="black"),
          strip.background=element_rect(fill="white", color="black"),
          legend.position="top",
          legend.direction="horizontal"
    ) +
    guides(fill=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded_Keyterms.Pieplots_Regardlesscelltypes.pdf", height=3, width=4)
plot_
dev.off()

#####################################



### Pie plot to show the proportions of the age-corr. genes in each supercluster per celltype at rough level
#####################################
###

library(clusterProfiler)
library(dplyr)

### Load the age-corr. genes and their enrichment in superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters.rds")
exist_genes_per_orientation=lapply(exist_genes_per_orientation, function(x) lapply(x, function(y) lapply(y, function(z) z %>% dplyr::select(gene) %>% tibble::deframe())))
# extract the genes not in any of the supercluster
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
not_there_geneNO=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  mutate(orientation=ifelse(rho<0,"downreg.","upreg.")) %>%
  split(.$orientation) %>%
  .[c("downreg.","upreg.")] %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% .[celltypes_] %>%
           lapply(., function(dff) dff %>% dplyr::select(gene) %>% tibble::deframe(.))) %>%
  lapply(., function(gene_per_orientation) lapply(gene_per_orientation, function(x) length(x[!(x %in% Reduce(c, Genes_in_SuperClusters))])))
# calculate the number of age-corr. genes in each supercluster per celltype
there_geneNO=exist_genes_per_orientation %>%
  lapply(., function(perorientation) lapply(perorientation, function(percelltype) lapply(percelltype, length)))

### Make the df
name_group_orientation=rep(c("downreg.","upreg."), each=length(celltypes_)*length(Genes_in_SuperClusters))
name_group_celltypes=rep(rep(celltypes_, each=length(Genes_in_SuperClusters)), times=2)
name_group_supercluster=rep(names(Genes_in_SuperClusters), times=2*length(celltypes_))
tempt=as.data.frame(there_geneNO)
df=data.frame(group=colnames(tempt), geneN=as.integer(tempt))
df$orientation=name_group_orientation
df$celltypes=name_group_celltypes
df$supercluster=name_group_supercluster
df=df %>% dplyr::select(-group)

name_group_orientation2=rep(c("downreg.","upreg."), each=length(celltypes_))
name_group_celltypes2=rep(celltypes_, times=2)
tempt2=as.data.frame(not_there_geneNO)
df2=data.frame(group=colnames(tempt2), geneN=as.integer(tempt2))
df2$orientation=name_group_orientation2
df2$celltypes=name_group_celltypes2
df2=df2 %>% dplyr::select(-group) %>%
  mutate(supercluster="other")

df_final=rbind(df, df2) %>%
  mutate(in_or_not=ifelse(supercluster!="other","in","not"))

plot_=
  ggplot(df_final, aes(x="", y=geneN, fill=in_or_not)) +
  facet_grid(orientation~celltypes) +
  geom_bar(stat="identity", width=0.5, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=c("grey20","grey80"), labels=c("Key terms","Others")) +
  theme_light() +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=9, angle=0),
        axis.title.y=element_blank(),
        strip.text=element_text(size=10, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="top",
        legend.direction="horizontal"
  ) +
  scale_y_continuous(breaks=c(0,2000,4000,6000,8000), labels={function(x) paste0(x/1000, "k")}) +
  guides(fill=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded_Keyterms.Pieplots_Rough.pdf", height=3.5, width=10)
plot_
dev.off()

### Plot the detailed proportions of the key enriched terms of age-corr. genes in lymphocytes at rough level
# df_subset=subset(df, celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT"))
supecluster_=names(table(df$supercluster))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
supecluster_reorder=supecluster_[supecluster_reorder_idx]
df$supercluster=forcats::fct_relevel(df$supercluster, supecluster_reorder)
plot_2=
  ggplot(df, aes(x=celltypes, fill=supercluster)) +
  geom_bar(data=.%>% subset(orientation=="upreg."), stat="identity",
           aes(y=geneN), position=position_dodge(width=0.7), width=0.3) +
  geom_bar(data=.%>% subset(orientation=="downreg."), stat="identity",
           aes(y=-geneN), position=position_dodge(width=0.7), width=0.3) +
  geom_hline(yintercept=0, linewidth=0.5, linetype="dotted", color="black") +
  scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Red_Blue_Brown"), labels=function(x){gsub("\\|.*","",x)}) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=9, angle=0),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="top",
        legend.direction="horizontal") +
  scale_y_continuous(labels={function(x) abs(x)}) +
  guides(fill=guide_legend(title=NULL, ncol=3))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded_Keyterms.Barplots_Rough.pdf", height=3.5, width=7)
plot_2
dev.off()


### Plot the proportions of age-correlated genes in each supercluster
# get the age-corr. genes in each supercluster
exist_genes_per_orientation=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters.rds")
exist_genes_per_orientation=lapply(exist_genes_per_orientation, function(x) lapply(x, function(y) lapply(y, function(z) z %>% dplyr::select(gene) %>% tibble::deframe())))
# get all the age-corr. genes in each celltype at different directions
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
all_cor_genes=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  mutate(orientation=ifelse(rho<0,"downreg.","upreg.")) %>%
  split(.$orientation) %>%
  .[c("downreg.","upreg.")] %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% .[celltypes_] %>%
           lapply(., function(dff) dff %>% dplyr::select(gene) %>% tibble::deframe(.)))
# calculate the proportions of genes in each supercluster
gene_proportion=lapply(names(exist_genes_per_orientation),
                       function(orientation) lapply(names(exist_genes_per_orientation[[orientation]]),
                                                    function(celltype) lapply(1:length(exist_genes_per_orientation[[orientation]][[celltype]]),
                                                                              function(idx) (length(exist_genes_per_orientation[[orientation]][[celltype]][[idx]] %>%
                                                                                                      .[. %in% all_cor_genes[[orientation]][[celltype]]]))/length(all_cor_genes[[orientation]][[celltype]]))))
names(gene_proportion)=c("downreg.","upreg.")
gene_proportion=lapply(gene_proportion, function(x) {names(x)=names(exist_genes_per_orientation[[1]]); x})
gene_proportion=lapply(gene_proportion, function(x) lapply(x, function(y) {names(y)=names(exist_genes_per_orientation[[1]][[1]]); y}))
# clean the df
name_group_orientation3=rep(c("downreg.","upreg."), each=length(celltypes_)*length(Genes_in_SuperClusters))
name_group_celltypes3=rep(rep(celltypes_, each=length(Genes_in_SuperClusters)), times=2)
name_group_supercluster3=rep(names(Genes_in_SuperClusters), times=2*length(celltypes_))
tempt=as.data.frame(gene_proportion)
df3=data.frame(group=colnames(tempt), geneportion=as.numeric(tempt))
df3$orientation=name_group_orientation3
df3$celltypes=name_group_celltypes3
df3$supercluster=name_group_supercluster3
df3=df3 %>% dplyr::select(-group)
supecluster_=names(table(df3$supercluster))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
supecluster_reorder=supecluster_[supecluster_reorder_idx]
df3$supercluster=forcats::fct_relevel(df3$supercluster, supecluster_reorder)

ggplot(df3, aes(x=supercluster, y=geneportion, fill=celltypes)) +
  geom_bar(stat="identity", position=position_dodge(width=0.5), width=0.5) +
  scale_x_discrete(guide=guide_axis(n.dodge=2), labels=function(x){gsub("\\|.*","",x)})

#####################################



### Get the genes related to positive/negative regulation of the key terms enriched based on the GO gaf file
#####################################
###

library(dplyr)

### Setup superclusters
Super_cluster=list(`energy metab.|'nucleoside triphosphate metabolic process_cellular respiration_electron transport chain_proton transmembrane transport_reactive nitrogen species metabolic process'`=
                     c("GO:0009141","GO:0045333","GO:0022900","GO:1902600","GO:2001057"),
                   `nucl. metab.|'nucleic acid metabolic process_nucleoside phosphate metabolic process'`=c("GO:0090304","GO:0006753"),
                   `prot. metab.|'cytoplasmic translation_protein maturation_vesicle-mediated transport_intracellular transport_protein modification process_regulation of protein modification process_proteolysis_protein catabolic process_regulation of protein stability'`=
                     c("GO:0002181","GO:0051604","GO:0016192","GO:0046907","GO:0036211","GO:0031399","GO:0006508","GO:0030163","GO:0031647"),
                   `ribosome syn.|'ribosome biogenesis'`=c("GO:0042254"),
                   `cellular org.|'organelle localization_cellular component organization_cellular macromolecule localization_maintenance of location'`=c("GO:0051640","GO:0016043","GO:0070727","GO:0051235"),
                   `cell division|'cell division_cell cycle'`=c("GO:0051301","GO:0007049"),
                   `autophagy|'autophagy'`=c("GO:0006914"),
                   `prog. death|'programmed cell death'`=c("GO:0012501"),
                   `anti-virus|'viral process_biological process involved in symbiotic interaction'`=c("GO:0016032","GO:0044403"),
                   `cell response|'response to stimulus'`=c("GO:0050896"),
                   `leuk. devel.|'hemopoiesis_regulation of cell differentiation'`=c("GO:0030097","GO:0045595"),
                   `immune proc.|'cell-cell adhesion_negative regulation of cell activation_platelet activation_follicular dendritic cell activation_leukocyte activation_cell activation involved in immune response_thrombocyte activation_positive regulation of cell activation_regulation of cell activation_leukocyte proliferation_MHC protein complex assembly_cytokine production_cell killing'`=
                     c("GO:0098609","GO:0050866","GO:0030168","GO:0045321","GO:0002263","GO:0071892","GO:0050867","GO:0050865","GO:0070661","GO:0002376","GO:0002396","GO:0001816","GO:0001906"))

### Extract all the GO terms related to regulation of each supercluster
GOSim::setOntology(ont="BP", loadIC=FALSE, DIR=NULL)
allGOs_pos=allGOs_neg=list()
for (i in 1:length(Super_cluster)) {
  other_offsprings=GOSim::getOffsprings()[Super_cluster[[i]]] %>% Reduce(c,.)
  tb=BiocGenerics::toTable(GO.db::GOBPCHILDREN); colnames(tb)[c(1,2)]=c("child","parent")
  other_offsprings=c(other_offsprings, subset(tb, parent %in% Super_cluster[i]) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]

  j=1
  for (j in 1:10) {
    other_offsprings=c(other_offsprings, subset(tb, parent %in% other_offsprings) %>% dplyr::select(child) %>% tibble::deframe(.)) %>% .[!duplicated(.)]
    j=j+1
  }

  other_offsprings_pos=tb %>% subset(child %in% other_offsprings) %>%
    subset(RelationshipType=="positively regulates") %>%
    select(child) %>% tibble::deframe()
  allGOs_pos=c(allGOs_pos, list(other_offsprings_pos))

  other_offsprings_neg=tb %>% subset(child %in% other_offsprings) %>%
    subset(RelationshipType=="negatively regulates") %>%
    select(child) %>% tibble::deframe()
  allGOs_neg=c(allGOs_neg, list(other_offsprings_neg))
}
names(allGOs_pos)=names(Super_cluster)
names(allGOs_neg)=names(Super_cluster)
saveRDS(list(allGOs_pos, allGOs_neg), "~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")

### Collect the genes from each supercluster
# filter with the obj features to make sure that the genes exist
seuratobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
features_=rownames(seuratobj)
gaf_file=mgsa::readGAF("~/enrichment/goa_human.gaf.gz", evidence=NULL)
Genes_in_SuperClusters_pos=Genes_in_SuperClusters_neg=list()
for (i in 1:length(allGOs_pos)) {
  go_list=allGOs_pos[[i]]
  geneIndices=Reduce(c, gaf_file@sets[go_list]) %>% .[!duplicated(.)] %>% .[!is.na(.)]
  genes_=gaf_file@itemAnnotations[geneIndices,]$symbol %>% .[. %in% features_]
  Genes_in_SuperClusters_pos[[i]]=genes_
}
names(Genes_in_SuperClusters_pos)=names(Super_cluster)

for (i in 1:length(allGOs_neg)) {
  go_list=allGOs_neg[[i]]
  geneIndices=Reduce(c, gaf_file@sets[go_list]) %>% .[!duplicated(.)] %>% .[!is.na(.)]
  genes_=gaf_file@itemAnnotations[geneIndices,]$symbol %>% .[. %in% features_]
  Genes_in_SuperClusters_neg[[i]]=genes_
}
names(Genes_in_SuperClusters_neg)=names(Super_cluster)
# save
saveRDS(list(Genes_in_SuperClusters_pos, Genes_in_SuperClusters_neg),
        "~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster_posAndnegRegulate.rds")

### Map age-corr. genes in each celltypes to the SuperCluster genes
# get the downreg./upreg. age-correlated genes in each celltype
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
Cor_analysis_DF_subset=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  arrange(desc(abs(rho))) %>%
  split(.$celltypes) %>%
  .[celltypes_] %>%
  lapply(., function(dff) dff %>% dplyr::select(gene, rho, rho_pval, celltypes))
# map to the pos-reg genes in superclusters
Genes_in_SuperClusters_pos=
  readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster_posAndnegRegulate.rds")[[1]]
exist_genes_per_celltype_pos=list()
for (i in 1:length(Cor_analysis_DF_subset)) {
  celltype_subset=Cor_analysis_DF_subset[[i]]
  exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters_pos),
                                      function(idx) celltype_subset %>% subset(gene %in% Genes_in_SuperClusters_pos[[idx]]) %>%
                                        mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters_pos)[idx])) %>%
                                        mutate(term_orientation="positively")) %>%
    data.table::rbindlist()
  exist_genes_per_celltype_pos[[i]]=exist_genes_per_supercluster
}
exist_genes_per_celltype_pos=data.table::rbindlist(exist_genes_per_celltype_pos)
# map to the neg-reg genes in superclusters
Genes_in_SuperClusters_neg=
  readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster_posAndnegRegulate.rds")[[2]]
exist_genes_per_celltype_neg=list()
for (i in 1:length(Cor_analysis_DF_subset)) {
  celltype_subset=Cor_analysis_DF_subset[[i]]
  exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters_neg),
                                      function(idx) celltype_subset %>% subset(gene %in% Genes_in_SuperClusters_neg[[idx]]) %>%
                                        mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters_neg)[idx])) %>%
                                        mutate(term_orientation="negatively")) %>%
    data.table::rbindlist()
  exist_genes_per_celltype_neg[[i]]=exist_genes_per_supercluster
}
exist_genes_per_celltype_neg=data.table::rbindlist(exist_genes_per_celltype_neg)

exist_genes_per_celltype_merged=data.table::rbindlist(list(exist_genes_per_celltype_pos, exist_genes_per_celltype_neg))

saveRDS(exist_genes_per_celltype_merged, "~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters_posAndnegRegulate.rds")

#####################################



### Enrichment term scoring of the 12 superclusters in each celltype at the rough level
#####################################
###

### Load the GO terms for each supercluster
GOs_both=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
GOs_both_df=
  as.data.frame(unlist(GOs_both)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_both)') %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_pos=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[1]]
names(GOs_pos)=gsub("\\|.*","",names(GOs_pos))
GOs_pos_df=
  as.data.frame(unlist(GOs_pos)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_pos)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_neg=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[2]]
names(GOs_neg)=gsub("\\|.*","",names(GOs_neg))
GOs_neg_df=
  as.data.frame(unlist(GOs_neg)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_neg)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))
  
GOs_all=data.table::rbindlist(list(GOs_pos_df, GOs_neg_df, GOs_both_df))

### Load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Rough.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)

### Map the GO IDs to the 12 supercluster
merged_df=full_join(ego_all, GOs_all, by=c("GO_ID")) %>% 
  subset(!is.na(Cluster) & !is.na(supercluster)) %>%
  mutate(p.adjust_wSign=ifelse(orientation=="upreg",p.adjust,-p.adjust)) %>%
  group_by(Cluster, supercluster, orientation) %>%
  slice_max(p.adjust, n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(Cluster, supercluster) %>%
  summarize_at("p.adjust_wSign", sum) %>%
  left_join(data.frame(Cluster='B cells', supercluster=NA, p.adjust_wSign=NA)) %>%
  tidyr::pivot_wider(names_from="supercluster", id_cols="Cluster", values_from="p.adjust_wSign") %>%
  tibble::column_to_rownames("Cluster") %>%
  replace(is.na(.), 0)

mtx_=scale(merged_df, scale=T, center=F)
# add back the missing celltypes
mtx_=rbind(mtx_, matrix(0, ncol=ncol(mtx_), nrow=8-nrow(mtx_)))
rownames(mtx_)=c(rownames(mtx_)[1:6], "B cells", "Other cells")
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]

range(mtx_, na.rm=T) # check
col_fun_sub=circlize::colorRamp2(c(-3, 0, 3), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(mtx_,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Scaled enrichment index",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Terekhova et al. (2023) dataset", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(6, "mm"),
          height=nrow(mtx_)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed.pdf", height=3, width=4.5)
draw(ht_final)
dev.off()

#####################################



### Ucell scoring of the pos/neg. regulated key-terms in each celltype at the rough level
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(ComplexHeatmap)

### Extract the age-corr. genes related to the superclusters
exist_genes_in_terms=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters_posAndnegRegulate.rds")
# get the gene list
orientation_=unique(exist_genes_in_terms$term_orientation)
superclusters_=unique(exist_genes_in_terms$superclusters)
genes_perorientation_persupercl=list()
for (i in 1:length(superclusters_)) {
  # since there're too many signatures for Ucell, we take top30
  exist_genes_in_terms_subset=exist_genes_in_terms %>% subset(superclusters==superclusters_[i]) %>%
    split(.$term_orientation) %>%
    lapply(., function(df) df %>% 
             mutate(rank=abs(rho)*(-log10(rho_pval))) %>%
             slice_max(rank, n=30) %>% dplyr::select(gene) %>% tibble::deframe())
  names(exist_genes_in_terms_subset)=paste0(superclusters_[i], "__",names(exist_genes_in_terms_subset))
  genes_perorientation_persupercl=c(genes_perorientation_persupercl, exist_genes_in_terms_subset)
}

### Score the sc obj
PseudobulkdObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
subset_Ucell=AddModuleScore_UCell(PseudobulkdObj, features=genes_perorientation_persupercl)
SCORE_DF=subset_Ucell[[]]
saveRDS(SCORE_DF, "~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_OneSamplePerDonorSingleCell.rds")
# process the score
score_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_OneSamplePerDonorSingleCell.rds")
score_df=score_df %>%
  dplyr::select(colnames(.)[grepl("Donor_id|Sex|^Age$|Annot\\.|_UCell",colnames(.))]) %>%
  mutate(Sex=ifelse(Sex=="Female","F","M")) %>%
  dplyr::rename(donor_id=Donor_id, sex=Sex, age=Age)
colnames(score_df)=gsub("_UCell$","",colnames(score_df))
data.table::fwrite(score_df, "~/Project_PBMCage/Results/Immunity_results/KeyTerms_Ucellw.Top30Genes_OneSamplePerDonorSingleCell.csv.gz", sep="\t")

### Score the pseudobulk_detailed obj
PseudobulkdObj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
subset_Ucell=AddModuleScore_UCell(PseudobulkdObj, features=genes_perorientation_persupercl)
SCORE_DF=subset_Ucell[[]]
saveRDS(SCORE_DF, "~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_pseudobulkDetailed.rds")
# process the score
score_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_pseudobulkDetailed.rds")
score_df=score_df %>%
  dplyr::select(-orig.ident) %>%
  mutate(sex=ifelse(sex=="Female","F","M"))
colnames(score_df)=gsub("_UCell$","",colnames(score_df))
data.table::fwrite(score_df, "~/Project_PBMCage/Results/Immunity_results/KeyTerms_top30genes_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")

#####################################



# ### Analyze the results of Ucell scoring of the 12 Key terms based on sc obj
# #####################################
# ###
# 
# ### Make the df with the analysis results on sc obj
# CommonAlteredGenes_score=data.table::fread("~/Project_PBMCage/Results/Immunity_results/KeyTerms_Ucellw.Top30Genes_OneSamplePerDonorSingleCell.csv.gz", sep="\t")
# CommonAlteredGenes_score=CommonAlteredGenes_score %>%
#   mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
#   mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
#   dplyr::rename(agecut=agecut.half)
# 
# # calculate the age correlation of the terms in each rough celltype, split sex
# cor_terms=CommonAlteredGenes_score %>% dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","cell_id","agecut")))
# cor_terms=cor_terms %>% split(.$Annot.rough) %>% lapply(., function(df) df %>% split(.$sex)) %>%
#   lapply(., function(x) lapply(x, function(y) y %>% dplyr::select(-c(sex, Annot.rough))))
# term_cor_results=lapply(cor_terms, 
#                         function(x) lapply(x,
#                                            function(y) lapply(2:ncol(y), function(idx) cor.test(y[["age"]], y[[idx]], method="spearman", exact=F))))
# 
# rho=lapply(term_cor_results, function(x) lapply(x, function(y) lapply(y, function(z) z$estimate))) %>% unlist()
# names_celltypes=names(term_cor_results)
# names_sex=names(term_cor_results[[1]])
# names_superclusters=colnames(cor_terms[[1]][[1]]) %>% .[.!="age"]
# names_celltypes_all=rep(names_celltypes, each=length(names_sex)*length(names_superclusters))
# names_sex_all=rep(rep(names_sex, each=length(names_superclusters)), times=length(names_celltypes))
# names_superclusters_all=rep(names_superclusters, times=length(names_sex)*length(names_celltypes))
# rho_df=data.frame(group=names(rho), rho=rho)
# rho_df$celltypes=names_celltypes_all; rho_df$sex=names_sex_all; rho_df$superclusters=names_superclusters_all
# rho_df=rho_df %>% dplyr::select(-group)
# 
# rho_pval=lapply(term_cor_results, function(x) lapply(x, function(y) lapply(y, function(z) z$p.value))) %>% unlist()
# rho_pval_df=data.frame(group=names(rho_pval), rho_pval=rho_pval)
# rho_pval_df$celltypes=names_celltypes_all; rho_pval_df$sex=names_sex_all; rho_pval_df$superclusters=names_superclusters_all
# rho_pval_df=rho_pval_df %>% dplyr::select(-group)
# 
# rho_value_pval_merged=inner_join(rho_df, rho_pval_df, by=c("celltypes","sex","superclusters")) %>%
#   dplyr::relocate(c("celltypes","sex","superclusters","rho","rho_pval"))
# 
# data.table::fwrite(rho_value_pval_merged, "~/Project_PBMCage/Results/PBMC_results/KeyTerms_Ucellw.Top30Genes_scobj_corResults.csv.gz", sep="\t")
# 
# # calculate the age correlation of the terms in each rough celltype, regardless of sex
# cor_terms=CommonAlteredGenes_score %>% dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","cell_id","agecut","sex")))
# cor_terms=cor_terms %>% split(.$Annot.rough) %>%
#   lapply(., function(x) x %>% dplyr::select(-c(Annot.rough)))
# term_cor_results=lapply(cor_terms, 
#                         function(y) lapply(2:ncol(y), function(idx) cor.test(y[["age"]], y[[idx]], method="spearman", exact=F)))
# 
# rho=lapply(term_cor_results, function(y) lapply(y, function(z) z$estimate)) %>% unlist()
# names_celltypes=names(term_cor_results)
# names_superclusters=colnames(cor_terms[[1]]) %>% .[.!="age"]
# names_celltypes_all=rep(names_celltypes, each=length(names_superclusters))
# names_superclusters_all=rep(names_superclusters, times=length(names_celltypes))
# rho_df=data.frame(group=names(rho), rho=rho)
# rho_df$celltypes=names_celltypes_all; rho_df$superclusters=names_superclusters_all
# rho_df=rho_df %>% dplyr::select(-group)
# 
# rho_pval=lapply(term_cor_results, function(y) lapply(y, function(z) z$p.value)) %>% unlist()
# rho_pval_df=data.frame(group=names(rho_pval), rho_pval=rho_pval)
# rho_pval_df$celltypes=names_celltypes_all; rho_pval_df$superclusters=names_superclusters_all
# rho_pval_df=rho_pval_df %>% dplyr::select(-group)
# 
# rho_value_pval_merged=inner_join(rho_df, rho_pval_df, by=c("celltypes","superclusters")) %>%
#   dplyr::relocate(c("celltypes","superclusters","rho","rho_pval"))
# 
# data.table::fwrite(rho_value_pval_merged, "~/Project_PBMCage/Results/PBMC_results/KeyTerms_Ucellw.Top30Genes_scobj_corResults_regardlessSex.csv.gz", sep="\t")
# 
# ### Plot the heatmap split sex
# 
# #####################################



### Analyze the results of Ucell scoring of the 12 Key terms based on pseudobulk_detailed obj
#####################################
###

library(ComplexHeatmap)

### Make the df with the analysis results on pseudobulk_detailed obj
CommonAlteredGenes_score=data.table::fread("~/Project_PBMCage/Results/Immunity_results/KeyTerms_top30genes_UcellScore_pseudobulkDetailed.csv.gz", sep="\t")
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
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
  dplyr::rename(agecut=agecut.half)

# # calculate the age correlation of the terms in each rough celltype, split sex

# # calculate the age correlation of the terms in each rough celltype, regardless of sex

### Plot the substract of positively-negatively
colnames(CommonAlteredGenes_score)[c(which(colnames(CommonAlteredGenes_score)=="anti-virus__negatively"),
                                     which(colnames(CommonAlteredGenes_score)=="anti-virus__positively"))]=
  c("anti-virus__positively", "anti-virus__negatively") # switch the anti-virus orientation because the the original meaning of the terms are "regulation by host of virus life cycle"

rho_value_pval_merged_subset=CommonAlteredGenes_score %>%
  tidyr::pivot_longer(cols=colnames(.)[grepl("__",colnames(.))], names_to="superclusters", values_to="score") %>%
  mutate(orientation=gsub(".*__","",superclusters)) %>%
  mutate(superclusters=gsub("__.*","",superclusters)) %>%
  tidyr::pivot_wider(names_from="orientation", values_from="score") %>%
  mutate(score_sub=positively-negatively) %>%
  group_by(donor_id, age, Annot.rough, sex, superclusters) %>%
  summarize_at("score_sub", mean) %>%
  tidyr::pivot_wider(names_from="superclusters", values_from="score_sub")
# calculate the age correlation of the terms in each rough celltype, regardless of sex
cor_terms=rho_value_pval_merged_subset %>% ungroup() %>%
  dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","cell_id","agecut","sex")))
cor_terms=cor_terms %>% split(.$Annot.rough) %>%
  lapply(., function(x) x %>% dplyr::select(-c(Annot.rough)))
term_cor_results=lapply(cor_terms, 
                        function(y) lapply(2:ncol(y), function(idx) cor.test(y[["age"]], y[[idx]], method="spearman", exact=F)))

rho=lapply(term_cor_results, function(y) lapply(y, function(z) z$estimate)) %>% unlist()
names_celltypes=names(term_cor_results)
names_superclusters=colnames(cor_terms[[1]]) %>% .[.!="age"]
names_celltypes_all=rep(names_celltypes, each=length(names_superclusters))
names_superclusters_all=rep(names_superclusters, times=length(names_celltypes))
rho_df=data.frame(group=names(rho), rho=rho)
rho_df$celltypes=names_celltypes_all; rho_df$superclusters=names_superclusters_all
rho_df=rho_df %>% dplyr::select(-group)

rho_pval=lapply(term_cor_results, function(y) lapply(y, function(z) z$p.value)) %>% unlist()
rho_pval_df=data.frame(group=names(rho_pval), rho_pval=rho_pval)
rho_pval_df$celltypes=names_celltypes_all; rho_pval_df$superclusters=names_superclusters_all
rho_pval_df=rho_pval_df %>% dplyr::select(-group)

rho_value_pval_merged=inner_join(rho_df, rho_pval_df, by=c("celltypes","superclusters")) %>%
  dplyr::relocate(c("celltypes","superclusters","rho","rho_pval")) %>%
  mutate(rho_pval=ifelse(rho_pval<0.0001,"****",
                         ifelse(rho_pval>=0.0001 & rho_pval<0.001, "***",
                                ifelse(rho_pval>=0.001 & rho_pval<0.01, "**",
                                       ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", "")))))
# make the df
columns_reorder=c("nucl. metab.","prot. metab.","ribosome syn.","cellular org.","energy metab.","cell division","autophagy","prog. death",
                  "cell response","anti-virus","leuk. devel.","immune proc.")
rows_reorder=names(table(rho_value_pval_merged$celltypes))[c(1,2,3,6,8,4,5,7)]
all_cell_cor=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
  dplyr::select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="superclusters", values_from="rho") %>%
  tibble::column_to_rownames("celltypes") %>%
  .[rows_reorder,columns_reorder[!grepl("ribosome syn\\.",columns_reorder)]]
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
  dplyr::select(-rho) %>%
  tidyr::pivot_wider(names_from="superclusters", values_from="rho_pval") %>%
  tibble::column_to_rownames("celltypes") %>%
  .[rows_reorder,columns_reorder[!grepl("ribosome syn\\.",columns_reorder)]] %>%
  as.matrix()
# plot
range(all_cell_cor) # check
col_fun_sub=circlize::colorRamp2(c(-0.45, 0, 0.45), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title=expression(Spearman~rho),
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Terekhova et al. (2023) dataset", column_title_gp=gpar(fontsize=11),
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(all_cell_fdr, i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_pos.neg.Sub.pdf", height=2, width=4.5)
draw(ht_final)
dev.off()


# ### Plot the relative b value of lm of pos/neg genes
# colnames(CommonAlteredGenes_score)[c(which(colnames(CommonAlteredGenes_score)=="anti-virus__negatively"),
#                                      which(colnames(CommonAlteredGenes_score)=="anti-virus__positively"))]=
#   c("anti-virus__positively", "anti-virus__negatively") # switch the anti-virus orientation because the the original meaning of the terms are "regulation by host of virus life cycle"
# 
# rho_value_pval_merged_subset_positive=CommonAlteredGenes_score %>%
#   tidyr::pivot_longer(cols=colnames(.)[grepl("__",colnames(.))], names_to="superclusters", values_to="score") %>%
#   mutate(orientation=gsub(".*__","",superclusters)) %>%
#   mutate(superclusters=gsub("__.*","",superclusters)) %>%
#   tidyr::pivot_wider(names_from="orientation", values_from="score") %>%
#   mutate(score_sub=positively) %>%
#   group_by(donor_id, age, Annot.rough, sex, superclusters) %>%
#   summarize_at("score_sub", mean) %>%
#   tidyr::pivot_wider(names_from="superclusters", values_from="score_sub")
# score_pos_bvalue=rho_value_pval_merged_subset_positive %>% ungroup() %>% split(.$Annot.rough) %>% 
#   lapply(., function(df) df %>% dplyr::select(-c(donor_id, sex, Annot.rough)) %>% arrange(age)) %>%
#   lapply(., function(df) lapply(2:ncol(df), function(x) lm(df[[x]]~df[[1]]) %>% .$coefficients %>% .[[2]])) %>%
#   lapply(., function(x) {df_=as.data.frame(x); colnames(df_)=colnames(rho_value_pval_merged_subset_positive) %>% .[5:ncol(rho_value_pval_merged_subset_positive)]; df_})
# score_pos_bvalue=lapply(1:length(score_pos_bvalue), function(x) score_pos_bvalue[[x]] %>% mutate(Annot.rough=names(score_pos_bvalue)[x])) %>%
#   data.table::rbindlist(.)
# 
# rho_value_pval_merged_subset_negative=CommonAlteredGenes_score %>%
#   tidyr::pivot_longer(cols=colnames(.)[grepl("__",colnames(.))], names_to="superclusters", values_to="score") %>%
#   mutate(orientation=gsub(".*__","",superclusters)) %>%
#   mutate(superclusters=gsub("__.*","",superclusters)) %>%
#   tidyr::pivot_wider(names_from="orientation", values_from="score") %>%
#   mutate(score_sub=negatively) %>%
#   group_by(donor_id, age, Annot.rough, sex, superclusters) %>%
#   summarize_at("score_sub", mean) %>%
#   tidyr::pivot_wider(names_from="superclusters", values_from="score_sub")
# score_neg_bvalue=rho_value_pval_merged_subset_negative %>% ungroup() %>% split(.$Annot.rough) %>% 
#   lapply(., function(df) df %>% dplyr::select(-c(donor_id, sex, Annot.rough)) %>% arrange(age)) %>%
#   lapply(., function(df) lapply(2:ncol(df), function(x) lm(df[[x]]~df[[1]]) %>% .$coefficients %>% .[[2]])) %>%
#   lapply(., function(x) {df_=as.data.frame(x); colnames(df_)=colnames(rho_value_pval_merged_subset_negative) %>% .[5:ncol(rho_value_pval_merged_subset_negative)]; df_})
# score_neg_bvalue=lapply(1:length(score_neg_bvalue), function(x) score_neg_bvalue[[x]] %>% mutate(Annot.rough=names(score_neg_bvalue)[x])) %>%
#   data.table::rbindlist(.)
# score_relative_bvalue=score_pos_bvalue[,1:(ncol(score_pos_bvalue)-1)]/score_neg_bvalue[,1:(ncol(score_neg_bvalue)-1)]
# score_relative_bvalue=as.matrix(score_relative_bvalue)
# rownames(score_relative_bvalue)=score_neg_bvalue$Annot.rough
# # make the df
# columns_reorder=c("nucl. metab.","prot. metab.","ribosome syn.","cellular org.","energy metab.","cell division","autophagy","prog. death",
#                   "cell response","anti-virus","leuk. devel.","immune proc.")
# rows_reorder=names(table(rownames(score_relative_bvalue)))[c(1,2,3,6,8,4,5,7)]
# all_cell_cor=score_relative_bvalue %>%
#   .[rows_reorder,columns_reorder[!grepl("ribosome syn\\.",columns_reorder)]]
# all_cell_cor[is.na(all_cell_cor)]=0
# # plot
# range(all_cell_cor) # check
# col_fun_sub=circlize::colorRamp2(c(-10, 0, 10), c("dodgerblue3", "white", "brown3"))
# # ht_=
#   Heatmap(all_cell_cor,
#           name="mat", show_heatmap_legend=T,
#           heatmap_legend_param=list(
#             title=expression(Spearman~rho),
#             legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
#             direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
#           ),
#           rect_gp=gpar(col="white", lwd=0.5),
#           col=col_fun_sub,
#           cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
#           column_names_rot=90,
#           show_column_names=T,
#           column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
#           column_title="Terekhova et al. (2023) dataset", column_title_gp=gpar(fontsize=11),
#           # layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
#           #   v=pindex(all_cell_fdr, i, j)
#           #   grid.text(v, x, y, rot=0, just="center",
#           #             gp=gpar(fontsize=10, fontface="bold"))
#           # },
#           width=ncol(all_cell_cor)*unit(7, "mm"),
#           height=nrow(all_cell_cor)*unit(5.5, "mm")
#   )
# ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
# ht_final=draw(ht_, heatmap_legend_side="left")

####################################



### Ucell scoring of the top30 age-cor.genes from 12 superclusters in each celltype at the rough level
# ... (UcellScore_pseudobulkDetailed_negposUcell)
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(ComplexHeatmap)

### Extract the age-corr. genes related to the superclusters
exist_genes_in_terms=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters_posAndnegRegulate.rds")
# get the gene list
superclusters_=unique(exist_genes_in_terms$superclusters)
genes_persupercl=list()
for (i in 1:length(superclusters_)) {
  # switch the anti-virus__negative and __positive, as the term is actually "virus control by host"
  if (superclusters_[i]=="anti-virus") {
    exist_genes_in_terms_subset=exist_genes_in_terms %>% 
      subset(superclusters==superclusters_[i]) %>%
      mutate(term_orientation=ifelse(term_orientation=="positively","negatively","positively")) %>%
      subset(rho_pval<0.05) %>% slice_max(abs(rho), n=50) %>%
      split(.$term_orientation) %>%
      lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe())
  } else {
    exist_genes_in_terms_subset=exist_genes_in_terms %>% 
      subset(superclusters==superclusters_[i]) %>%
      subset(rho_pval<0.05) %>% slice_max(abs(rho), n=50) %>%
      split(.$term_orientation) %>%
      lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe())
  }
  
  names(exist_genes_in_terms_subset)=paste0(superclusters_[i], "__",names(exist_genes_in_terms_subset))
  genes_persupercl=c(genes_persupercl, exist_genes_in_terms_subset)
}

### Add positive and negative signs for UCell scoring
genes_persupercl_bothorientation=list()
for (i in 1:length(genes_persupercl)) {
  if (grepl("negatively$",names(genes_persupercl)[i])) {
    genes_persupercl[[i]]=paste0(genes_persupercl[[i]],"-")
  }
  if (grepl("positively$",names(genes_persupercl)[i])) {
    genes_persupercl[[i]]=paste0(genes_persupercl[[i]],"+")
  }
  # merge the pos and neg for a supercluster
  if (i%%2==0) genes_persupercl_bothorientation[[i/2]]=c(genes_persupercl[[i-1]],genes_persupercl[[i]])
}
names(genes_persupercl_bothorientation)=superclusters_

### Score the pseudobulk_detailed obj
PseudobulkdObj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
subset_Ucell=AddModuleScore_UCell(PseudobulkdObj, features=genes_persupercl_bothorientation)
SCORE_DF=subset_Ucell[[]]
saveRDS(SCORE_DF, "~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_pseudobulkDetailed_negposUcell.rds")
# process the score
score_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/KeyTerms_Ucellw.Top30Genes_pseudobulkDetailed_negposUcell.rds")
score_df=score_df %>%
  dplyr::select(-orig.ident) %>%
  mutate(sex=ifelse(sex=="Female","F","M"))
colnames(score_df)=gsub("_UCell$","",colnames(score_df))
data.table::fwrite(score_df, "~/Project_PBMCage/Results/Immunity_results/KeyTerms_top30genes_UcellScore_pseudobulkDetailed_negposUcell.csv.gz", sep="\t")

#####################################



### Analyze the results (UcellScore_pseudobulkDetailed_negposUcell) of Ucell scoring of the 12 Key terms based on pseudobulk_detailed obj at rough
#####################################
###

library(ComplexHeatmap)

### Make the df with the analysis results on pseudobulk_detailed obj
CommonAlteredGenes_score=data.table::fread("~/Project_PBMCage/Results/Immunity_results/KeyTerms_top30genes_UcellScore_pseudobulkDetailed_negposUcell.csv.gz", sep="\t")
CommonAlteredGenes_score=CommonAlteredGenes_score %>%
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
  dplyr::rename(agecut=agecut.half)

# calculate the age correlation of the terms in each rough celltype, regardless of sex
cor_terms=CommonAlteredGenes_score %>% ungroup() %>%
  dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","cell_id","agecut","sex")))
cor_terms=cor_terms %>% split(.$Annot.rough) %>%
  lapply(., function(x) x %>% dplyr::select(-c(Annot.rough)))
term_cor_results=lapply(cor_terms, 
                        function(y) lapply(2:ncol(y), function(idx) cor.test(y[["age"]], y[[idx]], method="spearman", exact=F)))

rho=lapply(term_cor_results, function(y) lapply(y, function(z) z$estimate)) %>% unlist()
names_celltypes=names(term_cor_results)
names_superclusters=colnames(cor_terms[[1]]) %>% .[.!="age"]
names_celltypes_all=rep(names_celltypes, each=length(names_superclusters))
names_superclusters_all=rep(names_superclusters, times=length(names_celltypes))
rho_df=data.frame(group=names(rho), rho=rho)
rho_df$celltypes=names_celltypes_all; rho_df$superclusters=names_superclusters_all
rho_df=rho_df %>% dplyr::select(-group)

rho_pval=lapply(term_cor_results, function(y) lapply(y, function(z) z$p.value)) %>% unlist()
rho_pval_df=data.frame(group=names(rho_pval), rho_pval=rho_pval)
rho_pval_df$celltypes=names_celltypes_all; rho_pval_df$superclusters=names_superclusters_all
rho_pval_df=rho_pval_df %>% dplyr::select(-group)

rho_value_pval_merged=inner_join(rho_df, rho_pval_df, by=c("celltypes","superclusters")) %>%
  dplyr::relocate(c("celltypes","superclusters","rho","rho_pval")) %>%
  mutate(rho_pval=ifelse(rho_pval<0.0001,"****",
                         ifelse(rho_pval>=0.0001 & rho_pval<0.001, "***",
                                ifelse(rho_pval>=0.001 & rho_pval<0.01, "**",
                                       ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", "")))))

### Plot the heatmap regardless sex
columns_=names(table(gsub("__.*","",rho_value_pval_merged$superclusters)))
columns_reorder=columns_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
rows_=names(table(rho_value_pval_merged$celltypes))
rows_reorder=rows_[c(1,2,3,6,8,4,5,7)]
# make the df
all_cell_cor=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
  dplyr::select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="superclusters", values_from="rho") %>%
  tibble::column_to_rownames("celltypes") %>%
  .[rows_reorder,columns_reorder]
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
  dplyr::select(-rho) %>%
  tidyr::pivot_wider(names_from="superclusters", values_from="rho_pval") %>%
  tibble::column_to_rownames("celltypes") %>%
  .[rows_reorder,columns_reorder] %>%
  apply(., c(1,2), function(x) ifelse(!is.na(x), x, "")) %>%
  as.matrix()
# plot
range(all_cell_cor) # check
col_fun=circlize::colorRamp2(c(-0.2, 0, 0.2), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(4.5, "mm"),
          height=nrow(all_cell_cor)*unit(7, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(all_cell_fdr, i, j)
            grid.text(v, x, y, rot=90, just="top",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          column_gap=unit(1, "mm"), column_title=NULL,
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_final=draw(ht_, annotation_legend_list=lgd, annotation_legend_side="left")


# 
# ### Plot the substract of positively-negatively
# colnames(CommonAlteredGenes_score)[c(which(colnames(CommonAlteredGenes_score)=="anti-virus__negatively"),
#                                      which(colnames(CommonAlteredGenes_score)=="anti-virus__positively"))]=
#   c("anti-virus__positively", "anti-virus__negatively") # switch the anti-virus orientation because the the original meaning of the terms are "regulation by host of virus life cycle"
# 
# rho_value_pval_merged_subset=CommonAlteredGenes_score %>%
#   tidyr::pivot_longer(cols=colnames(.)[grepl("__",colnames(.))], names_to="superclusters", values_to="score") %>%
#   mutate(orientation=gsub(".*__","",superclusters)) %>%
#   mutate(superclusters=gsub("__.*","",superclusters)) %>%
#   tidyr::pivot_wider(names_from="orientation", values_from="score") %>%
#   mutate(score_sub=positively-negatively) %>%
#   group_by(donor_id, age, Annot.rough, sex, superclusters) %>%
#   summarize_at("score_sub", mean) %>%
#   tidyr::pivot_wider(names_from="superclusters", values_from="score_sub")
# # calculate the age correlation of the terms in each rough celltype, regardless of sex
# cor_terms=rho_value_pval_merged_subset %>% ungroup() %>%
#   dplyr::select(-any_of(c("donor_id","Annot.detailed","Annot.inter","cell_id","agecut","sex")))
# cor_terms=cor_terms %>% split(.$Annot.rough) %>%
#   lapply(., function(x) x %>% dplyr::select(-c(Annot.rough)))
# term_cor_results=lapply(cor_terms, 
#                         function(y) lapply(2:ncol(y), function(idx) cor.test(y[["age"]], y[[idx]], method="spearman", exact=F)))
# 
# rho=lapply(term_cor_results, function(y) lapply(y, function(z) z$estimate)) %>% unlist()
# names_celltypes=names(term_cor_results)
# names_superclusters=colnames(cor_terms[[1]]) %>% .[.!="age"]
# names_celltypes_all=rep(names_celltypes, each=length(names_superclusters))
# names_superclusters_all=rep(names_superclusters, times=length(names_celltypes))
# rho_df=data.frame(group=names(rho), rho=rho)
# rho_df$celltypes=names_celltypes_all; rho_df$superclusters=names_superclusters_all
# rho_df=rho_df %>% dplyr::select(-group)
# 
# rho_pval=lapply(term_cor_results, function(y) lapply(y, function(z) z$p.value)) %>% unlist()
# rho_pval_df=data.frame(group=names(rho_pval), rho_pval=rho_pval)
# rho_pval_df$celltypes=names_celltypes_all; rho_pval_df$superclusters=names_superclusters_all
# rho_pval_df=rho_pval_df %>% dplyr::select(-group)
# 
# rho_value_pval_merged=inner_join(rho_df, rho_pval_df, by=c("celltypes","superclusters")) %>%
#   dplyr::relocate(c("celltypes","superclusters","rho","rho_pval")) %>%
#   mutate(rho_pval=ifelse(rho_pval<0.0001,"****",
#                          ifelse(rho_pval>=0.0001 & rho_pval<0.001, "***",
#                                 ifelse(rho_pval>=0.001 & rho_pval<0.01, "**",
#                                        ifelse(rho_pval>=0.01 & rho_pval<0.05, "*", "")))))
# # make the df
# columns_reorder=c("nucl. metab.","prot. metab.","ribosome syn.","cellular org.","energy metab.","cell division","autophagy","prog. death",
#                   "cell response","anti-virus","leuk. devel.","immune proc.")
# rows_reorder=names(table(rho_value_pval_merged$celltypes))[c(1,2,3,6,8,4,5,7)]
# all_cell_cor=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
#   dplyr::select(-rho_pval) %>%
#   tidyr::pivot_wider(names_from="superclusters", values_from="rho") %>%
#   tibble::column_to_rownames("celltypes") %>%
#   .[rows_reorder,columns_reorder[!grepl("ribosome syn\\.",columns_reorder)]]
# all_cell_cor[is.na(all_cell_cor)]=0
# 
# all_cell_fdr=rho_value_pval_merged %>% dplyr::select(superclusters, celltypes, rho, rho_pval) %>%
#   dplyr::select(-rho) %>%
#   tidyr::pivot_wider(names_from="superclusters", values_from="rho_pval") %>%
#   tibble::column_to_rownames("celltypes") %>%
#   .[rows_reorder,columns_reorder[!grepl("ribosome syn\\.",columns_reorder)]] %>%
#   as.matrix()
# # plot
# range(all_cell_cor) # check
# col_fun_sub=circlize::colorRamp2(c(-0.45, 0, 0.45), c("dodgerblue3", "white", "brown3"))
# ht_=
#   Heatmap(all_cell_cor,
#           name="mat", show_heatmap_legend=T,
#           heatmap_legend_param=list(
#             title=expression(Spearman~rho),
#             legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
#             direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
#           ),
#           rect_gp=gpar(col="white", lwd=0.5),
#           col=col_fun_sub,
#           cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
#           column_names_rot=90,
#           show_column_names=T,
#           column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
#           column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
#           layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
#             v=pindex(all_cell_fdr, i, j)
#             grid.text(v, x, y, rot=0, just="center",
#                       gp=gpar(fontsize=10, fontface="bold"))
#           },
#           width=ncol(all_cell_cor)*unit(7, "mm"),
#           height=nrow(all_cell_cor)*unit(5.5, "mm")
#   )
# ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
# ht_final=draw(ht_, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_pos.neg.Sub.pdf", height=3, width=4.5)
draw(ht_final)
dev.off()

####################################



### Analyze the overlapping of age-cor. genes among different celltypes at Annot.rough level
#####################################
###

### Load the marker genes and the correlation df
merged_results=data.table::fread("~/Project_PBMCage/Results/Immunity_results/Cell_Markers_atAll3Levels.csv.gz", sep="\t")
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# take marker genes
marker_genes=merged_results %>%
  subset(p_val_adj<0.05 & avg_log2FC>0 & celltype.level=="Annot.rough") %>%
  split(.$Annot.rough) %>%
  .[names(table(merged_results$Annot.rough))] %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))

### Arrange pos.correlation df
Cor_analysis_subset=Cor_analysis_DF %>%
  subset(analysis=="all" & rho_pval<0.05) %>%
  split(.$celltypes) %>%
  .[names(table(Cor_analysis_DF$celltypes))] %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))

### Make the data for plotting
names(Cor_analysis_subset)
all_genes=Reduce(c,Cor_analysis_subset[c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")])
all_genes=all_genes[!duplicated(all_genes)]
data_df=data.frame(gene=all_genes,
                   `B cells`=as.numeric(all_genes %in% Cor_analysis_subset[["B cells"]]),
                   `CD4T cells`=as.numeric(all_genes %in% Cor_analysis_subset[["CD4T cells"]]),
                   `CD8T cells`=as.numeric(all_genes %in% Cor_analysis_subset[["CD8T cells"]]),
                   `NK cells`=as.numeric(all_genes %in% Cor_analysis_subset[["NK cells"]]),
                   `OtherT`=as.numeric(all_genes %in% Cor_analysis_subset[["OtherT"]]))
group=colnames(data_df)[!grepl("gene",colnames(data_df))]

### Plot
library(ComplexUpset)
# upset_plot=
  upset(
    data_df, group,
    sort_intersections_by=c('degree', 'cardinality'),
    # min_size=10,
    base_annotations=list(
      'Intersection size'=intersection_size(
        text_mapping=aes(
          label=
            ifelse(!!get_size_mode('intersect')/!!get_size_mode('union')==1,
                   " ",
                   paste0(' ',round(!!get_size_mode('intersect')/!!get_size_mode('union')*100, 1), '%'))
        ),
        color="black", fill="transparent", width=0.75, size=0.25,
        text=list(size=3, angle=90, hjust=0, vjust=0.5)) +
        theme(axis.text.y=element_text(size=9),
              axis.title.y=element_text(size=10))
    ),
    matrix=(
      intersection_matrix(
        geom=geom_point(shape='square', size=2),
        segment=geom_segment(linetype='dotted', linewidth=0.2)
      ) +
        scale_y_discrete(position='right')
    ),
    themes=upset_default_themes(axis.title.x=element_blank(),
                                axis.text.y=element_text(size=10),
                                panel.grid=element_blank()),
    set_sizes=(FALSE)
  )

# ... turns out that CD8T and CD4T shared ~1/3 of age-correlated genes

#####################################



# ### Analyze the overlap between age-correlated genes and highly expressed genes in each celltype level
# #####################################
# ###
#
# library(dplyr)
# library(ggplot2)
#
# ### For Rough level
# genelist_all=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/HighlyExpressedGenes_All.3.CelltypeLevels.rds")$Annot.rough
# celltype_names=names(genelist_all)
# markers_rough=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/Cell_Markers_rough.rds") %>%
#   subset(p_val_adj<0.05 & avg_log2FC>0) %>%
#   split(.$cluster) %>%
#   lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.)) %>%
#   .[celltype_names]
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# # both sex
# upreg_genes=subset(Cor_analysis_DF, rho_pval<0.05 & analysis=="all" & rho>0) %>% split(.$celltypes) %>%
#   lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.)) %>%
#   .[celltype_names]
# downreg_genes=subset(Cor_analysis_DF, rho_pval<0.05 & analysis=="all" & rho<0) %>% split(.$celltypes) %>%
#   lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.)) %>%
#   .[celltype_names]
#
# library(ggVennDiagram)
# p=list()
# for (i in 1:length(upreg_genes)) {
#   x=list(`upreg.`=upreg_genes[[i]], `downreg.`=downreg_genes[[i]], `highly expr.`=genelist_all[[i]])
#
#   p[[i]]=
#     ggVennDiagram(x, label_alpha=0, color="grey50", label_percent_digit=1, edge_size=0.5,
#                   label_size=3.5, set_size=4) +
#     scale_fill_gradient(low="white", high="coral4") +
#     labs(title=names(upreg_genes)[i]) +
#     theme(title=element_text(size=11),
#           legend.position="right",
#           # legend.box="vertical",
#           # legend.direction="vertical",
#           legend.key.width=unit(0.3,"cm")) +
#     guides(fill=guide_colorbar(label.theme=element_text(size=8),
#                                title.theme=element_text(size=9)))
# }
#
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Venn_Rough.pdf", height=3, width=3)
#
#
#
#
# ### Upset, 上调、下调、高表达基因的重合
#
#
#
# ### Upset analysis of the sig.cor.genes in all/females/males in each celltype
# #####################################
# ###
#
# library(dplyr)
#
# ### At the rough level
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# all_celltypes=names(table(Cor_analysis_DF$celltypes))
# upset_plot=list()
# PERCENTs_rough=list()
# for (i in 1:length(all_celltypes)) {
#   Subset_=subset(Cor_analysis_DF, celltypes==all_celltypes[i])
#   all_up=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   all_down=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_up=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_down=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_up=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_down=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#
#   # make the UpSet data
#   library(ComplexHeatmap)
#   library(ggplot2)
#   datalist=list(upreg=upreg_genes[[2]], highexpr=genelist_all[[2]], specific=markers_rough[[2]])
#   m=make_comb_mat(datalist, mode="intersect")
#   m=normalize_comb_mat(m, full_comb_sets=T)
#   set_name(m) # check
#   combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
#   m=m[, combinations_]
#   percents_=c(comb_size(m)[1:3]/set_size(m)[1],
#               comb_size(m)[4:6]/set_size(m)[4])
#   # save the precentages
#   PERCENTs_rough[[i]]=c(percents_, all_celltypes[i])
#   percents_=sprintf("%0.2f%%", percents_*100) # for plotting
#
#   # plot
#   # upset_plot_=
#     HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
#     UpSet(m,
#           top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
#           right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
#           row_names_gp=gpar(fontsize=9))
#
#     # %v%
#     # HeatmapAnnotation(
#     #   empty=anno_empty(border=F, height=unit(2,"mm")),
#     #   text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
#     # )
#
#       library(ggvenn)
#     ggvenn(
#       datalist,
#       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#       stroke_size = 0.5, set_name_size = 4
#     )
#
#
#
#     venn <- Venn(x)
#     data <- process_data(venn)
#     ggplot() +
#       # 1. region count layer
#       geom_sf(data = venn_region(data), aes(fill = count)) +
#       # 2. set edge layer
#       geom_sf(data = venn_setedge(data), aes(color = X),  show.legend = TRUE, size = 2) +
#       # 3. set label layer
#       geom_sf_text(data = venn_setlabel(data), aes(label = name)) +
#       # 4. region label layer
#       geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")),
#                     data = venn_region(data),
#                     size = 3) +
#       scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
#       scale_color_manual(values = c("A" = "yellow","B" ="steelblue",'C' = 'red', 'D' = 'black'),
#                          labels = c('D' = 'D = bdiv_human'))+
#       theme_void()
#
#
#
#     # , specific=markers_rough[[2]]#   ,
#     # set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
#     # comb_order=1:6
#
