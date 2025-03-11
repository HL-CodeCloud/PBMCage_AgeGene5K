
### Plot the alterations of the key TFs with age
#####################################
### 

library(Seurat)
library(dplyr)
library(ComplexHeatmap)

### Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj=NormalizeData(THEObj)

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
TF_expr=FetchData(THEObj, vars=c("Donor_id","Age","Sex","Annot.detailed","Annot.inter","Annot.rough",
                                 acts_list_total), layer="data")
TF_expr=TF_expr %>% tibble::rownames_to_column("cell_id")
data.table::fwrite(TF_expr, "~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")

### Arrange the df
TF_expr=data.table::fread("~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|Age|Sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  dplyr::rename(age=Age) %>%
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
  mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("Donor_id|age|Annot\\.|Sex|cell_id", colnames(.))], mean)

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|Sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
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
            right_annotation=
              rowAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                            col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                            annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                            show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
            cluster_columns=T, cluster_rows=F, show_column_dend=F, show_row_dend=F,
            show_row_names=F,
            row_names_gp=gpar(fontsize=9),
            row_title=names(SCORE_DFs)[i],
            row_title_gp=gpar(fontsize=10), row_title_side="right", row_title_rot=0,
            column_names_gp=gpar(fontsize=9),
            # column_title=names(SCORE_DFs)[i],
            # column_title_gp=gpar(fontsize=10), column_title_side="bottom", column_title_rot=0,
            height=unit(2.5,"mm")*nrow(score_df_t),
            width=unit(5,"mm")*ncol(score_df_t))
}
ht_list_=Reduce(`%v%`, ht_list)
lgd_age=Legend(col_fun=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50")),
               at=1:length(age_shared), 
               labels=rownames(score_df_t)[1:length(age_shared)],
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title=NULL, title_position="topcenter", title_gp=gpar(fontface="bold", fontsize=10, fontfamily="mono"),
               direction="vertical", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, heatmap_legend_side="top", merge_legend=FALSE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_TF_KeyTFs.WholeObj_Ucell_heatmap.pdf", height=5, width=4.5)
draw(whole_ht)
dev.off()



### Plot all the 8 Annot.rough celltypes for NFKB and MYC only
acts_list=list()
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
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
# # OneSamplerPerDonor
# THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
# THEObj=NormalizeData(THEObj)
# TF_expr=FetchData(THEObj, vars=c("Donor_id","Age","Sex","Annot.detailed","Annot.inter","Annot.rough",
#                                  acts_list_total), layer="data")
# TF_expr=TF_expr %>% tibble::rownames_to_column("cell_id")
# # PerTube
# THEObj2=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
# THEObj2=NormalizeData(THEObj2)
# TF_expr2=FetchData(THEObj2, vars=c("Donor_id","Age","Sex","Annot.detailed","Annot.inter","Annot.rough",
#                                  acts_list_total), layer="data")
# TF_expr2=TF_expr2 %>% tibble::rownames_to_column("cell_id")
# # make the df
# data.table::fwrite(data.table::rbindlist(list(TF_expr, TF_expr2)), 
#                    "~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData_allTubes.csv.gz", sep="\t")
TF_expr=data.table::fread("~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData_allTubes.csv.gz", sep="\t")
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|Age|Sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  dplyr::rename(age=Age) %>%
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
#   mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% 
  group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("Donor_id|age|Annot\\.|Sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(agecut, Annot.rough, MYC, NFKB2)

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|Sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
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
              title="Terekhova et al. (2023)\nZ-score",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"),
              direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
              title_gp=gpar(fontface="plain", fontsize=10)
            ), 
            top_annotation=
              HeatmapAnnotation(age=anno_simple(1:length(age_shared), height=unit(0.1, "cm"), width=unit(0.1, "cm"),
                                                col=circlize::colorRamp2(c(1,length(age_shared)), c("white","grey50"))),
                                annotation_name_gp=gpar(fontsize=10, fontfamily="mono"),
                                show_annotation_name=c(rep(F,(length(SCORE_DFs)-1)),T)[i]),
            col=rev(paletteer::paletteer_c("ggthemes::Red-Blue-White Diverging", 50)),
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
            width=unit(0.85,"mm")*ncol(t(score_df_t)))
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_TF_KeyTFs.WholeObj_Ucell_heatmap_NFKBandMYC.pdf", height=2.1, width=5.5)
draw(whole_ht)
dev.off()



### Plot all the 8 Annot.rough celltypes for other TFs than NFKB or MYC
TF_expr=data.table::fread("~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData_allTubes.csv.gz", sep="\t")
Fullnames=c("anti-virus","autophagy","cell division","cell response","cellular org.","energy metab.","immune proc.","leuk. devel.",
            "nucl. metab.","prog. death","prot. metab.","ribosome syn.")
short_names=c("antivirus","autophagy","celldivision","cellresponse","cellularorg","energymetab","immuneproc","leukdevel","nuclmetab","progdeath","protmetab","ribosomesyn")
chosen_TFs=list(TP53=c("TP53"),ESR1=c("ESR1"),E2F=c("E2F2","E2F3","E2F4"),HMGA2=c("HMGA2"),
                `CIITA/RFX`=c("RFXAP","RFXANK","RFX5","CIITA"),AR=c("AR"),
                NFKB1=c("NFKB1"),NFKB2=c("NFKB2"),
                # NFKB=c("NFKB1","NFKB2"),
                KMT2A=c("KMT2A"),
                NR2F1=c("NR2F1"),GABPA=c("GABPA"),ZBTB4=c("ZBTB4"),MYC=c("MYC"))
for (i in 1:length(chosen_TFs)) {
  tempt_value=TF_expr %>% dplyr::select(chosen_TFs[[i]]) %>% rowMeans(.)
  TF_expr$tempt=tempt_value
  colnames(TF_expr)[ncol(TF_expr)]=names(chosen_TFs)[i]
}
TF_expr=TF_expr %>% dplyr::select(c(colnames(.)[grepl("_id|Age|Sex|Annot\\.",colnames(.))], names(chosen_TFs)))
TF_expr=TF_expr %>% 
  dplyr::rename(age=Age) %>%
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
  mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)

# get the ages that occur in all the groups (Annot.rough*gene)
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.rough) %>% 
  summarize_at(colnames(.)[!grepl("Donor_id|age|Annot\\.|Sex|cell_id", colnames(.))], mean) %>%
  dplyr::select(-all_of(c("MYC", "NFKB2")))

age_shared=TF_expr_sumByagecut %>% group_by(Annot.rough, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|Sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
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
              legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
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
               legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"),
               title="age", title_position="topcenter", title_gp=gpar(fontface="plain", fontsize=10),
               direction="horizontal", labels_gp=gpar(fontsize=9)
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(2,"mm")
ht_opt$legend_gap=unit(0.5,"cm")
ht_opt$ANNOTATION_LEGEND_PADDING=unit(5,"mm")
ht_opt$ROW_ANNO_PADDING=unit(3,"mm")
whole_ht=draw(ht_list_, annotation_legend_list=lgd_age, annotation_legend_side="top", 
              heatmap_legend_side="top", merge_legend=TRUE, ht_gap=unit(0.2, "cm"))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_TF_KeyTFs.WholeObj_Ucell_heatmap_Other.than.NFKBorMYC.pdf", height=4, width=2.6)
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
TF_expr=data.table::fread("~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData_allTubes.csv.gz", sep="\t")
TF_expr=subset(TF_expr, Annot.rough %in% c("CD4T cells","CD8T cells","NK cells","OtherT")) %>%
  dplyr::rename(age=Age, sex=Sex)
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
TF_expr_sumByagecut=TF_expr %>% group_by(agecut, Annot.inter) %>% 
  summarize_at(colnames(.)[!grepl("donor_id|age|Annot\\.|sex|cell_id", colnames(.))], mean)
age_shared=TF_expr_sumByagecut %>% group_by(Annot.inter, agecut) %>% count() %>%
  tidyr::pivot_wider(names_from="agecut", values_from="n")
age_shared=age_shared[, colSums(is.na(age_shared))==0]
age_shared=colnames(age_shared)[!grepl("Annot|sex", colnames(age_shared))]
age_shared=age_shared[order(age_shared)]
SCORE_DFs=TF_expr_sumByagecut %>% subset(agecut %in% age_shared) %>% 
  tidyr::pivot_wider(id_cols="agecut", names_from="Annot.inter", values_from="MYC") %>%
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



### Plot MYC across ages at Annot.rough
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

### Plot all the CD4T and CD8T Annot.inter celltypes and all the key TFs
TF_expr=data.table::fread("~/Project_PBMCage/Results/Immunity_results/TF_KeyTFs_PBMCageFetchData_allTubes.csv.gz", sep="\t")
TF_expr=TF_expr %>% 
  mutate(cell_group=ifelse(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells","OtherT"), 
                           "Effector lymphocytes", 
                           "Supportive cells")) %>%
  dplyr::rename(age=Age, sex=Sex)
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

### Plot all the Annot.rough
TF_expr_sumByagecut=TF_expr %>% group_by(Donor_id, age, sex, Annot.rough, cell_group) %>% 
  summarize_at(colnames(.)[!grepl("age|Annot\\.|sex|id|cell_group", colnames(.))], mean)

MYC_Yazar=
ggplot(TF_expr_sumByagecut, aes(x=age, y=MYC, color=Annot.rough)) +
  facet_wrap(~cell_group) +
  coord_cartesian(ylim=c(0,0.5)) +
  ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=Annot.rough), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  ggpubr::stat_cor(aes(group=Annot.rough), method="spearman", 
                   label.y=c(0.5,0.5,0.45,0.45,0.4,0.4,0.35,0.35)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title=NULL, y="Expr. of MYC", x="age") +
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
  coord_cartesian(ylim=c(0,0.5)) +
  ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=Annot.rough), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  ggpubr::stat_cor(aes(group=Annot.rough), method="pearson", 
                   label.y=c(0.5,0.5,0.45,0.45,0.4,0.4,0.35,0.35)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  labs(title=NULL, y="Expr. of NFKB2", x="age") +
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

### Plot all the Annot.inter
TF_expr_sum.inter=TF_expr %>% group_by(Donor_id, age, sex, Annot.rough, Annot.inter, cell_group) %>% 
  summarize_at(colnames(.)[!grepl("age|Annot\\.|sex|_id|cell_group", colnames(.))], mean)
# MYC_Yazar_inter=
  ggplot(TF_expr_sum.inter %>% subset(cell_group=="Effector lymphocytes"), aes(x=age, y=MYC, color=Annot.inter)) +
  # coord_cartesian(ylim=c(0,0.5)) +
  ggrastr::geom_point_rast(size=0.05, shape=20, alpha=0.05) +
  geom_smooth(aes(group=Annot.inter), method="loess", fill="grey93", show.legend=F, linewidth=0.5) +
  ggpubr::stat_cor(aes(group=Annot.inter), method="spearman") +
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

# NFKB2_Yazar=
ggplot(TF_expr_sum.inter %>% subset(Annot.rough=="CD8T cells"), aes(x=age, y=NFKB2, color=Annot.inter)) +
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

### Confirm the age correlation of the kinases predicted at KEA3 based on PBMCage dataset
#####################################
###

library(dplyr)

### Load the predicted kinases
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
cor_df=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
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
col_fun_sub=circlize::colorRamp2(c(-0.4, 0, 0.4), c("dodgerblue3", "white", "brown3"))
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
            width=ncol(all_cell_cor)*unit(4, "mm"),
            height=nrow(all_cell_cor)*unit(5.5, "mm"))
}

ht_list=Reduce(`+`, ht_list)
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
whole_ht=draw(ht_list, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Bulkanalysis_RNAexpr.Cor.wAges_Kinase_in_subprocesses_ExprHeatmap.pdf", height=2.8, width=16)
draw(whole_ht)
dev.off()

#####################################