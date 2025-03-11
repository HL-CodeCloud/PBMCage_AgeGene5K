##########################################################################
################ Rb, Mt, cell freq across datasets #####################
##########################################################################

### Regardless of celltypes, Rb change in 5 datasets
#####################################
###

### Collect all 5 datasets
meta_info_list=list()

meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned" & Annot.inter!="Other.Eryth") %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  select(donor_id, sex, age, percent.rb)
meta_info$Dataset="Yazar et al. (2022) dataset"
meta_info_list[[1]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")
meta_info=meta_info %>% mutate(Sex=ifelse(Sex=="Male","M","F")) %>%
  dplyr::rename(percent.rb=percent.ribo, donor_id=Tube_id, sex=Sex, age=Age) %>%
  select(donor_id, sex, age, percent.rb)
meta_info$Dataset="Terekhova et al. (2023) dataset"
meta_info_list[[2]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/17PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned" & Annot.inter!="Other.Eryth") %>%
  select(donor_id, sex, age, percent.rb)
meta_info$Dataset="GSE213516"
meta_info_list[[3]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned" & Annot.inter!="Other.Eryth") %>%
  select(donor_id, sex, age, percent.rb)
meta_info$Dataset="GSE158055"
meta_info_list[[4]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/JapanSC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned" & Annot.inter!="Other.Eryth") %>%
  select(donor_id, sex, age, percent.rb)
meta_info$Dataset="SC cohort"
meta_info_list[[5]]=meta_info

meta_info_list=data.table::rbindlist(meta_info_list)
# modify the dataset names to save space when plotting
meta_info_list$Dataset=as.factor(meta_info_list$Dataset)
levels(meta_info_list$Dataset)
levels(meta_info_list$Dataset)=c("GSE158055","GSE213516","SC cohort","Terekhova et al.","Yazar et al.")

# plot all the celltypes in general
meta_info_subset=meta_info_list %>% group_by(Dataset, donor_id, sex, age) %>% summarize_at("percent.rb", mean)
# meta_info_subset$Dataset=forcats::fct_relevel(meta_info_subset$Dataset, c("Yazar et al. (2022) dataset","Terekhova et al. (2023) dataset","GSE213516","GSE158055","SC cohort"))
meta_info_subset$Dataset=forcats::fct_relevel(meta_info_subset$Dataset, c("Yazar et al.","Terekhova et al.","GSE213516","GSE158055","SC cohort"))
plot_rb_all=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  facet_wrap(~Dataset, nrow=2) +
  geom_point(alpha=0.1, size=0.1) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title=" ") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position=c(0.85,0.15),
        strip.text=element_text(size=10),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson", label.sep="\n", label.x.npc=0.1) +
  ggpubr::stat_cor(aes(label=after_stat(r.label), group=sex), method="pearson", show.legend=FALSE, cor.coef.name="r",
                   label.y.npc=0.2, size=3) +
  ggpubr::stat_cor(aes(label=after_stat(p.label), group=sex), method="pearson", show.legend=FALSE, cor.coef.name="r",
                   label.y.npc=0.6, size=3)

pdf("~/Project_PBMCage/Plots/Confirm_5_datasets_General.rb.pdf", height=3.5, width=3.75)
plot(plot_rb_all)
dev.off()

#####################################



### Confirm rb changes in CD8T Annot.detailed in PBMCage and Immunity datasets
#####################################
###

### PBMCage
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned") %>%
  mutate(sex=ifelse(sex=="female","F","M"))
# split Annot.detailed of CD8T
meta_info_subset=meta_info %>% subset(Annot.rough=="CD8T cells") %>% group_by(donor_id, sex, age, Annot.detailed) %>% summarize_at("percent.rb", mean)
plot_rb_sep=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.1, shape=20, size=0.1) +
  facet_wrap(~Annot.detailed, ncol=4) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="pearson", label.y.npc=c(0.25,0.22))

### Immunity
meta_data1=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
meta_data2=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
meta_data_merged=data.table::rbindlist(list(meta_data1, meta_data2))
meta_data_merged=meta_data_merged %>% 
  mutate(Sex=ifelse(Sex=="Male","M","F"))
# split Annot.detailed of CD8T
meta_info_subset=meta_data_merged %>% subset(Annot.rough=="CD8T cells") %>% 
  group_by(Tube_id, Sex, Age, Annot.detailed) %>% summarize_at("percent.ribo", mean) %>%
  dplyr::rename(sex=Sex, age=Age, percent.rb=percent.ribo)
plot_rb_sep2=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.1, shape=20, size=0.1) +
  facet_wrap(~Annot.detailed, ncol=4) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="pearson", label.y.npc=c(0.25,0.22))

plot_rb_sep + plot_rb_sep2

# ... turns out to be inconsistent at the Annot.detailed level

#####################################



### Rb change confirmation in 17_PBMC + CovidCtrl + JapanSC datasets
#####################################
###

### Collect all 3 datasets
meta_info_list=list()

meta_info=data.table::fread("~/Project_PBMCage/Results/17PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned")
meta_info$Dataset="GSE213516"
meta_info_list[[1]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
meta_info$Dataset="GSE158055"
meta_info_list[[2]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/JapanSC_results/Key_meta_info.csv.gz", sep="\t")
meta_info$Dataset="SC cohort"
meta_info_list[[3]]=meta_info

meta_info_list=data.table::rbindlist(meta_info_list)

# Plot CD8T
meta_info_subset=meta_info_list %>% 
  subset(Annot.rough=="CD8T cells") %>%
  group_by(donor_id, sex, age, Dataset) %>% summarize_at("percent.rb", mean)
meta_info_subset$Dataset=forcats::fct_relevel(meta_info_subset$Dataset, c("GSE213516","GSE158055","SC cohort"))
plot_rb_sep=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  facet_wrap(~Dataset, scales="free", nrow=1) +
  geom_point(alpha=0.25) +
  theme_classic() +
  labs(x="age", y="Percent of\nribosomal genes (%)", title="GSE213516, GSE158055, SC cohort", subtitle="<CD8T cells>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.text=element_text(size=9),
        strip.text=element_text(size=10),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson", label.y.npc=c(0.05,0.02)) +
  ggpmisc::stat_correlation(ggpmisc::use_label(c("R","p")), small.p=T, show.legend=FALSE, size=3.5, method="pearson", vstep=0.2)

pdf("~/Project_PBMCage/Plots/Confirm_3_datasets_CD8T_rb.pdf", height=2, width=6.5)
plot(plot_rb_sep)
dev.off()

#####################################



### Freq of OtherT.NKT confirmation in 17_PBMC + CovidCtrl + JapanSC datasets
#####################################
###

### Collect all 3 datasets
meta_info_list=list()

meta_info=data.table::fread("~/Project_PBMCage/Results/17PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned")
meta_info$Dataset="GSE213516"
meta_info_list[[1]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
meta_info$Dataset="GSE158055"
meta_info_list[[2]]=meta_info

meta_info=data.table::fread("~/Project_PBMCage/Results/JapanSC_results/Key_meta_info.csv.gz", sep="\t")
meta_info$Dataset="SC cohort"
meta_info_list[[3]]=meta_info

meta_info_list=data.table::rbindlist(meta_info_list)

# calculate the freq
meta_info_per.rough=meta_info_list %>% group_by(Dataset, donor_id, sex, age, Annot.rough) %>% dplyr::count() %>% mutate(N_per_rough=n) %>% select(-n)
meta_info_per.inter=meta_info_list %>% group_by(Dataset, donor_id, sex, age, Annot.rough, Annot.inter) %>% dplyr::count() %>% mutate(N_per_inter=n) %>% select(-n)
meta_info_per.detailed=meta_info_list %>% group_by(Dataset, donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% dplyr::count() %>% mutate(N_per_detailed=n) %>% select(-n)
meta_info_all=meta_info_list %>% group_by(Dataset, donor_id, sex, age) %>% dplyr::count() %>% mutate(N=n) %>% select(-n)
meta_info_combined=meta_info_per.detailed %>% 
  left_join(., meta_info_per.inter, by=intersect(colnames(.), colnames(meta_info_per.inter))) %>%
  left_join(., meta_info_per.rough, by=intersect(colnames(.), colnames(meta_info_per.rough))) %>%
  left_join(., meta_info_all, by=intersect(colnames(.), colnames(meta_info_all))) %>%
  mutate(percent.detailed=N_per_detailed/N,
         percent.inter=N_per_inter/N,
         percent.rough=N_per_rough/N) %>%
  ungroup()

# combine GSE158055 and GSE213516 as they are ~<70 year-old, and SC cohort is above
meta_info_combined=meta_info_combined %>%
  mutate(combined_dataset=ifelse(Dataset %in% c("GSE158055","GSE213516"), "GSE158055, GSE213516", "SC cohort"))

library(ggpmisc)
plot_detailed=
  ggplot(subset(meta_info_combined, Annot.inter=="OtherT.NKT"), aes(x=age, y=percent.inter, color=sex)) +
  ggrastr::geom_point_rast(size=1, shape=20, alpha=1) +
  facet_wrap(~combined_dataset, nrow=1, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="GSE213516, GSE158055, SC cohort", subtitle="<OtherT.NKT>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="none",
        legend.text=element_text(size=9),
        plot.subtitle=element_text(size=10),
        title=element_text(size=10), 
        plot.title.position="plot") +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpmisc::stat_correlation(aes(group=sex, 
                                label=paste(after_stat(r.label),
                                            after_stat(p.value.label),
                                            sep = "*\", \"*")),
                            show.legend=FALSE, size=3.5, method="spearman", vstep=0.1,
                            small.p=TRUE)

pdf("~/Project_PBMCage/Plots/Confirm_3_datasets_frequency_OtherT.NKT.pdf", width=4.6, height=2.5)
plot(plot_detailed)
dev.off()

#####################################



### Plot the peak and turning point of CD4T.naive and CD8T.naive
#####################################
###

# PBMCage
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% subset(Annot.inter %in% c("CD4T.naive","CD8T.naive")) %>% # have checked CD4T.Tcm and CD4T.Tem in Immunity dataset, inconsistent pattern, so most likely an artifact
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|inter",colnames(.))])) %>%
  subset(!duplicated(.))
plot_PBMCage=
  ggplot(meta_info_combined, aes(x=age, y=percent.inter, color=sex)) +
  facet_wrap(~Annot.inter, scales="free_y", nrow=1) +
  geom_point(size=0.1, shape=20, alpha=0.5) +
  labs(x="age", y="Cell percent (%)", title="Yazar et al. (2022) dataset", subtitle="<T.naive>") +
  geom_smooth(method="loess", fill="grey93", linewidth=0.5) +
  geom_vline(data=meta_info_combined %>% subset(Annot.inter=="CD4T.naive"), aes(xintercept=69), linewidth=0.5, linetype="dotted") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        title=element_text(size=10),
        legend.position="none",
        plot.title.position="plot",
        # plot.margin=unit(c(0.75, 0.25, 0.2, 0.25), "cm"),
        plot.subtitle=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19, linewidth=1, fill="transparent"), title=NULL))
# comfirmed by Immunity
meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% subset(Annot.inter %in% c("CD4T.naive","CD8T.naive")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|inter",colnames(.))])) %>%
  subset(!duplicated(.))
plot_immunity=
  ggplot(meta_info_combined, aes(x=age, y=percent.inter, color=sex)) +
  facet_wrap(~Annot.inter, scales="free_y", nrow=1) +
  geom_point(size=0.1, shape=20, alpha=0.5) +
  labs(x="age", y="Cell percent (%)", title="Terekhova et al. (2023) dataset", subtitle=" ") +
  geom_smooth(method="loess", fill="grey93", linewidth=0.5) +
  geom_vline(data=meta_info_combined %>% subset(Annot.inter=="CD4T.naive"), aes(xintercept=69), linewidth=0.5, linetype="dotted") +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        plot.title.position="plot",
        # plot.margin=unit(c(0.75, 0.25, 0.2, 0.25), "cm"),
        plot.subtitle=element_text(size=10),
        title=element_text(size=10),
        legend.position="none",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19, linewidth=1, fill="transparent")))
plot_both=cowplot::plot_grid(plotlist=list(plot_PBMCage, plot_immunity), axis="lr", align="h", nrow=1)

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_CD4Tnaive_TurningAt69.pdf", height=2.5, width=8.5)
plot(plot_both)
dev.off()

#####################################



### Plot the reproducable celltype freq at the Annot.detailed level in PBMCage vs. Immunity
#####################################
###

library(dplyr)
library(ggplot2)

### Determine the reproducable celltypes
reproducable_celltypes=c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm","CD4T.Tcm.HLA","CD4T.Treg.FOXP3hi","CD8T.naive.LINC02446","CD8T.Tcm.GATA3","CD8T.Tem.GZMK_XCL1","NK.CD56dim.FCER1Glo_KLRC2pos")

### Load
metadata_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% reproducable_celltypes) %>%
  mutate(Dataset="Yazar et al. (2022) dataset")
metadata_Immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
metadata_Immunity=metadata_Immunity %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% reproducable_celltypes) %>%
  mutate(Dataset="Terekhova et al. (2023) dataset")
# combine
meta_info_combined_subset=data.table::rbindlist(list(metadata_PBMCage, metadata_Immunity))
# wrap the celltype labels
meta_info_combined_subset=meta_info_combined_subset %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 10))

### Plot
plot_detailed=
  ggplot(subset(meta_info_combined_subset, Dataset=="Yazar et al. (2022) dataset"), aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, nrow=1, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman") +
  # ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman", label.y.npc=0.5) +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_detailed2=
  ggplot(subset(meta_info_combined_subset, Dataset=="Terekhova et al. (2023) dataset"), aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, nrow=1, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="Terekhova et al. (2023) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman") +
  # ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman", label.y.npc=0.5) +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_detailed_both=cowplot::plot_grid(plotlist=list(plot_detailed, plot_detailed2), axis="lr", align="hv", nrow=2)

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_CD4.CD8.NK_reproducable.celltypes_detailed.pdf", width=16.5, height=4.5)
plot(plot_detailed_both)
dev.off()

#####################################



### Plot the reproducable celltype freq at the Combo of Annot.detailed level in PBMCage vs. Immunity
#####################################
###

library(dplyr)
library(ggplot2)

### Determine the reproducable celltypes
reproducable_celltypes=list(CD8T.Tem.GZMK=c("CD8T.Tem.GZMK_NKG7","CD8T.Tem.GZMKhi"),
                            OtherT.gdT.CMC1=c("OtherT.gdT.CMC1"),
                            CD8T.Tem.GZMB_FCGR3A=c("CD8T.Tem.GZMB_FCGR3A")) # Tem GZMB+, Temra, NKT-like

### Load
metadata_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% unlist(reproducable_celltypes)) %>%
  mutate(Dataset="Yazar et al. (2022) dataset")
metadata_Immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
metadata_Immunity=metadata_Immunity %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% unlist(reproducable_celltypes)) %>%
  mutate(Dataset="Terekhova et al. (2023) dataset")
# combine
meta_info_combined_subset=data.table::rbindlist(list(metadata_PBMCage, metadata_Immunity))

### Plot
plot_CD8T.Tem.GZMK=
  ggplot(meta_info_combined_subset %>%
           subset(Annot.detailed %in% reproducable_celltypes[["CD8T.Tem.GZMK"]]) %>%
           group_by(donor_id, sex, age, Dataset) %>%
           summarize_at("percent.detailed", mean), 
         aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Dataset, nrow=1) +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="<CD8T.Tem.GZMK>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="bottom",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=10)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=" ", title.position="top")) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_OtherT.gdT.CMC1=
  ggplot(subset(meta_info_combined_subset, Annot.detailed %in% reproducable_celltypes[["OtherT.gdT.CMC1"]]), aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Dataset, nrow=2) +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL, subtitle="<OtherT.gdT.CMC1>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="bottom",
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_detailed_both=cowplot::plot_grid(plotlist=list(plot_CD8T.Tem.GZMK, plot_OtherT.gdT.CMC1), axis="lr", align="hv", nrow=2, rel_heights=c(1,1.15))

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_CD8TemGZMKorGDT_reproducable.celltypes_detailed.pdf", width=5, height=3.5)
plot(plot_CD8T.Tem.GZMK)
dev.off()

# check CD8T.Tem.GZMB_FCGR3A, which shows clear change in PBMCage but not Immunity dataset
plot_CD8T.Tem.GZMB_FCGR3A=
  ggplot(meta_info_combined_subset %>%
           subset(Annot.detailed %in% reproducable_celltypes[["CD8T.Tem.GZMB_FCGR3A"]]) %>%
           group_by(donor_id, sex, age, Dataset) %>%
           summarize_at("percent.detailed", mean), 
         aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Dataset, nrow=1) +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL, subtitle="<CD8T.Tem.GZMB_FCGR3A>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

#####################################



### Plot the insignificant/unreproduced celltype freq at the Annot.detailed level in PBMCage vs. Immunity
#####################################
###

library(dplyr)
library(ggplot2)

### Load
metadata_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% 
  subset(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  mutate(Dataset="Yazar et al. (2022) dataset")
metadata_PBMCage_types=unique(metadata_PBMCage$Annot.detailed)

metadata_Immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
metadata_Immunity=metadata_Immunity %>% 
  subset(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  mutate(Dataset="Terekhova et al. (2023) dataset")
metadata_Immunity_types=unique(metadata_Immunity$Annot.detailed)

### Filter the celltypes
reproducable_celltypes=c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm","CD4T.Tcm.HLA","CD4T.Treg.FOXP3hi","CD8T.naive.LINC02446","CD8T.Tcm.GATA3","CD8T.Tem.GZMK_XCL1","NK.CD56dim.FCER1Glo_KLRC2pos")
metadata_shared_celltypes=intersect(metadata_PBMCage_types, metadata_Immunity_types)
metadata_shared_celltypes=metadata_shared_celltypes[!(metadata_shared_celltypes %in% reproducable_celltypes)]

### Combine
meta_info_combined_subset=data.table::rbindlist(list(metadata_PBMCage, metadata_Immunity)) %>%
  subset(Annot.detailed %in% metadata_shared_celltypes) %>%
  mutate(Dataset=stringr::str_wrap(gsub(" dataset","",Dataset), 8))
# subsets
CD4_subset=meta_info_combined_subset %>% subset(grepl("CD4T\\.",Annot.detailed)) %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 10))
CD8_subset=meta_info_combined_subset %>% subset(grepl("CD8T\\.",Annot.detailed)) %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 10))
NK_subset=meta_info_combined_subset %>% subset(grepl("NK\\.",Annot.detailed)) %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 10))

### Plot
plot_CD4T=
  ggplot(CD4_subset, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_grid(Dataset ~ Annot.detailed, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL, subtitle="<CD4T cells>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_CD8T=
  ggplot(CD8_subset, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_grid(Dataset ~ Annot.detailed, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL, subtitle="<CD8T cells>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_NK=
  ggplot(NK_subset, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_grid(Dataset ~ Annot.detailed, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL, subtitle="<NK cells>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11)) +
  xlim(c(20,100)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman") +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

plot_detailed_both=cowplot::plot_grid(plotlist=list(plot_CD4T, plot_CD8T, plot_NK), axis="t", align="v", nrow=3)

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_CD4.CD8.NK_insignif.celltypes_detailed.pdf", width=14, height=8)
plot(plot_detailed_both)
dev.off()

#####################################



### Plot the CD8T.Tem.GZMK_NKG7 and CD8T.Tem.GZMKhi celltype freq in PBMCage vs. Immunity
#####################################
###

library(dplyr)
library(ggplot2)

### Load
metadata_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% 
  subset(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  mutate(Dataset="Yazar et al. (2022) dataset")
metadata_PBMCage_types=unique(metadata_PBMCage$Annot.detailed)

metadata_Immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
metadata_Immunity=metadata_Immunity %>% 
  subset(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  mutate(Dataset="Terekhova et al. (2023) dataset")
metadata_Immunity_types=unique(metadata_Immunity$Annot.detailed)

### Filter the celltypes
CD8T_selected_celltypes=c("CD8T.Tem.GZMK_NKG7","CD8T.Tem.GZMKhi")

### Combine
meta_info_combined_subset=data.table::rbindlist(list(metadata_PBMCage, metadata_Immunity)) %>%
  subset(Annot.detailed %in% CD8T_selected_celltypes) %>%
  mutate(Dataset=stringr::str_wrap(gsub(" dataset","",Dataset), 8)) %>%
  mutate(Dataset.sex=paste0(Dataset,", ",sex))
# subsets
CD8_subset=meta_info_combined_subset %>% subset(grepl("CD8T\\.",Annot.detailed)) %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 10))

### Plot
plot_CD8T_CD8Tem.GZMK.NKG7=
  ggplot(CD8_subset, 
         aes(x=age, y=percent.detailed, color=Dataset.sex)) +
    ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.05) +
    facet_wrap(~Annot.detailed, scales="free_y", nrow=2) +
    theme_classic() +
    labs(x="age", y="Cell percent (%)", title=NULL) +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position="bottom",
          legend.justification="left", 
          legend.box.margin=margin(0, 0, 0, -40),
          legend.spacing.x=unit(0.08, 'cm'),
          plot.subtitle=element_text(size=11),
          strip.text.x=element_text(size=10),
          strip.text.y=element_text(size=10),
          title=element_text(size=11)) +
    xlim(c(20,100)) +
    scale_color_manual(values=c("#fbafaa","#6adadd","#67312d","#005052","#f8766d","#00BFC4")) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    ggpmisc::stat_correlation(vstep=0.1, show.legend=FALSE, size=3.5, method="spearman",
                              ggpmisc::use_label(c("R", "p")), small.p=T) +
    geom_smooth(aes(group=Dataset.sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_CD8T.Tem.GZMK.NKG7_vs_GZMKhi.celltypes_detailed.pdf", width=3.5, height=5)
plot(plot_CD8T_CD8Tem.GZMK.NKG7)
dev.off()

### Analyze the difference in signatures and biological functions/properties between the two subsets
# ... CD8T.Tem.GZMK_NKG7 - Tte(terminal effector T)-like; CD8T.Tem.GZMKhi - Tem-like
library(Seurat)
library(clusterProfiler)
PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
PBMCage_subset=subset(PBMCage_obj, Annot.detailed %in% CD8T_selected_celltypes)
Idents(PBMCage_subset)="Annot.detailed"
CD8Tem_subpop_markers=FindAllMarkers(PBMCage_subset)
CD8Tem_subpop_markers_sorted=CD8Tem_subpop_markers %>% split(.$cluster) %>%
  lapply(., function(df) df %>% arrange(desc(avg_log2FC)))
gene_list=CD8Tem_subpop_markers_sorted %>% lapply(., function(df) df %>% select(gene, avg_log2FC) %>% tibble::deframe())
gse_results=compareCluster(gene_list, fun="gseGO", ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
go_genes=CD8Tem_subpop_markers_sorted %>% lapply(., function(df) df %>% subset(avg_log2FC>1) %>% select(gene) %>% tibble::deframe())
go_results=compareCluster(go_genes, fun="enrichGO", ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")

#####################################



### Plot the B.mem.CD27neg freq in PBMCage vs. Immunity
#####################################
###

library(dplyr)
library(ggplot2)

### Define the cell types to be checked
checked_celltypes=c("B.mem.CD27neg")

### Load
metadata_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
metadata_PBMCage=metadata_PBMCage %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% unlist(checked_celltypes)) %>%
  mutate(Dataset="Yazar et al. (2022) dataset")
metadata_Immunity=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
metadata_Immunity=metadata_Immunity %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  subset(Annot.detailed %in% unlist(checked_celltypes)) %>%
  mutate(Dataset="Terekhova et al. (2023) dataset")
# combine
meta_info_combined_subset=data.table::rbindlist(list(metadata_PBMCage, metadata_Immunity))

### Plot
plot_B.mem.CD27neg=
  ggplot(meta_info_combined_subset %>%
           group_by(donor_id, sex, age, Annot.detailed, Dataset) %>%
           summarize_at("percent.detailed", mean), 
         aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Dataset, nrow=1) +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="<B.mem.CD27neg>\n") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="bottom",
        legend.box.spacing=unit(1,"mm"),
        legend.key.height=unit(0.1,"mm"),
        legend.margin=margin(t=-1,unit="mm"),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=10)) +
  xlim(c(20,100)) +
  coord_cartesian(y=c(0,0.03)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=" ", title.position="top")) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="pearson", label.y=0.025) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25)

pdf("~/Project_PBMCage/Plots/Confirm_2_datasets_B.mem.CD27neg_freq.pdf", width=5.2, height=3.8)
plot(plot_B.mem.CD27neg)
dev.off()

#####################################






##########################################################################
################ Age distribution across datasets #####################
##########################################################################


### Collect all 5 datasets
meta_info_pbmc=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t") %>%
  dplyr::select(donor_id, age) %>%
  subset(!duplicated(.))

meta_info_immunity=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t") %>%
  dplyr::select(Tube_id, Age) %>%
  subset(!duplicated(.))

meta_info_GSE213516=data.table::fread("~/Project_PBMCage/Results/17PBMC_results/Key_meta_info.csv.gz", sep="\t") %>%
  dplyr::select(donor_id, age) %>%
  subset(!duplicated(.))

meta_info_GSE158055=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t") %>%
  dplyr::select(donor_id, age) %>%
  subset(!duplicated(.))

meta_info_SCcohort=data.table::fread("~/Project_PBMCage/Results/JapanSC_results/Key_meta_info.csv.gz", sep="\t") %>%
  dplyr::select(donor_id, age) %>%
  subset(!duplicated(.))

meta_info_3sets=data.table::rbindlist(list(meta_info_GSE213516 %>% mutate(dataset="GSE213516"),
                                           meta_info_GSE158055 %>% mutate(dataset="GSE158055"),
                                           meta_info_SCcohort %>% mutate(dataset="SC cohort")))

### Plot the distribution of ages
plot1=
  ggplot(meta_info_3sets, aes(x=age, fill=dataset)) +
  geom_histogram(bins=100, alpha=0.7) +
  scale_x_continuous(limits=c(18,111), n.breaks=10) +
  scale_fill_manual(values=paletteer::paletteer_d("ggthemes::few_Medium")[2:4]) +
  theme_minimal() +
  labs(subtitle="GSE213516, GSE158055, SC cohort") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11))

plot2=
  ggplot(meta_info_immunity, aes(x=Age)) +
  geom_histogram(bins=100, alpha=0.7) +
  scale_x_continuous(limits=c(18,111), n.breaks=10) +
  theme_minimal() +
  labs(x="age", subtitle="Terekhova et al. (2023) dataset") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11))

plot3=
  ggplot(meta_info_pbmc, aes(x=age)) +
  geom_histogram(bins=100, alpha=0.7) +
  scale_x_continuous(limits=c(18,111), n.breaks=10) +
  theme_minimal() +
  labs(subtitle="Yazar et al. (2022) dataset") +
  theme(plot.background=element_rect(color="transparent", fill="white"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        plot.subtitle=element_text(size=11),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        title=element_text(size=11))

plot_all=cowplot::plot_grid(plotlist=list(plot1, plot2, plot3), align="hv", ncol=1, rel_heights=c(1.3,1,1))

pdf("~/Project_PBMCage/Plots/AgeDistribution_5_datasets.pdf", height=5, width=5)
plot(plot_all)
dev.off()

