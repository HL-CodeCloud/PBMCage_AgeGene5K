
### Process PBMCage seurat obj
#####################################
###

object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

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
saveRDS(object, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")

# arrange the metadata
tempt=object[[]]
tempt=tempt[,c(1,2,3,4,5,6,8,10,12,49,50,51,52)]
object@meta.data=tempt
saveRDS(object, "~/Project_PBMCage/raw_data/PBMCage_annotated.rds")

# skip unassigned cells
object=subset(object, Annot.rough!="Unassigned")
saveRDS(object, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

#####################################



### Add cell cycle information for all the cells including unassigned
#####################################
###

library(Seurat)

###
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
object[[]] # percent.rb and percent.mt has been added
object # make sure the seurat_obj has been normalize, HVGs found, and scaled

### Cellcycle scoring
object=CellCycleScoring(object,
                        s.features=cc.genes.updated.2019$s.genes,
                        g2m.features=cc.genes.updated.2019$g2m.genes,
                        set.ident=TRUE)

### Save the mt, rb, and cell-cycle info
meta_info=object@meta.data %>% select(any_of(c("nCount_RNA","nFeature_RNA","dataset","donor_id","sex","age",
                                               "percent.rb","percent.mt","Annot.detailed","Annot.inter","Annot.rough",
                                               "agecut",
                                               "S.Score","G2M.Score","Phase"))) %>%
  mutate(cellid=rownames(.))
data.table::fwrite(meta_info, "~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")

### Plot the mt
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned") %>%
  mutate(sex=ifelse(sex=="female","F","M"))
# plot all the celltypes in general
meta_info_subset=meta_info %>% group_by(donor_id, sex, age) %>% summarize_at("percent.mt", mean)
plot_mt_all=
  ggplot(meta_info_subset, aes(x=age, y=percent.mt, color=sex)) +
  geom_point(alpha=0.25) +
  theme_classic() +
  labs(x="age", y="Percent of mitochondrial genes (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson")
# split Annot.rough
meta_info_subset=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% summarize_at("percent.mt", mean)
plot_mt_sep=
  ggplot(meta_info_subset, aes(x=age, y=percent.mt, color=sex)) +
  geom_point(alpha=0.25) +
  facet_wrap(~Annot.rough, ncol=4) +
  theme_classic() +
  labs(x="age", y="Percent of mitochondrial genes (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson")


### Plot the rb
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned" & Annot.inter!="Other.Eryth") %>%
  mutate(sex=ifelse(sex=="female","F","M"))
# plot all the celltypes in general
meta_info_subset=meta_info %>% group_by(donor_id, sex, age) %>% summarize_at("percent.rb", mean)
plot_rb_all=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.25) +
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
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson")

# split Annot.rough
meta_info_subset=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% summarize_at("percent.rb", mean)
plot_rb_sep=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.1, shape=20, size=0.1) +
  facet_wrap(~Annot.rough, ncol=4) +
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
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson", label.y.npc=c(0.25,0.23)) +
  ggpmisc::stat_correlation(ggpmisc::use_label(c("R","p")), small.p=T, show.legend=FALSE, size=3.5, method="pearson", vstep=0.1, label.y=c(0.25,0.1))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_rb_rough.pdf", height=4, width=6)
plot(plot_rb_sep)
dev.off()

# split Annot.inter of Other cells and OtherT
meta_info_subset=meta_info %>% subset(Annot.rough %in% c("Other cells","OtherT") & Annot.inter!="Other.Mast") %>% group_by(donor_id, sex, age, Annot.inter) %>% summarize_at("percent.rb", mean)
plot_rb_interOther=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.1, shape=20, size=0.1) +
  facet_wrap(~Annot.inter, ncol=3) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="Yazar et al. (2022) dataset", subtitle="<Other cells and OtherT>") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_x_continuous(breaks=c(30,60,90)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson", label.y.npc=c(0.25,0.22)) +
  ggpmisc::stat_correlation(ggpmisc::use_label(c("R","p")), small.p=T, show.legend=FALSE, size=3.5, method="pearson", vstep=0.1, label.y=c(0.25,0.1))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/PBMCage_rb_OtherTandOtherCells_inter.pdf", height=3.5, width=5)
plot(plot_rb_interOther)
dev.off()

# split Annot.detailed of CD8T
# ... in the 3_5_datasets_confirm.R


### Plot the cell-cycle phase
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned") %>%
  mutate(sex=ifelse(sex=="female","F","M"))
# plot all the celltypes in general
meta_info_phase=meta_info %>% group_by(donor_id, sex, age, Phase) %>% count() %>% mutate(PhaseN=n) %>% select(-n)
meta_info_allphase=meta_info %>% group_by(donor_id, sex, age) %>% count() %>% mutate(N=n) %>% select(-n)
meta_info_merge=right_join(meta_info_phase, meta_info_allphase, by=c("donor_id","sex","age")) %>%
  mutate(percent.phase=PhaseN/N)
plot_phase_all=
  ggplot(meta_info_merge, aes(x=age, y=percent.phase, color=sex)) +
  geom_point(alpha=0.25) +
  facet_wrap(~Phase) +
  theme_classic() +
  labs(x="age", y="Percent of cell phase (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson")
# split Annot.rough
meta_info_phase=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Phase) %>% count() %>% mutate(PhaseN=n) %>% select(-n)
meta_info_allphase=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% count() %>% mutate(N=n) %>% select(-n)
meta_info_merge=right_join(meta_info_phase, meta_info_allphase, by=c("donor_id","sex","age","Annot.rough")) %>%
  mutate(percent.phase=PhaseN/N)
plot_phase_sep=
  ggplot(meta_info_merge, aes(x=age, y=percent.phase, color=sex)) +
  geom_point(alpha=0.25) +
  facet_wrap(~Annot.rough+Phase, ncol=3) +
  theme_classic() +
  labs(x="age", y="Percent of cell phase (%)", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="pearson")

#####################################



### Make pseudobulk obj
#####################################
###

library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
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
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
# -- add sex metadata
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("male","female")
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")

# pseudobulk at the inter level
celltypes=names(table(Seurat_obj$Annot.inter))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.inter==celltypes[i])
  pseudo_inter=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("donor_id","age"))
  pseudo_inter$donor_id=gsub("_.*","",colnames(pseudo_inter))
  pseudo_inter$age=as.numeric(gsub(".*_","",colnames(pseudo_inter)))
  pseudo_inter$Annot.inter=celltypes[i]
  PseudoList[[i]]=pseudo_inter
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
tempt_annot_df=pbmc.merged[[]]
celltypes_=unique(tempt_annot_df$Annot.inter)
tempt_annot_df=tempt_annot_df %>%
  mutate(Annot.rough=ifelse(Annot.inter %in% celltypes_[grepl("^B\\.",celltypes_)], "B cells", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD4T\\.",celltypes_)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD8T\\.",celltypes_)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^DC\\.",celltypes_)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Mono\\.",celltypes_)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^NK\\.",celltypes_)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Other\\.",celltypes_)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^OtherT\\.",celltypes_)], "OtherT", Annot.rough))
table(rownames(tempt_annot_df)==colnames(pbmc.merged)) # check
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$Annot.rough, col.name="Annot.rough")
# -- add sex metadata
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("male","female")
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")

# pseudobulk at the detailed level
library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
sex_df=data.frame(donor_id=Seurat_obj$donor_id,
                  sex=Seurat_obj$sex)
sex_df=sex_df[!duplicated(sex_df$donor_id),]
celltypes=names(table(Seurat_obj$Annot.detailed))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.detailed==celltypes[i])
  pseudo_detailed=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("donor_id","age"))
  pseudo_detailed$donor_id=gsub("_.*","",colnames(pseudo_detailed))
  pseudo_detailed$age=as.numeric(gsub(".*_","",colnames(pseudo_detailed)))
  pseudo_detailed$Annot.detailed=celltypes[i]
  PseudoList[[i]]=pseudo_detailed
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
pbmc.merged$Annot.inter=sapply(pbmc.merged$Annot.detailed, function(x) paste0(strsplit(x, "\\.")[[1]][c(1,2)], collapse="."))
table(pbmc.merged$Annot.inter, useNA="ifany")
tempt_annot_df=pbmc.merged[[]]
celltypes_=unique(tempt_annot_df$Annot.inter)
tempt_annot_df=tempt_annot_df %>%
  mutate(Annot.rough=ifelse(Annot.inter %in% celltypes_[grepl("^B\\.",celltypes_)], "B cells", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD4T\\.",celltypes_)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD8T\\.",celltypes_)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^DC\\.",celltypes_)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Mono\\.",celltypes_)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^NK\\.",celltypes_)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Other\\.",celltypes_)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^OtherT\\.",celltypes_)], "OtherT", Annot.rough))
table(rownames(tempt_annot_df)==colnames(pbmc.merged)) # check
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$Annot.rough, col.name="Annot.rough")
# -- add sex metadata
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("male","female")
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")

# pseudobulk at the detailed level (updated for publication)
Seurat_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
sex_df=data.frame(donor_id=Seurat_obj$donor_id,
                  sex=Seurat_obj$sex)
sex_df=sex_df[!duplicated(sex_df$donor_id),]
Seurat_obj[[]]=Seurat_obj[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.prolif.TSHZ2"), "CD4T.prolif", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("CD4T.Tcm.ISG.TSHZ2"), "CD4T.Tcm.ISG", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))
# remove Other.Eryth and Other.Mast as we currently are not interested in the previous one and the second one is too rare
Seurat_obj=subset(Seurat_obj, Annot.detailed!="Other.Eryth" & Annot.detailed!="Other.Mast")
celltypes=names(table(Seurat_obj$Annot.detailed))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.detailed==celltypes[i])
  pseudo_detailed=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("donor_id","age"))
  pseudo_detailed$donor_id=gsub("_.*","",colnames(pseudo_detailed))
  pseudo_detailed$age=as.numeric(gsub(".*_","",colnames(pseudo_detailed)))
  pseudo_detailed$Annot.detailed=celltypes[i]
  PseudoList[[i]]=pseudo_detailed
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
pbmc.merged$Annot.inter=sapply(pbmc.merged$Annot.detailed, function(x) paste0(strsplit(x, "\\.")[[1]][c(1,2)], collapse="."))
table(pbmc.merged$Annot.inter, useNA="ifany")
tempt_annot_df=pbmc.merged[[]]
celltypes_=unique(tempt_annot_df$Annot.inter)
tempt_annot_df=tempt_annot_df %>%
  mutate(Annot.rough=ifelse(Annot.inter %in% celltypes_[grepl("^B\\.",celltypes_)], "B cells", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD4T\\.",celltypes_)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^CD8T\\.",celltypes_)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^DC\\.",celltypes_)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Mono\\.",celltypes_)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^NK\\.",celltypes_)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^Other\\.",celltypes_)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_[grepl("^OtherT\\.",celltypes_)], "OtherT", Annot.rough))
table(rownames(tempt_annot_df)==colnames(pbmc.merged)) # check
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$Annot.rough, col.name="Annot.rough")
# -- add sex metadata
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("male","female")
saveRDS(pbmc.merged, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")

# pseudobulk at the donor_id level
pbmc.merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned") %>%
  dplyr::select(donor_id, sex, age) %>% subset(!duplicated(.))
pseudo_donor=Seurat::AggregateExpression(pbmc.merged, assays="RNA", return.seurat=T, group.by=c("donor_id"))
pseudo_donor[[]]=pseudo_donor[[]] %>%
  tibble::rownames_to_column("donor_id") %>%
  mutate(donor_id=gsub("-","_",donor_id)) %>%
  left_join(., meta_info, by="donor_id")
saveRDS(pseudo_donor, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_donor.rds")

#####################################



# ### Draw boxplot of freq vs ages
# #####################################
# ###
#
# ###
# pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# ID_=pbmc.seu_merged$donor_id
# Age_=pbmc.seu_merged$age
# CellType_=pbmc.seu_merged$Annot.detailed
# cellNum=table(CellType_, ID_)
# cellProp=t(t(cellNum)/rowSums(t(cellNum)))
# data=t(cellProp*100)
# data=as.data.frame(data)
# colnames(data)=c("ID","CellType","Freq")
# ID_Age_match=data.frame(ID=ID_, Age=Age_)
# ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
# data=merge(data, ID_Age_match, by="ID")
# data$CellType=as.factor(data$CellType)
#
# # make a function
# library(ggpubr)
# draw_proportion_perAge=function(df, celltype.drawn) {
#   df.used=subset(data, CellType==celltype.drawn)
#   df.used=df.used %>%
#     dplyr::add_count(Age, name="Age_n") %>%
#     as.data.frame()
#   df.used$Age=as.factor(df.used$Age)
#   p=
#     ggplot(df.used, aes(x=Age, y=Freq)) +
#     geom_boxplot(color="gray80", outlier.size=0.5) +
#     # geom_point(df.used %>% subset(Age_n==1), aes(x=Age, y=Freq)) +
#     theme_classic() +
#     labs(x=NULL, y="Cell percent (%)") +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.y=element_text(size=10),
#           title=element_text(size=10),
#           legend.position="none") +
#     labs(title=paste0("Percent of ", celltype.drawn)) +
#     # geom_smooth(aes(x=as.numeric(Age), y=Freq), method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
#     stat_cor(aes(x=as.numeric(Age), y=Freq), method="kendall",show.legend=FALSE, size=3.5)
#
#   median_1=unlist(lapply(split(df.used$Freq, df.used$Age), median))
#   p1_l=data.frame(Age=names(median_1), Freq=median_1)
#   p=p + geom_line(data=p1_l, aes(x=Age, y=Freq, group=1), linetype="dashed", color="black")
#   plot(p)
# }
#
# # plot cell freq vs. ages for each celltype
# celltypes_ordered=names(table(data$CellType))
# celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B\\.", celltypes_ordered)],
#                       CD4_T=celltypes_ordered[grepl("^CD4T\\.", celltypes_ordered)],
#                       CD8_T=celltypes_ordered[grepl("^CD8T\\.", celltypes_ordered)],
#                       DC=celltypes_ordered[grepl("^DC\\.", celltypes_ordered)],
#                       Monocyte=celltypes_ordered[grepl("^Mono\\.", celltypes_ordered)],
#                       NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
#                       Other_cells=celltypes_ordered[grepl("^Other\\.", celltypes_ordered)],
#                       Other_T_cells=celltypes_ordered[grepl("^OtherT\\.", celltypes_ordered)])
# length_max=max(sapply(celltypes_groups, length))
# plistcowplot=list()
# for (i in 1:length(celltypes_groups)) {
#   current_types=celltypes_groups[[i]]
#   plist=list()
#   for (j in 1:length(current_types)) {
#     current_type=current_types[j]
#     plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
#   }
#   length_of_plist=length(plist)
#   codes=""
#   for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
#   if (length_of_plist<length_max) {for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")}
#   codes=gsub("^, ","",codes)
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=4, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Proportion_vs._age.pdf", width=25, height=16)
# for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
# dev.off()
#
# #####################################
#
#
#
# ### Draw boxplot of freq vs agecuts
# #####################################
# ###
#
# ###
# pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# ID_=pbmc.seu_merged$donor_id
# Age_=pbmc.seu_merged$agecut
# CellType_=pbmc.seu_merged$Annot.detailed
# cellNum=table(CellType_, ID_)
# cellProp=t(t(cellNum)/rowSums(t(cellNum)))
# data=t(cellProp*100)
# data=as.data.frame(data)
# colnames(data)=c("ID","CellType","Freq")
# ID_Age_match=data.frame(ID=ID_, Agecut=Age_)
# ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
# data=merge(data, ID_Age_match, by="ID")
# data$CellType=as.factor(data$CellType)
#
# # make a function
# library(ggpubr)
# draw_proportion_perAge=function(df, celltype.drawn) {
#   df.used=subset(data, CellType==celltype.drawn)
#   df.used=df.used %>%
#     as.data.frame()
#   df.used$Agecut=as.factor(df.used$Agecut)
#   p=
#     ggplot(df.used, aes(x=Agecut, y=Freq)) +
#     geom_boxplot(color="gray80", outlier.size=0.5) +
#     theme_classic() +
#     labs(x=NULL, y="Cell percent (%)") +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.y=element_text(size=10),
#           title=element_text(size=10),
#           legend.position="none") +
#     labs(title=paste0("Percent of ", celltype.drawn)) +
#     # geom_smooth(aes(x=as.numeric(Agecut), y=Freq), method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
#     stat_cor(aes(x=as.numeric(Agecut), y=Freq), method="kendall",show.legend=FALSE, size=3.5)
#
#   median_1=unlist(lapply(split(df.used$Freq, df.used$Agecut), median))
#   p1_l=data.frame(Agecut=names(median_1), Freq=median_1)
#   p=p + geom_line(data=p1_l, aes(x=Agecut, y=Freq, group=1), linetype="dashed", color="black")
#   plot(p)
# }
#
# # plot cell freq vs. ages for each celltype
# celltypes_ordered=names(table(data$CellType))
# celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B\\.", celltypes_ordered)],
#                       CD4_T=celltypes_ordered[grepl("^CD4T\\.", celltypes_ordered)],
#                       CD8_T=celltypes_ordered[grepl("^CD8T\\.", celltypes_ordered)],
#                       DC=celltypes_ordered[grepl("^DC\\.", celltypes_ordered)],
#                       Monocyte=celltypes_ordered[grepl("^Mono\\.", celltypes_ordered)],
#                       NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
#                       Other_cells=celltypes_ordered[grepl("^Other\\.", celltypes_ordered)],
#                       Other_T_cells=celltypes_ordered[grepl("^OtherT\\.", celltypes_ordered)])
# length_max=max(sapply(celltypes_groups, length))
# plistcowplot=list()
# for (i in 1:length(celltypes_groups)) {
#   current_types=celltypes_groups[[i]]
#   plist=list()
#   for (j in 1:length(current_types)) {
#     current_type=current_types[j]
#     plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
#   }
#   length_of_plist=length(plist)
#   codes=""
#   for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
#   if (length_of_plist<length_max) {for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")}
#   codes=gsub("^, ","",codes)
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=4, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Proportion_vs._agecut.pdf", width=25, height=16)
# for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
# dev.off()
#
# #####################################
#
#
#
# ### Draw boxplot of counts vs ages
# #####################################
# ###
#
# ###
# pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# ID_=pbmc.seu_merged$donor_id
# Age_=pbmc.seu_merged$age
# CellType_=pbmc.seu_merged$Annot.detailed
# cellNum=table(CellType_, ID_)
# data=as.data.frame(t(cellNum))
# colnames(data)=c("ID","CellType","Counts")
# ID_Age_match=data.frame(ID=ID_, Age=Age_)
# ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
# data=merge(data, ID_Age_match, by="ID")
# data$CellType=as.factor(data$CellType)
#
# # make a function
# library(ggpubr)
# draw_proportion_perAge=function(df, celltype.drawn) {
#   df.used=subset(data, CellType==celltype.drawn)
#   df.used=df.used %>%
#     dplyr::add_count(Age, name="Age_n") %>%
#     as.data.frame()
#   df.used$Age=as.factor(df.used$Age)
#   p=
#     ggplot(df.used, aes(x=Age, y=Counts)) +
#     geom_boxplot(color="gray80", outlier.size=0.5) +
#     # geom_point(df.used %>% subset(Age_n==1), aes(x=Age, y=Freq)) +
#     theme_classic() +
#     labs(x=NULL, y="Cell counts") +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.y=element_text(size=10),
#           title=element_text(size=10),
#           legend.position="none") +
#     labs(title=paste0("Counts of ", celltype.drawn)) +
#     # geom_smooth(aes(x=as.numeric(Age), y=Counts), method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
#     stat_cor(aes(x=as.numeric(Age), y=Counts), method="kendall",show.legend=FALSE, size=3.5)
#
#   median_1=unlist(lapply(split(df.used$Counts, df.used$Age), median))
#   p1_l=data.frame(Age=names(median_1), Counts=median_1)
#   p=p + geom_line(data=p1_l, aes(x=Age, y=Counts, group=1), linetype="dashed", color="black")
#   plot(p)
# }
#
# # plot cell counts vs. ages for each celltype
# celltypes_ordered=names(table(data$CellType))
# celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B\\.", celltypes_ordered)],
#                       CD4_T=celltypes_ordered[grepl("^CD4T\\.", celltypes_ordered)],
#                       CD8_T=celltypes_ordered[grepl("^CD8T\\.", celltypes_ordered)],
#                       DC=celltypes_ordered[grepl("^DC\\.", celltypes_ordered)],
#                       Monocyte=celltypes_ordered[grepl("^Mono\\.", celltypes_ordered)],
#                       NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
#                       Other_cells=celltypes_ordered[grepl("^Other\\.", celltypes_ordered)],
#                       Other_T_cells=celltypes_ordered[grepl("^OtherT\\.", celltypes_ordered)])
# length_max=max(sapply(celltypes_groups, length))
# plistcowplot=list()
# for (i in 1:length(celltypes_groups)) {
#   current_types=celltypes_groups[[i]]
#   plist=list()
#   for (j in 1:length(current_types)) {
#     current_type=current_types[j]
#     plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
#   }
#   length_of_plist=length(plist)
#   codes=""
#   for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
#   if (length_of_plist<length_max) {for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")}
#   codes=gsub("^, ","",codes)
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=4, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Counts_vs._age.pdf", width=25, height=16)
# for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
# dev.off()
#
# #####################################
#
#
#
# ### Draw boxplot of counts vs agecuts
# #####################################
# ###
#
# ###
# pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# ID_=pbmc.seu_merged$donor_id
# Age_=pbmc.seu_merged$agecut
# CellType_=pbmc.seu_merged$Annot.detailed
# cellNum=table(CellType_, ID_)
# data=as.data.frame(t(cellNum))
# colnames(data)=c("ID","CellType","Counts")
# ID_Age_match=data.frame(ID=ID_, Agecut=Age_)
# ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
# data=merge(data, ID_Age_match, by="ID")
# data$CellType=as.factor(data$CellType)
#
# # make a function
# library(ggpubr)
# draw_proportion_perAge=function(df, celltype.drawn) {
#   df.used=subset(data, CellType==celltype.drawn)
#   df.used=df.used %>%
#     as.data.frame()
#   df.used$Agecut=as.factor(df.used$Agecut)
#   p=
#     ggplot(df.used, aes(x=Agecut, y=Counts)) +
#     geom_boxplot(color="gray80", outlier.size=0.5) +
#     theme_classic() +
#     labs(x=NULL, y="Cell counts") +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.y=element_text(size=10),
#           title=element_text(size=10),
#           legend.position="none") +
#     labs(title=paste0("Counts of ", celltype.drawn)) +
#     # geom_smooth(aes(x=as.numeric(Agecut), y=Counts), method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
#     stat_cor(aes(x=as.numeric(Agecut), y=Counts), method="kendall",show.legend=FALSE, size=3.5)
#
#   median_1=unlist(lapply(split(df.used$Counts, df.used$Agecut), median))
#   p1_l=data.frame(Agecut=names(median_1), Counts=median_1)
#   p=p + geom_line(data=p1_l, aes(x=Agecut, y=Counts, group=1), linetype="dashed", color="black")
#   plot(p)
# }
#
# # plot cell counts vs. ages for each celltype
# celltypes_ordered=names(table(data$CellType))
# celltypes_groups=list(B_cell=celltypes_ordered[grepl("^B\\.", celltypes_ordered)],
#                       CD4_T=celltypes_ordered[grepl("^CD4T\\.", celltypes_ordered)],
#                       CD8_T=celltypes_ordered[grepl("^CD8T\\.", celltypes_ordered)],
#                       DC=celltypes_ordered[grepl("^DC\\.", celltypes_ordered)],
#                       Monocyte=celltypes_ordered[grepl("^Mono\\.", celltypes_ordered)],
#                       NK=celltypes_ordered[grepl("^NK\\.", celltypes_ordered)],
#                       Other_cells=celltypes_ordered[grepl("^Other\\.", celltypes_ordered)],
#                       Other_T_cells=celltypes_ordered[grepl("^OtherT\\.", celltypes_ordered)])
# length_max=max(sapply(celltypes_groups, length))
# plistcowplot=list()
# for (i in 1:length(celltypes_groups)) {
#   current_types=celltypes_groups[[i]]
#   plist=list()
#   for (j in 1:length(current_types)) {
#     current_type=current_types[j]
#     plist[[j]]=draw_proportion_perAge(df=data, celltype.drawn=current_type)
#   }
#   length_of_plist=length(plist)
#   codes=""
#   for (k in 1:length_of_plist) codes=paste0(codes, ", plist[[",k,"]]")
#   if (length_of_plist<length_max) {for (k in (length_of_plist+1):length_max) codes=paste0(codes, ", NULL")}
#   codes=gsub("^, ","",codes)
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=4, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Counts_vs._agecut.pdf", width=25, height=16)
# for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
# dev.off()
#
# #####################################



### Analyze the celltype freq and their correlation with ages, similar to the above "Draw boxplot" steps
#####################################
###

### Load the meta data with celltypes
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned")
# # plot the freq at the Annot.inter level
# meta_info_perinter=meta_info %>% group_by(donor_id, sex, age, Annot.inter) %>% count() %>% mutate(N_per_inter=n) %>% select(-n)
# meta_info_all=meta_info %>% group_by(donor_id, sex, age) %>% count() %>% mutate(N=n) %>% select(-n)
# meta_info_merge=right_join(meta_info_perinter, meta_info_all, by=c("donor_id","sex","age")) %>%
#   mutate(percent.inter=N_per_inter/N) %>%
#   mutate(sex=ifelse(sex=="female","F","M"))
# plot_inter=
#   ggplot(meta_info_merge, aes(x=age, y=percent.inter, color=sex)) +
#   geom_point(alpha=0.25) +
#   facet_wrap(~Annot.inter, ncol=5, scales="free_y") +
#   theme_classic() +
#   labs(x="age", y="Cell percent (%)", title=NULL) +
#   theme(axis.text.x=element_text(size=8),
#         axis.text.y=element_text(size=8),
#         axis.title.x=element_text(size=9),
#         axis.title.y=element_text(size=9),
#         legend.position="top",
#         # legend.title=element_blank(),
#         title=element_text(size=10)) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#   geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
#   ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

### Calculate the freq at the rough, inter, detailed level, respectively; save the result
meta_info=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned")
# calculate the freq
meta_info_per.rough=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% count() %>% mutate(N_per_rough=n) %>% select(-n)
meta_info_per.inter=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter) %>% count() %>% mutate(N_per_inter=n) %>% select(-n)
meta_info_per.detailed=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% count() %>% mutate(N_per_detailed=n) %>% select(-n)
meta_info_all=meta_info %>% group_by(donor_id, sex, age) %>% count() %>% mutate(N=n) %>% select(-n)
meta_info_combined=meta_info_per.detailed %>%
  left_join(., meta_info_per.inter, by=intersect(colnames(.), colnames(meta_info_per.inter))) %>%
  left_join(., meta_info_per.rough, by=intersect(colnames(.), colnames(meta_info_per.rough))) %>%
  left_join(., meta_info_all, by=intersect(colnames(.), colnames(meta_info_all))) %>%
  mutate(percent.detailed=N_per_detailed/N,
         percent.inter=N_per_inter/N,
         percent.rough=N_per_rough/N) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  ungroup()

### Calculate the correlation with ages at the detailed level
# regardless of sex
meta_info_combined_detailed=meta_info_combined %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.detailed")
# females
meta_info_combined_detailed=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_F=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.detailed")
# males
meta_info_combined_detailed=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_M=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.detailed")
# combine all
COR_detailed=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Calculate the correlation with ages at the inter level
# regardless of sex
meta_info_combined_inter=meta_info_combined %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.inter")
# females
meta_info_combined_inter=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_F=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.inter")
# males
meta_info_combined_inter=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_M=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.inter")

# combine all
COR_inter=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Calculate the correlation with ages at the rough level
# regardless of sex
meta_info_combined_rough=meta_info_combined %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.rough")
# females
meta_info_combined_rough=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_F=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.rough")
# males
meta_info_combined_rough=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_M=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.rough")
# combine all
COR_rough=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Combine the results from all the levels
COR_combined=data.table::rbindlist(list(COR_rough, COR_inter, COR_detailed))
write.table(COR_combined, "~/Project_PBMCage/Results/PBMC_results/Cell_frequency_correlation.txt", sep="\t")
write.table(meta_info_combined, "~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")

#####################################



### Resave the celltype freq and their correlation with ages (Moidfy Monocytes annotation according to the published version)
#####################################
###

### Modified Monocyte annotation
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final_122subsets.rds")
object[[]]=object[[]] %>%
  mutate(Annot.detailed=ifelse(Annot.detailed=="Mono.inter.C1Q","Mono.CD64",Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed=="Mono.nonclassical.GNLY","Mono.GNLY_FCGR3A",Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed=="Mono.inter.GNLY","Mono.GNLY_CD14",Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed=="Mono.CD64","Mono.CD64",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed=="Mono.GNLY_FCGR3A","Mono.GNLY",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed=="Mono.GNLY_CD14","Mono.GNLY",Annot.inter))
object # make sure the seurat_obj has been normalize, HVGs found, and scaled
# cellcycle scoring
object=CellCycleScoring(object,
                        s.features=cc.genes.updated.2019$s.genes,
                        g2m.features=cc.genes.updated.2019$g2m.genes,
                        set.ident=TRUE)
# save the mt, rb, and cell-cycle info
meta_info=object@meta.data %>% select(any_of(c("nCount_RNA","nFeature_RNA","dataset","donor_id","sex","age",
                                               "percent.rb","percent.mt","Annot.detailed","Annot.inter","Annot.rough",
                                               "agecut",
                                               "S.Score","G2M.Score","Phase"))) %>%
  mutate(cellid=rownames(.))
meta_info$age=as.integer(meta_info$age)
data.table::fwrite(meta_info, "~/Project_PBMCage/Results/PBMC_results/Key_meta_info_updated.csv.gz", sep="\t")

meta_info=subset(meta_info, Annot.rough!="Unassigned")
# calculate the freq
meta_info_per.rough=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% count() %>% mutate(N_per_rough=n) %>% select(-n)
meta_info_per.inter=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter) %>% count() %>% mutate(N_per_inter=n) %>% select(-n)
meta_info_per.detailed=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% count() %>% mutate(N_per_detailed=n) %>% select(-n)
meta_info_all=meta_info %>% group_by(donor_id, sex, age) %>% count() %>% mutate(N=n) %>% select(-n)
meta_info_combined=meta_info_per.detailed %>%
  left_join(., meta_info_per.inter, by=intersect(colnames(.), colnames(meta_info_per.inter))) %>%
  left_join(., meta_info_per.rough, by=intersect(colnames(.), colnames(meta_info_per.rough))) %>%
  left_join(., meta_info_all, by=intersect(colnames(.), colnames(meta_info_all))) %>%
  mutate(percent.detailed=N_per_detailed/N,
         percent.inter=N_per_inter/N,
         percent.rough=N_per_rough/N) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  ungroup()

### Calculate the correlation with ages at the detailed level
# regardless of sex
meta_info_combined_detailed=meta_info_combined %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.detailed")
# females
meta_info_combined_detailed=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_F=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.detailed")
# males
meta_info_combined_detailed=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(., names_from="Annot.detailed", values_from="percent.detailed")
nonzero_=apply(meta_info_combined_detailed[! colnames(meta_info_combined_detailed) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_detailed=meta_info_combined_detailed[, filter]
celltypes_M=colnames(meta_info_combined_detailed)[! colnames(meta_info_combined_detailed) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_detailed[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_detailed))], 2, function(y) cor.test(meta_info_combined_detailed[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.detailed")
# combine all
COR_detailed=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Calculate the correlation with ages at the inter level
# regardless of sex
meta_info_combined_inter=meta_info_combined %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.inter")
# females
meta_info_combined_inter=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_F=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.inter")
# males
meta_info_combined_inter=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.inter, percent.inter) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.inter", values_from="percent.inter")
nonzero_=apply(meta_info_combined_inter[! colnames(meta_info_combined_inter) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_inter=meta_info_combined_inter[, filter]
celltypes_M=colnames(meta_info_combined_inter)[! colnames(meta_info_combined_inter) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_inter[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_inter))], 2, function(y) cor.test(meta_info_combined_inter[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.inter")

# combine all
COR_inter=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Calculate the correlation with ages at the rough level
# regardless of sex
meta_info_combined_rough=meta_info_combined %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
ALL_COR_DF=data.frame(celltypes=celltypes_,
                      cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                      analysis="all", celltype.level="Annot.rough")
# females
meta_info_combined_rough=meta_info_combined %>%
  subset(sex=="F") %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_F=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
F_COR_DF=data.frame(celltypes=celltypes_F,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="F", celltype.level="Annot.rough")
# males
meta_info_combined_rough=meta_info_combined %>%
  subset(sex=="M") %>%
  select(donor_id, sex, age, Annot.rough, percent.rough) %>%
  subset(!duplicated(.)) %>%
  tidyr::pivot_wider(., names_from="Annot.rough", values_from="percent.rough")
nonzero_=apply(meta_info_combined_rough[! colnames(meta_info_combined_rough) %in% c("age","sex","donor_id")], 2, function(x) length(x[!is.na(x)])>=3)
filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
meta_info_combined_rough=meta_info_combined_rough[, filter]
celltypes_M=colnames(meta_info_combined_rough)[! colnames(meta_info_combined_rough) %in% c("age","sex", "donor_id")]

pearson_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="pearson"))
all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
spearman_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="spearman", exact=F))
all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
kendall_results=apply(meta_info_combined_rough[,!grepl("^age$|^sex$|^donor_id$",colnames(meta_info_combined_rough))], 2, function(y) cor.test(meta_info_combined_rough[["age"]], y, method="kendall"))
all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
M_COR_DF=data.frame(celltypes=celltypes_M,
                    cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                    analysis="M", celltype.level="Annot.rough")
# combine all
COR_rough=data.table::rbindlist(list(ALL_COR_DF, F_COR_DF, M_COR_DF))

### Combine the results from all the levels
COR_combined=data.table::rbindlist(list(COR_rough, COR_inter, COR_detailed))
write.table(COR_combined, "~/Project_PBMCage/Results/PBMC_results/Cell_frequency_correlation_updated.txt", sep="\t")
write.table(meta_info_combined, "~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values_updated.txt", sep="\t")

#####################################



### Plot the celltype freq at the Annot.rough level
#####################################
###

### Plot the cellratio
object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

library(scRNAtoolVis)
plot_rough=
  cellRatioPlot(object=object,
                sample.name="agecut",
                celltype.name="Annot.rough",
                col.width=0.5,
                flow.alpha=0.75,
                flow.curve=0.25,
                fill.col=ggsci::pal_d3("category10")(10)) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=8))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_rough.pdf", width=6, height=3)
plot(plot_rough)
dev.off()

### Plot the overall CD4T and overall CD8T
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|rough",colnames(.))])) %>%
  subset(!duplicated(.))
plot_rough=
  ggplot(meta_info_combined, aes(x=age, y=percent.rough, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.rough, ncol=4, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpmisc::stat_correlation(vstep=0.1, show.legend=FALSE, size=3, method="spearman", small.p=TRUE, ggpmisc::use_label(c("R", "p")))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_rough.CD4_CD8_respective.pdf", width=7, height=4.5)
plot(plot_rough)
dev.off()

### Plot CD4/CD8 ratio
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined_CD4.CD8=meta_info_combined %>% subset(Annot.rough %in% c("CD4T cells", "CD8T cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|rough",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  select(-N_per_rough) %>%
  tidyr::pivot_wider(names_from="Annot.rough", values_from="percent.rough") %>%
  mutate(CD4_CD8.ratio=`CD4T cells`/`CD8T cells`)
# add the agecut metadata
AGECUT=list(c(19:30),c(31:47),c(48:68),c(69:97))
AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}
meta_info_combined_CD4.CD8$agecut=meta_info_combined_CD4.CD8$age
for (i in 1:length(AGECUT)) {
  meta_info_combined_CD4.CD8=meta_info_combined_CD4.CD8 %>%
    mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
}
# plot with agecut, but can also change to age
plot_CD4CD8ratio_agecut=
  ggplot(meta_info_combined_CD4.CD8, aes(x=agecut, y=CD4_CD8.ratio)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="CD4T/CD8T ratio", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_compare_means(method="anova",aes(label=paste0("p = ", after_stat(p.format))), size=3.5)

plot_CD4CD8ratio_agecut_perSex=
  ggplot(meta_info_combined_CD4.CD8, aes(x=agecut, y=CD4_CD8.ratio, color=sex)) +
  geom_boxplot() +
  theme_classic() +
  labs(x=NULL, y="CD4T/CD8T ratio", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  ggpubr::stat_compare_means(method="anova", aes(label=paste0("p = ", after_stat(p.format))), size=3.5)

plot_both=cowplot::plot_grid(plotlist=list(plot_CD4CD8ratio_agecut, plot_CD4CD8ratio_agecut_perSex), axis="tb", align="h", ncol=2, rel_widths=c(0.6,1))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_rough.CD4_CD8_ratio.pdf", width=8, height=2.5)
plot(plot_both)
dev.off()

#####################################



### Plot the celltype freq at the Annot.inter level
#####################################
###

### Plot the Annot.inter freq
# load the cell proportion (use the monocyte-modified version to be consistent with the annotation to be published)
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values_updated.txt", sep="\t")
meta_info_combined=meta_info_combined %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|inter",colnames(.))])) %>%
  subset(!duplicated(.))
# remove Other.Eryth
meta_info_combined=meta_info_combined %>% subset(Annot.inter!="Other.Eryth")
# remove the rare populations
meta_info_combined=meta_info_combined %>% subset(Annot.inter!="Other.Mast") # mean<0.001

plot_inter=
  ggplot(meta_info_combined, aes(x=age, y=percent.inter, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.inter, ncol=4, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.9,0.05),
        legend.text=element_text(size=9),
        # legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_inter.pdf", width=8.5, height=10)
plot(plot_inter)
dev.off()

# plot the variation in freq of proliferating cells
meta_info_combined_sub=meta_info_combined %>% subset(Annot.inter %in% c("CD4T.prolif","CD8T.CTL","NK.prolif")) %>%
  mutate(Annot.inter.new=ifelse(Annot.inter=="CD8T.CTL"|Annot.inter=="CD4T.prolif","T.prolif","NK.prolif")) %>%
  group_by(donor_id, sex, age, Annot.inter.new) %>%
  summarize_at("percent.inter",sum)
# plot_inter=
ggplot(meta_info_combined_sub, aes(x=age, y=percent.inter, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.inter.new, ncol=1, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="top",
        legend.text=element_text(size=9),
        # legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_line(data=. %>% group_by(age, Annot.inter.new, sex) %>% summarize_at("percent.inter",mean),
            aes(group=1))

#####################################



### Plot the celltype freq at the Annot.detailed level
#####################################
###

library(dplyr)
library(ggplot2)

### Plot the Annot.detailed freq
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.))
# according to the Annot.inte plot, we focus on CD4T, CD8T, NK
celltype_selected=list(CD4T=meta_info_combined$Annot.detailed[grepl("^CD4T\\.", meta_info_combined$Annot.detailed)],
                       CD8T=meta_info_combined$Annot.detailed[grepl("^CD8T\\.", meta_info_combined$Annot.detailed)],
                       NK=meta_info_combined$Annot.detailed[grepl("^NK\\.", meta_info_combined$Annot.detailed)])
meta_info_combined=meta_info_combined %>% subset(Annot.detailed %in% unlist(celltype_selected))

plot_detailed=
  ggplot(meta_info_combined, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, ncol=8, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="top",
        legend.text=element_text(size=9),
        # legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_detailed.pdf", width=22, height=8)
plot(plot_detailed)
dev.off()

### the remaining: B cells, DCs, Monocytes, Other cells, OtherT
# ... by comparing with "~/Project_PBMCage/Plots/Immunity/Immunity_Cell_frequency_detailed_theRest.pdf",
# ... no change can be found/confirmed in these celltypes
celltype_selected=list(B=meta_info_combined$Annot.detailed[grepl("^B\\.", meta_info_combined$Annot.detailed)],
                       DC=meta_info_combined$Annot.detailed[grepl("^DC\\.", meta_info_combined$Annot.detailed)],
                       Mono=meta_info_combined$Annot.detailed[grepl("^Mono\\.", meta_info_combined$Annot.detailed)],
                       Other_cells=meta_info_combined$Annot.detailed[grepl("^Other\\.", meta_info_combined$Annot.detailed)],
                       OtherT=meta_info_combined$Annot.detailed[grepl("^OtherT\\.", meta_info_combined$Annot.detailed)])
meta_info_combined=meta_info_combined %>% subset(Annot.detailed %in% unlist(celltype_selected))

plot_detailed=
  ggplot(meta_info_combined, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, ncol=8, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="top",
        legend.text=element_text(size=9),
        # legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_detailed_theRest.pdf", width=22, height=8)
plot(plot_detailed)
dev.off()

#####################################



### Function inference of the celltypes with freq changing with age, as proven by cross-validation in PBMCage and Immunity datasets
#####################################
###
# CD8T.Tem.GZMK_NKG7, CD8T.Tem.GZMKhi

library(Seurat)
library(clusterProfiler)
library(GseaVis)

PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Analyze CD4T.naive.COTL1pos_SOX4pos
# GSEA analysis
PBMCage_subset=subset(PBMCage_obj, Annot.inter=="CD4T.naive")
Idents(PBMCage_subset)="Annot.detailed"
CD4T_subpop_markers=FindMarkers(PBMCage_subset, ident.1="CD4T.naive.COTL1pos_SOX4pos")
CD4T_subpop_markers_sorted=CD4T_subpop_markers %>% arrange(desc(avg_log2FC))

gene_list=CD4T_subpop_markers_sorted$avg_log2FC
names(gene_list)=rownames(CD4T_subpop_markers_sorted)

ego=gseGO(gene_list, ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
# ego2=gseGO(gene_list, ont="CC", OrgDb="org.Hs.eg.db", keyType="SYMBOL") # No result
# ego3=gseGO(gene_list, ont="MF", OrgDb="org.Hs.eg.db", keyType="SYMBOL") # No result
dotplot(ego)
# dotplot(ego2)
# dotplot(ego3)
View(ego@result) # check the result

# plot the meaningful enriched terms
gene_of_interest=ego@result %>% subset(ID %in% c("GO:0019221","GO:0031347")) %>% select(core_enrichment)
gene_of_interest=gene_of_interest %>% lapply(., function(x) strsplit(x, "/"))
gene_of_interest_top=Reduce(intersect, gene_of_interest) %>% lapply(., function(x) x[1:3])
gene_of_interest_top=unlist(gene_of_interest_top)

source("~/Rscripts/gseaNb2.R")
gse_CD4=
  gseaNb2(object=ego,
         geneSetID=c("GO:0019221","GO:0031347"),
         base_size=11,
         segCol="grey",
         lineSize=0.2,
         subPlot=2,
         termWidth=35,
         legend.position=c(0.75,0.76),
         addGene=gene_of_interest_top, geneSize=3,
         addPval=T,
         pvalX=0.05, pvalY=0.05)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/FunctionInference_CD4TnaiveSOX4pos_gsea.pdf", width=6, height=4)
plot(gse_CD4)
dev.off()

# plot the correlation between CD4T.naive.COTL1pos_SOX4pos and Th2 cells (CD4T.Tcm), as SOX4 is reported to inhibit Th2 diff
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% subset(Annot.detailed %in% c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm"))
meta_info_combined=meta_info_combined %>% select(donor_id, sex, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(names_from="Annot.detailed", values_from="percent.detailed")

plot_=
  ggplot(meta_info_combined, aes(x=CD4T.naive.COTL1pos_SOX4pos, y=CD4T.Tcm, color=sex)) +
  geom_point(size=0.5, shape=20, alpha=0.5) +
  theme_classic() +
  labs(x="Percent of \nCD4T.naive.COTL1pos_SOX4pos (%)", y="Percent of CD4T.Tcm (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="top",
        legend.text=element_text(size=9),
        plot.subtitle=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot") +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/FunctionInference_CD4TnaiveSOX4pos_cor.W.Th2cells.pdf", height=4, width=3)
plot(plot_)
dev.off()

# confirmed by Terekhova et al.'s dataset
meta_info_combined2=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined2=meta_info_combined2 %>% subset(Annot.detailed %in% c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm"))
tempt=meta_info_combined2 %>% select(donor_id, sex, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(names_from="Annot.detailed", values_from="percent.detailed")
ggplot(tempt, aes(x=CD4T.naive.COTL1pos_SOX4pos, y=CD4T.Tcm, color=sex)) +
  geom_point() +
  theme_classic() +
  labs(x="Percent of CD4T.naive.COTL1pos_SOX4pos (%)", y="Percent of CD4T.Tcm (%)", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position="top",
        legend.text=element_text(size=9),
        plot.subtitle=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot") +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=1), show.legend=FALSE, size=3.5, method="spearman")


### Analyze NK.CD56dim.FCER1Glo_KLRC2pos
# heatmap of the core definitive genes
PBMCage_subset=subset(PBMCage_obj, Annot.inter=="NK.CD56dim")
Idents(PBMCage_subset)="Annot.detailed"
core_gene_df=FetchData(PBMCage_subset, vars=c("donor_id","sex","age",
                                              "KLRC2","KLRC3","GZMH","KLRB1","KLRC1","KLRF1","FCER1G",
                                              "Annot.detailed"), layer="data")
core_gene_df_rearranged=core_gene_df %>%
  group_by(Annot.detailed) %>%
  summarize_at(c("KLRC2","KLRC3","GZMH","KLRB1","KLRC1","KLRF1","FCER1G"),mean) %>%
  ungroup() %>% tibble::column_to_rownames("Annot.detailed")
core_gene_df_rearranged=t(scale(core_gene_df_rearranged))

library(ComplexHeatmap)
col_fun=circlize::colorRamp2(c(-2, 0, 2), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(core_gene_df_rearranged.t,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=T, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(matrix_)*unit(4.5, "mm"),
          height=nrow(matrix_)*unit(4.5, "mm"),
  )
lgd=Legend(col_fun=col_fun,
           title="Z-score",
           legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_final=draw(ht_, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/FunctionInference_NK.CD56dim.FCER1Glo_KLRC2pos_coreGenesHeatmap.pdf", height=4, width=2.5)
ht_final
dev.off()

# gsea comparison with other NK.CD56dim subpopulations
PBMCage_subset=subset(PBMCage_obj, Annot.inter=="NK.CD56dim")
Idents(PBMCage_subset)="Annot.detailed"
NKCD56dim_subpop_markers=FindAllMarkers(PBMCage_subset)
NKCD56dim_subpop_markers_sorted=NKCD56dim_subpop_markers %>% split(.$cluster) %>%
  lapply(., function(df) df %>% arrange(desc(avg_log2FC)))
gene_list=NKCD56dim_subpop_markers_sorted %>% lapply(., function(df) df %>% select(gene, avg_log2FC) %>% tibble::deframe())
gse_results=compareCluster(gene_list, fun="gseGO", ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")

module_=gse_results@compareClusterResult %>%
  tibble::rowid_to_column("index") %>%
  select(index, Cluster, Description, NES, p.adjust) %>%
  mutate(orientation=ifelse(NES>0,"highly expr.","lowly expr."))
module_$Cluster=as.character(module_$Cluster)
module_=module_ %>%
  mutate(Cluster=ifelse(Cluster %in% c("NK.CD56dim.ISG","NK.CD56dim.FCER1G","NK.CD56dim.FOS","NK.CD56dim.HNRNPH1"), "NK.CD56dim.FCER1G+", Cluster)) %>%
  group_by(Cluster, orientation) %>%
  slice_max(index, n=3) %>%
  ungroup()

# module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
module_$Description=forcats::fct_relevel(module_$Description, unique(module_$Description))
module_$Cluster=forcats::fct_relevel(module_$Cluster, c("NK.CD56dim.FCER1G+",
                                                        "NK.CD56dim.FCER1Glo_KLRC2pos","NK.CD56dim.FCER1Glo_KLRC2neg"))
module_$Cluster=gsub("NK\\.CD56dim\\.","",module_$Cluster)
module_$Cluster=forcats::fct_relevel(module_$Cluster, gsub("NK\\.CD56dim\\.","",levels(module_$Cluster)))

plot_gse=
  ggplot(module_, aes(x=Cluster, y=Description, color=orientation, size=abs(NES))) +
  facet_wrap(~orientation, ncol=2, scales="free_x") +
  geom_point() +
  scale_size_continuous(range=c(1,4)) +
  scale_color_manual(values=c("brown3", "dodgerblue3")) +
  guides(size=guide_legend(order=1)) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10, lineheight=0.75),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.key.width=unit(dev.size()[1] / 50, "inches"),
        legend.direction="vertical",
        legend.box="vertical",
        # legend.box.spacing=margin(0.25),
        legend.position="right") +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  guides(color=FALSE,
         size=guide_legend(title="abs.NES", title.theme=element_text(size=9)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/FunctionInference_NK.CD56dim.FCER1Glo_KLRC2pos_gseaCompare.pdf", width=12, height=4)
plot(plot_gse)
dev.off()

#####################################



### Causual analysis between freq of CD4T.naive.COTL1pos_SOX4pos and Th2 CD4 cells
#####################################
###

if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
SOX4cell_on_CD4Tcm=meta_info_combined %>%
  subset(Annot.detailed %in% c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm")) %>%
  select(donor_id, sex, age, Annot.detailed, percent.detailed) %>%
  tidyr::pivot_wider(names_from="Annot.detailed", values_from="percent.detailed")

# # remove outliers
# outliers1=boxplot.stats(SOX4cell_on_CD4Tcm$CD4T.naive.COTL1pos_SOX4pos)$out
# outliers2=boxplot.stats(SOX4cell_on_CD4Tcm$CD4T.Tcm)$out
# SOX4cell_on_CD4Tcm=SOX4cell_on_CD4Tcm %>%
#   subset(!(CD4T.naive.COTL1pos_SOX4pos %in% outliers1)) %>%
#   subset(!(CD4T.Tcm %in% outliers2))

# set threshold
# ... based on normal distribution of the cell frequencies (cutoff=80%, i.e., the ones>=80% of the frequencies are considered highly present)
# ggplot(SOX4cell_on_CD4Tcm, aes(x=CD4T.naive.COTL1pos_SOX4pos)) +
#   geom_density()
# ggplot(SOX4cell_on_CD4Tcm, aes(x=CD4T.Tcm)) +
#   geom_density()
CD4T.naive.COTL1pos_SOX4pos_par=MASS::fitdistr(SOX4cell_on_CD4Tcm$CD4T.naive.COTL1pos_SOX4pos %>% .[!is.na(.)], densfun="normal")
CD4T.Tcm_par=MASS::fitdistr(SOX4cell_on_CD4Tcm$CD4T.Tcm %>% .[!is.na(.)], densfun="normal")
thr1=qnorm(0.68, mean=CD4T.naive.COTL1pos_SOX4pos_par$estimate[1], sd=CD4T.naive.COTL1pos_SOX4pos_par$estimate[2])
thr2=qnorm(0.68, mean=CD4T.Tcm_par$estimate[1], sd=CD4T.Tcm_par$estimate[2])
SOX4cell_on_CD4Tcm=SOX4cell_on_CD4Tcm %>%
  mutate(SOX4pos_or_neg=ifelse(CD4T.naive.COTL1pos_SOX4pos<=thr1 | is.na(CD4T.naive.COTL1pos_SOX4pos), "neg", "pos")) %>%
  mutate(Th2pos_or_neg=ifelse(CD4T.Tcm<=thr2 | is.na(CD4T.Tcm), "neg", "pos"))
# check the assumption1: no effect of the intervention underlying the control series
ggplot(SOX4cell_on_CD4Tcm %>% subset(SOX4pos_or_neg=="pos"), aes(x=CD4T.naive.COTL1pos_SOX4pos, y=CD4T.Tcm)) +
  geom_point() +
  geom_smooth(show.legend=FALSE, method="loess") +
  ggpubr::stat_cor(show.legend=FALSE, method="spearman")
ggplot(SOX4cell_on_CD4Tcm %>% subset(Th2pos_or_neg=="neg"), aes(x=CD4T.naive.COTL1pos_SOX4pos, y=CD4T.Tcm)) +
  geom_point() +
  geom_smooth(show.legend=FALSE, method="loess") +
  ggpubr::stat_cor(show.legend=FALSE, method="spearman")
# check the assumption2: the relation between y and covariate (CD4T.Tcm ~ age) is stable with/without the intervention (SOX4pos_or_neg)
ggplot(SOX4cell_on_CD4Tcm, aes(x=age, y=CD4T.Tcm, color=SOX4pos_or_neg)) +
  geom_point() +
  geom_smooth(aes(group=SOX4pos_or_neg), show.legend=FALSE, method="loess") +
  ggpubr::stat_cor(aes(group=SOX4pos_or_neg), show.legend=FALSE, method="spearman")
# check the assumption2: the relation between y and covariate (CD4T.naive.COTL1pos_SOX4pos ~ age) is stable with/without the intervention (Th2pos_or_neg)
ggplot(SOX4cell_on_CD4Tcm, aes(x=age, y=CD4T.naive.COTL1pos_SOX4pos, color=Th2pos_or_neg)) +
  geom_point() +
  geom_smooth(aes(group=Th2pos_or_neg), show.legend=FALSE, method="loess") +
  ggpubr::stat_cor(aes(group=Th2pos_or_neg), show.legend=FALSE, method="spearman")
# check the assumpiton3: prior.level.sd (the variability in the smoothing trend)
prior.level.sd_CD4T.Tcm=sd(predict(loess(CD4T.Tcm~age,SOX4cell_on_CD4Tcm)))
prior.level.sd_CD4T.Tcm<0.01
prior.level.sd_CD4T.naive.COTL1pos_SOX4pos=sd(predict(loess(CD4T.naive.COTL1pos_SOX4pos~age,SOX4cell_on_CD4Tcm)))
prior.level.sd_CD4T.naive.COTL1pos_SOX4pos<0.01

# check if CD4T.naive.COTL1pos_SOX4pos is the cause
mat_pre=SOX4cell_on_CD4Tcm %>% subset(SOX4pos_or_neg=="neg") %>% select(CD4T.Tcm, age) %>% as.matrix()
mat_post=SOX4cell_on_CD4Tcm %>% subset(SOX4pos_or_neg=="pos") %>% select(CD4T.Tcm, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
SOX4_res=list() # since the modeling faces some randomality, we here repeated the process 100 times
for (i in 1:100) {
  impact_SOX4=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
  SOX4_res[[i]]=impact_SOX4$summary
}
table(Reduce(c, lapply(SOX4_res, function(x) x$p))<0.05) # turns out to be >99% for p<0.05
plot_impact_SOX4=plot(impact_SOX4, "cumulative") # plot the last time as the representative result

# check if CD4T.Tcm is the cause
mat_pre=SOX4cell_on_CD4Tcm %>% subset(Th2pos_or_neg=="neg") %>% select(CD4T.naive.COTL1pos_SOX4pos, age) %>% as.matrix()
mat_post=SOX4cell_on_CD4Tcm %>% subset(Th2pos_or_neg=="pos" & !is.na(CD4T.naive.COTL1pos_SOX4pos)) %>% select(CD4T.naive.COTL1pos_SOX4pos, age) %>% as.matrix()
pre.period=c(1, nrow(mat_pre))
post.period=c(nrow(mat_pre)+1, nrow(mat_pre)+nrow(mat_post))
Th2_res=list() # since the modeling faces some randomality, we here repeated the process 100 times
for (i in 1:100) {
  impact_Th2=CausalImpact(rbind(mat_pre, mat_post), pre.period, post.period)
  Th2_res[[i]]=impact_Th2$summary
}
table(Reduce(c, lapply(Th2_res, function(x) x$p))<0.05) # turns out to be >99% for p>0.05
plot_impact_Th2=plot(impact_Th2, "cumulative") # plot the last time as the representative result

# merge the results for plotting
merged_df=data.table::rbindlist(list(plot_impact_SOX4$data %>% mutate(intervention="CD4T.naive.COTL1pos_SOX4pos"),
                                     plot_impact_Th2$data %>% mutate(intervention="CD4T.Tcm")))
summary_SOX4=data.frame(AbsEffect=sprintf("%.2e",impact_SOX4$summary$AbsEffect[1]),
                        RelEffect=sprintf("%.2f%%",impact_SOX4$summary$RelEffect[1]),
                        AbsEffect.lower=sprintf("%.2e",impact_SOX4$summary$AbsEffect.lower[1]),
                        AbsEffect.upper=sprintf("%.2e",impact_SOX4$summary$AbsEffect.upper[1]),
                        RelEffect.lower=sprintf("%.2f%%",impact_SOX4$summary$RelEffect.lower[1]),
                        RelEffect.upper=sprintf("%.2f%%",impact_SOX4$summary$RelEffect.upper[1]),
                        pval=sprintf("%.3f",impact_SOX4$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="CD4T.naive.COTL1pos_SOX4pos")

summary_Th2=data.frame(AbsEffect=sprintf("%.2e",impact_Th2$summary$AbsEffect[1]),
                        RelEffect=sprintf("%.2f%%",impact_Th2$summary$RelEffect[1]),
                        AbsEffect.lower=sprintf("%.2e",impact_Th2$summary$AbsEffect.lower[1]),
                        AbsEffect.upper=sprintf("%.2e",impact_Th2$summary$AbsEffect.upper[1]),
                        RelEffect.lower=sprintf("%.2f%%",impact_Th2$summary$RelEffect.lower[1]),
                        RelEffect.upper=sprintf("%.2f%%",impact_Th2$summary$RelEffect.upper[1]),
                        pval=sprintf("%.3f",impact_Th2$summary$p[1])) %>%
  mutate(`Effect`=paste0(AbsEffect," [", RelEffect, "]"),
         `95% CI`=paste0(AbsEffect.lower," ~ ",AbsEffect.upper, " [", RelEffect.lower," ~ ",RelEffect.upper, "]")) %>%
  select(Effect, `95% CI`, pval) %>%
  mutate(Intervention="CD4T.Tcm")

summary_df_merged=rbind(summary_SOX4, summary_Th2)

# plot
library(gridExtra)
library(ggplot2)
plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free") +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="donor #", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        plot.margin=margin(10,10,65,10)
  ) +
  annotation_custom(tableGrob(summary_df_merged %>% tibble::column_to_rownames("Intervention"),
                              theme=ttheme_minimal(base_size=10)),
                    xmin=1100, ymax=-5.5)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CausualInference_CD4TnaiveSOX4pos_cor.W.Th2cells.pdf", width=7, height=4)
plot_infer
dev.off()

#####################################



### Phenotype scoring
#####################################
###

### Loading
library(Seurat)
library(UCell)
library(dplyr)

object=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

### Score CD8T and NKT
CD8T.NKT_detailedcelltypes=names(table(object$Annot.detailed))[grepl("^CD8T\\.|^OtherT\\.NKT", names(table(object$Annot.detailed)))]
CD8T.NKT_subset=subset(object, Annot.detailed %in% CD8T.NKT_detailedcelltypes)
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
CD8T_Ucell=AddModuleScore_UCell(CD8T.NKT_subset, features=T_cell_scoring)
# extract the ucell scores
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
score_df=cbind(score_df, CD8T_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_CD8TandNKT.rds")

### Score CD4T
CD4T_detailedcelltypes=names(table(object$Annot.detailed))[grepl("^CD4T\\.", names(table(object$Annot.detailed)))]
CD4T_subset=subset(object, Annot.detailed %in% CD4T_detailedcelltypes)
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
CD4T_Ucell=AddModuleScore_UCell(CD4T_subset, features=T_cell_scoring)
# extract the ucell scores
score_df=CD4T_Ucell[[]][, colnames(CD4T_Ucell[[]])[grepl("_UCell$", colnames(CD4T_Ucell[[]]))]]
score_df=cbind(score_df, CD4T_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_CD4T.rds")

### Score B cells
Bcell_subset=subset(object, Annot.rough=="B cells")
B_cell_scoring=list(naive=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"),
                    `isotype switching`=c("AICDA","ATAD5","BATF","CD40LG","ERCC1","EXO1","EXOSC3","EXOSC6","ICOSLG","LIG4","MLH1","MSH2","MSH6","NBN","NFKBIZ","RNF8","RNF168","SANBR","SWAP70","XRCC4","UNG"),
                    activation=c("ABL1","CD19","CD180","GAPT","PLCL2","TLR4","IRF8","PLCG2","SPI1","ADA","CDH17","IL6","IL21","ITFG2","PHF14","POU2AF1","BCL3","CDH17","DLL1","DOCK10","DOCK11","LFNG","MFNG","PTK2B","NOTCH2","TNFAIP3","BCL6","ENPP1","IL2","IL10","ITM2A","LGALS1","LGALS8","NFKBIZ","NKX2-3","XBP1","ST3GAL1"),
                    migration=c("CH25H","CXCL13","HSD3B7","CYP7B1","PTK2B","GAS6","XCL1"),
                    maturation=c("ABL1","ATM","ATP11C","CD24","FNIP1","FOXP1","IGHM","IRF2BP2","KIT","LRRC8A","PRKDC","RAG1","RAG2","SPI1","SPIB","SYVN1","TNFSF13B","TRAF3IP2"),
                    `Ig production`=c("CRLF2","ENPP1","FGL2","GAPT","NFKBIZ","NOD2","POU2F2","TREX1"),
                    inhibitory=c("BCL6","LILRB4","ATM","BTK","BTLA","CASP3","CD24","CD300A","CDKN2A","CTLA4","FCGR2B","IL10","INPP5D","LYN","PAWR","PKN1","PRDM1","TNFRSF13B","PTEN","RC3H1","TNFRSF21","TSC2","TYROBP","NDFIP1","FOXP3","PARP3"),
                    `Ag presentation`=c("HLA-DOB","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DMA","HLA-DOA","HLA-DPA1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","FCER2","IGHE","HLA-E","HLA-G","HLA-F","CD70"),
                    proliferation=c("BAX","RAG2","ABL1","CD19","CD180","GAPT","PLCL2","TLR4"),
                    transcription=c("BCL2","PRKDC","RAG2","TCF3","TRP53"),
                    costimulation=c("CD320","IL4","TNFSF13B","TNFSF13C"))
Bcell_Ucell=AddModuleScore_UCell(Bcell_subset, features=B_cell_scoring)
# extract the ucell scores
score_df=Bcell_Ucell[[]][, colnames(Bcell_Ucell[[]])[grepl("_UCell$", colnames(Bcell_Ucell[[]]))]]
score_df=cbind(score_df, Bcell_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_Bcells.rds")

### Score monocytes
Myeloid_subset=subset(object, Annot.rough=="Monocytes")
Myeloid_scoring=list(
  activation=c("DYSF","FER1L5","BTK"),
  migration=c("ANXA1","CCL12","CCL26","CCR2","CTSG","FLT1","LGALS3","MSMP","MYO9B","PDGFB","PTPRO","RPS19","TNFSF11"),
  differentiation=c("BMP4","CSF1","CSF2","FASN","GMPR2","GPR68","IL3RA","IL3","IL31RA","JUN","MED1","MEF2C","MNDA","MYH9","PDE1B","PIR","PPARG","SP3","THOC5","VEGFA"),
  aggregation=c("BMP7","CD44","CD47","CFH","IL1B","THBS1"),
  adhesion=c("ADD2","CCR2","CHST4","CX3CR1","EXT1","GCNT1","GOLPH3","ITGA4","ITGAM","ITGB1","ITGB7","JAM2","LEP","LRG1","MADCAM1","PODXL2","ROCK1","SELE","SELL","SELP","SELPLG","SLC39A8","SPN","TNF","VCAM1"),
  transcription=c("MAFB","MAF","EGR1","IRF8"))
Mono_Ucell=AddModuleScore_UCell(Myeloid_subset, features=Myeloid_scoring)
# extract the ucell scores
score_df=Mono_Ucell[[]][, colnames(Mono_Ucell[[]])[grepl("_UCell$", colnames(Mono_Ucell[[]]))]]
score_df=cbind(score_df, Mono_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_Monocytes.rds")

### Score DCs
DC_subset=subset(object, Annot.rough=="DCs")
DC_scoring=list(maturation=c("CCL19","CCR7","F2RL1","PRTN3","BATF","BATF2","BATF3","CAMK4","CSF2","DCSTAMP","GIMAP3","GIMAP5","IL4","IRF4","LTBR","NOTCH2","PIRB","PSEN1","RELB","SPI1","TNFSF9","TRAF6","UBD","IRF8","SLAMF9","ZBTB46"),
                activation=c("DOCK2","PYCARD","SLAMF1"),
                migration=c("ALOX5","ANO6","ASB2","CCL19","CCL21","CDC42","DOCK8","EPS8","EXT1","GPR183","NLRP12","TRPM2","TRPM4"),
                `cytokine prod.`=c("CLEC7A","DDX1","DDX21","DHX36","KIT","MAVS","NOD2","PLCG2","RIGI","SCIMP","SLAMF9","TICAM1"),
                `Ag presentation`=c("HLA-DOB","HLA-DRA","HLA-DMB","HLA-DPB1","HLA-DMA","HLA-DOA","HLA-DPA1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5","CLEC4A"),
                inhibitory=c("THBS1","CD68","FGL2","JAK3","BST2","HAVCR2","CD37","IL10","TSPAN32","FCGR2B"),
                proliferation=c("TBK1","AZI2","BCL2L1"),
                transcription=c("BATF3","IRF8","ID2","IRF4","NFIL3","SPI1","TCF4"))
DC_Ucell=AddModuleScore_UCell(DC_subset, features=DC_scoring)
# extract the ucell scores
score_df=DC_Ucell[[]][, colnames(DC_Ucell[[]])[grepl("_UCell$", colnames(DC_Ucell[[]]))]]
score_df=cbind(score_df, DC_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_DCs.rds")

### Score NK
NK_detailedcelltypes=names(table(object$Annot.detailed))[grepl("NK\\.", names(table(object$Annot.detailed)))]
NK_subset=subset(object, Annot.detailed %in% NK_detailedcelltypes)
NK_cell_scoring=list(maturation=c("ITGAL","CXCR1","CX3CR1","CMKLR1","FCGR3A"),
                     adhesion=c("ITGAL","CD2","CD58","ITGAM","ICAM1","NCAM1","B3GAT1","SELL"),
                     cytokine=c("IL2RA","IL2RB","IL2RG","IL7R","IL18R1","KIT","CCR5","CCR7","CXCR1","CXCR3","CXCR4","CX3CR1","CMKLR1"),
                     activation=c("CD244","CD59","CD226","KLRF1","SLAMF6","SLAMF7","FCGR3A","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DL4","KIR2DS4","KIR2DS5","KIR3DS1","KLRC2","KLRK1","NCR3","NCR2","NCR1","CD96"),
                     cytotoxic=c("FAS","FASLG","CD40LG","TNFSF10","LAMP1","LAMP2","GNLY","NKG7","GZMB","PRF1"),
                     inhibition=c("KLRC1","KLRD1","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5","KIR3DL1","KIR3DL2","LILRB1","PDCD1","SIGLEC7","CD300A","TIGIT"))
NK_Ucell=AddModuleScore_UCell(NK_subset, features=NK_cell_scoring)
# extract the ucell scores
score_df=NK_Ucell[[]][, colnames(NK_Ucell[[]])[grepl("_UCell$", colnames(NK_Ucell[[]]))]]
score_df=cbind(score_df, NK_Ucell[[]][, c("age","sex","agecut","Annot.detailed","Annot.inter","Annot.rough")])
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/FunctionScore_NK.rds")

#####################################



### Analyze the B cell scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_Bcells.rds")
# combine extremely minor celltypes at the detailed level
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.inter.lambda","B.inter.lambda.IL4R"), "B.inter.IL4R", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("B.mem.IGHA","B.mem.IGHD","B.mem.IGHG"), "B.mem.IGHA", Annot.detailed))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_B=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
B_final=draw(ht_B, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_Bcells.pdf", height=5, width=4.5)
draw(B_final)
dev.off()

#####################################



### Analyze the CD4T scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_CD4T.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD4T cells")
# remove some celltypes that do not exist in the Immunity dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD4T.prolif","CD4T.prolif.TSHZ2","CD4T.Tcm.ISG","CD4T.Tcm.ISG.TSHZ2")))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_CD4=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
CD4T_final=draw(ht_CD4, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_CD4T.pdf", height=4.5, width=5)
draw(CD4T_final)
dev.off()

#####################################



### Analyze the CD8T scoring results and replot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_CD8TandNKT.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD8T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$age)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD8T.Tem.ISG")))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_CD8=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
CD8T_final=draw(ht_CD8, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_CD8T.pdf", height=5, width=5)
draw(CD8T_final)
dev.off()

#####################################



### Analyze the NK scoring results and replot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_NK.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="NK cells")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_NK=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
NK_final=draw(ht_NK, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_NK.pdf", height=4.5, width=4.5)
draw(NK_final)
dev.off()

#####################################



### Analyze the Monocytes scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_Monocytes.rds")
CD8T.NKT_df=CD8T.NKT_df %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.C1Q"), "Mono.CD64", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.nonclassical.GNLY"), "Mono.GNLY_FCGR3A", Annot.detailed)) %>%
  mutate(Annot.detailed=ifelse(Annot.detailed %in% c("Mono.inter.GNLY"), "Mono.GNLY_CD14", Annot.detailed))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_Mono=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
Mono_final=draw(ht_Mono, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_Monocytes.pdf", height=3, width=3.5)
draw(Mono_final)
dev.off()

#####################################



### Analyze the DC scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_DCs.rds")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean)
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,4:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed"))

# check which UCell_group is changing in which celltype
correlation_results %>% subset(rho_pval<0.05 & abs(rho)>0.4) %>% select(Ucell_group) %>% subset(!duplicated(.))

### Make the df
correlation_results_subset=correlation_results %>% select(Ucell_group, Annot.detailed, rho, rho_pval) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group))
all_cell_cor=correlation_results_subset %>% select(-rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results_subset %>% select(-rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
ht_DC=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=1),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(7, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y,
                      gp=gpar(fontsize=10, fontface="bold"))}
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
DC_final=draw(ht_DC, annotation_legend_list=lgd, annotation_legend_side="top")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_DCs.pdf", height=3.5, width=4)
draw(DC_final)
dev.off()

#####################################



### Analyze the CD4T scoring results in females vs. males seperately
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_CD4T.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD4T cells")
# remove some celltypes that do not exist in the Immunity dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD4T.prolif","CD4T.prolif.TSHZ2","CD4T.Tcm.ISG","CD4T.Tcm.ISG.TSHZ2")))

### Correlation analysis on females
df=CD8T.NKT_df %>% group_by(sex, age, agecut, Annot.detailed) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean) %>%
  subset(sex=="F")
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results_F=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed")) %>%
  mutate(sex="F")

### Correlation analysis on males
df=CD8T.NKT_df %>% group_by(sex, age, agecut, Annot.detailed) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean) %>%
  subset(sex=="M")
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results_M=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed")) %>%
  mutate(sex="M")

# merge results of females and males
correlation_results=data.table::rbindlist(list(correlation_results_F, correlation_results_M)) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group)) %>%
  mutate(Ucell_group_sex=paste0(Ucell_group,",",sex))

### Make the df
all_cell_cor=correlation_results %>% select(Ucell_group_sex, Annot.detailed, rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group_sex", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
colnames_=gsub(",[A-Z]","",colnames(all_cell_cor)); colnames_=colnames_[!duplicated(colnames_)]
all_cell_cor=all_cell_cor %>%
  dplyr::relocate(any_of(paste0(rep(colnames_, each=2),",",c("F","M"))))
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results %>% select(Ucell_group_sex, Annot.detailed, rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group_sex", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr=all_cell_fdr %>%
  dplyr::relocate(any_of(paste0(rep(colnames_, each=2),",",c("F","M"))))
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
split=rep(1:length(colnames_), each=2)
features_=HeatmapAnnotation(
  anno=anno_block(gp=gpar(fill="transparent", col="white"), labels=colnames_, labels_rot=90, labels_just=1,
                  labels_offset=1, labels_gp=gpar(fontsize=10.5))
)
sex_anno=HeatmapAnnotation(
  `sex: F/M`=rep(c(1,2), ncol(all_cell_cor)/2), 
  col=list(`sex: F/M`=c("1"=scales::hue_pal()(2)[1], "2"=scales::hue_pal()(2)[2])),
  simple_anno_size=unit(1, "mm"), show_legend=F, annotation_name_gp=gpar(fontsize=10, fontfamily="mono")
)
ht_CD4=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(3.5, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y, rot=90, just="top",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          column_split=split, column_gap=unit(1, "mm"), column_title=NULL,
          bottom_annotation=features_,
          top_annotation=sex_anno
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
CD4T_final=draw(ht_CD4, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_CD4T_Fvs.M.pdf", height=4.5, width=6)
draw(CD4T_final)
dev.off()

#####################################



### Analyze the CD8T scoring results in females vs. males seperately
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Tempt_RDS/FunctionScore_CD8TandNKT.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD8T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD8T.Tem.ISG")))

### Correlation analysis on females
df=CD8T.NKT_df %>% group_by(sex, age, agecut, Annot.detailed) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean) %>%
  subset(sex=="F")
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results_F=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed")) %>%
  mutate(sex="F")

### Correlation analysis on males
df=CD8T.NKT_df %>% group_by(sex, age, agecut, Annot.detailed) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  summarize_at(colnames(.)[!grepl("sex|age|Annot\\.",colnames(.))], mean) %>%
  subset(sex=="M")
df$age=as.integer(as.character(df$age))
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$estimate))
cor=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="pearson"))$p.value))
cor_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="cor_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$estimate))
rho=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="spearman", exact=F))$p.value))
rho_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="rho_pval")

df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$estimate))
tau=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau")
df_celltype=df %>% split(~Annot.detailed) %>% lapply(., function(df) apply(df[,5:ncol(df)], 2, function(y) (cor.test(df[["age"]], y, method="kendall"))$p.value))
tau_pval=as.data.frame(df_celltype) %>% mutate(Ucell_group=rownames(.)) %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("Ucell",colnames(.))], names_to="Annot.detailed", values_to="tau_pval")

correlation_results_M=cor %>%
  left_join(cor_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(rho_pval, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau, by=c("Ucell_group","Annot.detailed")) %>%
  left_join(tau_pval, by=c("Ucell_group","Annot.detailed")) %>%
  mutate(sex="M")

# merge results of females and males
correlation_results=data.table::rbindlist(list(correlation_results_F, correlation_results_M)) %>%
  mutate(Ucell_group=gsub("_UCell","",Ucell_group)) %>%
  mutate(Ucell_group_sex=paste0(Ucell_group,",",sex))

### Make the df
all_cell_cor=correlation_results %>% select(Ucell_group_sex, Annot.detailed, rho) %>%
  tidyr::pivot_wider(names_from="Ucell_group_sex", values_from="rho") %>%
  tibble::column_to_rownames("Annot.detailed")
colnames_=gsub(",[A-Z]","",colnames(all_cell_cor)); colnames_=colnames_[!duplicated(colnames_)]
all_cell_cor=all_cell_cor %>%
  dplyr::relocate(any_of(paste0(rep(colnames_, each=2),",",c("F","M"))))
all_cell_cor[is.na(all_cell_cor)]=0

all_cell_fdr=correlation_results %>% select(Ucell_group_sex, Annot.detailed, rho_pval) %>%
  tidyr::pivot_wider(names_from="Ucell_group_sex", values_from="rho_pval") %>%
  tibble::column_to_rownames("Annot.detailed")
all_cell_fdr=all_cell_fdr %>%
  dplyr::relocate(any_of(paste0(rep(colnames_, each=2),",",c("F","M"))))
all_cell_fdr[is.na(all_cell_fdr)]=1
all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001,"****",
                                                            ifelse(x>=0.0001 & x<0.001, "***",
                                                                   ifelse(x>=0.001 & x<0.01, "**",
                                                                          ifelse(x>=0.01 & x<0.05, "*", "")))))

### Plot the correlation matrix
col_fun=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
split=rep(1:length(colnames_), each=2)
features_=HeatmapAnnotation(
  anno=anno_block(gp=gpar(fill="transparent", col="white"), labels=colnames_, labels_rot=90, labels_just=1,
                  labels_offset=1, labels_gp=gpar(fontsize=10.5))
)
sex_anno=HeatmapAnnotation(
  `sex: F/M`=rep(c(1,2), ncol(all_cell_cor)/2), 
  col=list(`sex: F/M`=c("1"=scales::hue_pal()(2)[1], "2"=scales::hue_pal()(2)[2])),
  simple_anno_size=unit(1, "mm"), show_legend=F, annotation_name_gp=gpar(fontsize=10, fontfamily="mono")
)
ht_CD4=
  Heatmap(all_cell_cor,
          name="mat", show_heatmap_legend=F,
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0),
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          column_names_gp=gpar(fontsize=10.5), row_names_gp=gpar(fontsize=10),
          width=ncol(all_cell_cor)*unit(3.5, "mm"),
          height=nrow(all_cell_cor)*unit(6, "mm"),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(all_cell_fdr), x, y, rot=90, just="top",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          column_split=split, column_gap=unit(1, "mm"), column_title=NULL,
          bottom_annotation=features_,
          top_annotation=sex_anno
  )
lgd=Legend(col_fun=col_fun,
           title=expression(Spearman~rho),
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
CD4T_final=draw(ht_CD4, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_CD8T_Fvs.M.pdf", height=4.5, width=6)
draw(CD4T_final)
dev.off()

#####################################



### Try determining the functions of CD4T and CD8T subsets by analyzing the top age-cor.genes expressed by them
#####################################
###

### Load the age correlation results
cor_result=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
CD4_CD8T_subsets=cor_result$celltypes %>% unique(.) %>% .[grepl("^CD4T\\.|^CD8T\\.",.)]
cor_result_subset=cor_result %>% 
  subset(celltypes %in% CD4_CD8T_subsets & analysis=="all" & rho_pval<0.05) %>%
  split(.$celltypes)

### Enrich the age-correlated genes for each subset
gene_lists=lapply(cor_result_subset, function(df) df$gene)
enrich_result=clusterProfiler::compareCluster(gene_lists, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(enrich_result, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_EnrichResults_Detailed_CD4TandCD8T.rds")
clusterProfiler::dotplot(enrich_result)


#####################################



# ### Analyze the VDJ recombination in T and B cells
# #####################################
# ###
# 
# library(dplyr)
# library(ggplot2)
# 
# ### Load the correlation data
# donor_cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# donor_cor_mixedsex=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes))
# donor_cor_female=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes) & analysis=="females")
# donor_cor_male=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes) & analysis=="males")
# 
# ### Check the VDJ-related genes
# IgRecom_genes=c("ATM","CYREN","DCAF1","DCLRE1C","EZH2","FOXP1","HMGB1","LIG4","NHEJ1","POLB","PRKDC","RAG1","RAG2","TCF3","XRCC4","YY1","XRCC6")
# TRecom_genes=c("ATM","BCL11B","DCAF1","DCLRE1C","HMGB1","LEF1","LIG4","PRKDC","RAG1","RAG2","TCF7","XRCC6")
# Ig_gene_exist=IgRecom_genes[IgRecom_genes %in% donor_cor_mixedsex$gene]
# Tcell_gene_exist=TRecom_genes[TRecom_genes %in% donor_cor_mixedsex$gene]
# 
# ### Extract the expression data
# pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# TB_expr_data=Seurat::FetchData(pseudobulk_obj, vars=c("donor_id","sex","age","Annot.detailed","Annot.inter","Annot.rough",Ig_gene_exist,Tcell_gene_exist),
#                                layer="data")
# 
# ### Plot the VDJ-related gene expr in T naive cells
# TB_expr_data_df=TB_expr_data %>%
#   tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot", colnames(.))], names_to="gene", values_to="expr") %>%
#   subset(gene %in% Tcell_gene_exist) %>%
#   group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% summarize_at("expr", mean) %>%
#   mutate(sex=ifelse(sex=="female","F","M"))
# TB_expr_data_df$sex=forcats::fct_relevel(TB_expr_data_df$sex, c("F","M"))
# 
# VDJ_T_plot=
#   ggplot(subset(TB_expr_data_df, Annot.inter %in% c("CD4T.naive","CD8T.naive")), aes(x=age, y=expr, color=sex)) +
#   ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.05) +
#   facet_wrap(~Annot.inter, ncol=1) +
#   theme_classic() +
#   labs(x="age", y="Expression", title=NULL) +
#   theme(axis.text.x=element_text(size=10),
#         axis.text.y=element_text(size=10),
#         axis.title.x=element_text(size=11),
#         axis.title.y=element_text(size=11),
#         legend.position="top",
#         legend.title=element_blank(),
#         legend.text=element_text(size=10),
#         title=element_text(size=10),
#         strip.text=element_text(size=11)) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#   geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
#   ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=4)
# 
# ### Plot the VDJ-related gene expr in B.naive and B.inter cells
# TB_expr_data_df=TB_expr_data %>%
#   tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot", colnames(.))], names_to="gene", values_to="expr") %>%
#   subset(gene %in% Ig_gene_exist) %>%
#   group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% summarize_at("expr", mean) %>%
#   mutate(sex=ifelse(sex=="female","F","M"))
# TB_expr_data_df$sex=forcats::fct_relevel(TB_expr_data_df$sex, c("F","M"))
# 
# VDJ_Bcell_plot=
#   ggplot(subset(TB_expr_data_df, Annot.inter %in% c("B.naive","B.inter") & Annot.detailed!="B.inter.lambda.IL4R"), # checked, B.inter.lambda.IL4R has too few observations
#          aes(x=age, y=expr, color=Annot.detailed)) +
#   ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.05) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
#   # facet_wrap(~Annot.detailed, ncol=1) +
#   theme_classic() +
#   labs(x="age", y="Expression", title=NULL) +
#   theme(axis.text.x=element_text(size=10),
#         axis.text.y=element_text(size=10),
#         axis.title.x=element_text(size=11),
#         axis.title.y=element_text(size=11),
#         legend.position="top",
#         legend.title=element_blank(),
#         legend.text=element_text(size=10),
#         title=element_text(size=10)) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#   geom_smooth(aes(group=Annot.detailed), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
#   ggpubr::stat_cor(aes(group=Annot.detailed), show.legend=FALSE, size=4, label.x.npc=0.35, label.y.npc=1)
# 
# ### Plot both
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_VDJrecombination_TandB.pdf", height=5, width=7)
# cowplot::plot_grid(plotlist=list(VDJ_T_plot, VDJ_Bcell_plot), ncol=2, align="hv", axis="tb", rel_widths=c(1,2))
# dev.off()
# 
# #####################################



### Apoptosis analysis
#####################################
###

library(dplyr)
library(ggplot2)

### Apoptosis-related genes
anti_apoptosis_genes=c("ADA","ARG2","AURKB","AXL","BCL2","BCL2A1","BCL2L1","BCL3","BCL10","BCL11B","BLM","BMP4","CCL5","CCL19","CCL21",
                       "CCR5","CCR7","CD27","CD44","CD74","CXCL12","CXCR2","DOCK8","EFNA1","FADD","FCER1G","FCGR2B","FCMR","FOXP1","GAS6",
                       "GHSR","GPAM","HCLS1","HIF1A","HSH2D","IDO1","IL2","IL3","IL7R","IL18","IRS2","ITPKB","JAK3","KIFAP3","KITLG","MERTK",
                       "MIF","NOC2L","NOD2","ORMDL3","PDCD1","PIP","PNP","PRKCQ","PTCRA","RAG1","RORC","SELENOS","SERPINB9","SLC39A10","SLC46A2",
                       "ST3GAL1","ST6GAL1","STAT5A","TNFRSF4","TNFSF4","TSC22D3","VHL")
pro_apoptosis_genes=c("ADAM8","ANXA1","BAX","BBC3","BCL2L11","CCL5","CD24","CD44","CD274","CDKN2A","CEACAM1","FNIP1","HCAR2","IDO1","IL10","JAK3","LYN",
                      "MEF2C","MYC","NF1","NFKBID","NR4A3","P2RX7","PDCD1","PDCD7","PIK3CB","PIK3CD","PRELID1","RAPGEF2","SIGLEC1","SIRT1","TGFB2","TP53",
                      "WNT5A","ZC3H8")

### Load the expr data
# take only the subpopulations that exist in Immunity dataset
Immunity_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
confirm_celltypes=unique(Immunity_obj$Annot.detailed)
pseudobulk_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
ap_expr_data=Seurat::FetchData(pseudobulk_obj, vars=c("donor_id","sex","age","Annot.detailed","Annot.inter","Annot.rough",
                                                      anti_apoptosis_genes, pro_apoptosis_genes),
                               layer="data")

### Clean the data
genes_deselected=names(colSums(ap_expr_data[,7:ncol(ap_expr_data)])!=0)[colSums(ap_expr_data[,7:ncol(ap_expr_data)])==0]
ap_expr_data_df=ap_expr_data %>% select(-any_of(genes_deselected)) %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot",colnames(.))],
                                                     names_to="gene", values_to="expr") %>%
  mutate(pro_or_anti=ifelse(gene %in% intersect(anti_apoptosis_genes, pro_apoptosis_genes), "regulate",
                            ifelse(gene %in% anti_apoptosis_genes, "anti-apoptotic", "pro-apoptotic")))

### Plot
ap_expr_data_df_sum=ap_expr_data_df %>% group_by(across(-c(gene,expr))) %>%
  summarize_at("expr", mean) %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  subset(pro_or_anti!="regulate") %>%
  subset(Annot.rough %in% c("CD4T cells","CD8T cells") & Annot.detailed %in% confirm_celltypes) %>%
  tidyr::pivot_wider(names_from="pro_or_anti", values_from="expr") %>%
  mutate(pro_div_anti=`pro-apoptotic`/`anti-apoptotic`) %>%
  mutate(Annot.detailed=stringr::str_wrap(gsub("\\."," ",Annot.detailed), 15))

all_plot=
  ggplot(ap_expr_data_df_sum, aes(x=age, y=pro_div_anti, color=sex)) +
  facet_wrap(~Annot.detailed, ncol=4) +
  ggrastr::geom_point_rast(size=0.1, alpha=0.1, shape=20) +
  scale_x_continuous(breaks=c(30,60,90), limits=c(15,100)) +
  labs(title=NULL, y="Expression", x="age") +
  theme_classic() +
  geom_hline(yintercept=1, linetype="dotted", linewidth=0.5) +
  geom_smooth(aes(group=sex), se=F, linewidth=0.5, method="loess") +
  theme(title=element_text(size=11),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        panel.grid=element_line(colour="grey80", size=0.05),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text=element_text(size=10),
        legend.direction="horizontal",
        legend.box="horizontal",
        legend.position="top") +
  guides(color=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Apoptosis_CD4TandCD8T.pdf", height=5.5, width=6.5)
plot(all_plot)
dev.off()

#####################################



