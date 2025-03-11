

### Combine all 3 verification datasets for performance analysis of Xu's method
#####################################
###

library(Seurat)
library(dplyr)

CovidCtrl=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
CovidCtrl[["RNA"]]=split(CovidCtrl[["RNA"]], f=CovidCtrl$donor_id)
PBMC17=readRDS("~/Project_PBMCage/GSE213516_17PBMC/17_pbmc_annotated_AssignedOnly.rds")
PBMC17[["RNA"]]=split(PBMC17[["RNA"]], f=PBMC17$donor_id)
JapanSC=readRDS("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated.rds")
JapanSC[["RNA"]]=split(JapanSC[["RNA"]], f=JapanSC$donor_id)
# merge
merged_all=merge(CovidCtrl, c(PBMC17, JapanSC), add.cell.ids=c("CovidCtrl","PBMC17","JapanSC"))
saveRDS(merged_all, "~/Project_PBMCage/tempt.rds")
# integrate
merged_all=merged_all %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  FindClusters(resolution=0.5)
merged_all_integrated=IntegrateLayers(merged_all, method=RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca")
merged_all_integrated=merged_all_integrated %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  FindClusters(resolution=0.5, cluster.name="rpca_clusters") %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30, reduction.name="umap.rpca")
DimPlot(merged_all_integrated, reduction="umap.rpca", group.by="donor_id", label.size=2)
merged_all_integrated=JoinLayers(merged_all_integrated)
merged_all_integrated=ScaleData(merged_all_integrated)
DimPlot(merged_all_integrated, group.by="rpca_clusters", reduction="umap.rpca", label=T, label.size=2) + NoLegend()
saveRDS(merged_all_integrated, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Integrated.rds")

#####################################



### Make pseudobulk obj
#####################################
###

library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Integrated.rds")
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
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")
# -- add sex metadata
pbmc.merged=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(tempt_annot_df$donor_id, as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("M","F")
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_rough.rds")

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
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_inter.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_inter.rds")
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
tempt_annot_df$sex=sex_df$sex[match(tempt_annot_df$donor_id, as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("M","F")
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_inter.rds")

# pseudobulk at the detailed level
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
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_detailed.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_detailed.rds")
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
tempt_annot_df$sex=sex_df$sex[match(tempt_annot_df$donor_id, as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("M","F")
saveRDS(pbmc.merged, "~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_detailed.rds")

#####################################