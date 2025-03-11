
### Process GSE158055 (the COVID dataset) - Frozen_pbmc and Fresh_pbmc
#####################################
###

### Process Frozen_pbmc

library(SeuratDisk)
filepath="~/Project_PBMCage/GSE158055/COVID19_frozen_PBMC.h5ad"
tryCatch({Convert(filepath, dest="h5seurat", overwrite=FALSE)}, error=function(msg) {print("Pass.")})
frozen_pbmc=LoadH5Seurat("~/Project_PBMCage/GSE158055/COVID19_frozen_PBMC.h5seurat",
                         meta.data=FALSE, misc=FALSE)
frozen_pbmc
# add back the metadata
library(rhdf5)
info=h5ls(filepath)
table(colnames(frozen_pbmc)==h5read(filepath, "/obs/_index")) # check
metadata_category=info$group[grepl("\\/obs\\/",info$group) & info$dclass=="STRING"]
name_for_metadata=gsub("\\/obs\\/","",metadata_category)
name_for_metadata=gsub("-| ","_",name_for_metadata)
name_for_metadata=gsub("\\[.*","level",name_for_metadata)
for (i in 1:length(metadata_category)) {
  codes=h5read(filepath, paste0(metadata_category[i],"/codes"))+1
  cat=h5read(filepath, paste0(metadata_category[i],"/categories"))
  temp=cat[codes]
  frozen_pbmc[[name_for_metadata[i]]]=temp
}
# set the age variable as an ordered factor
frozen_pbmc$Age=as.factor(frozen_pbmc$Age)
ordered_levels=as.numeric(levels(frozen_pbmc$Age))
frozen_pbmc$Age=factor(frozen_pbmc$Age, levels=ordered_levels[order(ordered_levels)])
# have a look at the embedding and metadata
library(Seurat)
DimPlot(frozen_pbmc, group.by="majorType")
DimPlot(frozen_pbmc, group.by="celltype", label=F) + NoLegend()
head(frozen_pbmc[[]])
table(frozen_pbmc$majorType)
table(frozen_pbmc$celltype)
table(frozen_pbmc$Age)
table(frozen_pbmc$CoVID_19_severity)
table(frozen_pbmc$SARS_CoV_2)
# save the frozen_pbmc seuratobj
saveRDS(frozen_pbmc, "~/Project_PBMCage/GSE158055/frozen_pbmc.rds")

### Process Fresh_pbmc

library(SeuratDisk)
filepath="~/Project_PBMCage/GSE158055/COVID19_fresh_PBMC.h5ad"
tryCatch({Convert(filepath, dest="h5seurat", overwrite=FALSE)}, error=function(msg) {print("Pass.")})
fresh_pbmc=LoadH5Seurat("~/Project_PBMCage/GSE158055/COVID19_fresh_PBMC.h5seurat",
                        meta.data=FALSE, misc=FALSE)
fresh_pbmc
# add back the metadata
library(rhdf5)
info=h5ls(filepath)
table(colnames(fresh_pbmc)==h5read(filepath, "/obs/_index")) # check
metadata_category=info$group[grepl("\\/obs\\/",info$group) & info$dclass=="STRING"]
name_for_metadata=gsub("\\/obs\\/","",metadata_category)
name_for_metadata=gsub("-| ","_",name_for_metadata)
name_for_metadata=gsub("\\[.*","level",name_for_metadata)
for (i in 1:length(metadata_category)) {
  codes=h5read(filepath, paste0(metadata_category[i],"/codes"))+1
  cat=h5read(filepath, paste0(metadata_category[i],"/categories"))
  temp=cat[codes]
  fresh_pbmc[[name_for_metadata[i]]]=temp
}
# set the age variable as an ordered factor
fresh_pbmc$Age=as.factor(fresh_pbmc$Age)
ordered_levels=as.numeric(levels(fresh_pbmc$Age))
fresh_pbmc$Age=factor(fresh_pbmc$Age, levels=ordered_levels[order(ordered_levels)])
# have a look at the embedding and metadata
library(Seurat)
DimPlot(fresh_pbmc, group.by="majorType")
DimPlot(fresh_pbmc, group.by="celltype", label=F) + NoLegend()
head(fresh_pbmc[[]])
table(fresh_pbmc$majorType)
table(fresh_pbmc$celltype)
table(fresh_pbmc$Age)
table(fresh_pbmc$CoVID_19_severity)
table(fresh_pbmc$SARS_CoV_2)
# save the fresh_pbmc seuratobj
saveRDS(fresh_pbmc, "~/Project_PBMCage/GSE158055/fresh_pbmc.rds")

### Merge Frozen_pbmc and Fresh_pbmc

frozen_fresh_pbmc=merge(frozen_pbmc, fresh_pbmc, add.cell.ids=c("frozen", "fresh"))
frozen_fresh_counts=LayerData(frozen_fresh_pbmc, assay="RNA", layer="counts", fast=TRUE)

PBMC_meta=frozen_fresh_pbmc[[]]
colnames(PBMC_meta)
PBMC_meta_chosen=PBMC_meta[, c(14,11,16,1)]
PBMC_meta_chosen$Sex=factor(PBMC_meta_chosen$Sex, levels=c("M","F"))
PBMC_meta_chosen$Age=as.factor(PBMC_meta_chosen$Age)
colnames(PBMC_meta_chosen)=c("dataset","donor_id","sex","age")

pbmc.frozen_fresh=CreateSeuratObject(counts=frozen_fresh_counts, meta.data=PBMC_meta_chosen)

# add percent.rb and percent.mt
pbmc.frozen_fresh[["percent.rb"]]=PercentageFeatureSet(pbmc.frozen_fresh, pattern="^RP[SL]")
pbmc.frozen_fresh[["percent.mt"]]=PercentageFeatureSet(pbmc.frozen_fresh, pattern="^MT-")

saveRDS(pbmc.frozen_fresh, "~/Project_PBMCage/GSE158055/frozen_fresh_pbmc.rds")

# try merge without integration
pbmc.frozen_fresh_try=pbmc.frozen_fresh %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30) %>%
  FindClusters(resolution=0.5)
DimPlot(pbmc.frozen_fresh_try, group.by=c("dataset", "donor_id", "seurat_clusters"))
frozen.pbmc=subset(pbmc.frozen_fresh_try, dataset=="frozen PBMC") # or fresh PBMC
DimPlot(frozen.pbmc, group.by=c("donor_id", "seurat_clusters"))

# have to integrate
pbmc.frozen_fresh_try[["RNA"]]=split(pbmc.frozen_fresh[["RNA"]], f=pbmc.frozen_fresh$donor_id)
pbmc.frozen_fresh_try
pbmc.frozen_fresh_try=pbmc.frozen_fresh_try %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  FindClusters(resolution=0.5)
pbmc.frozen_fresh_try2=IntegrateLayers(pbmc.frozen_fresh_try, method=RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca")
pbmc.frozen_fresh_try3=pbmc.frozen_fresh_try2 %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  FindClusters(resolution=0.5, cluster.name="rpca_clusters") %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30, reduction.name="umap.rpca")
DimPlot(pbmc.frozen_fresh_try3, reduction="umap.rpca", group.by=c("dataset", "donor_id"), label.size=2)
saveRDS(pbmc.frozen_fresh_try3, "~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

### Annotate the cells
pbmc.frozen_fresh_try3=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
pbmc.seu_merged=JoinLayers(pbmc.frozen_fresh_try3)
DimPlot(pbmc.seu_merged, group.by="rpca_clusters", reduction="integrated.rpca", label=T, label.size=2) + NoLegend()
markers=FindAllMarkers(pbmc.seu_merged, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
saveRDS(markers, "~/Project_PBMCage/Tempt_RDS/markersTempt.rds")
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
markers2=FindMarkers(pbmc.seu_merged, ident.1=4, ident.2=5)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)
#0: CD14 Mono
#1: CD4 Tnaive
#2: CD4 Tcm
#3: CD8 Tem GZMB+
#4：NK
#5: NK
#6: CD8 naive
#7: B.naive
#8: CD8 Tem
#9: gdT
#10: CD14 Mono
#11: Bmem
#12: CD16 Mono
#13: MAIT
#14: CD4 Tcm(Trm?)
#15: Platelet+DC?
#16: cDC2
#17: CD4 Tnaive (?)
#18: Platelet
#19: NK.prolif
#20: pDC
#21: PC
#22: cDC1
saveRDS(pbmc.seu_merged, "~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

# tempt annot
tempt_annot_df=pbmc.seu_merged[[]]
tempt_annot_df$cell=rownames(tempt_annot_df)
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(1,2,14,17), "CD4T cells", seurat_clusters)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(3,6,8), "CD8T cells", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(7,11,21), "B cells", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(16,20,22), "DCs", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(0,10,12), "Monocytes", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(4,5,19), "NK cells", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(15,18), "Other cells", temptannot)) %>%
  mutate(temptannot=ifelse(seurat_clusters %in% c(9,13), "OtherT", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")
table(rownames(tempt_annot_df)==colnames(pbmc.seu_merged)) # check
pbmc.seu_merged=AddMetaData(pbmc.seu_merged, metadata=tempt_annot_df$temptannot, col.name="temptannot")
DimPlot(pbmc.seu_merged, reduction="umap.rpca", group.by="temptannot", label=T) + NoLegend()

B_subset=subset(pbmc.seu_merged, temptannot=="B cells")
B_subset=B_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
B_subset_try=FindClusters2(B_subset, cluster.range=c(3,4), by=0.1, res=0.2, verbose=T)
DimPlot(B_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(B_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(B_subset_try, idents=3)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(B_subset_try, idents=0), "B.naive.kappa", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(B_subset_try, idents=1), "B.mem.IGHA", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(B_subset_try, idents=2), "B.mem.IGHM", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(B_subset_try, idents=3), "B.pc", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")

CD4_subset=subset(pbmc.seu_merged, temptannot=="CD4T cells")
CD4_subset=CD4_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
CD4_subset_try=FindClusters2(CD4_subset, cluster.range=c(6,7), by=0.1, res=0.2, verbose=T)
DimPlot(CD4_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(CD4_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
markers2=FindMarkers(CD4_subset_try, ident.1=0, ident.2=4)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)
highlighted=WhichCells(CD4_subset_try, idents=0)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=0), "CD4T.naive.COTL1pos_SOX4neg", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=1), "CD4T.Tcm.COTL1", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=2), "CD4T.Tem.KLRB1", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=3), "CD4T.Treg", temptannot)) %>% #3:CD4 Treg (Tcm/Trm)
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=4), "CD4T.naive.HNRNPH1", temptannot)) %>% #4:CD4T.naive->Tcm trans
  mutate(temptannot=ifelse(cell %in% WhichCells(CD4_subset_try, idents=5), "CD4T.Tcm.ISG", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")

CD8_subset=subset(pbmc.seu_merged, temptannot=="CD8T cells")
CD8_subset=CD8_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
CD8_subset_try=FindClusters2(CD8_subset, cluster.range=c(7,8), by=0.1, res=0.5, verbose=T)
DimPlot(CD8_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(CD8_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(CD8_subset_try, idents=3)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=0), "CD8T.Tem.GZMB_FCGR3A", temptannot)) %>% #0:CD8 Tem
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=1), "CD8T.naive.LINC02446", temptannot)) %>% #1:CD8 Tnaive
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=2), "CD8T.Tem.GZMKhi", temptannot)) %>% #2:CD8 Tem.GZMKhi
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=3), "CD4T.CTL", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=4), "CD8T.Tem.GZMK_NKG7", temptannot)) %>% #2:CD8 Tem?
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=5), "CD8T.naive.LRRN3", temptannot)) %>% #1:CD8 Tnaive
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=6), "CD8T.Tcm.LYAR", temptannot)) %>% #3:CD8 Tcm/Trm
  mutate(temptannot=ifelse(cell %in% WhichCells(CD8_subset_try, idents=7), "CD8T.Tem.GZMB_ITGB1", temptannot)) #0:CD8 Tem
table(tempt_annot_df$temptannot, useNA="ifany")

DCs_subset=subset(pbmc.seu_merged, temptannot=="DCs")
DCs_subset=DCs_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
DCs_subset_try=FindClusters2(DCs_subset, cluster.range=c(3,5), by=0.1, res=0.2, verbose=T)
DimPlot(DCs_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(DCs_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(DCs_subset_try, idents=4)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(DCs_subset_try, idents=0), "DC.cDC2_2", temptannot)) %>% #0:cDC2_2
  mutate(temptannot=ifelse(cell %in% WhichCells(DCs_subset_try, idents=1), "DC.cDC2_1", temptannot)) %>% #1:cDC2_1
  mutate(temptannot=ifelse(cell %in% WhichCells(DCs_subset_try, idents=2), "DC.pDC", temptannot)) %>% #2:pDC
  mutate(temptannot=ifelse(cell %in% WhichCells(DCs_subset_try, idents=3), "DC.cDC1", temptannot)) %>% #3:cDC1
  mutate(temptannot=ifelse(cell %in% WhichCells(DCs_subset_try, idents=4), "DC.AXL_mDC", temptannot)) #4:ASDC
table(tempt_annot_df$temptannot, useNA="ifany")

Monocytes_subset=subset(pbmc.seu_merged, temptannot=="Monocytes")
Monocytes_subset=Monocytes_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
Monocytes_subset_try=FindClusters2(Monocytes_subset, cluster.range=c(4,5), by=0.1, res=0.1, verbose=T)
DimPlot(Monocytes_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(Monocytes_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(Monocytes_subset_try, idents=3)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Monocytes_subset_try, idents=0), "Mono.classical.HLA", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Monocytes_subset_try, idents=1), "Mono.classical", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Monocytes_subset_try, idents=2), "Mono.nonclassical", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Monocytes_subset_try, idents=3), "CD4T.naive.COTL1pos_SOX4neg", temptannot)) #4:CD4T.naive->Tcm trans
table(tempt_annot_df$temptannot, useNA="ifany")

NK_subset=subset(pbmc.seu_merged, temptannot=="NK cells")
NK_subset=NK_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
NK_subset_try=FindClusters2(NK_subset, cluster.range=c(5,6), by=0.1, res=0.2, verbose=T)
DimPlot(NK_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(NK_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(NK_subset_try, idents=3)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(NK_subset_try, idents=0), "NK.CD56dim.FCER1G", temptannot)) %>% #0:NK_2
  mutate(temptannot=ifelse(cell %in% WhichCells(NK_subset_try, idents=1), "NK.CD56dim.FCER1Glo_KLRC2pos", temptannot)) %>% #1:NK_1
  mutate(temptannot=ifelse(cell %in% WhichCells(NK_subset_try, idents=2), "NK.CD56hi.FCGR3A", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(NK_subset_try, idents=3), "OtherT.NKT", temptannot)) %>% #3: CD8 NKT
  mutate(temptannot=ifelse(cell %in% WhichCells(NK_subset_try, idents=4), "NK.prolif", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")

Othercells_subset=subset(pbmc.seu_merged, temptannot=="Other cells")
Othercells_subset=Othercells_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
Othercells_subset_try=FindClusters2(Othercells_subset, cluster.range=c(4,5), by=0.1, res=0.1, verbose=T)
DimPlot(Othercells_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(Othercells_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(Othercells_subset_try, idents=0)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Othercells_subset_try, idents=0), "CD4T.naive.COTL1pos_SOX4neg", temptannot)) %>% # CD4T.naive->Tcm trans
  mutate(temptannot=ifelse(cell %in% WhichCells(Othercells_subset_try, idents=1), "Other.Plt", temptannot)) %>% # Platelet
  mutate(temptannot=ifelse(cell %in% WhichCells(Othercells_subset_try, idents=2), "NK.CD56dim.FCER1Glo_KLRC2neg", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(Othercells_subset_try, idents=3), "CD8T.Tem.GZMKhi", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")

OtherT_subset=subset(pbmc.seu_merged, temptannot=="OtherT")
OtherT_subset=OtherT_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(reduction="integrated.rpca", dims=1:30) %>%
  RunUMAP(reduction="integrated.rpca", dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
OtherT_subset_try=FindClusters2(OtherT_subset, cluster.range=c(2,3), by=0.1, res=0.2, verbose=T)
DimPlot(OtherT_subset_try, label=T, label.size=5, repel=T) + NoLegend()
markers=FindAllMarkers(OtherT_subset_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
highlighted=WhichCells(OtherT_subset_try, idents=1)
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, repel=T, group.by="temptannot", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
tempt_annot_df=tempt_annot_df %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(OtherT_subset_try, idents=0), "OtherT.MAIT", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(OtherT_subset_try, idents=1), "OtherT.gdT.CMC1", temptannot)) %>%
  mutate(temptannot=ifelse(cell %in% WhichCells(OtherT_subset_try, idents=2), "OtherT.gdT.KLRC1", temptannot))
table(tempt_annot_df$temptannot, useNA="ifany")

### Add the metadata of celltype annotation
celltypes_detailed=unique(tempt_annot_df$temptannot)
tempt_annot_df=tempt_annot_df %>%
  mutate(Annot.inter=ifelse(temptannot %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", temptannot)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.inter %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter))
table(tempt_annot_df$Annot.inter, useNA="ifany")
tempt_annot_df=tempt_annot_df %>%
  mutate(Annot.rough=ifelse(temptannot %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", temptannot)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
table(tempt_annot_df$Annot.rough, useNA="ifany")

pbmc.seu_merged=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
table(rownames(tempt_annot_df)==colnames(pbmc.seu_merged)) # check
pbmc.seu_merged=AddMetaData(pbmc.seu_merged, metadata=tempt_annot_df$temptannot, col.name="Annot.detailed")
pbmc.seu_merged=AddMetaData(pbmc.seu_merged, metadata=tempt_annot_df$Annot.inter, col.name="Annot.inter")
pbmc.seu_merged=AddMetaData(pbmc.seu_merged, metadata=tempt_annot_df$Annot.rough, col.name="Annot.rough")
DimPlot(pbmc.seu_merged, reduction="umap.rpca", label=T, label.size=2, group.by="Annot.detailed") + NoLegend()
saveRDS(pbmc.seu_merged, "~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

# add the agecut metadata
AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}
agecut_df=pbmc.seu_merged[[]]
agecut_df$age=as.factor(agecut_df$age)
agecut_df$agecut=agecut_df$age
for (i in 1:length(AGECUT)) {
  agecut_df=agecut_df %>%
    mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
}
table(rownames(agecut_df)==colnames(pbmc.seu_merged)) # check
pbmc.seu_merged=AddMetaData(pbmc.seu_merged, metadata=agecut_df$agecut, col.name="agecut")
pbmc.seu_merged$agecut=as.factor(pbmc.seu_merged$agecut)
saveRDS(pbmc.seu_merged, "~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

#####################################



### Add cell cycle information for all the cells including unassigned
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

###
covidctrl=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
covidctrl[[]] # percent.rb and percent.mt has been added
covidctrl # make sure the seurat_obj has been normalize, HVGs found, and scaled

### Cellcycle scoring
covidctrl=CellCycleScoring(covidctrl, 
                           s.features=cc.genes.updated.2019$s.genes, 
                           g2m.features=cc.genes.updated.2019$g2m.genes,
                           set.ident=TRUE)

### Save the mt, rb, and cell-cycle info
meta_info=covidctrl@meta.data %>% select(any_of(c("nCount_RNA","nFeature_RNA","dataset","donor_id","sex","age",
                                                  "percent.rb","percent.mt","Annot.detailed","Annot.inter","Annot.rough",
                                                  "agecut",
                                                  "S.Score","G2M.Score","Phase"))) %>%
  mutate(cellid=rownames(.))
data.table::fwrite(meta_info, "~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")

### Plot the mt
meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
# plot all the celltypes in general
meta_info_subset=meta_info %>% group_by(donor_id, sex, age) %>% summarize_at("percent.mt", mean)
plot_mt_all=
  ggplot(meta_info_subset, aes(x=age, y=percent.mt, color=sex)) +
  geom_point(alpha=0.25) +
  theme_classic() +
  labs(x="age", y="Percent of mitochondrial genes (%)", title="GSE158055") +
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
  labs(x="age", y="Percent of mitochondrial genes (%)", title="GSE158055") +
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
meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
# plot all the celltypes in general
meta_info_subset=meta_info %>% group_by(donor_id, sex, age) %>% summarize_at("percent.rb", mean)
plot_rb_all=
  ggplot(meta_info_subset, aes(x=age, y=percent.rb, color=sex)) +
  geom_point(alpha=0.25) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="GSE158055") +
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
  geom_point(alpha=0.25) +
  facet_wrap(~Annot.rough, ncol=4) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="GSE158055") +
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

pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Mt.Rb.Cellcycle2.pdf", height=4.5, width=7)
plot(plot_rb_sep)
dev.off()

### Plot the cell-cycle phase
meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
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
  labs(x="age", y="Percent of cell phase (%)", title="GSE158055") +
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
  labs(x="age", y="Percent of cell phase (%)", title="GSE158055") +
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

pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Mt.Rb.Cellcycle1.pdf", height=4.5, width=2.5)
# plot(plot_mt_all)
# plot(plot_mt_sep)
plot(plot_rb_all)
# plot(plot_rb_sep)
# plot(plot_phase_all)
# plot(plot_phase_sep)
dev.off()

#####################################



### Make pseudobulk obj
#####################################
###

library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
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
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_rough.rds")
# -- add sex metadata
pbmc.merged=readRDS("~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_rough.rds")
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(tempt_annot_df$donor_id, as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("M","F")
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_rough.rds")

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
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_inter.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_inter.rds")
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
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_inter.rds")

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
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_detailed.rds")
# -- add celltypes metadata
pbmc.merged=readRDS("~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_detailed.rds")
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
saveRDS(pbmc.merged, "~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_detailed.rds")

#####################################



# ### Draw boxplot of freq vs ages
# #####################################
# ###
# 
# ###
# pbmc.seu_merged=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
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
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=1, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Proportion_vs._age.pdf", width=4.5, height=15)
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
# pbmc.seu_merged=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
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
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=1, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Proportion_vs._agecut.pdf", width=4, height=15)
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
# pbmc.seu_merged=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
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
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=1, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Counts_vs._age.pdf", width=4.5, height=15)
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
# pbmc.seu_merged=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
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
#   eval(parse(text=paste0("plistcowplot[[i]]=cowplot::plot_grid(",codes,", ncol=1, align='hv')"))) # add labels='AUTO': draw "A","B"...labels; "auto"="a","b"...
# }
# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Counts_vs._agecut.pdf", width=4, height=15)
# for (i in 1:length(plistcowplot)) plot(plistcowplot[[i]])
# dev.off()
# 
# #####################################



### Analyze the celltype freq and their correlation with ages, similar to the above "Draw boxplot" steps
#####################################
###

### Load the meta data with celltypes
meta_info=data.table::fread("~/Project_PBMCage/Results/CovidCtrl_results/Key_meta_info.csv.gz", sep="\t")
meta_info=subset(meta_info, Annot.rough!="Unassigned")

### Calculate the freq at the rough, inter, detailed level, respectively; save the result
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
write.table(COR_combined, "~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_correlation.txt", sep="\t")
write.table(meta_info_combined, "~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_values.txt", sep="\t")

#####################################



### Plot the celltype freq at the Annot.rough level
#####################################
###

### Plot the cellratio
object=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

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

pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Cell_frequency_rough.pdf", width=6, height=3)
plot(plot_rough)
dev.off()

### Plot CD4/CD8 ratio
meta_info_combined=read.delim("~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined_CD4.CD8=meta_info_combined %>% subset(Annot.rough %in% c("CD4T cells", "CD8T cells")) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|rough",colnames(.))])) %>%
  subset(!duplicated(.)) %>%
  select(-N_per_rough) %>%
  tidyr::pivot_wider(names_from="Annot.rough", values_from="percent.rough") %>%
  mutate(CD4_CD8.ratio=`CD4T cells`/`CD8T cells`)
# add the agecut metadata
AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
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
  geom_point(aes(color=sex)) +
  scale_color_manual(values=c("coral3","skyblue3")) +
  theme_classic() +
  labs(x=NULL, y="CD4T/CD8T ratio", title=NULL) +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=1), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25, color="black") +
  # ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3, method="spearman") +
  ggpubr::stat_compare_means(method="anova", aes(label=paste0("p = ", after_stat(p.format))), color="black", label.y.npc=0.9)

pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Cell_frequency_rough.CD4_CD8_ratio.pdf", width=6, height=3)
plot(plot_CD4CD8ratio_agecut)
dev.off()

#####################################



### Plot the celltype freq at the Annot.inter level
#####################################
###

### Load the data
meta_info_combined=read.delim("~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|inter",colnames(.))])) %>%
  subset(!duplicated(.))
plot_inter=
  ggplot(meta_info_combined, aes(x=age, y=percent.inter, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.inter, ncol=8, scales="free_y") +
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
pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Cell_frequency_inter.pdf", width=15, height=6)
plot(plot_inter)
dev.off()

#####################################



### Plot the celltype freq at the Annot.detailed level
#####################################
###

### Load the value
meta_info_combined=read.delim("~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_values.txt", sep="\t")

### !: According to the Annot.inte plot, we focus on CD4T, CD8T, NK, OtherT
COR_combined=read.delim("~/Project_PBMCage/Results/CovidCtrl_results/Cell_frequency_correlation.txt", sep="\t")
celltype_selected=list(CD4T=COR_combined$celltypes[grepl("^CD4T\\.", COR_combined$celltypes)],
                       CD8T=COR_combined$celltypes[grepl("^CD8T\\.", COR_combined$celltypes)],
                       NK=COR_combined$celltypes[grepl("^NK\\.", COR_combined$celltypes)],
                       OtherT=COR_combined$celltypes[grepl("^OtherT\\.", COR_combined$celltypes)])

### Take only those with significant correlation with ages
plot_detailed=list()
for (i in 1:length(celltype_selected)) {
  # take the sig. celltype
  COR_combined_subset=COR_combined %>% subset(analysis %in% c("F","M") & celltype.level=="Annot.detailed" & celltypes %in% celltype_selected[[i]])
  celltypes_sig=COR_combined_subset %>% select(celltypes) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  if (length(celltypes_sig)==0) break
  
  # subset the data
  meta_info_combined_subset=meta_info_combined %>% select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
    subset(!duplicated(.)) %>%
    subset(Annot.detailed %in% celltypes_sig)
  
  # plot
  plot_detailed[[i]]=
    ggplot(meta_info_combined_subset, aes(x=age, y=percent.detailed, color=sex)) +
    ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
    facet_wrap(~Annot.detailed, nrow=2, scales="free_y") +
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
    ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3, method="spearman")
}

# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/Cell_frequency_detailed.pdf", width=15, height=6)
# for (i in 1:length(plot_detailed)) plot(plot_detailed[[i]])
# dev.off()

#####################################



### Take CD8T, NK, and NKT for cytotoxicity scoring
#####################################
###
library(UCell)
library(Seurat)
library(dplyr)

object=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")

### Score CD8T and NKT
CD8T.NKT_detailedcelltypes=names(table(object$Annot.detailed))[grepl("^CD8T\\.|^Other\\.NKT", names(table(object$Annot.detailed)))]
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
score_df=cbind(score_df, CD8T_Ucell[[]][, c("age","agecut","Annot.detailed","Annot.inter","Annot.rough")])
score_df$age=as.factor(score_df$age)
levels(score_df$age)=names(table(score_df$age))
score_df$agecut=as.factor(score_df$agecut)
levels(score_df$agecut)=names(table(score_df$agecut))
saveRDS(score_df,"~/Project_PBMCage/GSE158055/Tempt_RDS/FunctionScore_CD8TandNKT.rds")

### Score NK
NK_detailedcelltypes=names(table(object$Annot.detailed))[grepl("^NK\\.", names(table(object$Annot.detailed)))]
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
score_df=cbind(score_df, NK_Ucell[[]][, c("age","agecut","Annot.detailed","Annot.inter","Annot.rough")])
score_df$age=as.factor(score_df$age)
levels(score_df$age)=names(table(score_df$age))
score_df$agecut=as.factor(score_df$agecut)
levels(score_df$agecut)=names(table(score_df$agecut))
saveRDS(score_df,"~/Project_PBMCage/GSE158055/Tempt_RDS/FunctionScore_NK.rds")

# Make a function for plotting on ages
geneset_scoring_forAge=function(dataframe, celltype_as_LegendName, scale) {
  score_df=dataframe
  library(ComplexHeatmap)
  score_df_agemean=aggregate(score_df[,grepl("_UCell",colnames(score_df))], list(Age=score_df$age, Celltype=score_df$Annot.detailed), FUN=mean)
  Celltypes_in_df=unique(score_df_agemean$Celltype)
  df_full_list=list()
  for (i in 1:length(Celltypes_in_df)) {
    df_=score_df_agemean[score_df_agemean$Celltype==Celltypes_in_df[i],]
    age_not_there=as.integer(levels(score_df$age)[sapply(levels(score_df$age), function(x) !(x %in% df_$Age))])
    df_empty=data.frame(matrix(nrow=length(age_not_there), ncol=ncol(df_)))
    colnames(df_empty)=colnames(df_)
    df_empty$Age=age_not_there
    tryCatch({df_empty$Celltype=Celltypes_in_df[i]}, error=function(msg) {print("Pass.")})
    df_full=data.table::rbindlist(list(df_, df_empty))
    df_full=df_full[order(df_full$Age),]
    df_full_list[[i]]=df_full
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  df_full_combined=df_full_combined[, 3:ncol(df_full_combined)]
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges*1.5, 0, ranges*1.5), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-colMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/4, 0, ranges/4), c("#00007E", "white", "#7E0000"))}

  p_age=list()
  codes=""
  for (i in 1:length(df_full_list)) {
    df_=as.data.frame(df_full_list[[i]])
    rownames(df_)=df_$Age
    df_to_draw=df_[,3:ncol(df_)]
    colnames(df_to_draw)=gsub("_UCell", "", colnames(df_to_draw))
    celltype_of_df=names(table(df_$Celltype))
    df_to_draw=t(df_to_draw)
    if (scale==T) {df_to_draw=t(scale(t(df_to_draw)))} else {df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)}

    if (i==1) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                name=paste0("UCell of ", celltype_as_LegendName),
                heatmap_legend_param=list(title_position="leftcenter-rot",
                                          title_gp=gpar(fontsize=8, fontface="bold")),
                column_title=NULL,
                row_title=celltype_of_df,
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
                row_title=celltype_of_df,
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
# a function for plotting on AgeCuts
geneset_scoring_forAgeCuts=function(dataframe, celltype_as_LegendName, scale) {
  score_df=dataframe
  library(ComplexHeatmap)
  score_df_agemean=aggregate(score_df[,grepl("_UCell",colnames(score_df))], list(Age_cut=score_df$agecut, Celltype=score_df$Annot.detailed), FUN=mean)
  Celltypes_in_df=unique(score_df_agemean$Celltype)
  df_full_list=list()
  for (i in 1:length(Celltypes_in_df)) {
    df_=score_df_agemean[score_df_agemean$Celltype==Celltypes_in_df[i],]
    AgeCut_not_there=levels(score_df$agecut)[sapply(levels(score_df$agecut), function(x) !(x %in% df_$Age_cut))]
    df_empty=data.frame(matrix(nrow=length(AgeCut_not_there), ncol=ncol(df_)))
    colnames(df_empty)=colnames(df_)
    df_empty$Age_cut=AgeCut_not_there
    tryCatch({df_empty$Celltype=Celltypes_in_df[i]}, error=function(msg) {print("Pass.")})
    df_full=data.table::rbindlist(list(df_, df_empty))
    df_full=df_full[match(levels(score_df$agecut), df_full$Age_cut), ]
    df_full_list[[i]]=df_full
  }
  df_full_combined=data.table::rbindlist(df_full_list)
  df_full_combined=df_full_combined[, 3:ncol(df_full_combined)]
  if (scale==T) {
    df_full_combined=t(scale(t(df_full_combined)))
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges*1.5, 0, ranges*1.5), c("#00007E", "white", "#7E0000"))}
  if (scale==F) {
    df_full_combined=df_full_combined-colMeans(df_full_combined, na.rm=T)
    quantiles=quantile(df_full_combined, na.rm=TRUE)
    ranges=min(abs(ceiling(quantiles[1]*100)/100), abs(floor(quantiles[5]*100)/100))
    col_fun=circlize::colorRamp2(c(-ranges/4, 0, ranges/4), c("#00007E", "white", "#7E0000"))}

  p_age=list()
  codes=""
  for (i in 1:length(df_full_list)) {
    df_=as.data.frame(df_full_list[[i]])
    rownames(df_)=df_$Age_cut
    df_to_draw=df_[,3:ncol(df_)]
    colnames(df_to_draw)=gsub("_UCell", "", colnames(df_to_draw))
    celltype_of_df=names(table(df_$Celltype))
    df_to_draw=t(df_to_draw)
    if (scale==T) {df_to_draw=t(scale(t(df_to_draw)))} else {df_to_draw=df_to_draw-rowMeans(df_to_draw, na.rm=T)}

    if (i==1) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                name=paste0("UCell of ", celltype_as_LegendName),
                heatmap_legend_param=list(title_position="leftcenter-rot",
                                          title_gp=gpar(fontsize=8, fontface="bold")),
                column_title=NULL,
                column_names_gp=gpar(fontsize=8),
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=T,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
                width=ncol(df_to_draw)*unit(3,"mm"),
                height=nrow(df_to_draw)*unit(3,"mm"))}
    if (i>=2) {
      p_age[[i]]=
        Heatmap(df_to_draw,
                na_col="gray99",
                col=col_fun,
                show_heatmap_legend=FALSE,
                column_title=NULL,
                column_names_gp=gpar(fontsize=8),
                row_title=celltype_of_df,
                row_title_gp=gpar(fontsize=8),
                row_title_rot=0,
                show_column_names=T,
                cluster_columns=F,
                cluster_rows=F,
                show_row_dend=F,
                row_names_gp=gpar(fontsize=8),
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

### Plot
# NK
NK_df=readRDS("~/Project_PBMCage/GSE158055/Tempt_RDS/FunctionScore_NK.rds")
NK_scoring_plots_forAge=geneset_scoring_forAge(dataframe=NK_df, celltype_as_LegendName="NK cells", scale=F)
NK_scoring_plots_forAgeCuts=geneset_scoring_forAgeCuts(dataframe=NK_df, celltype_as_LegendName="NK cells", scale=F)
# CD8T and NKT
CD8T.NKT_df=readRDS("~/Project_PBMCage/GSE158055/Tempt_RDS/FunctionScore_CD8TandNKT.rds")
CD8T.NKT_scoring_plots_forAge=geneset_scoring_forAge(dataframe=CD8T.NKT_df, celltype_as_LegendName="CD8 T and NKT cells", scale=F)
CD8T.NKT_scoring_plots_forAgeCuts=geneset_scoring_forAgeCuts(dataframe=CD8T.NKT_df, celltype_as_LegendName="CD8 T and NKT cells", scale=F)

pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Scoring_CD8_and_NK_rowmeans.pdf", width=8, height=12)
NK_scoring_plots_forAge
NK_scoring_plots_forAgeCuts
CD8T.NKT_scoring_plots_forAge
CD8T.NKT_scoring_plots_forAgeCuts
dev.off()

# pdf("~/Project_PBMCage/Plots/CovidCtrl_Plots/CovidCtrl_Scoring_CD8_and_NK_scale.pdf", width=8, height=12)
# NK_scoring_plots_forAge
# NK_scoring_plots_forAgeCuts
# CD8T.NKT_scoring_plots_forAge
# CD8T.NKT_scoring_plots_forAgeCuts
# dev.off()

#####################################



