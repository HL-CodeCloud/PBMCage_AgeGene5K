
### Take the rest of the Tube_ids to annotate the whole obj
#####################################
###

Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
Tubes_used=unique(Seurat_whole_obj$Tube_id)
Tubes_used=paste0("TubeID_",Tubes_used,".rds")

### Take the selected objs to merge
Seurat_OBJs=list()
files=list.files("~/Project_PBMCage/Immunity/raw_counts_h5ad/")
files=files[grepl("\\.rds",files)]
files=files[!(files %in% Tubes_used)]
for (i in 1:length(files)) {
  Seurat_OBJs[[i]]=readRDS(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
}
# truncate to speed up merging
chuncks=ceiling(length(files)/10)
Seurat_OBJs_CHUNKs=list()
for (j in 1:(chuncks-1)) {
  CHUNK_=c(1:chuncks)[j]
  if (j<15) {
    Seurat_OBJs_chunk=Seurat_OBJs[(10*(CHUNK_-1)+1):(10*(CHUNK_-1)+10)]
    Seurat_OBJs_CHUNKs[[j]]=merge(Seurat_OBJs_chunk[[1]], Seurat_OBJs_chunk[-1], add.cell.ids=NULL)
  } else {
    Seurat_OBJs_chunk=Seurat_OBJs[(10*(CHUNK_-1)+1):length(Seurat_OBJs)]
    Seurat_OBJs_CHUNKs[[j]]=merge(Seurat_OBJs_chunk[[1]], Seurat_OBJs_chunk[-1], add.cell.ids=NULL)
  }
  print(paste0("CHUNK: ", j, " DONE!"))
}
Seurat_whole_obj=merge(Seurat_OBJs_CHUNKs[[1]], Seurat_OBJs_CHUNKs[-1], add.cell.ids=NULL)
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")

#####################################



### Reannotate CD4T
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="CD4+ T cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/cd4_t_cells/cd4_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset naive
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Naive")
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures()
CD4T_obj_subset=SketchData(CD4T_obj_subset, ncells=50000, method="LeverageScore", sketched.assay="sketch")
DefaultAssay(CD4T_obj_subset)="sketch"
DimPlot(CD4T_obj_subset, group.by="Batch") # no batch effect, thus no need to integrate

CD4T_obj_subset=CD4T_obj_subset %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset)
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T) + NoLegend()

# cluster0: CD4 Tnaive with proliferative property (under transformation to Tcm as CD44+ SELL+)
VlnPlot(CD4T_obj_subset_try, features=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0) # 1 3
# cluster1: CD4 Tnaive under cell regulation of events including activation and proliferation
VlnPlot(CD4T_obj_subset_try, features=c("CD44","SELL","IL7R","COTL1","SOCS3","ANXA1","TAGLN2"), pt.size=0) # 0 2 5
# cluster3: CD4 Tnaive under cell regulation of quiescence (repressing differentiation), reported to be downreg during aging
VlnPlot(CD4T_obj_subset_try, features=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"), pt.size=0) #4

CD4T_obj_subset_try=ProjectData(
  CD4T_obj_subset_try,
  assay="RNA",
  full.reduction="pca.full",
  sketched.assay="sketch",
  sketched.reduction="pca",
  umap.model="umap",
  dims=1:50,
  refdata=list(seurat_clusters_full="seurat_clusters")
)

DefaultAssay(CD4T_obj_subset_try)="RNA"
Idents(CD4T_obj_subset_try)="seurat_clusters_full"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=0), "CD4T.naive.COTL1pos_SOX4neg", seurat_clusters_full)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=1), "CD4T.naive.HNRNPH1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=2), "CD4T.naive.COTL1pos_SOX4neg", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=3), "CD4T.naive.HNRNPH1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=4), "CD4T.naive.COTL1pos_SOX4pos", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents=5), "CD4T.naive.COTL1pos_SOX4neg", annotation))

table(obj_annotation$annotation, useNA="ifany")
data.table::fwrite(obj_annotation %>% select(cell, annotation), "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Subset Treg
Treg_celltypes=unique(CD4T_obj$Cluster_names_detailed)[grepl("Treg",unique(CD4T_obj$Cluster_names_detailed))]
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% Treg_celltypes)
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, group.by="Batch") # no batch effect, thus no need to integrate

Idents(CD4T_obj_subset)="Cluster_names_detailed"
# cluster0: subsets of activated/effector Tregs:
# (HLA-DRhi subset + FOXP3hi effector subset which is FOXP3hi PTPRChi DDX17hi)
# Path I is defined as the developmen path from naïve Tregs ending with the FOXP3hi subset
VlnPlot(CD4T_obj_subset, features=c("FOXP3","IL2RA","HLA-DRA","PTPRC","DDX17"), pt.size=0)
# cluster1: naive Tregs (at naïve status with CCR7hi TCF7hi HLA-DRlow FOXP3low profile)
VlnPlot(CD4T_obj_subset, features=c("CCR7","TCF7"), pt.size=0)

obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Treg naive"), "CD4T.Treg", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Treg memory"), "CD4T.Treg.FOXP3hi", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Treg KLRB1+RORC+"), "CD4T.Treg.KLRB1_RORC", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Treg cytotoxic"), "CD4T.Treg.cytotoxic", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Check Naive-IFN
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Naive-IFN")
CD4T_obj_subset=NormalizeData(CD4T_obj_subset)
VlnPlot(CD4T_obj_subset, features=c("TCF7","SELL","LEF1","CCR7","IL7R","CD44","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","COTL1","SOCS3","ANXA1","TAGLN2"), pt.size=0) # This one.
VlnPlot(CD4T_obj_subset, features=c("TCF7","SELL","LEF1","SOX4","ID2","SATB1"), pt.size=0)

Idents(CD4T_obj_subset)="Cluster_names_detailed"
obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Naive-IFN"), "CD4T.naive.COTL1pos_SOX4neg", Cluster_numbers))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Check Tcm/Tem
Tmem_celltypes=unique(CD4T_obj$Cluster_names_detailed)[!grepl("Treg|Naive",unique(CD4T_obj$Cluster_names_detailed))]
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% Tmem_celltypes)
CD4T_obj_subset=NormalizeData(CD4T_obj_subset)
Idents(CD4T_obj_subset)="Cluster_names_detailed"
# distinguish Tcm and Tem
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R"), pt.size=0) # Tem: Terminal effector, Temra, Th22; Tcm: others

# cluster0: CD4 Tcm, CD4T.Tcm
VlnPlot(CD4T_obj_subset, features=c("CCR7","SELL","GPR183","PIK3IP1"), pt.size=0) # Th2, Tfh, Th1
# cluster1: CD4 Tcm with cell proliferative property, CD4T.Tcm.HNRNPH1
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","C1orf56","HNRNPH1","CDC42SE1"), pt.size=0) # No.
# cluster2: CD4 Tcm TH17-like, CD4T.Tcm.KLRB1
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","KLRB1","CTSH","AQP3","TIMP1"), pt.size=0) # Th17, Th1/Th17
# cluster3: CD4 Tcm with upregulated HLA-DR, CD4T.Tcm.HLA
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","HLA-DRA","CD38","IL2RA"), pt.size=0) # HLA-DR+ memory
# cluster1: CD4 Tcm under cell regulation (events eg. proliferation/cytokine production/apoptosis...), CD4T.Tcm.COTL1
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","COTL1","ITGB1","LGALS1","S100A11","ANXA1"), pt.size=0) # Exhausted-like memory
# cluster1: IFN-I T at the status of Tcm, CD4T.Tcm.ISG
VlnPlot(CD4T_obj_subset, features=c("LIMS1","MAL","LTB","AQP3","MX1","ISG15","ISG20"), pt.size=0) # No.
# cluster2: IFN-I T at the status of transcriptional repression, CD4T.Tcm.ISG.TSHZ2
VlnPlot(ISAsigT_NKT_try, features=c("MX1","ISG15","ISG20","TSHZ2","FHIT","RNASET2"), pt.size=0) # No.

# cluster0: CD4 Tem with Th17 transcriptional signature (KLRB1+ GZMK+), CD4T.Tem.KLRB1
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R","GZMK","KLRB1"), pt.size=0) # Terminal effector
# cluster2: CD4 Tem with transcriptional regulation property, CD4T.Tem.MALAT1
VlnPlot(CD4T_obj_subset, features=c("JUND","MALAT1","DDX17"), pt.size=0) # Temra, Th22

Idents(CD4T_obj_subset)="Cluster_names_detailed"
obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Th22"), "CD4T.Tem.MALAT1", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Th17"), "CD4T.Tcm.KLRB1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Th2"), "CD4T.Tcm", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Th1/Th17"), "CD4T.Tcm.KLRB1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tfh"), "CD4T.Tcm", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Terminal effector"), "CD4T.Tem.KLRB1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="HLA-DR+ memory"), "CD4T.Tcm.HLA", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Th1"), "CD4T.Tcm", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Temra"), "CD4T.Tem.MALAT1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Exhausted-like memory"), "CD4T.Tcm.COTL1", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

#####################################



### Reannotate CD8T
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="TRAV1-2- CD8+ T cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/conventional_cd8_t_cells/conventional_cd8_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset naive
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Naive")
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(7,9), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster1: homeostatic proliferation (spontaneous proliferation and in partial activation), CD8T.naive.HNRNPH1
VlnPlot(CD4T_obj_subset_try, features=c("LEF1","SELL","C1orf56","HNRNPH1"), pt.size=0) # No.
# cluster2: cell regulation (events eg. proliferation/cytokine production/migration/apoptosis...), CD8T.naive.LINC02446
VlnPlot(CD4T_obj_subset_try, features=c("CD8B","LDHB","LEF1","LINC02446","CD8A","S100B","ID2","TCF7","VIM","CCR7"), pt.size=0) # all...
# cluster3: transcriptional repressors, essential for maintaining cellular homeostasis and regulating, CD8T.naive.TSHZ2
VlnPlot(CD4T_obj_subset_try, features=c("TSHZ2","FHIT","RNASET2")) # No.

DefaultAssay(CD4T_obj_subset_try)="RNA"
Idents(CD4T_obj_subset_try)="Cluster_names_detailed"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="Naive"), "CD8T.naive.LINC02446", seurat_clusters_full))
table(obj_annotation$annotation, useNA="ifany")
OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Check Naive-IFN
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Naive-IFN")
CD4T_obj_subset=NormalizeData(CD4T_obj_subset)
VlnPlot(CD4T_obj_subset, features=c("LEF1","SELL","C1orf56","HNRNPH1"), pt.size=0)
VlnPlot(CD4T_obj_subset, features=c("CD8B","LDHB","LEF1","LINC02446","CD8A","S100B","ID2","TCF7","VIM","CCR7"), pt.size=0) # This one.
VlnPlot(CD4T_obj_subset, features=c("TSHZ2","FHIT","RNASET2"))

Idents(CD4T_obj_subset)="Cluster_names_detailed"
obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Naive-IFN"), "CD8T.naive.LINC02446", Cluster_numbers))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Check Tcm/Tem
Tmem_celltypes=unique(CD4T_obj$Cluster_names_detailed)[!grepl("Naive",unique(CD4T_obj$Cluster_names_detailed))]
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% Tmem_celltypes)
CD4T_obj_subset=NormalizeData(CD4T_obj_subset)
Idents(CD4T_obj_subset)="Cluster_names_detailed"
# distinguish Tcm and Tem
VlnPlot(CD4T_obj_subset, features=c("CD44","SELL","IL7R"), pt.size=0) 
# Tem: Tem GZMB+, Temra, Tem GZMK+, Trm, Tmem KLRC2+ (Tcm+Tem), NKT-like (Teff)
# Tcm: Tcm CCR4-, Tcm CCR4+, HLA-DR+
# CTL: Proliferative

# cluster0: mature cytotoxic, CD8T.Tem.GZMB_FCGR3A
VlnPlot(CD4T_obj_subset, features=c("FCGR3A","PRF1","GZMB"), pt.size=0) # Tem GZMB+, Temra, NKT-like
# cluster1: less mature, but more proliferative, CD8T.Tem.GZMB_HNRNPH1
VlnPlot(CD4T_obj_subset, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0)
# cluster0: terminal effector T (Tte)-like cells (also a subset of conventional CD8 Tmem), CD8T.Tem.GZMK_NKG7
VlnPlot(CD4T_obj_subset, features=c("GZMB","GNLY","NKG7","ZEB2","GZMA"), pt.size=0) # No.
# cluster1: CCR7- GZMKhi, the conventional Tem cells (cytolytic and generally lacked the memory markers LEF1, CD27, CD28, and CD127), CD8T.Tem.GZMKhi
VlnPlot(CD4T_obj_subset, features=c("CCR7","GZMK","LEF1","CD27","CD28"), pt.size=0) # Tem GZMK+
# cluster2: NKG2C+GZMB–XCL1+, a unique CD8 Tmem that decreases during aging, CD8T.Tem.GZMK_XCL1
VlnPlot(CD4T_obj_subset, features=c("ZNF683","KLRC2","GZMB","LEF1","XCL1"), pt.size=0) # Tmem KLRC2+

# cluster0: early differentiation, CD8T.Tcm.LYAR
VlnPlot(CD4T_obj_subset, features=c("GZMK","LYAR","NELL2"), pt.size=0) # Tcm CCR4-
# cluster1: cytotoxic Tcm (late differentiation maybe), CD8T.Tcm.GZMB
VlnPlot(CD4T_obj_subset, features=c("NKG7","GNLY","GZMH","GZMB"), pt.size=0) # HLA-DR+
# cluster2: CD8T.Tcm.GATA3
VlnPlot(CD4T_obj_subset, features=c("CD8B", "C1orf162", "IL7R", "GATA3", "YBX3", "KRT1", "CD8A", "CTSW", "INPP4B", "LTB"), pt.size=0) # Tcm CCR4+
# cluster3: transcriptional regulation, CD8T.Tcm.MALAT1
VlnPlot(CD4T_obj_subset, features=c("MALAT1","NEAT1","DDX17"), pt.size=0)

# cluster2: highly proliferating CD8+ effector CTL, most are Tc1, CD8T.CTL
VlnPlot(CD4T_obj_subset, features=c("GZMH","NKG7","ID2","IFNG", "MKI67"), pt.size=0) # Proliferative

Idents(CD4T_obj_subset)="Cluster_names_detailed"
obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tem GZMB+"), "CD8T.Tem.GZMB_FCGR3A", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Temra"), "CD8T.Tem.GZMB_FCGR3A", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="NKT-like"), "CD8T.Tem.GZMB_FCGR3A", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tem GZMK+"), "CD8T.Tem.GZMKhi", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tmem KLRC2+"), "CD8T.Tem.GZMK_XCL1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tcm CCR4-"), "CD8T.Tcm.LYAR", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Tcm CCR4+"), "CD8T.Tcm.GATA3", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="HLA-DR+"), "CD8T.Tcm.GZMB", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Proliferative"), "CD8T.CTL", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Trm"), "CD8T.Trm", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# further subset Tem GZMK+
further_subset=subset(CD4T_obj_subset, Cluster_names_detailed=="Tem GZMK+")
further_subset=further_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(further_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
further_subset_try=FindClusters2(further_subset, cluster.range=c(2,4), by=0.1, res=0.05, verbose=T)
DimPlot(further_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: terminal effector T (Tte)-like cells (also a subset of conventional CD8 Tmem), CD8T.Tem.GZMK_NKG7
VlnPlot(further_subset_try, features=c("GZMB","GNLY","NKG7","ZEB2","GZMA","GZMK"), pt.size=0) # 2
# cluster1: CCR7- GZMKhi, the conventional Tem cells (cytolytic and generally lacked the memory markers LEF1, CD27, CD28, and CD127), CD8T.Tem.GZMKhi
VlnPlot(further_subset_try, features=c("CCR7","GZMK","LEF1","CD27","CD28"), pt.size=0) # 0, 1, 3

Idents(further_subset_try)="seurat_clusters"
obj_annotation=further_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="0"), "CD8T.Tem.GZMKhi", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="1"), "CD8T.Tem.GZMKhi", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="2"), "CD8T.Tem.GZMK_NKG7", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="3"), "CD8T.Tem.GZMKhi", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if CD8 is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Reannotate NK
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="NK cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/nk_cells/nk_cells_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset Proliferative
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Proliferative")
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(2), by=0.1, res=0.1, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: less proliferating but highly active CD56dim CD16+ NK cells, NK.prolif.MCM
VlnPlot(CD4T_obj_subset_try, features=c("FCGR3A","NCAM1","IL2RB","SPON2","FCER1G","GZMB","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MKI67"), pt.size=0) # 0
# cluster2: highly proliferating CD56dim CD16+ NK cells, NK.prolif.MKI67
VlnPlot(CD4T_obj_subset_try, features=c("FCGR3A","NCAM1","IL2RB","SPON2","FCER1G","GZMB","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MKI67"), pt.size=0) # 1

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "NK.prolif.MCM", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "NK.prolif.MKI67", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Subset CD56dim
Tmem_celltypes=unique(CD4T_obj$Cluster_names_detailed)[!grepl("Proliferative|CD56bright",unique(CD4T_obj$Cluster_names_detailed))]
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% Tmem_celltypes)
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(5,7), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: cytotoxic, NK.CD56dim.FCER1G
VlnPlot(CD4T_obj_subset_try, features=c("MYOM2","FCER1G","PRF1","SPON2"), pt.size=0) # CD56dim CD57int, # 0
# cluster1: less activated but more proliferative, NK.CD56dim.HNRNPH1
VlnPlot(CD4T_obj_subset_try, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0) # No.
# cluster2: adaptive NK (different from adaptive-like because KLRC2+, which is a marker of human cytomegalovirus/HCMV infection), NK.CD56dim.FCER1Glo_KLRC2pos
VlnPlot(CD4T_obj_subset_try, features=c("KLRC2","CD3E","CCL5","VIM"), pt.size=0) # CD56dim CD57+, # 3
# cluster4: adaptive-like NK (because KLRC2-; this is possibly a subset derived from infections), NK.CD56dim.FCER1Glo_KLRC2neg
VlnPlot(CD4T_obj_subset_try, features=c("CD3E","VIM","IL32","CD52","KLRC2","FCER1G"), pt.size=0) # 1
# cluster3: transcriptionally active (this is a homeostatic activation; apart from below also express "NFKBIA","CD69","CXCR4","ZFP36"), NK.CD56dim.FOS
VlnPlot(CD4T_obj_subset_try, features=c("FOS","FOSB","JUN","JUNB"), pt.size=0) # CD56dim CD57low, # 4
# cluster5: IFN-I+ NK (reported to have been found in early severe COVID-19 NK cells, so it's related to infections), NK.CD56dim.ISG
VlnPlot(CD4T_obj_subset_try, features=c("MX1","ISG15","ISG20"), pt.size=0) # No.
# cluster1: transitional NK between CD56dim and CD56hi, NK.CD56hi.FCGR3A
VlnPlot(CD4T_obj_subset_try, features=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"), pt.size=0) # CD56dim CD57-, # 2
# cluster0: NKT
VlnPlot(CD4T_obj_subset_try, features=c("CD3D","CD3G","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "NK.CD56dim.FCER1G", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "NK.CD56dim.FCER1Glo_KLRC2neg", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "NK.CD56hi.FCGR3A", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "NK.CD56dim.FCER1Glo_KLRC2pos", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "NK.CD56dim.FOS", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Subset CD56hi
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="CD56bright")
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(2,4), by=0.1, res=0.3, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: CD56bright NK with responses to cytokines and regulation on cytokine production, NK.CD56hi.COTL1
VlnPlot(CD4T_obj_subset_try, features=c("COTL1","XCL1","LTB","GZMK"), pt.size=0) # CD56bright
# cluster1: transitional NK between CD56dim and CD56hi, NK.CD56hi.FCGR3A
VlnPlot(CD4T_obj_subset_try, features=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"), pt.size=0) # No.
# cluster2: CD56bright NK with chemokine-secreting property (may play a critical role in recruitment of T cells and other immune cells), NK.CD56hi.CCL3
VlnPlot(CD4T_obj_subset_try, features=c("CCL4L2","CCL3","CCL4"), pt.size=0) # No.
# cluster3: CD56bright NK which is S100B+ cytotoxic (in the CD56hiCD16-CD49a-KIR-NK cluster, one of subsets of conventional CD56hi NK), NK.CD56hi.S100B
VlnPlot(CD4T_obj_subset_try, features=c("S100B","MYOM2","BEX3","SOX4","ITGA1","KIR2DL4"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="Cluster_names_detailed"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="CD56bright"), "NK.CD56hi.COTL1", Cluster_numbers))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if CD8 is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Reannotate gdT
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="gd T cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/gd_t_cells/gd_t_cells_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Check gdT
CD4T_obj=NormalizeData(CD4T_obj)
Idents(CD4T_obj)="Cluster_names_detailed"
# OtherT.gdT.LEF1
VlnPlot(CD4T_obj, features=c("RTKN2", "TRDC", "TRGC2", "LEF1", "IKZF2", "SOX4", "ZNF331", "ARID5B", "NUCB2", "CRTAM"), pt.size=0) # gd naive
# OtherT.gdT.KLRC1
VlnPlot(CD4T_obj, features=c("TRDC", "TRGC1", "TRGV9", "TRDV2", "KLRD1", "IL7R", "KLRC1", "DUSP2", "GNLY", "KLRG1"), pt.size=0) # Vd2 GZMK+, Vd2 GZMB+
# OtherT.gdT.CMC1
VlnPlot(CD4T_obj, features=c("TRDC", "TIGIT", "KLRC2", "TRGC2", "IKZF2", "GCSAM", "FCRL6", "TRDV1", "CST7", "CMC1"), pt.size=0) # Vd1 GZMK+, Vd1 GZMB+

obj_annotation=CD4T_obj[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="gd naive"), "OtherT.gdT.LEF1", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="Vd2 GZMK+"), "OtherT.gdT.KLRC1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="Vd2 GZMB+"), "OtherT.gdT.KLRC1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="Vd1 GZMK+"), "OtherT.gdT.CMC1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="Vd1 GZMB+"), "OtherT.gdT.CMC1", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if gdT is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Reannotate Progenitors
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="Progenitor cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/progenitor_cells/progenitor_cells_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset Progenitors
CD4T_obj_subset=CD4T_obj %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(3,4), by=0.1, res=0.3, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: CLP, Other.progenitor.CLP
VlnPlot(CD4T_obj_subset_try, features=c("DNTT","PRSS2","EGFL7","CYTL1","CD34"), pt.size=0) # 1
# cluster1: HSC, Other.progenitor.HSC
VlnPlot(CD4T_obj_subset_try, features=c("SAMHD1","LINC00299","TESC","ITGB7","COTL1"), pt.size=0) # 3
# cluster0: previously reported as presumptive multi-lineage progenitor, may actually be MPP, Other.progenitor.MPP
VlnPlot(CD4T_obj_subset_try, features=c("EGFL7","CYTL1","CD34","CSF3R","SMIM24","C1QTNF4"), pt.size=0) # 0, 2
# cluster1: MEP (with genes reported to be Ery-lineage progenitor gene signatures/the first row; and those to be megakaryocytic-erythroid progenitors/the 2nd row), Other.progenitor.MEP
VlnPlot(CD4T_obj_subset_try, features=c("EPOR","KLF1","CSF2RB","APOC1","CNRIP1","CD44","TFRC","ITGA2B"), pt.size=0) # No.
# cluster2: NK cell progenitors (NKPs in short)
VlnPlot(CD4T_obj_subset_try, features=c("IL7R","ITGA4","CD27","CD44","IL2RB","CD7"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "Other.progenitor.MPP", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "Other.progenitor.CLP", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "Other.progenitor.MPP", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "Other.progenitor.HSC", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if Progenitors is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Reannotate Myeloid cells
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="Myeloid cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/myeloid_cells/myeloid_cells_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset mono
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% c("Classical monocytes","Non-classical monocytes"))
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures()
CD4T_obj_subset=SketchData(CD4T_obj_subset, ncells=50000, method="LeverageScore", sketched.assay="sketch")
DefaultAssay(CD4T_obj_subset)="sketch"
DimPlot(CD4T_obj_subset, group.by="Batch") # no batch effect, thus no need to integrate

CD4T_obj_subset=CD4T_obj_subset %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T) + NoLegend()

# cluster0: CD14 classical Mono, Mono.classical
VlnPlot(CD4T_obj_subset_try, features=c("S100A9","VCAN","LYZ","S100A12","TREM1","CSF3R","CD14","FCGR3A"), pt.size=0) # Classical monocytes, 0 2 3 4 6
# cluster1: CD14 classical Mono with some mature DC phenotypes, Mono.classical.HLA
VlnPlot(CD4T_obj_subset_try, features=c("CD83","CD14","CCL17","CCL13","HLA-DPA1","HLA-DQA1","HLA-DRB1"), pt.size=0) # 5
# cluster3: CD14++CD16+ Intermediate monocytes with cytotoxic properties, Mono.inter.GNLY
VlnPlot(CD4T_obj_subset_try, features=c("S100A9","VCAN","LYZ","FCGR3A","RHOC","KLRD1","KLRF1","GNLY"), pt.size=0) # No.
# cluster0: CD14-CD16+ nonclassical Mono, Mono.nonclassical
VlnPlot(CD4T_obj_subset_try, features=c("CD14","FCGR3A","LST1","AIF1","CTSS","CTSL"), pt.size=0) # Non-classical monocytes, 1
# cluster1: CD14-CD16+ nonclassical Mono with cell proliferation property (EIF5A and C1orf56)
VlnPlot(CD4T_obj_subset_try, features=c("CD14","FCGR3A","EIF5A","C1orf56"), pt.size=0) # No.
# cluster2: CD14-CD16+ nonclassical monocytes with cytotoxic properties, Mono.nonclassical.GNLY
VlnPlot(CD4T_obj_subset_try, features=c("CD14","FCGR3A","KLRD1","KLRF1","GNLY"), pt.size=0) # No.
# cluster4: CD14+CD16+ Intermediate monocytes with C1qhi as a proinflammatory cytokine-secreting and enhanced phagocytotic marker, Mono.inter.C1Q
VlnPlot(CD4T_obj_subset_try, features=c("CD14","FCGR3A","C1QA","C1QB","C1QC"), pt.size=0) # No.

CD4T_obj_subset_try=ProjectData(
  CD4T_obj_subset_try,
  assay="RNA",
  full.reduction="pca.full",
  sketched.assay="sketch",
  sketched.reduction="pca",
  umap.model="umap",
  dims=1:50,
  refdata=list(seurat_clusters_full="seurat_clusters")
)

DefaultAssay(CD4T_obj_subset_try)="RNA"
Idents(CD4T_obj_subset_try)="seurat_clusters_full"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "Mono.classical", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "Mono.nonclassical", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "Mono.classical", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "Mono.classical", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "Mono.classical", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="5"), "Mono.classical.HLA", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="6"), "Mono.classical", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Subset DCs
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% c("cDCs","pDCs"))
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: AXL+ mDC, DC.AXL_mDC
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","FCER1A","CLEC10A"), pt.size=0) # No.
# cluster1: cDC1, DC.cDC1
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","THBD","CLEC9A","CADM1","XCR1","BATF3"), pt.size=0) # 5
# cluster2: AXL+ pDC, DC.AXL_pDC
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","TPM2","LRRC26"), pt.size=0) # No.
# cDC2_1, DC.cDC2_1
VlnPlot(CD4T_obj_subset_try, features=c("FCER1A", "CD14", "CLEC10A", "CTSS", "ENHO", "CD1C", "MRC1", "FCGR2B", "PID1", "IL13RA1"), pt.size=0) # 1, cDCs
# cDC2_2, DC.cDC2_2
VlnPlot(CD4T_obj_subset_try, features=c("FCER1A", "BASP1", "CD1C", "CD74", "CLEC10A", "HLA-DPA1", "ENHO", "HLA-DPB1", "PLD4", "HLA-DQA1"), pt.size=0) # 2, cDCs
# cluster1: pDC, DC.pDC
VlnPlot(CD4T_obj_subset_try, features=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","CLEC4C"), pt.size=0) # 0 3 4, pDCs
# cluster5: cDC2 with cell proliferative property, DC.cDC2_2.MKI67
VlnPlot(CD4T_obj_subset_try, features=c("PCLAF","TYMS","MKI67"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "DC.pDC", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "DC.cDC2_1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "DC.cDC2_2", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "DC.pDC", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "DC.pDC", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="5"), "DC.cDC1", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if Mono/DC is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Reannotate B cells
#####################################
###

library(Seurat)
library(dplyr)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="B cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/b_cells/b_cells_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj) + NoLegend()

### Subset Naive
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% c("Naive","Naive-IFN"))
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, reduction="umap") + NoLegend()

# B naive, B.naive.kappa
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"), pt.size=0) # 0 2 3 4 5 7, Naive+Naive_IFN
# lambda B naive, B.naive.lambda
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","MS4A1","TCL1A","IL4R","IGLC6","IGLC7","IGLC3"), pt.size=0) # 6
# B memory
VlnPlot(CD4T_obj_subset_try, features=c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781"), pt.size=0) # No.
# cluster1: IgMhi transitional B cells with MZB (marginal zone B) genes which will later develop into MZB cells, B.trans?
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","PLD4","MZB1"), pt.size=0) # No.
# B intermediate
VlnPlot(CD4T_obj_subset_try, features=c("MS4A1","TNFRSF13B","IGHM","IGHD","AIM2","CD79A","LINC01857","RALGPS2","BANK1","CD79B"), pt.size=0)
# cluster2: B intermediate (expressing naive markers - IL4R+TCL1A+, memory markers - CD24+CD27+, unswithced - IgMhi, switched - IgG+IgA+), B.inter.IL4R
VlnPlot(CD4T_obj_subset_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","IGHG3"), pt.size=0) # No.
# cluster1: B cells at the intermediate status between naive (IGHD+IL4R+TCL1A+) and Bmem (CD27+IGHA1+IGHG2+), slightly closer to Bnaive compared to Cluster2, B.inter.lambda.IL4R
VlnPlot(CD4T_obj_subset_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"), pt.size=0) # No.
# cluster2: B cells at the intermediate status between naive (IGHD+TCL1A+) and Bmem (CD27+IGHA1+IGHG2+CD24+) but closer to Bmem, B.inter.lambda
VlnPlot(CD4T_obj_subset_try, features=c("IGHD","IL4R","TCL1A","IGHM","CD27","IGHA1","IGHG2","CD24"), pt.size=0) # No.
# cluster2: B intermediate cells transforming from B naive to IgM+IgD+ Bmem (closer to IgM+IgD+ Bmem since CD27+ IGHMhi IGHD+ PLD4+ MZB1+), B.inter.IGHM
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","CD27","PLD4","MZB1","CD1C"), pt.size=0) # No.
# cluster0: B intermediate cells transforming from B naive to IgM+IgD+CD5-CD27- Bmem (close to Bmem since IL4R-TCL1A-CD27-IGHM+IGHD+), B.inter.CD27neg
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R","CD5"), pt.size=0) # 1

Idents(CD4T_obj_subset_try)="seurat_clusters"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "B.naive.kappa", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "B.inter.CD27neg", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="5"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="6"), "B.naive.lambda", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="7"), "B.naive.kappa", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Subset Plasma cells
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed=="Plasma cells")
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(2), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, reduction="umap") + NoLegend()

# Plasmablast, B.pb
VlnPlot(CD4T_obj_subset_try, features=c("CXCR3","HLA-DRA","KLF4","MKI67","MS4A1"), pt.size=0) # 1
# Plasma, B.pc
VlnPlot(CD4T_obj_subset_try, features=c("PTPRC","CXCR4","IGHA1","IGHG1","IGHA2","IGHG2","IGHG3","IGHM","TNFRSF17","SDC1","JCHAIN"), pt.size=0) # 0

Idents(CD4T_obj_subset_try)="seurat_clusters"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "B.pc", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "B.pb", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

### Check the rest
Tmem_celltypes=unique(CD4T_obj$Cluster_names_detailed)[!grepl("Naive|Plasma",unique(CD4T_obj$Cluster_names_detailed))]
CD4T_obj_subset=subset(CD4T_obj, Cluster_names_detailed %in% Tmem_celltypes)
CD4T_obj_subset=CD4T_obj_subset %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30, return.model=T)
DimPlot(CD4T_obj_subset, reduction="umap")
source("~/Rscripts/FindCluster2_Functions.R")
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(2), by=0.1, res=0.2, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, reduction="umap") + NoLegend()

Idents(CD4T_obj_subset)="Cluster_names_detailed"
# B memory, B.mem.IGHA
VlnPlot(CD4T_obj_subset, features=c("MS4A1","COCH","AIM2","BANK1","SSPN","CD79A","TEX9","RALGPS2","TNFRSF13C","LINC01781","IGLC2","IGLC3","IGHG3","IGKC","IGHG1","IGHG2","IGHE","IGHA1","IGHA2","IGHD","IGHM"), pt.size=0) # Switched memory
# cluster1: IgMhi transitional B cells with MZB (marginal zone B) genes which will later develop into MZB cells, B.trans
VlnPlot(CD4T_obj_subset, features=c("IGHM","PLD4","MZB1"), pt.size=0) # Transitional
# cluster0: B intermediate cells transforming from B naive to classical Bmem (close to classical Bmem since IL4R-TCL1A-CD27dim IGHMdim IGHDdim), B.inter.IGHM
VlnPlot(CD4T_obj_subset, features=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R"), pt.size=0) # CD5+ B cells
# DN Bmem - CD19+, MS4A1+, CD27-, IGHD-, IGHMlo, B.mem.CD27neg
VlnPlot(CD4T_obj_subset, features=c("CD19","MS4A1","CR2","CD27","IGHM","IGHD","MME","CD38", "POU2AF1","PAX5"), pt.size=0) # Atypical memory
# IgM Memory (CD19+IgD+CD27+), B.mem.IGHM
VlnPlot(CD4T_obj_subset, features=c("TCL1A","CD19","IGHD","CD27","MME","IGHM"), pt.size=0) # Non-switched memory
# Regulatory B - MS4A1neg
VlnPlot(CD4T_obj_subset, features=c("CD1D","CD5","CD19","MS4A1","CR2","CD24","CD38","CD40"), pt.size=0) # No.
# B naive, B.naive.kappa
VlnPlot(CD4T_obj_subset, features=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"), pt.size=0) # Activated

obj_annotation=CD4T_obj_subset[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Non-switched memory"), "B.mem.IGHM", Cluster_names_detailed)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Switched memory"), "B.mem.IGHA", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Transitional"), "B.trans", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Atypical memory"), "B.mem.CD27neg", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="Activated"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset, idents="CD5+ B cells"), "B.inter.IGHM", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# check if B is done
table(colnames(CD4T_obj) %in% OBJ_Annotation$cell)

#####################################



### Add DN T cell and MAIT cell annotation
#####################################
###

library(Seurat)
library(dplyr)

### Load
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="DN T cells" | Cluster_names=="MAIT cells")

Idents(CD4T_obj)="Cluster_names_detailed"
obj_annotation=CD4T_obj[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="MAIT cells"), "OtherT.MAIT", Cluster_names_detailed)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj, idents="DN T cells"), "OtherT.dnT", annotation))
table(obj_annotation$annotation, useNA="ifany")

OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=data.table::rbindlist(list(OBJ_Annotation, obj_annotation %>% select(cell, annotation)))
data.table::fwrite(OBJ_Annotation, "~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")

# Check if all celltypes are done
table(colnames(Seurat_whole_obj) %in% OBJ_Annotation$cell)

#####################################



### Add new annotation
#####################################
###

library(Seurat)
library(dplyr)

Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
OBJ_Annotation=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz", sep="\t")
OBJ_Annotation=OBJ_Annotation[match(colnames(Seurat_whole_obj), OBJ_Annotation$cell),]
# check
table(OBJ_Annotation$cell==colnames(Seurat_whole_obj))
Seurat_whole_obj$Annot.detailed=OBJ_Annotation$annotation

### Add the metadata of celltype annotation
celltypes_detailed=unique(OBJ_Annotation$annotation)
OBJ_Annotation=OBJ_Annotation %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", annotation)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(annotation %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter))
table(OBJ_Annotation$Annot.inter, useNA="ifany")
OBJ_Annotation=OBJ_Annotation %>%
  mutate(Annot.rough=ifelse(annotation %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", annotation)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
table(OBJ_Annotation$Annot.rough, useNA="ifany")

table(OBJ_Annotation$cell==colnames(Seurat_whole_obj)) # check
Seurat_whole_obj=AddMetaData(Seurat_whole_obj, metadata=OBJ_Annotation$Annot.inter, col.name="Annot.inter")
Seurat_whole_obj=AddMetaData(Seurat_whole_obj, metadata=OBJ_Annotation$Annot.rough, col.name="Annot.rough")
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")

# add the agecut metadata
df=Seurat_whole_obj[[]] %>%
  mutate(agecut=ifelse(Age %in% c(25:30), "25~30",
                       ifelse(Age %in% c(31:47), "31~47",
                              ifelse(Age %in% c(48:68), "48~68",
                                     ifelse(Age %in% c(69:81), "69~81", NA)))))
table(df$agecut, useNA='ifany')
table(rownames(df)==colnames(Seurat_whole_obj)) # check
Seurat_whole_obj=AddMetaData(Seurat_whole_obj, metadata=df$agecut, col.name="agecut")
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")

### Save metadata for RawCountObj_OneSamplePerDonor.rds
meta_data=Seurat_whole_obj[[]] %>% tibble::rownames_to_column("cell_id")
data.table::fwrite(meta_data, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
unlink("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz")

#####################################
