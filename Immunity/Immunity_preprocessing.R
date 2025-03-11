
### Prepare raw counts and meta of 10.1016/j.immuni.2023.10.013 (the Immunity dataset)
#####################################
###

library(Seurat)
library(dplyr)

### Download from Synapse:syn50542388
# b_cells.tar.gz
# cd4_t_cells.tar.tgz
# cd4_t_helper_memory_cells.tar.tgz (the subset of cd4_t_cells.tar.tgz, so skip this one)
# conventional_cd8_t_cells.tar.tgz
# gd_t_cells.tar.gz
# mait_cells.tar.tgz
# myeloid_cells.tar.gz
# nk_cells.tar.gz
# progenitor_cells.tar.tgz

#!: compared with the metadata in all_pbmcs.tar.gz, 1490 DN T cells are missing.

#!: 166 donors, splitted to 317 tubes, ran in 14 batches, generate 106 files
#!: various ages in each donor??? We use tube_id then
#!: each tube_id goes into only 1 batch, and by checking DimPlot(group.by="Batch"), limited Batch-effect was found, so skip harmonize

#####################################



# ### Modify the Cluster_names to fit Annot.rough/inter/detailed (THIS IS NOT USEFUL, DISCARD)
# #####################################
# ###
#
# ### Load the total meta (no detailed annotation here)
# raw_meta=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/all_pbmcs_metadata.csv", sep=",", row.names="X")
# raw_meta_toupdate=raw_meta
# raw_meta_toupdate$cell_id=rownames(raw_meta_toupdate)
# colnames(raw_meta_toupdate)
#
# ### Load the meta from each celltype (containing detailed anno)
# meta_perCelltype=list()
# folders=list.dirs("/home/jo_79/Project_PBMCage/Immunity/all_pbmcs")
# folders=folders[2:9]
# for (i in 1:length(folders)) {
#   files=list.files(folders[i])
#   meta.path=paste0(folders[i], "/", files[grepl("metadata\\.csv$",files)])
#   meta=read.csv(meta.path, row.names='X')
#   meta$cell_id=rownames(meta)
#
#   meta_from_raw=raw_meta[rownames(meta),] # because some metadata are incorrect, so here we directly extract from the original all_pbmcs_metadata
#   meta_from_raw$cell_id=rownames(meta_from_raw)
#
#   meta_merged=full_join(meta_from_raw %>% select(-c(Cluster_numbers, Cluster_names)), meta, by="cell_id")
#   shared_cols=gsub("\\.x|\\.y","",colnames(meta_merged))
#   shared_cols=shared_cols[duplicated(shared_cols)]
#   for (j in 1:length(shared_cols)) {
#     if (class(meta_merged[[paste0(shared_cols[j],".y")]])=="numeric") {
#       logic_=unique(sapply(meta_merged[[paste0(shared_cols[j],".x")]], function(x) signif(x,3))==
#                       sapply(meta_merged[[paste0(shared_cols[j],".y")]], function(x) signif(x,3)))
#     } else {
#       logic_=unique(meta_merged[[paste0(shared_cols[j],".x")]]==meta_merged[[paste0(shared_cols[j],".y")]])
#     }
#
#     if (length(logic_)==1) {
#       if (logic_==TRUE) {
#         meta_merged=select(meta_merged, -paste0(shared_cols[j],".x"))
#       }
#     }
#   }
#
#   meta_merged=meta_merged[,colSums(is.na(meta_merged))!=nrow(meta_merged)]
#   colnames(meta_merged)=gsub("\\.y$","",colnames(meta_merged))
#   meta_merged$Annot.rough.tempt=gsub(".*_pbmcs/","",folders[i])
#   meta_perCelltype[[i]]=meta_merged
# }
#
# ### Check
# for (i in 1:length(meta_perCelltype)) {
#   meta_merged=meta_perCelltype[[i]]
#   print(colnames(meta_merged)[grepl("\\.x",colnames(meta_merged))])
# }
# tempt_meta=meta_perCelltype[[3]]
# head(tempt_meta[tempt_meta$log2_nCount.x!=tempt_meta$log2_nCount,])
# meta_perCelltype[[3]]=meta_perCelltype[[3]] %>% select(-c(log2_nCount.x, log2_mt.x))
#
# tempt_meta=meta_perCelltype[[5]]
# head(tempt_meta[tempt_meta$Age.x!=tempt_meta$Age,])
# meta_perCelltype[[5]]=meta_perCelltype[[5]] %>% select(-Age) %>% dplyr::rename(Age=Age.x)
#
# for (i in 1:length(meta_perCelltype)) {
#   meta_merged=meta_perCelltype[[i]]
#   if ("Gender" %in% colnames(meta_merged)) print(i)
# }
# tempt_meta=meta_perCelltype[[5]]
# table(tempt_meta$Sex==tempt_meta$Gender)
# meta_perCelltype[[5]]=meta_perCelltype[[5]] %>% select(-Gender)
#
# ### Merge
# meta_perCelltype_all=data.table::rbindlist(meta_perCelltype, fill=TRUE)
#
# ### Arrange subpopulation names
# merged_meta=meta_perCelltype_all %>%
#   mutate(Cluster_names=ifelse(is.na(Cluster_names), Annot.rough.tempt, Cluster_names))
# table(merged_meta$Cluster_names, useNA="ifany") # check
# merged_meta$celltype=paste0(merged_meta$Annot.rough.tempt, ",", merged_meta$Cluster_names)
# # Annotate the detailed celltypes according to the celltypes in the original paper
# table(merged_meta$celltype, useNA="ifany")
# merged_meta$celltype=as.factor(merged_meta$celltype)
# levels(merged_meta$celltype)
# levels(merged_meta$celltype)=c("B.inter.IGHM","B.mem.CD27neg","B.naive.CD5pos",
#                                "B.naive","B.naive.ISG","B.mem.IGHM",
#                                "B.pc","B.mem.IGHG","B.trans",
#                                "CD4T.Tem.LAG3","CD4T.Tem.HLA","CD4T.naive",
#                                "CD4T.naive.ISG","CD4T.Temra","CD4T.CTL",
#                                "CD4T.Th.Tfh","CD4T.Th.Th1","CD4T.Th.Th1_Th17",
#                                "CD4T.Th.Th17","CD4T.Th.Th2","CD4T.Th.Th22",
#                                "CD4T.Treg.GZMB","CD4T.Treg.KLRB1","CD4T.Treg.HLA",
#                                "CD4T.Treg","CD8T.prolif.HLA","CD8T.naive",
#                                "CD8T.naive.ISG","CD8T.Tem.GZMB_FCGR3A","CD8T.prolif",
#                                "CD8T.Tcm.CCR4neg","CD8T.Tcm.CCR4pos","CD8T.Tem.GZMB_HNRNPH1",
#                                "CD8T.Tem.GZMKhi","CD8T.Temra","CD8T.Tem.GZMK_XCL1",
#                                "CD8T.Trm","OtherT.gdT.LEF1","OtherT.gdT.GZMB_Vd1",
#                                "OtherT.gdT.GZMK_Vd1","OtherT.gdT.GZMB_Vd2","OtherT.gdT.GZMK_Vd2",
#                                "OtherT.MAIT","DC.cDC2_2","Mono.classical",
#                                "Mono.nonclassical","DC.pDC","NK.CD56hi.FCGR3A",
#                                "NK.CD56dim.FCER1Glo_KLRC2neg","NK.CD56dim.FCER1Glo_KLRC2pos","NK.CD56dim.HNRNPH1",
#                                "NK.CD56dim.FOS","NK.prolif.MKI67","Other.progenitor")
# merged_meta$Annot.detailed=as.character(merged_meta$celltype)
# # generate the Annot.inter and Annot.rough accordingly
# detailed_celltypes=levels(merged_meta$celltype)
# merged_meta=merged_meta %>%
#   mutate(Annot.inter=ifelse(Annot.detailed %in% detailed_celltypes[grepl("^B\\.inter",detailed_celltypes)], "B.inter", Annot.detailed)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^B\\.mem",detailed_celltypes)], "B.mem", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^B\\.naive",detailed_celltypes)], "B.naive", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^B\\.pc|pb",detailed_celltypes)], "B.pbpc", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD4T\\.Tem",detailed_celltypes)], "CD4T.Tem", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD4T\\.naive",detailed_celltypes)], "CD4T.naive", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD4T\\.Th",detailed_celltypes)], "CD4T.Th", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD4T\\.Treg",detailed_celltypes)], "CD4T.Treg", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD8T\\.prolif",detailed_celltypes)], "CD8T.prolif", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD8T\\.naive",detailed_celltypes)], "CD8T.naive", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD8T\\.Tem",detailed_celltypes)], "CD8T.Tem", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^CD8T\\.Tcm",detailed_celltypes)], "CD8T.Tcm", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^OtherT\\.gdT",detailed_celltypes)], "OtherT.gdT", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^DC\\.cDC2",detailed_celltypes)], "DC.cDC2", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^NK\\.CD56hi",detailed_celltypes)], "NK.CD56hi", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^NK\\.CD56dim",detailed_celltypes)], "NK.CD56dim", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^NK\\.prolif",detailed_celltypes)], "NK.prolif", Annot.inter)) %>%
#   mutate(Annot.inter=ifelse(Annot.inter %in% detailed_celltypes[grepl("^Other\\.progenitor",detailed_celltypes)], "Other.progenitor", Annot.inter))
# table(merged_meta$Annot.inter, useNA="ifany")
# merged_meta=merged_meta %>%
#   mutate(Annot.rough=ifelse(Annot.detailed %in% detailed_celltypes[grepl("^B\\.",detailed_celltypes)], "B cells", Annot.detailed)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^CD4T\\.",detailed_celltypes)], "CD4T cells", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^CD8T\\.",detailed_celltypes)], "CD8T cells", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^DC\\.",detailed_celltypes)], "DCs", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^Mono\\.",detailed_celltypes)], "Monocytes", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^NK\\.",detailed_celltypes)], "NK cells", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^Other\\.",detailed_celltypes)], "Other cells", Annot.rough)) %>%
#   mutate(Annot.rough=ifelse(Annot.rough %in% detailed_celltypes[grepl("^OtherT\\.",detailed_celltypes)], "OtherT", Annot.rough))
# table(merged_meta$Annot.rough, useNA="ifany")
#
# ### Save the annotated metadata
# merged_meta=merged_meta %>% select(-c(celltype, Annot.rough.tempt))
# data.table::fwrite(merged_meta, "~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")
#
# #####################################
#
#
#
# ### Reconstruct the pbmc with the files above (THIS IS NOT USEFUL, DISCARD)
# #####################################
# ###
#
# ### Load meta
# merged_meta=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")
#
# PseudoList=list()
#
# folders=list.dirs("/home/jo_79/Project_PBMCage/Immunity/all_pbmcs")
# folders=folders[2:9]
# for (i in 1:length(folders)) {
#   files=list.files(folders[i])
#   rna.path=paste0(folders[i], "/", files[grepl("rna\\.rds$",files)])
#   meta.path=paste0(folders[i], "/", files[grepl("metadata\\.csv$",files)])
#
#   norm_counts=readRDS(rna.path)
#   meta=merged_meta %>% subset(cell_id %in% colnames(norm_counts))
#   meta=meta[match(colnames(norm_counts), meta$cell_id),]
#   meta=meta %>% tibble::column_to_rownames("cell_id") %>%
#     dplyr::rename(percent.rb=percent.ribo, donor_id=Donor_id, sex=Sex, age=Age) %>%
#     mutate(sex=ifelse(sex=="Male","M","F")) %>%
#     mutate(agecut=ifelse(age %in% c(25:30), "25~30",
#                          ifelse(age %in% c(31:47),"31~47",
#                                 ifelse(age %in% c(48:68), "48~68",
#                                        ifelse(age %in% c(69:81), "69~81", NA)))))
#   print(table(meta$age, useNA="ifany")) # check
#   seuratobj=CreateSeuratObject(counts=norm_counts, meta.data=meta)
#
#   # check if each tube_id goes into only 1 batch
#   tempt1=length(unique(seuratobj$Batch))
#   tempt2=length(unique(seuratobj$Tube_id))
#   tempt3=length(unique(paste0(seuratobj$Batch, "_", seuratobj$Tube_id)))
#   print(tempt1); print(tempt2); print(tempt3)
#   if (tempt3!=tempt2) {print("Batch > Tube_id"); break}
#
#   # get objs on each Annot.detailed to avoid memory overloading
#   detailed_celltypes=unique(seuratobj$Annot.detailed)
#   for (j in 1:length(detailed_celltypes)) {
#     detailed_obj=subset(seuratobj, Annot.detailed==detailed_celltypes[j])
#
#     # create pseudobulk on Annot.detailed
#     pseudo_detailed_=AggregateExpression(detailed_obj, assays="RNA", return.seurat=T, group.by=c("Batch","Tube_id"))
#     pseudo_detailed_$Batch=pseudo_detailed_$orig.ident
#     pseudo_detailed_$Tube_id=gsub(".*_","",colnames(pseudo_detailed_))
#     sex_match=data.frame(Tube_id=detailed_obj$Tube_id,
#                          donor_id=detailed_obj$donor_id,
#                          age=detailed_obj$age,
#                          sex=detailed_obj$sex,
#                          agecut=detailed_obj$agecut,
#                          Age_group=detailed_obj$Age_group) %>%
#       subset(!duplicated(.))
#
#     pseudo_detailed_$donor_id=sex_match$donor_id[match(pseudo_detailed_$Tube_id, sex_match$Tube_id)]
#     pseudo_detailed_$sex=sex_match$sex[match(pseudo_detailed_$Tube_id, sex_match$Tube_id)]
#     pseudo_detailed_$age=sex_match$age[match(pseudo_detailed_$Tube_id, sex_match$Tube_id)]
#     pseudo_detailed_$agecut=sex_match$agecut[match(pseudo_detailed_$Tube_id, sex_match$Tube_id)]
#     pseudo_detailed_$Age_group=sex_match$Age_group[match(pseudo_detailed_$Tube_id, sex_match$Tube_id)]
#     pseudo_detailed_$Annot.detailed=detailed_celltypes[j]
#     pseudo_detailed_$Annot.inter=unique(subset(seuratobj[[]], Annot.detailed==detailed_celltypes[j])$Annot.inter)
#     pseudo_detailed_$Annot.rough=unique(subset(seuratobj[[]], Annot.detailed==detailed_celltypes[j])$Annot.rough)
#     PseudoList=c(PseudoList, list(pseudo_detailed_))
#   }
# }
#
# # join all the pseudobulk_detailed
# pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=paste0("D",c(1:length(PseudoList))), merge.data=FALSE, project="Imm_pseudobulk_detailed")
# pbmc.merged=JoinLayers(pbmc.merged)
#
# # save the pseudobulk_detailed
# saveRDS(pbmc.merged, "~/Project_PBMCage/Immunity/all_pbmcs/Pseudobulk_detailed.rds")
#
# #####################################
#
#
#
# ### Generate pseudobulk_inter and pseudobulk_rough (THIS IS NOT USEFUL, DISCARD)
# #####################################
# ###
#
# ### Pseudobulk_inter
# pbmc.merged=readRDS("~/Project_PBMCage/Immunity/all_pbmcs/Pseudobulk_detailed.rds")
# match_df=pbmc.merged[[]] %>% select(-Annot.detailed) %>% subset(!duplicated(.))
# pbmc.merged_inter=AggregateExpression(pbmc.merged, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id","Annot.inter"))
#
# age=as.numeric(unlist(lapply(strsplit(colnames(pbmc.merged_inter), split="_"), function(x) x[1])))
# sex=unlist(lapply(strsplit(colnames(pbmc.merged_inter), split="_"), function(x) x[2]))
# donor_id=unlist(lapply(strsplit(colnames(pbmc.merged_inter), split="_"), function(x) x[3]))
# Annot.inter=unlist(lapply(strsplit(colnames(pbmc.merged_inter), split="_"), function(x) x[4]))
# Batch=match_df$Batch[match(donor_id, match_df$donor_id)]
# Annot.rough=match_df$Annot.rough[match(Annot.inter, match_df$Annot.inter)]
#
# pbmc.merged_inter$Batch=Batch
# pbmc.merged_inter$donor_id=donor_id
# pbmc.merged_inter$sex=sex
# pbmc.merged_inter$age=age
# pbmc.merged_inter$Annot.inter=Annot.inter
# pbmc.merged_inter$Annot.rough=Annot.rough
#
# saveRDS(pbmc.merged_inter, "~/Project_PBMCage/Immunity/all_pbmcs/Pseudobulk_inter.rds")
#
# ### Pseudobulk_rough
# pbmc.merged=readRDS("~/Project_PBMCage/Immunity/all_pbmcs/Pseudobulk_inter.rds")
# match_df=pbmc.merged[[]] %>% select(-Annot.inter) %>% subset(!duplicated(.))
# pbmc.merged_rough=AggregateExpression(pbmc.merged, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id","Annot.rough"))
#
# age=as.numeric(unlist(lapply(strsplit(colnames(pbmc.merged_rough), split="_"), function(x) x[1])))
# sex=unlist(lapply(strsplit(colnames(pbmc.merged_rough), split="_"), function(x) x[2]))
# donor_id=unlist(lapply(strsplit(colnames(pbmc.merged_rough), split="_"), function(x) x[3]))
# Annot.rough=unlist(lapply(strsplit(colnames(pbmc.merged_rough), split="_"), function(x) x[4]))
# Batch=match_df$Batch[match(donor_id, match_df$donor_id)]
#
# pbmc.merged_rough$Batch=Batch
# pbmc.merged_rough$donor_id=donor_id
# pbmc.merged_rough$sex=sex
# pbmc.merged_rough$age=age
# pbmc.merged_rough$Annot.rough=Annot.rough
#
# saveRDS(pbmc.merged_rough, "~/Project_PBMCage/Immunity/all_pbmcs/Pseudobulk_rough.rds")
#
# #####################################



### Reconstruct seurat objs with raw counts
#####################################
###

#... because the seurat objs given are with log data instead of raw counts, we cannot do pseudobulk analysis
#... here we reconstruct the objs with raw_counts_h5ad.tar.gz and all_pbmcs.tar.gz/all_pbmcs_metadata.csv
#... split the h5ad file based on Tube_id and convert to seurat, as so many cells in one file will lead to error

### Convert to seurat
# # in Python
# import scanpy as sc
# import anndata as ad
# import pandas as pd
#
# data = sc.read_h5ad('/home/jo_79/Project_PBMCage/Immunity/raw_counts_h5ad/pbmc_gex_raw_with_var_obs.h5ad')
# meta = pd.read_csv('/home/jo_79/Project_PBMCage/Immunity/all_pbmcs/all_pbmcs_metadata.csv', index_col=0)
# meta=meta.drop(columns=['nCount_RNA', 'nFeature_RNA', 'Cluster_names'])
# data.obs = data.obs.join(meta)
#
# tubes=set(meta['Tube_id'])
# for tube_id in tubes:
#   print(tube_id)
#   x=data[data.obs["Tube_id"] == tube_id]
#   del x.uns
#   del x.obsm
#   results_file = "/home/jo_79/Project_PBMCage/Immunity/raw_counts_h5ad/TubeID_" + str(tube_id) + ".h5ad"
#   x.write(results_file)

library(Seurat)
library(SeuratDisk)
files=list.files("~/Project_PBMCage/Immunity/raw_counts_h5ad/")
files=files[!grepl("pbmc_gex_raw",files)]
files=files[grepl("\\.h5ad",files)]

for (i in 1:length(files)) {
  Convert(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]), "h5seurat",
          overwrite=TRUE, assay="RNA")
  unlink(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
  seurat.obj=LoadH5Seurat(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",gsub("h5ad","h5seurat",files[i])), meta.data=FALSE, misc=FALSE)
  saveRDS(seurat.obj, paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",gsub("h5ad","rds",files[i])))
  unlink(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",gsub("h5ad","h5seurat",files[i])))

  print(paste0("=== Done with ",files[i]))
}

### Add metadata
library(dplyr)
raw_meta=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/all_pbmcs_metadata.csv", sep=",", row.names="X")
files=list.files("~/Project_PBMCage/Immunity/raw_counts_h5ad/")
files=files[grepl("\\.rds",files)]
for (i in 1:length(files)) {
  seurat.obj=readRDS(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
  seurat.obj[[]]=raw_meta[match(colnames(seurat.obj), rownames(raw_meta)),]
  saveRDS(seurat.obj, paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
}

#####################################



### Add the detailed clustering annotation from the original paper
#####################################
###

### Fetch all the detailed annotations
metadata_cluster=list()
folders=list.dirs("~/Project_PBMCage/Immunity/all_pbmcs")
folders=folders[grepl("_cell",folders)]
folders_correspondingCellAnnotation=c("B cells","CD4+ T cells","TRAV1-2- CD8+ T cells","gd T cells","MAIT cells","Myeloid cells","NK cells","Progenitor cells")
for (i in 1:length(folders)) {
  meta_files=list.files(folders[i])
  meta_file=meta_files[grepl("metadata\\.",meta_files)]
  metadata_csv=read.csv(paste0(folders[i],"/",meta_file), sep=",", row.names="X")
  metadata_cluster[[i]]=metadata_csv %>% tibble::rownames_to_column("cell_id") %>% select(any_of(c("cell_id","Cluster_names")))
  if (!("Cluster_names" %in% colnames(metadata_cluster[[i]]))) {
    metadata_cluster[[i]]=metadata_cluster[[i]] %>% mutate(Cluster_names=folders_correspondingCellAnnotation[i])
  }
}
metadata_cluster_all=data.table::rbindlist(metadata_cluster)
metadata_cluster_all$cell_id[duplicated(metadata_cluster_all$cell_id)] # check
# add 1490 DN T metadata
raw_meta=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/all_pbmcs_metadata.csv", sep=",", row.names="X")
DNT_cellid=raw_meta[!(rownames(raw_meta) %in% metadata_cluster_all$cell_id),]
DNT_meta=data.frame(cell_id=rownames(DNT_cellid), Cluster_names="DN T cells")
metadata_cluster_all=rbind(metadata_cluster_all, DNT_meta)
raw_meta$Cluster_names_detailed=metadata_cluster_all$Cluster_names[match(rownames(raw_meta), metadata_cluster_all$cell_id)]
raw_meta=raw_meta %>% tibble::rownames_to_column("cell_id")
data.table::fwrite(raw_meta, "~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")

### Add metadata detailed annotation
library(dplyr)
raw_meta=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")
files=list.files("~/Project_PBMCage/Immunity/raw_counts_h5ad/")
files=files[grepl("\\.rds",files)]
for (i in 1:length(files)) {
  seurat.obj=readRDS(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
  seurat.obj$Cluster_names_detailed=raw_meta$Cluster_names_detailed[match(colnames(seurat.obj), raw_meta$cell_id)]
  saveRDS(seurat.obj, paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
}

#####################################



### Take the representative Tube_ids to slim the whole obj
#####################################
###

###
raw_meta=data.table::fread("~/Project_PBMCage/Immunity/all_pbmcs/Metadata_annotated.csv.gz", sep="\t")
tempt=data.frame(raw_meta$Age, raw_meta$Donor_id) %>%
  subset(!duplicated(.))
split(tempt$raw_meta.Age, tempt$raw_meta.Donor_id) # seems that a donor visited 1-3 times for blood collection at different ages
age_donor_distribution=split(tempt$raw_meta.Age, tempt$raw_meta.Donor_id) %>% unlist() %>% unname()
ggplot2::ggplot(as.data.frame(age_donor_distribution), ggplot2::aes(x=age_donor_distribution)) +
  ggplot2::geom_density() # seems that the ages at visit of a donor distributed normally, thus randomly choose 1 age from a donor
# select the TubeIDs
tempt=data.frame(raw_meta$Age, raw_meta$Tube_id, raw_meta$Donor_id) %>%
  subset(!duplicated(.))
split(tempt$raw_meta.Tube_id, tempt$raw_meta.Donor_id)
set.seed(20240602)
TubeID_perDonor=split(tempt$raw_meta.Tube_id, tempt$raw_meta.Donor_id)
TubeID_selected=lapply(TubeID_perDonor, function(x) if (length(x)>1) {x[sample(1:length(x), size=1)]} else {x})
TubeID_selected_unlist=unlist(TubeID_selected) %>% unname()

### Take the selected objs to merge
Seurat_OBJs=list()
files=list.files("~/Project_PBMCage/Immunity/raw_counts_h5ad/")
files=files[grepl("\\.rds",files)]
files=files[gsub("TubeID_|\\.rds","",files) %in% TubeID_selected_unlist]
for (i in 1:length(files)) {
  Seurat_OBJs[[i]]=readRDS(paste0("~/Project_PBMCage/Immunity/raw_counts_h5ad/",files[i]))
}
# truncate to speed up merging
chuncks=ceiling(length(files)/10)
Seurat_OBJs_CHUNKs=list()
for (j in 1:chuncks) {
  CHUNK_=c(1:chuncks)[j]
  if (j<17) {
    Seurat_OBJs_chunk=Seurat_OBJs[(10*(CHUNK_-1)+1):(10*(CHUNK_-1)+10)]
    Seurat_OBJs_CHUNKs[[j]]=merge(Seurat_OBJs_chunk[[1]], Seurat_OBJs_chunk[-1], add.cell.ids=NULL)
  } else {
    Seurat_OBJs_chunk=Seurat_OBJs[(10*(CHUNK_-1)+1):length(Seurat_OBJs)]
    Seurat_OBJs_CHUNKs[[j]]=merge(Seurat_OBJs_chunk[[1]], Seurat_OBJs_chunk[-1], add.cell.ids=NULL)
  }
  print(paste0("CHUNK: ", j, " DONE!"))
}
Seurat_whole_obj=merge(Seurat_OBJs_CHUNKs[[1]], Seurat_OBJs_CHUNKs[-1], add.cell.ids=NULL)
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")

#####################################



### Reannotate CD4T
#####################################
###

library(Seurat)

### Add UMAP coordinates
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
CD4T_obj=subset(Seurat_whole_obj, Cluster_names=="CD4+ T cells")
CD4T_umap=read.csv("~/Project_PBMCage/Immunity/all_pbmcs/cd4_t_cells/cd4_umap.csv", sep=",", row.names="X", header=T)
CD4T_umap=CD4T_umap %>% tibble::rownames_to_column("cell_id") %>% subset(cell_id %in% colnames(CD4T_obj)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("cell_id")
UMAP_coordinates=as(CD4T_umap, "matrix")
UMAP_coordinates=UMAP_coordinates[match(colnames(CD4T_obj),rownames(UMAP_coordinates)),]
CD4T_obj[["UMAP"]]=CreateDimReducObject(embeddings=UMAP_coordinates, key="UMAP_", global=T, assay="RNA")
DimPlot(CD4T_obj_subset) + NoLegend()

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
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.3, verbose=T)
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
VlnPlot(CD4T_obj_subset, features=c("GZMB","GNLY","NKG7","ZEB2","GZMA","GZMK"), pt.size=0) # No.
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
further_subset_try=FindClusters2(further_subset, cluster.range=c(2,3), by=0.1, res=0.05, verbose=T)
DimPlot(further_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: terminal effector T (Tte)-like cells (also a subset of conventional CD8 Tmem), CD8T.Tem.GZMK_NKG7
VlnPlot(further_subset_try, features=c("GZMB","GNLY","NKG7","ZEB2","GZMA","GZMK"), pt.size=0) # 1
# cluster1: CCR7- GZMKhi, the conventional Tem cells (cytolytic and generally lacked the memory markers LEF1, CD27, CD28, and CD127), CD8T.Tem.GZMKhi
VlnPlot(further_subset_try, features=c("CCR7","GZMK","LEF1","CD27","CD28"), pt.size=0) # 0, 2

Idents(further_subset_try)="seurat_clusters"
obj_annotation=further_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="0"), "CD8T.Tem.GZMKhi", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="1"), "CD8T.Tem.GZMK_NKG7", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(further_subset_try, idents="2"), "CD8T.Tem.GZMKhi", annotation))
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(5,7), by=0.1, res=0.3, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: cytotoxic, NK.CD56dim.FCER1G
VlnPlot(CD4T_obj_subset_try, features=c("MYOM2","FCER1G","PRF1","SPON2"), pt.size=0) # CD56dim CD57int, # 0 2
# cluster1: less activated but more proliferative, NK.CD56dim.HNRNPH1
VlnPlot(CD4T_obj_subset_try, features=c("C1orf56","HNRNPH1","CDC42SE1"), pt.size=0) # No.
# cluster2: adaptive NK (different from adaptive-like because KLRC2+, which is a marker of human cytomegalovirus/HCMV infection), NK.CD56dim.FCER1Glo_KLRC2pos
VlnPlot(CD4T_obj_subset_try, features=c("KLRC2","CD3E","CCL5","VIM"), pt.size=0) # CD56dim CD57+, # 4
# cluster4: adaptive-like NK (because KLRC2-; this is possibly a subset derived from infections), NK.CD56dim.FCER1Glo_KLRC2neg
VlnPlot(CD4T_obj_subset_try, features=c("CD3E","VIM","IL32","CD52","KLRC2","FCER1G"), pt.size=0) # 1
# cluster3: transcriptionally active (this is a homeostatic activation; apart from below also express "NFKBIA","CD69","CXCR4","ZFP36"), NK.CD56dim.FOS
VlnPlot(CD4T_obj_subset_try, features=c("FOS","FOSB","JUN","JUNB"), pt.size=0) # CD56dim CD57low, No.
# cluster5: IFN-I+ NK (reported to have been found in early severe COVID-19 NK cells, so it's related to infections), NK.CD56dim.ISG
VlnPlot(CD4T_obj_subset_try, features=c("MX1","ISG15","ISG20"), pt.size=0) # No.
# cluster1: transitional NK between CD56dim and CD56hi, NK.CD56hi.FCGR3A
VlnPlot(CD4T_obj_subset_try, features=c("IL7R","SELL","KLRC1","CD44","XCL1","FCGR3A","GZMK","XCL2","CD160"), pt.size=0) # CD56dim CD57-, # 3
# cluster0: NKT
VlnPlot(CD4T_obj_subset_try, features=c("CD3D","CD3G","CD8A","CD8B","NKG7","KLRD1","TYROBP","GNLY","FCER1G"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "NK.CD56dim.FCER1G", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "NK.CD56dim.FCER1Glo_KLRC2neg", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "NK.CD56dim.FCER1G", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "NK.CD56hi.FCGR3A", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "NK.CD56dim.FCER1Glo_KLRC2pos", annotation))
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
CD4T_obj_subset_try=FindClusters2(CD4T_obj_subset, cluster.range=c(6,8), by=0.1, res=0.3, verbose=T)
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()

# cluster0: AXL+ mDC, DC.AXL_mDC
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","FCER1A","CLEC10A"), pt.size=0) # 6
# cluster1: cDC1, DC.cDC1
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","THBD","CLEC9A","CADM1","XCR1","BATF3"), pt.size=0) # 5
# cluster2: AXL+ pDC, DC.AXL_pDC
VlnPlot(CD4T_obj_subset_try, features=c("HLA-DRB1","HLA-DRA","PPP1R14A","AXL","SCT","TPM2","LRRC26"), pt.size=0) # No.
# cDC2_1, DC.cDC2_1
VlnPlot(CD4T_obj_subset_try, features=c("FCER1A", "CD14", "CLEC10A", "CTSS", "ENHO", "CD1C", "MRC1", "FCGR2B", "PID1", "IL13RA1"), pt.size=0) # 2, cDCs
# cDC2_2, DC.cDC2_2
VlnPlot(CD4T_obj_subset_try, features=c("FCER1A", "BASP1", "CD1C", "CD74", "CLEC10A", "HLA-DPA1", "ENHO", "HLA-DPB1", "PLD4", "HLA-DQA1"), pt.size=0) # 1, cDCs
# cluster1: pDC, DC.pDC
VlnPlot(CD4T_obj_subset_try, features=c("ITM2C","PLD4","SERPINF1","LILRA4","IL3RA","TPM2","CLEC4C"), pt.size=0) # 0 3 4, pDCs
# cluster5: cDC2 with cell proliferative property, DC.cDC2_2.MKI67
VlnPlot(CD4T_obj_subset_try, features=c("PCLAF","TYMS","MKI67"), pt.size=0) # No.

Idents(CD4T_obj_subset_try)="seurat_clusters"
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "DC.pDC", Cluster_numbers)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "DC.cDC2_2", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "DC.cDC2_1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="3"), "DC.pDC", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="4"), "DC.pDC", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="5"), "DC.cDC1", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="6"), "DC.AXL_mDC", annotation))
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"), pt.size=0) # 0 1 3 4 5 7, Naive+Naive_IFN
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
VlnPlot(CD4T_obj_subset_try, features=c("IGHM","IGHD","CD27","LINC01857","RALGPS2","TCL1A","IL4R","CD5"), pt.size=0) # 2

Idents(CD4T_obj_subset_try)="seurat_clusters"
DimPlot(CD4T_obj_subset_try, label=T, label.size=5, repel=T, reduction="umap") + NoLegend()
obj_annotation=CD4T_obj_subset_try[[]] %>% tibble::rownames_to_column("cell") %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="0"), "B.naive.kappa", seurat_clusters)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="1"), "B.naive.kappa", annotation)) %>%
  mutate(annotation=ifelse(cell %in% WhichCells(CD4T_obj_subset_try, idents="2"), "B.inter.CD27neg", annotation)) %>%
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
Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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

Seurat_whole_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
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
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")

# add the agecut metadata
df=Seurat_whole_obj[[]] %>%
  mutate(agecut=ifelse(Age %in% c(25:30), "25~30",
                       ifelse(Age %in% c(31:47), "31~47",
                              ifelse(Age %in% c(48:68), "48~68",
                                     ifelse(Age %in% c(69:81), "69~81", NA)))))
table(df$agecut, useNA='ifany')
table(rownames(df)==colnames(Seurat_whole_obj)) # check
Seurat_whole_obj=AddMetaData(Seurat_whole_obj, metadata=df$agecut, col.name="agecut")
# normalize and scale data
Seurat_whole_obj=Seurat::NormalizeData(Seurat_whole_obj) %>% Seurat::FindVariableFeatures(.) %>% Seurat::ScaleData(.)
saveRDS(Seurat_whole_obj, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")

### Save metadata for RawCountObj_OneSamplePerDonor.rds
meta_data=Seurat_whole_obj[[]] %>% tibble::rownames_to_column("cell_id")
data.table::fwrite(meta_data, "~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
unlink("~/Project_PBMCage/Immunity/all_pbmcs/Reannotation.csv.gz")

#####################################



### Plot the ribosomal percentage
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

meta_data1=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
meta_data2=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
meta_data_merged=data.table::rbindlist(list(meta_data1, meta_data2))
meta_data_merged=meta_data_merged %>%
  mutate(Sex=ifelse(Sex=="Male","M","F"))

### Plot the rb
# plot all the celltypes in general
meta_info_all=meta_data_merged %>% group_by(Donor_id, Sex, Age) %>% summarize_at("percent.ribo", mean)
# plot_rb_all=
  ggplot(meta_info_all, aes(x=Age, y=percent.ribo, color=Sex)) +
  geom_point(alpha=0.25) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="Terekhova et al. (2023) dataset") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=Sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=Sex), show.legend=FALSE, size=3.5, method="pearson")

# split Annot.rough
meta_info_subset=meta_data_merged %>% group_by(Tube_id, Sex, Age, Annot.rough) %>% summarize_at("percent.ribo", mean)
plot_rb_sep=
  ggplot(meta_info_subset, aes(x=Age, y=percent.ribo, color=Sex)) +
  geom_point(alpha=0.1, shape=20, size=0.1) +
  facet_wrap(~Annot.rough, ncol=4) +
  theme_classic() +
  labs(x="age", y="Percent of ribosomal genes (%)", title="Terekhova et al. (2023) dataset", subtitle=" ") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=Sex), method="lm", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  # ggpubr::stat_cor(aes(group=Sex), show.legend=FALSE, size=3.5, method="pearson") +
  ggpmisc::stat_correlation(ggpmisc::use_label(c("R","p")), small.p=T, show.legend=FALSE, size=3.5, method="pearson", vstep=0.2)

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_rb_rough.pdf", height=3.5, width=6)
plot(plot_rb_sep)
dev.off()

# split Annot.detailed of CD8T
# ... in the 3_5_datasets_confirm.R

#####################################



### Analyze the celltype freq and their correlation with ages (representative cells)
#####################################
###

### Load the meta data with celltypes
merged_meta=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
meta_info=merged_meta %>% mutate(Sex=ifelse(Sex=="Male","M","F")) %>%
  dplyr::rename(donor_id=Donor_id, sex=Sex, age=Age)
# calculate the freq
meta_info_per.rough=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% dplyr::count() %>% mutate(N_per_rough=n) %>% select(-n)
meta_info_per.inter=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter) %>% dplyr::count() %>% mutate(N_per_inter=n) %>% select(-n)
meta_info_per.detailed=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% dplyr::count() %>% mutate(N_per_detailed=n) %>% select(-n)
meta_info_all=meta_info %>% group_by(donor_id, sex, age) %>% dplyr::count() %>% mutate(N=n) %>% select(-n)
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
write.table(COR_combined, "~/Project_PBMCage/Results/Immunity_results/Cell_frequency_correlation.txt", sep="\t")
write.table(meta_info_combined, "~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")

#####################################



### Plot the celltype freq at the Annot.inter level (use representative data)
#####################################
###

library(ggplot2)

### Plot the Annot.inter freq
meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
meta_info_subset=meta_info_combined %>%
  # subset(grepl("CD4T\\.|CD8T\\.|NK\\.|OtherT\\.", Annot.inter)) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|inter",colnames(.))])) %>%
  subset(!duplicated(.))
# select the populations coexisted in PBMCage
meta_info_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
detailed_celltypes_existed=unique(meta_info_PBMCage$Annot.inter)
meta_info_subset=meta_info_subset %>% subset(Annot.inter %in% detailed_celltypes_existed)

plot_inter=
  ggplot(meta_info_subset, aes(x=age, y=percent.inter, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.inter, ncol=7, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="Terekhova et al. (2023) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10),
        legend.position=c(0.975,0.1),
        legend.text=element_text(size=9),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.25, fill="gray93", alpha=0.25) +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cell_frequency_inter.pdf", width=14, height=5)
plot(plot_inter)
dev.off()

#####################################



### Analyze the celltype freq and their correlation with ages (all cells)
#####################################
###

library(dplyr)

### Load the meta data with celltypes
metadata_1=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
metadata_2=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
merged_meta=data.table::rbindlist(list(metadata_1, metadata_2))
meta_info=merged_meta %>% mutate(Sex=ifelse(Sex=="Male","M","F")) %>%
  dplyr::rename(donor_id=Tube_id, sex=Sex, age=Age)
# calculate the freq
meta_info_per.rough=meta_info %>% group_by(donor_id, sex, age, Annot.rough) %>% dplyr::count() %>% mutate(N_per_rough=n) %>% select(-n)
meta_info_per.inter=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter) %>% dplyr::count() %>% mutate(N_per_inter=n) %>% select(-n)
meta_info_per.detailed=meta_info %>% group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% dplyr::count() %>% mutate(N_per_detailed=n) %>% select(-n)
meta_info_all=meta_info %>% group_by(donor_id, sex, age) %>% dplyr::count() %>% mutate(N=n) %>% select(-n)
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
write.table(COR_combined, "~/Project_PBMCage/Results/Immunity_results/AllCells_Cell_frequency_correlation.txt", sep="\t")
write.table(meta_info_combined, "~/Project_PBMCage/Results/Immunity_results/AllCells_Cell_frequency_values.txt", sep="\t")

#####################################



### Plot the celltype freq at the Annot.detailed level
#####################################
###

library(ggplot2)

### Load
meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/AllCells_Cell_frequency_values.txt", sep="\t")
meta_info_subset=meta_info_combined %>%
  subset(grepl("CD4T\\.|CD8T\\.|NK\\.", Annot.detailed)) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.))
# select the populations coexisted in PBMCage
meta_info_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
detailed_celltypes_existed=unique(meta_info_PBMCage$Annot.detailed)
meta_info_subset=meta_info_subset %>% subset(Annot.detailed %in% detailed_celltypes_existed)

plot_detailed=
  ggplot(meta_info_subset, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, ncol=8, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="Terekhova et al. (2023) dataset") +
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cell_frequency_detailed.pdf", width=22, height=6)
plot(plot_detailed)
dev.off()

### the remaining: B cells, DCs, Monocytes, Other cells, OtherT
# ... by comparing with "~/Project_PBMCage/Plots/PBMCage_Plots/Cell_frequency_detailed_theRest.pdf",
# ... no change can be found/confirmed in these celltypes
meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/AllCells_Cell_frequency_values.txt", sep="\t")
celltype_selected=list(B=meta_info_combined$Annot.detailed[grepl("^B\\.", meta_info_combined$Annot.detailed)],
                       DC=meta_info_combined$Annot.detailed[grepl("^DC\\.", meta_info_combined$Annot.detailed)],
                       Mono=meta_info_combined$Annot.detailed[grepl("^Mono\\.", meta_info_combined$Annot.detailed)],
                       Other_cells=meta_info_combined$Annot.detailed[grepl("^Other\\.", meta_info_combined$Annot.detailed)],
                       OtherT=meta_info_combined$Annot.detailed[grepl("^OtherT\\.", meta_info_combined$Annot.detailed)])
meta_info_combined=meta_info_combined %>% subset(Annot.detailed %in% unlist(celltype_selected)) %>%
  select(any_of(colnames(.)[grepl("^donor_id$|^sex$|^age$|detailed",colnames(.))])) %>%
  subset(!duplicated(.))
# select the populations coexisted in PBMCage
meta_info_PBMCage=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
detailed_celltypes_existed=unique(meta_info_PBMCage$Annot.detailed)
meta_info_subset=meta_info_combined %>% subset(Annot.detailed %in% detailed_celltypes_existed)

plot_detailed=
  ggplot(meta_info_subset, aes(x=age, y=percent.detailed, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
  facet_wrap(~Annot.detailed, ncol=8, scales="free_y") +
  theme_classic() +
  labs(x="age", y="Cell percent (%)", title="Terekhova et al. (2023) dataset") +
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Cell_frequency_detailed_theRest.pdf", width=22, height=6)
plot(plot_detailed)
dev.off()

#####################################



### Causual analysis between freq of CD4T.naive.COTL1pos_SOX4pos and Th2 CD4 cells
#####################################
###

if (!require("CausalImpact")) install.packages("CausalImpact")
library(CausalImpact)

meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
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
table(Reduce(c, lapply(Th2_res, function(x) x$p))>0.05) # turns out to be >99% for p>0.05
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

# plot the results
library(gridExtra)
library(ggplot2)
plot_infer=
  ggplot(merged_df, aes(x=time, y=response)) +
  facet_wrap(~intervention, scales="free") +
  geom_line(inherit.aes=F, aes(x=time, y=mean), linewidth=0.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, xmin=time, xmax=time), fill="gray80", alpha=0.5) +
  theme_classic() +
  labs(y="Cumulative effect (%)", x="donor #", title="Terekhova et al. (2023) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        strip.text=element_text(size=10)
  )

plot_infer_result.table=
  ggplot(data.frame(x=1, y=1), aes(x=x, y=y)) +
  theme_void() +
  theme(plot.margin=margin(0.5,0.1,0.5,0.1, "cm")) +
  annotation_custom(tableGrob(summary_df_merged %>% tibble::column_to_rownames("Intervention"),
                              theme=ttheme_minimal(base_size=8.2)),
                    xmin=-0.05, ymax=0.5)

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_CausualInference_CD4TnaiveSOX4pos_cor.W.Th2cells.pdf", width=6, height=3)
plot_infer
plot_infer_result.table
dev.off()

#####################################



### Make pseudobulk obj (representative tubes only)
#####################################
###

library(dplyr)
library(Seurat)
Seurat_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
sex_df=data.frame(donor_id=Seurat_obj$Donor_id,
                  sex=Seurat_obj$Sex)
sex_df=sex_df[!duplicated(sex_df$donor_id),]

# pseudobulk at the rough level
celltypes=names(table(Seurat_obj$Annot.rough))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.rough==celltypes[i])
  pseudo_rough=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("Donor_id","Age"))
  pseudo_rough$donor_id=gsub("_.*","",colnames(pseudo_rough))
  pseudo_rough$age=as.numeric(gsub(".*_","",colnames(pseudo_rough)))
  pseudo_rough$Annot.rough=celltypes[i]
  PseudoList[[i]]=pseudo_rough
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
# -- add sex metadata
tempt_annot_df=pbmc.merged[[]]
tempt_annot_df$sex=sex_df$sex[match(gsub("-","_",tempt_annot_df$donor_id), as.character(sex_df$donor_id))]
pbmc.merged=AddMetaData(pbmc.merged, metadata=tempt_annot_df$sex, col.name="sex")
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("Male","Female")
saveRDS(pbmc.merged, "~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")

# pseudobulk at the inter level
celltypes=names(table(Seurat_obj$Annot.inter))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.inter==celltypes[i])
  pseudo_inter=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("Donor_id","Age"))
  pseudo_inter$donor_id=gsub("_.*","",colnames(pseudo_inter))
  pseudo_inter$age=as.numeric(gsub(".*_","",colnames(pseudo_inter)))
  pseudo_inter$Annot.inter=celltypes[i]
  PseudoList[[i]]=pseudo_inter
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
# -- add celltypes metadata
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
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("Male","Female")
saveRDS(pbmc.merged, "~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")

# pseudobulk at the detailed level
celltypes=names(table(Seurat_obj$Annot.detailed))
PseudoList=list()
for (i in 1:length(celltypes)) {
  print(paste0("Pseudobulk processing on ",celltypes[i]," starts..."))
  Seurat_obj_subset=subset(Seurat_obj, Annot.detailed==celltypes[i])
  pseudo_detailed=AggregateExpression(Seurat_obj_subset, assays="RNA", return.seurat=T, group.by=c("Donor_id","Age"))
  pseudo_detailed$donor_id=gsub("_.*","",colnames(pseudo_detailed))
  pseudo_detailed$age=as.numeric(gsub(".*_","",colnames(pseudo_detailed)))
  pseudo_detailed$Annot.detailed=celltypes[i]
  PseudoList[[i]]=pseudo_detailed
  print(paste0("--- Pseudobulk processing on ",celltypes[i]," Done! ---"))
}
pbmc.merged=merge(PseudoList[[1]], PseudoList[-1], add.cell.ids=celltypes, merge.data=FALSE, project="PBMCage_pseudobulk")
pbmc.merged=JoinLayers(pbmc.merged)
# -- add celltypes metadata
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
pbmc.merged$sex=as.factor(pbmc.merged$sex); levels(pbmc.merged$sex)=c("Male","Female")
saveRDS(pbmc.merged, "~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")

#####################################



### Make pseudobulk obj on donors
#####################################
###

library(dplyr)
library(Seurat)

### Load the meta data
meta_data1=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
meta_data1=meta_data1 %>% dplyr::select(Donor_id, Tube_id, Sex, Age) %>% subset(!duplicated(.))
meta_data2=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
meta_data2=meta_data2 %>% dplyr::select(Donor_id, Tube_id, Sex, Age) %>% subset(!duplicated(.))
meta_data=rbind(meta_data1, meta_data2)

### Load the count data
# pseudbulk on the OneSamplePerDonor part per donor
pbmc.merged=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
pseudo_donor=Seurat::AggregateExpression(pbmc.merged, assays="RNA", return.seurat=T, group.by=c("donor_id"))
pseudo_donor[[]]=pseudo_donor[[]] %>%
  tibble::rownames_to_column("Donor_id") %>%
  left_join(., meta_data1, by="Donor_id") %>%
  rename(donor_id=Donor_id, tube_id=Tube_id, sex=Sex, age=Age)
saveRDS(pseudo_donor, "~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_onDonor.rds")
# pseudbulk on the rest part per tube
pbmc.merged2=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest.rds")
pseudo_donor2=Seurat::AggregateExpression(pbmc.merged2, assays="RNA", return.seurat=T, group.by=c("Tube_id"))
pseudo_donor2[[]]=pseudo_donor2[[]] %>%
  tibble::rownames_to_column("Tube_id") %>%
  left_join(., meta_data2, by="Tube_id") %>%
  rename(donor_id=Donor_id, tube_id=Tube_id, sex=Sex, age=Age)
saveRDS(pseudo_donor2, "~/Project_PBMCage/Immunity/Immunity_rest_Pseudobulk_onTube.rds")
# pseudbulk on the OneSamplePerDonor+rest parts per tube
pbmc.merged1=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
pseudo_donor1=Seurat::AggregateExpression(pbmc.merged1, assays="RNA", return.seurat=T, group.by=c("Tube_id"))
pseudo_donor1[[]]=pseudo_donor1[[]] %>%
  tibble::rownames_to_column("Tube_id") %>%
  left_join(., meta_data1, by="Tube_id") %>%
  rename(donor_id=Donor_id, tube_id=Tube_id, sex=Sex, age=Age)
pseudo_donors=merge(pseudo_donor1, pseudo_donor2)
pseudo_donors=JoinLayers(pseudo_donors)
saveRDS(pseudo_donors, "~/Project_PBMCage/Immunity/Immunity_all_Pseudobulk_onTube.rds")

#####################################



### Phenotype scoring
#####################################
###

### Loading
library(Seurat)
library(UCell)
library(dplyr)
object=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")

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
score_df=cbind(score_df, CD8T_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD8TandNKT.rds")

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
score_df=cbind(score_df, CD4T_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD4T.rds")

### Score B cells
Bcell_subset=subset(object, Annot.rough=="B cells")
B_cell_scoring=list(naive=c("IGHM","IGHD","CD79A","IL4R","MS4A1","CXCR4","BTG1","TCL1A","CD79B","YBX3"),
                    `isotype switching`=c("AICDA","ATAD5","BATF","CD40LG","ERCC1","EXO1","EXOSC3","EXOSC6","ICOSLG","LIG4","MLH1","MSH2","MSH6","NBN","NFKBIZ","RNF8","RNF168","SANBR","SWAP70","XRCC4","UNG"),
                    activation=c("ABL1","CD19","CD180","GAPT","PLCL2","TLR4","IRF8","PLCG2","SPI1","ADA","CDH17","IL6","IL21","ITFG2","PHF14","POU2AF1","BCL3","CDH17","DLL1","DOCK10","DOCK11","LFNG","MFNG","PTK2B","NOTCH2","TNFAIP3","BCL6","ENPP1","IL2","IL10","ITM2A","LGALS1","LGALS8","NFKBIZ","NKX2-3","XBP1","ST3GAL1"),
                    migration=c("CH25H","CXCL13","HSD3B7","CYP7B1","PTK2B","GAS6","XCL1"),
                    maturation=c("ABL1","ATM","ATP11C","CD24","FNIP1","FOXP1","IGHM","IRF2BP2","KIT","LRRC8A","PRKDC","RAG1","RAG2","SPI1","SPIB","SYVN1","TNFSF13B","TRAF3IP2"),
                    `Ig production`=c("CRLF2","ENPP1","FGL2","GAPT","NFKBIZ","NOD2","POU2F2","TREX1"),
                    inhibitory=c("BCL6","LILRB4","ATM","BTK","BTLA","CASP3","CD24","CD300A","CDKN2A","CTLA4","FCGR2B","IL10","INPP5D","LYN","PAWR","PKN1","PRDM1","TNFRSF13B","PTEN","RC3H1","TNFRSF21","TSC2","TYROBP","NDFIP1","FOXP3","PARP3"),
                    `Ag presentation`=c("FCER2","IGHE"),
                    proliferation=c("BAX","RAG2","ABL1","CD19","CD180","GAPT","PLCL2","TLR4"),
                    transcription=c("BCL2","PRKDC","RAG2","TCF3","TRP53"),
                    costimulation=c("CD320","IL4","TNFSF13B","TNFSF13C"))
Bcell_Ucell=AddModuleScore_UCell(Bcell_subset, features=B_cell_scoring)
# extract the ucell scores
score_df=Bcell_Ucell[[]][, colnames(Bcell_Ucell[[]])[grepl("_UCell$", colnames(Bcell_Ucell[[]]))]]
score_df=cbind(score_df, Bcell_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_Bcells.rds")

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
score_df=cbind(score_df, Mono_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_Monocytes.rds")

### Score DCs
DC_subset=subset(object, Annot.rough=="DCs")
DC_scoring=list(maturation=c("CCL19","CCR7","F2RL1","PRTN3","BATF","BATF2","BATF3","CAMK4","CSF2","DCSTAMP","GIMAP3","GIMAP5","IL4","IRF4","LTBR","NOTCH2","PIRB","PSEN1","RELB","SPI1","TNFSF9","TRAF6","UBD","IRF8","SLAMF9","ZBTB46"),
                activation=c("DOCK2","PYCARD","SLAMF1"),
                migration=c("ALOX5","ANO6","ASB2","CCL19","CCL21","CDC42","DOCK8","EPS8","EXT1","GPR183","NLRP12","TRPM2","TRPM4"),
                `cytokine prod.`=c("CLEC7A","DDX1","DDX21","DHX36","KIT","MAVS","NOD2","PLCG2","RIGI","SCIMP","SLAMF9","TICAM1"),
                `Ag presentation`=c("HLA-DRB1","HLA-DRB3","CLEC4A","HLA-DRA"),
                inhibitory=c("THBS1","CD68","FGL2","JAK3","BST2","HAVCR2","CD37","IL10","TSPAN32","FCGR2B"),
                proliferation=c("TBK1","AZI2","BCL2L1"),
                transcription=c("BATF3","IRF8","ID2","IRF4","NFIL3","SPI1","TCF4"))
DC_Ucell=AddModuleScore_UCell(DC_subset, features=DC_scoring)
# extract the ucell scores
score_df=DC_Ucell[[]][, colnames(DC_Ucell[[]])[grepl("_UCell$", colnames(DC_Ucell[[]]))]]
score_df=cbind(score_df, DC_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_DCs.rds")

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
score_df=cbind(score_df, NK_Ucell[[]][, c("Age","Sex","agecut","Annot.detailed","Annot.inter","Annot.rough")]) %>%
  dplyr::rename(age=Age)
saveRDS(score_df,"~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_NK.rds")

#####################################



### Analyze the B cell scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_Bcells.rds")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
B_final=draw(ht_B, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_Bcells.pdf", height=4, width=4.5)
draw(ht_B)
dev.off()

#####################################



### Analyze the CD4T scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD4T.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD4T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD4T.Treg.cytotoxic","CD4T.Treg.KLRB1_RORC")))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
CD4T_final=draw(ht_CD4, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_CD4T.pdf", height=4, width=5)
draw(ht_CD4)
dev.off()

#####################################



### Analyze the CD8T scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD8TandNKT.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD8T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD8T.Trm")))

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
CD8T_final=draw(ht_CD8, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_CD8T.pdf", height=3.5, width=5)
draw(ht_CD8)
dev.off()

#####################################



### Analyze the NK scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_NK.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="NK cells")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
NK_final=draw(ht_NK, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_NK.pdf", height=2.5, width=4.5)
draw(ht_NK)
dev.off()

#####################################



### Analyze the Monocytes scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_Monocytes.rds")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
Mono_final=draw(ht_Mono, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_Monocytes.pdf", height=2, width=4.5)
draw(ht_Mono)
dev.off()

#####################################



### Analyze the DC scoring results and plot
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_DCs.rds")

### Correlation analysis
df=CD8T.NKT_df %>% group_by(age, agecut, Annot.detailed) %>%
  summarize_at(colnames(.)[!grepl("Sex|age|Annot\\.",colnames(.))], mean)
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_DCs.pdf", height=2.5, width=4)
draw(ht_DC)
dev.off()

#####################################



### Analyze the CD4T scoring results in females vs. males seperately
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD4T.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD4T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD4T.Treg.cytotoxic","CD4T.Treg.KLRB1_RORC")))

### Correlation analysis on females
df=CD8T.NKT_df %>% group_by(Sex, age, agecut, Annot.detailed) %>% dplyr::rename(sex=Sex) %>%
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
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
df=CD8T.NKT_df %>% group_by(Sex, age, agecut, Annot.detailed) %>% dplyr::rename(sex=Sex) %>%
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_CD4T_Fvs.M.pdf", height=4, width=5)
draw(ht_CD4)
dev.off()

#####################################



### Analyze the CD8T scoring results in females vs. males seperately
#####################################
###

library(Seurat)
library(UCell)
library(ComplexHeatmap)

### Loading the score data
CD8T.NKT_df=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/FunctionScore_CD8TandNKT.rds")
CD8T.NKT_df=subset(CD8T.NKT_df, Annot.rough=="CD8T cells")
# remove some celltypes that do not exist in the PBMCage dataset
table(CD8T.NKT_df$Annot.detailed, CD8T.NKT_df$agecut)
CD8T.NKT_df=CD8T.NKT_df %>%
  subset(!(Annot.detailed %in% c("CD8T.Trm")))

### Correlation analysis on females
df=CD8T.NKT_df %>% group_by(Sex, age, agecut, Annot.detailed) %>% dplyr::rename(sex=Sex) %>%
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
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
df=CD8T.NKT_df %>% group_by(Sex, age, agecut, Annot.detailed) %>% dplyr::rename(sex=Sex) %>%
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_Scoring_Phenotye_CD8T_Fvs.M.pdf", height=4, width=5)
draw(ht_CD4)
dev.off()

#####################################



### Analyze the VDJ recombination in T and B cells
#####################################
###

library(dplyr)
library(ggplot2)

### Load the correlation data
donor_cor=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
donor_cor_mixedsex=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes))
donor_cor_female=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes) & analysis=="females")
donor_cor_male=donor_cor %>% subset(grepl("^B\\.|^CD4T\\.|^CD8T\\.", celltypes) & analysis=="males")

### Check the VDJ-related genes
IgRecom_genes=c("ATM","CYREN","DCAF1","DCLRE1C","EZH2","FOXP1","HMGB1","LIG4","NHEJ1","POLB","PRKDC","RAG1","RAG2","TCF3","XRCC4","YY1","XRCC6")
TRecom_genes=c("ATM","BCL11B","DCAF1","DCLRE1C","HMGB1","LEF1","LIG4","PRKDC","RAG1","RAG2","TCF7","XRCC6")
Ig_gene_exist=IgRecom_genes[IgRecom_genes %in% donor_cor_mixedsex$gene]
Tcell_gene_exist=TRecom_genes[TRecom_genes %in% donor_cor_mixedsex$gene]

### Extract the expression data
pseudobulk_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
TB_expr_data=Seurat::FetchData(pseudobulk_obj, vars=c("donor_id","sex","age","Annot.detailed","Annot.inter","Annot.rough",Ig_gene_exist,Tcell_gene_exist),
                               layer="data")

### Plot the VDJ-related gene expr in T naive cells
TB_expr_data_df=TB_expr_data %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot", colnames(.))], names_to="gene", values_to="expr") %>%
  subset(gene %in% Tcell_gene_exist) %>%
  group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% summarize_at("expr", mean) %>%
  mutate(sex=ifelse(sex=="Female","F","M"))
TB_expr_data_df$sex=forcats::fct_relevel(TB_expr_data_df$sex, c("F","M"))

VDJ_T_plot=
  ggplot(subset(TB_expr_data_df, Annot.inter %in% c("CD4T.naive","CD8T.naive")), aes(x=age, y=expr, color=sex)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.05) +
  facet_wrap(~Annot.inter, ncol=1) +
  theme_classic() +
  labs(x="age", y="Expression", title=NULL) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        title=element_text(size=10),
        strip.text=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
  ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=4)

### Plot the VDJ-related gene expr in B.naive and B.inter cells
TB_expr_data_df=TB_expr_data %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot", colnames(.))], names_to="gene", values_to="expr") %>%
  subset(gene %in% Ig_gene_exist) %>%
  group_by(donor_id, sex, age, Annot.rough, Annot.inter, Annot.detailed) %>% summarize_at("expr", mean) %>%
  mutate(sex=ifelse(sex=="Female","F","M"))
TB_expr_data_df$sex=forcats::fct_relevel(TB_expr_data_df$sex, c("F","M"))

VDJ_Bcell_plot=
  ggplot(subset(TB_expr_data_df, Annot.inter %in% c("B.naive","B.inter") & Annot.detailed!="B.inter.lambda.IL4R"), # checked, B.inter.lambda.IL4R has too few observations
         aes(x=age, y=expr, color=Annot.detailed)) +
  ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.05) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")) +
  # facet_wrap(~Annot.detailed, ncol=1) +
  theme_classic() +
  labs(x="age", y="Expression", title=NULL) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=11),
        axis.title.y=element_text(size=11),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  geom_smooth(aes(group=Annot.detailed), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
  ggpubr::stat_cor(aes(group=Annot.detailed), show.legend=FALSE, size=4, label.x.npc=0.35, label.y.npc=1)

### Plot both
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_VDJrecombination_TandB.pdf", height=5, width=7)
cowplot::plot_grid(plotlist=list(VDJ_T_plot, VDJ_Bcell_plot), ncol=2, align="hv", axis="tb", rel_widths=c(1,2))
dev.off()

#####################################



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
Immunity_obj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
Immunity_obj_expr_data=Seurat::FetchData(Immunity_obj, vars=c("donor_id","sex","age","Annot.detailed","Annot.inter","Annot.rough",
                                                              anti_apoptosis_genes, pro_apoptosis_genes),
                               layer="data")

### Clean the data
genes_deselected=names(colSums(Immunity_obj_expr_data[,7:ncol(Immunity_obj_expr_data)])!=0)[colSums(Immunity_obj_expr_data[,7:ncol(Immunity_obj_expr_data)])==0]
ap_expr_data_df=Immunity_obj_expr_data %>% select(-any_of(genes_deselected)) %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("donor_id|sex|age|Annot",colnames(.))],
                      names_to="gene", values_to="expr") %>%
  mutate(pro_or_anti=ifelse(gene %in% intersect(anti_apoptosis_genes, pro_apoptosis_genes), "regulate",
                            ifelse(gene %in% anti_apoptosis_genes, "anti-apoptotic", "pro-apoptotic")))

### Plot
ap_expr_data_df_sum=ap_expr_data_df %>% group_by(across(-c(gene,expr))) %>%
  summarize_at("expr", mean) %>%
  mutate(sex=ifelse(sex=="Female","F","M")) %>%
  subset(pro_or_anti!="regulate") %>%
  subset(Annot.rough %in% c("CD4T cells","CD8T cells")) %>%
  tidyr::pivot_wider(names_from="pro_or_anti", values_from="expr") %>%
  mutate(pro_div_anti=`pro-apoptotic`/`anti-apoptotic`)

# all_plot=
ggplot(ap_expr_data_df_sum, aes(x=age, y=pro_div_anti, color=sex)) +
  facet_wrap(~Annot.detailed, ncol=5) +
  scale_x_continuous(breaks=c(30,60,90), limits=c(15,100)) +
  labs(title=NULL, y="Expression", x="age") +
  theme_minimal() +
  geom_smooth(aes(group=sex), se=F, linewidth=0.5, method="loess") +
  theme(title=element_text(size=11),
        plot.subtitle=element_text(size=11, face="italic"),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=12),
        axis.title.x=element_text(size=12),
        panel.grid=element_line(colour="grey80", size=0.05),
        legend.title=element_text(size=10),
        legend.text=element_text(size=11),
        strip.text=element_text(size=11),
        legend.direction="vertical",
        legend.box="vertical",
        legend.position="top") +
  guides(color=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Apoptosis_inter.pdf", height=10.5, width=3)
plot(all_plot)
dev.off()

