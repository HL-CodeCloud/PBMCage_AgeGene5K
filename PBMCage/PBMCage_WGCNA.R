gc()
library(dplyr)
library(WGCNA)
library(hdWGCNA)
library(Seurat)

# ### Make slim object for all the assigned cells, [did not use this!]
# #####################################
# ###
# ### * Too large dataset for WGCNA analysis, so to slim it down
#
# library(dplyr)
# library(Seurat)
# THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
# print("Obj slim starts...")
# subset_=THEObj %>%
#   NormalizeData() %>%
#   FindVariableFeatures()
# subset_=SketchData(
#   subset_,
#   ncells=100000L,
#   sketched.assay="sketch",
#   method="LeverageScore",
#   var.name="leverage.score",
#   seed=2024L,
#   cast="dgCMatrix",
#   verbose=TRUE)
# assay.v5=GetAssay(subset_, assay="sketch")
# subset_=CreateSeuratObject(assay.v5)
# subset_[[]]=THEObj[[]][match(colnames(subset_), colnames(THEObj)),]
# saveRDS(subset_, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_5wCellsForWGCNA.rds")
# print("Obj slim done!")
#
# #####################################



# ### Take representative donor_ids to slim the obj
# #####################################
# ###
# ### * Too large dataset for WGCNA analysis, so to slim it down
#
# library(dplyr)
# library(Seurat)
# THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
#
# # select the donor_ids
# sex_age_df=data.frame(sex=as.character(THEObj$sex), age=as.numeric(as.character(THEObj$age)), donor_id=as.character(THEObj$donor_id)) %>%
#   mutate(sex_w_age=paste0(sex, ", ", age)) %>%
#   subset(!duplicated(.))
# donorid_per_sex.age=split(sex_age_df$donor_id, sex_age_df$sex_w_age)
# set.seed(202406)
# donorid_selected=lapply(donorid_per_sex.age, function(x) if (length(x)>4) {sample(x, size=4)} else {x})
#
# # slim the obj
# slim_obj=subset(THEObj, donor_id %in% unlist(donorid_selected))
# assay.v5=LayerData(slim_obj, assay="RNA", layer="counts")
# cor_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
# cor_data_subset=cor_data %>%
#   subset(analysis=="all" & rho_pval<0.05 & abs(rho)>0.1) %>%
#   split(.$celltypes) %>%
#   lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe()) %>%
#   Reduce(c, .) %>% .[!duplicated(.)]
# assay.v5=assay.v5[cor_data_subset, ]
# subset_=CreateSeuratObject(assay.v5)
# subset_[[]]=THEObj[[]][match(colnames(subset_), colnames(THEObj)),]
#
# # modify metadata
# tempt_df=subset_[[]] %>%
#   mutate(sex=ifelse(sex=="female","F","M")) %>%
#   mutate(agecut.half=ifelse(age %in% c(19:25), "19~30.1:19~25", age)) %>%
#   mutate(agecut.half=ifelse(age %in% c(26:30), "19~30.2:26~30", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(31:35), "31~47.1:31~35", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(36:40), "31~47.2:36~40", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(41:44), "31~47.3:41~44", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(45:47), "31~47.4:45~47", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(48:50), "48~68.1:48~50", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(51:55), "48~68.2:51~55", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(56:60), "48~68.3:56~60", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(61:65), "48~68.4:61~65", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(66:68), "48~68.5:66~68", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(69:70), "69~99.1:69~70", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(71:75), "69~99.2:71~75", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(76:80), "69~99.3:76~80", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(81:85), "69~99.4:81~85", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(86:90), "69~99.5:86~90", agecut.half)) %>%
#   mutate(agecut.half=ifelse(age %in% c(91:99), "69~99.6:91~99", agecut.half)) %>%
#   mutate(sex_w_agecut.half=paste0(sex, ", ", agecut.half))
#
# table(tempt_df$Annot.inter, tempt_df$agecut.half) # check
#
# subset_[[]]=tempt_df
# subset_$age=as.numeric(as.character(subset_$age))
# subset_$sex=as.character(subset_$sex)
# subset_$agecut.half=as.character(subset_$agecut.half)
# subset_$sex_w_agecut.half=as.character(subset_$sex_w_agecut.half)
#
# # determine the threshold and k value for metacell creation
# rowSums(table(subset_$Annot.inter, subset_$sex_w_agecut.half))
# rowMeans(table(subset_$Annot.inter, subset_$sex_w_agecut.half))
# table(subset_$Annot.inter, subset_$sex_w_agecut.half)
# # threshold: sum>250
# # k=20, min_cells=50, max_shared=5, target_metacells=100
#
# # scale the expression
# subset_=subset_ %>%
#   Seurat::NormalizeData() %>%
#   Seurat::FindVariableFeatures() %>%
#   Seurat::ScaleData()
#
# saveRDS(subset_, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_RepresentativeDonors_ForWGCNA.rds")
# print("Obj slim done!")
#
# #####################################



### Prepare metacells for WGCNA analysis (* generate at the inter level so as to analyze both rough and inter)
#####################################
###

### Load
slim_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# arrange the obj
inter_celltypes=names(table(slim_obj$Annot.inter))
slim_obj[["RNA"]]=as(slim_obj[["RNA"]], Class="Assay")
celltypes_detailed=unique(slim_obj$Annot.detailed)
slim_obj[[]]=slim_obj[[]] %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.inter",celltypes_detailed)], "B.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.naive",celltypes_detailed)], "B.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.mem",celltypes_detailed)], "B.mem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("B\\.pc|B\\.pb",celltypes_detailed)], "B.pbpc", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tcm",celltypes_detailed)], "CD4T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Tem",celltypes_detailed)], "CD4T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.prolif",celltypes_detailed)], "CD4T.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.Treg",celltypes_detailed)], "CD4T.Treg", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD4T\\.naive",celltypes_detailed)], "CD4T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tcm",celltypes_detailed)], "CD8T.Tcm", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.Tem",celltypes_detailed)], "CD8T.Tem", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("CD8T\\.naive",celltypes_detailed)], "CD8T.naive", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.AXL",celltypes_detailed)], "DC.AXL", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("DC\\.cDC2",celltypes_detailed)], "DC.cDC2", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.classical",celltypes_detailed)], "Mono.classical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.inter",celltypes_detailed)], "Mono.inter", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Mono\\.nonclassical",celltypes_detailed)], "Mono.nonclassical", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56dim",celltypes_detailed)], "NK.CD56dim", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.CD56hi",celltypes_detailed)], "NK.CD56hi", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("NK\\.prolif",celltypes_detailed)], "NK.prolif", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.progenitor",celltypes_detailed)], "Other.progenitor", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.gdT",celltypes_detailed)], "OtherT.gdT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("OtherT\\.NKT",celltypes_detailed)], "OtherT.NKT", Annot.inter)) %>%
  mutate(Annot.inter=ifelse(Annot.detailed %in% celltypes_detailed[grepl("Other\\.ILC",celltypes_detailed)], "Other.ILC", Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.detailed %in% celltypes_detailed[grepl("pltcontamin|lin_infidel",celltypes_detailed)], "Unassigned", Annot.detailed)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^B\\.",celltypes_detailed)], "B cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD4T\\.",celltypes_detailed)], "CD4T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^CD8T\\.",celltypes_detailed)], "CD8T cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^DC\\.",celltypes_detailed)], "DCs", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Mono\\.",celltypes_detailed)], "Monocytes", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^NK\\.",celltypes_detailed)], "NK cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^Other\\.",celltypes_detailed)], "Other cells", Annot.rough)) %>%
  mutate(Annot.rough=ifelse(Annot.rough %in% celltypes_detailed[grepl("^OtherT\\.",celltypes_detailed)], "OtherT", Annot.rough))
# remove the celltype with only 0 or 1 cell
More_than_oneCell_Celltype=names(table((slim_obj[[]]$Annot.inter)))[!table(slim_obj[[]]$Annot.inter) %in% c(0,1)]
slim_obj=subset(slim_obj, Annot.inter %in% More_than_oneCell_Celltype)

### Take the features
WGCNAobj=SetupForWGCNA(
  slim_obj,
  group.by="Annot.inter",
  gene_select="fraction",
  fraction=0.2,
  wgcna_name="WGCNA")
# save the setup WGCNAobj for later WGCNA analysis on those not using metacells
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")

### Metacell creation
WGCNAobj=WGCNAobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:18) %>%
  RunUMAP(dims=1:18)

### Add metadata
WGCNAobj$age=as.numeric(as.character(WGCNAobj$age))
WGCNAobj[[]]=WGCNAobj[[]] %>%
  mutate(sex=ifelse(sex=="female","F","M")) %>%
  mutate(agecut.half=ifelse(age %in% c(19:25), "19~30.1:19~25", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(26:30), "19~30.2:26~30", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:35), "31~47.1:31~35", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(36:40), "31~47.2:36~40", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(41:44), "31~47.3:41~44", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(45:47), "31~47.4:45~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:50), "48~68.1:48~50", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(51:55), "48~68.2:51~55", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(56:60), "48~68.3:56~60", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(61:65), "48~68.4:61~65", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(66:68), "48~68.5:66~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:70), "69~99.1:69~70", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(71:75), "69~99.2:71~75", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(76:80), "69~99.3:76~80", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(81:85), "69~99.4:81~85", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(86:90), "69~99.5:86~90", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(91:99), "69~99.6:91~99", agecut.half)) %>%
  mutate(sex_w_agecut.half=paste0(sex, ", ", agecut.half))
WGCNAobj$agecut.half=as.character(WGCNAobj$agecut.half)
WGCNAobj$sex_w_agecut.half=as.character(WGCNAobj$sex_w_agecut.half)

# determine the threshold and k value for metacell creation
rowSums(table(WGCNAobj$Annot.inter, WGCNAobj$sex_w_agecut.half))
rowMeans(table(WGCNAobj$Annot.inter, WGCNAobj$sex_w_agecut.half))
table(WGCNAobj$Annot.inter, WGCNAobj$sex_w_agecut.half)

# use determined appropriate k value and min_cells
WGCNAobj=MetacellsByGroups(
  WGCNAobj,
  group.by=c("Annot.inter","sex_w_agecut.half"),
  reduction="umap",
  k=5,
  min_cells=10,
  max_shared=2,
  target_metacells=20,
  mode="sum",
  verbose=TRUE,
  ident.group="Annot.inter",
  wgcna_name="WGCNA")

### Get the metacell obj, as the SetDatExpr step later will face errors if using the seuratobj
tempt_metacellobj=GetMetacellObject(WGCNAobj, wgcna_name="WGCNA")
tempt_metacellobj=JoinLayers(tempt_metacellobj)
tempt_metacellobj=tempt_metacellobj %>%
  NormalizeData(tempt_metacellobj) %>%
  FindVariableFeatures() %>%
  ScaleData()
table(tempt_metacellobj$Annot.inter, tempt_metacellobj$sex_w_agecut.half)
# add metadata
tempt_metacellobj$sex=gsub(", .*","",tempt_metacellobj$sex_w_agecut.half)
tempt_metacellobj$sex=forcats::fct_relevel(tempt_metacellobj$sex, c("M","F"))
tempt_metacellobj$agecut.updated=gsub("\\.[0-9]:.*","",gsub(".*, ","",tempt_metacellobj$sex_w_agecut.half))
tempt_metacellobj$agecut.updated=forcats::fct_relevel(tempt_metacellobj$agecut.updated, c("19~30","31~47","48~68","69~99"))
tempt_metacellobj[[]]=tempt_metacellobj[[]] %>%
  mutate(sex_w_agecut.updated=paste0(sex,", ",agecut.updated))
tempt_metacellobj$sex_w_agecut.updated=forcats::fct_relevel(tempt_metacellobj$sex_w_agecut.updated, c("M, 19~30","F, 19~30","M, 31~47","F, 31~47","M, 48~68","F, 48~68","M, 69~99","F, 69~99"))
# save the WGCNA obj
saveRDS(tempt_metacellobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")

### For the celltypes with unenough cell numbers, WGCNA analysis is not conducted on MetaCells
slim_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_RepresentativeDonors_ForWGCNA.rds")
WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
# determine those to use metacells and those not to
enough_cell_types=table(WGCNAobj$Annot.inter, WGCNAobj$sex_w_agecut.updated) %>% apply(., 1, function(x) min(x)) %>% .[.!=0] %>% names(.)
not_enough_cell_types=names(table(slim_obj$Annot.inter))[!(names(table(slim_obj$Annot.inter)) %in% enough_cell_types)]

### WGCNA on MetaCells
# Resetup on the features as the WGCNAobj is renewed
WGCNAobj=SetupForWGCNA(
  WGCNAobj,
  group.by="Annot.inter",
  gene_select="fraction",
  fraction=0.2,
  wgcna_name=enough_cell_types[1])
print("Resetup for MetaCell analysis is done!")

### Duplicate wgcna slots for all enough_cell_types
full_list=list()
full_list[[1]]=WGCNAobj@misc[[1]]
for (i in 2:(length(enough_cell_types)+1)) {
  full_list[[i]]=WGCNAobj@misc[[2]]
}
names(full_list)=c("active_wgcna", enough_cell_types)
WGCNAobj@misc=full_list
names(WGCNAobj@misc)
length(WGCNAobj@misc[[2]]$wgcna_genes)

### Set the expression data for each wgcna slots
for (i in 1:length(enough_cell_types)) {
  WGCNAobj@misc$active_wgcna=enough_cell_types[i]
  WGCNAobj=SetDatExpr(
    WGCNAobj,
    group_name=WGCNAobj@misc$active_wgcna,
    group.by="Annot.inter",
    assay="RNA",
    use_metacells=TRUE,
    slot="data")
}
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
print("Set expression for MetaCell analysis is done!")

### WGCNA analysis pipeline for those using metacells
for (i in 10:length(enough_cell_types)) {
  WGCNAobj@misc$active_wgcna=enough_cell_types[i]

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

  ### Construct the network
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
    modules=levels(WGCNAobj@misc[[enough_cell_types[i]]][["wgcna_modules"]][["module"]])
    saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
  } else {print("No selected_power found.")}
}

### WGCNA on Seuratobj
WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")
not_enough_cell_types=unique(WGCNAobj$Annot.inter)

### Duplicate wgcna slots for all not_enough_cell_types
full_list=list()
full_list[[1]]=WGCNAobj@misc[[1]]
for (i in 2:(length(not_enough_cell_types)+1)) {
  full_list[[i]]=WGCNAobj@misc[[2]]
}
names(full_list)=c("active_wgcna", not_enough_cell_types)
WGCNAobj@misc=full_list
names(WGCNAobj@misc)
length(WGCNAobj@misc[[2]]$wgcna_genes)

### Set the expression data for each wgcna slots
for (i in 1:length(not_enough_cell_types)) {
  WGCNAobj@misc$active_wgcna=not_enough_cell_types[i]
  WGCNAobj=SetDatExpr(
    WGCNAobj,
    group_name=WGCNAobj@misc$active_wgcna,
    group.by="Annot.inter",
    assay="RNA",
    use_metacells=FALSE,
    slot="data")
}
saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")
print("Set expression for MetaCell analysis is done!")

### WGCNA analysis pipeline for those not using metacells
for (i in 1:length(not_enough_cell_types)) {
  WGCNAobj@misc$active_wgcna=not_enough_cell_types[i]

  ### Choose the power
  selected_power=0
  WGCNAobj=TestSoftPowers(
    WGCNAobj,
    powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
  plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
  possible_power=c()
  for (i in 1:length(plot_list)) {
    tryCatch({
      possible_power[i]=(plot_list[[i]][[1]] %>% subset(text_color=="white"))$Power
    }, error=function(msg) {print("No Power found.")})
  }
  if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
  print(selected_power)

  ### Construct the network
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
    modules=levels(WGCNAobj@misc[[not_enough_cell_types[i]]][["wgcna_modules"]][["module"]])
    saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")
  } else {print("No selected_power found.")}
}

#####################################



# ### Extract the results from the WGCNA analysis at the inter level
# #####################################
# ###
# 
# ### Save genes in each module in the good WGCNA slots
# good_wgcna_slots_names=c()
# good_wgcna_slots=list()
# WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
# # take only good wgcna slots
# names(WGCNAobj@misc)
# for (i in 2:length(WGCNAobj@misc)) {
#   slot_=WGCNAobj@misc[[i]]
#   if (length(slot_)==2) {
#     print(paste0(names(WGCNAobj@misc)[i], ": No good power."))
#   } else if (length(slot_)==3) {
#     print(paste0(names(WGCNAobj@misc)[i], ": No module found."))
#   } else if (length(levels(WGCNAobj@misc[[i]][["wgcna_modules"]]$module))==1) {
#     print(paste0(names(WGCNAobj@misc)[i], ": Only grey module found."))
#   } else {
#     good_wgcna_slots_names=c(good_wgcna_slots_names, names(WGCNAobj@misc)[i])
#     good_wgcna_slots=c(good_wgcna_slots, list(WGCNAobj@misc[[i]][["wgcna_modules"]] %>% mutate(wgcna_slot=names(WGCNAobj@misc)[i])))
#   }
# }
# WGCNAobj_notenoughobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")
# names(WGCNAobj_notenoughobj@misc)
# for (i in 2:length(WGCNAobj_notenoughobj@misc)) {
#   slot_=WGCNAobj_notenoughobj@misc[[i]]
#   if (length(slot_)==2) {
#     print(paste0(names(WGCNAobj_notenoughobj@misc)[i], ": No good power."))
#   } else if (length(slot_)==3) {
#     print(paste0(names(WGCNAobj_notenoughobj@misc)[i], ": No module found."))
#   } else if (length(levels(WGCNAobj_notenoughobj@misc[[i]][["wgcna_modules"]]$module))==1) {
#     print(paste0(names(WGCNAobj_notenoughobj@misc)[i], ": Only grey module found."))
#   } else {
#     good_wgcna_slots_names=c(good_wgcna_slots_names, names(WGCNAobj_notenoughobj@misc)[i])
#     good_wgcna_slots=c(good_wgcna_slots, list(WGCNAobj_notenoughobj@misc[[i]][["wgcna_modules"]] %>% mutate(wgcna_slot=names(WGCNAobj_notenoughobj@misc)[i])))
#   }
# }
# GOOD_WGCNA.modules=data.table::rbindlist(good_wgcna_slots)
# data.table::fwrite(GOOD_WGCNA.modules, "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_modules.csv.gz", sep="\t")
# 
# ### Get the modules, trait-correlation, and kMEs from the good wgcna slots
# # for metacell obj
# WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
# GOOD_WGCNA.modules=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_modules.csv.gz", sep="\t")
# good_wgcna_slots_names=unique(GOOD_WGCNA.modules$wgcna_slot)
# WGCNAobj[["RNA"]]=as(object=WGCNAobj[["RNA"]], Class="Assay")
# WGCNAobj=Seurat::ScaleData(WGCNAobj)
# 
# ME_list=COR_PVAL_FDR_list=Module_kME_list=list()
# for (i in 1:length(good_wgcna_slots_names)) {
#   if (good_wgcna_slots_names[i] %in% names(WGCNAobj@misc)) {
#     WGCNAobj@misc$active_wgcna=good_wgcna_slots_names[i]
#     length(unique(WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$wgcna_modules$color))
# 
#     # MEs
#     module_table_minus=hdWGCNA::GetModules(WGCNAobj) %>% subset(color!="grey")
#     module_table_grey=hdWGCNA::GetModules(WGCNAobj) %>% subset(color=="grey") %>% .[1:100,]
#     module_table_back=rbind(module_table_minus, module_table_grey)
#     WGCNAobj=xpectr::suppress_mw(hdWGCNA::ModuleEigengenes(WGCNAobj, modules=module_table_back))
#     print(paste0(">>>>>> ModuleEigengenes calculation of No.",i,": ",WGCNAobj@misc$active_wgcna, " is done. <<<<<<"))
# 
#     MEs=hdWGCNA::GetMEs(WGCNAobj, wgcna_name=WGCNAobj@misc$active_wgcna) %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column("cell_id") %>%
#       tidyr::pivot_longer(cols=colnames(.)[!grepl("cell_id",colnames(.))], names_to="color", values_to="MEs") %>%
#       mutate(wgcna_slot=WGCNAobj@misc$active_wgcna) %>%
#       subset(color!="grey")
#     ME_list=c(ME_list, list(MEs))
# 
#     print(paste0(">>>>>> MEs of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
# 
#     # correlation
#     if (length(levels(WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]][["wgcna_modules"]]$module))>2) {
#       cur_traits=c("agecut.updated","sex","sex_w_agecut.updated")
#       WGCNAobj=hdWGCNA::ModuleTraitCorrelation(
#         WGCNAobj,
#         traits=cur_traits,
#         features="MEs",
#         cor_method="spearman",
#         group.by='sex'
#       )
#       cor_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["cor"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="cor")
#       pval_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["pval"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="pval")
#       fdr_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["fdr"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="fdr")
#       cor_pval_fdr=full_join(full_join(cor_values, pval_values, by=c("traits","color")), fdr_values, by=c("traits","color"))
#       cor_pval_fdr$wgcna_slot=WGCNAobj@misc$active_wgcna
#       COR_PVAL_FDR_list=c(COR_PVAL_FDR_list, list(cor_pval_fdr))
#     }
#     print(paste0(">>>>>> Correlation of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
# 
#     # kMEs
#     WGCNAobj=hdWGCNA::ModuleConnectivity(
#       WGCNAobj,
#       group.by='Annot.inter',
#       corFnc="spearman",
#       corOptions="use='p'",
#       harmonized=FALSE,
#       assay=NULL,
#       slot="data",
#       group_name=names(table(WGCNAobj$Annot.inter)))
#     Module_genes=hdWGCNA::GetModules(WGCNAobj) %>%
#       as.data.frame() %>%
#       tidyr::pivot_longer(cols=colnames(.)[grepl("kME_",colnames(.))], names_to="kME_group", values_to="kMEs") %>%
#       mutate(wgcna_slot=WGCNAobj@misc$active_wgcna) %>%
#       subset(color!="grey" & kME_group!="kME_grey")
#     Module_kME_list=c(Module_kME_list, list(Module_genes))
# 
#     print(paste0(">>>>>> kMEs of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
# 
#     gc()
#   }
# }
# ME_list_total=data.table::rbindlist(ME_list)
# COR_PVAL_FDR_list_total=data.table::rbindlist(COR_PVAL_FDR_list)
# Module_kME_list_total=data.table::rbindlist(Module_kME_list)
# 
# saveRDS(list(ME_list_total, COR_PVAL_FDR_list_total, Module_kME_list_total), "~/Project_PBMCage/tempt.rds")
# rm(list=ls());gc()
# 
# # for notenough obj
# WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter_NotEnoughCells.rds")
# GOOD_WGCNA.modules=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_modules.csv.gz", sep="\t")
# good_wgcna_slots_names=unique(GOOD_WGCNA.modules$wgcna_slot)
# WGCNAobj[["RNA"]]=as(object=WGCNAobj[["RNA"]], Class="Assay")
# # add metadata
# WGCNAobj$sex=forcats::fct_relevel(WGCNAobj$sex, c("M","F"))
# WGCNAobj$agecut.updated=gsub("-[0-9]","",gsub(".*, ","",WGCNAobj$sex_w_agecut.half))
# WGCNAobj$agecut.updated=forcats::fct_relevel(WGCNAobj$agecut.updated, c("<30","31~47","48~68",">69"))
# WGCNAobj[[]]=WGCNAobj[[]] %>%
#   mutate(sex_w_agecut.updated=paste0(sex,", ",agecut.updated))
# WGCNAobj$sex_w_agecut.updated=forcats::fct_relevel(WGCNAobj$sex_w_agecut.updated, c("M, <30","F, <30","M, 31~47","F, 31~47","M, 48~68","F, 48~68","M, >69","F, >69"))
# 
# for (i in 1:length(good_wgcna_slots_names)) {
#   if (good_wgcna_slots_names[i] %in% names(WGCNAobj@misc)) {
#     WGCNAobj@misc$active_wgcna=good_wgcna_slots_names[i]
#     length(unique(WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$wgcna_modules$color))
#     
#     # MEs
#     module_table_minus=hdWGCNA::GetModules(WGCNAobj) %>% subset(color!="grey")
#     module_table_grey=hdWGCNA::GetModules(WGCNAobj) %>% subset(color=="grey") %>% .[1:100,]
#     module_table_back=rbind(module_table_minus, module_table_grey)
#     WGCNAobj=xpectr::suppress_mw(hdWGCNA::ModuleEigengenes(WGCNAobj, modules=module_table_back))
#     print(paste0(">>>>>> ModuleEigengenes calculation of No.",i,": ",WGCNAobj@misc$active_wgcna, " is done. <<<<<<"))
#     
#     MEs=hdWGCNA::GetMEs(WGCNAobj, wgcna_name=WGCNAobj@misc$active_wgcna) %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column("cell_id") %>%
#       tidyr::pivot_longer(cols=colnames(.)[!grepl("cell_id",colnames(.))], names_to="color", values_to="MEs") %>%
#       mutate(wgcna_slot=WGCNAobj@misc$active_wgcna) %>%
#       subset(color!="grey")
# 
#     print(paste0(">>>>>> MEs of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
#     
#     # correlation
#     if (length(levels(WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]][["wgcna_modules"]]$module))>2) {
#       cur_traits=c("agecut.updated","sex","sex_w_agecut.updated")
#       WGCNAobj=hdWGCNA::ModuleTraitCorrelation(
#         WGCNAobj,
#         traits=cur_traits,
#         features="MEs",
#         cor_method="spearman",
#         group.by='sex'
#       )
#       cor_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["cor"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="cor")
#       pval_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["pval"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="pval")
#       fdr_values=
#         WGCNAobj@misc[[WGCNAobj@misc$active_wgcna]]$mt_cor[["fdr"]][["all_cells"]] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column("traits") %>%
#         tidyr::pivot_longer(cols=colnames(.)[!grepl("traits",colnames(.))], names_to="color", values_to="fdr")
#       cor_pval_fdr=full_join(full_join(cor_values, pval_values, by=c("traits","color")), fdr_values, by=c("traits","color"))
#       cor_pval_fdr$wgcna_slot=WGCNAobj@misc$active_wgcna
#     } else {
#       cor_pval_fdr=NULL
#     }
#     print(paste0(">>>>>> Correlation of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
#     
#     # kMEs
#     WGCNAobj=hdWGCNA::ModuleConnectivity(
#       WGCNAobj,
#       group.by='Annot.inter',
#       corFnc="spearman",
#       corOptions="use='p'",
#       harmonized=FALSE,
#       assay=NULL,
#       slot="data",
#       group_name=names(table(WGCNAobj$Annot.inter)))
#     Module_genes=hdWGCNA::GetModules(WGCNAobj) %>%
#       as.data.frame() %>%
#       tidyr::pivot_longer(cols=colnames(.)[grepl("kME_",colnames(.))], names_to="kME_group", values_to="kMEs") %>%
#       mutate(wgcna_slot=WGCNAobj@misc$active_wgcna) %>%
#       subset(color!="grey" & kME_group!="kME_grey")
# 
#     print(paste0(">>>>>> kMEs of No.",i,": ",WGCNAobj@misc$active_wgcna, " are taken. <<<<<<"))
#     
#     saveRDS(list(MEs, cor_pval_fdr, Module_genes), paste0("~/Project_PBMCage/tempt__o",i,".rds"))
#     rm(list=ls()[!grepl("^WGCNAobj$|^good_wgcna_slots_names$",ls())])
#     gc()
#   }
# }
# 
# # sort all the results
# ME_list=COR_PVAL_FDR_list=Module_kME_list=list()
# for (i in 1:length(c(15,16,17,18,19,21,22))) {
#   result_=readRDS(paste0("~/Project_PBMCage/tempt__o",c(15,16,17,18,19,21,22)[i],".rds"))
#   ME_list[[i]]=result_[[1]]
#   COR_PVAL_FDR_list[[i]]=result_[[2]]
#   Module_kME_list[[i]]=result_[[3]]
# }
# ME_list_total=data.table::rbindlist(ME_list)
# COR_PVAL_FDR_list_total=data.table::rbindlist(COR_PVAL_FDR_list)
# Module_kME_list_total=data.table::rbindlist(Module_kME_list)
# data.table::fwrite(ME_list_total, "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_notenough.wgcna.slots_MEs.csv.gz")
# data.table::fwrite(COR_PVAL_FDR_list_total, "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_notenough.wgcna.slots_Traitcor.csv.gz")
# data.table::fwrite(Module_kME_list_total, "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_notenough.wgcna.slots_kMEs.csv.gz")
# 
# all_tempt=readRDS("~/Project_PBMCage/tempt.rds")
# data.table::fwrite(all_tempt[[1]], "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_MEs.csv.gz")
# data.table::fwrite(all_tempt[[2]], "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_Traitcor.csv.gz")
# data.table::fwrite(all_tempt[[3]], "~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_kMEs.csv.gz")
# 
# #####################################
# 
# 
# 
# ### Enrich the WGCNA gene modules in each celltype
# #####################################
# ###
# 
# library(dplyr)
# library(ggplot2)
# 
# ### Load the modules
# module_gene_list1=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_kMEs.csv.gz", sep=",")
# module_gene_list2=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_notenough.wgcna.slots_kMEs.csv.gz", sep=",")
# module_gene_list=data.table::rbindlist(list(module_gene_list1, module_gene_list2))
# wgcna_slots_list=unique(module_gene_list$wgcna_slot)
# 
# ### Enrich
# compare_results=list()
# for (i in 1:length(wgcna_slots_list)) {
#   module_gene_list_subset=module_gene_list %>% subset(wgcna_slot==wgcna_slots_list[i])
#   module_gene_list_split=lapply(split(module_gene_list_subset$gene_name, module_gene_list_subset$module), function(x) x[!duplicated(x)])
#   compare_results[[i]]=
#     clusterProfiler::compareCluster(module_gene_list_split, 
#                                     fun=clusterProfiler::enrichGO, 
#                                     OrgDb="org.Hs.eg.db", 
#                                     keyType="SYMBOL", 
#                                     pvalueCutoff=0.05,
#                                     minGSSize=1,
#                                     maxGSSize=500,
#                                     ont="BP")
# }
# names(compare_results)=wgcna_slots_list
# saveRDS(compare_results, "~/Project_PBMCage/Tempt_RDS/wgcna_results/WGCNA_enrichGOresults_inter.rds")
# 
# #####################################
# 
# 
# 
# ### Plot the WGCNA gene module-trait correlations
# #####################################
# ###
# 
# library(dplyr)
# library(ggplot2)
# library(ComplexHeatmap)
# 
# ### Check the CD4T and CD8T cell-related slots
# GOOD_WGCNA.modules=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_modules.csv.gz", sep="\t")
# good_wgcna_slots_names=unique(GOOD_WGCNA.modules$wgcna_slot)
# T_related=good_wgcna_slots_names[grepl("CD4T|CD8T",good_wgcna_slots_names)]
# 
# ### Load the trait results
# trait_cor=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_good.wgcna.slots_Traitcor.csv.gz", sep=",")
# trait_cor_notenough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WGCNA_inter_notenough.wgcna.slots_Traitcor.csv.gz", sep=",")
# trait_cor_subset=rbind(trait_cor, trait_cor_notenough) %>%
#   subset(traits=="agecut.updated" & !grepl("Other\\.|CD4T\\.prolif",wgcna_slot))
# 
# ### Plot the heatmap list
# range(trait_cor_subset$cor)
# col_fun_up=circlize::colorRamp2(c(-0.2, 0, 0.2), c("dodgerblue3", "white", "brown3"))
# # arrange the data
# max.length=trait_cor_subset %>% split(.$wgcna_slot) %>% lapply(., nrow) %>% unlist() %>% max(.) +1
# order_by_length=sort(trait_cor_subset %>% split(.$wgcna_slot) %>% lapply(., nrow) %>% unlist()) %>% names()
# 
# ht_list=list()
# for (i in 1:length(order_by_length)) {
#   trait_cor_slot=subset(trait_cor_subset, wgcna_slot==order_by_length[i])
#   cor_df=trait_cor_slot %>% dplyr::select(color, cor, pval) %>% tibble::column_to_rownames("color") %>%
#     mutate(pval=ifelse(pval>0.05, "/", ""))
#   ht_empty=
#     Heatmap(rep(1, max.length-nrow(cor_df)) %>% t(.) %>% as.matrix(), name="empty", 
#             rect_gp=gpar(type="none"), 
#             cluster_columns=F, cluster_rows=F, show_heatmap_legend=F)
#   ht_list[[i]]=
#     Heatmap(cor_df %>% dplyr::select(-pval) %>% t(.),
#             name=" ",
#             show_heatmap_legend=F, 
#             # border_gp=gpar(col="black", lty=2, lwd=0.5),
#             col=col_fun_up,
#             rect_gp=gpar(col="black", lwd=0.5),
#             cluster_columns=T, cluster_rows=T, 
#             show_column_dend=F, show_row_dend=F, column_dend_side="top", row_dend_side="right",
#             row_names_gp=gpar(fontsize=10), show_row_names=F, 
#             row_title=order_by_length[i], row_title_gp=gpar(fontsize=11), row_title_rot=0, row_title_side="right",
#             show_column_names=F,
#             layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
#               v=pindex(cor_df %>% dplyr::select(-cor) %>% t(), i, j)
#               grid.text(v, x, y, rot=0, just="center",
#                         gp=gpar(fontsize=10, fontface="bold"))
#             },
#             height=unit(5, "mm")*(ncol(cor_df)-1), 
#             width=unit(3, "mm")*nrow(cor_df)) + ht_empty
#   ht_list[[i]]=draw(ht_list[[i]], ht_gap=unit(0, "cm"))
# }
# 
# gb_list=lapply(ht_list, function(x) grid.grabExpr(draw(x)))
# ht_all=cowplot::plot_grid(plotlist=gb_list, align="hv", axis="r", ncol=1)
#   
#   
#   ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
# ht_GOterm_sim=draw(ht_GOterm_sim)
# 
# 
# 
# # analyze age correlation only
# ggplot(trait_cor_subset, aes(x=wgcna_slot, y=cor, color=color, alpha=pval)) +
#   geom_point() +
#   scale_alpha_continuous(range=c(1,0)) +
#   guides(color="none")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# clusterProfiler::dotplot(compare_results)
# 
# WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_WGCNA_Inter.rds")
# WGCNAobj@misc$CD4T.Treg
# 
# df=WGCNAobj@assays$RNA$data[,grepl("CD4T\\.Treg",colnames(WGCNAobj@assays$RNA$data))]
# order_=c("<30-1","<30-2","31~47-1","31~47-2","31~47-3","31~47-4","48~68-1","48~68-2","48~68-3","48~68-4","48~68-5",">69-1",">69-2",">69-3",">69-4",">69-5",">69-6")
# match_colnames=colnames(df)[match(order_[sort(match(gsub(".*, |_.*","",colnames(df)), order_))], gsub(".*, |_.*","",colnames(df)))]
# df_reordered=df[,match_colnames]
# # take only the modules with sig. biological enrichment
# module_gene_list_split=module_gene_list_split[as.character(unique(compare_results@compareClusterResult$Cluster))]
# 
# slim_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_RepresentativeDonors_ForWGCNA.rds")
# slim_obj_subset=subset(slim_obj, Annot.inter=="CD4T.Treg")
# slim_obj_subset_percluster=list()
# for (i in 1:length(module_gene_list_split)) {
#   genes_=module_gene_list_split[[i]]
#   slim_obj_subset_percluster[[i]]=
#     df_reordered %>%
#     .[genes_, ] %>% 
#     as.data.frame()
#   slim_obj_subset_percluster[[i]]$gene=rownames(slim_obj_subset_percluster[[i]])
#   slim_obj_subset_percluster[[i]]=slim_obj_subset_percluster[[i]] %>%
#     tidyr::pivot_longer(cols=colnames(.)[!grepl("gene",colnames(.))], values_to="expr", names_to="cells") %>%
#     mutate(sex=gsub(".*\\#|, .*","",cells)) %>%
#     mutate(agecut=gsub(".*, |_.*","",cells)) %>%
#     group_by(agecut, gene) %>%
#     summarize_at("expr", mean) %>%
#     mutate(module=names(module_gene_list_split)[i])
# }
# slim_obj_subset_percluster_all=data.table::rbindlist(slim_obj_subset_percluster)
# mid_point=
#   slim_obj_subset_percluster_all %>%
#   group_by(agecut, module) %>%
#   summarize_at("expr", mean) %>%
#   rename(midpoint=expr)
# df_arranged=left_join(slim_obj_subset_percluster_all, mid_point, by=c("agecut","module")) %>%
#   mutate(expr.scaled=expr/midpoint)
# 
# df_arranged$agecut=forcats::fct_relevel(df_arranged$agecut, order_)
# ggplot(mid_point, aes(x=agecut, y=midpoint, color=module)) +
#   ggrastr::geom_point_rast() +
#   # facet_wrap(~module) +
#   geom_line(aes(group=module))
# 
# 
# 
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation_rough.pdf", width=9, height=5)
# plot(plot_)
# dev.off()
# 
# ### Plot with ggplot2
# all_cell_cor=WGCNAobj@misc$WGCNA$mt_cor$cor$all_cells; colnames(all_cell_cor)=paste0("M",c(1:ncol(all_cell_cor)))
# # switch sex and age
# all_cell_cor=all_cell_cor[c(2,1),]
# all_cell_fdr=WGCNAobj@misc$WGCNA$mt_cor$fdr$all_cells; colnames(all_cell_fdr)=paste0("M",c(1:ncol(all_cell_fdr)))
# all_cell_fdr=apply(all_cell_fdr, c(1,2), function(x) ifelse(x<0.0001, "****",
#                                                             ifelse(x>0.0001 & x<0.001, "***",
#                                                                    ifelse(x>0.001 & x<0.01, "**",
#                                                                           ifelse(x>0.01 & x<0.05, "*", NA)))))
# # switch sex and age
# all_cell_fdr=all_cell_fdr[c(2,1),]
# 
# library(ComplexHeatmap)
# col_fun=circlize::colorRamp2(c(-0.2, 0, 0.2), c("dodgerblue3", "white", "brown3"))
# ht=
# Heatmap(all_cell_cor,
#         name="mat", show_heatmap_legend=F,
#         col=col_fun,
#         cluster_rows=F, cluster_columns=F,
#         column_names_rot=0, column_names_centered=TRUE,
#         column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
# 
#         # layer_fun=function(j, i, x, y, width, height, fill) {
#         #   grid.text(sprintf("%.2f", pindex(as.matrix(all_cell_fdr), i, j)), x, y,
#         #             gp=gpar(fontsize=9, fontface="italic", fontfamily="HersheyScript"))}
# 
#         layer_fun=function(j, i, x, y, width, height, fill) {
#           grid.text(as.matrix(all_cell_fdr), x, y,
#                     gp=gpar(fontsize=10, fontface="italic"))}
# 
#         )
# lgd=Legend(col_fun=col_fun,
#            title=NULL,
#            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
#            direction="horizontal", title_position="lefttop", labels_gp=gpar(fontsize=9),
#            )
# ht_final=draw(ht, annotation_legend_list=lgd, annotation_legend_side="bottom")
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation_rough_ggplot.pdf", width=4, height=1)
# ht_final
# dev.off()
# 
# #####################################
# # 
# # 
# # 
# # ### Plot the WGCNA gene module-trait correlations of sex
# # #####################################
# # ###
# # 
# # library(dplyr)
# # library(hdWGCNA)
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # 
# # ### Get the correlation of sex
# # mt_cor=GetModuleTraitCorrelation(WGCNAobj)
# # Module_genes=GetModules(WGCNAobj)
# # 
# # cor_=lapply(1:length(mt_cor$cor), function(x) as.numeric(mt_cor$cor[[x]]["sex",]))
# # cor_df=as.data.frame(cor_)
# # colnames(cor_df)=names(mt_cor$cor)
# # rownames(cor_df)=Module_genes$module[match(colnames(mt_cor$cor[[1]]), Module_genes$color)]
# # cor_df$module=rownames(cor_df)
# # cor_df=tidyr::pivot_longer(cor_df, cols=-module, names_to="celltype", values_to="cor")
# # 
# # fdr_=lapply(1:length(mt_cor$fdr), function(x) as.numeric(mt_cor$fdr[[x]]["sex",]))
# # fdr_df=as.data.frame(fdr_)
# # colnames(fdr_df)=names(mt_cor$fdr)
# # rownames(fdr_df)=Module_genes$module[match(colnames(mt_cor$fdr[[1]]), Module_genes$color)]
# # fdr_df$module=rownames(fdr_df)
# # fdr_df=tidyr::pivot_longer(fdr_df, cols=-module, names_to="celltype", values_to="fdr")
# # 
# # pval_=lapply(1:length(mt_cor$pval), function(x) as.numeric(mt_cor$pval[[x]]["sex",]))
# # pval_df=as.data.frame(pval_)
# # colnames(pval_df)=names(mt_cor$pval)
# # rownames(pval_df)=Module_genes$module[match(colnames(mt_cor$pval[[1]]), Module_genes$color)]
# # pval_df$module=rownames(pval_df)
# # pval_df=tidyr::pivot_longer(pval_df, cols=-module, names_to="celltype", values_to="pval")
# # 
# # cor_fdr_pval_DF=merge(merge(cor_df, fdr_df, by=c("module","celltype")), pval_df, by=c("module","celltype"))
# # cor_fdr_pval_DF=cor_fdr_pval_DF %>% subset(celltype!="all_cells")
# # 
# # plot_sex=
# #   ggplot(cor_fdr_pval_DF, aes(x=cor, y=-log10(fdr), color=celltype)) +
# #   facet_wrap(~module, scales="free_y", ncol=2) +
# #   geom_point(size=2) +
# #   ggsci::scale_color_d3() +
# #   theme_classic() +
# #   labs(x="Pearson correlation with sex", y=expression(-log[10]~FDR), title=NULL) +
# #   theme(axis.text.x=element_text(size=9),
# #         axis.text.y=element_text(size=9),
# #         axis.title.x=element_text(size=10),
# #         axis.title.y=element_text(size=10),
# #         legend.position=c(0.75,0.15),
# #         legend.title=element_text(size=10),
# #         legend.text=element_text(size=9),
# #         title=element_text(size=11),
# #         strip.background=element_rect(fill="gray90")) +
# #   # coord_cartesian(ylim=c(1, 8)) +
# #   geom_hline(yintercept=-log10(0.05), linetype='dotted', linewidth=0.5) +
# #   geom_vline(xintercept=0, linetype='dotted', linewidth=0.5) +
# #   guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title="cell type", ncol=2))
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_rough_sex.pdf", height=4, width=5)
# # plot(plot_sex)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # 
# # ### Plot all the WGCNA gene module enrichments at the rough level
# # #####################################
# # ###
# # 
# # library(clusterProfiler)
# # library(ggplot2)
# # 
# # ### Get the module genes
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # gene_module=WGCNAobj@misc$WGCNA$wgcna_modules
# # gene_split=split(gene_module$gene_name, gene_module$module)
# # gene_split=gene_split[names(gene_split)!="grey"]
# # 
# # ### Enrich
# # MODULE_ENRICH_original=list()
# # for (i in 1:length(gene_split)) {
# #   genes_=gene_split[[i]]
# #   module_enrich=enrichGO(genes_,
# #                          keyType="SYMBOL",
# #                          OrgDb="org.Hs.eg.db",
# #                          ont="BP",
# #                          pvalueCutoff=0.05,
# #                          minGSSize=30,
# #                          maxGSSize=500)
# #   module_enrich_filter=gofilter(module_enrich, level=4)
# #   MODULE_ENRICH_original[[i]]=module_enrich
# # }
# # names(MODULE_ENRICH_original)=names(gene_split)
# # 
# # ### Plot
# # MODULE_ENRICH_Plots=list()
# # # get the top 10 terms for plotting
# # MODULE_ENRICH=lapply(MODULE_ENRICH_original, function(x) x@result[1:10,] %>% 
# #                        mutate(GeneRatio_div=as.numeric(strsplit(GeneRatio, split="/")[[1]][2])))
# # MODULE_ENRICH=lapply(MODULE_ENRICH, function(x) x %>% mutate(GeneRatio=Count/GeneRatio_div))
# # # plot with ggplot2
# # for (i in 1:length(MODULE_ENRICH)) {
# #   module_=MODULE_ENRICH[[i]] %>% arrange(GeneRatio) %>% select(Description, GeneRatio, p.adjust, Count)
# #   module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=50, collapse="\n"))
# #   module_$Description=forcats::fct_relevel(module_$Description, module_$Description)
# #   module_=rbind(module_, data.frame(Description=paste0(rep(" ",80), collapse=""), GeneRatio=NA, p.adjust=NA, Count=NA))
# #   x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10
# #   
# #   plot_=
# #     ggplot(module_, aes(x=GeneRatio, y=Description, color=p.adjust, size=Count)) +
# #     geom_point() +
# #     scale_colour_gradient(limits=quantile(module_$p.adjust, na.rm=T)[c(1,5)],
# #                           low="brown4", high="pink",
# #                           labels=sprintf(fmt="%0.01e", quantile(module_$p.adjust, na.rm=T)[c(1,5)]),
# #                           breaks=quantile(module_$p.adjust, na.rm=T)[c(1,5)]) +
# #     scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
# #                                    round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
# #                                    max(module_$Count, na.rm=T)),
# #                           labels=c(min(module_$Count, na.rm=T),
# #                                    round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
# #                                    max(module_$Count, na.rm=T)),
# #                           range=c(2,6)) +
# #     guides(size=guide_legend(order=1)) +
# #     labs(title=names(MODULE_ENRICH)[i], y=NULL) +
# #     theme_minimal() +
# #     theme(axis.text.x=element_text(size=9),
# #           axis.text.y=element_text(size=10),
# #           axis.title.x=element_text(size=10),
# #           axis.title.y=element_text(size=10),
# #           title=element_text(size=11),
# #           plot.title.position="plot",
# #           plot.title=element_text(hjust=0.5),
# #           legend.position="right",
# #           legend.box="vertical",
# #           legend.direction="vertical",
# #           legend.key.width=unit(0.5,"cm"),
# #           legend.margin=margin(t=20, b=20)
# #           ) +
# #       guides(color=guide_colorbar(label.theme=element_text(angle=90, size=9, vjust=0.5, hjust=0.5),
# #                                   title.theme=element_text(angle=90, size=10, vjust=1, hjust=0), order=1),
# #              size=guide_legend(order=2)) +
# #     scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(x_axis_expansion, x_axis_expansion)))
# #   
# #   MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
# # }
# # # # plot by clusterprofiler
# # # for (i in 1:length(gene_split)) {
# # #   if (MODULE_ENRICH[[i]]@result$p.adjust[1]<=0.05) {
# # #     x_axis_range=unlist(lapply(lapply(MODULE_ENRICH[[i]]@result$GeneRatio[1:10], function(x) as.numeric(strsplit(x, split="/")[[1]])), function(y) y[1]/y[2]))
# # #     x_axis_expansion=mean(range(x_axis_range))/10
# # #     
# # #     plot_=
# # #       dotplot(MODULE_ENRICH[[i]], showCategory=10,
# # #               label_format=50) +
# # #       scale_colour_gradient(limits=quantile(MODULE_ENRICH[[i]]@result$p.adjust[1:10])[c(1,5)],
# # #                             low="pink", high="brown4",
# # #                             labels=sprintf(fmt="%0.01e", quantile(MODULE_ENRICH[[i]]@result$p.adjust[1:10])[c(1,5)]),
# # #                             breaks=as.numeric(quantile(MODULE_ENRICH[[i]]@result$p.adjust[1:10])[c(1,5)])) +
# # #       scale_size_continuous(breaks=c(min(MODULE_ENRICH[[i]]@result$Count[1:10]),
# # #                                      round((min(MODULE_ENRICH[[i]]@result$Count[1:10])+max(MODULE_ENRICH[[i]]@result$Count[1:10]))/2),
# # #                                      max(MODULE_ENRICH[[i]]@result$Count[1:10])),
# # #                             labels=c(min(MODULE_ENRICH[[i]]@result$Count[1:10]),
# # #                                      round((min(MODULE_ENRICH[[i]]@result$Count[1:10])+max(MODULE_ENRICH[[i]]@result$Count[1:10]))/2),
# # #                                      max(MODULE_ENRICH[[i]]@result$Count[1:10])),
# # #                             range=c(2,6)) +
# # #       ggtitle(names(gene_split)[i]) +
# # #       theme(axis.text.x=element_text(size=9),
# # #             axis.text.y=element_text(size=10),
# # #             axis.title.x=element_text(size=10),
# # #             axis.title.y=element_text(size=10),
# # #             title=element_text(size=11),
# # #             legend.direction="horizontal",
# # #             legend.box="vertical",
# # #             legend.position="bottom") +
# # #       scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(x_axis_expansion, x_axis_expansion)))
# # #     
# # #     MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
# # #   }
# # # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_rough.pdf", height=4.5, width=5.5)
# # for (p in MODULE_ENRICH_Plots) plot(p)
# # dev.off()
# # 
# # ### Add ModuleTraitCorrelation info about celltypes for each module
# # library(hdWGCNA)
# # # get the info
# # mt_cor=GetModuleTraitCorrelation(WGCNAobj)
# # Module_genes=GetModules(WGCNAobj)
# # 
# # cor_=lapply(1:length(mt_cor$cor), function(x) as.numeric(mt_cor$cor[[x]]["age",]))
# # cor_df=as.data.frame(cor_)
# # colnames(cor_df)=names(mt_cor$cor)
# # rownames(cor_df)=Module_genes$module[match(colnames(mt_cor$cor[[1]]), Module_genes$color)]
# # cor_df$module=rownames(cor_df)
# # cor_df=tidyr::pivot_longer(cor_df, cols=-module, names_to="celltype", values_to="cor")
# # 
# # fdr_=lapply(1:length(mt_cor$fdr), function(x) as.numeric(mt_cor$fdr[[x]]["age",]))
# # fdr_df=as.data.frame(fdr_)
# # colnames(fdr_df)=names(mt_cor$fdr)
# # rownames(fdr_df)=Module_genes$module[match(colnames(mt_cor$fdr[[1]]), Module_genes$color)]
# # fdr_df$module=rownames(fdr_df)
# # fdr_df=tidyr::pivot_longer(fdr_df, cols=-module, names_to="celltype", values_to="fdr")
# # 
# # pval_=lapply(1:length(mt_cor$pval), function(x) as.numeric(mt_cor$pval[[x]]["age",]))
# # pval_df=as.data.frame(pval_)
# # colnames(pval_df)=names(mt_cor$pval)
# # rownames(pval_df)=Module_genes$module[match(colnames(mt_cor$pval[[1]]), Module_genes$color)]
# # pval_df$module=rownames(pval_df)
# # pval_df=tidyr::pivot_longer(pval_df, cols=-module, names_to="celltype", values_to="pval")
# # 
# # cor_fdr_pval_DF=merge(merge(cor_df, fdr_df, by=c("module","celltype")), pval_df, by=c("module","celltype"))
# # 
# # # clean the data
# # cor_fdr_pval_DF=cor_fdr_pval_DF %>% subset(celltype!="all_cells") %>% select(-pval)
# # plot_map_list=list()
# # for (i in 1:length(unique(cor_fdr_pval_DF$module))) {
# #   data_per_module=subset(cor_fdr_pval_DF, module==unique(cor_fdr_pval_DF$module)[i])
# #   plot_map_list[[i]]=
# #     ggplot(data_per_module, aes(x=celltype, y=cor)) +
# #     geom_segment(data=data_per_module,
# #                  aes(x=celltype, xend=celltype, y=0, yend=cor, alpha=fdr), linewidth=4, color="brown3") +
# #     theme_minimal() +
# #     labs(x=NULL, y="Pearson correlation \nwith age") +
# #     theme(axis.text.x=element_text(size=9),
# #           axis.text.y=element_text(size=9),
# #           axis.title.y=element_text(size=10),
# #           title=element_text(size=11),
# #           panel.grid.minor=element_blank(),
# #           panel.grid.major.x=element_blank(),
# #           legend.position="left",
# #           legend.title=element_text(size=9),
# #           legend.text=element_text(size=9)) +
# #     scale_y_continuous(breaks=seq(-1,1,0.2)) +
# #     scale_alpha_continuous(breaks=quantile(data_per_module$fdr)[c(1,3,5)],
# #                            
# #                            labels=sprintf(fmt="%0.1f", quantile(data_per_module$fdr)[c(1,3,5)]),
# #                            range=c(1,0.1)
# #                            ) +
# #     scale_x_discrete(guide=guide_axis(n.dodge=2)) +
# #     guides(alpha=guide_legend(title="FDR",
# #                               title.position="left", title.theme=element_text(angle=90, hjust=0.5)))
# # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_rough_MapToCelltypeRough.pdf", height=1.5, width=5.5)
# # for(i in 1:length(plot_map_list)) plot(plot_map_list[[i]])
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Plot the hub genes
# # #####################################
# # ###
# # 
# # ### Load
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # 
# # ### Mark the hub genes
# # marking=GetHubGenes(WGCNAobj, n_hubs=length(WGCNAobj@misc$WGCNA$wgcna_genes)) %>% subset(module!="grey") %>%
# #   group_by(module) %>% slice_max(kME, n=10) %>% as.data.frame() %>% mutate(marker=gene_name) %>% select(-c(module, kME))
# # hb=GetHubGenes(WGCNAobj, n_hubs=length(WGCNAobj@misc$WGCNA$wgcna_genes)) %>% subset(module!="grey") %>%
# #   dplyr::left_join(marking, by="gene_name")
# # hb$module=as.character(hb$module)
# # hb$module=forcats::fct_relevel(hb$module, unique(hb$module))
# # 
# # ### Plot each module
# # hb_list=
# #   hb %>% split(~module) %>%
# #   purrr::map(~ggplot(., aes(x=module, y=kME)) +
# #                geom_point(shape=1, size=1, color="brown4") +
# #                geom_label_repel(aes(label=marker), size=4, max.overlaps=20) +
# #                coord_flip() +
# #                labs(y="kME", x=NULL) +
# #                theme_minimal() +
# #                theme(panel.grid.major.x=element_blank(),
# #                      panel.grid.minor.x=element_blank(),
# #                      axis.title.x=element_text(size=12),
# #                      axis.text.y=element_text(size=12),
# #                      axis.text.x=element_text(size=10)))
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_rough_HubGenes.pdf", height=4, width=8)
# # for (i in 1:length(hb_list)) plot(hb_list[[i]])
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### EnrichmentNetwork with aPEAR : rough level
# # #####################################
# # ###
# # 
# # library(clusterProfiler)
# # library(hdWGCNA)
# # library(aPEAR)
# # library(dplyr)
# # library(ggplot2)
# # 
# # ### Extract the most sig. and correlated modules for each celltype at the rough level
# # WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # mt_cor=GetModuleTraitCorrelation(WGCNA_obj)
# # Module_genes=GetModules(WGCNA_obj)
# # 
# # # extract the cor, fdr, and pval
# # cor_=lapply(1:length(mt_cor$cor), function(x) as.numeric(mt_cor$cor[[x]]["age",]))
# # cor_df=as.data.frame(cor_)
# # colnames(cor_df)=names(mt_cor$cor)
# # rownames(cor_df)=Module_genes$module[match(colnames(mt_cor$cor[[1]]), Module_genes$color)]
# # cor_df$module=rownames(cor_df)
# # cor_df=tidyr::pivot_longer(cor_df, cols=-module, names_to="celltype", values_to="cor")
# # 
# # fdr_=lapply(1:length(mt_cor$fdr), function(x) as.numeric(mt_cor$fdr[[x]]["age",]))
# # fdr_df=as.data.frame(fdr_)
# # colnames(fdr_df)=names(mt_cor$fdr)
# # rownames(fdr_df)=Module_genes$module[match(colnames(mt_cor$fdr[[1]]), Module_genes$color)]
# # fdr_df$module=rownames(fdr_df)
# # fdr_df=tidyr::pivot_longer(fdr_df, cols=-module, names_to="celltype", values_to="fdr")
# # 
# # pval_=lapply(1:length(mt_cor$pval), function(x) as.numeric(mt_cor$pval[[x]]["age",]))
# # pval_df=as.data.frame(pval_)
# # colnames(pval_df)=names(mt_cor$pval)
# # rownames(pval_df)=Module_genes$module[match(colnames(mt_cor$pval[[1]]), Module_genes$color)]
# # pval_df$module=rownames(pval_df)
# # pval_df=tidyr::pivot_longer(pval_df, cols=-module, names_to="celltype", values_to="pval")
# # 
# # cor_fdr_pval_DF=merge(merge(cor_df, fdr_df, by=c("module","celltype")), pval_df, by=c("module","celltype"))
# # 
# # # sort sig. modules and the most cor. ones for each celltype
# # cor_fdr_pval_DF_pos=cor_fdr_pval_DF %>%
# #   subset(celltype!="all_cells" & fdr<=0.05 & cor>0) %>%
# #   group_by(celltype) %>%
# #   slice_max(cor, n=1)
# # module_celltype=split(cor_fdr_pval_DF_pos$celltype, cor_fdr_pval_DF_pos$module)
# # module_celltype_df_pos=data.frame(module=names(module_celltype),
# #                                   celltype=unlist(lapply(module_celltype, function(x) paste0(x, collapse=", "))))
# # celltype_m_pos=list()
# # for (i in 1:nrow(module_celltype_df_pos)) {
# #   celltype_m_pos[[i]]=subset(Module_genes, module==cor_fdr_pval_DF_pos$module[i])$gene_name
# # }
# # names(celltype_m_pos)=paste0(module_celltype_df_pos$module, " - ", module_celltype_df_pos$celltype)
# # # * CD8 T cell - pos: M3
# # 
# # cor_fdr_pval_DF_neg=cor_fdr_pval_DF %>%
# #   subset(celltype!="all_cells" & fdr<=0.05 & cor<0) %>%
# #   group_by(celltype) %>%
# #   slice_min(cor, n=1)
# # module_celltype=split(cor_fdr_pval_DF_neg$celltype, cor_fdr_pval_DF_neg$module)
# # module_celltype_df_neg=data.frame(module=names(module_celltype),
# #                                   celltype=unlist(lapply(module_celltype, function(x) paste0(x, collapse=", "))))
# # celltype_m_neg=list()
# # for (i in 1:nrow(module_celltype_df_neg)) {
# #   celltype_m_neg[[i]]=subset(Module_genes, module==cor_fdr_pval_DF_neg$module[i])$gene_name
# # }
# # names(celltype_m_neg)=paste0(module_celltype_df_neg$module, " - ", module_celltype_df_neg$celltype)
# # # * CD8 T cell - neg: M4
# # 
# # ### GO enrichment
# # EGO_pos=plotlist_pos=list()
# # for (i in 1:length(celltype_m_pos)) {
# #   genes_=celltype_m_pos[[i]]
# #   EGO_pos[[i]]=enrichGO(gene=genes_,
# #                         OrgDb="org.Hs.eg.db",
# #                         keyType="SYMBOL",
# #                         ont="BP",
# #                         pvalueCutoff=0.05,
# #                         readable=TRUE)
# #   result=gofilter(EGO_pos[[i]], level=3)@result
# #   message(paste0("EGO ",i," of ", length(celltype_m_pos), " has been enriched."))
# #   
# #   source("~/Rscripts/aPEAR2.R")
# #   theme_=aPEAR.theme
# #   theme_$colorBy="p.adjust"; theme_$nodeSize="Count"; theme_$fontSize=4; theme_$repelLabels=TRUE
# #   enrichs=
# #     enrichmentNetwork2(result,
# #                        title=names(celltype_m_neg),
# #                        minClusterSize=3,
# #                        theme=theme_,
# #                        plotOnly=FALSE,
# #                        palette="Blues")
# #   plotlist_pos[[i]]=enrichs
# # }
# # 
# # EGO_neg=plotlist_neg=list()
# # for (i in 1:length(celltype_m_neg)) {
# #   genes_=celltype_m_neg[[i]]
# #   EGO_neg[[i]]=enrichGO(gene=genes_,
# #                         OrgDb="org.Hs.eg.db",
# #                         keyType="SYMBOL",
# #                         ont="BP",
# #                         pvalueCutoff=0.05,
# #                         readable=TRUE)
# #   result=gofilter(EGO_neg[[i]], level=3)@result
# #   message(paste0("EGO ",i," of ", length(celltype_m_neg), " has been enriched."))
# #   
# #   source("~/Rscripts/aPEAR2.R")
# #   theme_=aPEAR.theme
# #   theme_$colorBy="p.adjust"; theme_$nodeSize="Count"; theme_$fontSize=4; theme_$repelLabels=TRUE
# #   enrichs=enrichmentNetwork2(result,
# #                              title=names(celltype_m_neg),
# #                              minClusterSize=3,
# #                              theme=theme_,
# #                              plotOnly=FALSE,
# #                              palette="Blues")
# #   plotlist_neg[[i]]=enrichs
# # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_rough.pdf", height=5, width=5)
# # for (p in plotlist_pos) plot(p)
# # for (p in plotlist_neg) plot(p)
# # dev.off()
# # 
# # # detailed terms in CD8 T cells
# # EGO_pos_CD8T=enrichGO(gene=celltype_m_pos[grepl("CD8T cells",names(celltype_m_pos))][[1]],
# #                       OrgDb="org.Hs.eg.db",
# #                       keyType="SYMBOL",
# #                       ont="BP",
# #                       pvalueCutoff=0.05,
# #                       readable=TRUE)
# # result_CD8T=gofilter(EGO_pos_CD8T, level=4)@result
# # enrichs_CD8T=enrichmentNetwork(result_CD8T,
# #                                minClusterSize=5,
# #                                colorBy="p.adjust",
# #                                nodeSize="Count",
# #                                fontsize=5,
# #                                repelLabels=TRUE,
# #                                plotOnly=FALSE)
# # plotlist_pos_CD8T=enrichs_CD8T+ggtitle("CD8T cells : pos.assoc.w/ ages")
# # 
# # EGO_neg_CD8T=enrichGO(gene=celltype_m_neg[grepl("CD8T cells",names(celltype_m_neg))][[1]],
# #                       OrgDb="org.Hs.eg.db",
# #                       keyType="SYMBOL",
# #                       ont="BP",
# #                       pvalueCutoff=0.05,
# #                       readable=TRUE)
# # result_CD8T=gofilter(EGO_neg_CD8T, level=4)@result
# # enrichs_CD8T=enrichmentNetwork(result_CD8T,
# #                                minClusterSize=5,
# #                                colorBy="p.adjust",
# #                                nodeSize="Count",
# #                                fontsize=5,
# #                                repelLabels=TRUE,
# #                                plotOnly=FALSE)
# # plotlist_neg_CD8T=enrichs_CD8T+ggtitle("CD8T cells : neg.assoc.w/ ages")
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_CD8T_rough.pdf")
# # plot(plotlist_pos_CD8T)
# # plot(plotlist_neg_CD8T)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Get the hub genes from each module : rough level
# # #####################################
# # ###
# # 
# # # * actually we extract all the genes but arrange them by KME so that the top_n will be the hub genes for each module
# # 
# # ### Load and extract
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # 
# # # get the list of modules:
# # modules=GetModules(WGCNAobj)
# # mods=levels(modules$module)
# # mods=mods[mods!='grey']
# # hb=GetHubGenes(WGCNAobj, n_hubs=length(WGCNAobj@misc$WGCNA$wgcna_genes))
# # 
# # # hubgene network
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/HubGene_network_Rough.pdf", height=12, width=12)
# # 
# # allModule_network=
# #   HubGeneNetworkPlot(
# #     WGCNAobj,
# #     n_hubs=5,
# #     n_other=10,
# #     edge_prop=1,
# #     mods="all"
# #   ) +
# #   title("All Modules\n")
# # 
# # CD8T_posM3_network=
# #   HubGeneNetworkPlot(
# #     WGCNAobj,
# #     n_hubs=10,
# #     n_other=50,
# #     edge_prop=1,
# #     mods=mods[3]
# #   ) +
# #   title("CD8T cells: M3 - pos.assoc.w/ ages\n")
# # 
# # CD8T_negM4_network=
# #   HubGeneNetworkPlot(
# #     WGCNAobj,
# #     n_hubs=10,
# #     n_other=50,
# #     edge_prop=1,
# #     mods=mods[4]
# #   ) +
# #   title("CD8T cells: M4 - neg.assoc.w/ ages\n")
# # 
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # 
# # 
# # 
# # ##########################################################################
# # #################### WGCNA and aPEAR at the inter level ####################
# # ##########################################################################
# # 
# # 
# # ### WGCNA analysis at the intermediate level
# # #####################################
# # ###
# # ###
# # library(dplyr)
# # library(WGCNA)
# # library(hdWGCNA)
# # library(Seurat)
# # 
# # # load
# # THE_WGCNA=readRDS("~/Project_PBMCage/Tempt_RDS/MetaCellObj_Inter.rds")
# # 
# # # remove the celltype with <=4 observations
# # celltype_rm=names(table(THE_WGCNA$Annot.inter))[table(THE_WGCNA$Annot.inter)<=4]
# # THE_WGCNA=subset(THE_WGCNA, (Annot.inter %in% celltype_rm)==FALSE)
# # 
# # # take the features
# # WGCNAobj=SetupForWGCNA(
# #   THE_WGCNA,
# #   group.by="Annot.inter",
# #   gene_select="fraction",
# #   fraction=0.05,
# #   wgcna_name="WGCNA")
# # length(WGCNAobj@misc$WGCNA$wgcna_genes)
# # 
# # # set the expression data
# # WGCNAobj=SetDatExpr(
# #   WGCNAobj,
# #   group_name=names(table(WGCNAobj$Annot.inter)),
# #   group.by="Annot.inter",
# #   assay="RNA",
# #   use_metacells=TRUE,
# #   slot="data"
# # )
# # saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # 
# # ### Choose the power
# # selected_power=0
# # WGCNAobj=TestSoftPowers(
# #   WGCNAobj,
# #   powers=c(seq(1, 10, by=1), seq(12, 30, by=2)))
# # plot_list=PlotSoftPowers(WGCNAobj, point_size=5, text_size=3)
# # possible_power=c()
# # for (i in 1:length(plot_list)) {
# #   tryCatch({
# #     possible_power[i]=(plot_list[[i]][[1]] %>% subset(text_color=="white"))$Power
# #   }, error=function(msg) {print("No Power found.")})
# # }
# # if (length(unique(possible_power))==1) {selected_power=unique(possible_power)} else {print("Error.")}
# # print(selected_power)
# # 
# # ### Construct the network
# # if (selected_power!=0) {
# #   WGCNAobj=ConstructNetwork(
# #     WGCNAobj,
# #     soft_power=selected_power,
# #     corType="pearson",
# #     networkType="signed",
# #     TOMType="signed",
# #     overwrite_tom=TRUE,
# #     tom_outdir="TOM",
# #     tom_name="inter"
# #   )
# #   modules=levels(WGCNAobj@misc[["WGCNA"]][["wgcna_modules"]][["module"]])
# #   saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # } else {print("No selected_power found.")}
# # 
# # ### Compute MEs
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # WGCNAobj[["RNA"]]=as(object=WGCNAobj[["RNA"]], Class="Assay")
# # modules=levels(WGCNAobj@misc[["WGCNA"]][["wgcna_modules"]][["module"]])
# # if (length(modules)!=1 & length(modules)!=2) {
# #   WGCNAobj=ScaleData(WGCNAobj)
# #   WGCNAobj=ModuleEigengenes(WGCNAobj, group.by.vars="donor_id")
# #   MEs=GetMEs(WGCNAobj)
# #   saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # } else {print("Fewer than 2 modules were found.")}
# # 
# # ### Associate gene modules with traits
# # if (length(modules)!=1 & length(modules)!=2) {
# #   WGCNAobj$sex=factor(WGCNAobj$sex, levels=c("male","female"))
# #   cur_traits=c("age","sex")
# #   WGCNAobj=ModuleTraitCorrelation(
# #     WGCNAobj,
# #     traits=cur_traits,
# #     features="hMEs",
# #     cor_method="pearson",
# #     group.by='Annot.inter'
# #   )
# #   saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # } else {print("Fewer than 2 modules were found.")}
# # 
# # ### Calculate module connectivity
# # if (length(modules)!=1 & length(modules)!=2) {
# #   WGCNAobj=ModuleConnectivity(
# #     WGCNAobj,
# #     group.by='Annot.inter',
# #     corFnc="bicor",
# #     corOptions="use='p'",
# #     harmonized=FALSE,
# #     assay=NULL,
# #     slot="data",
# #     group_name=names(table(WGCNAobj$Annot.inter)))
# #   # rename modules to replace color names
# #   WGCNAobj=ResetModuleNames(WGCNAobj, new_name="M")
# #   saveRDS(WGCNAobj, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # 
# #   Module_genes=GetModules(WGCNAobj)
# #   saveRDS(Module_genes, "~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter_ModuleGenes.rds")
# # }
# # 
# # #####################################
# # 
# # 
# # 
# # ### Plot the WGCNA gene module-trait correlations at the inter level
# # #####################################
# # ###
# # 
# # ###
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # 
# # plot_=PlotModuleTraitCorrelation(
# #   WGCNAobj,
# #   label='fdr', # or add pval label in each cell of the heatmap
# #   label_symbol='numeric', # labels as 'stars' or as 'numeric'
# #   text_size=3,
# #   text_digits=2,
# #   text_color='black',
# #   high_color='#fc9272',
# #   mid_color='#ffffbf',
# #   low_color='#9ecae1',
# #   plot_max=0.2,
# #   combine=T # 合并结果
# # )
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation_inter.pdf", width=9, height=12)
# # plot(plot_)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Compare the WGCNA gene modules at the rough and the inter level
# # #####################################
# # ###
# # 
# # ### Prepare the jaccard similarity function
# # jaccard=function(a, b) {
# #   intersection=length(intersect(a, b))
# #   union=length(a)+length(b)-intersection
# #   return (intersection/union)
# # }
# # 
# # ### Load the data
# # rough_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
# # inter_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # 
# # rough_module=rough_obj@misc$WGCNA$wgcna_modules %>% subset(module!="grey"); rough_module$module=as.character(rough_module$module)
# # rough_module_levels=sort(unique(rough_module$module))
# # inter_module=inter_obj@misc$WGCNA$wgcna_modules %>% subset(module!="grey"); inter_module$module=as.character(inter_module$module)
# # inter_module_levels=sort(unique(inter_module$module))
# # 
# # ### Compare each module
# # Jaccard_SIM=Jaccard_SIM_Name=c()
# # for (i in 1:length(rough_module_levels)) {
# #   the_rough_module=rough_module_levels[i]
# #   the_rough_genes_=rough_module %>% subset(module==the_rough_module) %>% select(gene_name) %>% unlist() %>% unname()
# #   
# #   for (j in 1:length(inter_module_levels)) {
# #     the_inter_module=inter_module_levels[j]
# #     the_inter_genes=inter_module %>% subset(module==the_inter_module) %>% select(gene_name) %>% unlist() %>% unname()
# #     
# #     jaccard_similarity=jaccard(the_rough_genes_, the_inter_genes)
# #     Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
# #     Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0("Rough_",i," Inter_",j))
# #   }
# # }
# # names(Jaccard_SIM)=Jaccard_SIM_Name
# # # jaccard similarity shows that the modules at the rough and inter levels are consistent
# # 
# # ### Plot the jaccard similarity
# # Jaccard_SIM_df=data.frame(Pair_1=gsub(" .*","",names(Jaccard_SIM)), Pair_2=gsub(".* ","",names(Jaccard_SIM)), jaccard=Jaccard_SIM)
# # Jaccard_SIM_df=Jaccard_SIM_df %>% subset(jaccard!=0) %>% mutate(comparison=paste0("M",gsub(".*_","",Pair_1)))
# # jaccard_plot=
# #   ggplot(Jaccard_SIM_df, aes(x=comparison, y=jaccard, fill=comparison, color=comparison)) +
# #   geom_bar(stat="identity", width=0.8, alpha=0.2) +
# #   coord_flip() +
# #   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20c_d3")[c(1,3,2,4,5)]) +
# #   scale_fill_manual(values=paletteer::paletteer_d("ggsci::category20c_d3")[c(1,3,2,4,5)]) +
# #   theme_classic() +
# #   labs(x=NULL, y="Jaccard similarity index", title=NULL) +
# #   theme(axis.text.x=element_text(size=10),
# #         axis.text.y=element_text(size=9),
# #         axis.title.x=element_text(size=10),
# #         axis.title.y=element_text(size=10),
# #         legend.position="none",
# #         # legend.title=element_blank(),
# #         title=element_text(size=10)) +
# #   geom_hline(yintercept=1, linetype="dashed", linewidth=0.4)
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_rough_vsInter_JaccardSimilarity.pdf", height=2, width=4.5)
# # plot(jaccard_plot)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Plot all the WGCNA gene module enrichments at the inter level
# # #####################################
# # ###
# # 
# # library(clusterProfiler)
# # 
# # # get the module genes
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # gene_module=WGCNAobj@misc$WGCNA$wgcna_modules
# # gene_split=split(gene_module$gene_name, gene_module$module)
# # gene_split=gene_split[names(gene_split)!="grey"]
# # 
# # # enrich
# # MODULE_ENRICH=list()
# # for (i in 1:length(gene_split)) {
# #   genes_=gene_split[[i]]
# #   module_enrich=enrichGO(genes_,
# #                          keyType="SYMBOL",
# #                          OrgDb="org.Hs.eg.db",
# #                          ont="BP",
# #                          pvalueCutoff=0.05,
# #                          minGSSize=30,
# #                          maxGSSize=500)
# #   module_enrich_filter=gofilter(module_enrich, level=4)
# #   MODULE_ENRICH[[i]]=module_enrich
# # }
# # names(MODULE_ENRICH)=names(gene_split)
# # 
# # # plot
# # MODULE_ENRICH_Plots=list()
# # for (i in 1:length(gene_split)) {
# #   if (MODULE_ENRICH[[i]]@result$p.adjust[1]<=0.05) {
# #     if (nrow(subset(MODULE_ENRICH[[i]]@result, p.adjust<=0.05))<3) {
# #       x_axis_range=unlist(lapply(lapply(MODULE_ENRICH[[i]]@result$GeneRatio, function(x) as.numeric(strsplit(x, split="/")[[1]])), function(y) y[1]/y[2]))
# #       x_axis_expansion=mean(range(x_axis_range))/10
# #       plot_=
# #         dotplot(MODULE_ENRICH[[i]], showCategory=10,
# #                 label_format=50) +
# #         scale_colour_gradient(low="pink", high="brown4",
# #                               labels= ~sprintf(fmt="%0.01e", .)) +
# #         scale_size_continuous(range=c(2,6)) +
# #         ggtitle(names(gene_split)[i]) +
# #         theme(axis.text.x=element_text(size=8),
# #               axis.text.y=element_text(size=8),
# #               axis.title.x=element_text(size=9),
# #               axis.title.y=element_text(size=9),
# #               title=element_text(size=11)) +
# #         scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(x_axis_expansion, x_axis_expansion)))
# #     } else {
# #       x_axis_range=unlist(lapply(lapply(MODULE_ENRICH[[i]]@result$GeneRatio[1:10], function(x) as.numeric(strsplit(x, split="/")[[1]])), function(y) y[1]/y[2]))
# #       x_axis_expansion=mean(range(x_axis_range))/10
# #       plot_=
# #         dotplot(MODULE_ENRICH[[i]], showCategory=10,
# #                 label_format=50) +
# #         scale_colour_gradient(limits=quantile(MODULE_ENRICH[[i]]@result$p.adjust[1:10])[c(1,5)],
# #                               low="pink", high="brown4",
# #                               labels= ~sprintf(fmt="%0.01e", .)) +
# #         scale_size_continuous(breaks=c(min(MODULE_ENRICH[[i]]@result$Count[1:10]),
# #                                        round((min(MODULE_ENRICH[[i]]@result$Count[1:10])+max(MODULE_ENRICH[[i]]@result$Count[1:10]))/2),
# #                                        max(MODULE_ENRICH[[i]]@result$Count[1:10])),
# #                               labels=c(min(MODULE_ENRICH[[i]]@result$Count[1:10]),
# #                                        round((min(MODULE_ENRICH[[i]]@result$Count[1:10])+max(MODULE_ENRICH[[i]]@result$Count[1:10]))/2),
# #                                        max(MODULE_ENRICH[[i]]@result$Count[1:10])),
# #                               range=c(2,6)) +
# #         ggtitle(names(gene_split)[i]) +
# #         theme(axis.text.x=element_text(size=8),
# #               axis.text.y=element_text(size=8),
# #               axis.title.x=element_text(size=9),
# #               axis.title.y=element_text(size=9),
# #               title=element_text(size=11)) +
# #         scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(x_axis_expansion, x_axis_expansion)))
# #     }
# # 
# #     MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
# #   }
# # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleEnrichment_inter.pdf", height=5, width=6)
# # for (p in MODULE_ENRICH_Plots) plot(p)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### EnrichmentNetwork with aPEAR : inter level
# # #####################################
# # ###
# # 
# # library(clusterProfiler)
# # library(hdWGCNA)
# # library(aPEAR)
# # library(dplyr)
# # library(ggplot2)
# # 
# # ### Extract the most sig. and correlated modules for each celltype at the inter level
# # WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # mt_cor=GetModuleTraitCorrelation(WGCNA_obj)
# # Module_genes=GetModules(WGCNA_obj)
# # 
# # # extract the cor, fdr, and pval
# # cor_=lapply(1:length(mt_cor$cor), function(x) as.numeric(mt_cor$cor[[x]]["age",]))
# # cor_df=as.data.frame(cor_)
# # colnames(cor_df)=names(mt_cor$cor)
# # rownames(cor_df)=Module_genes$module[match(colnames(mt_cor$cor[[1]]), Module_genes$color)]
# # cor_df$module=rownames(cor_df)
# # cor_df=tidyr::pivot_longer(cor_df, cols=-module, names_to="celltype", values_to="cor")
# # 
# # fdr_=lapply(1:length(mt_cor$fdr), function(x) as.numeric(mt_cor$fdr[[x]]["age",]))
# # fdr_df=as.data.frame(fdr_)
# # colnames(fdr_df)=names(mt_cor$fdr)
# # rownames(fdr_df)=Module_genes$module[match(colnames(mt_cor$fdr[[1]]), Module_genes$color)]
# # fdr_df$module=rownames(fdr_df)
# # fdr_df=tidyr::pivot_longer(fdr_df, cols=-module, names_to="celltype", values_to="fdr")
# # 
# # pval_=lapply(1:length(mt_cor$pval), function(x) as.numeric(mt_cor$pval[[x]]["age",]))
# # pval_df=as.data.frame(pval_)
# # colnames(pval_df)=names(mt_cor$pval)
# # rownames(pval_df)=Module_genes$module[match(colnames(mt_cor$pval[[1]]), Module_genes$color)]
# # pval_df$module=rownames(pval_df)
# # pval_df=tidyr::pivot_longer(pval_df, cols=-module, names_to="celltype", values_to="pval")
# # 
# # cor_fdr_pval_DF=merge(merge(cor_df, fdr_df, by=c("module","celltype")), pval_df, by=c("module","celltype"))
# # 
# # # sort sig. modules and the most cor. ones for each celltype
# # cor_fdr_pval_DF_pos=cor_fdr_pval_DF %>%
# #   subset(celltype!="all_cells" & fdr<=0.05 & cor>0) %>%
# #   group_by(celltype) %>%
# #   slice_max(cor, n=1)
# # module_celltype=split(cor_fdr_pval_DF_pos$celltype, cor_fdr_pval_DF_pos$module)
# # module_celltype_df_pos=data.frame(module=names(module_celltype),
# #                                   celltype=unlist(lapply(module_celltype, function(x) paste0(x, collapse=", "))))
# # celltype_m_pos=list()
# # for (i in 1:nrow(module_celltype_df_pos)) {
# #   celltype_m_pos[[i]]=subset(Module_genes, module==cor_fdr_pval_DF_pos$module[i])$gene_name
# # }
# # names(celltype_m_pos)=paste0(module_celltype_df_pos$module, " - ", module_celltype_df_pos$celltype)
# # # * CD8T.Tcm - pos: M4
# # 
# # cor_fdr_pval_DF_neg=cor_fdr_pval_DF %>%
# #   subset(celltype!="all_cells" & fdr<=0.05 & cor<0) %>%
# #   group_by(celltype) %>%
# #   slice_min(cor, n=1)
# # module_celltype=split(cor_fdr_pval_DF_neg$celltype, cor_fdr_pval_DF_neg$module)
# # module_celltype_df_neg=data.frame(module=names(module_celltype),
# #                                   celltype=unlist(lapply(module_celltype, function(x) paste0(x, collapse=", "))))
# # celltype_m_neg=list()
# # for (i in 1:nrow(module_celltype_df_neg)) {
# #   celltype_m_neg[[i]]=subset(Module_genes, module==cor_fdr_pval_DF_neg$module[i])$gene_name
# # }
# # names(celltype_m_neg)=paste0(module_celltype_df_neg$module, " - ", module_celltype_df_neg$celltype)
# # 
# # ### GO enrichment
# # EGO_pos=plotlist_pos=list()
# # for (i in 1:length(celltype_m_pos)) {
# #   genes_=celltype_m_pos[[i]]
# #   EGO_pos[[i]]=enrichGO(gene=genes_,
# #                         OrgDb="org.Hs.eg.db",
# #                         keyType="SYMBOL",
# #                         ont="BP",
# #                         pvalueCutoff=0.05,
# #                         readable=TRUE)
# #   result=gofilter(EGO_pos[[i]], level=3)@result
# #   message(paste0("EGO ",i," of ", length(celltype_m_pos), " has been enriched."))
# #   enrichs=enrichmentNetwork(result,
# #                             minClusterSize=5,
# #                             colorBy="p.adjust",
# #                             nodeSize="Count",
# #                             fontsize=5,
# #                             repelLabels=TRUE,
# #                             plotOnly=FALSE)
# #   plotlist_pos[[i]]=enrichs+ggtitle(paste0(gsub("^M[0-9] - ","",names(celltype_m_pos)[i]), " : pos.assoc.w/ ages\n\n"))
# # }
# # 
# # EGO_neg=plotlist_neg=list()
# # for (i in 1:length(celltype_m_neg)) {
# #   genes_=celltype_m_neg[[i]]
# #   EGO_neg[[i]]=enrichGO(gene=genes_,
# #                         OrgDb="org.Hs.eg.db",
# #                         keyType="SYMBOL",
# #                         ont="BP",
# #                         pvalueCutoff=0.05,
# #                         readable=TRUE)
# #   result=gofilter(EGO_neg[[i]], level=3)@result
# #   message(paste0("EGO ",i," of ", length(celltype_m_neg), " has been enriched."))
# #   enrichs=enrichmentNetwork(result,
# #                             minClusterSize=5,
# #                             colorBy="p.adjust",
# #                             nodeSize="Count",
# #                             fontsize=5,
# #                             repelLabels=TRUE,
# #                             plotOnly=FALSE)
# #   plotlist_neg[[i]]=enrichs+ggtitle(paste0(gsub("^M[0-9] - ","",names(celltype_m_neg)[i]), " : neg.assoc.w/ ages\n\n"))
# # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_inter.pdf", height=5, width=15)
# # for (p in plotlist_pos) plot(p)
# # for (p in plotlist_neg) plot(p)
# # dev.off()
# # 
# # # detailed terms in CD8 T cells
# # EGO_pos_CD8T=enrichGO(gene=celltype_m_pos[grepl("CD8T\\.",names(celltype_m_pos))][[1]],
# #                       OrgDb="org.Hs.eg.db",
# #                       keyType="SYMBOL",
# #                       ont="BP",
# #                       pvalueCutoff=0.05,
# #                       readable=TRUE)
# # result_CD8T=gofilter(EGO_pos_CD8T, level=4)@result
# # # since markov does not give results, so choose other cluster methods: "hier" or "spectral"
# # settings=aPEAR.methods
# # settings$cluster="hier"
# # enrichs_CD8T=enrichmentNetwork(result_CD8T,
# #                                minClusterSize=5,
# #                                colorBy="p.adjust",
# #                                nodeSize="Count",
# #                                fontsize=5,
# #                                repelLabels=TRUE,
# #                                plotOnly=FALSE, methods=settings) 
# # plotlist_pos_CD8T=enrichs_CD8T+ggtitle("CD8T.Tcm cells : pos.assoc.w/ ages")
# # 
# # # !! No module negatively associated with ages in CD8 T subtypes
# # # EGO_neg_CD8T=enrichGO(gene=celltype_m_neg[grepl("CD8T\\.",names(celltype_m_neg))][[1]],
# # #                       OrgDb="org.Hs.eg.db",
# # #                       keyType="SYMBOL",
# # #                       ont="BP",
# # #                       pvalueCutoff=0.05,
# # #                       readable=TRUE)
# # # result_CD8T=gofilter(EGO_neg_CD8T, level=4)@result
# # # enrichs_CD8T=enrichmentNetwork(result_CD8T,
# # #                                minClusterSize=5,
# # #                                colorBy="p.adjust",
# # #                                nodeSize="Count",
# # #                                fontsize=5,
# # #                                repelLabels=TRUE,
# # #                                plotOnly=FALSE)
# # # plotlist_neg_CD8T=enrichs_CD8T+ggtitle("CD8T? cells : neg.assoc.w/ ages")
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_CD8T_inter.pdf")
# # plot(plotlist_pos_CD8T)
# # # plot(plotlist_neg_CD8T)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### EnrichmentNetwork with aPEAR : inter level - for each module
# # #####################################
# # ###
# # 
# # library(clusterProfiler)
# # library(hdWGCNA)
# # library(aPEAR)
# # library(dplyr)
# # library(ggplot2)
# # 
# # ### Extract the modules
# # WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # mt_cor=GetModuleTraitCorrelation(WGCNA_obj)
# # Module_genes=GetModules(WGCNA_obj) %>% subset(module!="grey")
# # Module_genes$module=as.character(Module_genes$module)
# # genes_per_module=split(Module_genes$gene_name, Module_genes$module)
# # # set the palettes
# # palettes_=c("Blues","Greens","Oranges","Purples","Greys")
# # 
# # ### GO enrichment
# # EGO=plotlist=list()
# # for (i in 1:length(genes_per_module)) {
# #   genes_=genes_per_module[i][[1]]
# #   EGO[[i]]=enrichGO(gene=genes_,
# #                     OrgDb="org.Hs.eg.db",
# #                     keyType="SYMBOL",
# #                     ont="BP",
# #                     pvalueCutoff=0.05,
# #                     readable=TRUE)
# #   result=gofilter(EGO[[i]], level=3)@result
# #   message(paste0("EGO ",i," of ", length(genes_per_module), " has been enriched."))
# #   
# #   source("~/Rscripts/aPEAR2.R")
# #   theme_=aPEAR.theme
# #   theme_$colorBy="p.adjust"; theme_$nodeSize="Count"; theme_$fontSize=4; theme_$repelLabels=TRUE
# #   enrichs=
# #     enrichmentNetwork2(result,
# #                        title=paste0(names(genes_per_module)[i], " - "),
# #                        minClusterSize=3,
# #                        theme=theme_,
# #                        plotOnly=FALSE,
# #                        palette=palettes_[i])
# #   plotlist[[i]]=enrichs
# # }
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_moduleAPEAR_inter_AllModules.pdf", height=5, width=6)
# # for (p in plotlist) plot(p)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Connect the modules found at the inter level and the celltypes
# # #####################################
# # ###
# # 
# # library(hdWGCNA)
# # library(dplyr)
# # library(ggplot2)
# # 
# # ### Extract the most sig. and correlated modules for each celltype at the inter level
# # WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # mt_cor=GetModuleTraitCorrelation(WGCNA_obj)
# # Module_genes=GetModules(WGCNA_obj)
# # 
# # # extract the cor, fdr, and pval
# # cor_=lapply(1:length(mt_cor$cor), function(x) as.numeric(mt_cor$cor[[x]]["age",]))
# # cor_df=as.data.frame(cor_)
# # colnames(cor_df)=names(mt_cor$cor)
# # rownames(cor_df)=Module_genes$module[match(colnames(mt_cor$cor[[1]]), Module_genes$color)]
# # cor_df$module=rownames(cor_df)
# # cor_df=tidyr::pivot_longer(cor_df, cols=-module, names_to="celltype", values_to="cor")
# # 
# # fdr_=lapply(1:length(mt_cor$fdr), function(x) as.numeric(mt_cor$fdr[[x]]["age",]))
# # fdr_df=as.data.frame(fdr_)
# # colnames(fdr_df)=names(mt_cor$fdr)
# # rownames(fdr_df)=Module_genes$module[match(colnames(mt_cor$fdr[[1]]), Module_genes$color)]
# # fdr_df$module=rownames(fdr_df)
# # fdr_df=tidyr::pivot_longer(fdr_df, cols=-module, names_to="celltype", values_to="fdr")
# # 
# # pval_=lapply(1:length(mt_cor$pval), function(x) as.numeric(mt_cor$pval[[x]]["age",]))
# # pval_df=as.data.frame(pval_)
# # colnames(pval_df)=names(mt_cor$pval)
# # rownames(pval_df)=Module_genes$module[match(colnames(mt_cor$pval[[1]]), Module_genes$color)]
# # pval_df$module=rownames(pval_df)
# # pval_df=tidyr::pivot_longer(pval_df, cols=-module, names_to="celltype", values_to="pval")
# # 
# # cor_fdr_pval_DF=merge(merge(cor_df, fdr_df, by=c("module","celltype")), pval_df, by=c("module","celltype"))
# # cor_fdr_pval_DF=cor_fdr_pval_DF %>% subset(celltype!="all_cells" & fdr!=0 & !is.na(fdr))
# # 
# # ### Plot
# # library(ggalluvial)
# # 
# # cor_fdr_pval_DF_poscorr=cor_fdr_pval_DF %>% subset(cor>0) %>% mutate(correlation=abs(cor))
# # # take only those celltypes with a sum of correlation>0.15 with age in all the modules
# # cor_fdr_pval_DF_poscorr_filter=cor_fdr_pval_DF_poscorr %>% 
# #   group_by(celltype) %>% 
# #   summarise_at("correlation", sum) %>% 
# #   subset(correlation>=0.15) %>%
# #   select(celltype) %>% unlist() %>% unname()
# # cor_fdr_pval_DF_poscorr=cor_fdr_pval_DF_poscorr %>% subset(celltype %in% cor_fdr_pval_DF_poscorr_filter)
# # colnames(cor_fdr_pval_DF_poscorr)=c(colnames(cor_fdr_pval_DF_poscorr)[c(1,2,3)], "FDR", colnames(cor_fdr_pval_DF_poscorr)[c(5,6)])
# # cor_pos_plot=
# #   ggplot(cor_fdr_pval_DF_poscorr, aes(axis1=celltype, axis2=module, y=correlation)) +
# #   geom_alluvium(aes(fill=FDR), width=0.6, curve_type="sigmoid") +
# #   geom_stratum(color="black", linewidth=0.1, alpha=1, fill="white", width=0.6) +
# #   geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3.5) +
# #   scale_x_discrete(limits=c("celltype","module")) +
# #   theme_void() +
# #   theme(legend.position="bottom") +
# #   scale_fill_distiller(palette="Reds")
# # 
# # cor_fdr_pval_DF_negcorr=cor_fdr_pval_DF %>% subset(cor<0) %>% mutate(correlation=abs(cor))
# # # take only those celltypes with a sum of correlation>0.15 with age in all the modules
# # cor_fdr_pval_DF_negcorr_filter=cor_fdr_pval_DF_negcorr %>% 
# #   group_by(celltype) %>% 
# #   summarise_at("correlation", sum) %>% 
# #   subset(correlation>=0.15) %>%
# #   select(celltype) %>% unlist() %>% unname()
# # cor_fdr_pval_DF_negcorr=cor_fdr_pval_DF_negcorr %>% subset(celltype %in% cor_fdr_pval_DF_negcorr_filter)
# # colnames(cor_fdr_pval_DF_negcorr)=c(colnames(cor_fdr_pval_DF_negcorr)[c(1,2,3)], "FDR", colnames(cor_fdr_pval_DF_negcorr)[c(5,6)])
# # cor_neg_plot=
# #   ggplot(cor_fdr_pval_DF_negcorr, aes(axis1=celltype, axis2=module, y=correlation)) +
# #   geom_alluvium(aes(fill=FDR), width=0.6, curve_type="sigmoid") +
# #   geom_stratum(color="black", linewidth=0.1, alpha=1, fill="white", width=0.6) +
# #   geom_text(stat="stratum", aes(label=after_stat(stratum)), size=3.5) +
# #   scale_x_discrete(limits=c("celltype","module")) +
# #   theme_void() +
# #   theme(legend.position="bottom") +
# #   scale_fill_distiller(palette="Blues")
# # 
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WGCNA_correlation_inter_CelltypeModule_alluvial.pdf", height=7.2, width=4)
# # plot(cor_pos_plot)
# # plot(cor_neg_plot)
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # ### Get the hub genes from each module : inter level
# # #####################################
# # ###
# # 
# # # * actually we extract all the genes but arrange them by KME so that the top_n will be the hub genes for each module
# # 
# # ### Load and extract
# # WGCNAobj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Inter.rds")
# # 
# # # get the list of modules:
# # modules=GetModules(WGCNAobj)
# # mods=levels(modules$module)
# # mods=mods[mods!='grey']
# # hb=GetHubGenes(WGCNAobj, n_hubs=length(WGCNAobj@misc$WGCNA$wgcna_genes))
# # 
# # # hubgene network
# # pdf("~/Project_PBMCage/Plots/PBMCage_Plots/HubGene_network_Inter.pdf", height=12, width=12)
# # 
# # allModule_network=
# #   HubGeneNetworkPlot(
# #     WGCNAobj,
# #     n_hubs=5,
# #     n_other=10,
# #     edge_prop=1,
# #     mods="all"
# #   ) +
# #   title("All Modules\n")
# # 
# # CD8T.Tcm_posM4_network=
# #   HubGeneNetworkPlot(
# #     WGCNAobj,
# #     n_hubs=10,
# #     n_other=50,
# #     edge_prop=1,
# #     mods=mods[4]
# #   ) +
# #   title("CD8T.Tcm : M4 - pos.assoc.w/ ages\n")
# # 
# # dev.off()
# # 
# # #####################################
# # 
# # 
# # 
# # 
# # 
