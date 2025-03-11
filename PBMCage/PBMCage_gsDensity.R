### ====== / Installation / ======
# Turn off warning-error-conversion, because the tiniest warning stops installation
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"="true")
Sys.unsetenv("GITHUB_PAT") # remove my PAT otherwise github packages cannot be installed
remotes::install_github("https://github.com/RausellLab/CelliD", force=TRUE) # get from github which has been updated for Assay5
remotes::install_github("https://github.com/qingnanl/gsdensity", force=TRUE)





### ====== / GSDensity to score the pathway activities / ======
library(gsdensity)
library(ggplot2)
library(Seurat)
library(dplyr)

### Load the GO_BP genesets
mdb_c5=msigdbr::msigdbr(species="Homo sapiens", category="C5")
mdb_c5_bp=mdb_c5[mdb_c5$gs_subcat=="GO:BP",]
# convert to list
GO_BP_genelist=mdb_c5_bp %>% dplyr::select(gs_name, gene_symbol) %>% subset(!duplicated(.)) %>%
  split(.$gs_name) %>%
  lapply(., function(x) x$gene_symbol %>% subset(!duplicated(.)))

### Load the data and slim it down
# * Too large dataset for gsDensity
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
print("Obj slim starts...")
subset_=THEObj %>%
  NormalizeData() %>%
  FindVariableFeatures()
subset_=SketchData(
  subset_,
  ncells=10000L,
  sketched.assay="sketch",
  method="LeverageScore",
  var.name="leverage.score",
  seed=2024L,
  cast="dgCMatrix",
  verbose=TRUE)
assay.v5=GetAssay(subset_, assay="sketch")
subset_=CreateSeuratObject(assay.v5)
subset_[[]]=THEObj[[]][match(colnames(subset_), colnames(THEObj)),]
subset_[["RNA"]]=as(object=subset_[["RNA"]], Class="Assay") # convert to Assay
saveRDS(subset_, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_5wCellsForgsDensity.rds")

### MCA analysis
# compute cell and gene embedding
ce=compute.mca(subset_, assay="RNA", slot="data")
saveRDS(ce, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_ceOBJ.rds")
# find gene sets with differential density
source("~/Rscripts/compute.kld2.R")
res=compute.kld2(coembed=ce, 
                 genes.use=intersect(rownames(ce),rownames(subset_)),
                 n.grids=100,
                 gene.set.list=GO_BP_genelist[31:35],
                 gene.set.cutoff=3,
                 n.times=100)
saveRDS(res, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_DiffPathRes.rds")
# build nearest neighbor graph
subset_=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_5wCellsForgsDensity.rds")
ce=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_ceOBJ.rds")
res=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_DiffPathRes.rds")
el=compute.nn.edges(coembed=ce, nn.use=300)
# get label propagation probability im each cell for all the sig. genesets
sig.genesets=(res %>% subset(p.adj<0.05))$gene.set
CVs=list()
CVs_names=c()
for (i in 1:length(sig.genesets)) {
  cv=run.rwr(el=el, gene_set=GO_BP_genelist[[sig.genesets[i]]], cells=colnames(subset_))
  CVs=c(CVs, list(cv))
  CVs_names=c(CVs_names, sig.genesets[i])
  names(CVs)=CVs_names
  print(paste0(i," in ",length(sig.genesets)," is done."))
  saveRDS(CVs, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVs.rds")
}
# label the cells with positive/negative for all the genesets
CVs=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVs.rds")
CLs=lapply(CVs, function(cv) compute.cell.label(cv))
saveRDS(CLs, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CLs.rds")

### Arrange the MCA analysis results
# remove the genesets with no pos cell
CVs=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVs.rds")
CLs=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CLs.rds")
selected_genesets=names(CLs)[lapply(CLs, function(x) length(table(x))!=0) %>% unlist(.)]
CVs=CVs[selected_genesets]
CLs=CLs[selected_genesets]
# arrange the CVs and CLs as dataframes
table(names(CVs)==names(CLs)) # check
CV.CLs=lapply(1:length(CVs), function(idx) {
  data.frame(cellid=names(CVs[[idx]]), geneset=names(CVs)[idx], propagation.p=CVs[[idx]]) %>%
    left_join(data.frame(cellid=names(CLs[[idx]]), geneset=names(CLs)[idx], label=CLs[[idx]]))
})
names(CV.CLs)=names(CVs)
saveRDS(CV.CLs, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CV.CLs.rds")
CV.CL_df=CV.CLs %>% data.table::rbindlist()
data.table::fwrite(CV.CL_df, "~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVandCL_RESULTs.csv.gz")
# add the Annot info to the CV/CLs dfs
subset_=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_5wCellsForgsDensity.rds")
CV.CL_dfs=readRDS("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CV.CLs.rds")
CV.CL_dfs_annotated=lapply(CV.CL_dfs, function(df) {
  df %>% left_join(subset_[[]] %>% tibble::rownames_to_column("cellid") %>% 
                     dplyr::select(cellid, donor_id, age, agecut, sex, Annot.detailed, Annot.inter, Annot.rough))
})
# summarize the dfs to get the most sig. genesets
CV.CL_dfs_annotated_summarized=lapply(CV.CL_dfs_annotated, function(df) {
  df %>% 
    mutate(mean.prop.all=mean(propagation.p, na.rm=T)) %>%
    group_by(label) %>%
    mutate(count.neg.or.pos=n()) %>%
    ungroup() %>% 
    group_by(mean.prop.all, count.neg.or.pos, donor_id, age, agecut, sex, Annot.rough) %>%
    summarize(mean.prop=mean(propagation.p, na.rm=T))
})
CV.CL_dfs_annotated_summarized_all=lapply(1:length(CV.CL_dfs_annotated_summarized), function(idx) {
  CV.CL_dfs_annotated_summarized[[idx]] %>% 
    mutate(geneset=names(CV.CL_dfs_annotated_summarized)[idx])
}) %>% data.table::rbindlist()
CV.CL_dfs_annotated_sigGeneset=CV.CL_dfs_annotated_summarized_all %>%
  group_by(geneset) %>%
  summarize(mean_=mean(mean.prop, na.rm=T),
            median_=median(mean.prop, na.rm=T),
            min_=min(mean.prop[mean.prop!=0], na.rm=T),
            max_=max(mean.prop, na.rm=T)) %>%
  slice_max(mean_, n=20)
CV.CL_dfs_annotated_sigGeneset.perAnnotrough=CV.CL_dfs_annotated_summarized_all %>%
  group_by(geneset, Annot.rough) %>%
  summarize(mean_=mean(mean.prop, na.rm=T),
            median_=median(mean.prop, na.rm=T),
            min_=min(mean.prop[mean.prop!=0], na.rm=T),
            max_=max(mean.prop, na.rm=T)) %>%
  ungroup() %>%
  group_by(Annot.rough) %>%
  slice_max(mean_, n=20)
CV.CL_dfs_annotated_sigGeneset.perAnnotrough_mean=CV.CL_dfs_annotated_sigGeneset.perAnnotrough %>%
  dplyr::select(geneset, Annot.rough, mean_) %>%
  tidyr::pivot_wider(id_cols="geneset", names_from="Annot.rough", values_from="mean_")




### Plot
# rearrange the df to the wide format
CV.CL_df=data.table::fread("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVandCL_RESULTs.csv.gz")
CV_wide=CV.CL_df %>% 
  dplyr::select(-label) %>%
  tidyr::pivot_wider(id_cols="cellid", names_from="geneset", values_from="propagation.p") %>%
  tibble::column_to_rownames("cellid")
CL_wide=CV.CL_df %>% 
  dplyr::select(-propagation.p) %>%
  tidyr::pivot_wider(id_cols="cellid", names_from="geneset", values_from="label") %>%
  tibble::column_to_rownames("cellid")


### Run Umap for visualization
subset_=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed_5wCellsForgsDensity.rds")
subset_=subset_ %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
# check
DimPlot(subset_, group.by="Annot.rough")

### Add the label propagation probability and labels to the metadata
table(rownames(subset_[[]])==rownames(CV_wide)) # check
subset_[[]]=subset_[[]] %>% cbind(CV_wide)

subset_$GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS__bin=cl[colnames(subset_)]

### Plot the probabilities
CV.CL_res=data.table::fread("~/Project_PBMCage/Tempt_RDS/PBMCage_gsDensity_CVandCL_RESULTs.csv.gz")
FeaturePlot(subset_, features=CV.CL_dfs_annotated_sigGeneset.perAnnotrough$geneset[1], raster=T)


