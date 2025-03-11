
### Make slim object for all the assigned cells
#####################################
###
### * Too large dataset for CytoTRACE analysis, so to slim it down

library(dplyr)
library(Seurat)
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")

print("Obj slim starts...")
subset_=THEObj %>%
  NormalizeData() %>%
  FindVariableFeatures()
subset_=SketchData(
  subset_,
  ncells=200000L,
  sketched.assay="sketch",
  method="LeverageScore",
  var.name="leverage.score",
  seed=2024L,
  cast="dgCMatrix",
  verbose=TRUE)
assay.v5=GetAssay(subset_, assay="sketch")
subset_=CreateSeuratObject(assay.v5)
subset_[[]]=THEObj[[]][match(colnames(subset_), colnames(THEObj)),]
saveRDS(subset_, "~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed.rds")
print("Obj slim done!")

#####################################



### CytoTRACE analysis for all the assigned cells
#####################################
###

library(dplyr)
library(Seurat)
library(CytoTRACE)
subset_=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Slimmed.rds")

phe=subset_$Annot.detailed
names(phe)=colnames(subset_)
mat=as.matrix(subset_@assays$RNA@counts)

print("CytoTRACE starts...")
results=CytoTRACE(mat=mat, ncores=15)
plotGene=results$cytoGenes
print("CytoTRACE done!")

plotCytoTRACE(results, phenotype=phe, gene=names(plotGene[1]), outputDir="~/Project_PBMCage/CytoTrace_results/AllAssignedCells_")
cytotrace_table=read.table("~/Project_PBMCage/CytoTrace_results/AllAssignedCells_CytoTRACE_plot_table.txt", sep="\t")
cytotrace_table$agecut=subset_$agecut[match(rownames(cytotrace_table), colnames(subset_))]
cytotrace_table$Annot.inter=subset_$Annot.inter[match(rownames(cytotrace_table), colnames(subset_))]
cytotrace_table$Annot.rough=subset_$Annot.rough[match(rownames(cytotrace_table), colnames(subset_))]
write.table(cytotrace_table, "~/Project_PBMCage/Results/PBMC_Slim_CytoTRACE_table.txt", sep="\t")

#####################################



### Plot the CytoTRACE results that are generated with all the assigned cells
#####################################
###

library(ggplot2)
library(dplyr)
cytotrace_table=read.table("~/Project_PBMCage/Results/PBMC_Slim_CytoTRACE_table.txt")

### Umap For all the agecuts
cytoplot=
  ggplot(cytotrace_table, aes(x=Component1, y=Component2, color=CytoTRACE)) +
  geom_point(size=0.1) +
  scale_color_gradient2(midpoint=0.5, low="blue", mid="white", high="red", space="Lab") +
  labs(title="All ages") +
  theme_classic() +
  theme(title=element_text(size=11),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

namePos=cytotrace_table %>% group_by(Phenotype) %>% summarise(posMedia1=median(Component1), posMedia2=median(Component2)) %>%
  as.data.frame()
namePos$Celltype_inter=cytotrace_table$Annot.inter[match(namePos$Phenotype, cytotrace_table$Phenotype)]

phenoplot=
  ggplot(cytotrace_table, aes(x=Component1, y=Component2, color=Annot.inter)) +
  geom_point(size=0.1) +
  guides(color=guide_legend(override.aes=list(size=3), title="Celltype")) +
  labs(title="All ages") +
  theme_classic() +
  theme(title=element_text(size=11),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        axis.title.x=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  ggrepel::geom_text_repel(data=namePos, aes(x=posMedia1, y=posMedia2, label=Celltype_inter), color="black", size=3, show.legend=FALSE, max.overlaps=3)

allagecut=cowplot::plot_grid(cytoplot, phenoplot, align="hv", ncol=2, rel_widths=c(3,3))

png("~/Project_PBMCage/Plots/PBMCage_Plots/Cytotrace_umap_Overall.png", width=2000, height=800, res=120)
plot(allagecut)
dev.off()

### Umap For each agecuts
agecuts=unique(cytotrace_table$agecut)
agecuts=agecuts[order(as.numeric(sapply(strsplit(agecuts, "~"), function(x) x[1])))]
cytoplot=list()

for (agecut_i in 1:length(agecuts)) {
  cytotrace_table_=subset(cytotrace_table, agecut==agecuts[agecut_i])
  
  cytoplot[[agecut_i]]=
    ggplot(cytotrace_table_, aes(x=Component1, y=Component2, color=CytoTRACE)) +
    geom_point(size=0.1) +
    scale_color_gradient2(midpoint=0.5, low="blue", mid="white", high="red", space="Lab") +
    labs(title=agecuts[agecut_i]) +
    theme_classic() +
    theme(title=element_text(size=11),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9))
}
agecut_plots=cowplot::plot_grid(plotlist=cytoplot, align="hv", ncol=4)

png("~/Project_PBMCage/Plots/PBMCage_Plots/Cytotrace_umap_PerAgecuts.png", width=2000, height=1200, res=120)
plot(agecut_plots)
dev.off()

### Boxplot For each agecuts
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Cytotrace_boxplot_PerAgecuts.pdf", height=5, width=50)
ggplot(cytotrace_table, aes(x=Annot.inter, y=CytoTRACE, color=agecut)) +
  geom_boxplot(coef=5)
dev.off()

#####################################



### Plot the CytoTRACE results to show stemness of different celltypes along with age
#####################################
###

### Load the data
cytotrace_table=read.table("~/Project_PBMCage/Results/PBMC_Slim_CytoTRACE_table.txt", sep="\t")

### Clean the data
celltype_rough_original=names(table(cytotrace_table$Annot.rough))
# remove Other cells as there are not enough observations for a solid conclusion
celltype_rough=celltype_rough[celltype_rough!="Other cells"]
color_assigned=paletteer::paletteer_d("ggsci::category20_d3")[match(celltype_rough, celltype_rough_original)]
# remove the subpopulations without enough observations
selected_type=cytotrace_table %>% subset(Annot.rough %in% celltype_rough) %>% 
  group_by(Phenotype, agecut) %>% count() %>% subset(n>100) %>% 
  ungroup() %>% group_by(Phenotype) %>% count() %>% subset(n>5) %>%
  select(Phenotype) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
df=cytotrace_table %>% subset(Phenotype %in% selected_type) %>%
  mutate(Phenotype=gsub("\\."," ",Phenotype))

### Plot
plot_=
  ggplot(df, aes(x=agecut, y=CytoTRACE, color=Annot.rough)) +
  facet_wrap(~Phenotype, nrow=8, labeller=labeller(Phenotype=scales::label_wrap(10))) +
  # ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.01, color=paletteer::paletteer_d("ggsci::category20_d3")[i]) +
  scale_color_manual(values=color_assigned) +
  theme_classic() +
  labs(x=NULL, y="CytoTRACE", title=NULL) +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  scale_y_continuous(limits=c(0,1)) +
  geom_smooth(aes(group=1), method="gam", 
              show.legend=FALSE, linewidth=0.5, fill="gray93") +
  ggpubr::stat_cor(aes(group=1), size=3)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/CellTRACE_corplots.pdf", height=12, width=6.5)
plot(plot_)
dev.off()

#####################################




### Determine the root_celltype based on the CytoTRACE results that are generated with all the assigned cells
#####################################
###

###
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
cytotrace_table=read.table("~/Project_PBMCage/Results/PBMC_Slim_CytoTRACE_table.txt", sep="\t")
cytotrace_table$age=THEObj$age[match(rownames(cytotrace_table),colnames(THEObj))]

Root_celltype=list()
### root_celltype per intermediate celltypes
cytotrace_table_summarizeInter=cytotrace_table %>%
  group_by(age, Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  group_by(age, Annot.inter) %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeInter=cytotrace_table_summarizeInter[,c("age","Phenotype","Annot.inter","Annot.rough")]
colnames(cytotrace_table_summarizeInter)=c("age","root","from_which_inter","from_which_rough")
Root_celltype[[1]]=cytotrace_table_summarizeInter
# all age together
cytotrace_table_summarizeInter=cytotrace_table %>%
  group_by(Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  group_by(Annot.inter) %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeInter$age="all"
cytotrace_table_summarizeInter=cytotrace_table_summarizeInter[,c("age","Phenotype","Annot.inter","Annot.rough")]
colnames(cytotrace_table_summarizeInter)=c("age","root","from_which_inter","from_which_rough")
Root_celltype[[2]]=cytotrace_table_summarizeInter

### root_celltype per rough celltypes
cytotrace_table_summarizeRough=cytotrace_table %>%
  group_by(age, Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  group_by(age, Annot.rough) %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeRough=cytotrace_table_summarizeRough[,c("age","Annot.inter","Annot.rough")]
colnames(cytotrace_table_summarizeRough)=c("age","root","from_which_rough")
Root_celltype[[3]]=cytotrace_table_summarizeRough
# all age together
cytotrace_table_summarizeRough=cytotrace_table %>%
  group_by(Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  group_by(Annot.rough) %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeRough$age="all"
cytotrace_table_summarizeRough=cytotrace_table_summarizeRough[,c("age","Annot.inter","Annot.rough")]
colnames(cytotrace_table_summarizeRough)=c("age","root","from_which_rough")
Root_celltype[[4]]=cytotrace_table_summarizeRough

### root_celltype in the whole PBMC
cytotrace_table_summarizeAll=cytotrace_table %>%
  group_by(age, Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  group_by(age) %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeAll=cytotrace_table_summarizeAll[,c("age","Annot.rough")]
colnames(cytotrace_table_summarizeAll)=c("age","root")
Root_celltype[[5]]=cytotrace_table_summarizeAll
# all age together
cytotrace_table_summarizeAll=cytotrace_table %>%
  group_by(Phenotype, Annot.inter, Annot.rough) %>%
  summarize(mean=mean(CytoTRACE)) %>%
  ungroup() %>%
  slice_max(order_by=mean, n=1) %>%
  as.data.frame()
cytotrace_table_summarizeAll$age="all"
cytotrace_table_summarizeAll=cytotrace_table_summarizeAll[,c("age","Annot.rough","mean")]
colnames(cytotrace_table_summarizeAll)=c("age","root","CytoTrace")
Root_celltype[[6]]=cytotrace_table_summarizeAll

names(Root_celltype)=c("detailed_root_perAge","detailed_root_allAges",
                       "inter_root_perAge","inter_root_allAges",
                       "rough_root_perAge","rough_root_allAges")
saveRDS(Root_celltype, "~/Project_PBMCage/Tempt_RDS/Pseudotime_RootCelltype.rds")

#####################################



### Pseudotime analysis
#####################################
###

library(dplyr)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(Seurat)
library(scRNAtoolVis)

# load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T_subset=subset(THEObj, Annot.rough=="CD8T cells")

# set sim
sim=SingleCellExperiment(assays=List(counts=CD8T_subset@assays$RNA@counts))

# run umap
CD8T_subset=CD8T_subset %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
saveRDS(list(CD8T_subset, sim), "~/Project_PBMCage/Tempt_RDS/CD8T_subset_sim_metadata.rds")

# transfer metadata
umap=CD8T_subset@reductions$umap@cell.embeddings
colnames(umap)=c("UMAP_1", "UMAP_2")
reducedDims(sim)=SimpleList(UMAP=umap)
for (colname_ in colnames(CD8T_subset[[]])) {
  colData(sim)[[colname_]]=CD8T_subset[[colname_]][,1]
}

# sim analysis
sim=slingshot(sim,
              clusterLabels="Annot.inter",
              start.clus="CD8T.naive",
              reducedDim="UMAP",
              approx_points=FALSE)

# extract pseudotime
pseudotimes=slingPseudotime(sim, na=FALSE)
pseudotimes=as.data.frame(pseudotimes)
pseudotimes$cellid=rownames(pseudotimes)
metadata=CD8T_subset[[]]
metadata$cellid=rownames(metadata)
metadata_plus_pseudotime=merge(metadata, pseudotimes, by="cellid")

# save
saveRDS(list(CD8T_subset, sim, metadata_plus_pseudotime), "~/Project_PBMCage/Tempt_RDS/CD8T_subset_sim_metadata.rds")

#####################################



### Plot pseudotime results
#####################################
###

library(dplyr)
library(tidyr)

### Load
CD8T_subset_sim_metadata.plus.pseudotime=readRDS("~/Project_PBMCage/Tempt_RDS/CD8T_subset_sim_metadata.rds")
CD8T_subset=CD8T_subset_sim_metadata.plus.pseudotime[[1]]
sim=CD8T_subset_sim_metadata.plus.pseudotime[[2]]
metadata_plus_pseudotime=CD8T_subset_sim_metadata.plus.pseudotime[[3]]

### Make the dataframe
DF_LONG_List=list()
celltype_d=unique(metadata_plus_pseudotime$Annot.detailed)
for (i in 1:length(celltype_d)) {
  subset_=subset(metadata_plus_pseudotime, Annot.detailed==celltype_d[i])
  split_=split(subset_$Lineage1, subset_$age)
  split_df_=plyr::ldply(split_, rbind)
  rownames(split_df_)=split_df_$.id
  split_df_=select(split_df_, -.id)
  split_df_=t(split_df_)
  split_df_=as.data.frame(split_df_)
  df_long=pivot_longer(split_df_, cols=colnames(split_df_), names_to="age", values_to="pseudotime")
  df_long$Annot.detailed=celltype_d[i]
  DF_LONG_List[[i]]=df_long
}
DF_LONG=data.table::rbindlist(DF_LONG_List)
DF_LONG=na.omit(DF_LONG)

### Plot the CD8 T cell pseudotime
library(ggplot2)
library(ggsci)
plot_=
  ggplot(metadata_plus_pseudotime, aes(x=age, y=Lineage1, color=Annot.detailed)) +
  theme_classic() +
  labs(x="Age", y="Pseudotime") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        title=element_text(size=10)) +
  scale_x_discrete(breaks=seq(19, 97, 5)) +
  guides(color=guide_legend(title="Celltype")) +
  scale_color_d3("category20") +
  geom_smooth(aes(group=Annot.detailed), show.legend=TRUE, linewidth=0.5, fill="gray93")

### Plot the CD8 T cell freq (for comparison with pseudotime alterations)
pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T_cell_subset=subset(pbmc.seu_merged, Annot.rough=="CD8T cells")
ID_=CD8T_cell_subset$donor_id
Age_=CD8T_cell_subset$age
CellType_=CD8T_cell_subset$Annot.detailed
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","Annot.detailed","Freq")
ID_Age_match=data.frame(ID=ID_, age=Age_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$Annot.detailed=as.factor(data$Annot.detailed)
plot_2=
  ggplot(data, aes(x=age, y=Freq, color=Annot.detailed)) +
  theme_classic() +
  labs(x="Age", y="Cell percent (%)") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        title=element_text(size=10)) +
  scale_x_discrete(breaks=seq(19, 97, 5)) +
  guides(color=guide_legend(title="Celltype")) +
  scale_color_d3("category20") +
  geom_smooth(aes(group=Annot.detailed), show.legend=TRUE, linewidth=0.5, fill="gray93")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Pseudotime_CD8T.pdf", width=6, height=6)
plot(plot_)
plot(plot_2)
dev.off()

#####################################



### Pseudotime Umap plotting for each agecut'
#####################################
###
# Merge some agecuts according to the Pseudotime_CD8T.pdf

###

library(dplyr)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(Seurat)
library(scRNAtoolVis)

### Add the agecut metadata
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
CD8T=subset(THEObj, Annot.rough=="CD8T cells")

AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
AGECUT=list(c(19:24),c(25:30),c(31:51),c(52:58),c(59:73),c(74:81),c(82:91),c(92:97)) # merge

AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}
agecut_df=CD8T[[]]
agecut_df$age=as.factor(agecut_df$age)
agecut_df$agecut=agecut_df$age
for (i in 1:length(AGECUT)) {
  agecut_df=agecut_df %>%
    mutate(agecut=ifelse(age %in% AGECUT[[i]], AGECUT_name[i], agecut))
}
table(rownames(agecut_df)==colnames(CD8T)) # check
CD8T=AddMetaData(CD8T, metadata=agecut_df$agecut, col.name="agecut_merged")

### Analyze and plot the trajectory in each agecut_merged
agecut_merged_list=names(table(CD8T$agecut_merged))
SIMs=list()
for (i in 1:length(agecut_merged_list)) {
  subset_=subset(CD8T, agecut_merged==agecut_merged_list[i])
  
  # make slim obj to 30,000 cells if the orig. subset_ exceeds 30,000 to speed up the slingshot process
  if (ncol(subset_)>30000) {
    print("Orig. subset_ exceeds 30,000 cells. Start slimming...")
    subset_=subset_ %>%
      NormalizeData() %>%
      FindVariableFeatures()
    subset_=SketchData(
      subset_,
      ncells=30000L,
      sketched.assay="sketch",
      method="LeverageScore",
      var.name="leverage.score",
      seed=2024L,
      cast="dgCMatrix",
      verbose=TRUE)
    assay.v5=GetAssay(subset_, assay="sketch")
    subset_=CreateSeuratObject(assay.v5)
    subset_[[]]=CD8T[[]][match(colnames(subset_), colnames(CD8T)),]
    print("Obj slimming done!")
  }
  
  # # * if Error in shortest_paths, this's because some celltypes have only 1 single cell, then run: 
  # print(table(subset_$Annot.detailed))
  # celltype_containing_more_cells=names(table(subset_$Annot.detailed))[table(subset_$Annot.detailed)!=1]
  # subset_=subset(subset_, Annot.detailed %in% celltype_containing_more_cells)
  # print(table(subset_$Annot.detailed))
  
  sim0=SingleCellExperiment(assays=List(counts=subset_@assays$RNA@counts))
  subset_=subset_ %>%
    SCTransform() %>%
    RunPCA() %>%
    FindNeighbors(dims=1:50) %>%
    RunUMAP(dims=1:50)
  umap=subset_@reductions$umap@cell.embeddings
  colnames(umap)=c("UMAP_1", "UMAP_2")
  reducedDims(sim0)=SimpleList(UMAP=umap)
  for (colname_ in colnames(subset_[[]])) {
    colData(sim0)[[colname_]]=subset_[[colname_]][,1]
  }
  
  # trajectory analysis
  sim0=slingshot(sim0,
                 clusterLabels="Annot.detailed",
                 start.clus="CD8T.naive.LINC02446",
                 reducedDim="UMAP",
                 approx_points=FALSE
  )
  SIMs[[i]]=sim0
  pseudotimes=slingPseudotime(sim0, na=FALSE)
  # 
  # # set the color and the clusters
  # clustering=subset_$Annot.detailed
  # clustering=as.factor(clustering)
  # levels(clustering)=sort(levels(clustering))
  # pal=paletteer::paletteer_d("ggsci::category20_d3")
  # 
  # # plot
  # pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Trajectory_",i,".pdf"), height=5, width=8)
  # layout(matrix(c(1,1,2), nrow=1, ncol=3, byrow=TRUE))
  # par(mar=c(5.1, 4.1, 4.1, 2.1))
  # plot(
  #   reducedDims(sim0)[[1]],
  #   col=pal[match(clustering, levels(clustering))],
  #   pch=19,
  #   cex=0.5,
  #   main=paste0(agecut_merged_list[i])
  # )
  # lines(
  #   SlingshotDataSet(sim0),
  #   lwd=1,
  #   type='lineages',
  #   col='black',
  #   show.constraints=TRUE
  # )
  # par(mar=c(5.1, 0, 4.1, 0))
  # plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1, main="\n")
  # legend(
  #   "topleft",
  #   legend=levels(clustering),
  #   fill=pal[1:length(levels(clustering))]
  # )
  # dev.off()
}
# ### Merge the plots of all the agecut_merged
# qpdf::pdf_combine(input=paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Trajectory_",c(1:8),".pdf"),
#                   output="~/Project_PBMCage/Plots/PBMCage_Plots/Trajectory_CD8T_For.Each.Agecutmerged.pdf")

### Save the SIMs obj
saveRDS(SIMs, "~/Project_PBMCage/Tempt_RDS/Pseudotime_CD8T.rds")

### Plot by ggplot2
SIMs_CD8T=readRDS("~/Project_PBMCage/Tempt_RDS/Pseudotime_CD8T.rds")

sim_=SIMs_CD8T[[2]]
umapData=reducedDims(sim_)[[1]] %>% as.data.frame()
mst=slingMST(sim_, as.df=TRUE) %>% mutate(Cluster=gsub("CD8T\\.","",Cluster))

ggplot(umapData, aes(x=UMAP_1, y=UMAP_2)) +
  ggrastr::geom_point_rast(size=0.1, alpha=0.1) +
  geom_point(data=mst, size=1, alpha=0.5) +
  geom_path(data=mst %>% arrange(Order), aes(group=Lineage), linewidth=0.5) +
  ggrepel::geom_label_repel(data=mst %>% arrange(Order) %>% subset(!duplicated(.$Cluster)), 
                            aes(x=UMAP_1, y=UMAP_2, label=Cluster)) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        panel.grid=element_blank(),
        title=element_text(size=11))
  
SlingshotDataSet(SIMs_CD8T[[1]])

AGECUT=list(c(19:24),c(25:30),c(31:51),c(52:58),c(59:73),c(74:81),c(82:91),c(92:97)) # merge
AGECUT_name=c()
for (i in 1:length(AGECUT)) {
  AGECUT_name=c(AGECUT_name, paste0(AGECUT[[i]][1],"~",AGECUT[[i]][length(AGECUT[[i]])]))
}


mst=lapply(SIMs_CD8T, function(obj) slingMST(obj, as.df=TRUE))
mst=lapply(1:length(mst), function(i) mst[[i]] %>% mutate(agecut=AGECUT_name[i]))
mst_total=data.table::rbindlist(mst)

eee=list()

aaa=mst[[1]]
aaa=aaa %>% select(-c(UMAP_1,UMAP_2)) %>% tidyr::pivot_wider(names_from="Lineage", values_from="Order")
# for (i in 1:(ncol(aaa)-1)) {
#   aaa[[as.character(i)]]=sapply(c(1:14), function(j) ifelse(is.na(aaa[[as.character(i)]][j]), NA, j))
# }
colnames(aaa)=c("Cluster", paste0("Lin",colnames(aaa)[!grepl("Cluster",colnames(aaa))]))
celltypes=aaa$Cluster

aaa=aaa %>% 
  tidyr::pivot_longer(cols=colnames(.)[!grepl("Cluster",colnames(.))], names_to="Lineage", values_to="Order")
bbb=aaa %>% na.omit(.) %>% group_by(Cluster, Lineage) %>% dplyr::count() %>% 
  mutate(lin_to_mark=as.numeric(gsub("Lin","",Lineage))) %>%
  ungroup() %>% group_by(Cluster) %>%
  slice_min(lin_to_mark, n=1) %>%
  ungroup() %>% mutate(mark=Cluster) %>% select(-c(n, lin_to_mark))

ddd=left_join(aaa, bbb, by=c("Cluster","Lineage"))
eee[[1]]=ddd %>% mutate(AgeCut="one")

ggplot(ddd %>% arrange(Order), aes(x=Lineage, y=Order, alpha=Lineage)) +
  scale_alpha_manual(values=rep(1/length(unique(aaa$Lineage)), length(unique(aaa$Lineage)))) +
  geom_label(aes(label=mark)) +
  geom_point() +
  geom_path()




aaa=mst[[2]]
aaa=aaa %>% select(-c(UMAP_1,UMAP_2)) %>% tidyr::pivot_wider(names_from="Lineage", values_from="Order")
# for (i in 1:(ncol(aaa)-1)) {
#   aaa[[as.character(i)]]=sapply(c(1:14), function(j) ifelse(is.na(aaa[[as.character(i)]][j]), NA, j))
# }
colnames(aaa)=c("Cluster", paste0("Lin",colnames(aaa)[!grepl("Cluster",colnames(aaa))]))
celltypes=aaa$Cluster

aaa=aaa %>% 
  tidyr::pivot_longer(cols=colnames(.)[!grepl("Cluster",colnames(.))], names_to="Lineage", values_to="Order")
bbb=aaa %>% na.omit(.) %>% group_by(Cluster, Lineage) %>% dplyr::count() %>% 
  mutate(lin_to_mark=as.numeric(gsub("Lin","",Lineage))) %>%
  ungroup() %>% group_by(Cluster) %>%
  slice_min(lin_to_mark, n=1) %>%
  ungroup() %>% mutate(mark=Cluster) %>% select(-c(n, lin_to_mark))

ddd=left_join(aaa, bbb, by=c("Cluster","Lineage"))
eee[[2]]=ddd %>% mutate(AgeCut="two")

ggplot(ddd %>% arrange(Order), aes(x=Lineage, y=Order, alpha=Lineage)) +
  scale_alpha_manual(values=rep(1/length(unique(aaa$Lineage)), length(unique(aaa$Lineage)))) +
  geom_label(aes(label=mark)) +
  geom_point() +
  geom_path()

eee_2=data.table::rbindlist(eee)

ggplot(eee_2 %>% arrange(Order), aes(x=Lineage, y=Order, alpha=Lineage)) +
  facet_wrap(~AgeCut) +
  scale_alpha_manual(values=rep(1/length(unique(aaa$Lineage)), length(unique(aaa$Lineage)))) +
  geom_label(aes(label=mark)) +
  geom_point() +
  geom_path()


#####################################



