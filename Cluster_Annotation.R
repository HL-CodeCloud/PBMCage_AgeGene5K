setwd("~/Project_PBMCage")

library(SeuratDisk)
library(SeuratData)
library(Seurat)
library(edgeR)
library(vcd)
library(ggplot2)
library(scRNAtoolVis)
library(pheatmap)
library(dittoSeq)
library(ggsignif)
library(tidyr)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(harmony)
library(ggpubr)
library(patchwork)
library(cowplot)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(CellChat)
library(data.table)
library(paletteer)
library(ComplexHeatmap)
library(SCpubr)
library(RColorBrewer) 

col.p=col.pDark=scater:::.get_palette("tableau20")
col.pMedium=scater:::.get_palette("tableau10medium")
col.pDark=scater:::.get_palette("tableau20")[2*(1:10)-1]
col.pLight=scater:::.get_palette("tableau20")[2*(1:10)]
col.sex=scater:::.get_palette("tableau20")[c(1,11)]
col.age=paletteer::paletteer_c("ggthemes::Gray", 78)

### Check the data from the cellxgene db
#####################################
###

### According to the original paper, the data have been demultiplexed and doublet-removed
PBMC=readRDS("~/Project_PBMCage/cellxgene.rds")
str(PBMC)
# in the data slot is the raw UMI
head(PBMC@assays$RNA@counts)
head(PBMC@assays$RNA@data)
head(PBMC@assays$RNA@scale.data)
# gene names were ENSEMBL ID
head(PBMC@assays$RNA@meta.features$feature_name)
head(rownames(PBMC))
# the obj has been de-multiplexed
table(PBMC$pool_number)
options(max.print=6000); table(PBMC$donor_id, PBMC$pool_number)
# the obj has been doublet-removed and quality-controlled
# VlnPlot(PBMC, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=3) 

#####################################



### Reconstruct obj with Gene SYMBOL and useful metadata
#####################################
###

###
PBMC_counts=LayerData(PBMC, assay="RNA", layer="data", fast=TRUE)
dimnames(PBMC_counts)=list(PBMC@assays$RNA@meta.features$feature_name, colnames(PBMC_counts))

PBMC_meta=PBMC[[]]
PBMC_meta_chosen=PBMC_meta[, c(2,3,4,5,6,7,9,20,24,27)]
PBMC_meta_chosen$sex=factor(PBMC_meta_chosen$sex, levels=c("male","female"))
PBMC_meta_chosen$age=as.factor(PBMC_meta_chosen$age)

object=CreateSeuratObject(counts=PBMC_counts, meta.data=PBMC_meta_chosen)

# add percent.rb
object[["percent.rb"]]=PercentageFeatureSet(object, pattern="^RP[SL]")

saveRDS(object, "~/Project_PBMCage/raw_data/Seurats_all.rds")

# save seurt4 obj for scPred
query_counts=LayerData(object, assay="RNA", layer="counts", fast=TRUE)
dimnames(query_counts)=list(rownames(object@assays$RNA@features), rownames(object@assays$RNA@cells))
query_meta=object[[]]
save(query_counts, query_meta, file="~/Project_PBMCage/raw_data/Seurats_ElementsForScpred.RData")

#####################################



### Percent of Mt, Rb, and cellcyle evaluation
#####################################
###

### Add cellcyle info
# match ENSEMBL in the data to SYMBOL in the cc.genes
anno=rbind(bitr(cc.genes.updated.2019$s.genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db"),
           bitr(cc.genes.updated.2019$g2m.genes, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db"))
anno[duplicated(anno$SYMBOL),]
for (symbol_ in anno[duplicated(anno$SYMBOL), 1]) {
  right_one=anno[anno$SYMBOL==symbol_,] %>% subset(!(ENSEMBL %in% rownames(PBMC)))
  anno=anno %>% subset(!(ENSEMBL %in% right_one$ENSEMBL))
}
s.features=anno[match(cc.genes.updated.2019$s.genes, anno$SYMBOL), 2]
g2m.features=anno[match(cc.genes.updated.2019$g2m.genes, anno$SYMBOL), 2]
# normalize, find HVGs, and scale
PBMC=SetAssayData(PBMC, slot="counts", new.data=GetAssayData(PBMC, slot="data"))
PBMC=NormalizeData(PBMC)
# cellcycle scoring
PBMC=CellCycleScoring(PBMC, 
                      s.features=s.features, 
                      g2m.features=g2m.features,
                      set.ident=TRUE)
# create a dataframe for cellcyle info
list_=list()
pb=txtProgressBar(min=0, max=length(names(table(PBMC$donor_id))), style=3, width=60, char="=")
for (idx in 1:length(names(table(PBMC$donor_id)))) {
  ID=names(table(PBMC$donor_id))[idx]
  temp=subset(PBMC@"meta.data", donor_id==ID)
  Age=as.integer(temp$age[1])
  Sex=as.character(temp$sex[1])
  percent.G1=as.numeric(table(temp$Phase)/nrow(temp))[grepl("G1",names(table(temp$Phase)))]
  percent.G2M=as.numeric(table(temp$Phase)/nrow(temp))[grepl("G2M",names(table(temp$Phase)))]
  percent.S=as.numeric(table(temp$Phase)/nrow(temp))[grepl("S",names(table(temp$Phase)))]
  PerRow=data.table(ID, Age, Sex, percent.G1, percent.G2M, percent.S)
  list_=c(list_, list(PerRow))
  setTxtProgressBar(pb, idx)
}
close(pb)
df=rbindlist(list_)

### Create a dataframe for Mt and Rb info
Sample_files=gsub(".rds","",list.files("~/Project_PBMCage/raw_data/preprocessedData/"))
list2_=list()
pb=txtProgressBar(min=0, max=length(Sample_files), style=3, width=60, char="=")
for (i in 1:length(Sample_files)){
  temp=readRDS(file.path("~/Project_PBMCage/raw_data/preprocessedData", paste0(Sample_files[i],".rds")))
  ids=names(table(temp$donor_id))[table(temp$donor_id)!=0]
  for (idx in 1:length(ids)) {
    ID=ids[idx]
    temp2=subset(temp@"meta.data", donor_id==ID)
    percent.mt=mean(temp2$percent.mt)
    percent.rb=mean(temp2$percent.rb)
    PerRow=data.table(ID, percent.mt, percent.rb)
    list2_=c(list2_, list(PerRow))
  }
  setTxtProgressBar(pb, i)
}
close(pb)
df2=rbindlist(list2_)

### Merge the two info dataframe and save it
merge_=merge(df, df2, by="ID")
merge_$Age=as.factor(merge_$Age)
merge_$Sex=as.factor(merge_$Sex); merge_$Sex=factor(merge_$Sex, levels=c("male","female"))
write.table(merge_, "~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt", sep="\t")

### Plot
p1=
  ggplot(merge_, aes(x=Age, y=percent.mt, color=Sex)) +
  geom_point(size=1.5, shape="o") +
  scale_color_manual(values=col.sex) +
  theme_classic() +
  geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
  stat_cor(show.legend=FALSE, size=3.5) +
  labs(title="Percent of reads mapped to mitochondrial genes") +
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,1,0,1), "lines"),
    plot.title=element_text(size=12))

p2=
  ggplot(merge_, aes(x=Age, y=percent.rb, color=Sex)) +
  geom_point(size=1.5, shape="o") +
  scale_color_manual(values=col.sex) +
  theme_classic() +
  geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
  stat_cor(show.legend=FALSE, size=3.5) +
  labs(title="Percent of reads mapped to ribosomal genes") +
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,1,0,1), "lines"),
    plot.title=element_text(size=12))

p3=
  ggplot(merge_, aes(x=Age, y=percent.G1, color=Sex)) +
    geom_point(size=1.5, shape="o") +
    scale_color_manual(values=col.sex) +
    theme_classic() +
    geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
    stat_cor(show.legend=FALSE, size=3.5) +
    labs(title="Percent of G1 phase") +
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,1,0,1), "lines"),
    plot.title=element_text(size=12))

p4=  
  ggplot(merge_, aes(x=Age, y=percent.G2M, color=Sex)) +
    geom_point(size=1.5, shape="o") +
    scale_color_manual(values=col.sex) +
    theme_classic() +
    geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
    stat_cor(show.legend=FALSE, size=3.5)  +
    labs(title="Percent of G2M phase") +
  theme(
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,1,0,1), "lines"),
    plot.title=element_text(size=12))

p5=
  ggplot(merge_, aes(x=Age, y=percent.S, color=Sex)) +
    geom_point(size=1.5, shape="o") +
    scale_color_manual(values=col.sex) +
    theme_classic() +
    geom_smooth(method='lm', show.legend=FALSE, linewidth=0.5, fill="gray") +
    stat_cor(show.legend=FALSE, size=3.5) +
    labs(title="Percent of S phase") +
  theme(
    plot.margin=unit(c(1,1,1,1), "lines"),
    plot.title=element_text(size=12))

plots=p1/p2/p3/p4/p5

pdf("~/Project_PBMCage/Plots/MtRbCellcycle_Eval.pdf", width=5, height=15)
plot(plots)
dev.off()

#####################################



### Sample clustering of donor_id
#####################################
###

### For each age, rowsum on donor_id
Sample_files=gsub(".rds","",list.files("~/Project_PBMCage/raw_data/preprocessedData/"))
pb=txtProgressBar(min=0, max=length(Sample_files), style=3, width=60, char="=")
counts_sum_list=list()
for (i in 1:length(Sample_files)) {
  temp=readRDS(file.path("~/Project_PBMCage/raw_data/preprocessedData", paste0(Sample_files[i],".rds")))
  donor_id=temp$donor_id
  counts_matrix=as.matrix(temp@assays$RNA@counts)
  counts_sum=t(rowsum(t(counts_matrix), group=donor_id))
  counts_sum_list=c(counts_sum_list, list(counts_sum))
  setTxtProgressBar(pb, i)
}
close(pb)
names(counts_sum_list)=Sample_files

### cbind all the donors
df_by_donor=merge(counts_sum_list[[1]], counts_sum_list[-1], by="row.names")
rownames(df_by_donor)=df_by_donor[,1]
df_by_donor=df_by_donor[,2:ncol(df_by_donor)]
colnames(df_by_donor)=gsub("^age_[0-9]{2}\\.|X", "", colnames(df_by_donor))

### Create a pseudo-bulk DGEList obj
pseudoDGE=DGEList(df_by_donor)
metainfo=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt", sep="\t", stringsAsFactors=FALSE)
pseudoDGE$samples$Donor_ID=colnames(df_by_donor)
pseudoDGE$samples$Donor_ID=as.factor(pseudoDGE$samples$Donor_ID)
pseudoDGE$samples$Age=metainfo$Age[match(colnames(df_by_donor), metainfo$ID)]
pseudoDGE$samples$Age=as.factor(pseudoDGE$samples$Age)
pseudoDGE$samples$Sex=metainfo$Sex[match(colnames(df_by_donor), metainfo$ID)]
pseudoDGE$samples$Sex=as.factor(pseudoDGE$samples$Sex); pseudoDGE$samples$Sex=factor(pseudoDGE$samples$Sex, levels=c("male","female"))
saveRDS(pseudoDGE, "~/Project_PBMCage/pseudoDGE_by_donor.rds")

### Calculate cpm
CPM_Matrix=edgeR::cpm(pseudoDGE, log=TRUE, prior.count=1)

### Plot sample distance
sampleDists=dist(t(CPM_Matrix))
sampleDistMatrix=as.matrix(sampleDists)
colnames(sampleDistMatrix)=NULL
library("RColorBrewer") 
col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255) 
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=col)

### Plot the heatmap
# find the sig. genes across ages
PBMC=readRDS("~/Project_PBMCage/cellxgene.rds")
PBMC=SetAssayData(PBMC, slot="counts", new.data=GetAssayData(PBMC, slot="data"))
PBMC=PBMC %>%
  NormalizeData()
age_cut=data.frame(Age=PBMC$age) %>%
  mutate(age_cut=ifelse(Age %in% c(19:65), "19:65", 
                        ifelse(Age %in% c(66:100), "66:100", " ")))
PBMC$age_cut=age_cut$age_cut
Idents(PBMC)="age_cut"
DEGs=FindAllMarkers(PBMC, min.pct=0, logfc.threshold=0.1, only.pos=TRUE)
DEGs$SYMBOL=PBMC@assays$RNA@meta.features$feature_name[match(DEGs$gene, rownames(PBMC))]
saveRDS(DEGs, "~/Project_PBMCage/sigGenes_young_elder.rds")
# prepare the matrix
df=CPM_Matrix[match(DEGs$SYMBOL, rownames(CPM_Matrix)), ]
df_scaled=t(scale(t(df)))
# column annotation
sudocolor=rep("white", ncol(df_scaled)); names(sudocolor)=1:ncol(df_scaled)
group=HeatmapAnnotation(class=
                          anno_block(gp=gpar(fill=col.age, col="white")),
                        Age=1:ncol(df_scaled),
                        annotation_height=unit(c(2, 0.001), c("mm", "mm")),
                        col=list(Age=sudocolor),
                        annotation_label=gt_render(c("Age", "<sup>Age</sup>"), gp=gpar(fontsize=14)),
                        show_legend=c("Age"=FALSE))
quantile(df_scaled, c(0.01,0.05,0.95,0.99)) # for setting matrix colors
cols=circlize::colorRamp2(c(-3,0,3), c("#260F99","white","#990F0F"))
# plot
p=
Heatmap(df_scaled,
        col=cols,
        column_split=length(Sample_files),
        top_annotation=group,
        column_gap=unit(0.1, "mm"),
        
        cluster_columns=TRUE,
        show_column_dend=TRUE,
        column_title=NULL,
        show_column_names=FALSE,
        
        cluster_rows=TRUE,
        show_row_dend=FALSE,
        row_title=NULL,
        show_row_names=FALSE,
        
        name=" ")
pdf("~/Project_PBMCage/Plots/Heatmap_by_age.pdf", width=15, height=6)
plot(p)
dev.off()

#####################################



### Sample clustering of age
#####################################
###

### Reduce the per_donor_matrix to per_age_matrix
temp=readRDS("~/Project_PBMCage/pseudoDGE_by_donor.rds")
df_by_age=t(rowsum(t(temp$counts), group=temp$samples$Age))

### Create a pseudo-bulk DGEList obj
pseudoDGE=DGEList(df_by_age)
metainfo=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt", sep="\t", stringsAsFactors=FALSE)
pseudoDGE$samples$Age=colnames(df_by_age)
pseudoDGE$samples$Age=as.factor(pseudoDGE$samples$Age)
saveRDS(pseudoDGE, "~/Project_PBMCage/pseudoDGE_by_age.rds")

### Calculate cpm
CPM_Matrix=edgeR::cpm(pseudoDGE, log=TRUE, prior.count=1)

### Plot sample distance
sampleDists=dist(t(CPM_Matrix))
sampleDistMatrix=as.matrix(sampleDists)
sampleDistMatrix[lower.tri(sampleDistMatrix)]=NA
colnames(sampleDistMatrix)=NULL
col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255) 
ht=
Heatmap(sampleDistMatrix,
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        col=col,
        na_col="white",
        row_names_gp=gpar(fontsize=6),
        heatmap_legend_param=list(title=NULL)
        )
pdf("~/Project_PBMCage/Plots/SampleClustering_by_age.pdf")
draw(ht, background="transparent")
dev.off()

### Plot the heatmap
# find the sig. genes across ages
PBMC=readRDS("~/Project_PBMCage/cellxgene.rds")
PBMC=SetAssayData(PBMC, slot="counts", new.data=GetAssayData(PBMC, slot="data"))
PBMC=PBMC %>%
  NormalizeData()
age_cut=data.frame(Age=PBMC$age) %>%
  mutate(age_cut=ifelse(Age %in% c(19:65), "19:65", 
                        ifelse(Age %in% c(66:100), "66:100", " ")))
PBMC$age_cut=age_cut$age_cut
Idents(PBMC)="age_cut"
DEGs=FindAllMarkers(PBMC, min.pct=0, logfc.threshold=0.1, only.pos=TRUE)
DEGs$SYMBOL=PBMC@assays$RNA@meta.features$feature_name[match(DEGs$gene, rownames(PBMC))]
# get the meaningful genes among the sig. ones
meaningful=enrichGO(DEGs$SYMBOL, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="ALL")
as.data.frame(meaningful) %>%
  group_by(ONTOLOGY) %>%
  top_n(5) %>%
  View()
selected_genes=strsplit(as.data.frame(meaningful)$geneID[1], "/")[[1]]
# prepare the matrix
df=CPM_Matrix[match(DEGs$SYMBOL, rownames(CPM_Matrix)), ]
df_scaled=t(scale(t(df)))
# column annotation
sudocolor=paletteer_c("grDevices::Grays", 78); names(sudocolor)=rev(colnames(df_scaled))
group=HeatmapAnnotation(Age=anno_simple(colnames(df_scaled), col=sudocolor),
                        annotation_height=unit(2,"mm"),
                        annotation_label=gt_render(" Age", gp=gpar(fontsize=12)))
# gene annotation
index=which(rownames(df_scaled) %in% selected_genes)
rowannot=rowAnnotation(genemark=anno_mark(at=index,
                                          labels=selected_genes,
                                          labels_gp=gpar(fontsize=8),
                                          lines_gp=gpar()))
quantile(df_scaled, c(0.01,0.05,0.95,0.99)) # for setting matrix colors
cols=circlize::colorRamp2(c(-3,0,3), c("#5c4bb3","white","#b34b4b"))
# plot
p=
  Heatmap(df_scaled,
          col=cols,
          right_annotation=rowannot,
          top_annotation=group,
          
          cluster_columns=FALSE,
          show_column_dend=FALSE,
          column_title=NULL,
          show_column_names=FALSE,
          
          cluster_rows=TRUE,
          show_row_dend=FALSE,
          row_title=NULL,
          show_row_names=FALSE,
          
          name=" ")
pdf("~/Project_PBMCage/Plots/Heatmap_by_age_sum.pdf", width=0.05*nrow(df_scaled), height=0.05*ncol(df_scaled))
plot(p)
dev.off()

#####################################



### Sample clustering of AgeCut
#####################################
###

### Reduce the per_donor_matrix to per_age_matrix
temp=readRDS("~/Project_PBMCage/pseudoDGE_by_donor.rds")
df_by_age=t(rowsum(t(temp$counts), group=temp$samples$Age))

### Create a pseudo-bulk DGEList obj
pseudoDGE=DGEList(df_by_age)
metainfo=read.delim("~/Project_PBMCage/Results/MtRbCellcycle_Eval.txt", sep="\t", stringsAsFactors=FALSE)
pseudoDGE$samples$Age=colnames(df_by_age)
pseudoDGE$samples$Age=as.factor(pseudoDGE$samples$Age)
saveRDS(pseudoDGE, "~/Project_PBMCage/pseudoDGE_by_age.rds")

### Calculate cpm
CPM_Matrix=edgeR::cpm(pseudoDGE, log=TRUE, prior.count=1)

### Take only the aging-related genes based on DEMAGALHAES_AGING_UP Dataset
AgingDB=read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")
Aging_Matrix=CPM_Matrix[match(AgingDB$gene, rownames(CPM_Matrix)), ]

### Visualize the age effects
p=
Heatmap(t(scale(t(Aging_Matrix))),
        cluster_columns=F,
        clustering_method_rows="ward.D2",
        heatmap_legend_param=list(title=NULL),
        row_names_gp=gpar(fontsize=8),
        column_names_gp=gpar(fontsize=9))
# seems that age cuts can be set at 19-40, 41-70, 71-100
pdf("~/Project_PBMCage/Plots/Heatmap_for_AgeCut.pdf", height=0.15*nrow(Aging_Matrix), width=0.15*ncol(Aging_Matrix))
plot(p)
dev.off()

#####################################



### Evaluate the pool_number effect
#####################################
###

### Make the dataframe
Ages_DF=Donors_DF=list()
for (i in 1:length(levels(object$pool_number))) {
  ages_=names(table(split(as.character(object$age), as.character(object$pool_number))[i]))
  donors_=names(table(split(as.character(object$donor_id), as.character(object$pool_number))[i]))
  ages_df=data.table(batch=rep(levels(object$pool_number)[i], length(ages_)), age=ages_)
  donors_df=data.table(batch=rep(levels(object$pool_number)[i], length(donors_)), donor=donors_)
  Ages_DF=c(Ages_DF, list(ages_df))
  Donors_DF=c(Donors_DF, list(donors_df))
}
ages_df=rbindlist(Ages_DF)
donors_df=rbindlist(Donors_DF)

### Plot
ages_df$age=factor(ages_df$age, levels=sort(as.integer(names(table(ages_df$age)))))
p1=
  ggplot(ages_df, aes(x=batch, y=age)) +
  geom_point() +
  theme_classic()

# p2=
# ggplot(donors_df, aes(x=batch, y=donor)) +
#   geom_point() +
#   theme_bw()

pdf("~/Project_PBMCage/Plots/Age_Batch_Effect.pdf", width=12, height=10)
plot(p1)
dev.off()

#####################################



if (F) {
### Merge and Harmony based on AgeCuts, and a rough clustering
#####################################
###

### Harmony the samples per AgeCut
Sample_files=gsub(".rds","",list.files("~/Project_PBMCage/raw_data/preprocessedData/"))
Seurats_1940=Seurats_4170=Seurats_7197=list()
for (i in 1:length(Sample_files)) {
  temps=readRDS(file.path("~/Project_PBMCage/raw_data/preprocessedData", paste0(Sample_files[i],".rds")))
  age_=gsub("age_","",Sample_files[i])
  if (as.numeric(age_) %in% c(19:40)) {
    Seurats_1940=c(Seurats_1940, list(temps))
  } 
  else if (as.numeric(age_) %in% c(41:70)) {
    Seurats_4170=c(Seurats_4170, list(temps))
  }
  else if (as.numeric(age_) %in% c(71:100)) {
    Seurats_7197=c(Seurats_7197, list(temps))
  }
}
saveRDS(Seurats_1940, "~/Project_PBMCage/raw_data/Seurats_1940.rds")
saveRDS(Seurats_4170, "~/Project_PBMCage/raw_data/Seurats_4170.rds")
saveRDS(Seurats_7197, "~/Project_PBMCage/raw_data/Seurats_7197.rds")

### Split by pool_number, run SCT, and integrate
lapply(Seurats_1940, function(x) dim(x)) # check it is proper to run SCT instead of LogNorm
seurats=merge(Seurats_1940[[1]], Seurats_1940[-1])
seurats=SplitObject(seurats, split.by="pool_number")
seurats=lapply(seurats, SCTransform, return.only.var.genes=FALSE)
var.features=SelectIntegrationFeatures(seurats, nfeatures=3000)
seurats_=merge(seurats[[1]], seurats[-1])
VariableFeatures(seurats_)=var.features
seurats_=seurats_ %>%
  RunPCA(verbose=T) %>%
  RunHarmony(group.by.vars="pool_number", plot_convergence=TRUE)
saveRDS(seurats_, "~/Project_PBMCage/raw_data/Seurats_1940_harmony.rds")

### Pack as function and run on the other two Seurats List
Temp_function=
  function(Seurats_List) {
    lapply(Seurats_List, function(x) dim(x)) # check it is proper to run SCT instead of LogNorm
    seurats=merge(Seurats_List[[1]], Seurats_List[-1])
    seurats=SplitObject(seurats, split.by="pool_number")
    seurats=lapply(seurats, SCTransform, return.only.var.genes=FALSE)
    var.features=SelectIntegrationFeatures(seurats, nfeatures=3000)
    seurats_=merge(seurats[[1]], seurats[-1])
    VariableFeatures(seurats_)=var.features
    seurats_=seurats_ %>%
      RunPCA(verbose=T) %>%
      RunHarmony(group.by.vars="pool_number", plot_convergence=TRUE)
    saveRDS(seurats_, file.path("~/Project_PBMCage/raw_data", "_harmony.rds"))
  }

Temp_function(Seurats_4170)
file.rename("~/Project_PBMCage/raw_data/_harmony.rds","~/Project_PBMCage/raw_data/Seurats_4170_harmony.rds")
Temp_function(Seurats_7197)
file.rename("~/Project_PBMCage/raw_data/_harmony.rds","~/Project_PBMCage/raw_data/Seurats_7197_harmony.rds")

#####################################
}



### Sketch-based analysis
#####################################
###

### Create the sketch set
object=readRDS("~/Project_PBMCage/raw_data/Seurats_all.rds")
object=NormalizeData(object)
object=FindVariableFeatures(object, verbose=T)
# draw sketch
object=SketchData(object=object, ncells=50000, method="LeverageScore", sketched.assay="sketch")
saveRDS(object, "~/Project_PBMCage/raw_data/sketch.rds")

#####################################



### Download the reference for supervised annotation
#####################################
###

### 
library(SeuratData)
AvailableData()
# download
# wget http://seurat.nygenome.org/src/contrib/pbmc3k.SeuratData_3.1.4.tar.gz -O ~/SeuratData/pbmc3k.SeuratData.tar.gz
# install.packages("~/SeuratData/pbmc3k.SeuratData.tar.gz", repos=NULL, type="source")
# wget http://seurat.nygenome.org/src/contrib/pbmcMultiome.SeuratData_0.1.4.tar.gz -O ~/SeuratData/pbmcMultiome.SeuratData.tar.gz
# install.packages("~/SeuratData/pbmcMultiome.SeuratData.tar.gz", repos=NULL, type="source")
# wget http://seurat.nygenome.org/src/contrib/pbmcref.SeuratData_1.0.0.tar.gz -O ~/SeuratData/pbmcref.SeuratData.tar.gz
# install.packages("~/SeuratData/pbmcref.SeuratData.tar.gz", repos=NULL, type="source")

# pbmcsca.SeuratData_3.0.0

#####################################



### Annotate using supervised scPred
#####################################
###

# source("~/scPred/Build_ref.R")

### Load the scPred prediction result
query=readRDS("~/Project_PBMCage/raw_data/Seurats_all_scPredresults.rds")
# DimPlot(query, group.by="scpred_prediction", label=TRUE, label.size=3, repel=TRUE, raster=F) + NoLegend()

#####################################



### Supervised classification of the unannotated cells using TransferAnchors
#####################################
###

### Load the reference dataset
# have been processed
reference=readRDS("~/scPred/RefData/PBMC_multimodal_ref.rds")
reference
# add the preferable celltype list
table(reference$celltype.l1)
table(reference$celltype.l2)
table(reference$celltype.l3)
reference=subset(reference, celltype.l2!="Doublet")
cell_type=data.frame(cell=colnames(reference),
                     celltype.l1=reference$celltype.l1,
                     celltype.l2=reference$celltype.l2,
                     celltype.l3=reference$celltype.l3,
                     mycelltype=reference$celltype.l2)
rownames(cell_type)=cell_type$cell
cell_type=cell_type %>%
  mutate(mycelltype=ifelse(celltype.l2 %in% c("cDC1","cDC2"), "cDC", mycelltype)) %>%
  mutate(mycelltype=ifelse(celltype.l3=="Plasma", "Plasma", mycelltype)) %>%
  mutate(mycelltype=ifelse(celltype.l3=="Plasmablast", "Plasmablast", mycelltype)) %>%
  mutate(mycelltype=ifelse(celltype.l2 %in% c("CD4 CTL","CD4 Proliferating","CD4 TCM","CD4 TEM"), "Experienced CD4 T", mycelltype)) %>%
  mutate(mycelltype=ifelse(celltype.l2 %in% c("CD8 Proliferating","CD8 TCM","CD8 TEM"), "Experienced CD8 T", mycelltype))
table(cell_type$mycelltype) # check
table(names(reference$celltype.l1)==cell_type$cell) # check the order
reference[["mycelltype"]]=cell_type$mycelltype

### Process the query dataset
# have been normalized
query=readRDS("~/Project_PBMCage/raw_data/Seurats_all_scPredresults.rds")
object=subset(query, scpred_prediction=="unassigned")

### Transfer the annotation to the query dataset
anchors=FindTransferAnchors(reference=reference,
                            query=object,
                            dims=1:30,
                            reference.reduction="pca",
                            verbose=T)
predictions=TransferData(anchorset=anchors,
                         refdata=reference$mycelltype,
                         dims=1:30)
write.table(predictions, "~/Project_PBMCage/Results/remaining_cells_PredictResults.txt", sep="\t")

### Assign the cells
table(predictions$predicted.id)
query_annot=data.frame(cell=names(query$scpred_prediction), scpred_prediction=query$scpred_prediction)
re_annot=data.frame(cell=rownames(predictions), re_annotation=predictions$predicted.id)
com_annot=merge(query_annot, re_annot, by="cell", all=T)
com_annot=com_annot[match(query_annot$cell, com_annot$cell),]
com_annot=com_annot %>%
  mutate(final_celltype_annotation=ifelse(scpred_prediction=="unassigned", re_annotation, scpred_prediction))
# add metadata
table(com_annot$cell==colnames(query)) # check
query[["final_celltype_annotation"]]=com_annot$final_celltype_annotation
table(query$final_celltype_annotation)
write.table(com_annot, "~/Project_PBMCage/Results/all_cells_annoted.txt", sep="\t")
saveRDS(query, "~/Project_PBMCage/raw_data/Seurats_all_annotated.rds")

#####################################



### This is the final annotation for celltypes after subsets unsupervised clustering
# Steps that have been done:
# 1) scPred, end up with a portion of unassigned cells
# 2) Anchor-based Label Transfer for those unassigned cells
# Steps to be done:
# 3) For each celltype, run unsupervised clustering
# 4) FindAllMarkers to get cell markers
# 5) Updated annotation is then saved in another vector

### Read the obj after scPred + Label Transfer (a rough supervised annotation)
query=readRDS("~/Project_PBMCage/raw_data/Seurats_all_annotated.rds")
table(query$final_celltype_annotation)

### Reassign B naive + B memory + B intermediate
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("B memory","B naive","B intermediate"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(10,12), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("B_","B_","B_","B_","B_","dnT","B_","NKT","B_","B_","Megakaryocyte","CD14 Mono")

### Take away non-B cells and cluster again
sub_obj2=subset(sub_obj_try, reassign=="B_")
sub_obj2=sub_obj2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj2)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try2=FindClusters2(sub_obj2, cluster.range=c(5,8), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try2, label=T, label.size=5, repel=T) + NoLegend()

markers2=FindAllMarkers(sub_obj_try2, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_2)

sub_obj_try2[["reassign"]]=sub_obj_try2$seurat_clusters
levels(sub_obj_try2$reassign)=c("0 B naive","1 Switched B memory","2 Switched B memory","3 Non-classical B memory","4 Switched B memory","5 B naive lambda","6 B intermediate","7 ISG+ switched B memory")

# split "5 B naive lambda"
sub_obj3=subset(sub_obj_try2, reassign=="5 B naive lambda")
sub_obj3=sub_obj3 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj3)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try3=FindClusters2(sub_obj3, cluster.range=2, by=0.1, res=0.1, verbose=T)
DimPlot(sub_obj_try3, label=T, label.size=5, repel=T) + NoLegend()

markers3=FindAllMarkers(sub_obj_try3, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_3=markers3 %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_3)

sub_obj_try3[["reassign"]]=sub_obj_try3$seurat_clusters
levels(sub_obj_try3$reassign)=c("B naive lambda","B memory lambda")

### Combine the multi-degree clustering results
table(sub_obj_try$reassign)
table(sub_obj_try2$reassign)
table(sub_obj_try3$reassign)
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
reassignment=reassignment %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("dnT"))), "dnT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("NKT"))), "NKT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("Megakaryocyte"))), "Megakaryocyte", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("CD14 Mono"))), "CD14 Mono", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("0 B naive"))), "B naive", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("1 Switched B memory","2 Switched B memory","4 Switched B memory"))), "Switched B memory", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("3 Non-classical B memory"))), "Non-classical B memory", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("6 B intermediate"))), "B intermediate", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("7 ISG+ switched B memory"))), "ISG+ switched B memory", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try3, reassign %in% c("B naive lambda"))), "B naive lambda", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try3, reassign %in% c("B memory lambda"))), "B memory lambda", reassign))

table(reassignment$reassign, useNA="ifany") # check

# modify the names
reassignment=reassignment %>%
  mutate(reassign=ifelse(reassign %in% c("B memory lambda","Switched B memory"), "B memory", reassign)) %>%
  mutate(reassign=ifelse(reassign %in% c("B naive lambda","B naive"), "B naive", reassign)) %>%
  mutate(reassign=ifelse(reassign=="ISG+ switched B memory", "ISG+ B memory", reassign))
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_B_memory.txt", sep="\t")

#####################################



### Reassign Experienced CD4 T
#####################################
###

sub_obj=subset(query, final_celltype_annotation=="Experienced CD4 T")
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(6,8), by=0.1, res=0.25, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

VlnPlot(sub_obj_try, features=c("GZMK","CD3D","CD3E","CD3G","CD4","CD8A","CD8B"), pt.size=0)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("CD4 Tcm","CD4 Tem","CD4 Tem","CD4 CTL","ISAGhi T","Megakaryocyte","6 CD4 Tcm","CD8 Tcm")
# * Cluster 0: CD4 TCM_3

### Split Cluster 6: CD4 Tcm
sub_obj2=subset(sub_obj_try, reassign=="6 CD4 Tcm")
sub_obj2=sub_obj2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try2=FindClusters2(sub_obj2, cluster.range=3, by=0.1, res=0.25, verbose=T)
DimPlot(sub_obj_try2, label=T, label.size=5, repel=T) + NoLegend()

markers2=FindAllMarkers(sub_obj_try2, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_2)

sub_obj_try2[["reassign"]]=sub_obj_try2$seurat_clusters
levels(sub_obj_try2$reassign)=c("NKT","B naive lambda","CD4 Tcm")

### Split Cluster 1 + 2 : CD4 Tem
sub_obj3=subset(sub_obj_try, reassign=="CD4 Tem")
sub_obj3=sub_obj3 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try3=FindClusters2(sub_obj3, cluster.range=c(3,6), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try3, label=T, label.size=5, repel=T) + NoLegend()

markers3=FindAllMarkers(sub_obj_try3, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_3=markers3 %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_3)

sub_obj_try3[["reassign"]]=sub_obj_try3$seurat_clusters
levels(sub_obj_try3$reassign)=c("CD4 Tem","CD4 Tem","CD4 Tem","CD4 Tem","CD4 Tem","Cycling CD4 T")

### Combine the multi-degree clustering results
table(sub_obj_try$reassign)
table(sub_obj_try2$reassign)
table(sub_obj_try3$reassign)
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
reassignment=reassignment %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try3, reassign %in% c("CD4 Tem"))), "CD4 Tem", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try3, reassign %in% c("Cycling CD4 T"))), "Cycling CD4 T", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("NKT"))), "NKT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("B naive lambda"))), "B naive lambda", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("CD4 Tcm"))), "CD4 Tcm", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("CD4 Tcm"))), "CD4 Tcm", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("CD4 CTL"))), "CD4 CTL", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("ISAGhi T"))), "ISAGhi T", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("Megakaryocyte"))), "Megakaryocyte", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("CD8 Tcm"))), "CD8 Tcm", reassign))

table(reassignment$reassign, useNA="ifany") # check
View(reassignment)
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_Experienced_CD4_T.txt", sep="\t")

#####################################



### Reassign Experienced CD8 T
#####################################
###

sub_obj=subset(query, final_celltype_annotation=="Experienced CD8 T")
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(10,15), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

VlnPlot(sub_obj_try, features=c("GZMB"), pt.size=0)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("CD8 Trm","CD8 Temra","GZMB+ CD8 Tem","CD8 Tcm","CD8 NKT-like","Eryth","HLA-DR+ CD8 T","Megakaryocyte","B memory","Exhausted CD8 T","CD8 Tcm","dnT")
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_Experienced_CD8_T.txt", sep="\t")

#####################################



### Reassign NK
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("NK","NK Proliferating","NK_CD56bright"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(9,12), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("NKT","CD56dim NK","NKT","CD56hi NK","Eryth","dnT","NKT","Megakaryocyte","Immature NK","B naive")
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_NK.txt", sep="\t")

#####################################



### Reassign DC
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("ASDC","cDC","pDC"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(10,15), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

markers2=FindMarkers(sub_obj_try, ident.1="0",ident.2="5", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("cDC2","pDC","pre-cDC","pre-pDC","AXL+ DC","cDC1","Megakaryocyte","B naive lambda","Eryth","B naive lambda")
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_DC.txt", sep="\t")

#####################################



### Reassign Mono
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("CD14 Mono","CD16 Mono"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(5,8), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("CD14+ Mono","CD14+ Mono","CD16+ Mono","CD14+ Mono","Megakaryocyte","CD16+ Mono","B naive")

reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_Mono.txt", sep="\t")

#####################################



### Reassign Treg
#####################################
###

sub_obj=subset(query, final_celltype_annotation=="Treg")
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(5,6), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("CD4 naive","Naive Treg","CD4 CTL","CD4 CTL","B memory","Megakaryocyte")

reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_Treg.txt", sep="\t")

#####################################



### Reassign CD4 naive
#####################################
###

sub_obj=subset(query, final_celltype_annotation=="CD4 Naive")

sub_obj=sub_obj %>% 
  NormalizeData() %>%
  FindVariableFeatures()

sub_obj=SketchData(sub_obj, ncells=50000, method="LeverageScore", sketched.assay="sketch")
DefaultAssay(sub_obj)="sketch"
sub_obj=sub_obj %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50, return.model=T)
DimPlot(sub_obj)
sub_obj_try=FindClusters(sub_obj, resolution=0.2, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()
# ---take a look---
markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
if (F) {
  sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
  levels(sub_obj_try$reassign)=c("CD4 Tnaive","CD4 Tnaive","CD4 Tnaive","CD8 NKT-like","Eryth","ISAGhi T","Megakaryocyte")
}
# ---   end    ---
sub_obj_try=ProjectData(
  sub_obj_try,
  assay="RNA",
  full.reduction="pca.full",
  sketched.assay="sketch",
  sketched.reduction="pca",
  umap.model="umap",
  dims=1:50,
  refdata=list(seurat_clusters_full="seurat_clusters")
)

DefaultAssay(sub_obj_try)="RNA"
Idents(sub_obj_try)="seurat_clusters_full"
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters_full
sub_obj_try$reassign=as.factor(sub_obj_try$reassign)
levels(sub_obj_try$reassign)=c("0 CD4 Tnaive","1 CD4 Tnaive","10 CD4 Tnaive","2 NKT","3 NKT","4 Eryth","5 ISAGhi T","6 Megakaryocyte","7 CD4 Tnaive","8 NKT")
levels(sub_obj_try$reassign)=c("CD4 Tnaive","CD4 Tnaive","CD4 Tnaive","NKT","NKT","Eryth","ISAGhi T","Megakaryocyte","CD4 Tnaive","NKT")

reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_CD4Tnaive.txt", sep="\t")

#####################################



### Reassign CD8 Naive
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("CD8 Naive"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(5,8), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

markers2=FindMarkers(sub_obj_try, ident.1="0",ident.2="5", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("CD8 Tnaive","CD8 Tnaive","CD8 Tnaive","CD8 Tnaive","CD8 NKT-like","Eryth","ISAGhi T","B naive")
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_CD8Tnaive.txt", sep="\t")

#####################################



### Reassign gdT, dnT, HSPC, ILC, MAIT
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("gdT","dnT","HSPC","ILC","MAIT"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(6,10), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

markers2=FindMarkers(sub_obj_try, ident.1="0",ident.2=c("1","2"), only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("gdT","MAIT","gdT","3","HSPC","GZMK+ CD8 Tem","HSPC","ILC","Megakaryocyte")

### Split Cluster 3
sub_obj2=subset(sub_obj_try, seurat_clusters=="3")
sub_obj2=sub_obj2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)
DimPlot(sub_obj2)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try2=FindClusters2(sub_obj2, cluster.range=c(3,5), by=0.1, res=0.25, verbose=T)
DimPlot(sub_obj_try2, label=T, label.size=5, repel=T) + NoLegend()

markers2=FindAllMarkers(sub_obj_try2, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_2)

sub_obj_try2[["reassign"]]=sub_obj_try2$seurat_clusters
levels(sub_obj_try2$reassign)=c("CD4 Tnaive","CD4 Tcm","MAIT","CD4 Tcm")

### Combine the multi-degree clustering results
table(sub_obj_try$reassign)
table(sub_obj_try2$reassign)
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
reassignment=reassignment %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("gdT"))), "gdT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("MAIT"))), "MAIT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("HSPC"))), "HSPC", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("GZMK+ CD8 Tem"))), "GZMK+ CD8 Tem", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("ILC"))), "ILC", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("Megakaryocyte"))), "Megakaryocyte", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("CD4 Tnaive"))), "CD4 Tnaive", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("CD4 Tcm"))), "CD4 Tcm", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("MAIT"))), "MAIT", reassign))
table(reassignment$reassign, useNA="ifany") # check
View(reassignment)
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_gdT_dnT_HSPC_ILC_MAIT.txt", sep="\t")

#####################################



### Reassign Plasma, Plasmablast
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("Plasma","Plasmablast"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(6,10), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("Plasma cell","Plasma cell","CD4 Tem","Plasmablast","Plasma cell","B intermediate","CD56dim NK")

reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
View(reassignment)
table(reassignment$reassign, useNA="ifany")
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_PB_PC.txt", sep="\t")

#####################################



### Reassign Eryth, Platelet
#####################################
###

sub_obj=subset(query, final_celltype_annotation %in% c("Eryth","Platelet"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(6,10), by=0.1, res=0.3, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

markers2=FindMarkers(sub_obj_try, ident.1="4", ident.2="3", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("0","1 CD4 Tnaive","2 CD8 NKT-like","3 Megakaryocyte","4 Eryth","5 CD14+ Mono","6","7 NK")

### Split Cluster 0 + 6
sub_obj2=subset(sub_obj_try, seurat_clusters %in% c("0","6"))
sub_obj2=sub_obj2 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)
DimPlot(sub_obj2)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try2=FindClusters2(sub_obj2, cluster.range=c(3,5), by=0.1, res=0.25, verbose=T)
DimPlot(sub_obj_try2, label=T, label.size=5, repel=T) + NoLegend()

markers3=FindAllMarkers(sub_obj_try2, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_3=markers3 %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_3)

sub_obj_try2[["reassign"]]=sub_obj_try2$seurat_clusters
levels(sub_obj_try2$reassign)=c("B naive","Platelet","CD4 Tcm")

markers4=FindMarkers(sub_obj_try2, ident.1="2", ident.2="1", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_4=markers4 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_4)

### Combine the multi-degree clustering results
table(sub_obj_try$reassign)
table(sub_obj_try2$reassign)
reassignment=data.frame(cell=colnames(sub_obj_try), reassign=sub_obj_try$reassign)
reassignment=reassignment %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("1 CD4 Tnaive"))), "CD4 Tnaive", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("2 CD8 NKT-like"))), "CD8 NKT-like", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("3 Megakaryocyte"))), "Megakaryocyte", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("4 Eryth"))), "Eryth", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("5 CD14+ Mono"))), "CD14+ Mono", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign %in% c("7 NK"))), "NK", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("B naive"))), "B naive", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("Platelet"))), "Platelet", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try2, reassign %in% c("CD4 Tcm"))), "CD4 Tcm", reassign))
table(reassignment$reassign, useNA="ifany") # check
View(reassignment)
write.table(reassignment, "~/Project_PBMCage/Results/reassignment_Eryth_Platelet.txt", sep="\t")

#####################################



### Combine everything
#####################################
###

files=list.files("~/Project_PBMCage/Results/")[grepl("^reassignment_", list.files("~/Project_PBMCage/Results/"))]
fileLists=list()
for (i in 1:length(files)) {
  fileLists[[i]]=read.delim(file.path("~/Project_PBMCage/Results",files[i]))
}
names(fileLists)=gsub("reassignment_|\\.txt$","",files)
fileComb=data.table::rbindlist(fileLists)
nrow(fileComb) # check
query # check

# edit the names to remove duplicates
names(table(fileComb$reassign)) # check
fileComb=fileComb %>%
  mutate(reassign=ifelse(reassign=="CD14+ Mono","CD14 Mono", reassign)) %>%
  mutate(reassign=ifelse(reassign=="CD16+ Mono","CD16 Mono", reassign)) %>%
  mutate(reassign=ifelse(reassign=="CD4 naive","CD4 Tnaive", reassign)) %>%
  mutate(reassign=ifelse(reassign=="NK","CD56dim NK", reassign)) %>%
  mutate(reassign=ifelse(reassign=="NK","CD56dim NK", reassign)) %>%
  mutate(reassign=ifelse(reassign=="Megakaryocyte","Platelet", reassign)) %>%
  mutate(reassign=ifelse(reassign=="Eryth","Erythroid cell", reassign))

names(table(fileComb$reassign)) # check

# add metadata
fileComb_ordered=fileComb[match(colnames(query), fileComb$cell),]
head(fileComb_ordered) # check
head(colnames(query)) # check
reassigned_annotation=fileComb_ordered$reassign
names(reassigned_annotation)=fileComb_ordered$cell
query=AddMetaData(query, metadata=reassigned_annotation, col.name="reassigned")

names(table(query$reassigned)) # final check
DimPlot(query, label=T, label.size=3, repel=T, group.by="reassigned") + NoLegend()

#####################################



### Modify strange clustering
#####################################
###

###
sub_obj=subset(query, reassigned %in% c("ILC"))
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:30) %>%
  RunUMAP(dims=1:30)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=c(2,3), by=0.1, res=0.1, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)

markers2=FindMarkers(sub_obj_try, ident.1="4", ident.2="3", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("NK","NK")

###
sub_obj=subset(query, reassigned %in% c("NKT"))

sub_obj=sub_obj %>% 
  NormalizeData() %>%
  FindVariableFeatures()

sub_obj=SketchData(sub_obj, ncells=50000, method="LeverageScore", sketched.assay="sketch")
DefaultAssay(sub_obj)="sketch"
sub_obj=sub_obj %>% 
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50, return.model=T)
DimPlot(sub_obj)
sub_obj_try=FindClusters(sub_obj, resolution=0.2, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()
# ---take a look---
markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=20, order_by=avg_log2FC)
View(top10)
if (F) {
  sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
  levels(sub_obj_try$reassign)=c("0 NK","1 NK","2 CD4 Tcm","3 NKT","4 Proliferating NK","5 B naive","6 CD14 Mono","7 NK")
}
# ---   end    ---
sub_obj_try=ProjectData(
  sub_obj_try,
  assay="RNA",
  full.reduction="pca.full",
  sketched.assay="sketch",
  sketched.reduction="pca",
  umap.model="umap",
  dims=1:50,
  refdata=list(seurat_clusters_full="seurat_clusters")
)

DefaultAssay(sub_obj_try)="RNA"
Idents(sub_obj_try)="seurat_clusters_full"
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers2=FindMarkers(sub_obj_try, ident.1="6", ident.2=c("1","0"), only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10_2=markers2 %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10_2)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters_full
sub_obj_try$reassign=as.factor(sub_obj_try$reassign)
levels(sub_obj_try$reassign)=c("0 NK","1 NK","2 CD4 Tcm","3 NKT","4 Proliferating NK","5 B naive","6 CD14 Mono")
levels(sub_obj_try$reassign)=c("NK","NK","CD4 Tcm","NKT","Proliferating NK","B naive","CD14 Mono")

### modifty
fileComb=fileComb %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="NK")), "CD56dim NK", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="CD4 Tcm")), "CD4 Tcm", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="NKT")), "NKT", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="Proliferating NK")), "Proliferating NK", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="B naive")), "B naive", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="CD14 Mono")), "CD14 Mono", reassign))

names(table(fileComb$reassign)) # check

# add metadata
fileComb_ordered=fileComb[match(colnames(query), fileComb$cell),]
head(fileComb_ordered) # check
head(colnames(query)) # check
reassigned_annotation=fileComb_ordered$reassign
names(reassigned_annotation)=fileComb_ordered$cell
query=AddMetaData(query, metadata=reassigned_annotation, col.name="reassigned")

names(table(query$reassigned)) # final check
DimPlot(query, label=T, label.size=2, repel=T, group.by="reassigned", raster=F) + NoLegend()

write.table(fileComb_ordered, "~/Project_PBMCage/Results/reassignment_Combined.txt", sep="\t")
saveRDS(query, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned.rds")

### Check again in the umap for final modification
# Idents(query)=query$reassigned
clusters=names(table(query_modified$reassigned));clusters
cluster_to_choose=c("Bnaive")
highlighted=WhichCells(query_modified, idents=cluster_to_choose)
DimPlot(query_modified, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="transparent") + NoLegend()
sub_obj=subset(query_modified, reassigned %in% cluster_to_choose)
sub_obj=sub_obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors(dims=1:50) %>%
  RunUMAP(dims=1:50)
DimPlot(sub_obj)
source("~/Rscripts/FindCluster2_Functions.R")
sub_obj_try=FindClusters2(sub_obj, cluster.range=3, by=0.1, res=0.1, verbose=T)
DimPlot(sub_obj_try, label=T, label.size=5, repel=T) + NoLegend()

markers=FindAllMarkers(sub_obj_try, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
top10=markers %>%
  group_by(cluster) %>%
  slice_max(n=10, order_by=avg_log2FC)
View(top10)

sub_obj_try[["reassign"]]=sub_obj_try$seurat_clusters
levels(sub_obj_try$reassign)=c("pDC","pDC","CLP")
DimPlot(sub_obj_try, label=T, label.size=5, repel=T, group.by="reassign") + NoLegend()

# modifty
# fileComb=read.delim("~/Project_PBMCage/Results/reassignment_Combined.txt")
fileComb_modified=read.delim("~/Project_PBMCage/Results/all_cells_annoted_final.txt")
fileComb_modified=fileComb_modified %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="pDC")), "pDC", reassign)) %>%
  mutate(reassign=ifelse(cell %in% colnames(subset(sub_obj_try, reassign=="CLP")), "CLP", reassign))

names(table(fileComb_modified$reassign)) # check
write.table(fileComb_modified, "~/Project_PBMCage/Results/all_cells_annoted_final.txt", sep="\t")
# add metadata
fileComb_ordered=fileComb_modified[match(colnames(query_modified), fileComb_modified$cell),]
head(fileComb_ordered) # check
head(colnames(query_modified)) # check
reassigned_annotation=fileComb_ordered$reassign
names(reassigned_annotation)=fileComb_ordered$cell
query_modified=AddMetaData(query_modified, metadata=reassigned_annotation, col.name="reassigned")

names(table(query_modified$reassigned)) # final check
Idents(query_modified)="reassigned"
highlighted=WhichCells(query_modified, idents="dnT")
DimPlot(query_modified, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="gray80") + NoLegend()
subset_=subset(query_modified, reassigned %in% cluster_to_choose)
DimPlot(subset_, label=T, label.size=5, repel=T) + NoLegend()

saveRDS(query_modified, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")

plists=list()
clusters=names(table(query_modified$reassigned))
for (i in 1:length(clusters)) {
  highlighted=WhichCells(query_modified, idents=clusters[i])
  plists[[i]]=
    DimPlot(query_modified, label=T, label.size=2, repel=T, group.by="reassigned", cells.highlight=list(Cluster=highlighted), cols.highlight=c("darkred"), cols="transparent") + NoLegend()
}
for (i in 1:length(plists)) {
  pdf(file.path("~/Project_PBMCage/a_temp", paste0(i, ".pdf")))
  plot(plists[[i]])
  dev.off()  
}

#####################################



### Group the celltypes
#####################################
###

###
table(query_modified$reassigned)
fileComb_modified=read.delim("~/Project_PBMCage/Results/all_cells_annoted_final.txt")
fileComb_interm=fileComb_modified %>%
  mutate(intermediate=ifelse(reassign %in% c("cDC1","cDC2"), "cDC", reassign)) %>%
  mutate(intermediate=ifelse(reassign %in% c("Bnaive","lambda+ Bnaive"), "Bnaive", intermediate)) %>%
  mutate(intermediate=ifelse(reassign %in% c("Classical Bmem","DN Bmem","IgDdim IgG Bmem","IgDhi IgG Bmem"), "Bmem", intermediate)) %>%
  mutate(intermediate=ifelse(reassign %in% c("CLP","EMP"), "HSPC", intermediate)) %>%
  mutate(intermediate=ifelse(reassign %in% c("GZMB+ CD8 Tem","GZMK+ CD8 Tem"), "CD8 Tem", intermediate)) %>%
  mutate(intermediate=ifelse(reassign %in% c("IgA+ PC","IgG+ PC"), "PC", intermediate))
table(fileComb_interm$intermediate)
write.table(fileComb_interm, "~/Project_PBMCage/Results/all_cells_annoted_final_WithIntermAnnot.txt")

fileComb_rough=fileComb_interm %>%
  mutate(rough=ifelse(intermediate %in% c("AXL+ DC","cDC","pDC"), "DC", intermediate)) %>%
  mutate(rough=ifelse(intermediate %in% c("CD4 Tcm","CD4 Tem"), "CD4 Tmem", rough)) %>%
  mutate(rough=ifelse(intermediate %in% c("CD8 Tcm","CD8 Tem"), "CD8 Tmem", rough)) %>%
  mutate(rough=ifelse(intermediate %in% c("dnT","gdT","ISAGhi T","MAIT","NKT","Treg"), "Other T", rough)) %>%
  mutate(rough=ifelse(intermediate %in% c("Erythroid cell","Platelet","Mast cell","ILC"), "Other cells", rough)) %>%
  mutate(rough=ifelse(intermediate %in% c("Proliferating NK","Proliferating T"), "Proliferating NK/T", rough)) %>%
  mutate(rough=ifelse(intermediate %in% c("PB","PC"), "PB/PC", rough))
table(fileComb_rough$rough)
write.table(fileComb_rough, "~/Project_PBMCage/Results/all_cells_annoted_final_WithIntermAndRoughAnnot.txt")

head(fileComb_interm) # check
colnames(query_modified)[1:6] # check
query_modified[["reassigned_interm"]]=fileComb_interm$intermediate
head(fileComb_rough) # check
colnames(query_modified)[1:6] # check
query_modified[["reassigned_rough"]]=fileComb_rough$rough

saveRDS(query_modified, "~/Project_PBMCage/raw_data/Seurats_all_annotated_reassigned_final.rds")

#####################################


