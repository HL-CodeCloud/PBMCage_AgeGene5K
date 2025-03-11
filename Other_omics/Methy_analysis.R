
### Analyze age_methylation data from EWAS_Atlas
#####################################
###

# ### Prepare the age_methylation data from the EWAS-Atlas
# load("/home/jo_79/Project_PBMCage/Other_omics/age_methylation_v1.RData")
# Tissues=as.character(age_download["tissue",])
# table(Tissues)
# tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
# processbar=txtProgressBar(min=0, max=length(tissues_selected), style=3, width=60, char="=")
# for (i in 1:length(tissues_selected)) {
#   tissue_=tissues_selected[i]
#   EWASage_PBMC=age_download[, as.character(age_download["tissue",])==tissue_]
#   saveRDS(EWASage_PBMC, paste0("/home/jo_79/Project_PBMCage/Other_omics/age_methylation_",gsub("\\+| ","",tissue_),".rds"))
#   setTxtProgressBar(processbar, i)
# }
# close(processbar)

### DMP/DMR analysis on age_methylation in each tissue selected
library(dplyr)
library(ChAMP)
library(DMRcate)

### Load the whole MetaData
load("~/Project_PBMCage/Other_omics/age_methylation_v1_sample_age.RData")
head(sample_age)

###
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
file_name=gsub("\\+| ","",tissues_selected)

METAOBJs=list()
processbar=txtProgressBar(min=0, max=length(file_name), style=3, width=60, char="=")
for (f in 1:length(file_name)) {
  MethyData=readRDS(paste0("~/Project_PBMCage/Other_omics/age_methylation_",file_name[f],".rds"))
  message("--- MethyData loaded. ---")

  # process metadata
  meta=sample_age[match(colnames(MethyData),rownames(sample_age)), ]
  table(colnames(MethyData)==rownames(meta)) # check
  meta=subset(meta, is.na(infection)) # remove infected individuals
  meta=meta[, colSums(is.na(meta))!=nrow(meta)]
  table(meta$platform) # check the platform
  colnames(meta)=c("Sample_Name",colnames(meta)[2:ncol(meta)])

  # process beta value data
  BetaValue=MethyData[4:nrow(MethyData),]
  beta=BetaValue[, colnames(BetaValue) %in% rownames(meta)]

  # create the obj
  MetaAgeObj=list(beta=beta, pd=meta)
  METAOBJs[[f]]=MetaAgeObj

  # preprocess the obj
  MetaAgeObj$beta=mutate_all(MetaAgeObj$beta, function(x) as.numeric(as.character(x)))
  MetaAgeObj$beta=as.matrix(MetaAgeObj$beta)
  MetaAgeObj=champ.filter(beta=MetaAgeObj$beta, pd=MetaAgeObj$pd)
  MetaAgeObj=champ.impute(beta=MetaAgeObj$beta, pd=MetaAgeObj$pd, method="Combine", SampleCutoff=0.2) # original samplecutoff=0.1 but too strict threshold

  # Calculate DMPs
  myDMP=champ.DMP(beta=MetaAgeObj$beta,
                  pheno=MetaAgeObj$pd$age,
                  compare.group=NULL,
                  adjPVal=0.05,
                  adjust.method="BH",
                  arraytype="450K")

  # Calculate DMRs
  source("~/Rscripts/champ.DMR2.R")
  myDMR=champ.DMR2(beta=MetaAgeObj$beta, pheno=MetaAgeObj$pd$age, method="DMRcate")

  # GSEA analysis
  myGSEA=champ.GSEA(beta=MetaAgeObj$beta, pheno=MetaAgeObj$pd$age, DMP=myDMP[[1]], DMR=myDMR, arraytype="450K", adjPval=0.05, method="gometh", Rplot=FALSE)
  # myebayGSEA=champ.GSEA(beta=MetaAgeObj$beta, pheno=MetaAgeObj$pd$age, arraytype="450K", adjPval=0.05, method="ebayes")

  # save
  save(myDMP, myDMR, myGSEA, file=paste0("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_",file_name[f],"_results.RData"))


  message(paste0("RData till ", file_name[f], " has been saved."))
  setTxtProgressBar(processbar, f)
}
close(processbar)

names(METAOBJs)=file_name
saveRDS(METAOBJs, "~/Project_PBMCage/Other_omics/AgeMeth_data.rds")

#####################################



### Extract the top DMPs (all feat.cgi in general) and plot their cor with ages
#####################################
###

###
library(dplyr)
library(tidyverse)
library(ggplot2)

### Make a function to extract and arrange the methylation data based on the selected DMPs
Beta_LongDF=function(beta, meta, DMP) {
  beta_extract=as.data.frame(t(beta[rownames(DMP), ]))
  beta_extract$age=meta$age[match(rownames(beta_extract), rownames(meta))]
  beta_extract=mutate_all(beta_extract, function(x) as.numeric(as.character(x)))
  beta_extract=beta_extract %>% arrange(age)

  # take the mean of the same age
  longdf=tidyr::pivot_longer(beta_extract, cols=colnames(beta_extract)[1:ncol(beta_extract)-1], names_to="probeid", values_to="beta")
  longdf$age=as.numeric(longdf$age)
  longdf$beta=as.numeric(longdf$beta)

  return(longdf)
}

Beta_WideTable=function(longdf) {
  # longdf is the result from Beta_LongDF(beta, meta, DMP)
  df=longdf %>%
    group_by(age, probeid) %>%
    summarize_at("beta", mean, na.rm=T)

  # rearrange the matrix as a wide table
  widedf=tidyr::pivot_wider(df, names_from="probeid", values_from="beta") %>% column_to_rownames("age")
  widedf=t(widedf)
  widedf=mutate_all(as.data.frame(widedf), function(x) as.numeric(as.character(x)))

  return(widedf)
}

### Load the methylation data
Methylation_Data=readRDS("~/Project_PBMCage/Other_omics/AgeMeth_data.rds")

### Analyze the pos. and neg. age-associated DMPs (all feat.cgi in general) in each celltype
LongDF_up_list=WideTable_up_list=LongDF_down_list=WideTable_down_list=list()
for (i in 1:length(Methylation_Data)) {
  celltype=names(Methylation_Data)[i]

  # load the DMP result
  results_per_celltype=load(paste0("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_",celltype,"_results.RData"))
  DMP_=myDMP[[1]]

  # load the methylation levels and metadata
  beta_=Methylation_Data[[i]]$beta
  meta_=Methylation_Data[[i]]$pd

  # take the top sig. DMPs
  threshold_=quantile(DMP_$adj.P.Val, c(0,0.1))[2]
  myDMP_up=DMP_ %>% subset(adj.P.Val<threshold_ & logFC>0) %>% slice_max(logFC, n=50)
  myDMP_down=DMP_ %>% subset(adj.P.Val<threshold_ & logFC<0) %>% slice_min(logFC, n=50)

  # get the long dataframe and wide table with the sig. upregulated DMPs
  longdf=Beta_LongDF(beta=beta_, meta=meta_, DMP=myDMP_up)
  LongDF_up_list[[i]]=longdf
  WideTable_up_list[[i]]=Beta_WideTable(longdf)

  # get the long dataframe and wide table with the sig. downreg DMPs
  longdf=Beta_LongDF(beta=beta_, meta=meta_, DMP=myDMP_down)
  LongDF_down_list[[i]]=longdf
  WideTable_down_list[[i]]=Beta_WideTable(longdf)
}
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
names(LongDF_up_list)=names(WideTable_up_list)=names(LongDF_down_list)=names(WideTable_down_list)=tissues_selected

### Plot the pos. and neg. age-associated DMPs (all feat.cgi in general) in each celltype
for (i in 1:length(LongDF_up_list)) {
  LongDF_up=LongDF_up_list[[i]]
  LongDF_down=LongDF_down_list[[i]]
  celltype=names(LongDF_up_list)[i]

  plot_up_PerCpG=
    ggplot(LongDF_up, aes(x=age, y=beta, color=probeid)) +
    # geom_boxplot(alpha=0.05) +
    theme_classic() +
    labs(x="age", y="Beta-value") +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0(celltype, " : Pos. assoc. w/ age")) +
    geom_smooth(show.legend=FALSE, method='loess', linewidth=0.5, fill="gray93", alpha=0.5) +
    ggpubr::stat_cor(aes(x=age, y=beta), show.legend=FALSE, size=1.5, label.x.npc="left", label.y.npc="top", method="spearman")

  plot_down_PerCpG=
    ggplot(LongDF_down, aes(x=age, y=beta, color=probeid)) +
    # geom_point(alpha=0.05) +
    theme_classic() +
    labs(x="age", y="Beta-value") +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0(celltype, " : Neg. assoc. w/ age")) +
    geom_smooth(show.legend=FALSE, method='loess', linewidth=0.5, fill="gray93", alpha=0.5) +
    ggpubr::stat_cor(aes(x=age, y=beta), show.legend=FALSE, size=1.5, label.x.npc="left", label.y.npc="top", method="spearman")

  # plot all sig. CpGs in total
  LongDF_up_all=LongDF_up %>%
    group_by(age) %>%
    summarize_at("beta", mean, na.rm=T)
  LongDF_up_all$age=as.numeric(LongDF_up_all$age)
  LongDF_up_all$beta=as.numeric(LongDF_up_all$beta)

  LongDF_down_all=LongDF_down %>%
    group_by(age) %>%
    summarize_at("beta", mean, na.rm=T)
  LongDF_down_all$age=as.numeric(LongDF_down_all$age)
  LongDF_down_all$beta=as.numeric(LongDF_down_all$beta)

  plot_up_AllSig.CpG=
    ggplot(LongDF_up_all, aes(x=age, y=beta)) +
    geom_point(alpha=1) +
    theme_classic() +
    labs(x="age", y="Beta-value") +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0(celltype)) +
    geom_smooth(show.legend=FALSE, method='loess', linewidth=0.5, fill="gray") +
    ggpubr::stat_cor(aes(x=age, y=beta), show.legend=FALSE, size=3.5, method="spearman")

  plot_down_AllSig.CpG=
    ggplot(LongDF_down_all, aes(x=age, y=beta)) +
    geom_point(alpha=1) +
    theme_classic() +
    labs(x="age", y="Beta-value") +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=10),
          title=element_text(size=10),
          legend.position="none") +
    labs(title=paste0(celltype)) +
    geom_smooth(show.legend=FALSE, method='loess', linewidth=0.5, fill="gray") +
    ggpubr::stat_cor(aes(x=age, y=beta), show.legend=FALSE, size=3.5, method="spearman")

  plot_AllSig.CpG=cowplot::plot_grid(plotlist=list(plot_up_AllSig.CpG, plot_down_AllSig.CpG), align="hv", ncol=2)

  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/",i,".pdf"), width=10, height=5)
  plot(plot_up_PerCpG)
  plot(plot_down_PerCpG)
  plot(plot_AllSig.CpG)
  dev.off()
}

### Merge the plots of all the celltypes
qpdf::pdf_combine(input=paste0("~/Project_PBMCage/Other_omics/Methy_plots/",c(1:length(LongDF_up_list)),".pdf"),
                  output="~/Project_PBMCage/Other_omics/Methy_plots/DMP_all_SigAssocWithAges.pdf")

#####################################



### Correlation analysis of gene-wise TSS200/TSS1500 DMPs with ages
#####################################
###

library(dplyr)
library(tidyverse)
library(ggplot2)

### Load the methylation data
Methylation_Data=readRDS("~/Project_PBMCage/Other_omics/AgeMeth_data.rds")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")

### Analyze the correlation and save the analysis results
Correlation_analysis=function(i) { # make the function for multicore
  celltype=names(Methylation_Data)[i]
  
  message(paste0(tissues_selected[i], " Starts..."))
  # load the DMP result
  results_per_celltype=load(paste0("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_",celltype,"_results.RData"))
  DMP_=myDMP[[1]]
  DMP_=subset(DMP_, feat.cgi %in% c("TSS200-island", "TSS1500-island"))
  DMP_$probeid=rownames(DMP_)
  gene_cpg_pair_df=data.frame(probeid=rownames(DMP_),
                              gene=DMP_$gene)
  
  # load and clean the raw data
  beta=Methylation_Data[[i]]$beta
  meta=Methylation_Data[[i]]$pd
  beta=beta[rownames(DMP_), ]
  beta_t=data.table::transpose(beta)
  colnames(beta_t)=rownames(beta)
  rownames(beta_t)=colnames(beta)
  beta_meta_merged=cbind(meta, beta_t[rownames(meta),])
  beta_meta_merged_longDF=tidyr::pivot_longer(beta_meta_merged, 
                                              cols=colnames(beta_meta_merged)[grepl("^cg[0-9]", colnames(beta_meta_merged))],
                                              names_to="probeid", 
                                              values_to="beta")
  beta_meta_merged_longDF=dplyr::left_join(beta_meta_merged_longDF, gene_cpg_pair_df, by="probeid")
  beta_meta_merged_longDF$age=as.numeric(as.character(beta_meta_merged_longDF$age))
  beta_meta_merged_longDF$beta=as.numeric(as.character(beta_meta_merged_longDF$beta))
  
  # take the mean of beta value for each gene
  beta_meta_merged_longDF_merged=beta_meta_merged_longDF %>% 
    subset((!is.na(beta)) & (!is.na(age))) %>%
    group_by(across(-c(probeid, beta))) %>% 
    summarize_at("beta", mean)
  saveRDS(beta_meta_merged_longDF_merged, paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_",i,".rds"))
  
  # remove the genes with <3 observations:
  genes_=beta_meta_merged_longDF_merged %>%
    group_by(sex, gene) %>%
    count() %>%
    subset(n>=3) %>% 
    ungroup() %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  genes_f=beta_meta_merged_longDF_merged %>%
    subset(sex=="F") %>%
    group_by(gene) %>%
    count() %>%
    subset(n>=3) %>% 
    ungroup() %>% 
    select(gene) %>% unlist() %>% unname()
  genes_m=beta_meta_merged_longDF_merged %>%
    subset(sex=="M") %>%
    group_by(gene) %>%
    count() %>%
    subset(n>=3) %>% 
    select(gene) %>% unlist() %>% unname()

  # analyze the correlation of all/female/male samples
  # ..get all the expr data
  pb=txtProgressBar(min=0, max=length(genes_), style=3, width=60, char="-")
  all_cor=all_cor_pval=all_rho=all_rho_pval=all_tau=all_tau_pval=list()
  for (j in 1:length(genes_)) {
    subset_=subset(beta_meta_merged_longDF_merged, gene==genes_[j])
    cor_result=cor.test(subset_$age, subset_$beta, method="pearson")
    rho_result=cor.test(subset_$age, subset_$beta, method="spearman", exact=F)
    tau_result=cor.test(subset_$age, subset_$beta, method="kendall")
    all_cor[[j]]=cor_result$estimate; all_cor_pval[[j]]=cor_result$p.value
    all_rho[[j]]=rho_result$estimate; all_rho_pval[[j]]=rho_result$p.value
    all_tau[[j]]=tau_result$estimate; all_tau_pval[[j]]=tau_result$p.value
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # ..get the females' expr
  pb=txtProgressBar(min=0, max=length(genes_f), style=3, width=60, char="=")
  female_cor=female_cor_pval=female_rho=female_rho_pval=female_tau=female_tau_pval=list()
  for (j in 1:length(genes_f)) {
    subset_=subset(beta_meta_merged_longDF_merged, sex=="F" & gene==genes_f[j])
    cor_result=cor.test(subset_$age, subset_$beta, method="pearson")
    rho_result=cor.test(subset_$age, subset_$beta, method="spearman", exact=F)
    tau_result=cor.test(subset_$age, subset_$beta, method="kendall")
    female_cor[[j]]=cor_result$estimate; female_cor_pval[[j]]=cor_result$p.value
    female_rho[[j]]=rho_result$estimate; female_rho_pval[[j]]=rho_result$p.value
    female_tau[[j]]=tau_result$estimate; female_tau_pval[[j]]=tau_result$p.value
    setTxtProgressBar(pb, j)
  }
  close(pb)
  # ..get the males' expr
  pb=txtProgressBar(min=0, max=length(genes_m), style=3, width=60, char="+")
  male_cor=male_cor_pval=male_rho=male_rho_pval=male_tau=male_tau_pval=list()
  for (j in 1:length(genes_m)) {
    subset_=subset(beta_meta_merged_longDF_merged, sex=="M" & gene==genes_m[j])
    cor_result=cor.test(subset_$age, subset_$beta, method="pearson")
    rho_result=cor.test(subset_$age, subset_$beta, method="spearman", exact=F)
    tau_result=cor.test(subset_$age, subset_$beta, method="kendall")
    male_cor[[j]]=cor_result$estimate; male_cor_pval[[j]]=cor_result$p.value
    male_rho[[j]]=rho_result$estimate; male_rho_pval[[j]]=rho_result$p.value
    male_tau[[j]]=tau_result$estimate; male_tau_pval[[j]]=tau_result$p.value
    setTxtProgressBar(pb, j)
  }
  close(pb)
  
  # combine the correlation analysis results in dataframe
  ALL_COR_DF=Female_COR_DF=Male_COR_DF=data.frame()
  if (length(genes_)>=1) {
    ALL_COR_DF=data.frame(gene=genes_, 
                          cor=unlist(all_cor), cor_pval=unlist(all_cor_pval), rho=unlist(all_rho), rho_pval=unlist(all_rho_pval), tau=unlist(all_tau), tau_pval=unlist(all_tau_pval),
                          analysis="all", celltypes=tissues_selected[i])
  }
  if (length(genes_f)>=1) {
    Female_COR_DF=data.frame(gene=genes_f, 
                             cor=unlist(female_cor), cor_pval=unlist(female_cor_pval), rho=unlist(female_rho), rho_pval=unlist(female_rho_pval), tau=unlist(female_tau), tau_pval=unlist(female_tau_pval),
                             analysis="females", celltypes=tissues_selected[i])
  }
  if (length(genes_m)>=1) {
    Male_COR_DF=data.frame(gene=genes_m, 
                           cor=unlist(male_cor), cor_pval=unlist(male_cor_pval), rho=unlist(male_rho), rho_pval=unlist(male_rho_pval), tau=unlist(male_tau), tau_pval=unlist(male_tau_pval),
                           analysis="males", celltypes=tissues_selected[i])
  }
  correlation_DF=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))
  
  message(paste0("All the genes in ",tissues_selected[i]," have been gone through!"))
  saveRDS(correlation_DF, paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult_",i,".rds"))
}
# run with multicores
Correlation_analysis(1)
Correlation_analysis(2)
Correlation_analysis(3)
Correlation_analysis(4)
Correlation_analysis(5)
Correlation_analysis(6)
Correlation_analysis(7)

### Save the gene-wise beta value from 7 celltypes into the same file
Data_list=list()
for (i in 1:7) {
  Data_list[[i]]=readRDS(paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_",i,".rds"))
  Data_list[[i]]=Data_list[[i]] %>% ungroup() %>% select(Sample_Name, age, project_id, sex, gene, beta)
  Data_list[[i]]$celltypes=tissues_selected[i]
}
Data=data.table::rbindlist(Data_list)
data.table::fwrite(Data, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")
Data=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t") # check
unlink(paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_",c(1:7),".rds")) # delete the seven files

### Save the correlation_DFs from 7 celltypes into the same file
correlation_DF_list=list()
for (i in 1:7) {
  correlation_DF_list[[i]]=readRDS(paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult_",i,".rds"))
}
correlation_DF=data.table::rbindlist(correlation_DF_list)
data.table::fwrite(correlation_DF, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
correlation_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t") # check
unlink(paste0("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult_",c(1:7),".rds")) # delete the seven files

#####################################



### Upset analysis of the age-correlated gene-wise TSS DMPs
#####################################
###

library(dplyr)

Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
all_celltypes=names(table(Cor_analysis_DF$celltypes))
upset_plot=list()
for (i in 1:length(all_celltypes)) {
  Subset_=subset(Cor_analysis_DF, celltypes==all_celltypes[i])
  all_up=Subset_ %>%
    subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
    select(gene) %>%
    unlist() %>% unname()
  all_down=Subset_ %>%
    subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
    select(gene) %>%
    unlist() %>% unname()
  female_up=Subset_ %>%
    subset(analysis=="females" & rho_pval<0.05 & rho>0) %>%
    select(gene) %>%
    unlist() %>% unname()
  female_down=Subset_ %>%
    subset(analysis=="females" & rho_pval<0.05 & rho<0) %>%
    select(gene) %>%
    unlist() %>% unname()
  male_up=Subset_ %>%
    subset(analysis=="males" & rho_pval<0.05 & rho>0) %>%
    select(gene) %>%
    unlist() %>% unname()
  male_down=Subset_ %>%
    subset(analysis=="males" & rho_pval<0.05 & rho<0) %>%
    select(gene) %>%
    unlist() %>% unname()

  # make the UpSet data
  library(ComplexHeatmap)
  library(ggplot2)
  datalist=list(Hypo.both=all_down, Hypo.F=female_down, Hypo.M=male_down, Hyper.both=all_up, Hyper.F=female_up, Hyper.M=male_up)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[1:3]/set_size(m)[1],
              comb_size(m)[4:6]/set_size(m)[4])
  percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
  title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated gene-wise TSS DMPs in\n",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
  
  # plot
  upset_plot_=
    HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
    UpSet(m,
          top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
          right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
          row_names_gp=gpar(fontsize=9),
          set_order=c("Hypo.both","Hypo.F","Hypo.M","Hyper.both","Hyper.F","Hyper.M"),
          comb_order=1:6) %v%
    HeatmapAnnotation(
      empty=anno_empty(border=F, height=unit(2,"mm")),
      text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
    )
  upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
pdf("~/Project_PBMCage/Other_omics/Methy_plots/TSS.DMPs_Correlation_wAges_UpSet.pdf", width=5, height=4)
for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
dev.off()

#####################################



### Plot the correlation between the top30 age-correlated gene-wise TSS DMPs and ages
#####################################
###

library(dplyr)
library(tidyverse)
library(ggplot2)
library(Seurat)

### Filter the correlation results with p<0.05 and analysis=="females"|"males"
correlation_DF_ALL=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")

### Load the beta value data and the metadata
beta_meta_merged_ALL=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")

### Plot the top age-correlated TSS DMPs with ages
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
pos_plot=neg_plot=list()
pb=txtProgressBar(min=0, max=7, style=3, width=60, char="=")
for (i in 1:length(tissues_selected)) {

    # extract the top30 correlated genes
  beta_meta_merged_longDF=beta_meta_merged_ALL %>% subset(sex %in% c("M","F") & celltypes==tissues_selected[i])
  correlation_DF=subset(correlation_DF_ALL, celltypes==tissues_selected[i])
  
  Cor_df_neg=correlation_DF %>%
    subset(rho_pval<0.05 & rho<0) %>%
    slice_min(rho, n=200) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
  
  Cor_df_pos=Cor_analysis_DF %>%
    subset(rho_pval<0.05 & rho>0) %>%
    slice_max(rho, n=200) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]

  # make the df for plotting
  Top_betaDF_neg=beta_meta_merged_longDF %>%
    subset(gene %in% Cor_df_neg)
  Top_betaDF_pos=beta_meta_merged_longDF %>%
    subset(gene %in% Cor_df_pos)

  # plot the hypomethylated genes
  neg_plot[[i]]=
    ggplot(Top_betaDF_neg, aes(x=age, y=beta, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Beta-value", title=paste0("Hypomethylated genes in ", tissues_selected[i])) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          # legend.position="none",
          # legend.title=element_blank(),
          title=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
    ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

  # plot the hypermethylated genes
  pos_plot[[i]]=
    ggplot(Top_betaDF_pos, aes(x=age, y=beta, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Beta-value", title=paste0("Hypermethylated genes in ", tissues_selected[i])) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          # legend.position="none",
          # legend.title=element_blank(),
          title=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
    ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")

  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
  plot(neg_plot[[i]])
  dev.off()

  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
  plot(pos_plot[[i]])
  dev.off()

  setTxtProgressBar(pb, i)
}
close(pb)

### Combine the pdf files
qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",c(1:7),"_neg.pdf"),
                          paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",c(1:7),"_pos.pdf")),
                  output="~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots.pdf")
unlink(c(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",c(1:7),"_neg.pdf"),
         paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_CorPlots_",c(1:7),"_pos.pdf")))

#####################################



### Analyze the consistency between the age-correlated gene-wise TSS DMPs and transcripts
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

### Load the TSS-age correlation and the RNA-age correlation
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")

RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))])

### Extract the consistent genes in each celltype
All.hyper.neg_List=F.hyper.neg_List=M.hyper.neg_List=All.hypo.pos_List=F.hypo.pos_List=M.hypo.pos_List=list()
for (i in 1:length(tissues_selected)) {
  # for the genes at hypermethylation and downregulation in all samples
  TSS_subset_all.hyper=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho>0 & analysis=="all")
  RNA_subset_all.neg=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho<0 & analysis=="all")
  shared_subset_all.hyper.neg=dplyr::inner_join(TSS_subset_all.hyper, RNA_subset_all.neg, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_all.hyper.neg=shared_subset_all.hyper.neg %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_all.hyper.neg$analysis.orientation="Hyper_Neg"
  All.hyper.neg_List[[i]]=shared_subset_all.hyper.neg
  # for the genes at hypermethylation and downregulation in females
  TSS_subset_f.hyper=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho>0 & analysis=="females")
  RNA_subset_f.neg=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho<0 & analysis=="females")
  shared_subset_f.hyper.neg=dplyr::inner_join(TSS_subset_f.hyper, RNA_subset_f.neg, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_f.hyper.neg=shared_subset_f.hyper.neg %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_f.hyper.neg$analysis.orientation="Hyper_Neg"
  F.hyper.neg_List[[i]]=shared_subset_f.hyper.neg
  # for the genes at hypermethylation and downregulation in all males
  TSS_subset_m.hyper=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho>0 & analysis=="males")
  RNA_subset_m.neg=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho<0 & analysis=="males")
  shared_subset_m.hyper.neg=dplyr::inner_join(TSS_subset_m.hyper, RNA_subset_m.neg, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_m.hyper.neg=shared_subset_m.hyper.neg %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_m.hyper.neg$analysis.orientation="Hyper_Neg"
  M.hyper.neg_List[[i]]=shared_subset_m.hyper.neg

  # for the genes at hypomethylation and upregulation in all samples
  TSS_subset_all.hypo=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho<0 & analysis=="all")
  RNA_subset_all.pos=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho>0 & analysis=="all")
  shared_subset_all.hypo.pos=dplyr::inner_join(TSS_subset_all.hypo, RNA_subset_all.pos, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_all.hypo.pos=shared_subset_all.hypo.pos %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_all.hypo.pos$analysis.orientation="Hypo_Pos"
  All.hypo.pos_List[[i]]=shared_subset_all.hypo.pos
  # for the genes at hypomethylation and upregulation in females
  TSS_subset_f.hypo=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho<0 & analysis=="females")
  RNA_subset_f.pos=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho>0 & analysis=="females")
  shared_subset_f.hypo.pos=dplyr::inner_join(TSS_subset_f.hypo, RNA_subset_f.pos, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_f.hypo.pos=shared_subset_f.hypo.pos %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_f.hypo.pos$analysis.orientation="Hypo_Pos"
  F.hypo.pos_List[[i]]=shared_subset_f.hypo.pos
  # for the genes at hypomethylation and upregulation in males
  TSS_subset_m.hypo=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05) %>% subset(rho<0 & analysis=="males")
  RNA_subset_m.pos=RNA_Cor %>% subset(celltypes %in% shared_celltype[[i]] & rho_pval<0.05) %>% subset(rho>0 & analysis=="males")
  shared_subset_m.hypo.pos=dplyr::inner_join(TSS_subset_m.hypo, RNA_subset_m.pos, by=c("gene", "analysis"), suffix=c(".TSS",".RNA"))
  shared_subset_m.hypo.pos=shared_subset_m.hypo.pos %>%
    mutate(Cor_Consis=sqrt(abs(cor.TSS*cor.RNA)))
  shared_subset_m.hypo.pos$analysis.orientation="Hypo_Pos"
  M.hypo.pos_List[[i]]=shared_subset_m.hypo.pos
}

### Save the consistent correlation table
Cor_Consistency_DF=data.table::rbindlist(c(All.hyper.neg_List, F.hyper.neg_List, M.hyper.neg_List, All.hypo.pos_List, F.hypo.pos_List, M.hypo.pos_List))
data.table::fwrite(Cor_Consistency_DF, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")

#####################################



### Plot the correlation consistency volcano plots
#####################################
###

###
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
Cor_Consistency_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")
plot_hyper.neg_all=plot_hypo.pos_all=list()
for (i in 1:length(tissues_selected)) {
  # hypermethylated-downregulated
  shared_subset_hyper.neg=Cor_Consistency_DF %>% subset(analysis.orientation=="Hyper_Neg" & celltypes.TSS==tissues_selected[i]) %>%
    mutate(analysis=ifelse(analysis=="all","<both>",
                           ifelse(analysis=="females","<females>","<males>")),
           pval.TSS=rho_pval.TSS)
  plot_hyper.neg_all[[i]]=
    ggplot(shared_subset_hyper.neg, aes(x=Cor_Consis, y=-log10(rho_pval.RNA), color=pval.TSS)) +
    ggrastr::geom_point_rast() +
    labs(title=paste0("Hypermethylated genes in ",tissues_selected[i]), y=expression(-log[10]~pval.RNA), x="Correlation") +
    facet_wrap(~analysis) +
    theme_classic() +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
          strip.text.x=element_text(size=10, color="black", face="italic")) +
    scale_colour_gradientn(colours=c("brown4", "chocolate1")) +
    geom_text_repel(aes(label=gene), size=3, show.legend=FALSE)

  # hypomethylated-upregulated
  shared_subset_hypo.pos=Cor_Consistency_DF %>% subset(analysis.orientation=="Hypo_Pos" & celltypes.TSS==tissues_selected[i]) %>%
    mutate(analysis=ifelse(analysis=="all","<both>",
                           ifelse(analysis=="females","<females>","<males>")),
           pval.TSS=rho_pval.TSS)
  plot_hypo.pos_all[[i]]=
    ggplot(shared_subset_hypo.pos, aes(x=Cor_Consis, y=-log10(rho_pval.RNA), color=pval.TSS)) +
    ggrastr::geom_point_rast() +
    labs(title=paste0("Hypomethylated genes in ",tissues_selected[i]), y=expression(-log[10]~pval.RNA), x="Correlation") +
    facet_wrap(~analysis) +
    theme_classic() +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
          strip.text.x=element_text(size=10, color="black", face="italic")) +
    scale_colour_gradientn(colours=c("brown4", "chocolate1")) +
    geom_text_repel(aes(label=gene), size=3, show.legend=FALSE)
}
pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency.pdf", height=5, width=20)
for (i in 1:length(tissues_selected)) {plot(plot_hyper.neg_all[[i]]); plot(plot_hypo.pos_all[[i]])}
dev.off()

### Plot the consistent genes in the mixed sex with volcano plots
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
Cor_Consistency_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")

shared_subset=Cor_Consistency_DF %>% subset(analysis=="all" & celltypes.TSS %in% c("CD14+ monocyte","CD4+ T cell","CD8+ T cell")) %>%
  mutate(Correlation=ifelse(analysis.orientation=="Hyper_Neg", -Cor_Consis, Cor_Consis),
         pval.TSS=rho_pval.TSS)

plots_=
  ggplot(shared_subset, aes(x=Correlation, y=-log10(rho_pval.RNA), alpha=pval.TSS, color=celltypes.TSS)) +
  facet_wrap(~celltypes.TSS, nrow=1) +
  ggrastr::geom_point_rast() +
  labs(y=expression(-log[10]~pval.RNA), x="Correlation") +
  theme_classic() +
  theme(title=element_text(size=11),
        plot.title=element_text(hjust=0.5),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        axis.title.x=element_text(size=11),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10),
        strip.text=element_text(size=11)) +
  scale_alpha_continuous(range=c(1,0.1)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(5,2,3)]) +
  guides(color="none") +
  scale_x_continuous(labels=abs, limits=c(-0.6,0.6)) +
  ggrepel::geom_text_repel(aes(label=gene), size=3, show.legend=FALSE)

pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_VolcanoPlot.pdf", height=4.5, width=15)
plot(plots_)
dev.off()

#####################################



### Plot the RNA/TSS age-correlation consistency
#####################################
###

# ### Load the TSS-age correlation and the RNA-age correlation
# TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
# TSS_Cor=TSS_Cor %>%
#   mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
#                           ifelse(celltypes=="CD8+ T cell","CD8T cells",
#                                  ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
#   # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
#   subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
#   subset(analysis=="all") %>%
#   dplyr::select(gene, rho, rho_pval, celltypes)
# 
# RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# RNA_Cor=RNA_Cor %>%
#   # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
#   subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
#   subset(analysis=="all") %>%
#   dplyr::select(gene, rho, rho_pval, celltypes)
# 
# Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes"), suffix=c(".TSS",".RNA"))

### Load the TSS-age correlation and the RNA-age correlation, for F and M respectively
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
TSS_Cor=TSS_Cor %>%
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)

RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
RNA_Cor=RNA_Cor %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)

Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes","analysis"), suffix=c(".TSS",".RNA")) %>%
  mutate(analysis=ifelse(analysis=="females","F","M"))

### Plot
plot_=
  ggplot(Cor_Consistency_df, aes(x=rho.TSS, y=rho.RNA)) +
  facet_wrap(~celltypes) +
  ggrastr::geom_point_rast(size=0.2, alpha=0.2, shape=20, color="grey90") +
  geom_smooth(aes(group=analysis, color=analysis), method="lm", linewidth=0.5, se=F) +
  labs(y=expression(atop(Spearman~rho~"in",RNA~transcription)),
       x=expression(Spearman~rho~"in"~TSS~methylation)) +
  # labs(y=expression(Spearman~rho~"in"~RNA~transcription),
  #      x=expression(Spearman~rho~"in"~TSS~methylation)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="right",
        legend.direction="vertical") +
  ggpubr::stat_cor(aes(group=analysis, color=analysis), size=3.5, show.legend=F) +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
    
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_corrplot.pdf", height=2, width=4)
plot_
dev.off()

# ### Analyze in each supercluster of the 12 key terms
# TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
# TSS_Cor=TSS_Cor %>%
#   mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
#                           ifelse(celltypes=="CD8+ T cell","CD8T cells",
#                                  ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
#   # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
#   subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
#   subset(analysis!="all") %>%
#   dplyr::select(gene, rho, rho_pval, celltypes, analysis)
# RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# RNA_Cor=RNA_Cor %>%
#   # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
#   subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
#   subset(analysis!="all") %>%
#   dplyr::select(gene, rho, rho_pval, celltypes, analysis)
# Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes","analysis"), suffix=c(".TSS",".RNA")) %>%
#   mutate(analysis=ifelse(analysis=="females","F","M"))
# 
# Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
# subset_TSS_cordf=lapply(1:length(Genes_in_SuperClusters), 
#                         function(idx) Cor_Consistency_df %>% 
#                           subset(gene %in% Genes_in_SuperClusters[[idx]]) %>% 
#                           # subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05) %>%
#                           mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters))[idx])) %>%
#   data.table::rbindlist()
# # reorder the stripts
# columns_reorder=c("nucl. metab.","prot. metab.","ribosome syn.","cellular org.","energy metab.","cell division","autophagy","prog. death",
#                   "cell response","anti-virus","leuk. devel.","immune proc.")
# subset_TSS_cordf$superclusters=forcats::fct_relevel(subset_TSS_cordf$superclusters, columns_reorder)
# 
# plot_supercluster=
#   ggplot(subset_TSS_cordf, aes(x=rho.TSS, y=rho.RNA, color=celltypes)) +
#   ggrastr::geom_point_rast(alpha=0.01, size=0.1, shape=1) +
#   scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(2,3)]) +
#   facet_wrap(~superclusters) +
#   geom_smooth(data=. %>% subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05),
#               aes(group=celltypes), method="lm", fill="grey93", alpha=0.5, linewidth=0.5) +
#   ggpubr::stat_cor(data=. %>% subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05), 
#                    aes(group=celltypes), show.legend=F, label.y=c(-0.3,-0.5),
#                    size=3) +
#   theme_classic() +
#   labs(y=expression(Spearman~rho~"in"~RNA~transcription), x=expression(Spearman~rho~"in"~TSS~methylation)) +
#   theme(axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         strip.text=element_text(size=10),
#         legend.text=element_text(size=9),
#         title=element_text(size=11),
#         legend.position="top",
#         legend.direction="horizontal") +
#   guides(color=guide_legend(title=NULL, override.aes=list(size=2, shape=19, alpha=1, fill="transparent")))
# 
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_MapToSuperclusters.pdf", height=5.5, width=6.5)
# plot_supercluster
# dev.off()

#####################################



### Plot the heatmap with pearson correlation between rho.TSS and rho.RNA in each supercluster in CD4T and CD8T
#####################################
###

library(ComplexHeatmap)

### PBMCage dataset 
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
TSS_Cor=TSS_Cor %>%
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
RNA_Cor=RNA_Cor %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)
Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes","analysis"), suffix=c(".TSS",".RNA")) %>%
  mutate(analysis=ifelse(analysis=="females","F","M"))

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
subset_TSS_cordf=lapply(1:length(Genes_in_SuperClusters), 
                        function(idx) Cor_Consistency_df %>% 
                          subset(gene %in% Genes_in_SuperClusters[[idx]]) %>% 
                          # subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05) %>%
                          mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters))[idx])) %>%
  data.table::rbindlist()
# calculate the pearson correlation between rho.TSS and rho.RNA
subset_TSS_cordf_pearson=subset_TSS_cordf %>% 
  subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05) %>%
  split(.$superclusters) %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% lapply(., function(dff) cor.test(dff$rho.TSS, dff$rho.RNA, method="pearson")))
subset_TSS_cordf_pearson_cor=subset_TSS_cordf_pearson %>%
  lapply(., function(item) {lapply(item, function(item_sub) item_sub$estimate)}) %>%
  lapply(., function(x) unlist(x)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("celltype") %>%
  mutate(celltype=gsub("\\.cor","",celltype)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("celltype")
colnames(subset_TSS_cordf_pearson_cor)=names(subset_TSS_cordf_pearson)
subset_TSS_cordf_pearson_pval=subset_TSS_cordf_pearson %>%
  lapply(., function(item) {lapply(item, function(item_sub) item_sub$p.value)}) %>%
  lapply(., function(x) unlist(x)) %>%
  as.data.frame()
colnames(subset_TSS_cordf_pearson_pval)=names(subset_TSS_cordf_pearson)
subset_TSS_cordf_pearson_pval=
  apply(subset_TSS_cordf_pearson_pval, c(1,2), function(x)
        ifelse(x<0.0001,"****",
               ifelse(x>=0.0001 & x<0.001, "***",
                      ifelse(x>=0.001 & x<0.01, "**",
                             ifelse(x>=0.01 & x<0.05, "*", "")))))
# plot
range(subset_TSS_cordf_pearson_cor) # check
col_fun_sub=circlize::colorRamp2(c(-0.5, 0, 0.5), c("dodgerblue3", "white", "brown3"))
# ht_=
  Heatmap(subset_TSS_cordf_pearson_cor,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title=expression(Spearman~rho),
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
          layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
            v=pindex(subset_TSS_cordf_pearson_pval, i, j)
            grid.text(v, x, y, rot=0, just="center",
                      gp=gpar(fontsize=10, fontface="bold"))
          },
          width=ncol(subset_TSS_cordf_pearson_cor)*unit(7, "mm"),
          height=nrow(subset_TSS_cordf_pearson_cor)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

### Terekhova dataset 
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
TSS_Cor=TSS_Cor %>%
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
RNA_Cor=RNA_Cor %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)
Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes","analysis"), suffix=c(".TSS",".RNA")) %>%
  mutate(analysis=ifelse(analysis=="females","F","M"))

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
subset_TSS_cordf=lapply(1:length(Genes_in_SuperClusters), 
                        function(idx) Cor_Consistency_df %>% 
                          subset(gene %in% Genes_in_SuperClusters[[idx]]) %>% 
                          # subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05) %>%
                          mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters))[idx])) %>%
  data.table::rbindlist()
# calculate the pearson correlation between rho.TSS and rho.RNA
subset_TSS_cordf_pearson=subset_TSS_cordf %>% 
  subset(rho_pval.RNA<0.05 & rho_pval.TSS<0.05) %>%
  split(.$superclusters) %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% lapply(., function(dff) cor.test(dff$rho.TSS, dff$rho.RNA, method="pearson")))
subset_TSS_cordf_pearson_cor=subset_TSS_cordf_pearson %>%
  lapply(., function(item) {lapply(item, function(item_sub) item_sub$estimate)}) %>%
  lapply(., function(x) unlist(x)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("celltype") %>%
  mutate(celltype=gsub("\\.cor","",celltype)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("celltype")
colnames(subset_TSS_cordf_pearson_cor)=names(subset_TSS_cordf_pearson)
subset_TSS_cordf_pearson_pval=subset_TSS_cordf_pearson %>%
  lapply(., function(item) {lapply(item, function(item_sub) item_sub$p.value)}) %>%
  lapply(., function(x) unlist(x)) %>%
  as.data.frame()
colnames(subset_TSS_cordf_pearson_pval)=names(subset_TSS_cordf_pearson)
subset_TSS_cordf_pearson_pval=
  apply(subset_TSS_cordf_pearson_pval, c(1,2), function(x)
    ifelse(x<0.0001,"****",
           ifelse(x>=0.0001 & x<0.001, "***",
                  ifelse(x>=0.001 & x<0.01, "**",
                         ifelse(x>=0.01 & x<0.05, "*", "")))))
# plot
range(subset_TSS_cordf_pearson_cor) # check
col_fun_sub=circlize::colorRamp2(c(-0.5, 0, 0.5), c("dodgerblue3", "white", "brown3"))
# ht_=
Heatmap(subset_TSS_cordf_pearson_cor,
        name="mat", show_heatmap_legend=T,
        heatmap_legend_param=list(
          title=expression(Spearman~rho),
          legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
          direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
        ),
        rect_gp=gpar(col="white", lwd=0.5),
        col=col_fun_sub,
        cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
        column_names_rot=90,
        show_column_names=T,
        column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
        column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
        layer_fun=function(j, i, x, y, width, height, fill, slice_r, slice_c) {
          v=pindex(subset_TSS_cordf_pearson_pval, i, j)
          grid.text(v, x, y, rot=0, just="center",
                    gp=gpar(fontsize=10, fontface="bold"))
        },
        width=ncol(subset_TSS_cordf_pearson_cor)*unit(7, "mm"),
        height=nrow(subset_TSS_cordf_pearson_cor)*unit(5.5, "mm")
)
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

#####################################



### Plot the heatmap with pearson correlation between TSS and RNA in each supercluster in CD8T
#####################################
###

library(ComplexHeatmap)
library(Seurat)

### TSS
Data=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t") # check
Data_subset=Data %>% subset(celltypes=="CD8+ T cell")
Data_subset=Data_subset %>% 
  mutate(age=round(age)) %>%
  group_by(age, gene) %>% summarize_at("beta", mean)

### RNA
pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj=subset(pseudobulkobj, Annot.rough=="CD8T cells")
RNA_data=LayerData(pseudobulkobj, layer="data") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("gene",colnames(.))], names_to="group", values_to="data") %>%
  mutate(donor_id=gsub("CD8T cells_|_[0-9]+$","",group),
         age=as.integer(gsub(".*-[0-9]+_","",group))) %>%
  dplyr::select(gene, age, data) %>%
  group_by(gene, age) %>%
  summarize_at("data", mean) %>%
  rename(RNA.expr=data)

### Merge
Meth_RNA=inner_join(Data_subset, RNA_data, by=c("age","gene"))
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
Meth_RNA_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) Meth_RNA %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx]))) %>%
  data.table::rbindlist(.)

### Calculate pearson correlation
Meth_RNA_percluster_cortest=Meth_RNA_percluster %>% 
  arrange(age) %>% 
  split(.$superclusters) %>%
  lapply(. , function(df) {
    test=cor.test(df$beta, df$RNA.expr)
    data.frame(cor=test$estimate, pval=test$p.value)
  }) %>%
  unlist()
Meth_RNA_percluster_df=data.frame(group=names(Meth_RNA_percluster_cortest),
                                  value=Meth_RNA_percluster_cortest) %>%
  mutate(supercluster=gsub("\\.cor|\\.pval","",group)) %>%
  mutate(value.type=ifelse(grepl("cor$",group),"cor","pval")) %>%
  tibble::remove_rownames() %>%
  dplyr::select(-group) %>%
  dplyr::relocate(c("supercluster","value.type","value"))
# cor
Meth_RNA_percluster_cor.pval=Meth_RNA_percluster_df %>%
  tidyr::pivot_wider(names_from="value.type", values_from="value") %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("supercluster") %>%
  mutate(pval=-log10(pval))
supecluster_=names(table(rownames(Meth_RNA_percluster_cor.pval)))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
Meth_RNA_percluster_cor.pval=Meth_RNA_percluster_cor.pval[supecluster_[supecluster_reorder_idx],]

range(Meth_RNA_percluster_cor.pval[,1]) # check
col_fun_cor=circlize::colorRamp2(c(-0.15, 0), c("dodgerblue3", "white"))
ht_cor=
Heatmap(Meth_RNA_percluster_cor.pval %>% dplyr::select(cor),
        name="Pearson r", show_heatmap_legend=T,
        heatmap_legend_param=list(
          title="Pearson r", title_gp=gpar(fontface="plain"),
          legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
          direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
        ),
        rect_gp=gpar(col="white", lwd=0.5),
        border_gp=gpar(col="black", lty=1, lwd=0.25),
        col=col_fun_cor,
        cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
        column_names_rot=90,
        show_column_names=F,
        show_row_names=F,
        column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
        column_title=NULL, column_title_gp=gpar(fontsize=11),
        width=1*unit(5, "mm"),
        height=nrow(Meth_RNA_percluster_cor.pval)*unit(8.5, "mm")
)

range(Meth_RNA_percluster_cor.pval[,2]) # check
col_fun_pval=circlize::colorRamp2(c(0, 150), c("white","brown3"))
ht_pval=
Heatmap(Meth_RNA_percluster_cor.pval %>% dplyr::select(pval),
        name="p-value", show_heatmap_legend=T,
        heatmap_legend_param=list(
          title="p-value", title_gp=gpar(fontface="plain"),
          legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
          direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
        ),
        rect_gp=gpar(col="white", lwd=0.5),
        border_gp=gpar(col="black", lty=1, lwd=0.25),
        col=col_fun_pval,
        cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
        column_names_rot=90,
        show_column_names=F,
        show_row_names=T,
        column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
        column_title=NULL, column_title_gp=gpar(fontsize=11),
        width=1*unit(5, "mm"),
        height=nrow(Meth_RNA_percluster_cor.pval)*unit(8.5, "mm")
)

ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_cor + ht_pval, heatmap_legend_side="left", merge_legend=TRUE)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_MapToSuperclusters_heatmap.pdf", height=4.5, width=2.5)
ht_final
dev.off()

#####################################



### Plot the RNA/TSS age-correlation consistency with Immunity dataset
#####################################
###

### Load the TSS-age correlation and the RNA-age correlation, for F and M respectively
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
TSS_Cor=TSS_Cor %>%
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)

RNA_Cor=data.table::fread("~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
RNA_Cor=RNA_Cor %>%
  # subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(celltypes %in% c("CD8T cells")) %>%
  subset(analysis!="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes, analysis)

Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes","analysis"), suffix=c(".TSS",".RNA")) %>%
  mutate(analysis=ifelse(analysis=="females","F","M"))

### Plot
plot_=
  ggplot(Cor_Consistency_df, aes(x=rho.TSS, y=rho.RNA)) +
  facet_wrap(~celltypes) +
  ggrastr::geom_point_rast(size=0.2, alpha=0.2, shape=20, color="grey90") +
  geom_smooth(aes(group=analysis, color=analysis), method="lm", linewidth=0.5, se=F) +
  labs(y=expression(atop(Spearman~rho~"in",RNA~transcription)),
       x=expression(Spearman~rho~"in"~TSS~methylation)) +
  # labs(y=expression(Spearman~rho~"in"~RNA~transcription),
  #      x=expression(Spearman~rho~"in"~TSS~methylation)) +
  theme_classic() +
  theme(plot.background=element_rect(fill="transparent", color="transparent"),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="right",
        legend.direction="vertical") +
  ggpubr::stat_cor(aes(group=analysis, color=analysis), size=3.5, show.legend=F, label.y=c(0.4,-0.3)) +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_corrplot.pdf", height=1.5, width=3.5)
plot_
dev.off()

#####################################



### Plot the heatmap with pearson correlation between TSS and RNA in each supercluster in CD8T
#####################################
###

library(ComplexHeatmap)
library(Seurat)

### TSS
Data=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t") # check
Data_subset=Data %>% subset(celltypes=="CD8+ T cell")
Data_subset=Data_subset %>% 
  mutate(age=round(age)) %>%
  group_by(age, gene) %>% summarize_at("beta", mean)

### RNA
pseudobulkobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
pseudobulkobj=subset(pseudobulkobj, Annot.rough=="CD8T cells")
RNA_data=LayerData(pseudobulkobj, layer="data") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("gene",colnames(.))], names_to="group", values_to="data") %>%
  mutate(donor_id=gsub("CD8T cells_|_[0-9]+$","",group),
         age=as.integer(gsub(".*_[A-Z]+[0-9]+_","",group))) %>%
  dplyr::select(gene, age, data) %>%
  group_by(gene, age) %>%
  summarize_at("data", mean) %>%
  rename(RNA.expr=data)

### Merge
Meth_RNA=inner_join(Data_subset, RNA_data, by=c("age","gene"))
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
Meth_RNA_percluster=lapply(1:length(Genes_in_SuperClusters), 
                           function(idx) Meth_RNA %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                             mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx]))) %>%
  data.table::rbindlist(.)

### Calculate pearson correlation
Meth_RNA_percluster_cortest=Meth_RNA_percluster %>% 
  arrange(age) %>% 
  split(.$superclusters) %>%
  lapply(. , function(df) {
    test=cor.test(df$beta, df$RNA.expr)
    data.frame(cor=test$estimate, pval=test$p.value)
  }) %>%
  unlist()
Meth_RNA_percluster_df=data.frame(group=names(Meth_RNA_percluster_cortest),
                                  value=Meth_RNA_percluster_cortest) %>%
  mutate(supercluster=gsub("\\.cor|\\.pval","",group)) %>%
  mutate(value.type=ifelse(grepl("cor$",group),"cor","pval")) %>%
  tibble::remove_rownames() %>%
  dplyr::select(-group) %>%
  dplyr::relocate(c("supercluster","value.type","value"))
# cor
Meth_RNA_percluster_cor.pval=Meth_RNA_percluster_df %>%
  tidyr::pivot_wider(names_from="value.type", values_from="value") %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("supercluster") %>%
  mutate(pval=-log10(pval))
supecluster_=names(table(rownames(Meth_RNA_percluster_cor.pval)))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
Meth_RNA_percluster_cor.pval=Meth_RNA_percluster_cor.pval[supecluster_[supecluster_reorder_idx],]

range(Meth_RNA_percluster_cor.pval[,1]) # check
col_fun_cor=circlize::colorRamp2(c(-0.15, 0), c("dodgerblue3", "white"))
ht_cor=
  Heatmap(Meth_RNA_percluster_cor.pval %>% dplyr::select(cor),
          name="Pearson r", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Pearson r", title_gp=gpar(fontface="plain"),
            legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=1, lwd=0.25),
          col=col_fun_cor,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          show_row_names=F,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title=NULL, column_title_gp=gpar(fontsize=11),
          width=1*unit(5, "mm"),
          height=nrow(Meth_RNA_percluster_cor.pval)*unit(5, "mm")
  )

range(Meth_RNA_percluster_cor.pval[,2]) # check
col_fun_pval=circlize::colorRamp2(c(0, 150), c("white","brown3"))
ht_pval=
  Heatmap(Meth_RNA_percluster_cor.pval %>% dplyr::select(pval),
          name="p-value", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="p-value", title_gp=gpar(fontface="plain"),
            legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9)
          ),
          rect_gp=gpar(col="white", lwd=0.5),
          border_gp=gpar(col="black", lty=1, lwd=0.25),
          col=col_fun_pval,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=F,
          show_row_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title=NULL, column_title_gp=gpar(fontsize=11),
          width=1*unit(5, "mm"),
          height=nrow(Meth_RNA_percluster_cor.pval)*unit(5, "mm")
  )

ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_cor + ht_pval, heatmap_legend_side="left", merge_legend=TRUE)

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_MapToSuperclusters_heatmap.pdf", height=3, width=2.5)
ht_final
dev.off()

#####################################



### Analyze the percentage of the consistent genes
#####################################
###

### Load the data
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))])
Cor_Consistency_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")

### Plot the UpSet plot on hypermethylated-downregulated genes
PERCENTs_hyperNeg=list()
upset_plot=list()
for (i in 1:length(tissues_selected)) {
  Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho>0 & analysis=="all") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho<0 & analysis=="all") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  F_Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho>0 & analysis=="females") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  F_Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho<0 & analysis=="females") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  M_Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho>0 & analysis=="males") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  M_Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho<0 & analysis=="males") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()

  library(ComplexHeatmap)
  library(ggplot2)
  datalist=list(TSS.both=Hyper.TSS_genes, RNA.both=Neg.RNA_genes, TSS.F=F_Hyper.TSS_genes, RNA.F=F_Neg.RNA_genes, TSS.M=M_Hyper.TSS_genes, RNA.M=M_Neg.RNA_genes)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("100100","010010","001001") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[c(1:3)]/c(length(unique(c(Neg.RNA_genes, Hyper.TSS_genes))),
                                     length(unique(c(F_Hyper.TSS_genes, F_Neg.RNA_genes))),
                                     length(unique(c(M_Hyper.TSS_genes, M_Neg.RNA_genes)))))
  # save the precentages
  PERCENTs_hyperNeg[[i]]=c(percents_, tissues_selected[i])
  percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
  title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in\n",tissues_selected[i]), x=unit(0, "mm"), just="left"))}

  # plot
  upset_plot_=
    HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
    UpSet(m,
          top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
          right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
          row_names_gp=gpar(fontsize=9),
          set_order=c("TSS.both","RNA.both","TSS.F","RNA.F","TSS.M","RNA.M"),
          comb_order=1:3) %v%
    HeatmapAnnotation(
      empty=anno_empty(border=F, height=unit(2,"mm")),
      text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
    )
  upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_UpSet_hyper.neg.pdf", height=3.5, width=3.5)
for (i in 1:length(tissues_selected)) plot(upset_plot[[i]])
dev.off()

### Plot the UpSet plot on hypomethylated-upregulated genes
PERCENTs_hypoPos=list()
upset_plot=list()
for (i in 1:length(tissues_selected)) {
  Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho<0 & analysis=="all") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho>0 & analysis=="all") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  F_Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho<0 & analysis=="females") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  F_Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho>0 & analysis=="females") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  M_Hyper.TSS_genes=TSS_Cor %>% subset(celltypes==tissues_selected[i] & rho_pval<0.05 & rho<0 & analysis=="males") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  M_Neg.RNA_genes=RNA_Cor %>% subset(celltypes==shared_celltype[[i]] & rho_pval<0.05 & rho>0 & analysis=="males") %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  
  library(ComplexHeatmap)
  library(ggplot2)
  datalist=list(TSS.both=Hyper.TSS_genes, RNA.both=Neg.RNA_genes, TSS.F=F_Hyper.TSS_genes, RNA.F=F_Neg.RNA_genes, TSS.M=M_Hyper.TSS_genes, RNA.M=M_Neg.RNA_genes)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("100100","010010","001001") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[c(1:3)]/c(length(unique(c(Neg.RNA_genes, Hyper.TSS_genes))),
                                     length(unique(c(F_Hyper.TSS_genes, F_Neg.RNA_genes))),
                                     length(unique(c(M_Hyper.TSS_genes, M_Neg.RNA_genes)))))
  # save the precentages
  PERCENTs_hypoPos[[i]]=c(percents_, tissues_selected[i])
  percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
  title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in\n",tissues_selected[i]), x=unit(0, "mm"), just="left"))}
  
  # plot
  upset_plot_=
    HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
    UpSet(m,
          top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
          right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
          row_names_gp=gpar(fontsize=9),
          set_order=c("TSS.both","RNA.both","TSS.F","RNA.F","TSS.M","RNA.M"),
          comb_order=1:3) %v%
    HeatmapAnnotation(
      empty=anno_empty(border=F, height=unit(2,"mm")),
      text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
    )
  upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_UpSet_hypo.pos.pdf", height=3.5, width=3.5)
for (i in 1:length(tissues_selected)) plot(upset_plot[[i]])
dev.off()

### Save all the precentages for plotting the females vs. males
names_precents=c("(TSS.both & RNA.both)/(TSS.both | RNA.both)", "(TSS.F & RNA.F)/(TSS.F | RNA.F)","(TSS.M & RNA.M)/(TSS.M | RNA.M)")
hyperneg_df=data.frame(c1=unlist(lapply(PERCENTs_hyperNeg, function(x) x[1])),
                    c2=unlist(lapply(PERCENTs_hyperNeg, function(x) x[2])),
                    c3=unlist(lapply(PERCENTs_hyperNeg, function(x) x[3])),
                    celltype=unlist(lapply(PERCENTs_hyperNeg, function(x) x[4])),
                    analysis.orientation="Hyper_Neg")
hypopos_df=data.frame(c1=unlist(lapply(PERCENTs_hypoPos, function(x) x[1])),
                    c2=unlist(lapply(PERCENTs_hypoPos, function(x) x[2])),
                    c3=unlist(lapply(PERCENTs_hypoPos, function(x) x[3])),
                    celltype=unlist(lapply(PERCENTs_hypoPos, function(x) x[4])),
                    analysis.orientation="Hypo_Pos")
DFs=data.table::rbindlist(list(hyperneg_df, hypopos_df))
colnames(DFs)=c(names_precents, "celltype","analysis.orientation")
write.table(DFs, "~/Project_PBMCage/Results/PBMC_results/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_UpSet_percents.txt", sep="\t")

### Plot the percentages of overlapping between TSS. and RNA. in the Upset analysis for females and males
# Load the precent results
DFs=read.delim("~/Project_PBMCage/Results/PBMC_results/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_UpSet_percents.txt", sep="\t")
names_cols=c("(TSS.both & RNA.both)/(TSS.both | RNA.both)", "(TSS.F & RNA.F)/(TSS.F | RNA.F)","(TSS.M & RNA.M)/(TSS.M | RNA.M)",
             "celltype","orientation")
colnames(DFs)=names_cols
DFs=DFs %>% tidyr::pivot_longer(cols=names_cols[1:3], names_to="comparison", values_to="percent") %>%
  mutate(orientation=ifelse(orientation=="Hyper_Neg", "hyper._downreg.", "hypo._upreg.")) %>%
  subset(comparison %in% c("(TSS.F & RNA.F)/(TSS.F | RNA.F)","(TSS.M & RNA.M)/(TSS.M | RNA.M)")) %>%
  mutate(comparison=ifelse(comparison=="(TSS.M & RNA.M)/(TSS.M | RNA.M)", "M", "F")) %>%
  mutate(celltype=ifelse(celltype=="peripheral blood mononuclear cell","PBMC",celltype))
DFs$orientation=fct_relevel(DFs$orientation, "hypo._upreg.","hyper._downreg.")
DFs$percent=100*(DFs$percent)
# Plot
library(ggplot2)
# for rough
plot_=
  ggplot(DFs, aes(x=celltype, y=percent, color=comparison)) +
  geom_linerange(data=subset(DFs, orientation=="hypo._upreg."),
                 aes(ymin=0, ymax=percent, color=comparison), size=1, alpha=0.9, position=position_dodge(width=0.2)) +
  geom_linerange(data=subset(DFs, orientation=="hyper._downreg."),
                 aes(ymin=0, ymax=-percent, color=comparison), size=1, alpha=0.9, position=position_dodge(width=0.2)) +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=8, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(breaks=seq(0,16,2)) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1)))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_UpSet_TSSvsRNA_overlap_percents.pdf", height=3, width=3)
plot(plot_)
dev.off()

#####################################



### Plot the jaccard similarity of genes between TSS and RNA
#####################################
### 

### Load the data
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))])

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Calculate the jaccard similarity for hypo.pos genes
TSS_Cor_genes=subset(TSS_Cor, rho<0 & rho_pval<0.05 & analysis=="all")
TSS_Cor_genes=split(TSS_Cor_genes$gene, TSS_Cor_genes$celltypes)
names(TSS_Cor_genes) # check the order
RNA_Cor_genes_df=subset(RNA_Cor, rho>0 & rho_pval<0.05 & analysis=="all")
RNA_Cor_genes=list()
for (i in 1:length(shared_celltype)) {
  RNA_Cor_genes[[i]]=subset(RNA_Cor_genes_df, celltypes %in% shared_celltype[[i]])
  RNA_Cor_genes[[i]]=split(RNA_Cor_genes[[i]]$gene, RNA_Cor_genes[[i]]$celltypes)[[1]]
}
names(RNA_Cor_genes)=names(TSS_Cor_genes)

Jaccard_SIM=Jaccard_SIM_Name=c()
percent_of_cor.RNA.genes=c()
for (i in 1:length(TSS_Cor_genes)) {
  the_TSS_celltype=TSS_Cor_genes[[i]]
  the_RNA_celltype=RNA_Cor_genes[[i]]
  
  jaccard_similarity=jaccard(the_TSS_celltype, the_RNA_celltype)
  Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
  Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(TSS_Cor_genes)[i],"_vs_",names(RNA_Cor_genes)[i]))
  
  percent_of_cor.RNA.genes=c(percent_of_cor.RNA.genes, length(intersect(the_TSS_celltype, the_RNA_celltype))/length(the_RNA_celltype))
}
names(Jaccard_SIM)=names(percent_of_cor.RNA.genes)=Jaccard_SIM_Name
Jaccard_SIM_pos=data.frame(pair=names(Jaccard_SIM), hypo.pos.jaccard=Jaccard_SIM, hypo.pos.percent=percent_of_cor.RNA.genes)

### Calculate the jaccard similarity for hyper.neg genes
TSS_Cor_genes=subset(TSS_Cor, rho>0 & rho_pval<0.05 & analysis=="all")
TSS_Cor_genes=split(TSS_Cor_genes$gene, TSS_Cor_genes$celltypes)
names(TSS_Cor_genes) # check the order
RNA_Cor_genes_df=subset(RNA_Cor, rho<0 & rho_pval<0.05 & analysis=="all")
RNA_Cor_genes=list()
for (i in 1:length(shared_celltype)) {
  RNA_Cor_genes[[i]]=subset(RNA_Cor_genes_df, celltypes %in% shared_celltype[[i]])
  RNA_Cor_genes[[i]]=split(RNA_Cor_genes[[i]]$gene, RNA_Cor_genes[[i]]$celltypes)[[1]]
}
names(RNA_Cor_genes)=names(TSS_Cor_genes)

Jaccard_SIM=Jaccard_SIM_Name=percent_of_cor.RNA.genes=c()
for (i in 1:length(TSS_Cor_genes)) {
  the_TSS_celltype=TSS_Cor_genes[[i]]
  the_RNA_celltype=RNA_Cor_genes[[i]]
  
  jaccard_similarity=jaccard(the_TSS_celltype, the_RNA_celltype)
  Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
  Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(TSS_Cor_genes)[i],"_vs_",names(RNA_Cor_genes)[i]))
  
  percent_of_cor.RNA.genes=c(percent_of_cor.RNA.genes, length(intersect(the_TSS_celltype, the_RNA_celltype))/length(the_RNA_celltype))
}
names(Jaccard_SIM)=names(percent_of_cor.RNA.genes)=Jaccard_SIM_Name
Jaccard_SIM_neg=data.frame(pair=names(Jaccard_SIM), hyper.neg.jaccard=Jaccard_SIM, hyper.neg.percent=percent_of_cor.RNA.genes)

### Merge the pos and neg results
# jaccard
Jaccard_SIM_merged=merge(Jaccard_SIM_pos, Jaccard_SIM_neg, by="pair")
Jaccard_SIM_merged=Jaccard_SIM_merged %>% 
  tidyr::pivot_longer(cols=c("hyper.neg.jaccard","hypo.pos.jaccard"), names_to="orientation", values_to="jaccard") %>%
  mutate(Tss.celltype=gsub("_vs_.*","",pair), RNA.celltypecluster=gsub(".*_vs_","",pair))
Jaccard_SIM_merged$orientation=forcats::fct_relevel(Jaccard_SIM_merged$orientation, c("hypo.pos.jaccard","hyper.neg.jaccard"))
# percent in age-correlated transcripts
Percent_SIM_merged=merge(Jaccard_SIM_pos, Jaccard_SIM_neg, by="pair")
Percent_SIM_merged=Percent_SIM_merged %>% 
  tidyr::pivot_longer(cols=c("hyper.neg.percent","hypo.pos.percent"), names_to="orientation", values_to="percent") %>%
  mutate(Tss.celltype=gsub("_vs_.*","",pair), RNA.celltypecluster=gsub(".*_vs_","",pair))
Percent_SIM_merged$orientation=forcats::fct_relevel(Percent_SIM_merged$orientation, c("hypo.pos.percent","hyper.neg.percent"))

### Plot the results
plot_jaccard=
  ggplot(Jaccard_SIM_merged, aes(x=Tss.celltype, y=jaccard, alpha=orientation)) +
  geom_bar(stat="identity", width=0.5) +
  scale_alpha_manual(values=c(0.3,0.7), labels=c("hypo. upreg.","hyper. downreg."), guide_legend(title="")) +
  scale_x_discrete(labels=scales::label_wrap(20)) +
  theme_classic() +
  labs(x=NULL, y="Jaccard similarity\nindex", title=NULL) +
  theme(axis.text.x=element_text(size=11, angle=60, hjust=1, vjust=1),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        legend.position="right",
        title=element_text(size=10)) +
  scale_x_discrete(labels=scales::label_wrap(20))
  
# plot_percent=
  ggplot(Percent_SIM_merged, aes(x=Tss.celltype, y=percent, alpha=orientation)) +
  geom_bar(stat="identity", width=0.3) +
  scale_alpha_manual(values=c(0.3,0.7), labels=c("hypo. upreg.","hyper. downreg."), guide_legend(title="")) +
  scale_x_discrete(labels=scales::label_wrap(20)) +
  theme_classic() +
  labs(x=NULL, y="Jaccard similarity index", title=NULL) +
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        legend.position="right",
        title=element_text(size=10))
pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_JaccardSimilarity_TSSvsRNA.pdf", height=2.5, width=5)
plot(plot_jaccard)
dev.off()

#####################################



### Plot the percentages of overlapping between age-cor TSS in age-cor RNA
#####################################
### 

library(dplyr)

### Load the data
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# select only the celltypes of interest
TSS_Cor=TSS_Cor %>%
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(analysis=="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes)
RNA_Cor=RNA_Cor %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  subset(analysis=="all") %>%
  dplyr::select(gene, rho, rho_pval, celltypes)

Cor_Consistency_df=inner_join(TSS_Cor, RNA_Cor, by=c("gene","celltypes"), suffix=c(".TSS",".RNA"))

### Analyze CD4 and CD8T
TSS_subset=Cor_Consistency_df %>% 
  subset(rho_pval.RNA<0.05 & abs(rho.RNA)>0.1) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  mutate(TSS.orientation=ifelse(rho_pval.TSS>=0.05, "age-uncorr.", NA)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS<0 & rho.RNA<0 & rho_pval.TSS<0.05, "promiscuous", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS<0 & rho.RNA>0 & rho_pval.TSS<0.05, "hypo. & upreg.", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS>0 & rho.RNA<0 & rho_pval.TSS<0.05, "hyper. & downreg.", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS>0 & rho.RNA>0 & rho_pval.TSS<0.05, "promiscuous", TSS.orientation)) %>%
  subset(celltypes!="Monocytes") # remove Monocytes

plot_=
  ggplot(TSS_subset, aes(x=TSS.orientation, fill=TSS.orientation)) +
  facet_wrap(~celltypes) +
  geom_bar(width=1, stat="count", color="black", linewidth=0.2) +
  coord_polar("x", start=0) +
  theme_light() +
  scale_fill_manual(values=c("transparent","dodgerblue3","brown3","grey80")) +
  theme(plot.margin=margin(-3,5,-3,5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10, color="black"),
        strip.background=element_rect(fill="white", color="black")) +
  guides(fill=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_Proportion_CD4T.CD8T.pdf", height=3, width=6)
plot(plot_)
dev.off()

### Number of the right direction
TSS_subset=Cor_Consistency_df %>% 
  subset(rho_pval.RNA<0.05 & abs(rho.RNA)>0.1) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  mutate(TSS.orientation=ifelse(rho_pval.TSS>=0.05, "age-uncorr.", NA)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS<0 & rho.RNA<0 & rho_pval.TSS<0.05, "promiscuous", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS<0 & rho.RNA>0 & rho_pval.TSS<0.05, "hypo. & upreg.", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS>0 & rho.RNA<0 & rho_pval.TSS<0.05, "hyper. & downreg.", TSS.orientation)) %>%
  mutate(TSS.orientation=ifelse(rho.TSS>0 & rho.RNA>0 & rho_pval.TSS<0.05, "promiscuous", TSS.orientation))

plot_2=
  ggplot(subset(TSS_subset, TSS.orientation %in% c("hypo. & upreg.","hyper. & downreg.")), 
         aes(x=celltypes, color=TSS.orientation)) +
  geom_bar(stat="count", position=position_stack(), width=0.7, fill="transparent") +
  scale_color_manual(values=c("dodgerblue3","brown3")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=10, color="black"),
        strip.background=element_rect(fill="white", color="black")) +
  # scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  guides(color=guide_legend(title=NULL))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_Count_CD4T.CD8T.Mono.pdf", height=2, width=5)
plot(plot_2)
dev.off()

#####################################



### Plot the percentages of overlapping between females and males in TSS-cor genes
#####################################
###

### Load the data
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))])

### Clean the data
TSS_result_DF=list()
for (i in 1:length(unique(TSS_Cor$celltypes))) {
  TSS_Cor_genes_M_df=subset(TSS_Cor, rho_pval<0.05 & analysis=="males" & celltypes==unique(TSS_Cor$celltypes)[i]) %>% mutate(orientation=ifelse(rho>0,"hyper","hypo"))
  TSS_Cor_genes_M=split(TSS_Cor_genes_M_df$gene, TSS_Cor_genes_M_df$orientation)
  TSS_Cor_genes_F_df=subset(TSS_Cor, rho_pval<0.05 & analysis=="females" & celltypes==unique(TSS_Cor$celltypes)[i]) %>% mutate(orientation=ifelse(rho>0,"hyper","hypo"))
  TSS_Cor_genes_F=split(TSS_Cor_genes_F_df$gene, TSS_Cor_genes_F_df$orientation)
  TSS_result_DF[[i]]=data.table::rbindlist(list(data.frame(overlap.percent=length(intersect(TSS_Cor_genes_M[["hyper"]], TSS_Cor_genes_F[["hyper"]]))/length(unique(c(TSS_Cor_genes_M[["hyper"]], TSS_Cor_genes_F[["hyper"]]))),
                                                           celltypes=unique(TSS_Cor$celltypes)[i],
                                                           orientation="hyper"),
                                                data.frame(overlap.percent=length(intersect(TSS_Cor_genes_M[["hypo"]], TSS_Cor_genes_F[["hypo"]]))/length(unique(c(TSS_Cor_genes_M[["hypo"]], TSS_Cor_genes_F[["hypo"]]))),
                                                           celltypes=unique(TSS_Cor$celltypes)[i],
                                                           orientation="hypo")))
}
DFs_FvsM=data.table::rbindlist(TSS_result_DF)

### Plot
library(ggplot2)
# plot_rough=
  ggplot(DFs_FvsM, aes(x=celltypes, y=overlap.percent, color=orientation)) +
  geom_linerange(data=subset(DFs_FvsM, orientation=="hyper"),
                 aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=overlap.percent), size=1, alpha=0.9) +
  geom_linerange(data=subset(DFs_FvsM, orientation=="hypo"),
                 aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=-overlap.percent), size=1, alpha=0.9) +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=11, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        title=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(labels=abs) +
  scale_x_discrete(labels=scales::label_wrap(20))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/TSS.DMPs_Correlation_wAges_FvsM.pdf", height=3.5, width=4.5)
plot(plot_rough)
dev.off()

#####################################



### Plot the percentages of overlapping between females and males in TSS/RNA consistent genes
#####################################
###

### Load the data
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
RNA_Cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))],
                     unique(RNA_Cor$celltypes)[!grepl("Other cells", unique(RNA_Cor$celltypes))])
Cor_Consistency_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")

### Clean the data
Consis_result_DF=list()
for (i in 1:length(unique(Cor_Consistency_DF$celltypes.TSS))) {
  Consis_genes_M_df=subset(Cor_Consistency_DF, analysis=="males" & celltypes.TSS==unique(Cor_Consistency_DF$celltypes.TSS)[i])
  Consis_genes_M=split(Consis_genes_M_df$gene, Consis_genes_M_df$analysis.orientation)
  Consis_genes_F_df=subset(Cor_Consistency_DF, analysis=="females" & celltypes.TSS==unique(Cor_Consistency_DF$celltypes.TSS)[i])
  Consis_genes_F=split(Consis_genes_F_df$gene, Consis_genes_F_df$analysis.orientation)
  Consis_result_DF[[i]]=data.table::rbindlist(list(data.frame(overlap.percent=length(intersect(Consis_genes_M[["Hyper_Neg"]], Consis_genes_F[["Hyper_Neg"]]))/length(unique(c(Consis_genes_M[["Hyper_Neg"]], Consis_genes_F[["Hyper_Neg"]]))),
                                                              celltypes=unique(Cor_Consistency_DF$celltypes.TSS)[i],
                                                              orientation="hyper. downreg."),
                                                   data.frame(overlap.percent=length(intersect(Consis_genes_M[["Hypo_Pos"]], Consis_genes_F[["Hypo_Pos"]]))/length(unique(c(Consis_genes_M[["Hypo_Pos"]], Consis_genes_F[["Hypo_Pos"]]))),
                                                              celltypes=unique(Cor_Consistency_DF$celltypes.TSS)[i],
                                                              orientation="hypo. upreg.")))
}
DFs_FvsM=data.table::rbindlist(Consis_result_DF)

### Plot
library(ggplot2)
plot_rough=
  ggplot(DFs_FvsM, aes(x=celltypes, y=overlap.percent, color=orientation)) +
  geom_linerange(data=subset(DFs_FvsM, orientation=="hyper. downreg."),
                 aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=overlap.percent), size=1, alpha=0.9) +
  geom_linerange(data=subset(DFs_FvsM, orientation=="hypo. upreg."),
                 aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=-overlap.percent), size=1, alpha=0.9) +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of\noverlap (%)") +
  theme(axis.text.x=element_text(size=11, angle=60, hjust=1, vjust=1),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=11),
        title=element_text(size=12),
        legend.position="none") +
  scale_y_continuous(labels=abs) +
  scale_x_discrete(labels=scales::label_wrap(20))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_Cor.Consistency_FvsM.pdf", height=2.5, width=4)
plot(plot_rough)
dev.off()

#####################################



### Analyze the correlation between the TSS beta-values and RNA expr of consistently age-correlated genes
#####################################
###

library(dplyr)

### Get the RNA and Methylation data
pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
TSS.beta=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")

### Get the consistently age-correlated genes
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))],
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))])
Cor_Consistency_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Cor.Consistency.csv.gz", sep="\t")

### Analyze
Cor_analysis_DF_list=list()
for (i in 1:length(tissues_selected)) {
  Hyper_Neg.Consis.genes=Cor_Consistency_DF %>%
    subset(celltypes.TSS==tissues_selected[i] & analysis.orientation=="Hyper_Neg" & analysis %in% c("females","males")) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  # extract the expr
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% shared_celltype[[i]])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Hyper_Neg.Consis.genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("^age$|^sex$|^donor_id$",colnames(.))], names_to="gene", values_to="expression") %>%
    group_by(age, sex, gene) %>%
    summarize_at("expression", mean)
  # extract the beta
  TSS_subset=subset(TSS.beta, celltypes==tissues_selected[i] & (sex %in% c("F","M")) & (gene %in% Hyper_Neg.Consis.genes)) %>%
    group_by(age, sex, gene) %>%
    summarize_at("beta", mean)
  # merge the expr data and the beta value
  merged_data=dplyr::inner_join(TSS_subset, expr_data, by=c("age","sex","gene"), suffix=c(".TSS",".RNA"))
  # remove the genes with <3 observations
  genes_=merged_data %>% group_by(gene) %>% count() %>% subset(n>=3) %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  genes_f=merged_data %>% subset(sex=="F") %>% group_by(gene) %>% count() %>% subset(n>=3) %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  genes_m=merged_data %>% subset(sex=="M") %>% group_by(gene) %>% count() %>% subset(n>=3) %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()

  # analyze correlation with ages for all/females/males
  all_cor=all_cor_pval=all_rho=all_rho_pval=all_tau=all_tau_pval=list()
  if (length(genes_)>=1) {
    for (j in 1:length(genes_)) {
      df_=subset(merged_data, gene==genes_[j])
      cor_result=cor.test(df_$beta, df_$expression, method="pearson")
      all_cor[[j]]=cor_result$estimate
      all_cor_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="spearman", exact=F)
      all_rho[[j]]=cor_result$estimate
      all_rho_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="kendall")
      all_tau[[j]]=cor_result$estimate
      all_tau_pval[[j]]=cor_result$p.value
    }
    all_cor=unname(unlist(all_cor)); all_cor_pval=unname(unlist(all_cor_pval))
    all_rho=unname(unlist(all_rho)); all_rho_pval=unname(unlist(all_rho_pval))
    all_tau=unname(unlist(all_tau)); all_tau_pval=unname(unlist(all_tau_pval))
  }
  female_cor=female_cor_pval=female_rho=female_rho_pval=female_tau=female_tau_pval=list()
  if (length(genes_f)>=1) {
    for (j in 1:length(genes_f)) {
      df_=subset(merged_data, sex=="F" & gene==genes_f[j])
      cor_result=cor.test(df_$beta, df_$expression, method="pearson")
      female_cor[[j]]=cor_result$estimate
      female_cor_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="spearman", exact=F)
      female_rho[[j]]=cor_result$estimate
      female_rho_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="kendall")
      female_tau[[j]]=cor_result$estimate
      female_tau_pval[[j]]=cor_result$p.value
    }
    female_cor=unname(unlist(female_cor)); female_cor_pval=unname(unlist(female_cor_pval))
    female_rho=unname(unlist(female_rho)); female_rho_pval=unname(unlist(female_rho_pval))
    female_tau=unname(unlist(female_tau)); female_tau_pval=unname(unlist(female_tau_pval))
  }
  male_cor=male_cor_pval=male_rho=male_rho_pval=male_tau=male_tau_pval=list()
  if (length(genes_m)>=1) {
    for (j in 1:length(genes_m)) {
      df_=subset(merged_data, sex=="F" & gene==genes_f[j])
      cor_result=cor.test(df_$beta, df_$expression, method="pearson")
      male_cor[[j]]=cor_result$estimate
      male_cor_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="spearman", exact=F)
      male_rho[[j]]=cor_result$estimate
      male_rho_pval[[j]]=cor_result$p.value
      
      cor_result=cor.test(df_$beta, df_$expression, method="kendall")
      male_tau[[j]]=cor_result$estimate
      male_tau_pval[[j]]=cor_result$p.value
    }
    male_cor=unname(unlist(male_cor)); male_cor_pval=unname(unlist(male_cor_pval))
    male_rho=unname(unlist(male_rho)); male_rho_pval=unname(unlist(male_rho_pval))
    male_tau=unname(unlist(male_tau)); male_tau_pval=unname(unlist(male_tau_pval))
  }

  # combine the correlation analysis results in dataframe
  ALL_COR_DF=Female_COR_DF=Male_COR_DF=data.frame()
  if (length(genes_)>=1) {
    ALL_COR_DF=data.frame(gene=genes_, 
                          cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval,
                          analysis="all", celltypes=tissues_selected[i])
  }
  if (length(genes_f)>=1) {
    Female_COR_DF=data.frame(gene=genes_f, 
                             cor=female_cor, cor_pval=female_cor_pval, rho=female_rho, rho_pval=female_rho_pval, tau=female_tau, tau_pval=female_tau_pval,
                             analysis="females", celltypes=tissues_selected[i])
  }
  if (length(genes_m)>=1) {
    Male_COR_DF=data.frame(gene=genes_m, 
                           cor=male_cor, cor_pval=male_cor_pval, rho=male_rho, rho_pval=male_rho_pval, tau=male_tau, tau_pval=male_tau_pval,
                           analysis="males", celltypes=tissues_selected[i])
  }
  Cor_analysis_DF_list[[i]]=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))
}

Cor_analysis_DF=data.table::rbindlist(Cor_analysis_DF_list)
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Correlation.csv.gz", sep="\t")


###
### Plot the correlation results
# get the RNA and Methylation data
pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
TSS.beta=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")
# get the correlation results
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
shared_celltype=list(c("Monocytes"),
                     c("CD4T cells"),
                     c("CD8T cells"),
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))],
                     c("B cells", "CD4T cells", "CD8T cells", "NK cells", "OtherT"),
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))],
                     unique(pseudobulkobj$Annot.rough)[!grepl("Other cells", unique(pseudobulkobj$Annot.rough))])
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_RNA.PBMCage_Correlation.csv.gz", sep="\t")
Cor_analysis_DF=Cor_analysis_DF %>% subset(analysis %in% c("females","males") & cor_pval<0.05 & cor<0 & rho_pval<0.05 & rho<0)
# update the celltypes that exist in the filtered Correlation DF
select_=tissues_selected %in% unique(Cor_analysis_DF$celltypes)
tissues_selected=tissues_selected[select_]
shared_celltype=shared_celltype[select_]
cor_plots=list()
for (i in 1:length(tissues_selected)) {
  # determine the correlated genes
  genes_=Cor_analysis_DF %>% subset(celltypes==tissues_selected[i]) %>% arrange(cor) %>% select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unique() %>% unname()
  # extract the expr
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% shared_celltype[[i]])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", genes_), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!grepl("^age$|^sex$|^donor_id$",colnames(.))], names_to="gene", values_to="expression") %>%
    group_by(age, sex, gene) %>%
    summarize_at("expression", mean)
  # extract the beta
  TSS_subset=subset(TSS.beta, celltypes==tissues_selected[i] & (sex %in% c("F","M")) & (gene %in% genes_)) %>%
    group_by(age, sex, gene) %>%
    summarize_at("beta", mean)
  # merge the expr data and the beta value
  merged_data=dplyr::inner_join(TSS_subset, expr_data, by=c("age","sex","gene"), suffix=c(".TSS",".RNA"))
  # plot
  cor_plots[[i]]=
    ggplot(merged_data, aes(x=beta, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.5) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="Beta-value", y="Expression", title=paste0("Correlation between methylation and expression in ", tissues_selected[i])) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          title=element_text(size=10)) +
    guides(guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))) +
    geom_smooth(show.legend=FALSE, method='lm', linewidth=0.5, fill="gray93", alpha=0.5) +
    ggpubr::stat_cor(show.legend=FALSE, size=3.5, method="pearson")
  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",i,".pdf"), height=ceiling(length(unique(merged_data$gene))/3)*2, width=9)
  plot(cor_plots[[i]])
  dev.off()
}
qpdf::pdf_combine(input=paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",c(1:length(tissues_selected)),".pdf"),
                  output="~/Project_PBMCage/Other_omics/Methy_plots/AgeMeth.TSS_RNA.PBMCage_CorPlots.pdf")
unlink(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",c(1:length(tissues_selected)),".pdf"))

#####################################



### Extract and enrich the age-associated DMRs
#####################################
###

###
library(dplyr)
library(tidyverse)

### Load the methylation data
Methylation_Data=readRDS("~/Project_PBMCage/Other_omics/AgeMeth_data.rds")

### Extract the pos. and neg. age-associated DMRs in each celltype
TERM2GENE_list=list()
tissues_selected=c("CD14+ monocyte","CD4+ T cell","CD8+ T cell","leukocyte","lymphocyte","peripheral blood mononuclear cell","whole blood")
for (i in 1:length(Methylation_Data)) {
  celltype=names(Methylation_Data)[i]

  # load the DMR result
  results_per_celltype=load(paste0("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_",celltype,"_results.RData"))
  DMR_=myDMR[[1]]

  # take the sig. DMRs
  myDMR_sig_genes=DMR_ %>%
    mutate(abs_meandiff=abs(meandiff)) %>%
    subset(Fisher<0.05 & abs_meandiff!=0)
  myDMR_sig_genes=split(myDMR_sig_genes$overlapping.genes, myDMR_sig_genes$seqnames)
  myDMR_sig_genes=lapply(myDMR_sig_genes, function(x) {unlist(lapply(x, function(y) str_split(y, pattern=", ")[[1]]))})
  myDMR_sig_genes=lapply(myDMR_sig_genes, function(x) x[(!duplicated(x)) & (!is.na(x))])
  chr_names=names(myDMR_sig_genes)
  myDMR_sig_genes_df=lapply(1:length(chr_names), function(x) {
    data.frame(term=rep(paste0(tissues_selected[i], ".sig.", chr_names[x]), length(myDMR_sig_genes[[x]])),
               gene=myDMR_sig_genes[[x]])
  })
  myDMR_sig_genes_df=data.table::rbindlist(myDMR_sig_genes_df)

  # take the pos. sig. DMRs
  myDMR_pos_genes=DMR_ %>%
    subset(Fisher<0.05 & meandiff>0)
  myDMR_pos_genes=split(myDMR_pos_genes$overlapping.genes, myDMR_pos_genes$seqnames)
  myDMR_pos_genes=lapply(myDMR_pos_genes, function(x) {unlist(lapply(x, function(y) str_split(y, pattern=", ")[[1]]))})
  myDMR_pos_genes=lapply(myDMR_pos_genes, function(x) x[(!duplicated(x)) & (!is.na(x))])
  chr_names=names(myDMR_pos_genes)
  myDMR_pos_genes_df=lapply(1:length(chr_names), function(x) {
    data.frame(term=rep(paste0(tissues_selected[i], ".pos.", chr_names[x]), length(myDMR_pos_genes[[x]])),
               gene=myDMR_pos_genes[[x]])
  })
  myDMR_pos_genes_df=data.table::rbindlist(myDMR_pos_genes_df)

  # take the neg. sig. DMRs
  myDMR_neg_genes=DMR_ %>%
    subset(Fisher<0.05 & meandiff<0)
  myDMR_neg_genes=split(myDMR_neg_genes$overlapping.genes, myDMR_neg_genes$seqnames)
  myDMR_neg_genes=lapply(myDMR_neg_genes, function(x) {unlist(lapply(x, function(y) str_split(y, pattern=", ")[[1]]))})
  myDMR_neg_genes=lapply(myDMR_neg_genes, function(x) x[(!duplicated(x)) & (!is.na(x))])
  chr_names=names(myDMR_neg_genes)
  myDMR_neg_genes_df=lapply(1:length(chr_names), function(x) {
    data.frame(term=rep(paste0(tissues_selected[i], ".neg.", chr_names[x]), length(myDMR_neg_genes[[x]])),
               gene=myDMR_neg_genes[[x]])
  })
  myDMR_neg_genes_df=data.table::rbindlist(myDMR_neg_genes_df)

  # create the TermTOGene dataframe for later enrichment analysis
  TERM2GENE_list[[i]]=data.table::rbindlist(list(
    myDMR_sig_genes_df,
    myDMR_pos_genes_df,
    myDMR_neg_genes_df
  ))
}
TERM2GENE_df=data.table::rbindlist(TERM2GENE_list)
saveRDS(TERM2GENE_df, "~/Project_PBMCage/Other_omics/Methy_results/age_methylation_results_sigDMRs.rds")

### Add the meta info
TERM2GENE_df=readRDS("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_results_sigDMRs.rds")
TERM2GENE_df$celltype=sapply(TERM2GENE_df$term, function(x) strsplit(x, split="\\.")[[1]][1])
TERM2GENE_df$orientation=sapply(TERM2GENE_df$term, function(x) strsplit(x, split="\\.")[[1]][2])
TERM2GENE_df$chromosome=sapply(TERM2GENE_df$term, function(x) strsplit(x, split="\\.")[[1]][3])
saveRDS(TERM2GENE_df, "~/Project_PBMCage/Other_omics/Methy_results/age_methylation_results_sigDMRs.rds")

#####################################



### Enrich the age-associated DMRs
#####################################
###

### Get the all DMR sig. genes from each celltype
TERM2GENE_df=readRDS("~/Project_PBMCage/Other_omics/Methy_results/age_methylation_results_sigDMRs.rds")
TERM2GENE_df$celltype.orient=paste0(TERM2GENE_df$celltype, "_", TERM2GENE_df$orientation)
gene_split=split(TERM2GENE_df$gene, TERM2GENE_df$celltype.orient)

# enrich
MODULE_ENRICH=list()
celltypes=names(table(TERM2GENE_df$celltype))
for (i in 1:length(celltypes)) {
  celltype_=celltypes[i]
  genes_list=gene_split[grepl(gsub("\\+","\\\\+",celltype_), names(gene_split))]
  # run compareCluster to get the terms in which up.genes and down.genes are enriched
  comparecluster_res=compareCluster(genes_list,
                                    fun="enrichGO",
                                    keyType="SYMBOL",
                                    OrgDb="org.Hs.eg.db",
                                    ont="BP",
                                    pvalueCutoff=0.05,
                                    minGSSize=30,
                                    maxGSSize=500)
  pos_ids=subset(comparecluster_res@compareClusterResult, grepl("_pos", Cluster))$ID
  neg_ids=subset(comparecluster_res@compareClusterResult, grepl("_neg", Cluster))$ID
  sig_ids=subset(comparecluster_res@compareClusterResult, grepl("_sig", Cluster))$ID
  # run the enrichGO
  module_enrich=enrichGO(unname(unlist(genes_list)),
                         keyType="SYMBOL",
                         OrgDb="org.Hs.eg.db",
                         ont="BP",
                         pvalueCutoff=0.05,
                         minGSSize=30,
                         maxGSSize=500)
  # add the up and down info
  module_enrich@result=module_enrich@result %>%
    mutate(`.sign`=ifelse(ID %in% pos_ids, "hypermethylation",
                          ifelse(ID %in% neg_ids, "hypomethylation", NA)))
  module_enrich_filter=gofilter(module_enrich, level=4)

  MODULE_ENRICH[[i]]=module_enrich_filter
}
names(MODULE_ENRICH)=celltypes

# plot
MODULE_ENRICH_Plots=list()
for (i in 1:length(MODULE_ENRICH)) {
  if (MODULE_ENRICH[[i]]@result$p.adjust[1]<=0.05) {
    if (nrow(subset(MODULE_ENRICH[[i]]@result, p.adjust<=0.05))<3) {
      plot_=
        dotplot(MODULE_ENRICH[[i]], showCategory=10,
                label_format=50, split='.sign') +
        facet_grid(.~.sign) +
        scale_colour_gradient(low="pink", high="brown4",
                              labels= ~sprintf(fmt="%0.01e", .)) +
        scale_size_continuous(range=c(2,6)) +
        ggtitle(names(MODULE_ENRICH)[i]) +
        theme(axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title.x=element_text(size=9),
              axis.title.y=element_text(size=9),
              title=element_text(size=11)) +
        scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(0.005, 0.005)))
    } else {
      show_cat=rbind(MODULE_ENRICH[[i]]@result %>% subset(p.adjust<0.05 & `.sign`=="hypermethylation") %>% slice_min(pvalue, n=10),
                     MODULE_ENRICH[[i]]@result %>% subset(p.adjust<0.05 & `.sign`=="hypomethylation") %>% slice_min(pvalue, n=10))

      plot_=
        dotplot(MODULE_ENRICH[[i]], showCategory=10,
                label_format=50, split='.sign') +
        facet_grid(.~.sign) +
        scale_colour_gradient(limits=quantile(show_cat$p.adjust)[c(1,5)],
                              low="pink", high="brown4",
                              labels= ~sprintf(fmt="%0.01e", .)) +
        scale_size_continuous(breaks=c(min(show_cat$Count),
                                       round((min(show_cat$Count)+max(show_cat$Count))/2),
                                       max(show_cat$Count)),
                              labels=c(min(show_cat$Count),
                                       round((min(show_cat$Count)+max(show_cat$Count))/2),
                                       max(show_cat$Count)),
                              range=c(2,6)) +
        ggtitle(names(MODULE_ENRICH)[i]) +
        theme(axis.text.x=element_text(size=8),
              axis.text.y=element_text(size=8),
              axis.title.x=element_text(size=9),
              axis.title.y=element_text(size=9),
              title=element_text(size=11)) +
        scale_x_continuous(expand=expansion(mult=c(0, 0), add=c(0.005, 0.005)))
    }

    MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
  }
}

pdf("~/Project_PBMCage/Other_omics/Methy_plots/DMRs_enrichmentResults.pdf", height=5, width=7)
for (p in MODULE_ENRICH_Plots) plot(p)
dev.off()

#####################################


TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
key_regulators=c("MYC","NFKB1","NFKB2","TP53","ERBB3","JAK3","AKT1","FYN","LCK","JAK1","MAPK8","MAP3K7","MAPK1")
TSS_Cor_subset=TSS_Cor %>% subset(analysis=="all" & gene %in% key_regulators)
ggplot(TSS_Cor_subset, aes(x=rho, y=-log10(rho_pval), color=celltypes)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=gene)) +
  geom_hline(yintercept=1)



