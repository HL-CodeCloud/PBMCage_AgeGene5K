
### Analyze the correlation between the RNA_expr and ages in each celltype at the rough level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_rough.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.rough))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="Female", "F", "M"))
  nonzero_=apply(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data=expr_data[, filter]
  genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
  
  # get the females' expr
  expr_data_female=subset(expr_data, sex=="F")
  nonzero_=apply(expr_data_female[! colnames(expr_data_female) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_female=expr_data_female[, filter]
  genes_f=colnames(expr_data_female)[! colnames(expr_data_female) %in% c("age","sex", "donor_id")]
  
  # get the males' expr
  expr_data_male=subset(expr_data, sex=="M")
  nonzero_=apply(expr_data_male[! colnames(expr_data_male) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_male=expr_data_male[, filter]
  genes_m=colnames(expr_data_male)[! colnames(expr_data_male) %in% c("age","sex", "donor_id")]
  
  # analyze correlation with ages for all/females/males
  if (length(genes_)>=1) {
    pearson_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="pearson"))
    all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="spearman", exact=F))
    all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="kendall"))
    all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_f)>=1) {
    pearson_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="pearson"))
    female_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    female_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="spearman", exact=F))
    female_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    female_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="kendall"))
    female_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    female_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_m)>=1) {
    pearson_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="pearson"))
    male_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    male_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="spearman", exact=F))
    male_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    male_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="kendall"))
    male_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    male_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  
  # combine the correlation analysis results in dataframe
  ALL_COR_DF=Female_COR_DF=Male_COR_DF=data.frame()
  if (length(genes_)>=1) {
    ALL_COR_DF=data.frame(gene=genes_, 
                          cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval, 
                          analysis="all", celltypes=celltypes[i])
  }
  if (length(genes_f)>=1) {
    Female_COR_DF=data.frame(gene=genes_f, 
                             cor=female_cor, cor_pval=female_cor_pval, rho=female_rho, rho_pval=female_rho_pval, tau=female_tau, tau_pval=female_tau_pval, 
                             analysis="females", celltypes=celltypes[i])
  }
  if (length(genes_m)>=1) {
    Male_COR_DF=data.frame(gene=genes_m, 
                           cor=male_cor, cor_pval=male_cor_pval, rho=male_rho, rho_pval=male_rho_pval, tau=male_tau, tau_pval=male_tau_pval, 
                           analysis="males", celltypes=celltypes[i])
  }
  Cor_analysis_DF_list[[i]]=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))
}
Cor_analysis_DF=data.table::rbindlist(Cor_analysis_DF_list)
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the inter level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_inter.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.inter))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="Female", "F", "M"))
  nonzero_=apply(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data=expr_data[, filter]
  genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
  
  # get the females' expr
  expr_data_female=subset(expr_data, sex=="F")
  nonzero_=apply(expr_data_female[! colnames(expr_data_female) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_female=expr_data_female[, filter]
  genes_f=colnames(expr_data_female)[! colnames(expr_data_female) %in% c("age","sex", "donor_id")]
  
  # get the males' expr
  expr_data_male=subset(expr_data, sex=="M")
  nonzero_=apply(expr_data_male[! colnames(expr_data_male) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_male=expr_data_male[, filter]
  genes_m=colnames(expr_data_male)[! colnames(expr_data_male) %in% c("age","sex", "donor_id")]
  
  # analyze correlation with ages for all/females/males
  if (length(genes_)>=1) {
    pearson_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="pearson"))
    all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="spearman", exact=F))
    all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="kendall"))
    all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_f)>=1) {
    pearson_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="pearson"))
    female_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    female_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="spearman", exact=F))
    female_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    female_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="kendall"))
    female_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    female_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_m)>=1) {
    pearson_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="pearson"))
    male_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    male_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="spearman", exact=F))
    male_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    male_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="kendall"))
    male_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    male_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  
  # combine the correlation analysis results in dataframe
  ALL_COR_DF=Female_COR_DF=Male_COR_DF=data.frame()
  if (length(genes_)>=1) {
    ALL_COR_DF=data.frame(gene=genes_, 
                          cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval, 
                          analysis="all", celltypes=celltypes[i])
  }
  if (length(genes_f)>=1) {
    Female_COR_DF=data.frame(gene=genes_f, 
                             cor=female_cor, cor_pval=female_cor_pval, rho=female_rho, rho_pval=female_rho_pval, tau=female_tau, tau_pval=female_tau_pval, 
                             analysis="females", celltypes=celltypes[i])
  }
  if (length(genes_m)>=1) {
    Male_COR_DF=data.frame(gene=genes_m, 
                           cor=male_cor, cor_pval=male_cor_pval, rho=male_rho, rho_pval=male_rho_pval, tau=male_tau, tau_pval=male_tau_pval, 
                           analysis="males", celltypes=celltypes[i])
  }
  Cor_analysis_DF_list[[i]]=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))
}

Cor_analysis_DF=data.table::rbindlist(Cor_analysis_DF_list)
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the detailed level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/Immunity/Immunity_OneSamplePerDonor_Pseudobulk_detailed.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.detailed))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="Female", "F", "M"))
  nonzero_=apply(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data=expr_data[, filter]
  genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
  
  # get the females' expr
  expr_data_female=subset(expr_data, sex=="F")
  nonzero_=apply(expr_data_female[! colnames(expr_data_female) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_female=expr_data_female[, filter]
  genes_f=colnames(expr_data_female)[! colnames(expr_data_female) %in% c("age","sex", "donor_id")]
  
  # get the males' expr
  expr_data_male=subset(expr_data, sex=="M")
  nonzero_=apply(expr_data_male[! colnames(expr_data_male) %in% c("age","sex","donor_id")], 2, function(x) length(x[x!=0])>=3)
  filter=c("age", "sex", "donor_id", names(nonzero_)[nonzero_==TRUE])
  expr_data_male=expr_data_male[, filter]
  genes_m=colnames(expr_data_male)[! colnames(expr_data_male) %in% c("age","sex", "donor_id")]
  
  # analyze correlation with ages for all/females/males
  if (length(genes_)>=1) {
    pearson_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="pearson"))
    all_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    all_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="spearman", exact=F))
    all_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    all_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data))], 2, function(y) cor.test(expr_data[["age"]], y, method="kendall"))
    all_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    all_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_f)>=1) {
    pearson_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="pearson"))
    female_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    female_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="spearman", exact=F))
    female_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    female_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_female[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_female))], 2, function(y) cor.test(expr_data_female[["age"]], y, method="kendall"))
    female_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    female_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  if (length(genes_m)>=1) {
    pearson_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="pearson"))
    male_cor=unname(unlist(lapply(pearson_results, function(x) x$estimate)))
    male_cor_pval=unname(unlist(lapply(pearson_results, function(x) x$p.value)))
    
    spearman_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="spearman", exact=F))
    male_rho=unname(unlist(lapply(spearman_results, function(x) x$estimate)))
    male_rho_pval=unname(unlist(lapply(spearman_results, function(x) x$p.value)))
    
    kendall_results=apply(expr_data_male[,!grepl("^age$|^sex$|^donor_id$",colnames(expr_data_male))], 2, function(y) cor.test(expr_data_male[["age"]], y, method="kendall"))
    male_tau=unname(unlist(lapply(kendall_results, function(x) x$estimate)))
    male_tau_pval=unname(unlist(lapply(kendall_results, function(x) x$p.value)))
  }
  
  # combine the correlation analysis results in dataframe
  ALL_COR_DF=Female_COR_DF=Male_COR_DF=data.frame()
  if (length(genes_)>=1) {
    ALL_COR_DF=data.frame(gene=genes_, 
                          cor=all_cor, cor_pval=all_cor_pval, rho=all_rho, rho_pval=all_rho_pval, tau=all_tau, tau_pval=all_tau_pval, 
                          analysis="all", celltypes=celltypes[i])
  }
  if (length(genes_f)>=1) {
    Female_COR_DF=data.frame(gene=genes_f, 
                             cor=female_cor, cor_pval=female_cor_pval, rho=female_rho, rho_pval=female_rho_pval, tau=female_tau, tau_pval=female_tau_pval, 
                             analysis="females", celltypes=celltypes[i])
  }
  if (length(genes_m)>=1) {
    Male_COR_DF=data.frame(gene=genes_m, 
                           cor=male_cor, cor_pval=male_cor_pval, rho=male_rho, rho_pval=male_rho_pval, tau=male_tau, tau_pval=male_tau_pval, 
                           analysis="males", celltypes=celltypes[i])
  }
  Cor_analysis_DF_list[[i]]=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))
}

Cor_analysis_DF=data.table::rbindlist(Cor_analysis_DF_list)
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Immunity/Results/Bulkanalysis_OneSamplePerDonor_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")

#####################################



# ### Upset analysis of the sig.cor.genes in all/females/males in each celltype
# #####################################
# ###
# 
# library(dplyr)
# 
# ### At the rough level
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# all_celltypes=names(table(Cor_analysis_DF$celltypes))
# upset_plot=list()
# for (i in 1:length(all_celltypes)) {
#   Subset_=subset(Cor_analysis_DF, celltypes==all_celltypes[i])
#   all_up=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   all_down=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_up=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_down=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_up=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_down=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   
#   # make the UpSet data
#   library(ComplexHeatmap)
#   library(ggplot2)
#   datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
#   m=make_comb_mat(datalist, mode="intersect")
#   m=normalize_comb_mat(m, full_comb_sets=T)
#   set_name(m) # check
#   combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
#   m=m[, combinations_]
#   percents_=c(comb_size(m)[1:3]/set_size(m)[1],
#               comb_size(m)[4:6]/set_size(m)[4])
#   percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
#   title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in ",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
# 
#   # plot
#   upset_plot_=
#     HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
#     UpSet(m, 
#           top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
#           right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
#           row_names_gp=gpar(fontsize=9),
#           set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
#           comb_order=1:6) %v% 
#     HeatmapAnnotation(
#       empty=anno_empty(border=F, height=unit(2,"mm")),
#       text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
#     )
#   upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Rough.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()
# 
# ### At the inter level
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
# all_celltypes=names(table(Cor_analysis_DF$celltypes))
# upset_plot=list()
# for (i in 1:length(all_celltypes)) {
#   Subset_=subset(Cor_analysis_DF, celltypes==all_celltypes[i])
#   all_up=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   all_down=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_up=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_down=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_up=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_down=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   
#   # make the UpSet data
#   library(ComplexHeatmap)
#   library(ggplot2)
#   datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
#   m=make_comb_mat(datalist, mode="intersect")
#   m=normalize_comb_mat(m, full_comb_sets=T)
#   set_name(m) # check
#   combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
#   m=m[, combinations_]
#   percents_=c(comb_size(m)[1:3]/set_size(m)[1],
#               comb_size(m)[4:6]/set_size(m)[4])
#   percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
#   title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in ",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
#   
#   # plot
#   upset_plot_=
#     HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
#     UpSet(m, 
#           top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
#           right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
#           row_names_gp=gpar(fontsize=9),
#           set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
#           comb_order=1:6) %v% 
#     HeatmapAnnotation(
#       empty=anno_empty(border=F, height=unit(2,"mm")),
#       text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
#     )
#   upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Inter.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()
# 
# ### At the detailed level
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# all_celltypes=names(table(Cor_analysis_DF$celltypes))
# upset_plot=list()
# for (i in 1:length(all_celltypes)) {
#   Subset_=subset(Cor_analysis_DF, celltypes==all_celltypes[i])
#   all_up=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   all_down=Subset_ %>%
#     subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_up=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   female_down=Subset_ %>%
#     subset(analysis=="females" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_up=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho>0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
#   male_down=Subset_ %>%
#     subset(analysis=="males" & rho_pval<0.05 & rho<0) %>%
#     select(gene) %>%
#     unlist() %>% unname()
# 
#   # make the UpSet data
#   library(ComplexHeatmap)
#   library(ggplot2)
#   datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
#   m=make_comb_mat(datalist, mode="intersect")
#   m=normalize_comb_mat(m, full_comb_sets=T)
#   set_name(m) # check
#   combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
#   m=m[, combinations_]
#   percents_=c(comb_size(m)[1:3]/set_size(m)[1],
#               comb_size(m)[4:6]/set_size(m)[4])
#   percents_=sprintf("%0.2f%%", percents_*100) # calculate the percentages
#   title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in\n",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
# 
#   # plot
#   upset_plot_=
#     HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
#     UpSet(m,
#           top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
#           right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
#           row_names_gp=gpar(fontsize=9),
#           set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
#           comb_order=1:6) %v%
#     HeatmapAnnotation(
#       empty=anno_empty(border=F, height=unit(2,"mm")),
#       text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
#     )
#   upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Detailed.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()
# 
# ##################################
# 
# 
# 
# ### Volcano plots to show the sig. age-correlated genes in each celltype at the rough level
# # ...* The inter and detailed levels generate plots of too large size and limited information, so take only CD8T cells
# #####################################
# ###
# 
# library(dplyr)
# library(scRNAtoolVis)
# 
# ### For Rough level
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# # all
# df_all=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="all") %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_=
#   jjVolcano2(diffData=df_all, 
#             # myMarkers=shared_kinases,
#             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
#             aesCol=c("lightgrey","lightgrey"),
#             back.col="grey97",
#             log2FC.cutoff=0.1,
#             size=3,
#             tile.fontsize=3.5,
#             tile.linewidth=0.05,
#             pSize=0.5,
#             fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<both>")
# # females
# df_f=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="females") %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_f=
#   jjVolcano2(diffData=df_f, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
#              aesCol=c("lightpink","lightpink"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3,
#              tile.fontsize=3.5,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# # males
# df_m=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="males") %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_m=
#   jjVolcano2(diffData=df_m, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
#              aesCol=c("lightblue","lightblue"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3,
#              tile.fontsize=3.5,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# # save all/females/males
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Rough.pdf", width=length(unique(df_all$cluster)), height=5)
# plot(plot_); plot(plot_f); plot(plot_m)
# dev.off()
# 
# ### For Inter level CD8T cells
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
# # all
# df_all=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="all" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_=
#   jjVolcano2(diffData=df_all, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
#              aesCol=c("lightgrey","lightgrey"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3,
#              tile.fontsize=3.5,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<both>")
# # females
# df_f=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="females" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_f=
#   jjVolcano2(diffData=df_f, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
#              aesCol=c("lightpink","lightpink"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3,
#              tile.fontsize=3.5,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# # males
# df_m=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="males" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_m=
#   jjVolcano2(diffData=df_m, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
#              aesCol=c("lightblue","lightblue"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3,
#              tile.fontsize=3.5,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# # save all/females/males
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Inter_CD8T.pdf", width=5, height=5)
# plot(plot_); plot(plot_f); plot(plot_m)
# dev.off()
# 
# ### For Detailed level CD8T cells
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# # all
# df_all=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="all" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_=
#   jjVolcano2(diffData=df_all, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
#              aesCol=c("lightgrey","lightgrey"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3.5,
#              tile.fontsize=3,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<both>")
# # females
# df_f=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="females" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_f=
#   jjVolcano2(diffData=df_f, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
#              aesCol=c("lightpink","lightpink"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3.5,
#              tile.fontsize=3,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# # males
# df_m=Cor_analysis_DF %>% 
#   mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
#   subset(analysis=="males" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
#   mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
#   select(p_val, avg_log2FC, p_val_adj, cluster, gene)
# source("~/Rscripts/jjVolcano2.R")
# plot_m=
#   jjVolcano2(diffData=df_m, 
#              # myMarkers=shared_kinases,
#              tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
#              aesCol=c("lightblue","lightblue"),
#              back.col="grey97",
#              log2FC.cutoff=0.1,
#              size=3.5,
#              tile.fontsize=3,
#              tile.linewidth=0.05,
#              pSize=0.5,
#              fontface="italic") +
#   guides(color=guide_legend(size=3)) +
#   theme(title=element_text(size=10), 
#         legend.text=element_text(size=12),
#         legend.position="none",
#         axis.title.y=element_text(size=12)) +
#   labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# # save all/females/males
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Detailed_CD8T.pdf", width=20, height=5)
# plot(plot_); plot(plot_f); plot(plot_m)
# dev.off()
# 
# ##################################
# 
# 
# 
# ### Plot the correlation between the top30 age-correlated genes and ages in the celltypes at the rough level
# #####################################
# ###
# 
# pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# celltypes_all=unique(Cor_analysis_DF$celltypes)
# cor_plots_neg=cor_plots_pos=list()
# for (i in 1:length(celltypes_all)) {
#   ### Plot downregulated genes
#   Cor_df_neg=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
#     slice_min(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_neg[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Plot upregulated genes
#   Cor_df_pos=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
#     slice_max(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_pos[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Save all
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
#   plot(cor_plots_neg[[i]])
#   dev.off()
#   
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
#   plot(cor_plots_pos[[i]])
#   dev.off()
# }
# qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_neg.pdf"),
#                           paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_pos.pdf")),
#                   output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough.pdf")
# unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_neg.pdf"),
#          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_pos.pdf")))
# 
# ##################################
# 
# 
# 
# ### Plot the correlation between the top30 age-correlated genes and ages in the CD8T cells at the inter level
# #####################################
# ###
# 
# pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
# celltypes_all=unique(Cor_analysis_DF$celltypes)
# celltypes_all=celltypes_all[grepl("^CD8T\\.",celltypes_all)]
# cor_plots_neg=cor_plots_pos=list()
# for (i in 1:length(celltypes_all)) {
#   ### Plot downregulated genes
#   Cor_df_neg=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
#     slice_min(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_neg[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Plot upregulated genes
#   Cor_df_pos=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
#     slice_max(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_pos[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Save all
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
#   plot(cor_plots_neg[[i]])
#   dev.off()
#   
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
#   plot(cor_plots_pos[[i]])
#   dev.off()
# }
# qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_neg.pdf"),
#                           paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_pos.pdf")),
#                   output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_CD8T.pdf")
# unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_neg.pdf"),
#          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_pos.pdf")))
# 
# ##################################
# 
# 
# 
# ### Plot the correlation between the top30 age-correlated genes and ages in the CD8T cells at the detailed level
# #####################################
# ###
# 
# pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# celltypes_all=unique(Cor_analysis_DF$celltypes)
# celltypes_all=celltypes_all[grepl("^CD8T\\.",celltypes_all)]
# cor_plots_neg=cor_plots_pos=list()
# for (i in 1:length(celltypes_all)) {
#   ### Plot downregulated genes
#   Cor_df_neg=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
#     slice_min(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_neg[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Plot upregulated genes
#   Cor_df_pos=Cor_analysis_DF %>%
#     subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
#     slice_max(rho, n=100) %>%
#     select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
#   Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
#   # get the expr data
#   pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes_all[i])
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M")) %>%
#     tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
#   # plot
#   cor_plots_pos[[i]]=
#     ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93", alpha=0.5) +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5, method="spearman")
#   
#   ### Save all
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
#   plot(cor_plots_neg[[i]])
#   dev.off()
#   
#   pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
#   plot(cor_plots_pos[[i]])
#   dev.off()
# }
# qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_neg.pdf"),
#                           paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_pos.pdf")),
#                   output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_CD8T.pdf")
# unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_neg.pdf"),
#          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_pos.pdf")))
# 
# ##################################
# 
# 
# 
# ### Map the reported/known age-associated genes in the Pval-Cor plot
# #####################################
# ###
# 
# ### Take only the aging-related genes based on DEMAGALHAES_AGING_UP Dataset
# AgingDB=read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")
# marker_genes=AgingDB$gene
# 
# ### Load the Correlation df
# Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# 
# ### Plot the Pval-Cor volcano plots for each celltype
# celltypes_all=unique(Cor_analysis_DF$celltypes)
# plot_corr=list()
# for (i in 1:length(celltypes_all)) {
#   Cor_analysis_subset=Cor_analysis_DF %>% 
#     subset(celltypes==celltypes_all[i]) %>% 
#     mutate(mark=ifelse(gene %in% marker_genes, gene, NA),
#            analysis=ifelse(analysis=="all","<both>",ifelse(analysis=="females","<females>","<males>")))
#   plot_corr[[i]]=
#     ggplot(Cor_analysis_subset, aes(x=rho, y=-log10(rho_pval), color=analysis)) +
#     ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.5) +
#     labs(title=paste0("Gene-Age correlation in ",celltypes_all[i]), y=expression(-log[10]~pval), x=expression(Spearman~rho)) +
#     facet_wrap(~analysis) +
#     theme_classic() +
#     theme(title=element_text(size=11),
#           plot.subtitle=element_text(size=11, face="italic"),
#           axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.y=element_text(size=9),
#           axis.title.x=element_text(size=9),
#           legend.position="none",
#           legend.title=element_text(size=10),
#           legend.text=element_text(size=9),
#           strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
#           strip.text.x=element_text(size=10, color="black", face="italic")) +
#     scale_color_manual(values=c("grey","lightpink","lightblue")) +
#     geom_text_repel(aes(label=mark), size=3, show.legend=FALSE, color="black")
# }
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Rough_MarkedwKnownGenes.pdf", height=5, width=6)
# for(i in 1:length(plot_corr)) plot(plot_corr[[i]])
# dev.off()
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
# 
# 
# 
# RNA.sig.cor_genes_FROM_DMP.neg.cor_genes=RNA.sig.cor_genes_FROM_DMP.pos.cor_genes=list()
# neg_plot=pos_plot=list()
# for (i in 1:7) {
#   pseudo_subset=subset(pseudobulkobj, Annot.rough %in% rough_celltype[[i]])
# 
# 
#   ### Analyze the correlation between expr of the hypomethylated genes and ages
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", as.character(Neg.Cor.genes[[i]])), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M"))
#   colsum_=colSums(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")])
#   filter=c("age", "sex", "donor_id", names(colsum_)[colsum_!=0])
#   expr_data=expr_data[, filter]
#   genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
#   age_correlated_expr_genes=c()
#   for (j in 1:length(genes_)) {
#     cor_result=cor.test(expr_data$age, expr_data[[genes_[j]]])
#     if (cor_result$p.value<0.05) {
#       the_cor=cor_result$estimate
#       names(the_cor)=genes_[j]
#       age_correlated_expr_genes=c(age_correlated_expr_genes, the_cor)
#     }
#   }
#   # save all the (RNA-expr)sig.age-cor.genes from the hypomethylated genes and their correlation with ages
#   RNA.sig.cor_genes_FROM_DMP.neg.cor_genes[[i]]=
#     data.frame(gene=names(age_correlated_expr_genes),
#                cor=age_correlated_expr_genes,
#                DMP_type="Significantly hypomethylated",
#                celltype=tissues_selected[i])
# 
#   # make the df with the top30 abs(cor) genes
#   age_correlated_expr_genes=names(sort(abs(age_correlated_expr_genes), decreasing=T)[1:30])
#   age_correlated_expr_genes=age_correlated_expr_genes[!is.na(age_correlated_expr_genes)]
#   expr_data=expr_data %>%
#     select(c("age", "sex", "donor_id", age_correlated_expr_genes))
#   expr_data_longDF=tidyr::pivot_longer(expr_data, cols=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex","donor_id")], names_to="gene", values_to="expression")
# 
#   # plot
#   neg_plot[[i]]=
#     ggplot(expr_data_longDF, aes(x=age, y=expression, color=sex)) +
#     geom_point(size=0.1, shape=20, alpha=0.05) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Expression of the hypomethylated genes in ", tissues_selected[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5)
# 
#   pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",i,"_neg.pdf"), height=ceiling(length(age_correlated_expr_genes)/3)*2, width=9)
#   plot(neg_plot[[i]])
#   dev.off()
# 
# 
#   ### Analyze the correlation between expr of the hypermethylated genes and ages
#   expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", as.character(Pos.Cor.genes[[i]])), layer="data")
#   expr_data=expr_data %>%
#     mutate(sex=ifelse(sex=="female", "F", "M"))
#   colsum_=colSums(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")])
#   filter=c("age", "sex", "donor_id", names(colsum_)[colsum_!=0])
#   expr_data=expr_data[, filter]
#   genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
#   age_correlated_expr_genes=c()
#   for (j in 1:length(genes_)) {
#     cor_result=cor.test(expr_data$age, expr_data[[genes_[j]]])
#     if (cor_result$p.value<0.05) {
#       the_cor=cor_result$estimate
#       names(the_cor)=genes_[j]
#       age_correlated_expr_genes=c(age_correlated_expr_genes, the_cor)
#     }
#   }
#   # save all the (RNA-expr)sig.age-cor.genes from the hypermethylated genes and their correlation with ages
#   RNA.sig.cor_genes_FROM_DMP.pos.cor_genes[[i]]=
#     data.frame(gene=names(age_correlated_expr_genes),
#                cor=age_correlated_expr_genes,
#                DMP_type="Significantly hypermethylated",
#                celltype=tissues_selected[i])
# 
#   # make the df with the top30 abs(cor) genes
#   age_correlated_expr_genes=names(sort(abs(age_correlated_expr_genes), decreasing=T)[1:30])
#   age_correlated_expr_genes=age_correlated_expr_genes[!is.na(age_correlated_expr_genes)]
#   expr_data=expr_data %>%
#     select(c("age", "sex", "donor_id", age_correlated_expr_genes))
#   expr_data_longDF=tidyr::pivot_longer(expr_data, cols=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex","donor_id")], names_to="gene", values_to="expression")
# 
#   # plot
#   pos_plot[[i]]=
#     ggplot(expr_data_longDF, aes(x=age, y=expression, color=sex)) +
#     geom_point(size=0.1, shape=20, alpha=0.05) +
#     # ggsci::scale_color_d3() +
#     facet_wrap(~gene, ncol=3, scales="free") +
#     theme_classic() +
#     labs(x="age", y="Expression", title=paste0("Expression of the hypermethylated genes in ", tissues_selected[i])) +
#     theme(axis.text.x=element_text(size=8),
#           axis.text.y=element_text(size=8),
#           axis.title.x=element_text(size=9),
#           axis.title.y=element_text(size=9),
#           # legend.position="none",
#           # legend.title=element_blank(),
#           title=element_text(size=10)) +
#     guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
#     geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
#     ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5)
# 
#   pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",i,"_pos.pdf"), height=ceiling(length(age_correlated_expr_genes)/3)*2, width=9)
#   plot(pos_plot[[i]])
#   dev.off()
# 
#   message(paste0("Plotting of the expr of the RNA.sig.age-cor.genes in ",tissues_selected[i], " is done!"))
# }
# 
# ### Save all the (RNA-expr)sig.age-cor.genes from the (DMP-beta)sig.age-cor.genes and their correlation with ages
# RNA.sig.cor_genes_WITH_cor=data.table::rbindlist(c(RNA.sig.cor_genes_FROM_DMP.neg.cor_genes, RNA.sig.cor_genes_FROM_DMP.pos.cor_genes))
# saveRDS(RNA.sig.cor_genes_WITH_cor, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/RelatedGenes_RNAexpression_ageCorrelationResult.rds")
# 
# 
# 
# 
