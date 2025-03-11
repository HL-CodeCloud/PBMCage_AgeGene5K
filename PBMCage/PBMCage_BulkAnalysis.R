

### Cosine similarity between gene expr and ages in each celltype at the rough level
#####################################
###

library(dplyr)
library(Seurat)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")

### Calculate the cos<gene-expr, age> in both sex
celltypes=names(table(pseudobulkobj$Annot.rough))

cosine_df=list()
for (i in 1:length(celltypes)) {
  pseudobulkobj_celltype=subset(pseudobulkobj, Annot.rough==celltypes[i])
  
  pseudo_data=LayerData(subset(pseudobulkobj_celltype), layer="data")
  ages=as.numeric(unlist(lapply(strsplit(colnames(pseudo_data), split="_"), function(x) x[3])))
  pseudo_data.t=data.table::transpose(as.data.frame(pseudo_data))
  colnames(pseudo_data.t)=rownames(pseudo_data)
  rownames(pseudo_data.t)=colnames(pseudo_data)
  cosine_value=apply(pseudo_data.t, 2, function(x) lsa::cosine(x, ages))
  cosine_df[[i]]=data.frame(gene=names(cosine_value), cosine=cosine_value, analysis="all", celltypes=celltypes[i])
}
Cosine_DF_all=data.table::rbindlist(cosine_df)

### Calculate the cos<gene-expr, age> in females
pseudobulkobj_females=subset(pseudobulkobj, sex=="female")
celltypes=names(table(pseudobulkobj_females$Annot.rough))

cosine_df=list()
for (i in 1:length(celltypes)) {
  pseudobulkobj_celltype=subset(pseudobulkobj_females, Annot.rough==celltypes[i])
  
  pseudo_data=LayerData(subset(pseudobulkobj_celltype), layer="data")
  ages=as.numeric(unlist(lapply(strsplit(colnames(pseudo_data), split="_"), function(x) x[3])))
  pseudo_data.t=data.table::transpose(as.data.frame(pseudo_data))
  colnames(pseudo_data.t)=rownames(pseudo_data)
  rownames(pseudo_data.t)=colnames(pseudo_data)
  cosine_value=apply(pseudo_data.t, 2, function(x) lsa::cosine(x, ages))
  cosine_df[[i]]=data.frame(gene=names(cosine_value), cosine=cosine_value, analysis="females", celltypes=celltypes[i])
}
Cosine_DF_females=data.table::rbindlist(cosine_df)

### Calculate the cos<gene-expr, age> in males
pseudobulkobj_males=subset(pseudobulkobj, sex=="male")
celltypes=names(table(pseudobulkobj_males$Annot.rough))

cosine_df=list()
for (i in 1:length(celltypes)) {
  pseudobulkobj_celltype=subset(pseudobulkobj_males, Annot.rough==celltypes[i])
  
  pseudo_data=LayerData(subset(pseudobulkobj_celltype), layer="data")
  ages=as.numeric(unlist(lapply(strsplit(colnames(pseudo_data), split="_"), function(x) x[3])))
  pseudo_data.t=data.table::transpose(as.data.frame(pseudo_data))
  colnames(pseudo_data.t)=rownames(pseudo_data)
  rownames(pseudo_data.t)=colnames(pseudo_data)
  cosine_value=apply(pseudo_data.t, 2, function(x) lsa::cosine(x, ages))
  cosine_df[[i]]=data.frame(gene=names(cosine_value), cosine=cosine_value, analysis="males", celltypes=celltypes[i])
}
Cosine_DF_males=data.table::rbindlist(cosine_df)

Cosine_DF_everything=data.table::rbindlist(list(Cosine_DF_all, Cosine_DF_females, Cosine_DF_males))
data.table::fwrite(Cosine_DF_everything, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_CosineSimilarity_Rough.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in all celltypes
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj=AggregateExpression(pseudobulkobj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
pseudobulkobj$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[1])))
pseudobulkobj$sex=unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[2]))
pseudobulkobj$donor_id=unlist(lapply(strsplit(colnames(pseudobulkobj), split="_"), function(x) x[3]))
all_genes=rownames(pseudobulkobj)

pseudo_subset=pseudobulkobj

# get all the expr data
expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
expr_data=expr_data %>%
  mutate(sex=ifelse(sex=="female", "F", "M"))
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
                        analysis="all", celltypes="all")
}
if (length(genes_f)>=1) {
  Female_COR_DF=data.frame(gene=genes_f, 
                           cor=female_cor, cor_pval=female_cor_pval, rho=female_rho, rho_pval=female_rho_pval, tau=female_tau, tau_pval=female_tau_pval, 
                           analysis="females", celltypes="all")
}
if (length(genes_m)>=1) {
  Male_COR_DF=data.frame(gene=genes_m, 
                         cor=male_cor, cor_pval=male_cor_pval, rho=male_rho, rho_pval=male_rho_pval, tau=male_tau, tau_pval=male_tau_pval, 
                         analysis="males", celltypes="all")
}

Cor_analysis_DF_list=data.table::rbindlist(list(ALL_COR_DF, Female_COR_DF, Male_COR_DF))

data.table::fwrite(Cor_analysis_DF_list, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_donorid.csv.gz", sep="\t")

#####################################



### Analyze the overlapping and difference in the age-correlated genes at the donorid level between female and males
#####################################
###

library(dplyr)
library(ggplot2)

### Load the data
donor_cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_donorid.csv.gz", sep="\t")
donor_cor_female=donor_cor %>% subset(rho_pval<0.05 & analysis=="females")
donor_cor_male=donor_cor %>% subset(rho_pval<0.05 & analysis=="males")

### Take the genes by sex
female_only=(donor_cor_female %>% subset(! gene %in% donor_cor_male$gene))$gene
male_only=(donor_cor_male %>% subset(! gene %in% donor_cor_female$gene))$gene
shared=(donor_cor_female %>% subset(gene %in% donor_cor_male$gene))$gene

### Plot
df_for_plot=data.table::rbindlist(list(donor_cor %>% subset(gene %in% shared & analysis=="all") %>% mutate(sex="both"),
                                       donor_cor %>% subset(gene %in% female_only & analysis=="all") %>% mutate(sex="F"),
                                       donor_cor %>% subset(gene %in% male_only & analysis=="all") %>% mutate(sex="M")))

plot_=
  ggplot(subset(df_for_plot, sex!="both"), aes(x=rho, color=sex)) +
  geom_density() +
  coord_flip() +
  theme_classic() +
  scale_color_manual(values=c("coral3","deepskyblue3")) +
  labs(x=expression(Spearman~rho), y=NULL, title=NULL) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line.x=element_blank(),
        axis.title.y=element_text(size=10),
        legend.position="top",
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_density_FvsM.pdf", height=4, width=2)
plot(plot_)
dev.off()

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the rough level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.rough))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
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
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the rough level, <age81
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj=subset(pseudobulkobj, age<81)
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.rough))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
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
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough_LowerThanAge81.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the inter level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.inter))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
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
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the detailed level
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.detailed))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
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
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")

#####################################



### Analyze the correlation between the RNA_expr and ages in each celltype at the detailed level,
### ...based on the pseudobulk obj after being updated for publication
#####################################
###

library(dplyr)

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed_updated.rds")
all_genes=rownames(pseudobulkobj)
celltypes=names(table(pseudobulkobj$Annot.detailed))

Cor_analysis_DF_list=list()
for (i in 1:length(celltypes)) {
  pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes[i])
  
  # get all the expr data
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", all_genes), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
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
data.table::fwrite(Cor_analysis_DF, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed_updated.csv.gz", sep="\t")

#####################################



### Analyze the overlapping and difference in the age-correlated genes in reproducible freq-changing detailed celltypes between female and males
#####################################
###

library(dplyr)
library(ggplot2)

reproducable_celltypes=c("CD4T.naive.COTL1pos_SOX4pos","CD4T.Tcm","CD4T.Tcm.HLA","CD4T.Treg.FOXP3hi","CD8T.naive.LINC02446","CD8T.Tcm.GATA3","CD8T.Tem.GZMK_XCL1","NK.CD56dim.FCER1Glo_KLRC2pos")

### Load the data
donor_cor=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
donor_cor_female=donor_cor %>% subset(rho_pval<0.05 & analysis=="females" & celltypes %in% reproducable_celltypes)
donor_cor_male=donor_cor %>% subset(rho_pval<0.05 & analysis=="males" & celltypes %in% reproducable_celltypes)
donor_cor_both=donor_cor %>% subset(celltypes %in% reproducable_celltypes)

# plot_=
  ggplot(subset(donor_cor_both, analysis!="all"), aes(x=rho, color=analysis)) +
  facet_wrap(~celltypes, scales="fixed") +
  geom_density(aes(y=after_stat(density)-15), linewidth=0.5) +
  geom_density(inherit.aes=F, aes(x=rho_pval, color=analysis, y=15-after_stat(density)), alpha=0.5, linewidth=0.25) +
  geom_vline(xintercept=0, color="grey45", linewidth=0.25, linetype="dotted") +
  geom_vline(xintercept=c(0.25,-0.25), color="grey70", linewidth=0.25, linetype="dotted") +
  theme_classic() +
  scale_color_manual(values=c("coral3","deepskyblue3")) +
  labs(x=expression(Spearman~rho), y=NULL, title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        axis.line.y=element_blank(),
        axis.title.x=element_text(size=10),
        legend.position="top",
        title=element_text(size=10)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
  scale_x_continuous(sec.axis=dup_axis(name="pval"), limits=c(-1,1))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_density_FvsM_reproducableDetailedCelltypes.pdf", height=4, width=2)
plot(plot_)
dev.off()

#####################################



### Upset analysis of the sig.cor.genes in all/females/males in each celltype
#####################################
###

library(dplyr)

### At the rough level
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
all_celltypes=names(table(Cor_analysis_DF$celltypes))
upset_plot=list()
PERCENTs_rough=list()
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
  datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[1:3]/set_size(m)[1],
              comb_size(m)[4:6]/set_size(m)[4])
  # save the precentages
  PERCENTs_rough[[i]]=c(percents_, all_celltypes[i])
  percents_=sprintf("%0.2f%%", percents_*100) # for plotting
  
  # title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in ",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
  # 
  # # plot
  # upset_plot_=
  #   HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
  #   UpSet(m,
  #         top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
  #         right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
  #         row_names_gp=gpar(fontsize=9),
  #         set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
  #         comb_order=1:6) %v%
  #   HeatmapAnnotation(
  #     empty=anno_empty(border=F, height=unit(2,"mm")),
  #     text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
  #   )
  # upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Rough.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()

### At the inter level
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
all_celltypes=names(table(Cor_analysis_DF$celltypes))
upset_plot=list()
PERCENTs_inter=list()
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
  datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[1:3]/set_size(m)[1],
              comb_size(m)[4:6]/set_size(m)[4])
  # save the precentages
  PERCENTs_inter[[i]]=c(percents_, all_celltypes[i])
  percents_=sprintf("%0.2f%%", percents_*100) # for plotting
  
  # title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in ",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
  # 
  # # plot
  # upset_plot_=
  #   HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
  #   UpSet(m,
  #         top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
  #         right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
  #         row_names_gp=gpar(fontsize=9),
  #         set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
  #         comb_order=1:6) %v%
  #   HeatmapAnnotation(
  #     empty=anno_empty(border=F, height=unit(2,"mm")),
  #     text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
  #   )
  # upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Inter.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()

### At the detailed level
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
all_celltypes=names(table(Cor_analysis_DF$celltypes))
upset_plot=list()
PERCENTs_detailed=list()
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
  datalist=list(Neg.both=all_down, Neg.F=female_down, Neg.M=male_down, Pos.both=all_up, Pos.F=female_up, Pos.M=male_up)
  m=make_comb_mat(datalist, mode="intersect")
  m=normalize_comb_mat(m, full_comb_sets=T)
  set_name(m) # check
  combinations_=c("110000","101000","011000","000110","000101","000011") # set the combinations & order
  m=m[, combinations_]
  percents_=c(comb_size(m)[1:3]/set_size(m)[1],
              comb_size(m)[4:6]/set_size(m)[4])
  # save the precentages
  PERCENTs_detailed[[i]]=c(percents_, all_celltypes[i])
  percents_=sprintf("%0.2f%%", percents_*100) # for plotting
  
  # title_fun=function(ht) {decorate_annotation("title", grid.text(paste0("Age-correlated genes in\n",all_celltypes[i]), x=unit(0, "mm"), just="left"))}
  # 
  # # plot
  # upset_plot_=
  #   HeatmapAnnotation(title=anno_empty(border=FALSE)) %v%
  #   UpSet(m,
  #         top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
  #         right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
  #         row_names_gp=gpar(fontsize=9),
  #         set_order=c("Neg.both","Neg.F","Neg.M","Pos.both","Pos.F","Pos.M"),
  #         comb_order=1:6) %v%
  #   HeatmapAnnotation(
  #     empty=anno_empty(border=F, height=unit(2,"mm")),
  #     text=anno_text(percents_, rot=0, just="center", gp=gpar(fontsize=9))
  #   )
  # upset_plot[[i]]=draw(upset_plot_, post_fun=title_fun)
}
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_Detailed.pdf", width=5, height=4)
# for(i in 1:length(upset_plot)) plot(upset_plot[[i]])
# dev.off()

### Save all the precentages for plotting the females vs. males
names_precents=c("(Neg.both & Neg.F)/Neg.both", "(Neg.both & Neg.M)/Neg.both","(Neg.F & Neg.M)/Neg.both",
                 "(Pos.both & Pos.F)/Pos.both","(Pos.both & Pos.M)/Pos.both","(Pos.F & Pos.M)/Pos.both")
rough_df=data.frame(c1=unlist(lapply(PERCENTs_rough, function(x) x[1])),
                    c2=unlist(lapply(PERCENTs_rough, function(x) x[2])),
                    c3=unlist(lapply(PERCENTs_rough, function(x) x[3])),
                    c4=unlist(lapply(PERCENTs_rough, function(x) x[4])),
                    c5=unlist(lapply(PERCENTs_rough, function(x) x[5])),
                    c6=unlist(lapply(PERCENTs_rough, function(x) x[6])),
                    celltype=unlist(lapply(PERCENTs_rough, function(x) x[7])),
                    celltype.level="Annot.rough")
inter_df=data.frame(c1=unlist(lapply(PERCENTs_inter, function(x) x[1])),
                    c2=unlist(lapply(PERCENTs_inter, function(x) x[2])),
                    c3=unlist(lapply(PERCENTs_inter, function(x) x[3])),
                    c4=unlist(lapply(PERCENTs_inter, function(x) x[4])),
                    c5=unlist(lapply(PERCENTs_inter, function(x) x[5])),
                    c6=unlist(lapply(PERCENTs_inter, function(x) x[6])),
                    celltype=unlist(lapply(PERCENTs_inter, function(x) x[7])),
                    celltype.level="Annot.inter")
detailed_df=data.frame(c1=unlist(lapply(PERCENTs_detailed, function(x) x[1])),
                       c2=unlist(lapply(PERCENTs_detailed, function(x) x[2])),
                       c3=unlist(lapply(PERCENTs_detailed, function(x) x[3])),
                       c4=unlist(lapply(PERCENTs_detailed, function(x) x[4])),
                       c5=unlist(lapply(PERCENTs_detailed, function(x) x[5])),
                       c6=unlist(lapply(PERCENTs_detailed, function(x) x[6])),
                       celltype=unlist(lapply(PERCENTs_detailed, function(x) x[7])),
                       celltype.level="Annot.detailed")
DFs=data.table::rbindlist(list(rough_df, inter_df, detailed_df))
colnames(DFs)=c(names_precents, "celltype","celltype.level")
write.table(DFs, "~/Project_PBMCage/Results/PBMC_results/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_percents.txt", sep="\t")

#####################################



### Analyze the overlapping of age-cor genes in CD4T and CD8T Annot.inter
#####################################
###

###
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")

gene_subset=Cor_analysis_DF %>% 
  subset(analysis=="all" & rho_pval<0.05) %>%
  subset(grepl("CD4|CD8",celltypes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Calculate the jaccard similarity for pos genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(gene_subset)) {
  the_first_celltype=gene_subset[[i]]
  
  for (j in 1:length(gene_subset)) {
    the_second_celltype=gene_subset[[j]]
    
    jaccard_similarity=jaccard(the_first_celltype, the_second_celltype)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(gene_subset)[i],"_vs_",names(gene_subset)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name

library(ComplexHeatmap)
Jaccard_SIM_df=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)
Jaccard_SIM_df=Jaccard_SIM_df %>% mutate(celltype=gsub("_vs_.*","",pair), module=gsub(".*_vs_","",pair)) %>%
  select(-pair) %>%
  tidyr::pivot_wider(names_from="module", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype")
# reorder the matrix
reorder_=colnames(Jaccard_SIM_df); reorder_
reorder_=reorder_[c(1,3,4,5,2,7,8,9,6)]; reorder_
Jaccard_SIM_df=Jaccard_SIM_df[reorder_,reorder_]

range(Jaccard_SIM_df[Jaccard_SIM_df!=1]) # check

### Plot
col_fun=circlize::colorRamp2(c(0, 0.5), c("white", "brown3"))
ht=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F, show_column_names=F,
          row_names_gp=gpar(fontsize=9),
          row_names_side="right",
          row_names_rot=0,
          column_names_gp=gpar(fontsize=9),
                  layer_fun=function(j, i, x, y, width, height, fill) {
                    grid.text(sprintf("%.2f",as.matrix(Jaccard_SIM_df)), x, y,
                              gp=gpar(fontsize=9, fontface="plain"))},
          height=unit(8,"mm")*nrow(Jaccard_SIM_df), 
          width=unit(8,"mm")*ncol(Jaccard_SIM_df))
lgd=Legend(col_fun=col_fun, 
           title="Jaccard similarity index", title_gp=gpar(fontface="plain"),
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_opt(RESET=T)
ht_opt$ANNOTATION_LEGEND_PADDING=unit(0.75,"cm")
ht_final=draw(ht, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_JaccardSimilarity_CD4T.CD8T_Inter.pdf", height=3.5, width=4.5)
ht_final
dev.off()

#####################################



### Analyze the overlapping of age-cor genes in CD4T and CD8T Annot.inter
#####################################
###

###
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")

gene_subset=Cor_analysis_DF %>% 
  subset(analysis=="all" & rho_pval<0.05) %>%
  subset(grepl("CD4|CD8",celltypes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Calculate the jaccard similarity for pos genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(gene_subset)) {
  the_first_celltype=gene_subset[[i]]
  
  for (j in 1:length(gene_subset)) {
    the_second_celltype=gene_subset[[j]]
    
    jaccard_similarity=jaccard(the_first_celltype, the_second_celltype)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(gene_subset)[i],"_vs_",names(gene_subset)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name

library(ComplexHeatmap)
Jaccard_SIM_df=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)
Jaccard_SIM_df=Jaccard_SIM_df %>% mutate(celltype=gsub("_vs_.*","",pair), module=gsub(".*_vs_","",pair)) %>%
  select(-pair) %>%
  tidyr::pivot_wider(names_from="module", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype")
# # reorder the matrix
# reorder_=colnames(Jaccard_SIM_df); reorder_
# reorder_=reorder_[c(1,3,4,5,2,7,8,9,6)]; reorder_
# Jaccard_SIM_df=Jaccard_SIM_df[reorder_,reorder_]

range(Jaccard_SIM_df[Jaccard_SIM_df!=1]) # check

### Plot
col_fun=circlize::colorRamp2(c(0, 1), c("white", "brown3"))
# ht=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F, show_column_names=F,
          row_names_gp=gpar(fontsize=9),
          row_names_side="right",
          row_names_rot=0,
          column_names_gp=gpar(fontsize=9),
          layer_fun=function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f",as.matrix(Jaccard_SIM_df)), x, y,
                      gp=gpar(fontsize=9, fontface="plain"))},
          height=unit(8,"mm")*nrow(Jaccard_SIM_df), 
          width=unit(8,"mm")*ncol(Jaccard_SIM_df))
lgd=Legend(col_fun=col_fun, 
           title="Jaccard similarity index", title_gp=gpar(fontface="plain"),
           legend_height=unit(6, "cm"), grid_width=unit(0.2, "cm"),
           direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
)
ht_opt(RESET=T)
ht_opt$ANNOTATION_LEGEND_PADDING=unit(0.75,"cm")
ht_final=draw(ht, annotation_legend_list=lgd, annotation_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_JaccardSimilarity_CD4T.CD8T_Inter.pdf", height=3.5, width=4.5)
ht_final
dev.off()











### Plot the percentages of overlapping between females and males in the Upset analysis of the sig.cor.genes in each celltype
#####################################
###

### Load the precent results
DFs=read.delim("~/Project_PBMCage/Results/PBMC_results/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_percents.txt", sep="\t")
names_cols=c("(Neg.both & Neg.F)/Neg.both", "(Neg.both & Neg.M)/Neg.both","(Neg.F & Neg.M)/Neg.both",
             "(Pos.both & Pos.F)/Pos.both","(Pos.both & Pos.M)/Pos.both","(Pos.F & Pos.M)/Pos.both",
             "celltype","celltype.level")
colnames(DFs)=names_cols
DFs=DFs %>% tidyr::pivot_longer(cols=names_cols[1:6], names_to="comparison", values_to="percent")

### Prepare the Fvs.M data for plotting
DFs_FvsM=subset(DFs, comparison %in% c("(Neg.F & Neg.M)/Neg.both","(Pos.F & Pos.M)/Pos.both")) %>%
  mutate(orientation=ifelse(comparison=="(Pos.F & Pos.M)/Pos.both", "upreg.", "downreg."))
DFs_FvsM$orientation=forcats::fct_relevel(DFs_FvsM$orientation, "upreg.","downreg.")
DFs_FvsM$percent=100*(DFs_FvsM$percent)

### Plot
library(ggplot2)
# for rough
plot_rough=
  ggplot(subset(DFs_FvsM, celltype.level=="Annot.rough"), aes(x=celltype, y=percent, color=orientation)) +
  geom_linerange(data=subset(DFs_FvsM, celltype.level=="Annot.rough" & orientation=="upreg."),
                 aes(xmin=celltype, xmax=celltype, ymin=0, ymax=percent), size=1, alpha=0.9) +
  geom_linerange(data=subset(DFs_FvsM, celltype.level=="Annot.rough" & orientation=="downreg."),
                 aes(xmin=celltype, xmax=celltype, ymin=0, ymax=-percent), size=1, alpha=0.9) +
  scale_color_manual(values=c("lightskyblue3","coral3")) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(breaks=seq(0,100,20)) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_FvsM_overlap_percents_rough.pdf", height=3.5, width=4.5)
plot(plot_rough)
dev.off()

# for inter
plot_inter=
  ggplot(subset(DFs_FvsM, celltype.level=="Annot.inter"), aes(x=celltype, y=percent, color=orientation)) +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.inter" & orientation=="upreg."),
               aes(x=celltype, xend=celltype, y=0, yend=percent), size=1, alpha=0.9, color="coral3") +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.inter" & orientation=="downreg."),
               aes(x=celltype, xend=celltype, y=0, yend=-percent), size=1, alpha=0.9, color="lightskyblue3") +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  geom_line(alpha=0) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=8, angle=-90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(breaks=seq(0,100,20)) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_FvsM_overlap_percents_2.pdf", height=5, width=5)
plot(plot_inter)
dev.off()

# for detailed
DFs_FvsM_switched=DFs_FvsM
DFs_FvsM_switched$orientation=forcats::fct_relevel(DFs_FvsM_switched$orientation, c("downreg.","upreg."))
plot_detailed=
  ggplot(subset(DFs_FvsM_switched, celltype.level=="Annot.detailed"), aes(x=celltype, y=percent, color=orientation)) +
  geom_segment(data=subset(DFs_FvsM_switched, celltype.level=="Annot.detailed" & orientation=="upreg."),
               aes(x=celltype, xend=celltype, y=0, yend=percent), size=1, alpha=0.9, color="coral3") +
  geom_segment(data=subset(DFs_FvsM_switched, celltype.level=="Annot.detailed" & orientation=="downreg."),
               aes(x=celltype, xend=celltype, y=0, yend=-percent), size=1, alpha=0.9, color="lightskyblue3") +
  geom_line(alpha=0) +
  scale_color_manual(values=c("lightskyblue3","coral3")) +
  coord_flip() +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(labels=abs) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1, size=10)))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_FvsM_overlap_percents_detailed.pdf", height=16, width=4.5)
plot(plot_detailed)
dev.off()

# for CD8T-inter
CD8T_celltypes=DFs_FvsM$celltype[grepl("^CD8T\\.",DFs_FvsM$celltype)]
plot_inter_CD8=
  ggplot(subset(DFs_FvsM, celltype.level=="Annot.inter" & celltype %in% CD8T_celltypes), aes(x=celltype, y=percent, color=orientation)) +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.inter" & celltype %in% CD8T_celltypes & orientation=="upreg."),
               aes(x=celltype, xend=celltype, y=0, yend=percent), size=1, alpha=0.9, color="coral3") +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.inter" & celltype %in% CD8T_celltypes & orientation=="downreg."),
               aes(x=celltype, xend=celltype, y=0, yend=-percent), size=1, alpha=0.9, color="lightskyblue3") +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  geom_line(alpha=0) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(breaks=seq(0,100,20)) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1)))

# for CD8T-detailed
CD8T_celltypes=DFs_FvsM$celltype[grepl("^CD8T\\.",DFs_FvsM$celltype)]
plot_detailed_CD8=
  ggplot(subset(DFs_FvsM, celltype.level=="Annot.detailed" & celltype %in% CD8T_celltypes), aes(x=celltype, y=percent, color=orientation)) +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.detailed" & celltype %in% CD8T_celltypes & orientation=="upreg."),
               aes(x=celltype, xend=celltype, y=0, yend=percent), size=1, alpha=0.9, color="coral3") +
  geom_segment(data=subset(DFs_FvsM, celltype.level=="Annot.detailed" & celltype %in% CD8T_celltypes & orientation=="downreg."),
               aes(x=celltype, xend=celltype, y=0, yend=-percent), size=1, alpha=0.9, color="lightskyblue3") +
  scale_color_manual(values=c("coral3","lightskyblue3")) +
  geom_line(alpha=0) +
  theme_minimal() +
  labs(x=NULL, y="Precentage of overlap (%)") +
  theme(axis.text.x=element_text(size=8, angle=-90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=8),
        axis.title.y=element_text(size=9),
        title=element_text(size=10),
        legend.position="top") +
  scale_y_continuous(breaks=seq(0,100,20)) +
  guides(color=guide_legend(title=NULL, override.aes=list(linewidth=0.8, alpha=1)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_UpSet_FvsM_overlap_percents_3.pdf", height=5, width=3)
plot(plot_detailed_CD8)
dev.off()

##################################



### Count analysis of the sig.cor.genes across the celltypes at the detailed level
#####################################
###

library(dplyr)
library(ggplot2)

### Load the data
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
Cor_analysis_DF_rough=Cor_analysis_DF %>% mutate(rough=gsub("\\..*","",celltypes))
Cor_analysis_DF_rough=subset(Cor_analysis_DF_rough) %>% subset(analysis=="all" & rho_pval<0.05)

### Analyze the number of age-correlated genes in the mixed sex in each detailed celltypes
DF_up=Cor_analysis_DF_rough %>% subset(rho>0) %>% group_by(rough, celltypes) %>% count() %>% mutate(`upreg.`=n) %>% select(-n)
DF_down=Cor_analysis_DF_rough %>% subset(rho<0) %>% group_by(rough, celltypes) %>% count() %>% mutate(`downreg.`=n) %>% select(-n)
DF_merged=merge(DF_up, DF_down, by=c("rough","celltypes"))
DF_merged=DF_merged %>% tidyr::pivot_longer(cols=c("upreg.","downreg."), names_to="orientation", values_to="count")
counts_of_detailedCelltype_inEachRough=unlist(lapply(split(DF_merged$celltypes, DF_merged$rough), length))
rough_types=unique(DF_merged$rough)

plots_list=list()
for (i in 1:length(rough_types)) {
  the_rough_type=rough_types[i]
  DF_merged_subset=subset(DF_merged, rough==the_rough_type) %>% arrange(desc(count))
  DF_merged_subset$celltypes=forcats::fct_relevel(DF_merged_subset$celltypes, unique(DF_merged_subset$celltypes))
  
  plot_=
    ggplot(DF_merged_subset, aes(x=celltypes, y=count, alpha=orientation)) +
    geom_linerange(data=subset(DF_merged_subset, orientation=="upreg."),
                   aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=count), size=1, color=paletteer::paletteer_d("ggsci::category20_d3")[i]) +
    geom_linerange(data=subset(DF_merged_subset, orientation=="downreg."),
                   aes(xmin=celltypes, xmax=celltypes, ymin=0, ymax=-count), size=1, color=paletteer::paletteer_d("ggsci::category20_d3")[i]) +
    scale_alpha_manual(values=c(0.4,0.8)) +
    theme_minimal() +
    labs(x=NULL, y="Number of age-correlated genes") +
    theme(axis.text.x=element_text(size=11, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=11),
          axis.title.y=element_text(size=12),
          title=element_text(size=12),
          legend.position="none") +
    scale_y_continuous(labels=abs)
  
  plots_list[[i]]=plot_
}

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Counts_Detailed.pdf", height=10, width=15)
cowplot::plot_grid(cowplot::plot_grid(plotlist=plots_list[1:3], align="hv", nrow=1, ncol=3, rel_widths=counts_of_detailedCelltype_inEachRough[1:3]),
                   cowplot::plot_grid(plotlist=plots_list[4:8], align="hv", nrow=1, ncol=5, rel_widths=counts_of_detailedCelltype_inEachRough[4:8]),
                   align="hv", axis="lr", nrow=2)
dev.off()

#####################################



### Volcano plots to show the sig. age-correlated genes in each celltype at the rough level
# ...* The inter and detailed levels generate plots of too large size and limited information, so take only CD8T cells
#####################################
###

library(dplyr)
library(ggplot2)
library(scRNAtoolVis)

### For Rough level
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# all
df_all=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="all") %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_=
  jjVolcano2(diffData=df_all, 
             # myMarkers=specific_genes$gene,
             topGeneN=3,
             # tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
             aesCol=c("lightgrey","lightgrey"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=4,
             tile.fontsize=4.5,
             tile.linewidth=0.05,
             pSize=0.1,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title=NULL)
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Rough_1.pdf", width=9, height=6)
# plot(plot_)
# dev.off()
# females
df_f=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="females") %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_f=
  jjVolcano2(diffData=df_f, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
             aesCol=c("lightpink","lightpink"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3,
             tile.fontsize=3.5,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# males
df_m=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="males") %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_m=
  jjVolcano2(diffData=df_m, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20],
             aesCol=c("lightblue","lightblue"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3,
             tile.fontsize=3.5,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# save all/females/males
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Rough.pdf", width=length(unique(df_all$cluster)), height=5)
plot(plot_); plot(plot_f); plot(plot_m)
dev.off()

### For Inter level CD8T cells
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
# all
df_all=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="all" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_=
  jjVolcano2(diffData=df_all, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
             aesCol=c("lightgrey","lightgrey"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3,
             tile.fontsize=3.5,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<both>")
# females
df_f=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="females" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_f=
  jjVolcano2(diffData=df_f, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
             aesCol=c("lightpink","lightpink"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3,
             tile.fontsize=3.5,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# males
df_m=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="males" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_m=
  jjVolcano2(diffData=df_m, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[11:20][1:length(unique(df_all$cluster))],
             aesCol=c("lightblue","lightblue"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3,
             tile.fontsize=3.5,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# save all/females/males
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Inter_CD8T.pdf", width=5, height=5)
plot(plot_); plot(plot_f); plot(plot_m)
dev.off()

### For Detailed level CD8T cells
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# all
df_all=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="all" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_=
  jjVolcano2(diffData=df_all, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
             aesCol=c("lightgrey","lightgrey"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3.5,
             tile.fontsize=3,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<both>")
# females
df_f=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="females" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_f=
  jjVolcano2(diffData=df_f, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
             aesCol=c("lightpink","lightpink"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3.5,
             tile.fontsize=3,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<females>")
# males
df_m=Cor_analysis_DF %>% 
  mutate(p_val=rho_pval, avg_log2FC=rho, p_val_adj=rho_pval, cluster=celltypes) %>%
  subset(analysis=="males" & cluster %in% unique(Cor_analysis_DF$celltypes)[grepl("^CD8T\\.",unique(Cor_analysis_DF$celltypes))]) %>%
  mutate(cluster=gsub("^CD8T\\.","",cluster)) %>%
  select(p_val, avg_log2FC, p_val_adj, cluster, gene)
source("~/Rscripts/jjVolcano2.R")
plot_m=
  jjVolcano2(diffData=df_m, 
             # myMarkers=shared_kinases,
             tile.col=paletteer::paletteer_d("ggsci::category20c_d3")[7:20],
             aesCol=c("lightblue","lightblue"),
             back.col="grey97",
             log2FC.cutoff=0.1,
             size=3.5,
             tile.fontsize=3,
             tile.linewidth=0.05,
             pSize=0.5,
             fontface="italic") +
  guides(color=guide_legend(size=3)) +
  theme(title=element_text(size=10), 
        legend.text=element_text(size=12),
        legend.position="none",
        axis.title.y=element_text(size=12)) +
  labs(x=NULL, y=expression(Spearman~rho), title="<males>")
# save all/females/males
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Detailed_CD8T.pdf", width=20, height=5)
plot(plot_); plot(plot_f); plot(plot_m)
dev.off()

##################################



### Analyze the shared age-correlated genes across the celltypes at the rough level
#####################################
###

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

### Load the data
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
Cor_analysis_DF=Cor_analysis_DF %>% subset(analysis=="all")
pos_genes=Cor_analysis_DF %>% subset(rho_pval<0.05 & rho>0)
pos_genes=split(pos_genes$gene, pos_genes$celltypes)
neg_genes=Cor_analysis_DF %>% subset(rho_pval<0.05 & rho<0)
neg_genes=split(neg_genes$gene, neg_genes$celltypes)

### The jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}

### Calculate the jaccard similarity for pos genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(pos_genes)) {
  the_first_celltype=pos_genes[[i]]
  
  for (j in 1:length(pos_genes)) {
    the_second_celltype=pos_genes[[j]]
    
    jaccard_similarity=jaccard(the_first_celltype, the_second_celltype)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(pos_genes)[i],"_vs_",names(pos_genes)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name
Jaccard_SIM_pos=data.frame(pair=names(Jaccard_SIM), pos=Jaccard_SIM)

### Calculate the jaccard similarity for neg genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(neg_genes)) {
  the_first_celltype=neg_genes[[i]]
  
  for (j in 1:length(neg_genes)) {
    the_second_celltype=neg_genes[[j]]
    
    jaccard_similarity=jaccard(the_first_celltype, the_second_celltype)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(neg_genes)[i],"_vs_",names(neg_genes)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name
Jaccard_SIM_neg=data.frame(pair=names(Jaccard_SIM), neg=Jaccard_SIM)

### Merge the pos and neg results
Jaccard_SIM_merged=merge(Jaccard_SIM_pos, Jaccard_SIM_neg, by="pair")
Jaccard_SIM_merged=Jaccard_SIM_merged %>% 
  tidyr::pivot_longer(cols=c("pos","neg"), names_to="orientation", values_to="jaccard") %>%
  mutate(celltype_1=gsub("_vs_.*","",pair), celltype_2=gsub(".*_vs_","",pair))

### Plot the results
# neg
Jaccard_SIM_df=Jaccard_SIM_merged %>% subset(orientation=="pos") %>%
  select(-c(pair, orientation)) %>%
  tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype_1")
Jaccard_SIM_df=-Jaccard_SIM_df
range(Jaccard_SIM_df[Jaccard_SIM_df!=-1]) # check

col_fun=circlize::colorRamp2(c(-0.25, 0), c("dodgerblue3","white"))
ht1=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10),
          column_names_side="top",
          column_names_rot=45,
          column_names_gp=gpar(fontsize=10),
          height=0.1*nrow(Jaccard_SIM_df), 
          width=0.1*ncol(Jaccard_SIM_df))
lgd1=Legend(col_fun=col_fun, 
            title="Jaccard similarity of downreg.",
            legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)

# pos
Jaccard_SIM_df=Jaccard_SIM_merged %>% subset(orientation=="pos") %>%
  select(-c(pair, orientation)) %>%
  tidyr::pivot_wider(names_from="celltype_2", values_from="jaccard") %>%
  tibble::column_to_rownames("celltype_1")
range(Jaccard_SIM_df[Jaccard_SIM_df!=1]) # check

col_fun=circlize::colorRamp2(c(0, 0.25), c("white", "brown3"))
ht2=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10),
          column_names_side="top",
          column_names_rot=45,
          column_names_gp=gpar(fontsize=10),
          height=0.1*nrow(Jaccard_SIM_df), 
          width=0.1*ncol(Jaccard_SIM_df))
lgd2=Legend(col_fun=col_fun, 
           title="Jaccard similarity of upreg.",
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)

ht_final_merged=draw(ht1 + ht2, annotation_legend_list=list(lgd1, lgd2), annotation_legend_side="bottom")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_JaccardSimilarty_amongCelltypes_rough.pdf", height=4, width=6)
draw(ht_final_merged)
dev.off()

#####################################



### Analyze the overlapping between WGCNA modules at the rough and the age-asscoiated genes in each celltype
#####################################
###

library(dplyr)
library(ggplot2)
library(scRNAtoolVis)
library(hdWGCNA)

### Load the data
# age-associated genes
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
# WGCNA modules
WGCNA_obj=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Analysis_Rough.rds")
Module_genes=GetModules(WGCNA_obj); Module_genes$module=as.character(Module_genes$module)

### Map the data
merged_data=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05) %>% 
  right_join(., Module_genes %>% subset(module!="grey") %>% mutate(gene=gene_name) %>% select(gene, module), by="gene")
age_cor_genes=split(merged_data$gene, merged_data$celltypes)
wgcna_genes=split(merged_data$gene, merged_data$module)

### Calculate the jaccard similarity
# the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
# compare each module vs. each celltype of the age-correlated genes
Jaccard_SIM=Jaccard_SIM_Name=c()
for (i in 1:length(age_cor_genes)) {
  the_age_cor=age_cor_genes[[i]]

  for (j in 1:length(wgcna_genes)) {
    the_wgcna=wgcna_genes[[j]]

    jaccard_similarity=jaccard(the_age_cor, the_wgcna)
    Jaccard_SIM=c(Jaccard_SIM, jaccard_similarity)
    Jaccard_SIM_Name=c(Jaccard_SIM_Name, paste0(names(age_cor_genes)[i],"_vs_",names(wgcna_genes)[j]))
  }
}
names(Jaccard_SIM)=Jaccard_SIM_Name

### Plot the results
library(ComplexHeatmap)
Jaccard_SIM_df=data.frame(pair=names(Jaccard_SIM), jaccard=Jaccard_SIM)
Jaccard_SIM_df=Jaccard_SIM_df %>% mutate(celltype=gsub("_vs_.*","",pair), module=gsub(".*_vs_","",pair)) %>%
  select(-pair) %>%
  tidyr::pivot_wider(names_from="module", values_from="jaccard") %>%
  column_to_rownames("celltype")
range(Jaccard_SIM_df) # check

col_fun=circlize::colorRamp2(c(0, 0.1, 0.2), c("dodgerblue3", "white", "brown3"))
ht=
  Heatmap(Jaccard_SIM_df,
          name=" ",
          show_heatmap_legend=FALSE, 
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F,
          row_names_gp=gpar(fontsize=9),
          column_names_side="top",
          column_names_rot=0,
          column_names_centered=T,
          column_names_gp=gpar(fontsize=9),
          height=0.15*nrow(Jaccard_SIM_df), 
          width=0.1*ncol(Jaccard_SIM_df))
lgd=Legend(col_fun=col_fun, 
           title="Jaccard similarity index",
           legend_width=unit(6, "cm"), grid_height=unit(0.2, "cm"),
           direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9),
)
ht_final=draw(ht, annotation_legend_list=lgd, annotation_legend_side="bottom")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_JaccardSimilarty_vsWGCNAmodules_rough.pdf", height=0.5*nrow(Jaccard_SIM_df), width=0.6*ncol(Jaccard_SIM_df))
draw(ht_final)
dev.off()

#####################################



### Plot the correlation between the top30 age-correlated genes and ages in the celltypes at the rough level
#####################################
###

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")
celltypes_all=unique(Cor_analysis_DF$celltypes)
cor_plots_neg=cor_plots_pos=list()
for (i in 1:length(celltypes_all)) {
  ### Plot downregulated genes
  Cor_df_neg=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
    slice_min(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_neg[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
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
  
  ### Plot upregulated genes
  Cor_df_pos=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
    slice_max(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_pos[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
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
  
  ### Save all
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
  plot(cor_plots_neg[[i]])
  dev.off()
  
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
  plot(cor_plots_pos[[i]])
  dev.off()
}
qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_neg.pdf"),
                          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_pos.pdf")),
                  output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough.pdf")
unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_neg.pdf"),
         paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Rough_",c(1:length(celltypes_all)),"_pos.pdf")))

##################################



### Plot the correlation between the top30 age-correlated genes and ages in the CD8T cells at the inter level
#####################################
###

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Inter.csv.gz", sep="\t")
celltypes_all=unique(Cor_analysis_DF$celltypes)
celltypes_all=celltypes_all[grepl("^CD8T\\.",celltypes_all)]
cor_plots_neg=cor_plots_pos=list()
for (i in 1:length(celltypes_all)) {
  ### Plot downregulated genes
  Cor_df_neg=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
    slice_min(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_neg[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
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
  
  ### Plot upregulated genes
  Cor_df_pos=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
    slice_max(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.inter %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_pos[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
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
  
  ### Save all
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
  plot(cor_plots_neg[[i]])
  dev.off()
  
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
  plot(cor_plots_pos[[i]])
  dev.off()
}
qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_neg.pdf"),
                          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_pos.pdf")),
                  output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_CD8T.pdf")
unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_neg.pdf"),
         paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Inter_",c(1:length(celltypes_all)),"_pos.pdf")))

##################################



### Plot the correlation between the top30 age-correlated genes and ages in the CD8T cells at the detailed level
#####################################
###

pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
celltypes_all=unique(Cor_analysis_DF$celltypes)
celltypes_all=celltypes_all[grepl("^CD8T\\.",celltypes_all)]
cor_plots_neg=cor_plots_pos=list()
for (i in 1:length(celltypes_all)) {
  ### Plot downregulated genes
  Cor_df_neg=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho<0) %>%
    slice_min(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_neg=Cor_df_neg[1:30]; Cor_df_neg=Cor_df_neg[!is.na(Cor_df_neg)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_neg), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_neg[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Downregulated genes in ", celltypes_all[i])) +
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
  
  ### Plot upregulated genes
  Cor_df_pos=Cor_analysis_DF %>%
    subset(celltypes==celltypes_all[i] & rho_pval<0.05 & rho>0) %>%
    slice_max(rho, n=100) %>%
    select(gene) %>% subset(!duplicated(.)) %>% unlist() %>% unname()
  Cor_df_pos=Cor_df_pos[1:30]; Cor_df_pos=Cor_df_pos[!is.na(Cor_df_pos)]
  # get the expr data
  pseudo_subset=subset(pseudobulkobj, Annot.detailed %in% celltypes_all[i])
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", Cor_df_pos), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M")) %>%
    tidyr::pivot_longer(cols=colnames(.)[!(colnames(.) %in% c("age","sex","donor_id"))], names_to="gene", values_to="expression")
  # plot
  cor_plots_pos[[i]]=
    ggplot(expr_data, aes(x=age, y=expression, color=sex)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.1) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Upregulated genes in ", celltypes_all[i])) +
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
  
  ### Save all
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",i,"_neg.pdf"), height=ceiling(length(Cor_df_neg)/3)*2, width=9)
  plot(cor_plots_neg[[i]])
  dev.off()
  
  pdf(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",i,"_pos.pdf"), height=ceiling(length(Cor_df_pos)/3)*2, width=9)
  plot(cor_plots_pos[[i]])
  dev.off()
}
qpdf::pdf_combine(input=c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_neg.pdf"),
                          paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_pos.pdf")),
                  output="~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_CD8T.pdf")
unlink(c(paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_neg.pdf"),
         paste0("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Detailed_",c(1:length(celltypes_all)),"_pos.pdf")))

##################################



### Map the reported/known age-associated genes in the Pval-Cor plot
#####################################
###

### Take only the aging-related genes based on DEMAGALHAES_AGING_UP Dataset
AgingDB=read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")
marker_genes=AgingDB$gene

### Load the Correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough.csv.gz", sep="\t")

### Plot the Pval-Cor volcano plots for each celltype
celltypes_all=unique(Cor_analysis_DF$celltypes)
plot_corr=list()
for (i in 1:length(celltypes_all)) {
  Cor_analysis_subset=Cor_analysis_DF %>% 
    subset(celltypes==celltypes_all[i]) %>% 
    mutate(mark=ifelse(gene %in% marker_genes, gene, NA),
           analysis=ifelse(analysis=="all","<both>",ifelse(analysis=="females","<females>","<males>")))
  plot_corr[[i]]=
    ggplot(Cor_analysis_subset, aes(x=rho, y=-log10(rho_pval), color=analysis)) +
    ggrastr::geom_point_rast(size=0.5, shape=20, alpha=0.5) +
    labs(title=paste0("Gene-Age correlation in ",celltypes_all[i]), y=expression(-log[10]~pval), x=expression(Spearman~rho)) +
    facet_wrap(~analysis) +
    theme_classic() +
    theme(title=element_text(size=11),
          plot.subtitle=element_text(size=11, face="italic"),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.position="none",
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.background=element_rect(color="transparent", fill="transparent", linewidth=1.5, linetype=NULL),
          strip.text.x=element_text(size=10, color="black", face="italic")) +
    scale_color_manual(values=c("grey","lightpink","lightblue")) +
    geom_text_repel(aes(label=mark), size=3, show.legend=FALSE, color="black")
}
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_Volcano_Rough_MarkedwKnownGenes.pdf", height=5, width=6)
for(i in 1:length(plot_corr)) plot(plot_corr[[i]])
dev.off()




















RNA.sig.cor_genes_FROM_DMP.neg.cor_genes=RNA.sig.cor_genes_FROM_DMP.pos.cor_genes=list()
neg_plot=pos_plot=list()
for (i in 1:7) {
  pseudo_subset=subset(pseudobulkobj, Annot.rough %in% rough_celltype[[i]])


  ### Analyze the correlation between expr of the hypomethylated genes and ages
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", as.character(Neg.Cor.genes[[i]])), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
  colsum_=colSums(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")])
  filter=c("age", "sex", "donor_id", names(colsum_)[colsum_!=0])
  expr_data=expr_data[, filter]
  genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
  age_correlated_expr_genes=c()
  for (j in 1:length(genes_)) {
    cor_result=cor.test(expr_data$age, expr_data[[genes_[j]]])
    if (cor_result$p.value<0.05) {
      the_cor=cor_result$estimate
      names(the_cor)=genes_[j]
      age_correlated_expr_genes=c(age_correlated_expr_genes, the_cor)
    }
  }
  # save all the (RNA-expr)sig.age-cor.genes from the hypomethylated genes and their correlation with ages
  RNA.sig.cor_genes_FROM_DMP.neg.cor_genes[[i]]=
    data.frame(gene=names(age_correlated_expr_genes),
               cor=age_correlated_expr_genes,
               DMP_type="Significantly hypomethylated",
               celltype=tissues_selected[i])

  # make the df with the top30 abs(cor) genes
  age_correlated_expr_genes=names(sort(abs(age_correlated_expr_genes), decreasing=T)[1:30])
  age_correlated_expr_genes=age_correlated_expr_genes[!is.na(age_correlated_expr_genes)]
  expr_data=expr_data %>%
    select(c("age", "sex", "donor_id", age_correlated_expr_genes))
  expr_data_longDF=tidyr::pivot_longer(expr_data, cols=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex","donor_id")], names_to="gene", values_to="expression")

  # plot
  neg_plot[[i]]=
    ggplot(expr_data_longDF, aes(x=age, y=expression, color=sex)) +
    geom_point(size=0.1, shape=20, alpha=0.05) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Expression of the hypomethylated genes in ", tissues_selected[i])) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          # legend.position="none",
          # legend.title=element_blank(),
          title=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
    ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5)

  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",i,"_neg.pdf"), height=ceiling(length(age_correlated_expr_genes)/3)*2, width=9)
  plot(neg_plot[[i]])
  dev.off()


  ### Analyze the correlation between expr of the hypermethylated genes and ages
  expr_data=FetchData(pseudo_subset, vars=c("age", "sex", "donor_id", as.character(Pos.Cor.genes[[i]])), layer="data")
  expr_data=expr_data %>%
    mutate(sex=ifelse(sex=="female", "F", "M"))
  colsum_=colSums(expr_data[! colnames(expr_data) %in% c("age","sex","donor_id")])
  filter=c("age", "sex", "donor_id", names(colsum_)[colsum_!=0])
  expr_data=expr_data[, filter]
  genes_=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex", "donor_id")]
  age_correlated_expr_genes=c()
  for (j in 1:length(genes_)) {
    cor_result=cor.test(expr_data$age, expr_data[[genes_[j]]])
    if (cor_result$p.value<0.05) {
      the_cor=cor_result$estimate
      names(the_cor)=genes_[j]
      age_correlated_expr_genes=c(age_correlated_expr_genes, the_cor)
    }
  }
  # save all the (RNA-expr)sig.age-cor.genes from the hypermethylated genes and their correlation with ages
  RNA.sig.cor_genes_FROM_DMP.pos.cor_genes[[i]]=
    data.frame(gene=names(age_correlated_expr_genes),
               cor=age_correlated_expr_genes,
               DMP_type="Significantly hypermethylated",
               celltype=tissues_selected[i])

  # make the df with the top30 abs(cor) genes
  age_correlated_expr_genes=names(sort(abs(age_correlated_expr_genes), decreasing=T)[1:30])
  age_correlated_expr_genes=age_correlated_expr_genes[!is.na(age_correlated_expr_genes)]
  expr_data=expr_data %>%
    select(c("age", "sex", "donor_id", age_correlated_expr_genes))
  expr_data_longDF=tidyr::pivot_longer(expr_data, cols=colnames(expr_data)[! colnames(expr_data) %in% c("age","sex","donor_id")], names_to="gene", values_to="expression")

  # plot
  pos_plot[[i]]=
    ggplot(expr_data_longDF, aes(x=age, y=expression, color=sex)) +
    geom_point(size=0.1, shape=20, alpha=0.05) +
    # ggsci::scale_color_d3() +
    facet_wrap(~gene, ncol=3, scales="free") +
    theme_classic() +
    labs(x="age", y="Expression", title=paste0("Expression of the hypermethylated genes in ", tissues_selected[i])) +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          # legend.position="none",
          # legend.title=element_blank(),
          title=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL)) +
    geom_smooth(aes(group=sex), method="loess", show.legend=FALSE, linewidth=0.5, fill="gray93") +
    ggpubr::stat_cor(aes(group=sex), show.legend=FALSE, size=3.5)

  pdf(paste0("~/Project_PBMCage/Other_omics/Methy_plots/DMP_TSS_SigAssocWithAges_RNAexpr_CorPlots_",i,"_pos.pdf"), height=ceiling(length(age_correlated_expr_genes)/3)*2, width=9)
  plot(pos_plot[[i]])
  dev.off()

  message(paste0("Plotting of the expr of the RNA.sig.age-cor.genes in ",tissues_selected[i], " is done!"))
}

### Save all the (RNA-expr)sig.age-cor.genes from the (DMP-beta)sig.age-cor.genes and their correlation with ages
RNA.sig.cor_genes_WITH_cor=data.table::rbindlist(c(RNA.sig.cor_genes_FROM_DMP.neg.cor_genes, RNA.sig.cor_genes_FROM_DMP.pos.cor_genes))
saveRDS(RNA.sig.cor_genes_WITH_cor, "~/Project_PBMCage/Other_omics/MethyTSS_subsets/RelatedGenes_RNAexpression_ageCorrelationResult.rds")




