
### Prepare the published genesets for comparison
#####################################
###
### !: to compare, here scRNA_AgeGene geneset is ranking regardless the celltypes

### Load Xu's list
Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
Gene_weights_5000=Gene_weights_5000 %>% select(symbol, mean) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% arrange(desc(mean))
colnames(Gene_weights_5000)=c("gene","weight")
Xu_5000.ranked=Gene_weights_5000$gene

### To compare with the published results, prepare their genesets
bulkRNA_ageGenes=read.csv("~/Project_PBMCage/Publicated_bulkRNA_AgeGenes.csv", sep=",")
scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
scRNA_ageGenes=scRNA_ageGenes %>% group_by(Gene.name) %>% summarize_at("Vip", mean) %>% arrange(desc(Vip)) %>% as.data.frame()
colnames(scRNA_ageGenes)=c("gene","weight")
Demagalhaes_AgingUp_2023=clusterProfiler::read.gmt("~/Project_PBMCage/DEMAGALHAES_AGING_UP.v2023.2.Hs.gmt")

### Make the list to go through
gene_list=as.list(bulkRNA_ageGenes)
gene_list=lapply(gene_list, function(x) x[grepl("[A-Z]|[a-z]|[0-9]", x)])
gene_list=c(gene_list, list(Demagalhaes_2023.unranked=Demagalhaes_AgingUp_2023$gene,
                            scRNA_ageGenes.ranked=scRNA_ageGenes$gene))
# rearrange the unranked genesets according to the Xu_5000 set and cut to the same length as the Xu_5000
for (i in 1:length(gene_list)) {
  if (grepl("\\.unranked", names(gene_list)[i])) {
    gene_list[[i]]=gene_list[[i]][match(Xu_5000.ranked, gene_list[[i]])]
    gene_list[[i]]=gene_list[[i]][!is.na(gene_list[[i]])]
  }
}

### Save the genesets
gene_list=c(gene_list, list(Xu_5000.ranked=Xu_5000.ranked))
saveRDS(gene_list, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")

#####################################



### Compare the performance among all the published age-associated genesets on the PBMCage dataset
#####################################
###

library(Seurat)
library(glmnet)
library(dplyr)

### Make the pseudobulk on donor_id
pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj_donor=AggregateExpression(pseudobulkobj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
pseudobulkobj_donor$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[1])))
pseudobulkobj_donor$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[2]))
pseudobulkobj_donor$donor_id=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[3]))

### Load genesets for comparison
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
# cut the genesets all to the same length as Xu_5000 to compare
for (i in 1:length(gene_list)) {
  gene_list[[i]]=gene_list[[i]][1:min(length(gene_list[[i]]), length(Xu_5000.ranked))]
}
length(gene_list) # to check

### Analyze the prediction ability of 5000 genes with sliding windows
# set the steps according to the Xu_5000 geneset to compare across all genesets
steps_=seq(2, length(Xu_5000.ranked), 10)
gene_list_steps=lapply(gene_list, function(x) {
  ifelse(length(x)>max(steps_), 
         list(steps_), 
         ifelse(length(x)<100, list(seq(2, length(x), 1)), list(seq(2, length(x), 10)))
  )
})
gene_list_steps=lapply(gene_list_steps, function(x) x[[1]])

# manually run for each genesets (... to avoid dealing with the complicated env. viable problem)
RSQs_all=list()

for (geneset_idx in 1:length(gene_list)) {
  library(doParallel)
  # run on parallel
  cl=makeCluster(10)
  registerDoParallel(cl)
  steps_=gene_list_steps[[geneset_idx]]
  determine_if_shuffled=grepl("\\.unranked",names(gene_list_steps)[geneset_idx])
  
  Parame=list()
  Parame=
    foreach(i=1:length(steps_)) %dopar% {
      # have to reload libraries in each core
      library(Seurat)
      library(dplyr)
      library(glmnet)
      
      # take the genes
      topn=steps_[i]
      if (determine_if_shuffled) {
        genes_selected=gene_list[[geneset_idx]] %>% sample(topn)
      } else {
        genes_selected=gene_list[[geneset_idx]][1:topn]
      }
      
      pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex","donor_id",genes_selected), layer="data")
      
      # make the train and test sets
      x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age","donor_id")]) %>% as.matrix()
      y=pseudo_data %>% select("age") %>% unlist() %>% unname()
      set_trainset_prob=0.9
      ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
      x_train=x[ind, ]; y_train=y[ind]
      x_test=x[!ind, ]; y_test=y[!ind]
      
      rsq=MAE=MAD=RMSE=gene_n=0
      
      tryCatch({
        # fit
        cv_model=cv.glmnet(x_train, y_train, alpha=1)
        best_lambda=cv_model$lambda.min
        # plot(cv_model)
        Model=glmnet(x_train, y_train, alpha=1, lambda=best_lambda)
        
        # predict the test sets
        y_pred=predict(Model, s=best_lambda, newx=x_test)
        
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred[,1], y_test)
        MAD=mad(y_pred[,1], y_test)
        RMSE=Metrics::rmse(y_pred[,1], y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      }, error=function(msg) print("Skip."))
      
      Parame[[i]]=c(rsq, MAE, MAD, RMSE, gene_n)
    }
  stopCluster(cl)
  
  RSQs_df=data.frame(topN=steps_,
                     `R-squared`=unlist(lapply(Parame, function(x) x[1])),
                     MAE=unlist(lapply(Parame, function(x) x[2])),
                     MAD=unlist(lapply(Parame, function(x) x[3])),
                     RMSE=unlist(lapply(Parame, function(x) x[4])),
                     lasso_selectedN=unlist(lapply(Parame, function(x) x[5])))
  
  RSQs_all[[geneset_idx]]=RSQs_df
}

### Save all the results
names(RSQs_all)=names(gene_list)
saveRDS(RSQs_all, "~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Rearrange the results
RSQs_all=readRDS("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")
geneset_result=data.table::rbindlist(lapply(1:length(RSQs_all), function(i) as.data.frame(rep(names(RSQs_all)[i], nrow(RSQs_all[[i]])))))
topN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[1]])))
Rsquared_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[2]])))
MAE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[3]])))
MAD_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[4]])))
RMSE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[5]])))
lasso_selectedN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[6]])))
RSQs_DF_long=do.call(cbind,
                     list(geneset_result, topN_result, Rsquared_result, MAE_result, MAD_result, RMSE_result, lasso_selectedN_result))
colnames(RSQs_DF_long)=c("info", "topN", "R-squared", "MAE", "MAD", "RMSE", "lasso_selectedN")
# change the names of the genesets
RSQs_DF_long$geneset=as.factor(RSQs_DF_long$info)
levels(RSQs_DF_long$geneset) # check
levels(RSQs_DF_long$geneset)=c("de Magalhães et al. (2009)","GenAge Build 21","Dørum et al. (2024)","Yang et al. (2015)","Chuffa et al. (2022)",
                               "Peters et al. (2015)","Shokhirev et al. (2021); cubist-DG","Shokhirev et al. (2021); cubist-VG",
                               "Shokhirev et al. (2021); Rf-DG","Shokhirev et al. (2021); Rf-VG","Ren et al. (2020)","Zhu et al. (2023)",
                               "AgeGene5K")
write.table(RSQs_DF_long, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_perDonorID.txt", sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the R-squared
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_perDonorID.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")

RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)

plot_legend=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  geom_line(alpha=0) +
  ggsci::scale_color_d3("category20") +
  theme_void() +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position=c(0.5,0.5)) +
  guides(color=guide_legend(override.aes=list(linetype=c("solid",rep("dashed", 12)), alpha=1), 
                            nrow=5, 
                            title=NULL))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_donor_legend.pdf", width=6.5, height=1.2)
plot(plot_legend)
dev.off()

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_wrap(~parameter, ncol=1, scales="free") +
  ggsci::scale_color_d3("category20") +
  theme_classic() +
  labs(x="geneN", y="parameter", title="Yazar et al. (2022) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_donor.pdf", width=4, height=5)
plot(plot_)
dev.off()

#####################################



### Compare the performance among all the published age-associated genesets on the Combined_verification dataset
#####################################
###

library(Seurat)
library(glmnet)
library(dplyr)

### Make the pseudobulk
seurat_obj=readRDS("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Integrated.rds")
sex_df=data.frame(donor_id=seurat_obj$donor_id,
                  sex=seurat_obj$sex)
sex_df=sex_df[!duplicated(sex_df$donor_id),]
pseudobulk=AggregateExpression(seurat_obj, assays="RNA", return.seurat=T, group.by=c("donor_id","age"))
pseudobulk$donor_id=gsub("_.*","",colnames(pseudobulk))
pseudobulk$age=as.numeric(gsub(".*_","",colnames(pseudobulk)))
pseudobulk=AddMetaData(pseudobulk, metadata=sex_df$sex[match(pseudobulk$donor_id, as.character(sex_df$donor_id))], col.name="sex")
pseudobulk$sex=as.factor(pseudobulk$sex)
saveRDS(pseudobulk, "~/Project_PBMCage/For_or_From_Xu/Results/Verification_datasets_Pseudobulk_perDonorID.rds")

### Load the pseudobulk on donor_id
pseudobulkobj_donor=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Verification_datasets_Pseudobulk_perDonorID.rds")
unique(pseudobulkobj_donor$donor_id) # remove the CTs since they do not have accurate ages but only periods
donors_remained=unique(pseudobulkobj_donor$donor_id)[!grepl("^CT[0-9]", unique(pseudobulkobj_donor$donor_id))]
pseudobulkobj_donor=subset(pseudobulkobj_donor, donor_id %in% donors_remained)

### Load genesets for comparison
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
# cut the genesets all to the same length as Xu_5000 to compare
for (i in 1:length(gene_list)) {
  gene_list[[i]]=gene_list[[i]][1:min(length(gene_list[[i]]), length(Xu_5000.ranked))]
}
length(gene_list) # to check
lapply(gene_list, length)

### Analyze the prediction ability of 5000 genes with sliding windows
# set the steps according to the Xu_5000 geneset to compare across all genesets
steps_=seq(2, length(Xu_5000.ranked), 10)
gene_list_steps=lapply(gene_list, function(x) {
  ifelse(length(x)>max(steps_), 
         list(steps_), 
         ifelse(length(x)<100, list(seq(2, length(x), 1)), list(seq(2, length(x), 10)))
         )
})
gene_list_steps=lapply(gene_list_steps, function(x) x[[1]])

# manually run for each genesets (... to avoid dealing with the complicated env. viable problem)
RSQs_all=list()

for (geneset_idx in 1:length(gene_list)) {
  library(doParallel)
  # run on parallel
  cl=makeCluster(10)
  registerDoParallel(cl)
  steps_=gene_list_steps[[geneset_idx]]
  determine_if_shuffled=grepl("\\.unranked",names(gene_list_steps)[geneset_idx])
  
  Parame=list()
  Parame=
    foreach(i=1:length(steps_)) %dopar% {
      # have to reload libraries in each core
      library(Seurat)
      library(dplyr)
      library(glmnet)
      
      # take the genes
      topn=steps_[i]
      if (determine_if_shuffled) {
        genes_selected=gene_list[[geneset_idx]] %>% sample(topn)
      } else {
        genes_selected=gene_list[[geneset_idx]][1:topn]
      }
      
      pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex","donor_id",genes_selected), layer="data")
      
      # make the train and test sets
      x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age","donor_id")]) %>% as.matrix()
      y=pseudo_data %>% select("age") %>% unlist() %>% unname()
      set_trainset_prob=0.9
      ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
      x_train=x[ind, ]; y_train=y[ind]
      x_test=x[!ind, ]; y_test=y[!ind]
      
      rsq=MAE=MAD=RMSE=gene_n=0
      
      tryCatch({
        # fit
        cv_model=cv.glmnet(x_train, y_train, alpha=1)
        best_lambda=cv_model$lambda.min
        # plot(cv_model)
        Model=glmnet(x_train, y_train, alpha=1, lambda=best_lambda)
        
        # predict the test sets
        y_pred=predict(Model, s=best_lambda, newx=x_test)
        
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred[,1], y_test)
        MAD=mad(y_pred[,1], y_test)
        RMSE=Metrics::rmse(y_pred[,1], y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      }, error=function(msg) print("Skip."))
      
      Parame[[i]]=c(rsq, MAE, MAD, RMSE, gene_n)
    }
  stopCluster(cl)
  
  RSQs_df=data.frame(topN=steps_,
                     `R-squared`=unlist(lapply(Parame, function(x) x[1])),
                     MAE=unlist(lapply(Parame, function(x) x[2])),
                     MAD=unlist(lapply(Parame, function(x) x[3])),
                     RMSE=unlist(lapply(Parame, function(x) x[4])),
                     lasso_selectedN=unlist(lapply(Parame, function(x) x[5])))
  
  RSQs_all[[geneset_idx]]=RSQs_df
}

### Save all the results
names(RSQs_all)=names(gene_list)
saveRDS(RSQs_all, "~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Rearrange the results
RSQs_all=readRDS("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")
geneset_result=data.table::rbindlist(lapply(1:length(RSQs_all), function(i) as.data.frame(rep(names(RSQs_all)[i], nrow(RSQs_all[[i]])))))
topN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[1]])))
Rsquared_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[2]])))
MAE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[3]])))
MAD_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[4]])))
RMSE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[5]])))
lasso_selectedN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[6]])))
RSQs_DF_long=do.call(cbind,
                     list(geneset_result, topN_result, Rsquared_result, MAE_result, MAD_result, RMSE_result, lasso_selectedN_result))
colnames(RSQs_DF_long)=c("info", "topN", "R-squared", "MAE", "MAD", "RMSE", "lasso_selectedN")
# change the names of the genesets
RSQs_DF_long$geneset=as.factor(RSQs_DF_long$info)
levels(RSQs_DF_long$geneset) # check
levels(RSQs_DF_long$geneset)=c("de Magalhães et al. (2009)","GenAge Build 21","Dørum et al. (2024)","Yang et al. (2015)","Chuffa et al. (2022)",
                               "Peters et al. (2015)","Shokhirev et al. (2021); cubist-DG","Shokhirev et al. (2021); cubist-VG",
                               "Shokhirev et al. (2021); Rf-DG","Shokhirev et al. (2021); Rf-VG","Ren et al. (2020)","Zhu et al. (2023)",
                               "AgeGene5K")
write.table(RSQs_DF_long, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_perDonorID.txt", sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the R-squared
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_perDonorID.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")

RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_wrap(~parameter, ncol=1, scales="free") +
  ggsci::scale_color_d3("category20") +
  theme_classic() +
  labs(x="geneN", y="parameter", title="GSE213516, GSE158055, SC cohort") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_donor.pdf", width=4, height=5)
plot(plot_)
dev.off()

#####################################



### Compare the performance among all the published age-associated genesets on a bulkseq dataset
#####################################
###

library(Seurat)
library(glmnet)
library(dplyr)

### Clean the bulkseq data
counts.data=read.delim("~/Project_PBMCage/For_or_From_Xu/VerificationDatasets/Shokhirev et al. (2021)/raw_filtered.txt", sep="\t")
rownames(counts.data)=counts.data$X
counts.data=counts.data[,2:ncol(counts.data)]
meta.data=read.delim("~/Project_PBMCage/For_or_From_Xu/VerificationDatasets/Shokhirev et al. (2021)/meta_filtered.txt", sep="\t")
rownames(meta.data)=meta.data$SRR.ID
meta.data=meta.data[,3:ncol(meta.data)]
colnames(meta.data)=c("GEO","donor_id","age","sex", colnames(meta.data)[5:9])
meta.data$sex=as.factor(meta.data$sex)
levels(meta.data$sex)=c("F","M","Unknown")
bulkseq.seu=CreateSeuratObject(counts=counts.data, min.cells=0, min.features=0, meta.data=meta.data)
bulkseq.seu=subset(bulkseq.seu, Condition=="Healthy" & Tissue=="Blood;PBMC")
bulkseq.seu=NormalizeData(bulkseq.seu)
saveRDS(bulkseq.seu, "~/Project_PBMCage/For_or_From_Xu/VerificationDatasets/Shokhirev et al. (2021)/Bulkseq_Shokhirev.rds")


pseudobulkobj_donor=readRDS("~/Project_PBMCage/For_or_From_Xu/VerificationDatasets/Shokhirev et al. (2021)/Bulkseq_Shokhirev.rds")
### Load genesets for comparison
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
# cut the genesets all to the same length as Xu_5000 to compare
for (i in 1:length(gene_list)) {
  gene_list[[i]]=gene_list[[i]][1:min(length(gene_list[[i]]), length(Xu_5000.ranked))]
}
length(gene_list) # to check

### Analyze the prediction ability of 5000 genes with sliding windows
# set the steps according to the Xu_5000 geneset to compare across all genesets
steps_=seq(2, length(Xu_5000.ranked), 10)
gene_list_steps=lapply(gene_list, function(x) {
  ifelse(length(x)>max(steps_), 
         list(steps_), 
         ifelse(length(x)<100, list(seq(2, length(x), 1)), list(seq(2, length(x), 10)))
  )
})
gene_list_steps=lapply(gene_list_steps, function(x) x[[1]])

# manually run for each genesets (... to avoid dealing with the complicated env. viable problem)
RSQs_all=list()

for (geneset_idx in 1:length(gene_list)) {
  library(doParallel)
  # run on parallel
  cl=makeCluster(10)
  registerDoParallel(cl)
  steps_=gene_list_steps[[geneset_idx]]
  determine_if_shuffled=grepl("\\.unranked",names(gene_list_steps)[geneset_idx])
  
  Parame=list()
  Parame=
    foreach(i=1:length(steps_)) %dopar% {
      # have to reload libraries in each core
      library(Seurat)
      library(dplyr)
      library(glmnet)
      
      # take the genes
      topn=steps_[i]
      if (determine_if_shuffled) {
        genes_selected=gene_list[[geneset_idx]] %>% sample(topn)
      } else {
        genes_selected=gene_list[[geneset_idx]][1:topn]
      }
      
      pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex","donor_id",genes_selected), layer="data")
      
      # make the train and test sets
      x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age","sex","donor_id")]) %>% as.matrix()
      y=pseudo_data %>% select("age") %>% unlist() %>% unname()
      set_trainset_prob=0.9
      ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
      x_train=x[ind, ]; y_train=y[ind]
      x_test=x[!ind, ]; y_test=y[!ind]
      
      rsq=MAE=MAD=RMSE=gene_n=0
      
      tryCatch({
        # fit
        cv_model=cv.glmnet(x_train, y_train, alpha=1)
        best_lambda=cv_model$lambda.min
        # plot(cv_model)
        Model=glmnet(x_train, y_train, alpha=1, lambda=best_lambda)
        
        # predict the test sets
        y_pred=predict(Model, s=best_lambda, newx=x_test)
        
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred[,1], y_test)
        MAD=mad(y_pred[,1], y_test)
        RMSE=Metrics::rmse(y_pred[,1], y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      }, error=function(msg) print("Skip."))
      
      Parame[[i]]=c(rsq, MAE, MAD, RMSE, gene_n)
    }
  stopCluster(cl)
  
  RSQs_df=data.frame(topN=steps_,
                     `R-squared`=unlist(lapply(Parame, function(x) x[1])),
                     MAE=unlist(lapply(Parame, function(x) x[2])),
                     MAD=unlist(lapply(Parame, function(x) x[3])),
                     RMSE=unlist(lapply(Parame, function(x) x[4])),
                     lasso_selectedN=unlist(lapply(Parame, function(x) x[5])))
  
  RSQs_all[[geneset_idx]]=RSQs_df  # CHANGE THIS MANUALLY !!!
}

### Save all the results
names(RSQs_all)=names(gene_list)
saveRDS(RSQs_all, "~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Rearrange the results
RSQs_all=readRDS("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")
geneset_result=data.table::rbindlist(lapply(1:length(RSQs_all), function(i) as.data.frame(rep(names(RSQs_all)[i], nrow(RSQs_all[[i]])))))
topN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[1]])))
Rsquared_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[2]])))
MAE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[3]])))
MAD_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[4]])))
RMSE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[5]])))
lasso_selectedN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[6]])))
RSQs_DF_long=do.call(cbind,
                     list(geneset_result, topN_result, Rsquared_result, MAE_result, MAD_result, RMSE_result, lasso_selectedN_result))
colnames(RSQs_DF_long)=c("info", "topN", "R-squared", "MAE", "MAD", "RMSE", "lasso_selectedN")
# change the names of the genesets
RSQs_DF_long$geneset=as.factor(RSQs_DF_long$info)
levels(RSQs_DF_long$geneset) # check
levels(RSQs_DF_long$geneset)=c("de Magalhães et al. (2009)","GenAge Build 21","Dørum et al. (2024)","Yang et al. (2015)","Chuffa et al. (2022)",
                               "Peters et al. (2015)","Shokhirev et al. (2021); cubist-DG","Shokhirev et al. (2021); cubist-VG",
                               "Shokhirev et al. (2021); Rf-DG","Shokhirev et al. (2021); Rf-VG","Ren et al. (2020)","Zhu et al. (2023)",
                               "AgeGene5K")
write.table(RSQs_DF_long, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.txt", sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the R-squared
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")

RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_wrap(~parameter, ncol=1, scales="free") +
  ggsci::scale_color_d3("category20") +
  theme_classic() +
  labs(x="geneN", y="parameter", title="Shokhirev et al. (2021) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.pdf", width=4, height=5)
plot(plot_)
dev.off()

#####################################



### Compare the performance among all the published age-associated genesets on the immunity dataset
#####################################
###

library(Seurat)
library(glmnet)
library(dplyr)

### Make the pseudobulk on donor_id
pseudobulkobj_donor=readRDS("~/Project_PBMCage/Immunity/Immunity_all_Pseudobulk_onTube.rds")

### Load genesets for comparison
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
# cut the genesets all to the same length as Xu_5000 to compare
for (i in 1:length(gene_list)) {
  gene_list[[i]]=gene_list[[i]][1:min(length(gene_list[[i]]), length(Xu_5000.ranked))]
}
length(gene_list) # to check

### Analyze the prediction ability of 5000 genes with sliding windows
# set the steps according to the Xu_5000 geneset to compare across all genesets
steps_=seq(2, length(Xu_5000.ranked), 10)
gene_list_steps=lapply(gene_list, function(x) {
  ifelse(length(x)>max(steps_), 
         list(steps_), 
         ifelse(length(x)<100, list(seq(2, length(x), 1)), list(seq(2, length(x), 10)))
  )
})
gene_list_steps=lapply(gene_list_steps, function(x) x[[1]])

# manually run for each genesets (... to avoid dealing with the complicated env. viable problem)
RSQs_all=list()

for (geneset_idx in 1:length(gene_list)) {
  library(doParallel)
  # run on parallel
  cl=makeCluster(10)
  registerDoParallel(cl)
  steps_=gene_list_steps[[geneset_idx]]
  determine_if_shuffled=grepl("\\.unranked",names(gene_list_steps)[geneset_idx])
  
  Parame=list()
  Parame=
    foreach(i=1:length(steps_)) %dopar% {
      # have to reload libraries in each core
      library(Seurat)
      library(dplyr)
      library(glmnet)
      
      # take the genes
      topn=steps_[i]
      if (determine_if_shuffled) {
        genes_selected=gene_list[[geneset_idx]] %>% sample(topn)
      } else {
        genes_selected=gene_list[[geneset_idx]][1:topn]
      }
      
      pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex","donor_id",genes_selected), layer="data")
      
      # make the train and test sets
      x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age","donor_id")]) %>% as.matrix()
      y=pseudo_data %>% select("age") %>% unlist() %>% unname()
      set_trainset_prob=0.9
      ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
      x_train=x[ind, ]; y_train=y[ind]
      x_test=x[!ind, ]; y_test=y[!ind]
      
      rsq=MAE=MAD=RMSE=gene_n=0
      
      tryCatch({
        # fit
        cv_model=cv.glmnet(x_train, y_train, alpha=1)
        best_lambda=cv_model$lambda.min
        # plot(cv_model)
        Model=glmnet(x_train, y_train, alpha=1, lambda=best_lambda)
        
        # predict the test sets
        y_pred=predict(Model, s=best_lambda, newx=x_test)
        
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred[,1], y_test)
        MAD=mad(y_pred[,1], y_test)
        RMSE=Metrics::rmse(y_pred[,1], y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      }, error=function(msg) print("Skip."))
      
      Parame[[i]]=c(rsq, MAE, MAD, RMSE, gene_n)
    }
  stopCluster(cl)
  
  RSQs_df=data.frame(topN=steps_,
                     `R-squared`=unlist(lapply(Parame, function(x) x[1])),
                     MAE=unlist(lapply(Parame, function(x) x[2])),
                     MAD=unlist(lapply(Parame, function(x) x[3])),
                     RMSE=unlist(lapply(Parame, function(x) x[4])),
                     lasso_selectedN=unlist(lapply(Parame, function(x) x[5])))
  
  RSQs_all[[geneset_idx]]=RSQs_df
}

### Save all the results
names(RSQs_all)=names(gene_list)
saveRDS(RSQs_all, "~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Rearrange the results
RSQs_all=readRDS("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")
geneset_result=data.table::rbindlist(lapply(1:length(RSQs_all), function(i) as.data.frame(rep(names(RSQs_all)[i], nrow(RSQs_all[[i]])))))
topN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[1]])))
Rsquared_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[2]])))
MAE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[3]])))
MAD_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[4]])))
RMSE_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[5]])))
lasso_selectedN_result=data.table::rbindlist(lapply(RSQs_all, function(x) as.data.frame(x[[6]])))
RSQs_DF_long=do.call(cbind,
                     list(geneset_result, topN_result, Rsquared_result, MAE_result, MAD_result, RMSE_result, lasso_selectedN_result))
colnames(RSQs_DF_long)=c("info", "topN", "R-squared", "MAE", "MAD", "RMSE", "lasso_selectedN")
# change the names of the genesets
RSQs_DF_long$geneset=as.factor(RSQs_DF_long$info)
levels(RSQs_DF_long$geneset) # check
levels(RSQs_DF_long$geneset)=c("de Magalhães et al. (2009)","GenAge Build 21","Dørum et al. (2024)","Yang et al. (2015)","Chuffa et al. (2022)",
                               "Peters et al. (2015)","Shokhirev et al. (2021); cubist-DG","Shokhirev et al. (2021); cubist-VG",
                               "Shokhirev et al. (2021); Rf-DG","Shokhirev et al. (2021); Rf-VG","Ren et al. (2020)","Zhu et al. (2023)",
                               "AgeGene5K")
write.table(RSQs_DF_long, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inImmunity_perDonorID.txt", sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the R-squared
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inImmunity_perDonorID.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")

RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_wrap(~parameter, ncol=1, scales="free") +
  ggsci::scale_color_d3("category20") +
  theme_classic() +
  labs(x="geneN", y="parameter", title="Terekhova et al. (2023) dataset") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inImmunity_donor.pdf", width=4, height=5)
plot(plot_)
dev.off()

#####################################



### Combine all the results of performance comparison between the two genesets on bulk and pseudobulk-donor_id datasets
#####################################
###

qpdf::pdf_combine(input=c("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_donor.pdf",
                          "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_donor.pdf",
                          "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.pdf",
                          "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inImmunity_donor.pdf"),
                  output="~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_donor.pdf")
unlink(c("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_donor.pdf",
         "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_donor.pdf",
         "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.pdf",
         "~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inImmunity_donor.pdf"))

#####################################



### Plot the Combination of all the results of performance comparison on all the test datasets
#####################################
###

### Load the results
res_Yazar=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_perDonorID.txt", sep="\t") %>% 
  select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0) %>% 
  tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value") %>% 
  subset(!is.na(value) & value!=0) %>% 
  mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
res_Yazar$parameter=forcats::fct_relevel(res_Yazar$parameter, "MAD","MAE","RMSE","R-squared")
res_Yazar=res_Yazar %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000) %>%
  mutate(test.dataset="Yazar dataset")

res_Terekhova=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inImmunity_perDonorID.txt", sep="\t") %>% 
  select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0) %>% 
  tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value") %>% 
  subset(!is.na(value) & value!=0) %>% 
  mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
res_Terekhova$parameter=forcats::fct_relevel(res_Terekhova$parameter, "MAD","MAE","RMSE","R-squared")
res_Terekhova=res_Terekhova %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000) %>%
  mutate(test.dataset="Terekhova dataset")

res_combination=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_perDonorID.txt", sep="\t") %>% 
  select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0) %>% 
  tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value") %>% 
  subset(!is.na(value) & value!=0) %>% 
  mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
res_combination$parameter=forcats::fct_relevel(res_combination$parameter, "MAD","MAE","RMSE","R-squared")
res_combination=res_combination %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000) %>%
  mutate(
    # test.dataset="GSE213516, GSE158055, SC cohort",
    test.dataset="Combined dataset"
         )

res_bulk=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.txt", sep="\t") %>% 
  select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0) %>% 
  tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value") %>% 
  subset(!is.na(value) & value!=0) %>% 
  mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
res_bulk$parameter=forcats::fct_relevel(res_bulk$parameter, "MAD","MAE","RMSE","R-squared")
res_bulk=res_bulk %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000) %>%
  mutate(test.dataset="Shokhirev dataset")

### Combine the results
res_all=data.table::rbindlist(list(res_Yazar, res_Terekhova, res_combination, res_bulk))
res_all$test.dataset=forcats::fct_relevel(res_all$test.dataset,
                                          c("Yazar dataset", "Terekhova dataset", 
                                            # "GSE213516, GSE158055, SC cohort", 
                                            "Combined dataset",
                                            "Shokhirev dataset"))


plot_all=
  ggplot(res_all, aes(x=topN, y=value, color=geneset)) +
  facet_grid(parameter~test.dataset, scales="free") +
  ggsci::scale_color_d3("category20") +
  theme_classic() +
  labs(x="geneN (x100)", y="parameter", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_x_continuous(labels=function(x) x/100) +
  geom_smooth(data=subset(res_all, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(res_all, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inAllTestDatasets.pdf", width=6, height=2.5)
plot(plot_all)
dev.off()