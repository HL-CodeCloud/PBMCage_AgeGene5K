

### Compare the prediction ability between Xu_5000 geneset and the ScRNA_ageGeneset on the 17PBMC datasets
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="rough" # CHANGE MANUALLY!!!

pseudobulkobj=readRDS(paste0("~/Project_PBMCage/GSE213516_17PBMC/17PBMC_annotated_AssignedOnly_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj

### Set the sliding window
steps=seq(50, 4000, 50)

### Run with PLSR model
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
RSQs=list()
RSQs=
  foreach(i=1:length(steps)) %dopar% {
    # have to reload libraries in each core
    library(Seurat)
    library(dplyr)
    library(pls)
    
    geneN_chosen=steps[i]
    
    ### Take the genes from the genesets according to the geneN
    Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
    Gene_weights_5000=Gene_weights_5000 %>% select(symbol, mean) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% slice_max(mean, n=geneN_chosen)
    colnames(Gene_weights_5000)=c("gene","weight")
    scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
    scRNA_ageGenes=scRNA_ageGenes %>% 
      mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
      group_by(Gene.name) %>% summarize_at("Vip", mean) %>% slice_max(Vip, n=geneN_chosen) %>% as.data.frame()
    colnames(scRNA_ageGenes)=c("gene","weight")
    
    ### Take the genes from each list
    genelist=list(Xu_list=Gene_weights_5000$gene,
                  scRNA_ageGenes=scRNA_ageGenes$gene)
    
    ### Run on Xu_gene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[1]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    Xu_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    ### Run on scRNA_ageGene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[2]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    scRNA.AgeGenes_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    # return both
    RSQs[[i]]=list(Xu_5000=Xu_results,
                   scRNA.AgeGenes=scRNA.AgeGenes_results)
  }
stopCluster(cl)

### Rearrange the results
RSQs_df_Xu=data.frame(topN=steps,
                      `R-squared`=unlist(lapply(RSQs, function(x) x[[1]][1])),
                      MAE=unlist(lapply(RSQs, function(x) x[[1]][2])),
                      MAD=unlist(lapply(RSQs, function(x) x[[1]][3])),
                      RMSE=unlist(lapply(RSQs, function(x) x[[1]][4])),
                      pls_selectedN=unlist(lapply(RSQs, function(x) x[[1]][5])))
RSQs_df_Xu$geneset="Xu_5000.ranked"
RSQs_df_Xu=RSQs_df_Xu %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_df_scRNA_AgeGenes=data.frame(topN=steps,
                                  `R-squared`=unlist(lapply(RSQs, function(x) x[[2]][1])),
                                  MAE=unlist(lapply(RSQs, function(x) x[[2]][2])),
                                  MAD=unlist(lapply(RSQs, function(x) x[[2]][3])),
                                  RMSE=unlist(lapply(RSQs, function(x) x[[2]][4])),
                                  pls_selectedN=unlist(lapply(RSQs, function(x) x[[2]][5])))
RSQs_df_scRNA_AgeGenes$geneset="scRNA_ageGenes.ranked"
RSQs_df_scRNA_AgeGenes=RSQs_df_scRNA_AgeGenes %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_both_DF=data.table::rbindlist(list(RSQs_df_Xu, RSQs_df_scRNA_AgeGenes))
write.table(RSQs_both_DF, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters_in17PBMC_",celltype.level,".txt"), sep="\t")

#####################################



### Compare the prediction ability between Xu_5000 geneset and all other published genesets on the 17PBMC datasets
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="rough" # CHANGE MANUALLY!!!
pseudobulkobj=readRDS(paste0("~/Project_PBMCage/GSE213516_17PBMC/17PBMC_annotated_AssignedOnly_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj

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
# cut on 1000 -> we do not need so many steps
gene_list_steps=lapply(gene_list_steps, function(x) x[x<=1000])

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
write.table(RSQs_DF_long, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_in17PBMC_",celltype.level,".txt"), sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the results
RSQs_DF_long_extend_all=list()
# at the rough level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_in17PBMC_rough.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level1"
RSQs_DF_long_extend$celltype.level="coarse"
RSQs_DF_long_extend_all[["rough"]]=RSQs_DF_long_extend
# at the inter level
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_in17PBMC_inter.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level2"
RSQs_DF_long_extend$celltype.level="inter."
RSQs_DF_long_extend_all[["inter"]]=RSQs_DF_long_extend
# at the detailed level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_in17PBMC_detailed.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level3"
RSQs_DF_long_extend$celltype.level="fine"
RSQs_DF_long_extend_all[["detailed"]]=RSQs_DF_long_extend
# merged
RSQs_DF_long_extend=data.table::rbindlist(RSQs_DF_long_extend_all)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)
RSQs_DF_long_extend$celltype.level=forcats::fct_relevel(RSQs_DF_long_extend$celltype.level, 
                                                        # c("celltype.level1","celltype.level2","celltype.level3"),
                                                        c("coarse","inter.","fine"))

# plot_legend=
#   ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
#   geom_line(alpha=0) +
#   ggsci::scale_color_d3("category20") +
#   theme_void() +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         legend.position=c(0.5,0.5)) +
#   guides(color=guide_legend(override.aes=list(linetype=c("solid",rep("dashed", 12)), alpha=1), 
#                             ncol=5, 
#                             title=NULL))
# pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_donor_legend.pdf", width=10, height=1)
# plot(plot_legend)
# dev.off()

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_grid(celltype.level~parameter, scales="free", switch="y") +
  ggsci::scale_color_d3("category20") +
  theme_minimal() +
  labs(x="geneN (x100)", y="parameter", title="GSE213516") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_blank(),
        strip.text.y=element_text(size=10),
        panel.spacing.y=unit(1, "lines"),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_y_continuous(position="right", n.breaks=4) +
  scale_x_continuous(n.breaks=4, labels=function(x) x/100) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_in17PBMC_celltypeLevel.pdf", width=1.75, height=3)
plot(plot_)
dev.off()

#####################################



### Compare the prediction ability between Xu_5000 geneset and the ScRNA_ageGeneset on the CovidCtrl dataset
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="detailed" # CHANGE THIS MANUALLY!!!

pseudobulkobj=readRDS(paste0("~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj

### Set the sliding window
steps=seq(50, 4000, 50)

### Run with PLSR model
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
RSQs=list()
RSQs=
  foreach(i=1:length(steps)) %dopar% {
    # have to reload libraries in each core
    library(Seurat)
    library(dplyr)
    library(pls)
    
    geneN_chosen=steps[i]
    
    ### Take the genes from the genesets according to the geneN
    Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
    Gene_weights_5000=Gene_weights_5000 %>% select(symbol, mean) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% slice_max(mean, n=geneN_chosen)
    colnames(Gene_weights_5000)=c("gene","weight")
    scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
    scRNA_ageGenes=scRNA_ageGenes %>% 
      mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
      group_by(Gene.name) %>% summarize_at("Vip", mean) %>% slice_max(Vip, n=geneN_chosen) %>% as.data.frame()
    colnames(scRNA_ageGenes)=c("gene","weight")
    
    ### Take the genes from each list
    genelist=list(Xu_list=Gene_weights_5000$gene,
                  scRNA_ageGenes=scRNA_ageGenes$gene)
    
    ### Run on Xu_gene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[1]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    Xu_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    ### Run on scRNA_ageGene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[2]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    scRNA.AgeGenes_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    # return both
    RSQs[[i]]=list(Xu_5000=Xu_results,
                   scRNA.AgeGenes=scRNA.AgeGenes_results)
  }
stopCluster(cl)

### Rearrange the results
RSQs_df_Xu=data.frame(topN=steps,
                      `R-squared`=unlist(lapply(RSQs, function(x) x[[1]][1])),
                      MAE=unlist(lapply(RSQs, function(x) x[[1]][2])),
                      MAD=unlist(lapply(RSQs, function(x) x[[1]][3])),
                      RMSE=unlist(lapply(RSQs, function(x) x[[1]][4])),
                      pls_selectedN=unlist(lapply(RSQs, function(x) x[[1]][5])))
RSQs_df_Xu$geneset="Xu_5000.ranked"
RSQs_df_Xu=RSQs_df_Xu %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_df_scRNA_AgeGenes=data.frame(topN=steps,
                                  `R-squared`=unlist(lapply(RSQs, function(x) x[[2]][1])),
                                  MAE=unlist(lapply(RSQs, function(x) x[[2]][2])),
                                  MAD=unlist(lapply(RSQs, function(x) x[[2]][3])),
                                  RMSE=unlist(lapply(RSQs, function(x) x[[2]][4])),
                                  pls_selectedN=unlist(lapply(RSQs, function(x) x[[2]][5])))
RSQs_df_scRNA_AgeGenes$geneset="scRNA_ageGenes.ranked"
RSQs_df_scRNA_AgeGenes=RSQs_df_scRNA_AgeGenes %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_both_DF=data.table::rbindlist(list(RSQs_df_Xu, RSQs_df_scRNA_AgeGenes))
write.table(RSQs_both_DF, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters_inCovidCtrl_",celltype.level,".txt"), sep="\t")

#####################################



### Compare the prediction ability between Xu_5000 geneset and all other published genesets on the CovidCtrl datasets
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="rough" # CHANGE MANUALLY!!!
pseudobulkobj=readRDS(paste0("~/Project_PBMCage/GSE158055/CovidCtrl_annotated_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj

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
# cut on 1000 -> we do not need so many steps
gene_list_steps=lapply(gene_list_steps, function(x) x[x<=1000])

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
write.table(RSQs_DF_long, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCovidCtrl_",celltype.level,".txt"), sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the results
RSQs_DF_long_extend_all=list()
# at the rough level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCovidCtrl_rough.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level1"
RSQs_DF_long_extend$celltype.level="coarse"
RSQs_DF_long_extend_all[["rough"]]=RSQs_DF_long_extend
# at the inter level
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCovidCtrl_inter.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level2"
RSQs_DF_long_extend$celltype.level="inter."
RSQs_DF_long_extend_all[["inter"]]=RSQs_DF_long_extend
# at the detailed level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCovidCtrl_detailed.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level3"
RSQs_DF_long_extend$celltype.level="fine"
RSQs_DF_long_extend_all[["detailed"]]=RSQs_DF_long_extend
# merged
RSQs_DF_long_extend=data.table::rbindlist(RSQs_DF_long_extend_all)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)
RSQs_DF_long_extend$celltype.level=forcats::fct_relevel(RSQs_DF_long_extend$celltype.level, 
                                                        # c("celltype.level1","celltype.level2","celltype.level3"),
                                                        c("coarse","inter.","fine"))

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_grid(celltype.level~parameter, scales="free", switch="y") +
  ggsci::scale_color_d3("category20") +
  theme_minimal() +
  labs(x="geneN (x100)", y="parameter", title="GSE158055") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_blank(),
        strip.text.y=element_text(size=10),
        panel.spacing.y=unit(1, "lines"),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_y_continuous(position="right", n.breaks=4) +
  scale_x_continuous(n.breaks=4, labels=function(x) x/100) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inCovidCtrl_celltypeLevel.pdf", width=1.75, height=3)
plot(plot_)
dev.off()

#####################################



### Compare the prediction ability between Xu_5000 geneset and the ScRNA_ageGeneset on the SC dataset
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="detailed" # CHANGE THIS MANUALLY!!!

pseudobulkobj=readRDS(paste0("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj
ncol(pseudobulkobj_donor)

### Set the sliding window
steps=seq(50, 4000, 50)

### Run with PLSR model
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
RSQs=list()
RSQs=
  foreach(i=1:length(steps)) %dopar% {
    # have to reload libraries in each core
    library(Seurat)
    library(dplyr)
    library(pls)
    
    geneN_chosen=steps[i]
    
    ### Take the genes from the genesets according to the geneN
    Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
    Gene_weights_5000=Gene_weights_5000 %>% select(symbol, mean) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% slice_max(mean, n=geneN_chosen)
    colnames(Gene_weights_5000)=c("gene","weight")
    scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
    scRNA_ageGenes=scRNA_ageGenes %>% 
      mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
      group_by(Gene.name) %>% summarize_at("Vip", mean) %>% slice_max(Vip, n=geneN_chosen) %>% as.data.frame()
    colnames(scRNA_ageGenes)=c("gene","weight")
    
    ### Take the genes from each list
    genelist=list(Xu_list=Gene_weights_5000$gene,
                  scRNA_ageGenes=scRNA_ageGenes$gene)
    
    ### Run on Xu_gene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[1]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    Xu_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    ### Run on scRNA_ageGene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[2]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    scRNA.AgeGenes_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    # return both
    RSQs[[i]]=list(Xu_5000=Xu_results,
                   scRNA.AgeGenes=scRNA.AgeGenes_results)
  }
stopCluster(cl)

### Rearrange the results
RSQs_df_Xu=data.frame(topN=steps,
                      `R-squared`=unlist(lapply(RSQs, function(x) x[[1]][1])),
                      MAE=unlist(lapply(RSQs, function(x) x[[1]][2])),
                      MAD=unlist(lapply(RSQs, function(x) x[[1]][3])),
                      RMSE=unlist(lapply(RSQs, function(x) x[[1]][4])),
                      pls_selectedN=unlist(lapply(RSQs, function(x) x[[1]][5])))
RSQs_df_Xu$geneset="Xu_5000.ranked"
RSQs_df_Xu=RSQs_df_Xu %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_df_scRNA_AgeGenes=data.frame(topN=steps,
                                  `R-squared`=unlist(lapply(RSQs, function(x) x[[2]][1])),
                                  MAE=unlist(lapply(RSQs, function(x) x[[2]][2])),
                                  MAD=unlist(lapply(RSQs, function(x) x[[2]][3])),
                                  RMSE=unlist(lapply(RSQs, function(x) x[[2]][4])),
                                  pls_selectedN=unlist(lapply(RSQs, function(x) x[[2]][5])))
RSQs_df_scRNA_AgeGenes$geneset="scRNA_ageGenes.ranked"
RSQs_df_scRNA_AgeGenes=RSQs_df_scRNA_AgeGenes %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_both_DF=data.table::rbindlist(list(RSQs_df_Xu, RSQs_df_scRNA_AgeGenes))
write.table(RSQs_both_DF, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters_inJapanSC_",celltype.level,".txt"), sep="\t")

#####################################



### Compare the prediction ability between Xu_5000 geneset and all other published genesets on the SC datasets
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="rough" # CHANGE MANUALLY!!!
pseudobulkobj=readRDS(paste0("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj

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
# cut on 1000 -> we do not need so many steps
gene_list_steps=lapply(gene_list_steps, function(x) x[x<=1000])

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
write.table(RSQs_DF_long, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inJapanSC_",celltype.level,".txt"), sep="\t")
unlink("~/Project_PBMCage/Tempt_RDS/Xu_5000Genes_prediction.MAE.analysis.rds")

### Plot the results
RSQs_DF_long_extend_all=list()
# at the rough level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inJapanSC_rough.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level1"
RSQs_DF_long_extend$celltype.level="coarse"
RSQs_DF_long_extend_all[["rough"]]=RSQs_DF_long_extend
# at the inter level
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inJapanSC_inter.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level2"
RSQs_DF_long_extend$celltype.level="inter."
RSQs_DF_long_extend_all[["inter"]]=RSQs_DF_long_extend
# at the detailed level
unique(RSQs_DF_long$geneset) # check
RSQs_DF_long=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inJapanSC_detailed.txt", sep="\t")
RSQs_DF_long=RSQs_DF_long %>% select(-lasso_selectedN) %>% subset(!is.na(MAE)) %>% subset(MAE!=0)
RSQs_DF_long_extend=RSQs_DF_long %>% tidyr::pivot_longer(cols=3:6, names_to="parameter", values_to="value")
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(!is.na(value) & value!=0)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
RSQs_DF_long_extend$parameter=forcats::fct_relevel(RSQs_DF_long_extend$parameter, "MAD","MAE","RMSE","R-squared")
# RSQs_DF_long_extend$celltype.level="celltype.level3"
RSQs_DF_long_extend$celltype.level="fine"
RSQs_DF_long_extend_all[["detailed"]]=RSQs_DF_long_extend
# merged
RSQs_DF_long_extend=data.table::rbindlist(RSQs_DF_long_extend_all)
RSQs_DF_long_extend=RSQs_DF_long_extend %>% subset(parameter!="R-squared" & parameter!="MAD" & topN<1000)
RSQs_DF_long_extend$celltype.level=forcats::fct_relevel(RSQs_DF_long_extend$celltype.level, 
                                                        # c("celltype.level1","celltype.level2","celltype.level3"),
                                                        c("coarse","inter.","fine"))

plot_=
  ggplot(RSQs_DF_long_extend, aes(x=topN, y=value, color=geneset)) +
  facet_grid(celltype.level~parameter, scales="free", switch="y") +
  ggsci::scale_color_d3("category20") +
  theme_minimal() +
  labs(x="geneN (x100)", y="parameter", title="SC cohort") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_blank(),
        strip.text.y=element_text(size=10),
        panel.spacing.y=unit(1, "lines"),
        legend.position="top",
        # legend.title=element_blank(),
        title=element_text(size=10)) +
  scale_y_continuous(position="right", n.breaks=4) +
  scale_x_continuous(n.breaks=4, labels=function(x) x/100) +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset!="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="dashed") +
  geom_smooth(data=subset(RSQs_DF_long_extend, geneset=="AgeGene5K"), aes(group=geneset), 
              method="loess", show.legend=FALSE, linewidth=0.3, se=FALSE, alpha=1, linetype="solid")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.PublishedGenesets_PredictParameters_inJapanSC_celltypeLevel.pdf", width=1.75, height=3)
plot(plot_)
dev.off()

#####################################



### Compare the prediction ability between Xu_5000 geneset and the ScRNA_ageGeneset on all 3 verification datasets
#####################################
###

library(Seurat)
library(dplyr)

### Load the bulk data
celltype.level="inter" # CHANGE THIS MANUALLY!!!

pseudobulkobj=readRDS(paste0("~/Project_PBMCage/For_or_From_Xu/Verification_datasets_Pseudobulk_",celltype.level,".rds"))
pseudobulkobj_donor=pseudobulkobj
ncol(pseudobulkobj_donor)

### Set the sliding window
steps=seq(50, 4000, 50)

### Run with PLSR model
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl)
RSQs=list()
RSQs=
  foreach(i=1:length(steps)) %dopar% {
    # have to reload libraries in each core
    library(Seurat)
    library(dplyr)
    library(pls)
    
    geneN_chosen=steps[i]
    
    ### Take the genes from the genesets according to the geneN
    Gene_weights_5000=readRDS("~/Project_PBMCage/For_or_From_Xu/Gene_weights_5000.rds")
    Gene_weights_5000=Gene_weights_5000 %>% select(symbol, mean) %>% mutate(symbol=gsub("_ENSG","-ENSG",symbol)) %>% slice_max(mean, n=geneN_chosen)
    colnames(Gene_weights_5000)=c("gene","weight")
    scRNA_ageGenes=read.csv("~/Project_PBMCage/Zhu_2023_scRNA_AgeGenes.csv", sep=",")
    scRNA_ageGenes=scRNA_ageGenes %>% 
      mutate(Gene.name=gsub("HLA\\.","HLA-",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.AS1","-AS1",Gene.name)) %>% 
      mutate(Gene.name=gsub("\\.2HG","-2HG",Gene.name)) %>% 
      group_by(Gene.name) %>% summarize_at("Vip", mean) %>% slice_max(Vip, n=geneN_chosen) %>% as.data.frame()
    colnames(scRNA_ageGenes)=c("gene","weight")
    
    ### Take the genes from each list
    genelist=list(Xu_list=Gene_weights_5000$gene,
                  scRNA_ageGenes=scRNA_ageGenes$gene)
    
    ### Run on Xu_gene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[1]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    Xu_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    ### Run on scRNA_ageGene list
    pseudo_data=FetchData(pseudobulkobj_donor, vars=c("age","sex",genelist[[2]]), layer="data")
    colnames(pseudo_data)=gsub("-","_",colnames(pseudo_data))
    # make the train and test sets
    x=pseudo_data %>% select(colnames(.)[!colnames(.) %in% c("age")])
    y=pseudo_data %>% select("age") %>% unlist() %>% unname()
    set_trainset_prob=0.9
    ind=sample(c(TRUE, FALSE), nrow(x)*set_trainset_prob, replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
    xy_train=pseudo_data[ind, ]
    xy_test=pseudo_data[!ind, ]
    x_train=xy_train[,c(2:ncol(xy_train))]; y_train=xy_train[[1]]
    x_test=xy_test[,c(2:ncol(xy_test))]; y_test=xy_test[[1]]
    if (!all(x_train==0)) { # glmnet is only doable when x_train is not all zero
      # fit
      all_genes_codes=paste0(colnames(x_train), collapse="+")
      eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=xy_train, scale=TRUE, validation='CV')")))
      # determine the best number of comp
      # summary(Model)
      cverr=RMSEP(Model)$val[1,,]
      imin_cv=which.min(cverr)-1
      r2=R2(Model)$val[1,,]
      imax_r2=which.max(r2)-1
      best.ncomp=as.integer(mean(imin_cv, imax_r2))
      if (best.ncomp!=0) {
        # predict
        y_pred=predict(Model, x_test, ncomp=best.ncomp)[,,1]
        # calculate parameters
        SST=sum((y_train-mean(y_train))^2)
        SSE=sum((y_pred-y_test)^2)
        rsq=1-SSE/SST
        MAE=Metrics::mae(y_pred, y_test)
        MAD=mad(y_pred, y_test)
        RMSE=Metrics::rmse(y_pred, y_test)
        gene_coeff=unlist(coefficients(Model))[2:nrow(coefficients(Model))]
        gene_n=length(gene_coeff[gene_coeff!=0])
      } else {rsq=MAE=MAD=RMSE=gene_n=0}
    } else {
      rsq=MAE=MAD=RMSE=gene_n=0
    }
    scRNA.AgeGenes_results=c(rsq, MAE, MAD, RMSE, gene_n)
    
    # return both
    RSQs[[i]]=list(Xu_5000=Xu_results,
                   scRNA.AgeGenes=scRNA.AgeGenes_results)
  }
stopCluster(cl)

### Rearrange the results
RSQs_df_Xu=data.frame(topN=steps,
                      `R-squared`=unlist(lapply(RSQs, function(x) x[[1]][1])),
                      MAE=unlist(lapply(RSQs, function(x) x[[1]][2])),
                      MAD=unlist(lapply(RSQs, function(x) x[[1]][3])),
                      RMSE=unlist(lapply(RSQs, function(x) x[[1]][4])),
                      pls_selectedN=unlist(lapply(RSQs, function(x) x[[1]][5])))
RSQs_df_Xu$geneset="Xu_5000.ranked"
RSQs_df_Xu=RSQs_df_Xu %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_df_scRNA_AgeGenes=data.frame(topN=steps,
                                  `R-squared`=unlist(lapply(RSQs, function(x) x[[2]][1])),
                                  MAE=unlist(lapply(RSQs, function(x) x[[2]][2])),
                                  MAD=unlist(lapply(RSQs, function(x) x[[2]][3])),
                                  RMSE=unlist(lapply(RSQs, function(x) x[[2]][4])),
                                  pls_selectedN=unlist(lapply(RSQs, function(x) x[[2]][5])))
RSQs_df_scRNA_AgeGenes$geneset="scRNA_ageGenes.ranked"
RSQs_df_scRNA_AgeGenes=RSQs_df_scRNA_AgeGenes %>% tidyr::pivot_longer(cols=colnames(.)[!colnames(.) %in% c("topN","geneset")], names_to="parameter", values_to="value")

RSQs_both_DF=data.table::rbindlist(list(RSQs_df_Xu, RSQs_df_scRNA_AgeGenes))
write.table(RSQs_both_DF, paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters_inALL3_",celltype.level,".txt"), sep="\t")

#####################################



### Plot the prediction ability comparisons between the genesets in all datasets
#####################################
###

library(ggplot2)
library(dplyr)

### Combine all the results
file=c("17PBMC_","CovidCtrl_","JapanSC_","ALL3_")
dataset_name=c("GSE213516_","GSE158055_","SC cohort_","Combined_")
celltype.level=c("rough","inter","detailed")
file.celltypelevel=unlist(lapply(file, function(x) paste0(x, celltype.level)))
dataset.celltypelevel=unlist(lapply(dataset_name, function(x) paste0(x, celltype.level)))

RSQs_both_DF_list=list()
for (i in 1:length(file.celltypelevel)) {
  RSQs_both_DF=read.delim(paste0("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters_in",file.celltypelevel[i],".txt"), sep="\t")
  RSQs_both_DF=RSQs_both_DF %>%
    mutate(geneset=ifelse(geneset=="Xu_5000.ranked","AgeGene5K","Zhu et al. (2023)")) %>%
    mutate(parameter=ifelse(parameter=="R.squared","R-squared",parameter))
  RSQs_both_DF$verification_file=file.celltypelevel[i]
  RSQs_both_DF$verification_dataset=gsub("_.*","",dataset.celltypelevel[i])
  RSQs_both_DF$celltype.level=gsub(".*_","",dataset.celltypelevel[i])
  
  RSQs_both_DF_list[[i]]=RSQs_both_DF
}
RSQs_both_DF_all=data.table::rbindlist(RSQs_both_DF_list)
data.table::fwrite(RSQs_both_DF_all, "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters.csv.gz", sep="\t")

### Plot
RSQs_both_DF=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters.csv.gz", sep="\t")
RSQs_both_DF$celltype.level=as.factor(RSQs_both_DF$celltype.level)
levels(RSQs_both_DF$celltype.level)=c("celltype.l3","celltype.l2","celltype.l1")
RSQs_both_DF$celltype.level=forcats::fct_relevel(RSQs_both_DF$celltype.level, "celltype.l1","celltype.l2","celltype.l3")
RSQs_both_DF$geneset=forcats::fct_relevel(RSQs_both_DF$geneset, "AgeGene5K", "Zhu et al. (2023)")
RSQs_both_DF$parameter=forcats::fct_relevel(RSQs_both_DF$parameter, "MAD","MAE","RMSE","R-squared","pls_selectedN")

plot_all.datasets=list()
for (i in 1:length(unique(RSQs_both_DF$verification_dataset))) {
  RSQs_both_DF_dataset=subset(RSQs_both_DF, verification_dataset==unique(RSQs_both_DF$verification_dataset)[i])
  plot_all.celltype.levels=list()
  for (j in 1:length(unique(RSQs_both_DF_dataset$celltype.level))) {
    RSQs_both_DF_dataset.celllevel=subset(RSQs_both_DF_dataset, celltype.level==unique(RSQs_both_DF_dataset$celltype.level)[j])
    topN_to_remove=RSQs_both_DF_dataset.celllevel %>% subset(value==0) %>% select(topN) %>% unlist() %>% unname()
    RSQs_both_DF_dataset.celllevel=RSQs_both_DF_dataset.celllevel %>% subset(!(topN %in% topN_to_remove) & parameter!="pls_selectedN")
    plot_=
      ggplot(RSQs_both_DF_dataset.celllevel, aes(x=topN, y=value, color=geneset)) +
      ggrastr::geom_point_rast(size=0.1, shape=20, alpha=0.1) +
      facet_wrap(~parameter, ncol=4, scales="free") +
      # ggsci::scale_color_d3("category20c") +
      scale_color_manual(values=c("brown3","grey50")) +
      theme_classic() +
      labs(x="geneN", y="parameter", 
           title=unique(RSQs_both_DF_dataset.celllevel$verification_dataset), 
           subtitle=paste0("<",unique(RSQs_both_DF_dataset.celllevel$celltype.level),">")) +
      theme(axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=8),
            axis.title.x=element_text(size=9),
            axis.title.y=element_text(size=9),
            legend.position="top",
            # legend.title=element_blank(),
            title=element_text(size=10), 
            plot.subtitle=element_text(size=9, face="italic")) +
      guides(color=guide_legend(override.aes=list(size=2, alpha=1))) +
      geom_smooth(aes(group=geneset), method="loess", show.legend=FALSE, linewidth=0.5, se=FALSE, alpha=0.5)
    plot_all.celltype.levels[[j]]=plot_
  }
  plot_all.datasets[[i]]=cowplot::plot_grid(plotlist=plot_all.celltype.levels, align="hv", ncol=1)
}

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Compare.to.scRNAAgeGenes_PredictParameters.pdf", height=9, width=8)
for (i in 1:length(plot_all.datasets)) {plot(plot_all.datasets[[i]])}
dev.off()

#####################################


