# 
# ### Find the best geneN
# #####################################
# ###
# 
# ### Load the performance comparison results
# # on PBMCage_donor, Combo_donor, bulkseq
# Result_DF_path=c("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inPBMCage_perDonorID.txt",
#                  "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inCombinedVeri_perDonorID.txt",
#                  "~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.PublishedGenesets_PredictParameters_inBulkRNA.txt")
# 
# Selected_TopNs=list()
# for (i in 1:length(Result_DF_path)) {
#   RSQs=read.delim(Result_DF_path[i], sep="\t")
#   RSQs_Rsquared=RSQs %>% select(topN, R.squared, geneset) %>% group_by(topN) %>% filter(R.squared==max(R.squared, na.rm=T)) %>% arrange(topN) %>% 
#     subset(geneset=="AgeGene5K") %>%
#     select(topN) %>% unlist() %>% unname()
#   RSQs_MAE=RSQs %>% select(topN, MAE, geneset) %>% group_by(topN) %>% filter(MAE==min(MAE, na.rm=T)) %>% arrange(topN) %>% 
#     subset(geneset=="AgeGene5K") %>%
#     select(topN) %>% unlist() %>% unname()
#   RSQs_MAD=RSQs %>% select(topN, MAD, geneset) %>% group_by(topN) %>% filter(MAD==min(MAD, na.rm=T)) %>% arrange(topN) %>% 
#     subset(geneset=="AgeGene5K") %>%
#     select(topN) %>% unlist() %>% unname()
#   RSQs_RMSE=RSQs %>% select(topN, RMSE, geneset) %>% group_by(topN) %>% filter(RMSE==min(RMSE, na.rm=T)) %>% arrange(topN) %>% 
#     subset(geneset=="AgeGene5K") %>%
#     select(topN) %>% unlist() %>% unname()
#   RSQs_topN=Reduce(intersect, list(RSQs_Rsquared, RSQs_MAE, RSQs_MAD, RSQs_RMSE))
#   
#   Selected_TopNs[[i]]=RSQs_topN
# }
# # on 3 different scRNA datasets, per celltype.level
# RSQs=data.table::fread("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.scRNAAgeGenes_PredictParameters.csv.gz", sep="\t")
# parames=unique(RSQs$parameter)[!grepl("pls_selectedN", unique(RSQs$parameter))]
# 
# Selected_TopNs_para=list()
# for (i in 1:length(parames)) {
#   if (parames[i]=="R-squared") {
#     RSQs_topN=RSQs %>% subset(parameter==parames[i]) %>% group_by(topN, parameter, verification_file) %>% filter(value==max(value, na.rm=T)) %>%
#       subset(geneset=="AgeGene5K") %>% ungroup() %>% count(topN) %>%
#       subset(n>=length(unique(RSQs$verification_file))/3) %>%
#       arrange(topN) %>% select(topN) %>% unlist() %>% unname()
#   } else {
#     RSQs_topN=RSQs %>% subset(parameter==parames[i]) %>% group_by(topN, parameter, verification_file) %>% filter(value==min(value, na.rm=T)) %>%
#       subset(geneset=="AgeGene5K") %>% ungroup() %>% count(topN) %>%
#       subset(n>=length(unique(RSQs$verification_file))/3) %>%
#       arrange(topN) %>% select(topN) %>% unlist() %>% unname()
#   }
#   
#   Selected_TopNs_para[[i]]=RSQs_topN
# }
# # combine both
# Selected_TopNs=c(Selected_TopNs, list(Reduce(intersect, Selected_TopNs_para)))
# 
# ### Select the TopN that looks better:
# # TopN=50, 100, 200
# # !: I chose 300 at last.
# xxx=lapply(Selected_TopNs, function(x) round(x/100)*100)
# Reduce(intersect, xxx)
# 
# #####################################






###############################################################################################################
##################################### Train model based on the PBMCage RNA #####################################
###############################################################################################################

RESULT_DFs=list()
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")

#####################################



### Check the performance with VerifiCombo dataset
#####################################
###

library(dplyr)
library(UCell)
library(Seurat)
library(pls)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

### Load the PBMCage dataset
# pseudobulk on donor_id
training_PBMCage=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
training_PBMCage=AggregateExpression(training_PBMCage, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
training_PBMCage$age=as.numeric(unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[1])))
training_PBMCage$sex=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[2]))
training_PBMCage$donor_id=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[3]))
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(training_PBMCage, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
training_PBMCage$CD8T.naive=score_df
# extract the topN gene expression
training_PBMCage_expr=FetchData(training_PBMCage, vars=c("age","sex","CD8T.naive",topN_genes), layer="counts")

## The testing data
test_verificombo=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Verification_datasets_Pseudobulk_perDonorID.rds")
unique(test_verificombo$donor_id) # remove the CTs since they do not have accurate ages but only periods
donors_remained=unique(test_verificombo$donor_id)[!grepl("^CT[0-9]|^SC[0-9]", unique(test_verificombo$donor_id))]
test_verificombo=subset(test_verificombo, donor_id %in% donors_remained)
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_verificombo, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_verificombo$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_verificombo_expr=FetchData(test_verificombo, vars=c("age","sex","CD8T.naive", topN_genes), layer="counts")

### Synchronize the genes occurred in both the training and the testing dataset
# in case gene symbols here do not match in the training and the testing dataset
feature.exist.In_train=c(colnames(training_PBMCage_expr), gsub("-","_",colnames(training_PBMCage_expr)), gsub("_","-",colnames(training_PBMCage_expr)))
feature.exist.In_test=c(colnames(test_verificombo_expr), gsub("-","_",colnames(test_verificombo_expr)), gsub("_","-",colnames(test_verificombo_expr)))
coexist.genes=intersect(feature.exist.In_train, feature.exist.In_test)
# update the training and the testing dataset
training_PBMCage_expr=training_PBMCage_expr %>% select(any_of(coexist.genes))
test_verificombo_expr=test_verificombo_expr %>% select(any_of(coexist.genes))
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr.t=data.table::transpose(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
colnames(training_PBMCage_expr.t)=rownames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
rownames(training_PBMCage_expr.t)=colnames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
test_verificombo_expr.t=data.table::transpose(test_verificombo_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
colnames(test_verificombo_expr.t)=rownames(test_verificombo_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
rownames(test_verificombo_expr.t)=colnames(test_verificombo_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))

### Integrate the training and the testing datasets
training_obj=CreateSeuratObject(counts=training_PBMCage_expr.t, min.cells=0, min.features=0, meta.data=training_PBMCage_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive"))))
training_obj$dataset="train"
testing_obj=CreateSeuratObject(counts=test_verificombo_expr.t, min.cells=0, min.features=0, meta.data=test_verificombo_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive"))))
testing_obj$dataset="test"
train_test.mergedobj=merge(training_obj, testing_obj, add.cell.ids=c("train","test"))
# take the highest acceptable npcs_dims and kweight
npcs_dims=ifelse(ncol(training_obj)>1000, 100, min(ncol(training_obj),nrow(training_obj)))
kweight=min(ncol(training_obj), ncol(testing_obj))
# integrate
train_test.mergedobj_integrated=train_test.mergedobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs=npcs_dims) %>%
  FindNeighbors(dims=1:(npcs_dims-1)) %>%
  FindClusters(resolution=0.5)
train_test.mergedobj_integrated=train_test.mergedobj_integrated %>% 
  IntegrateLayers(method=RPCAIntegration, k.weight=30, orig.reduction="pca", new.reduction="integrated.rpca")

### Extract the integrated.rpca for model training
# training data
training_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="train")@reductions$integrated.rpca@cell.embeddings)
training_rpca.df=t(scale(t(as.matrix(training_rpca.df)))) %>% as.data.frame()
training_rpca.df$age=training_obj$age
training_rpca.df$CD8T.naive=training_obj$CD8T.naive
training_rpca.df$sex=training_obj$sex
training_rpca.df$sex=as.factor(training_rpca.df$sex)
levels(training_rpca.df$sex)
levels(training_rpca.df$sex)=c("female","male")
# testing data
testing_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="test")@reductions$integrated.rpca@cell.embeddings)
testing_rpca.df=t(scale(t(as.matrix(testing_rpca.df)))) %>% as.data.frame()
testing_rpca.df$age=testing_obj$age
testing_rpca.df$CD8T.naive=testing_obj$CD8T.naive
testing_rpca.df$sex=testing_obj$sex
testing_rpca.df$sex=as.factor(testing_rpca.df$sex)
levels(testing_rpca.df$sex)
levels(testing_rpca.df$sex)=c("female","male")

### Test only on the ages that have occurred more than 10 times in the training data
ggplot2::ggplot(data.frame(no=1:length(training_rpca.df$age), age=training_rpca.df$age), ggplot2::aes(x=age)) +
  ggplot2::geom_density()
# testing_rpca.df=testing_rpca.df %>% subset(age %in% c(60:80))
# possible_prediction=as.numeric(names(table(training_rpca.df$age)[table(training_rpca.df$age)>15]))

### Create the final train and test sets
# train set
train.x=training_rpca.df %>% select(-age)
train.y=training_rpca.df %>% select(age) %>% unlist() %>% unname()
train.xANDy=cbind(train.x, age=train.y)
# test set
test.x=testing_rpca.df %>% select(-age)
test.y=testing_rpca.df %>% select(age) %>% unlist() %>% unname()

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]
# calculate parameters
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE) # check

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
RESULT_DFs[[1]]=result_df
names(RESULT_DFs)[1]="VerifiCombo"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")

#####################################



### Check the performance with SLE dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

### Load the PBMCage dataset
# pseudobulk on donor_id
training_PBMCage=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
training_PBMCage=AggregateExpression(training_PBMCage, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
training_PBMCage$age=as.numeric(unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[1])))
training_PBMCage$sex=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[2]))
training_PBMCage$donor_id=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[3]))
training_PBMCage$condition="healthy" # since test data have this info, we add it also to the train data
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(training_PBMCage, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
training_PBMCage$CD8T.naive=score_df
# extract the topN gene expression
training_PBMCage_expr=FetchData(training_PBMCage, vars=c("age","sex","CD8T.naive","condition",topN_genes), layer="counts")

## The testing data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
test_sle$donor_id=test_sle$subject_id
test_sle=AggregateExpression(test_sle, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id","classification"))
test_sle$age=as.numeric(unlist(lapply(strsplit(colnames(test_sle), split="_"), function(x) x[1])))
test_sle$donor_id=unlist(lapply(strsplit(colnames(test_sle), split="_"), function(x) x[2]))
test_sle$sex="unknown" # since train data have this info, we add it also to the test data
test_sle$condition=unlist(lapply(strsplit(colnames(test_sle), split="_"), function(x) x[3]))
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr=FetchData(test_sle, vars=c("age","sex","CD8T.naive","condition", topN_genes), layer="counts")

### Synchronize the genes occurred in both the training and the testing dataset
# in case gene symbols here do not match in the training and the testing dataset
feature.exist.In_train=c(colnames(training_PBMCage_expr), gsub("-","_",colnames(training_PBMCage_expr)), gsub("_","-",colnames(training_PBMCage_expr)))
feature.exist.In_test=c(colnames(test_sle_expr), gsub("-","_",colnames(test_sle_expr)), gsub("_","-",colnames(test_sle_expr)))
coexist.genes=intersect(feature.exist.In_train, feature.exist.In_test)
# update the training and the testing dataset
training_PBMCage_expr=training_PBMCage_expr %>% select(any_of(coexist.genes))
test_sle_expr=test_sle_expr %>% select(any_of(coexist.genes))
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr.t=data.table::transpose(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
colnames(training_PBMCage_expr.t)=rownames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
rownames(training_PBMCage_expr.t)=colnames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
test_sle_expr.t=data.table::transpose(test_sle_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
colnames(test_sle_expr.t)=rownames(test_sle_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
rownames(test_sle_expr.t)=colnames(test_sle_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive","condition"))))

### Integrate the training and the testing datasets
training_obj=CreateSeuratObject(counts=training_PBMCage_expr.t, min.cells=0, min.features=0, meta.data=training_PBMCage_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
training_obj$dataset="train"
testing_obj=CreateSeuratObject(counts=test_sle_expr.t, min.cells=0, min.features=0, meta.data=test_sle_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive","condition"))))
testing_obj$dataset="test"
train_test.mergedobj=merge(training_obj, testing_obj, add.cell.ids=c("train","test"))
# take the highest acceptable npcs_dims and kweight
npcs_dims=ifelse(ncol(training_obj)>1000, 100, min(ncol(training_obj),nrow(training_obj)))
kweight=min(ncol(training_obj), ncol(testing_obj))
# integrate
train_test.mergedobj_integrated=train_test.mergedobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs=npcs_dims) %>%
  FindNeighbors(dims=1:(npcs_dims-1)) %>%
  FindClusters(resolution=0.5)
train_test.mergedobj_integrated=train_test.mergedobj_integrated %>% 
  IntegrateLayers(method=RPCAIntegration, k.weight=30, orig.reduction="pca", new.reduction="integrated.rpca")

### Extract the integrated.rpca for model training
# !: split the test data to get healthy and SLE
# !: because the test data do not have sex info, so remove this feature
# training data
training_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="train")@reductions$integrated.rpca@cell.embeddings)
training_rpca.df=t(scale(t(as.matrix(training_rpca.df)))) %>% as.data.frame()
training_rpca.df$age=training_obj$age
training_rpca.df$CD8T.naive=training_obj$CD8T.naive
# testing data
testing_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="test" & condition=="Control")@reductions$integrated.rpca@cell.embeddings)
testing_rpca.df=t(scale(t(as.matrix(testing_rpca.df)))) %>% as.data.frame()
testing_rpca.df$age=subset(testing_obj, condition=="Control")$age
testing_rpca.df$CD8T.naive=subset(testing_obj, condition=="Control")$CD8T.naive

### Test only on the ages that have occurred more than 10 times in the training data
ggplot2::ggplot(data.frame(no=1:length(training_rpca.df$age), age=training_rpca.df$age), ggplot2::aes(x=age)) +
  ggplot2::geom_density()
# testing_rpca.df=testing_rpca.df %>% subset(age %in% c(60:80))
# possible_prediction=as.numeric(names(table(training_rpca.df$age)[table(training_rpca.df$age)>15]))

### Create the final train and test sets
# train set
train.x=training_rpca.df %>% select(-age)
train.y=training_rpca.df %>% select(age) %>% unlist() %>% unname()
train.xANDy=cbind(train.x, age=train.y)
# test set
test.x=testing_rpca.df %>% select(-age)
test.y=testing_rpca.df %>% select(age) %>% unlist() %>% unname()

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]
# calculate parameters
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE) # check

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
RESULT_DFs[[2]]=result_df
names(RESULT_DFs)[2]="SLECtrl"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")

#####################################



### Check the performance with Immunity dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

### Load the PBMCage dataset
# pseudobulk on donor_id
training_PBMCage=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
training_PBMCage=AggregateExpression(training_PBMCage, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
training_PBMCage$age=as.numeric(unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[1])))
training_PBMCage$sex=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[2]))
training_PBMCage$donor_id=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[3]))
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(training_PBMCage, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
training_PBMCage$CD8T.naive=score_df
# extract the topN gene expression
training_PBMCage_expr=FetchData(training_PBMCage, vars=c("age","sex","CD8T.naive",topN_genes), layer="counts")

## The testing data
test_dataset=readRDS("~/Project_PBMCage/Immunity/Immunity_all_Pseudobulk_onTube.rds")
head(test_dataset[[]]) # check
test_dataset[[]]=test_dataset[[]] %>% mutate(sex=ifelse(sex=="Male","male","female"))
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_dataset, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_dataset$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_dataset_expr=FetchData(test_dataset, vars=c("age","sex","CD8T.naive", topN_genes), layer="counts")

### Synchronize the genes occurred in both the training and the testing dataset
# in case gene symbols here do not match in the training and the testing dataset
feature.exist.In_train=c(colnames(training_PBMCage_expr), gsub("-","_",colnames(training_PBMCage_expr)), gsub("_","-",colnames(training_PBMCage_expr)))
feature.exist.In_test=c(colnames(test_dataset_expr), gsub("-","_",colnames(test_dataset_expr)), gsub("_","-",colnames(test_dataset_expr)))
coexist.genes=intersect(feature.exist.In_train, feature.exist.In_test)
# update the training and the testing dataset
training_PBMCage_expr=training_PBMCage_expr %>% select(any_of(coexist.genes))
test_dataset_expr=test_dataset_expr %>% select(any_of(coexist.genes))
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr.t=data.table::transpose(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
colnames(training_PBMCage_expr.t)=rownames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
rownames(training_PBMCage_expr.t)=colnames(training_PBMCage_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
test_dataset_expr.t=data.table::transpose(test_dataset_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
colnames(test_dataset_expr.t)=rownames(test_dataset_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))
rownames(test_dataset_expr.t)=colnames(test_dataset_expr %>% select(which(! colnames(.) %in% c("age","sex","CD8T.naive"))))

### Integrate the training and the testing datasets
training_obj=CreateSeuratObject(counts=training_PBMCage_expr.t, min.cells=0, min.features=0, meta.data=training_PBMCage_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive"))))
training_obj$dataset="train"
testing_obj=CreateSeuratObject(counts=test_dataset_expr.t, min.cells=0, min.features=0, meta.data=test_dataset_expr %>% select(which(colnames(.) %in% c("age","sex","CD8T.naive"))))
testing_obj$dataset="test"
train_test.mergedobj=merge(training_obj, testing_obj, add.cell.ids=c("train","test"))
# take the highest acceptable npcs_dims and kweight
npcs_dims=ifelse(ncol(training_obj)>1000, 100, min(ncol(training_obj),nrow(training_obj)))
npcs_dims=100
kweight=min(ncol(training_obj), ncol(testing_obj))
kweight=30
# integrate
train_test.mergedobj_integrated=train_test.mergedobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs=npcs_dims) %>%
  FindNeighbors(dims=1:(npcs_dims-1)) %>%
  FindClusters(resolution=0.5)
train_test.mergedobj_integrated=train_test.mergedobj_integrated %>% 
  IntegrateLayers(method=RPCAIntegration, k.weight=kweight, orig.reduction="pca", new.reduction="integrated.rpca")

### Extract the integrated.rpca for model training
# training data
training_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="train")@reductions$integrated.rpca@cell.embeddings)
training_rpca.df=t(scale(t(as.matrix(training_rpca.df)))) %>% as.data.frame()
training_rpca.df$age=training_obj$age
training_rpca.df$CD8T.naive=training_obj$CD8T.naive
training_rpca.df$sex=training_obj$sex
training_rpca.df$sex=as.factor(training_rpca.df$sex)
levels(training_rpca.df$sex)
levels(training_rpca.df$sex)=c("female","male")
# testing data
testing_rpca.df=as.data.frame(subset(train_test.mergedobj_integrated, dataset=="test")@reductions$integrated.rpca@cell.embeddings)
testing_rpca.df=t(scale(t(as.matrix(testing_rpca.df)))) %>% as.data.frame()
testing_rpca.df$age=testing_obj$age
testing_rpca.df$CD8T.naive=testing_obj$CD8T.naive
testing_rpca.df$sex=testing_obj$sex
testing_rpca.df$sex=as.factor(testing_rpca.df$sex)
levels(testing_rpca.df$sex)
levels(testing_rpca.df$sex)=c("female","male")

### Test only on the ages that have occurred more than 10 times in the training data
ggplot2::ggplot(data.frame(no=1:length(training_rpca.df$age), age=training_rpca.df$age), ggplot2::aes(x=age)) +
  ggplot2::geom_density()
# testing_rpca.df=testing_rpca.df %>% subset(age %in% c(60:80))
# possible_prediction=as.numeric(names(table(training_rpca.df$age)[table(training_rpca.df$age)>15]))

### Create the final train and test sets
# train set
train.x=training_rpca.df %>% select(-age)
train.y=training_rpca.df %>% select(age) %>% unlist() %>% unname()
train.xANDy=cbind(train.x, age=train.y)
# test set
test.x=testing_rpca.df %>% select(-age)
test.y=testing_rpca.df %>% select(age) %>% unlist() %>% unname()

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]
# calculate parameters
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE) # check

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
RESULT_DFs[[3]]=result_df
names(RESULT_DFs)[3]="Immunity"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")

#####################################



### Check the performance with PBMCage dataset (cross-validation)
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

### Load the PBMCage dataset
# pseudobulk on donor_id
training_PBMCage=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
training_PBMCage=AggregateExpression(training_PBMCage, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
training_PBMCage$age=as.numeric(unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[1])))
training_PBMCage$sex=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[2]))
training_PBMCage$donor_id=unlist(lapply(strsplit(colnames(training_PBMCage), split="_"), function(x) x[3]))
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(training_PBMCage, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
training_PBMCage$CD8T.naive=score_df
# extract the topN gene expression
training_PBMCage_expr=FetchData(training_PBMCage, vars=c("age","sex","CD8T.naive",topN_genes), layer="counts")
training_PBMCage_expr=cbind(training_PBMCage_expr[,1:2], t(scale(t(training_PBMCage_expr[,3:ncol(training_PBMCage_expr)]))))

### Prepare the training and testing data
Data.wide=training_PBMCage_expr %>% dplyr::rename(Age=age, Sex=sex)
# make sure that the colnames are ok for modeling
colnames(Data.wide)=gsub("-","_",colnames(Data.wide)) # for genes
# split the dataset to train/test
set_trainset_prob=0.8
ind=sample(c(TRUE, FALSE), nrow(Data.wide), replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
train.xANDy=Data.wide[ind, ]
test.xANDy=Data.wide[!ind, ]
# prepare X and y
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
# calculate parameters
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
RESULT_DFs[[4]]=result_df
names(RESULT_DFs)[4]="PBMCage"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")

#####################################



### Save the performance of the PBMCage_RNA-based model with all the predictive range
#####################################
### 

### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
file_names=names(RESULT_DFs)
dataset_names=c("Combination of GSE213516, GSE158055, SC cohort",
                "GSE189050",
                "Terekhova et al. (2023) dataset",
                "Yazar et al. (2022) dataset")

plot_list_=RSQs=MAEs=MADs=RMSEs=list()
for (i in 1:length(RESULT_DFs)) {
  df_=RESULT_DFs[[i]]
  train.y=df_$train.y
  pred=df_$pred
  test.y=df_$test.y
  ggplot_data=data.frame(Predicted.Age=pred, Chronological.Age=test.y)
  
  # attach MAE, MAD and RMSE
  SST=sum((train.y-mean(train.y))^2)
  SSE=sum((pred-test.y)^2)
  rsq=1-SSE/SST
  MAE=Metrics::mae(pred, test.y)
  MAD=mad(pred, test.y)
  RMSE=Metrics::rmse(pred, test.y)
  
  RSQs[[i]]=rsq; MAEs[[i]]=MAE; MADs[[i]]=MAD; RMSEs[[i]]=RMSE
}

RESULT_=data.frame(rsq=unlist(RSQs), MAE=unlist(MAEs), MAD=unlist(MADs), RMSE=unlist(RMSEs), dataset=dataset_names, geneset="AgeGene5K", info="Xu_5000.ranked")
write.table(RESULT_, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.txt", sep="\t")

#####################################



### Save the performance of the PBMCage_RNA-based model with range 60-80
#####################################
### 

### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.rds")
file_names=names(RESULT_DFs)
dataset_names=c("Combination of GSE213516, GSE158055, SC cohort",
                "GSE189050",
                "Terekhova et al. (2023) dataset",
                "Yazar et al. (2022) dataset")

plot_list_=RSQs=MAEs=MADs=RMSEs=list()
for (i in 1:length(RESULT_DFs)) {
  df_=RESULT_DFs[[i]]
  train.y=df_$train.y
  pred=df_$pred
  test.y=df_$test.y
  ggplot_data=data.frame(Predicted.Age=pred, Chronological.Age=test.y)
  
  # choose the predictive range
  filter_60_80=(test.y>=60 & test.y<=80)
  pred=pred[filter_60_80]
  test.y=test.y[filter_60_80]
  ggplot_data_limited=ggplot_data[filter_60_80,]
  
  # attach MAE, MAD and RMSE
  SST=sum((train.y-mean(train.y))^2)
  SSE=sum((pred-test.y)^2)
  rsq=1-SSE/SST
  MAE=Metrics::mae(pred, test.y)
  MAD=mad(pred, test.y)
  RMSE=Metrics::rmse(pred, test.y)
  
  RSQs[[i]]=rsq; MAEs[[i]]=MAE; MADs[[i]]=MAD; RMSEs[[i]]=RMSE
}

RESULT_=data.frame(rsq=unlist(RSQs), MAE=unlist(MAEs), MAD=unlist(MADs), RMSE=unlist(RMSEs), dataset=dataset_names, geneset="AgeGene5K", info="Xu_5000.ranked")
write.table(RESULT_, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results_60To80.txt", sep="\t")

#####################################






###############################################################################################################
##################################### Train model based on the PBMCage celltype frequency #####################################
###############################################################################################################

RESULT_DFs=list()
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")

#####################################



### Freq+GeneExpr model cross-validation on PBMCage dataset (cross-validation)
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Get the celltype frequency
pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=pbmc.seu_merged$donor_id
Age_=pbmc.seu_merged$age
Sex_=pbmc.seu_merged$sex
CellType_=pbmc.seu_merged$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
Data.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression
# use the pseudobulk obj
pseudobulkobj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
pseudobulkobj_donor=AggregateExpression(pseudobulkobj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
pseudobulkobj_donor$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[1])))
pseudobulkobj_donor$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[2]))
pseudobulkobj_donor$donor_id=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
pseudo_data=FetchData(pseudobulkobj_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
pseudo_data$ID=gsub("-","_",pseudobulkobj_donor$donor_id)
pseudo_data=pseudo_data %>% select(-donor_id)
# scale the data
pseudo_data[,!grepl("ID", colnames(pseudo_data))]=t(scale(t(pseudo_data[,!grepl("ID", colnames(pseudo_data))])))

### Combine freq and expression
Data.wide=dplyr::left_join(Data.wide, pseudo_data, by="ID")

### Prepare the training and testing data
# remove ID
Data.wide=Data.wide %>% select(-ID)
# make sure that the colnames are ok for modeling
colnames(Data.wide)=gsub(" ","\\.",colnames(Data.wide)) # for Annot.rough
colnames(Data.wide)=gsub("-","_",colnames(Data.wide)) # for genes
# convert age to int and sex to char
Data.wide$Age=as.integer(as.character(Data.wide$Age))
Data.wide$Sex=as.character(Data.wide$Sex)
# split the dataset to train/test
set_trainset_prob=0.8
ind=sample(c(TRUE, FALSE), nrow(Data.wide), replace=TRUE, prob=c(set_trainset_prob, 1-set_trainset_prob))
train.xANDy=Data.wide[ind, ]
test.xANDy=Data.wide[!ind, ]
# prepare X and y
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
# calculate parameters
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
RESULT_DFs[[1]]=result_df
names(RESULT_DFs)[1]="PBMCage"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")

#####################################



### Check the performance with the 17_PBMC dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Get the celltype frequency of the PBMCage
PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=PBMCage_obj$donor_id
Age_=PBMCage_obj$age
Sex_=as.factor(as.character(PBMCage_obj$sex)); levels(Sex_)=c("female","male")
CellType_=PBMCage_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
PBMCAGE.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the PBMCage
# use the pseudobulk obj
PBMCage_obj_pseudo=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
PBMCage_obj_pseudo_donor=AggregateExpression(PBMCage_obj_pseudo, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
PBMCage_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[1])))
PBMCage_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[2]))
PBMCage_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
PBMCage_expr_data=FetchData(PBMCage_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
PBMCage_expr_data$ID=gsub("-","_",PBMCage_obj_pseudo_donor$donor_id)
PBMCage_expr_data=PBMCage_expr_data %>% select(-donor_id)
# scale the data
PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))]=t(scale(t(PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))])))
# combine freq and expression
PBMCAGE.wide=dplyr::left_join(PBMCAGE.wide, PBMCage_expr_data, by="ID")

### Get the celltype frequency of the testing dataset
test_obj=readRDS("~/Project_PBMCage/GSE213516_17PBMC/17_pbmc_annotated_AssignedOnly.rds")
ID_=test_obj$donor_id
Age_=test_obj$age
Sex_=as.factor(as.character(test_obj$sex)); levels(Sex_)=c("female","male")
CellType_=test_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
TEST.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the Test dataset
# make the pseudobulk obj
test_obj_pseudo_donor=AggregateExpression(test_obj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
test_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[1])))
test_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[2]))
test_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
test_expr_data=FetchData(test_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
test_expr_data$ID=gsub("-","_",test_obj_pseudo_donor$donor_id)
test_expr_data=test_expr_data %>% select(-donor_id)
# scale the data
test_expr_data[,!grepl("ID", colnames(test_expr_data))]=t(scale(t(test_expr_data[,!grepl("ID", colnames(test_expr_data))])))
# combine freq and expression
TEST.wide=dplyr::left_join(TEST.wide, test_expr_data, by="ID")

### Prepare the training and testing data
# remove ID
PBMCAGE.wide=PBMCAGE.wide %>% select(-ID)
TEST.wide=TEST.wide %>% select(-ID)
# make sure that the colnames are ok for modeling
colnames(PBMCAGE.wide)=gsub(" ","\\.",colnames(PBMCAGE.wide)) # for Annot.rough
colnames(PBMCAGE.wide)=gsub("-","_",colnames(PBMCAGE.wide)) # for genes
colnames(TEST.wide)=gsub(" ","\\.",colnames(TEST.wide)) # for Annot.rough
colnames(TEST.wide)=gsub("-","_",colnames(TEST.wide)) # for genes
# convert age to int and sex to char
PBMCAGE.wide$Age=as.integer(as.character(PBMCAGE.wide$Age))
PBMCAGE.wide$Sex=as.character(PBMCAGE.wide$Sex)
TEST.wide$Age=as.integer(as.character(TEST.wide$Age))
TEST.wide$Sex=as.character(TEST.wide$Sex)
# synchronize the training and testing dataset features
coexist.features=intersect(colnames(PBMCAGE.wide), colnames(TEST.wide))
PBMCAGE.wide=PBMCAGE.wide %>% select(any_of(coexist.features))
TEST.wide=TEST.wide %>% select(any_of(coexist.features))
# prepare X and y
train.xANDy=PBMCAGE.wide
test.xANDy=TEST.wide
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
RESULT_DFs[[2]]=result_df
names(RESULT_DFs)[2]="17_PBMC"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")

#####################################



### Check the performance with the CovidCtrl dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Get the celltype frequency of the PBMCage
# PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=PBMCage_obj$donor_id
Age_=PBMCage_obj$age
Sex_=as.factor(as.character(PBMCage_obj$sex)); levels(Sex_)=c("female","male")
CellType_=PBMCage_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
PBMCAGE.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the PBMCage
# use the pseudobulk obj
PBMCage_obj_pseudo=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
PBMCage_obj_pseudo_donor=AggregateExpression(PBMCage_obj_pseudo, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
PBMCage_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[1])))
PBMCage_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[2]))
PBMCage_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
PBMCage_expr_data=FetchData(PBMCage_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
PBMCage_expr_data$ID=gsub("-","_",PBMCage_obj_pseudo_donor$donor_id)
PBMCage_expr_data=PBMCage_expr_data %>% select(-donor_id)
# scale the data
PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))]=t(scale(t(PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))])))
# combine freq and expression
PBMCAGE.wide=dplyr::left_join(PBMCAGE.wide, PBMCage_expr_data, by="ID")

### Get the celltype frequency of the testing dataset
test_obj=readRDS("~/Project_PBMCage/GSE158055/frozen_fresh_pbmc_INTEGRATED.rds")
ID_=test_obj$donor_id
Age_=test_obj$age
Sex_=as.factor(as.character(test_obj$sex)); levels(Sex_)=c("female","male")
CellType_=test_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
TEST.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the Test dataset
# make the pseudobulk obj
test_obj_pseudo_donor=AggregateExpression(test_obj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
test_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[1])))
test_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[2]))
test_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
test_expr_data=FetchData(test_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
test_expr_data$ID=test_obj_pseudo_donor$donor_id
test_expr_data=test_expr_data %>% select(-donor_id)
# scale the data
test_expr_data[,!grepl("ID", colnames(test_expr_data))]=t(scale(t(test_expr_data[,!grepl("ID", colnames(test_expr_data))])))
# combine freq and expression
TEST.wide=dplyr::left_join(TEST.wide, test_expr_data, by="ID")

### Prepare the training and testing data
# remove ID
PBMCAGE.wide=PBMCAGE.wide %>% select(-ID)
TEST.wide=TEST.wide %>% select(-ID)
# make sure that the colnames are ok for modeling
colnames(PBMCAGE.wide)=gsub(" ","\\.",colnames(PBMCAGE.wide)) # for Annot.rough
colnames(PBMCAGE.wide)=gsub("-","_",colnames(PBMCAGE.wide)) # for genes
colnames(TEST.wide)=gsub(" ","\\.",colnames(TEST.wide)) # for Annot.rough
colnames(TEST.wide)=gsub("-","_",colnames(TEST.wide)) # for genes
# convert age to int and sex to char
PBMCAGE.wide$Age=as.integer(as.character(PBMCAGE.wide$Age))
PBMCAGE.wide$Sex=as.character(PBMCAGE.wide$Sex)
TEST.wide$Age=as.integer(as.character(TEST.wide$Age))
TEST.wide$Sex=as.character(TEST.wide$Sex)
# synchronize the training and testing dataset features
coexist.features=intersect(colnames(PBMCAGE.wide), colnames(TEST.wide))
PBMCAGE.wide=PBMCAGE.wide %>% select(any_of(coexist.features))
TEST.wide=TEST.wide %>% select(any_of(coexist.features))
# prepare X and y
train.xANDy=PBMCAGE.wide
test.xANDy=TEST.wide
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
RESULT_DFs[[3]]=result_df
names(RESULT_DFs)[3]="CovidCtrl"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")

#####################################



### Check the performance with the JapanSC dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Get the celltype frequency of the PBMCage
# PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=PBMCage_obj$donor_id
Age_=PBMCage_obj$age
Sex_=as.factor(as.character(PBMCage_obj$sex)); levels(Sex_)=c("female","male")
CellType_=PBMCage_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
PBMCAGE.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the PBMCage
# use the pseudobulk obj
PBMCage_obj_pseudo=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
PBMCage_obj_pseudo_donor=AggregateExpression(PBMCage_obj_pseudo, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
PBMCage_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[1])))
PBMCage_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[2]))
PBMCage_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
PBMCage_expr_data=FetchData(PBMCage_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
PBMCage_expr_data$ID=gsub("-","_",PBMCage_obj_pseudo_donor$donor_id)
PBMCage_expr_data=PBMCage_expr_data %>% select(-donor_id)
# scale the data
PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))]=t(scale(t(PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))])))
# combine freq and expression
PBMCAGE.wide=dplyr::left_join(PBMCAGE.wide, PBMCage_expr_data, by="ID")

### Get the celltype frequency of the testing dataset
test_obj=readRDS("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated.rds")
ID_=test_obj$donor_id
Age_=test_obj$age
Sex_=as.factor(as.character(test_obj$sex)); levels(Sex_)=c("female","male")
CellType_=test_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
TEST.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the Test dataset
# make the pseudobulk obj
test_obj_pseudo_donor=AggregateExpression(test_obj, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
test_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[1])))
test_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[2]))
test_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(test_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
test_expr_data=FetchData(test_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
test_expr_data$ID=gsub("-","_",test_obj_pseudo_donor$donor_id)
test_expr_data=test_expr_data %>% select(-donor_id)
# scale the data
test_expr_data[,!grepl("ID", colnames(test_expr_data))]=t(scale(t(test_expr_data[,!grepl("ID", colnames(test_expr_data))])))
# combine freq and expression
TEST.wide=dplyr::left_join(TEST.wide, test_expr_data, by="ID")

### Prepare the training and testing data
# remove ID
PBMCAGE.wide=PBMCAGE.wide %>% select(-ID)
TEST.wide=TEST.wide %>% select(-ID)
# make sure that the colnames are ok for modeling
colnames(PBMCAGE.wide)=gsub(" ","\\.",colnames(PBMCAGE.wide)) # for Annot.rough
colnames(PBMCAGE.wide)=gsub("-","_",colnames(PBMCAGE.wide)) # for genes
colnames(TEST.wide)=gsub(" ","\\.",colnames(TEST.wide)) # for Annot.rough
colnames(TEST.wide)=gsub("-","_",colnames(TEST.wide)) # for genes
# convert age to int and sex to char
PBMCAGE.wide$Age=as.integer(as.character(PBMCAGE.wide$Age))
PBMCAGE.wide$Sex=as.character(PBMCAGE.wide$Sex)
TEST.wide$Age=as.integer(as.character(TEST.wide$Age))
TEST.wide$Sex=as.character(TEST.wide$Sex)
# synchronize the training and testing dataset features
coexist.features=intersect(colnames(PBMCAGE.wide), colnames(TEST.wide))
PBMCAGE.wide=PBMCAGE.wide %>% select(any_of(coexist.features))
TEST.wide=TEST.wide %>% select(any_of(coexist.features))
# prepare X and y
train.xANDy=PBMCAGE.wide
test.xANDy=TEST.wide
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
RESULT_DFs[[4]]=result_df
names(RESULT_DFs)[4]="JapanSC"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")


#####################################



### Check the performance with the Immunity dataset
#####################################
###

library(UCell)
library(Seurat)
library(pls)

### Get the celltype frequency of the PBMCage
# PBMCage_obj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=PBMCage_obj$donor_id
Age_=PBMCage_obj$age
Sex_=as.factor(as.character(PBMCage_obj$sex)); levels(Sex_)=c("female","male")
CellType_=PBMCage_obj$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
PBMCAGE.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the PBMCage
# use the pseudobulk obj
PBMCage_obj_pseudo=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
PBMCage_obj_pseudo_donor=AggregateExpression(PBMCage_obj_pseudo, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id"))
PBMCage_obj_pseudo_donor$age=as.numeric(unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[1])))
PBMCage_obj_pseudo_donor$sex=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[2]))
PBMCage_obj_pseudo_donor$donor_id=unlist(lapply(strsplit(colnames(PBMCage_obj_pseudo_donor), split="_"), function(x) x[3]))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
PBMCage_expr_data=FetchData(PBMCage_obj_pseudo_donor, vars=c("donor_id",Xu_5000.ranked), layer="data")
PBMCage_expr_data$ID=gsub("-","_",PBMCage_obj_pseudo_donor$donor_id)
PBMCage_expr_data=PBMCage_expr_data %>% select(-donor_id)
# scale the data
PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))]=t(scale(t(PBMCage_expr_data[,!grepl("ID", colnames(PBMCage_expr_data))])))
# combine freq and expression
PBMCAGE.wide=dplyr::left_join(PBMCAGE.wide, PBMCage_expr_data, by="ID")

### Get the celltype frequency of the testing dataset
merged_meta1=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor_metadata.csv.gz", sep="\t")
merged_meta2=data.table::fread("~/Project_PBMCage/Immunity/Immunity_RawCountObj_rest_metadata.csv.gz", sep="\t")
merged_meta=rbind(merged_meta1, merged_meta2)
ID_=merged_meta$Tube_id
Age_=merged_meta$Age
Sex_=as.factor(as.character(merged_meta$Sex)); levels(Sex_)=c("female","male"); Sex_=as.character(Sex_)
CellType_=merged_meta$Annot.rough
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_, Sex=Sex_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data$Age=ID_Age_match$Age[match(data$ID, ID_Age_match$ID)]
data$Sex=ID_Age_match$Sex[match(data$ID, ID_Age_match$ID)]
data$CellType=as.factor(data$CellType)
TEST.wide=tidyr::pivot_wider(data, names_from="CellType", values_from="Freq")

### Add expression to the Test dataset
# make the pseudobulk obj
test_obj_pseudo_donor=readRDS("~/Project_PBMCage/Immunity/Immunity_all_Pseudobulk_onTube.rds")
head(test_obj_pseudo_donor[[]]) # check
test_obj_pseudo_donor[[]]=test_obj_pseudo_donor[[]] %>% mutate(sex=ifelse(sex=="Male","male","female"))
# extract the genes
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.ranked=gene_list$Xu_5000.ranked
Xu_5000.ranked=Xu_5000.ranked[1:300] # take topN genes
# fetch the data
test_expr_data=FetchData(test_obj_pseudo_donor, vars=c("tube_id",Xu_5000.ranked), layer="data")
test_expr_data$ID=gsub("-","_",test_obj_pseudo_donor$tube_id)
test_expr_data=test_expr_data %>% select(-tube_id)
# scale the data
test_expr_data[,!grepl("ID", colnames(test_expr_data))]=t(scale(t(test_expr_data[,!grepl("ID", colnames(test_expr_data))])))
# combine freq and expression
TEST.wide=dplyr::left_join(TEST.wide, test_expr_data, by="ID")

### Prepare the training and testing data
# remove ID
PBMCAGE.wide=PBMCAGE.wide %>% select(-ID)
TEST.wide=TEST.wide %>% select(-ID)
# make sure that the colnames are ok for modeling
colnames(PBMCAGE.wide)=gsub(" ","\\.",colnames(PBMCAGE.wide)) # for Annot.rough
colnames(PBMCAGE.wide)=gsub("-","_",colnames(PBMCAGE.wide)) # for genes
colnames(TEST.wide)=gsub(" ","\\.",colnames(TEST.wide)) # for Annot.rough
colnames(TEST.wide)=gsub("-","_",colnames(TEST.wide)) # for genes
# convert age to int and sex to char
PBMCAGE.wide$Age=as.integer(as.character(PBMCAGE.wide$Age))
PBMCAGE.wide$Sex=as.character(PBMCAGE.wide$Sex)
TEST.wide$Age=as.integer(as.character(TEST.wide$Age))
TEST.wide$Sex=as.character(TEST.wide$Sex)
# synchronize the training and testing dataset features
coexist.features=intersect(colnames(PBMCAGE.wide), colnames(TEST.wide))
PBMCAGE.wide=PBMCAGE.wide %>% select(any_of(coexist.features))
TEST.wide=TEST.wide %>% select(any_of(coexist.features))
# prepare X and y
train.xANDy=PBMCAGE.wide
test.xANDy=TEST.wide
train.x=train.xANDy %>% select(-Age)
train.y=train.xANDy$Age
test.x=test.xANDy %>% select(-Age)
test.y=test.xANDy$Age

### Fit the model
all_genes_codes=paste0(colnames(train.x), collapse="+")
eval(parse(text=paste0("Model=plsr(Age~",all_genes_codes,", data=train.xANDy, scale=TRUE, validation='CV')")))
# determine the best number of comp
cverr=RMSEP(Model)$val[1,,]
imin_cv=which.min(cverr)-1
r2=R2(Model)$val[1,,]
imax_r2=which.max(r2)-1
best.ncomp=as.integer(mean(imin_cv, imax_r2))

### Predict
pred=(predict(Model, test.x, ncomp=best.ncomp)[,,1])
SST=sum((train.y-mean(train.y))^2)
SSE=sum((pred-test.y)^2)
rsq=1-SSE/SST
MAE=Metrics::mae(pred, test.y)
MAD=mad(pred, test.y)
RMSE=Metrics::rmse(pred, test.y)
cor(pred, test.y)
c(best.ncomp, rsq, MAE, MAD, RMSE)

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
result_df=list(train.y=train.y,
               test.y=test.y,
               pred=unname(pred))
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
RESULT_DFs[[5]]=result_df
names(RESULT_DFs)[5]="Immunity"
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")

#####################################



### Plot the performance of the PBMCage_RNA freq+GeneExpr model
#####################################
### 

### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
file_names=names(RESULT_DFs)
dataset_names=c("Yazar et al. (2022) dataset",
                "GSE213516",
                "GSE158055",
                "SC cohort",
                "Immunity")
file_names=names(RESULT_DFs)

plot_list_=RSQs=MAEs=MADs=RMSEs=list()
for (i in 1:length(RESULT_DFs)) {
  df_=RESULT_DFs[[i]]
  train.y=df_$train.y
  pred=df_$pred
  test.y=df_$test.y
  ggplot_data=data.frame(Predicted.Age=pred, Chronological.Age=test.y)
  
  # attach MAE, MAD and RMSE
  SST=sum((train.y-mean(train.y))^2)
  SSE=sum((pred-test.y)^2)
  rsq=1-SSE/SST
  MAE=Metrics::mae(pred, test.y)
  MAD=mad(pred, test.y)
  RMSE=Metrics::rmse(pred, test.y)
  
  RSQs[[i]]=rsq; MAEs[[i]]=MAE; MADs[[i]]=MAD; RMSEs[[i]]=RMSE
}

RESULT_=data.frame(rsq=unlist(RSQs), MAE=unlist(MAEs), MAD=unlist(MADs), RMSE=unlist(RMSEs), dataset=dataset_names, geneset="AgeGene5K", info="Xu_5000.ranked")
write.table(RESULT_, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.txt", sep="\t")

#####################################



### Save the performance of the PBMCage_RNA freq+GeneExpr model on range 60-80
#####################################
### 


### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.rds")
file_names=names(RESULT_DFs)
dataset_names=c("Yazar et al. (2022) dataset",
                "GSE213516",
                "GSE158055",
                "SC cohort",
                "Immunity")
file_names=names(RESULT_DFs)

plot_list_=RSQs=MAEs=MADs=RMSEs=list()
for (i in 1:length(RESULT_DFs)) {
  df_=RESULT_DFs[[i]]
  train.y=df_$train.y
  pred=df_$pred
  test.y=df_$test.y
  ggplot_data=data.frame(Predicted.Age=pred, Chronological.Age=test.y)
  
  # choose the predictive range
  filter_60_80=(test.y>=60 & test.y<=80)
  pred=pred[filter_60_80]
  test.y=test.y[filter_60_80]
  ggplot_data_limited=ggplot_data[filter_60_80,]
  
  # attach MAE, MAD and RMSE
  SST=sum((train.y-mean(train.y))^2)
  SSE=sum((pred-test.y)^2)
  rsq=1-SSE/SST
  MAE=Metrics::mae(pred, test.y)
  MAD=mad(pred, test.y)
  RMSE=Metrics::rmse(pred, test.y)
  
  RSQs[[i]]=rsq; MAEs[[i]]=MAE; MADs[[i]]=MAD; RMSEs[[i]]=RMSE
}

RESULT_=data.frame(rsq=unlist(RSQs), MAE=unlist(MAEs), MAD=unlist(MADs), RMSE=unlist(RMSEs), dataset=dataset_names, geneset="AgeGene5K", info="Xu_5000.ranked")
write.table(RESULT_, "~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results_60TO80.txt", sep="\t")

#####################################



### Combine the performance of the RNAexpr-based model and the PBMCage_RNA freq+GeneExpr model
#####################################
### 

rna_based=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results.txt", sep="\t")
rna_based$model="RNA-based"
rna_freq_based=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results.txt", sep="\t")
rna_freq_based$model="RNA.Freq-based"
both_models=data.table::rbindlist(list(rna_based, rna_freq_based))
write.table(both_models, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_AgeGene5K_results_20TO100.txt", sep="\t")

# age range 60-80 only
rna_based=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_results_60To80.txt", sep="\t")
rna_based$model="RNA-based"
rna_freq_based=read.delim("~/Project_PBMCage/For_or_From_Xu/Results/Model.freq_prediction_results_60TO80.txt", sep="\t")
rna_freq_based$model="RNA.Freq-based"
both_models=data.table::rbindlist(list(rna_based, rna_freq_based))
write.table(both_models, "~/Project_PBMCage/For_or_From_Xu/Results/Model_prediction_AgeGene5K_results_60TO80.txt", sep="\t")

#####################################





