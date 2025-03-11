
##########################################################################

### Apply AgeGene5K-based model on SLE dataset

##########################################################################


RESULT_DFs=list()
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")

### RNA-based training with the SLE Ctrl at the donor level
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

## Process the SLE data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
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
test_sle_expr_total=FetchData(test_sle, vars=c("age","sex","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("Control",rownames(test_sle_expr_total)),]
test_sle_expr=test_sle_expr_total[grepl("SLE",rownames(test_sle_expr_total)),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, sex, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, sex, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=training_PBMCage_expr$condition
test.y_withname=test.y; names(test.y_withname)=test_sle_expr$condition
pred_withname=pred; names(pred_withname)=test_sle_expr$condition
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")
RESULT_DFs$`RNA-based`[["SLE_donorid"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")

#####################################



### RNA-based training with the SLE Ctrl at the rough level
#####################################
###

### Prepare the annot.rough
test_obj=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
unique(test_obj$fine_cell_type)
df=data.frame(donor_id=test_obj$subject_id, condition=test_obj$classification, age=test_obj$age, celltypes=test_obj$fine_cell_type) %>%
  mutate(celltypes=as.character(celltypes)) %>%
  mutate(Annot.inter=ifelse(celltypes=="Progenitors","Other.progenitor",
                            ifelse(celltypes=="CD14+ monocytes","Mono.classical",
                                   ifelse(celltypes=="CD4+ T cells","CD4T cells",
                                          ifelse(celltypes=="cDCs","DC.cDC",
                                                 ifelse(celltypes=="pDCs","DC.pDC",
                                                        ifelse(celltypes=="CD16+ monocytes","Mono.nonclassical",
                                                               ifelse(celltypes=="Naive B cells","B.naive",
                                                                      ifelse(celltypes=="Transitional B cells","B.trans",
                                                                             ifelse(celltypes=="Memory B cells","B.mem",
                                                                                    ifelse(celltypes=="Plasmablasts","B.pbpc",
                                                                                           ifelse(celltypes=="CD8+ cytotoxic T cells","CD8T.CTL",
                                                                                                  ifelse(celltypes=="CD8+ T cells","CD8T.naive/mem",
                                                                                                         ifelse(celltypes=="ABCs","B.ABC", celltypes))))))))))))))
# check
table(test_obj$fine_cell_type)
table(df$Annot.inter, useNA="ifany")
df=df %>% 
  mutate(Annot.rough=ifelse(grepl("Other\\.",Annot.inter),"Other cells",
                            ifelse(grepl("Mono\\.",Annot.inter),"Monocytes",
                                   ifelse(grepl("CD4T",Annot.inter),"CD4T cells",
                                          ifelse(grepl("DC\\.",Annot.inter),"DCs",
                                                 ifelse(grepl("B\\.",Annot.inter),"B cells",
                                                        ifelse(grepl("CD8T\\.",Annot.inter),"CD8T cells", Annot.inter)))))))
# check
table(df$Annot.rough, useNA="ifany")
table(rownames(test_obj[[]])==rownames(df))
# add the celltype annotation
test_obj[[]]$donor_id=test_obj[[]]$subject_id
test_obj[[]]$Annot.inter=df$Annot.inter
test_obj[[]]$Annot.rough=df$Annot.rough

### Make pseudobulk at the rough level
test_obj_pseudo_rough=AggregateExpression(test_obj, assays="RNA", return.seurat=T, group.by=c("donor_id","Annot.rough"))
df_to_match=test_obj_pseudo_rough[[]]
df_to_match=dplyr::left_join(test_obj_pseudo_rough[[]] %>% mutate(donor_id=gsub("_.*","",rownames(.)),
                                                                  Annot.rough=gsub(".*_","",rownames(.))), 
                             df %>% select(-c(celltypes, Annot.inter)) %>% subset(!duplicated(.)), by=c("donor_id","Annot.rough"))
test_obj_pseudo_rough[[]]$donor_id=as.character(df_to_match$donor_id)
test_obj_pseudo_rough[[]]$age=as.numeric(as.character(df_to_match$age))
test_obj_pseudo_rough[[]]$condition=as.character(df_to_match$condition)
test_obj_pseudo_rough[[]]$Annot.rough=as.character(df_to_match$Annot.rough)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

## Process the SLE data
test_sle=test_obj_pseudo_rough
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("age","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("Control",test_sle_expr_total$condition),]
test_sle_expr=test_sle_expr_total[grepl("SLE",test_sle_expr_total$condition),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=paste0(gsub(".*_","",rownames(training_rpca.df)),"***",training_PBMCage_expr$condition)
test.y_withname=test.y; names(test.y_withname)=paste0(gsub(".*_","",rownames(testing_rpca.df)),"***",test_sle_expr$condition)
pred_withname=pred; names(pred_withname)=paste0(gsub(".*_","",names(pred)),"***",test_sle_expr$condition)
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")
RESULT_DFs$`RNA-based`[["SLE_rough"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")

#####################################



### RNA-based training with the SLE Ctrl at the inter level
#####################################
###

### Prepare the annot.rough
test_obj=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
unique(test_obj$fine_cell_type)
df=data.frame(donor_id=test_obj$subject_id, condition=test_obj$classification, age=test_obj$age, celltypes=test_obj$fine_cell_type) %>%
  mutate(celltypes=as.character(celltypes)) %>%
  mutate(Annot.inter=ifelse(celltypes=="Progenitors","Other.progenitor",
                            ifelse(celltypes=="CD14+ monocytes","Mono.classical",
                                   ifelse(celltypes=="CD4+ T cells","CD4T cells",
                                          ifelse(celltypes=="cDCs","DC.cDC",
                                                 ifelse(celltypes=="pDCs","DC.pDC",
                                                        ifelse(celltypes=="CD16+ monocytes","Mono.nonclassical",
                                                               ifelse(celltypes=="Naive B cells","B.naive",
                                                                      ifelse(celltypes=="Transitional B cells","B.trans",
                                                                             ifelse(celltypes=="Memory B cells","B.mem",
                                                                                    ifelse(celltypes=="Plasmablasts","B.pbpc",
                                                                                           ifelse(celltypes=="CD8+ cytotoxic T cells","CD8T.CTL",
                                                                                                  ifelse(celltypes=="CD8+ T cells","CD8T.naive/mem",
                                                                                                         ifelse(celltypes=="ABCs","B.ABC", celltypes))))))))))))))
# check
table(test_obj$fine_cell_type)
table(df$Annot.inter, useNA="ifany")
df=df %>% 
  mutate(Annot.rough=ifelse(grepl("Other\\.",Annot.inter),"Other cells",
                            ifelse(grepl("Mono\\.",Annot.inter),"Monocytes",
                                   ifelse(grepl("CD4T",Annot.inter),"CD4T cells",
                                          ifelse(grepl("DC\\.",Annot.inter),"DCs",
                                                 ifelse(grepl("B\\.",Annot.inter),"B cells",
                                                        ifelse(grepl("CD8T\\.",Annot.inter),"CD8T cells", Annot.inter)))))))
# check
table(df$Annot.rough, useNA="ifany")
table(rownames(test_obj[[]])==rownames(df))
# add the celltype annotation
test_obj[[]]$donor_id=test_obj[[]]$subject_id
test_obj[[]]$Annot.inter=df$Annot.inter
test_obj[[]]$Annot.rough=df$Annot.rough

### Make pseudobulk at the rough level
test_obj_pseudo_inter=AggregateExpression(test_obj, assays="RNA", return.seurat=T, group.by=c("donor_id","Annot.inter"))
df_to_match=test_obj_pseudo_rough[[]]
df_to_match=dplyr::left_join(test_obj_pseudo_inter[[]] %>% mutate(donor_id=gsub("_.*","",rownames(.)),
                                                                  Annot.inter=gsub(".*_","",rownames(.))), 
                             df %>% select(-celltypes) %>% subset(!duplicated(.)), by=c("donor_id","Annot.inter"))
test_obj_pseudo_inter[[]]$donor_id=as.character(df_to_match$donor_id)
test_obj_pseudo_inter[[]]$age=as.numeric(as.character(df_to_match$age))
test_obj_pseudo_inter[[]]$condition=as.character(df_to_match$condition)
test_obj_pseudo_inter[[]]$Annot.rough=as.character(df_to_match$Annot.rough)
test_obj_pseudo_inter[[]]$Annot.inter=as.character(df_to_match$Annot.inter)

### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

## Process the SLE data
test_sle=test_obj_pseudo_inter
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("age","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("Control",test_sle_expr_total$condition),]
test_sle_expr=test_sle_expr_total[grepl("SLE",test_sle_expr_total$condition),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=paste0(gsub(".*_","",rownames(training_rpca.df)),"***",training_PBMCage_expr$condition)
test.y_withname=test.y; names(test.y_withname)=paste0(gsub(".*_","",rownames(testing_rpca.df)),"***",test_sle_expr$condition)
pred_withname=pred; names(pred_withname)=paste0(gsub(".*_","",names(pred)),"***",test_sle_expr$condition)
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")
RESULT_DFs$`RNA-based`[["SLE_inter"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")

#####################################


### Analyze the results
#####################################
###

### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_SLE.rds")
RESULT_DFs=RESULT_DFs$`RNA-based`

### At the donor_id level
RESULT_donor=data.frame(condition=names(RESULT_DFs$SLE_donorid$test.y),
                        CA=RESULT_DFs$SLE_donorid$test.y, 
                        pred.condition=names(RESULT_DFs$SLE_donorid$pred),
                        BA=RESULT_DFs$SLE_donorid$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_donor$condition==RESULT_donor$pred.condition)

RESULT_donor=RESULT_donor %>% select(-pred.condition) %>% 
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(gsub("SLE ","",condition), ".", CA_or_BA))
RESULT_donor$condition_train=forcats::fct_relevel(RESULT_donor$condition_train, c("INACT.CA","INACT.BA","ACT.CA","ACT.BA"))
RESULT_donor$condition=forcats::fct_relevel(RESULT_donor$condition, c("SLE INACT","SLE ACT"))

plot_=
  ggplot(RESULT_donor, aes(x=condition_train, y=age, color=condition)) +
  geom_point(size=2) +
  # facet_wrap(~celltypes) +
  scale_color_manual(values=c("lightblue3","coral3")) +
  geom_line(aes(group=pair), linewidth=0.5) +
  ggpubr::stat_compare_means(aes(label=paste0()),
                             comparisons=list(c("INACT.CA","INACT.BA"), c("ACT.CA","ACT.BA")),
                             tip.length=0,
                             paired=T, size=3.5, bracket.size=0, 
                             label.y=65) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9, angle=45, hjust=1, vjust=1),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none")

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inSLE_donor.pdf", height=5, width=2)
plot(plot_)
dev.off()

### At the Annot.rough level
RESULT_rough=data.frame(condition=names(RESULT_DFs$SLE_rough$test.y),
                        CA=RESULT_DFs$SLE_rough$test.y, 
                        pred.condition=names(RESULT_DFs$SLE_rough$pred),
                        BA=RESULT_DFs$SLE_rough$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_rough$condition==RESULT_rough$pred.condition)
RESULT_rough=RESULT_rough %>% select(-pred.condition) %>% 
  mutate(celltypes=gsub("\\*\\*\\*.*","",condition),
         condition=gsub(".*\\*\\*\\*","",condition)) %>%
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(gsub("SLE ","",condition), ".", CA_or_BA))
RESULT_rough$condition_train=forcats::fct_relevel(RESULT_rough$condition_train, c("INACT.CA","INACT.BA","ACT.CA","ACT.BA"))
RESULT_rough$condition=forcats::fct_relevel(RESULT_rough$condition, c("SLE INACT","SLE ACT"))

plot_rough=
  ggplot(RESULT_rough, aes(x=condition_train, y=age, color=condition)) +
  geom_point(size=1) +
  facet_wrap(~celltypes) +
  scale_color_manual(values=c("lightblue3","coral3")) +
  geom_line(aes(group=pair), linewidth=0.25) +
  ggpubr::stat_compare_means(comparisons=list(c("INACT.CA","INACT.BA"), c("ACT.CA","ACT.BA")), label.y=65,
                             paired=T, size=3.5, bracket.size=0, label=after_stat("p.format")) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none",
        strip.background=element_rect(fill="transparent"),
        strip.text.x=element_text(size=10, color="black")) +
  scale_x_discrete(guide=guide_axis(n.dodge=2))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inSLE_rough.pdf", height=5, width=8)
plot(plot_rough)
dev.off()

### At the Annot.inter level
RESULT_inter=data.frame(condition=names(RESULT_DFs$SLE_inter$test.y),
                        CA=RESULT_DFs$SLE_inter$test.y, 
                        pred.condition=names(RESULT_DFs$SLE_inter$pred),
                        BA=RESULT_DFs$SLE_inter$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_inter$condition==RESULT_inter$pred.condition)
RESULT_inter=RESULT_inter %>% select(-pred.condition) %>% 
  mutate(celltypes=gsub("\\*\\*\\*.*","",condition),
         condition=gsub(".*\\*\\*\\*","",condition)) %>%
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(gsub("SLE ","",condition), ".", CA_or_BA))
RESULT_inter$condition_train=forcats::fct_relevel(RESULT_inter$condition_train, c("INACT.CA","INACT.BA","ACT.CA","ACT.BA"))
RESULT_inter$condition=forcats::fct_relevel(RESULT_inter$condition, c("SLE INACT","SLE ACT"))

plot_inter=
  ggplot(RESULT_inter, aes(x=condition_train, y=age, color=condition)) +
    geom_point(size=1) +
    facet_wrap(~celltypes, ncol=3) +
    scale_color_manual(values=c("lightblue3","coral3")) +
    geom_line(aes(group=pair), linewidth=0.25) +
    ggpubr::stat_compare_means(comparisons=list(c("INACT.CA","INACT.BA"), c("ACT.CA","ACT.BA")), label.y=65,
                               paired=T, size=3.5, bracket.size=0, label=after_stat("p.signif")) +
    theme_light() +
    labs(x=NULL, y="age") +
    # scale_x_discrete(guide=guide_axis(n.dodge=2)) +
    theme(axis.text.x=element_text(size=9, angle=90, vjust=0.5, hjust=1),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          panel.grid.major=element_blank(),
          title=element_text(size=10),
          legend.position="none",
          strip.background=element_rect(fill="transparent"),
          strip.text.x=element_text(size=10, color="black"))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inSLE_inter.pdf", height=7, width=4)
plot(plot_inter)
dev.off()

#####################################



### Analyze the CD4T, CD8T, and NK proportions
#####################################
### 

library(Seurat)
library(UCell)
library(dplyr)

### Process the SLE data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
test_sle[[]]=test_sle[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))

# ### Check if integration is needed before clustering
# # turns out to be not needed
# test_sle_try=test_sle %>%
#   NormalizeData() %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA() %>%
#   FindNeighbors(dims=1:30) %>%
#   RunUMAP(dims=1:30) %>%
#   FindClusters(resolution=0.5)
# DimPlot(test_sle_try, group.by=c("subject_id", "fine_cell_type"), reduction="umap")

# ### Cluster and annotate CD4T.naive, CD4T.mem, CD8T.naive, CD8T.mem, and NK.CD56dim
# turns out to be bad clustering due to too low seq depth

### Determine CD4T.naive, CD4T.mem, CD8T.naive, CD8T.mem, and NK.CD56dim by Ucell scoring
# CD4T
CD4T_subset=subset(test_sle, fine_cell_type %in% c("CD4+ T cells"))
CD4T_scoring=list(CD4T.naive=c("TCF7", "CD4", "CCR7", "IL7R", "FHIT","LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1"),
                  CD4T.mem=c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL", "CCL5", "FYB1", "GZMK", "GZMA", "KLRB1"))
CD4T_Ucell=AddModuleScore_UCell(CD4T_subset, features=CD4T_scoring)
score_df=CD4T_Ucell[[]][, colnames(CD4T_Ucell[[]])[grepl("_UCell$", colnames(CD4T_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df=list(CD4T=score_df)
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")
# CD8T
CD8T_subset=subset(test_sle, fine_cell_type %in% c("CD8+ T cells"))
CD8T_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP","LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"),
                  CD8T.mem=c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB","CCR4","CCR7","CD27","CD3D","CD3G","CD3E","SELL","CD69",
                             "KLRC1","KLRC2","KLRC3","IKZF2","CCL5","GZMH","KLRD1","TRGC2","NKG7","GZMK","CST7","GZMB","FCGR3A","HNRNPH1","XCL1","GATA3",
                             "LYAR","MALAT1","TSHZ2"))
CD8T_Ucell=AddModuleScore_UCell(CD8T_subset, features=CD8T_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df_orig=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")
score_df=c(score_df_orig, list(CD8T=score_df))
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")
# NK
NK_subset=subset(test_sle, fine_cell_type %in% c("NK cells"))
NK_scoring=list(NK.CD56dim=c("GNLY", "TYROBP", "NKG7", "FCER1G", "GZMB", "TRDC", "PRF1", "FGFBP2", "SPON2", "KLRF1"),
                  NK.CD56hi=c("XCL2", "FCER1G", "SPINK2", "TRDC", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A"))
NK_Ucell=AddModuleScore_UCell(NK_subset, features=NK_scoring)
score_df=NK_Ucell[[]][, colnames(NK_Ucell[[]])[grepl("_UCell$", colnames(NK_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df_orig=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")
score_df=c(score_df_orig, list(NK=score_df))
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")

### Arrange the scores for CD4T, CD8T, and NK
score_df=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.SLE_CD4TCD8TNK_scores.rds")
df_arranged=test_sle[[]] %>% 
  subset(fine_cell_type %in% c("CD4+ T cells","CD8+ T cells","NK cells")) %>%
  dplyr::select(subject_id, classification, age, fine_cell_type, age_cut) %>%
  tibble::rownames_to_column("cellid") %>%
  left_join(., score_df[[1]]) %>%
  left_join(., score_df[[2]]) %>%
  left_join(., score_df[[3]])
df_arranged=df_arranged %>% subset(age_cut!="69~97") # 69~97 has only SLE INACT group but no SLE ACT or Control, so remove this
  
### Check the scores across age_cut and celltypes
# turns out that only CD8T.mem scores (in 48~68) showed a diff across disease status
df_selected=subset(df_arranged, fine_cell_type %in% c("CD4+ T cells", "CD8+ T cells", "NK cells")) %>% 
  # dplyr::select(-c(CD4T.naive_UCell,CD4T.mem_UCell,NK.CD56dim_UCell,NK.CD56hi_UCell)) %>%
  tidyr::pivot_longer(cols=colnames(.)[grepl("UCell",colnames(.))], values_to="Score", names_to="group") %>%
  mutate(group=gsub("_UCell","",group)) %>%
  mutate(classification=ifelse(classification=="Control","CTRL",ifelse(classification=="SLE ACT","ACT","INACT")))
df_selected$group=forcats::fct_relevel(df_selected$group, c("CD4T.naive","CD4T.mem","CD8T.naive","CD8T.mem","NK.CD56dim","NK.CD56hi"))
df_selected$classification=forcats::fct_relevel(df_selected$classification, c("CTRL","INACT","ACT"))

plot_all=
  ggplot(df_selected, aes(x=Score, color=classification)) +
  facet_grid(group~age_cut) +
  geom_density(linewidth=0.25) +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  scale_x_continuous(breaks=c(0.4,0.8), labels=c(0.4,0.8)) +
  theme_classic() +
  labs(x="UCell score", y="density", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(linewidth=0.5)))

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_CellProportions_CD4T.CD8T.NK.pdf", height=6.5, width=3.5)
plot_all
dev.off()

### Plot CD8T.mem scores across disease status
df_CD8T.mem=subset(df_selected, fine_cell_type=="CD8+ T cells" & group=="CD8T.mem")
plot_=
  ggplot(df_CD8T.mem, aes(x=Score, color=classification)) +
  facet_wrap(~age_cut) +
  geom_density(linewidth=0.25) +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  theme_classic() +
  labs(x="UCell score", y="density", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(linewidth=0.5)))

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_CellProportions_CD8Tmem.pdf", height=1.5, width=5)
plot_
dev.off()

#####################################



### Analyze the expression of metabolic- and immune-related hub genes in the CD4T & CD8T from the SLE dataset
#####################################
### 

### Load the metabolic- and immune-related hub genes from the PBMCage WGCNA modules
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe()
  )
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="regulation of T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)
names(Module_genes)=ordered_terms
Merged_modules=
  list(
    metabolism=c("cytoplasmic translation","tRNA modification","chromosome localization/mitotic checkpoint",
                 "oxidative phosphorylation","regulation of hemopoiesis"),
    immunity=c("defense response","humoral immunity","regulation of T cell differentiation",
               "leukocyte-mediated cytotoxicity","antigen processing and presentation",
               "hemostasis")
  )
Module_genes=lapply(Merged_modules, function(x) Module_genes[x] %>% Reduce(union, .))

### Load the expr data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
test_sle[[]]=test_sle[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
test_sle_subset=subset(test_sle, fine_cell_type %in% c("CD4+ T cells","CD8+ cytotoxic T cells","CD8+ T cells"))
test_sle_subset=subset(test_sle_subset, age_cut!="69~97")
DefaultAssay(test_sle_subset)="RNA"

### Arrange the metabolic- and immunity-related df
DF_all=list()
for (i in 1:length(Merged_modules)) {
  genes.selected=Module_genes[[i]]
  expr_data_meta=FetchData(test_sle_subset, vars=c("age", "classification", "age_cut", "subject_id", "fine_cell_type"), layer="data")
  expr_data_all=FetchData(test_sle_subset, vars=c(genes.selected), layer="data")
  unique(rownames(expr_data_meta)==rownames(expr_data_all)) # check
  expr_data_all_filter=names(expr_data_all)[colSums(expr_data_all)>0]
  expr_data_all_arranged=expr_data_all %>% dplyr::select(all_of(expr_data_all_filter)) %>% rowMeans(.)
  expr_data_all_arranged=data.frame(cellid=names(expr_data_all_arranged), mean.expr=expr_data_all_arranged) %>%
    right_join(., expr_data_meta %>% tibble::rownames_to_column("cellid"), by="cellid") %>%
    mutate(classification=ifelse(classification=="Control","CTRL",ifelse(classification=="SLE ACT","ACT","INACT")))
  expr_data_all_arranged$classification=forcats::fct_relevel(expr_data_all_arranged$classification, c("CTRL","INACT","ACT"))
  DF_all[[i]]=expr_data_all_arranged
}
names(DF_all)=names(Merged_modules)

### Plot
# plot metabolism
plot_1=
  ggplot(DF_all[[1]], aes(x=classification, y=mean.expr, color=classification)) +
  facet_wrap(~age_cut) +
  geom_boxplot(outlier.size=0.5, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(data=. %>% subset(age_cut!="19~30"), aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","INACT"),c("INACT","ACT"),c("CTRL","ACT"))) +
  ggpubr::stat_compare_means(data=. %>% subset(age_cut=="19~30"), aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","ACT"))) +
  theme_classic() +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  labs(x=NULL, y="Avg. expression", subtitle="metabolism") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
  scale_y_continuous(expand=c(0, 0.005))
# plot immunity
plot_2=
  ggplot(DF_all[[2]], aes(x=classification, y=mean.expr, color=classification)) +
  facet_wrap(~age_cut) +
  geom_boxplot(outlier.size=0.5, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(data=. %>% subset(age_cut!="19~30"), aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","INACT"),c("INACT","ACT"),c("CTRL","ACT"))) +
  ggpubr::stat_compare_means(data=. %>% subset(age_cut=="19~30"), aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","ACT"))) +
  theme_classic() +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  labs(x=NULL, y="Avg. expression", subtitle="immunity") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
    scale_y_continuous(expand=c(0, 0.1))
  
plot_both=cowplot::plot_grid(plotlist=list(plot_1, plot_2), ncol=1)

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_Metabolism.Immune_hub.genes.expr.pdf", height=5.5, width=5)
plot_both
dev.off()

### If not split age_cut when comparing across disease status
plot_notsplit_1=
  ggplot(DF_all[[1]], aes(x=classification, y=mean.expr, color=classification)) +
  # facet_wrap(~fine_cell_type) +
  geom_boxplot(outlier.size=1, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","INACT"),c("INACT","ACT"),c("CTRL","ACT"))) +
  theme_classic() +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  labs(x=NULL, y="Avg. expression", subtitle="metabolism") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
  scale_y_continuous(expand=c(0,0.005))
  
plot_notsplit_2=
  ggplot(DF_all[[2]], aes(x=classification, y=mean.expr, color=classification)) +
  # facet_wrap(~fine_cell_type) +
  geom_boxplot(outlier.size=1, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","INACT"),c("INACT","ACT"),c("CTRL","ACT"))) +
  theme_classic() +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  labs(x=NULL, y="Avg. expression", subtitle="immunity") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
    scale_y_continuous(expand=c(0,0.1))

plot_both_notsplit=cowplot::plot_grid(plotlist=list(plot_notsplit_1, plot_notsplit_2), ncol=1)

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_Metabolism.Immune_hub.genes.expr_notsplitagecut.pdf", height=5.5, width=2)
plot_both_notsplit
dev.off()

#####################################



### Analyze the expression of the key kinases (JAK1/3, MAPK1/8/MAP3K7, LCK, AKT1) and TFs (NFKB familiy, MYC)
#####################################
###

library(Seurat)

## Process the SLE data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
unique(test_sle$fine_cell_type)
test_sle_subset=subset(test_sle, fine_cell_type %in% c("CD4+ T cells","CD8+ cytotoxic T cells","CD8+ T cells"))
test_sle_subset$donor_id=test_sle_subset$subject_id
test_sle_subset_agg=AggregateExpression(test_sle_subset, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id","classification"))
test_sle_subset_agg$age=as.numeric(unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[1])))
test_sle_subset_agg$donor_id=unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[2]))
test_sle_subset_agg$condition=unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[3]))

### Extract the gene expression
genes_of_interest=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","AKT1","LCK","MYC","NFKB1","NFKB2")
expr_data=FetchData(test_sle_subset_agg, vars=c(genes_of_interest, colnames(test_sle_subset_agg[[]])), layer="data")
expr_data_arranged=expr_data %>%
  dplyr::select(-c(orig.ident)) %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("age|donor_id|condition",colnames(.))], names_to="gene", values_to="data")
# sort the gene expr
genes_of_interest_list=list(MYC=c("MYC"), 
                            `Metabolism-related kinases`=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","AKT1"),
                            NFKB2=c("NFKB2"),
                            `Immunity-related kinases`=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","LCK"))
mean_expr_data=
  lapply(1:length(genes_of_interest_list), function(idx) expr_data_arranged %>% 
           subset(gene %in% genes_of_interest_list[[idx]]) %>%
           group_by(donor_id, age, condition) %>%
           summarize_at("data", mean) %>%
           mutate(group=names(genes_of_interest_list)[idx])) %>%
  data.table::rbindlist(.) %>%
  mutate(condition=ifelse(condition=="Control","CTRL",ifelse(condition=="SLE ACT","ACT","INACT"))) %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
mean_expr_data$group=forcats::fct_relevel(mean_expr_data$group, c("MYC","NFKB2","Metabolism-related kinases","Immunity-related kinases"))
mean_expr_data$condition=forcats::fct_relevel(mean_expr_data$condition, c("CTRL","INACT","ACT"))

### Plot
plot_TFKinase=
  ggplot(mean_expr_data, aes(x=age, y=data, color=condition)) +
  facet_wrap(~group) +
  geom_point(size=0.5, shape=1) +
  geom_smooth(aes(group=condition), method="lm", se=F, linewidth=0.25, show.legend=F) +
  theme_classic() +
  scale_color_manual(values=c("grey50","lightblue3","coral3")) +
  labs(x="age", y="Avg. expression", subtitle=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical",
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.subtitle=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, shape=19)))

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_Metabolism.Immune_TFsAndKinases.expr.pdf", height=3.5, width=5)
plot_TFKinase
dev.off()

#####################################



### Analyze the communication among cell types
#####################################
###

library(dplyr)
library(ggplot2)
library(Seurat)

### Turns out that CellChat doesn't work as the netP is all ==1

### Load the expr data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
test_sle[[]]=test_sle[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
test_sle_subset=subset(test_sle, age_cut!="69~97")

### Extract the key interactions
unique(test_sle_subset$fine_cell_type)
B_CD8T=list(`Naive B cells,Memory B cells,ABCs`=c("HLA-A","HLA-B","HLA-C","HLA-E"), 
            `CD8+ T cells`=c("CD8B"))
allother_NK=list(`Naive B cells,Memory B cells,ABCs,CD4+ T cells,CD14+ monocytes,CD16+ monocytes,CD8+ T cells,NK cells`=c("HLA-E"), 
                 `NK cells`=c("KLRK1","KLRD1","KLRC1"))
CD8T_NK=list(`CD8+ T cells`=c("CD99"), `NK cells`=c("CD99"))
NK_CD8T=list(`NK cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"), `CD8+ T cells`=c("CD8A"))
CD4T_CD8T=list(`CD4+ T cells`=c("MIF"), `CD8+ T cells`=c("CD44","CD74","CXCR4"))

all_to_CD8T=list(`Naive B cells,Memory B cells,ABCs,NK cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"), `CD8+ T cells`=c("CD8B"))
all_to_NK=list(`Naive B cells,Memory B cells,ABCs,CD4+ T cells,CD14+ monocytes,CD16+ monocytes,CD8+ T cells,NK cells`=c("HLA-E"), 
                 `NK cells`=c("KLRK1","KLRD1","KLRC1"))

### Plot the signals to CD8T and NK, respectively
molecule_list=list(all_to_CD8T, all_to_NK)
plot_lists=list()
for (j in 1:length(molecule_list)) {
  DF_all=lapply(1:length(molecule_list[[j]]), function(idx) FetchData(test_sle_subset %>% 
                                                                   subset(fine_cell_type %in% strsplit(names(molecule_list[[j]][idx]), split=",")[[1]]), 
                                                                 vars=c(molecule_list[[j]][[idx]], "age", "age_cut", "classification", "subject_id"), 
                                                                 layer="data"))
  tempt=lapply(1:length(DF_all), function(idx) DF_all[[idx]] %>% dplyr::select(any_of(molecule_list[[j]][[idx]])) %>% rowSums())
  tempt=lapply(1:length(tempt), function(idx) data.frame(cellid=names(tempt[[idx]]), expr=tempt[[idx]]))
  DF_all_mean=lapply(1:length(DF_all), function(idx) {
    tempt_=DF_all[[idx]] %>% 
      dplyr::select(all_of(c("age", "age_cut", "classification", "subject_id"))) %>%
      tibble::rownames_to_column("cellid") %>%
      left_join(., tempt[[idx]]) %>%
      group_by(subject_id, age, age_cut, classification) %>%
      summarize_at("expr", mean)
    colnames(tempt_)[ncol(tempt_)]=names(molecule_list[[j]][idx])
    tempt_
  }) %>%
    Reduce(full_join, .) %>%
    mutate(classification=ifelse(classification=="Control","CTRL",ifelse(classification=="SLE ACT","ACT","INACT")))
  colnames(DF_all_mean)[c(5,6)]=c("interaction_A","interaction_B")
  DF_all_mean$classification=forcats::fct_relevel(DF_all_mean$classification, c("CTRL","INACT"))
  
  plot_lists[[j]]=
    ggplot(DF_all_mean, aes(x=interaction_A, y=interaction_B, color=classification)) +
    geom_point(size=0.5, shape=1) +
    geom_smooth(aes(group=classification), show.legend=F, se=F, method="lm", linewidth=0.5) +
    ggpmisc::stat_correlation(mapping=ggpmisc::use_label(c("R", "P"), aes(group=classification)), 
                              small.p=T,
                              show.legend=F, size=3, vstep=0.05) +
    # scale_color_discrete(labels=c("CTRL","INACT","ACT")) +
    theme_classic() +
    scale_color_manual(values=c("grey50","lightblue3","coral3")) +
    labs(x="Avg. expr. of MHC-I", y=c("Avg. expr. of CD8B","Avg. expr. of NK receptor")[j], 
         subtitle=c("Signals to CD8T cells","Signals to NK cells")[j]) +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position=c("none","right")[j],
          legend.title=element_blank(),
          legend.text=element_text(size=9),
          plot.subtitle=element_text(size=11)) +
    guides(color=guide_legend(override.aes=list(size=2, shape=19)))
}

plot_both=cowplot::plot_grid(plotlist=plot_lists, align="h", axis="tb", rel_widths=c(0.65,1))

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_CellchatInteractions.expr.pdf", height=5.5, width=5)
plot_both
dev.off()

#####################################



### Analyze the expression patterns
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

### Load the expr data
test_sle=readRDS("~/Project_PBMCage/Lupus_RA/SLEdataset_seuobj.RDS")
DefaultAssay(test_sle)="RNA"
test_sle_aggr=AggregateExpression(test_sle, return.seurat=T, group.by=c("subject_id","fine_cell_type"))
updated_meta=test_sle_aggr[[]] %>% 
  mutate(subject_id=orig.ident, fine_cell_type=gsub(".*_","",rownames(.))) %>%
  tibble::rownames_to_column("id") %>% 
  left_join(., test_sle[[]] %>% dplyr::select(age, subject_id, classification, fine_cell_type) %>% subset(!duplicated(.))) %>%
  tibble::column_to_rownames("id")
test_sle_aggr[[]]=updated_meta
test_sle_aggr[[]]=test_sle_aggr[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
unique(test_sle$fine_cell_type)
test_sle_subset=subset(test_sle_aggr, age_cut!="69~97")

### Load the 5 expression patterns from the PBMCage project
pattern_trend_genes_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.txt.gz", sep="\t")
celltype_map=unique(pattern_trend_genes_df$celltype)
pattern_df_modified=pattern_trend_genes_df %>%
  mutate(celltype=ifelse(celltype=="B","Naive B cells,Transitional B cells,Memory B cells,Plasmablasts,ABCs",
                         ifelse(celltype=="CD4T","CD4+ T cells",
                                ifelse(celltype=="CD8T","CD8+ cytotoxic T cells,CD8+ T cells",
                                       ifelse(celltype=="Mono","CD14+ monocytes,CD16+ monocytes",
                                              ifelse(celltype=="NK","NK cells",
                                                     ifelse(celltype=="DC","cDCs,pDCs",celltype))))))) %>%
  subset(celltype!="OtherT" & celltype!="Other") %>%
  split(.$pattern)

### Extract the mean.expr for each pattern
expr_all_patterns=list()
for (i in 1:length(pattern_df_modified)) {
  celltype_list=unique(pattern_df_modified[[i]]$celltype)
  expr_for_celltype_list=list()
  for (j in 1:length(celltype_list)) {
    gene_per_cell_per_pattern=pattern_df_modified[[i]] %>% 
      subset(celltype==celltype_list[j]) %>% dplyr::select(gene) %>% tibble::deframe()
    cellname=strsplit(celltype_list[j], split=",")[[1]]
    expr_data=
      FetchData(test_sle_subset %>% subset(fine_cell_type %in% cellname), vars=c(gene_per_cell_per_pattern, "age","age_cut","classification","subject_id"),
                layer="data")
    # remove the genes with all zeros
    expr_data_filter=colSums(expr_data %>% .[,colnames(.)[!grepl("age|classification|subject",colnames(.))]]) %>% .[.>0] %>% names(.)
    expr_data_cleaned=cbind(expr_data[, c(expr_data_filter)], 
                            expr_data[, c("age","age_cut","classification","subject_id")])
    expr_data_cleaned=expr_data_cleaned %>% group_by(age, age_cut, classification, subject_id) %>%
      summarize_at(colnames(.)[!grepl("age|classification|subject",colnames(.))], mean) %>%
      tidyr::pivot_longer(cols=colnames(.)[!grepl("age|classification|subject",colnames(.))],
                          values_to="mean.expr", names_to="gene") %>%
      mutate(celltype=celltype_list[j])
    expr_for_celltype_list[[j]]=expr_data_cleaned
  }
  all_expr_for_pattern=data.table::rbindlist(expr_for_celltype_list) %>%
    mutate(classification=ifelse(classification=="Control","CTRL", ifelse(classification=="SLE INACT","INACT","ACT")))
  all_expr_for_pattern$classification=forcats::fct_relevel(all_expr_for_pattern$classification, c("CTRL","INACT","ACT"))
  
  expr_all_patterns[[i]]=all_expr_for_pattern
}
saveRDS(expr_all_patterns, "~/Project_PBMCage/Map_to_SLE/ExpressionPatterns_mean.expr.per.genes.rds")

### Plot
expr_all_patterns=readRDS("~/Project_PBMCage/Map_to_SLE/ExpressionPatterns_mean.expr.per.genes.rds")
plot_lists=list()
for (i in 1:length(expr_all_patterns)) {
  quantile_lim=expr_all_patterns[[i]] %>% group_by(age_cut, classification) %>% 
    mutate(median=quantile(mean.expr, 0.5),
           upper_hinge=quantile(mean.expr, 0.75),
           lower_hinge=quantile(mean.expr, 0.25)) %>%
    mutate(IQR=upper_hinge-lower_hinge) %>%
    mutate(upper_whisker=upper_hinge+1.5*IQR,
           lower_whisker=lower_hinge-1.5*IQR) %>%
    dplyr::select(age_cut, classification, upper_hinge, lower_hinge, IQR, upper_whisker, lower_whisker) %>%
    subset(!duplicated(.))
  
  plot_lists[[i]]=
    ggplot(expr_all_patterns[[i]], aes(x=age_cut, y=mean.expr, color=classification)) +
    facet_wrap(~classification) +
    geom_boxplot(outlier.shape=NA, outlier.size=0.5, width=0.5, fill="transparent", coef=0, linewidth=0.25) +
    scale_color_manual(values=c("grey50","lightblue3","coral3")) +
    theme_classic() +
    labs(x=NULL, y="Avg. expression", subtitle=paste0("Pattern",i)) +
    theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position="none",
          legend.title=element_blank(),
          plot.subtitle=element_text(size=11)) +
    coord_cartesian(ylim=c(min(quantile_lim$lower_hinge), max(quantile_lim$upper_hinge)))
}
plot_patterns_all=cowplot::plot_grid(plotlist=plot_lists, axis="tb", align="h", nrow=1)

pdf("~/Project_PBMCage/Plots/Map_to_SLE/MapToSLE_ExprPatterns.pdf", height=1.75, width=18)
plot_patterns_all
dev.off()

#####################################



### Analyze the timing patterns
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

#####################################



##########################################################################

### Apply AgeGene5K-based model on Covid_flu dataset

##########################################################################


### Reconstruct the Covid-Flu dataset
#####################################
###

###
Covid_flu=readRDS("/home/jo_79/Project_PBMCage/GSE158055/GSE206265_covid_flu.CITEseq.pseudobulk.gene.expr.RDS")
Covid_flu_celltypes=list()
for (i in 1:length(Covid_flu)) {
  meta_data=Covid_flu[[i]]$samples
  count_data=Covid_flu[[i]]$counts
  Covid_flu_celltypes[[i]]=Seurat::CreateSeuratObject(counts=count_data, meta.data=meta_data, min.cells=0, min.features=0)
  Covid_flu_celltypes[[i]]$celltype=names(Covid_flu)[i]
}
Covid_flu_celltypes_total=merge(Covid_flu_celltypes[[1]], Covid_flu_celltypes[-1], add.cell.ids=names(Covid_flu))
Covid_flu_celltypes_merged=JoinLayers(Covid_flu_celltypes_total)

# select rough level and save pseudobulk obj
lapply(Covid_flu_celltypes, dim)
names(Covid_flu)
Covid_flu_celltypes_rough=subset(Covid_flu_celltypes_merged, celltype %in% c("B","CD4","CD8","cDC","gdT-Vd2","HSPC","ILC","Mono","MAIT","Neut","NK","pDC","Plasmablast","Platelet"))
rough_df=Covid_flu_celltypes_rough[[]]
rough_df=rough_df %>% 
  mutate(condition=ifelse(is.na(long.covid.symptoms), group, paste0(group, ".", long.covid.symptoms))) %>%
  mutate(Annot.rough=ifelse(celltype=="B" | celltype=="Plasmablast","B cells",celltype)) %>%
  mutate(Annot.rough=ifelse(celltype=="CD4","CD4T cells",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype=="CD8","CD8T cells",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype=="cDC" | celltype=="pDC","DCs",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype=="gdT-Vd2" | celltype=="MAIT","OtherT",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype %in% c("HSPC","ILC","Neut","Platelet"),"Other cells",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype=="Mono","Monocytes",Annot.rough)) %>%
  mutate(Annot.rough=ifelse(celltype=="NK","NK cells",Annot.rough))
table(rough_df$Annot.rough, useNA="ifany")
table(rough_df$condition, useNA="ifany")

Covid_flu_celltypes_rough$Annot.rough=rough_df$Annot.rough
Covid_flu_celltypes_rough$condition=rough_df$condition
table(Covid_flu_celltypes_rough$Annot.rough, useNA="ifany")
table(Covid_flu_celltypes_rough$condition, useNA="ifany")

saveRDS(Covid_flu_celltypes_rough, "~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_rough.rds")

# select rough level and save pseudobulk obj
Covid_flu_celltypes_inter=subset(Covid_flu_celltypes_merged, celltype %in% c("B_Mem","B_Naive","B_Naive_Intermediate",
                                                                             "CD4_CM","CD4_EM","CD4_Naive","CD4_Tfh","CD4_Treg",
                                                                             "CD8_CM","CD8_EM","CD8_Naive","CD8_proliferating","CD8_TEMRA","CD8_TRM",
                                                                             "cDC","gdT-Vd2","HSPC","ILC","Mono_Classical","Mono_Intermediate","Mono_NonClassical",
                                                                             "MAIT","Neut","NK_CD16hi","NK_CD56hiCD16lo","NK_proliferating","pDC","Plasmablast","Platelet"))

inter_df=Covid_flu_celltypes_inter[[]]
inter_df=inter_df %>% 
  mutate(condition=ifelse(is.na(long.covid.symptoms), group, paste0(group, ".", long.covid.symptoms))) %>%
  mutate(Annot.inter=ifelse(celltype=="B_Mem","B.mem",celltype)) %>%
  mutate(Annot.inter=ifelse(celltype=="B_Naive","B.naive",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="B_Naive_Intermediate","B.inter",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD4_CM","CD4T.Tcm",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD4_EM","CD4T.Tem",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD4_Naive","CD4T.naive",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD4_Tfh","CD4T.Tfh",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD4_Treg","CD4T.Treg",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_CM","CD8T.Tcm",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_EM","CD8T.Tem",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_Naive","CD8T.naive",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_proliferating","CD8T.prolif",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_TEMRA","CD8T.Temra",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="CD8_TRM","CD8T.Trm",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="cDC","DC.cDC",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="gdT-Vd2","OtherT.gdT",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="HSPC","Other.HSPC",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="ILC","Other.ILC",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Mono_Classical","Mono.classical",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Mono_Intermediate","Mono.inter",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Mono_NonClassical","Mono.nonclassical",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="MAIT","OtherT.MAIT",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Neut","Other.Neutro",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="NK_CD16hi","NK.CD56dim",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="NK_CD56hiCD16lo","NK.CD56hi",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="NK_proliferating","NK.prolif",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="pDC","DC.pDC",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Plasmablast","B.pbpc",Annot.inter)) %>%
  mutate(Annot.inter=ifelse(celltype=="Platelet","Other.plt",Annot.inter))
table(inter_df$Annot.inter, useNA="ifany")
table(inter_df$condition, useNA="ifany")

inter_df=inter_df %>%
  mutate(Annot.rough=gsub("\\..*","",Annot.inter)) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="B","B cells",
                            ifelse(Annot.rough=="CD4T","CD4T cells",
                                   ifelse(Annot.rough=="CD8T","CD8T cells",
                                          ifelse(Annot.rough=="DC","DCs",
                                                 ifelse(Annot.rough=="NK","NK cells",
                                                        ifelse(Annot.rough=="Other","Other cells","OtherT")))))))
table(inter_df$Annot.rough, useNA="ifany")

Covid_flu_celltypes_inter$Annot.inter=inter_df$Annot.inter
Covid_flu_celltypes_inter$Annot.rough=inter_df$Annot.rough
Covid_flu_celltypes_inter$condition=inter_df$condition
table(Covid_flu_celltypes_inter$Annot.inter, useNA='ifany')
table(Covid_flu_celltypes_inter$Annot.rough, useNA='ifany')
table(Covid_flu_celltypes_inter$condition, useNA='ifany')

saveRDS(Covid_flu_celltypes_inter, "~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")

# save the pseudobulk obj at the donorid level
pseudobulkobj=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_rough.rds")
pseudobulkobj$subset.id_visit.id=gsub(".*_","",rownames(pseudobulkobj[[]]))
pseudobulkobj_donor=AggregateExpression(pseudobulkobj, assays="RNA", return.seurat=T, group.by=c("age","sex","subset.id_visit.id","condition"))
pseudobulkobj_donor$age=as.numeric(unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[1])))
pseudobulkobj_donor$sex=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[2]))
pseudobulkobj_donor$subset.id_visit.id=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[3]))
pseudobulkobj_donor$condition=unlist(lapply(strsplit(colnames(pseudobulkobj_donor), split="_"), function(x) x[4]))

saveRDS(pseudobulkobj_donor, "~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_donorid.rds")

#####################################



#####################################

RESULT_DFs=list()
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")

### RNA-based training with the CovidFlue Ctrl at the donor level
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

## Process the SLE data
test_sle=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_donorid.rds")

# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("age","sex","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("HC$",rownames(test_sle_expr_total)),]
test_sle_expr=test_sle_expr_total[grepl("COVR\\.",rownames(test_sle_expr_total)),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, sex, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, sex, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=training_PBMCage_expr$condition
test.y_withname=test.y; names(test.y_withname)=test_sle_expr$condition
pred_withname=pred; names(pred_withname)=test_sle_expr$condition
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs$`RNA-based`[["CovidFlu_donorid"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")

#####################################



### RNA-based training with the CovidFlu Ctrl at the rough level
#####################################
###

### Load the pseudobulk_rough data
test_obj_pseudo_rough=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_rough.rds")
test_obj_pseudo_rough=subset(test_obj_pseudo_rough, condition!="NA")
test_obj_pseudo_rough=test_obj_pseudo_rough %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

## Process the covid data
test_sle=test_obj_pseudo_rough
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("Annot.rough","age","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))
rownames(test_sle_expr_total)=paste0(test_sle_expr_total$Annot.rough, ">>", rownames(test_sle_expr_total))
test_sle_expr_total=test_sle_expr_total %>% select(-Annot.rough)

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("HC",test_sle_expr_total$condition),]
test_sle_expr=test_sle_expr_total[grepl("COVR",test_sle_expr_total$condition),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=paste0(gsub(">>.*","",rownames(training_rpca.df)),"***",training_PBMCage_expr$condition)
test.y_withname=test.y; names(test.y_withname)=paste0(gsub(">>.*","",rownames(test_sle_expr)),"***",test_sle_expr$condition)
pred_withname=pred; names(pred_withname)=paste0(gsub(">>.*","",rownames(test_sle_expr)),"***",test_sle_expr$condition)
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs$`RNA-based`[["CovidFlu_rough"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")

#####################################



### RNA-based training with the CovidFlu Ctrl at the rough level and then merged on donor
#####################################
###
library(UCell)
library(Seurat)
library(pls)

### Load the pseudobulk_rough data
test_obj_pseudo_rough=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_rough.rds")
test_obj_pseudo_rough=subset(test_obj_pseudo_rough, condition!="NA")
test_obj_pseudo_rough=test_obj_pseudo_rough %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

## Process the covid data
test_sle=test_obj_pseudo_rough
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("Annot.rough","age","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))
rownames(test_sle_expr_total)=paste0(test_sle_expr_total$Annot.rough, ">>", rownames(test_sle_expr_total))
test_sle_expr_total=test_sle_expr_total %>% select(-Annot.rough)

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("HC",test_sle_expr_total$condition),]
test_sle_expr=test_sle_expr_total[grepl("COVR",test_sle_expr_total$condition),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=rownames(training_rpca.df)
test.y_withname=test.y; names(test.y_withname)=rownames(test_sle_expr)
pred_withname=pred; names(pred_withname)=rownames(test_sle_expr)
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs$`RNA-based`[["CovidFlu_rough_onDonor"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")

#####################################



### RNA-based training with the CovidFlu Ctrl at the inter level
#####################################
###

### Load the pseudobulk_rough data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter=test_obj_pseudo_inter %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()
### Load the Xu_5000 geneset
gene_list=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Compare.to.Published_Genesets_GenesetsToCompare.rds")
Xu_5000.genes=gene_list$Xu_5000.ranked
topN_genes=Xu_5000.genes[1:300] # take topN genes
# in case gene symbols here do not match those in the training/testing dataset
topN_genes=c(topN_genes, gsub("-","_",topN_genes))

## Process the covid data
test_sle=test_obj_pseudo_inter
# add CD8T scores as another feature
T_cell_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP", "LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"))
CD8T_Ucell=AddModuleScore_UCell(test_sle, features=T_cell_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]]
test_sle$CD8T.naive=score_df
# extract the expression of the genes used in the training data
test_sle_expr_total=FetchData(test_sle, vars=c("Annot.inter","age","CD8T.naive","condition",topN_genes), layer="counts")
colnames(test_sle_expr_total)=gsub("-","_",colnames(test_sle_expr_total))
rownames(test_sle_expr_total)=paste0(test_sle_expr_total$Annot.inter, ">>", rownames(test_sle_expr_total))
test_sle_expr_total=test_sle_expr_total %>% select(-Annot.inter)

### Prepare Ctrl as the training set and patients as the applied set
# scale the gene expression but leave age, sex, CD8T.naive as origin
training_PBMCage_expr=test_sle_expr_total[grepl("HC",test_sle_expr_total$condition),]
test_sle_expr=test_sle_expr_total[grepl("COVR",test_sle_expr_total$condition),]

training_rpca.df=t(scale(t(as.matrix(training_PBMCage_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
training_rpca.df$age=training_PBMCage_expr$age
training_rpca.df$CD8T.naive=training_PBMCage_expr$CD8T.naive

testing_rpca.df=t(scale(t(as.matrix(test_sle_expr %>% select(-c(age, CD8T.naive, condition)))))) %>% as.data.frame()
testing_rpca.df$age=test_sle_expr$age
testing_rpca.df$CD8T.naive=test_sle_expr$CD8T.naive

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
if (best.ncomp==0) best.ncomp=1

### Predict
pred=predict(Model, test.x, ncomp=best.ncomp)[,,1]

### Plot the predicted_vs._true
ggplot2::ggplot(data.frame(pred_val=pred, true_val=test.y), ggplot2::aes(x=true_val, y=pred_val)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline() +
  ggplot2::xlim(1,100) +
  ggplot2::ylim(1,100) +
  ggpubr::stat_cor(method="pearson")

### Save the results for later analysis
train.y_withname=train.y; names(train.y_withname)=paste0(gsub(">>.*","",rownames(training_rpca.df)),"***",training_PBMCage_expr$condition)
test.y_withname=test.y; names(test.y_withname)=paste0(gsub(">>.*","",rownames(test_sle_expr)),"***",test_sle_expr$condition)
pred_withname=pred; names(pred_withname)=paste0(gsub(">>.*","",rownames(test_sle_expr)),"***",test_sle_expr$condition)
result_df=list(train.y=train.y_withname,
               test.y=test.y_withname,
               pred=pred_withname)

RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs$`RNA-based`[["CovidFlu_inter"]]=result_df
saveRDS(RESULT_DFs, "~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")

#####################################



### Analyze the results
#####################################
###

### Load the results
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs=RESULT_DFs$`RNA-based`

### At the donor_id level
RESULT_donor=data.frame(condition=names(RESULT_DFs$CovidFlu_donorid$test.y),
                        CA=RESULT_DFs$CovidFlu_donorid$test.y, 
                        pred.condition=names(RESULT_DFs$CovidFlu_donorid$pred),
                        BA=RESULT_DFs$CovidFlu_donorid$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_donor$condition==RESULT_donor$pred.condition)

RESULT_donor=RESULT_donor %>% select(-pred.condition) %>% 
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(condition, "_", CA_or_BA))
RESULT_donor$condition_train=forcats::fct_relevel(RESULT_donor$condition_train, c("COVR.N_CA","COVR.N_BA","COVR.Y_CA","COVR.Y_BA"))
RESULT_donor$condition=forcats::fct_relevel(RESULT_donor$condition, c("COVR.N","COVR.Y"))

plot_=
  ggplot(RESULT_donor, aes(x=condition_train, y=age, color=condition)) +
  geom_point(size=2) +
  # facet_wrap(~celltypes) +
  scale_color_manual(values=c("lightblue3","coral3")) +
  geom_line(aes(group=pair), linewidth=0.5) +
  ggpubr::stat_compare_means(aes(label=paste0()),
                             comparisons=list(c("COVR.N_CA","COVR.N_BA"), c("COVR.Y_CA","COVR.Y_BA")),
                             tip.length=0,
                             paired=T, size=3.5, bracket.size=0, 
                             label.y=65) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inCovidFlu_donor.pdf", height=5, width=2)
plot(plot_)
dev.off()

# put Covery.non-longcovid and Covery.longcovid together
RESULT_donor_modified=RESULT_donor %>% mutate(condition_train=gsub("\\.[A-Z]_","\\.",condition_train))
RESULT_donor_modified$condition_train=forcats::fct_relevel(RESULT_donor_modified$condition_train, c("COVR.CA","COVR.BA"))
plot_=
  ggplot(RESULT_donor_modified, aes(x=condition_train, y=age)) +
  geom_point(size=1, color="grey50") +
  geom_line(aes(group=pair), linewidth=0.25, color="grey50") +
  ggpubr::stat_compare_means(aes(label=paste0()),
                             comparisons=list(c("COVR.CA","COVR.BA")),
                             tip.length=0,
                             paired=T, size=3.5, bracket.size=0, 
                             label.y=65) +
  theme_light() +
  labs(x=NULL, y="age", title=" ") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none")
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inCovidFlu_donor_LongCovidandNot_mixed.pdf", height=5, width=2)
plot(plot_)
dev.off()

### At the Annot.rough level
RESULT_rough=data.frame(condition=names(RESULT_DFs$CovidFlu_rough$test.y),
                        CA=RESULT_DFs$CovidFlu_rough$test.y, 
                        pred.condition=names(RESULT_DFs$CovidFlu_rough$pred),
                        BA=RESULT_DFs$CovidFlu_rough$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_rough$condition==RESULT_rough$pred.condition)
RESULT_rough=RESULT_rough %>% select(-pred.condition) %>% 
  mutate(celltypes=gsub("\\*\\*\\*.*","",condition),
         condition=gsub(".*\\*\\*\\*","",condition)) %>%
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(condition, "_", CA_or_BA))
RESULT_rough$condition_train=forcats::fct_relevel(RESULT_rough$condition_train, c("COVR.N_CA","COVR.N_BA","COVR.Y_CA","COVR.Y_BA"))
RESULT_rough$condition=forcats::fct_relevel(RESULT_rough$condition, c("COVR.N","COVR.Y"))

plot_rough=
  ggplot(RESULT_rough, aes(x=condition_train, y=age, color=condition)) +
  geom_point(size=1) +
  facet_wrap(~celltypes, nrow=2) +
  scale_color_manual(values=c("lightblue3","coral3")) +
  geom_line(aes(group=pair), linewidth=0.25) +
  ggpubr::stat_compare_means(comparisons=list(c("COVR.N_CA","COVR.N_BA"), c("COVR.Y_CA","COVR.Y_BA")), label.y=65,
                             paired=T, size=3.5, bracket.size=0, label=after_stat("p.format")) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none",
        strip.background=element_rect(fill="transparent"),
        strip.text.x=element_text(size=10, color="black")) +
  scale_x_discrete(guide=guide_axis(n.dodge=2))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inCovidFlu_rough.pdf", height=5, width=12)
plot(plot_rough)
dev.off()

# show donors instead of celltypes and put Covery.non-longcovid and Covery.longcovid together
RESULT_DFs=readRDS("~/Project_PBMCage/For_or_From_Xu/Results/Analysis.with.AgeGene5K.model_CovidFlu.rds")
RESULT_DFs=RESULT_DFs$`RNA-based`
RESULT_rough=data.frame(condition=names(RESULT_DFs$CovidFlu_rough_onDonor$test.y),
                        CA=RESULT_DFs$CovidFlu_rough_onDonor$test.y, 
                        pred.condition=names(RESULT_DFs$CovidFlu_rough_onDonor$pred),
                        BA=RESULT_DFs$CovidFlu_rough_onDonor$pred) %>%
  mutate(donor_id=gsub(".*_","",condition),
         celltypes=gsub(">>.*","",condition)) %>%
  group_by(donor_id, celltypes) %>% # because the pseudobulk was created on the original celltype annotation, the no of samples were more than no.of.donorid*8
  summarize_at(c("CA","BA"), mean) %>%
  ungroup() %>%
  mutate(pair=rownames(.)) %>%
  mutate(condition="COVR") %>%
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(condition, ".", CA_or_BA))
RESULT_rough$condition_train=forcats::fct_relevel(RESULT_rough$condition_train, c("COVR.CA","COVR.BA"))
plot_=
  ggplot(RESULT_rough, aes(x=condition_train, y=age)) +
  geom_point(size=1, color="grey50") +
  facet_wrap(~celltypes, nrow=2) +
  geom_line(aes(group=pair), linewidth=0.25, color="grey50") +
  ggpubr::stat_compare_means(comparisons=list(c("COVR.CA","COVR.BA")), label.y=60,
                             paired=T, size=3.5, bracket.size=0, label=after_stat("p.format")) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none",
        strip.background=element_rect(fill="transparent"),
        strip.text.x=element_text(size=10, color="black"))
pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inCovidFlu_rough_LongCovidandNot_mixed.pdf", height=5, width=6)
plot(plot_)
dev.off()

### At the Annot.inter level
RESULT_inter=data.frame(condition=names(RESULT_DFs$CovidFlu_inter$test.y),
                        CA=RESULT_DFs$CovidFlu_inter$test.y, 
                        pred.condition=names(RESULT_DFs$CovidFlu_inter$pred),
                        BA=RESULT_DFs$CovidFlu_inter$pred) %>%
  mutate(pair=rownames(.))
# check
table(RESULT_inter$condition==RESULT_inter$pred.condition)
RESULT_inter=RESULT_inter %>% select(-pred.condition) %>% 
  mutate(celltypes=gsub("\\*\\*\\*.*","",condition),
         condition=gsub(".*\\*\\*\\*","",condition)) %>%
  tidyr::pivot_longer(cols=c("CA","BA"), names_to="CA_or_BA", values_to="age") %>%
  mutate(condition_train=paste0(condition, "_", CA_or_BA))
RESULT_inter$condition_train=forcats::fct_relevel(RESULT_inter$condition_train, c("COVR.N_CA","COVR.N_BA","COVR.Y_CA","COVR.Y_BA"))
RESULT_inter$condition=forcats::fct_relevel(RESULT_inter$condition, c("COVR.N","COVR.Y"))

plot_inter=
  ggplot(RESULT_inter, aes(x=condition_train, y=age, color=condition)) +
  geom_point(size=1) +
  facet_wrap(~celltypes) +
  scale_color_manual(values=c("lightblue3","coral3")) +
  geom_line(aes(group=pair), linewidth=0.25) +
  ggpubr::stat_compare_means(comparisons=list(c("COVR.N_CA","COVR.N_BA"), c("COVR.Y_CA","COVR.Y_BA")), label.y=65,
                             paired=T, size=3.5, bracket.size=0, label=after_stat("p.signif")) +
  theme_light() +
  labs(x=NULL, y="age") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_blank(),
        title=element_text(size=10),
        legend.position="none",
        strip.background=element_rect(fill="transparent"),
        strip.text.x=element_text(size=10, color="black")) +
  scale_x_discrete(guide=guide_axis(n.dodge=2))

pdf("~/Project_PBMCage/For_or_From_Xu/Xu_Plots/Model.use_inCovidFlu_inter.pdf", height=12, width=12.5)
plot(plot_inter)
dev.off()

#####################################



### Analyze the CD4T, CD8T, and NK proportions
#####################################
### 

library(Seurat)
library(UCell)
library(dplyr)

## Process the Covid data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=Seurat::NormalizeData(test_obj_pseudo_inter)
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter[[]]=test_obj_pseudo_inter[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))

### Determine CD4T.naive, CD4T.mem, CD8T.naive, CD8T.mem, and NK.CD56dim by Ucell scoring
# CD4T
CD4T_subset=subset(test_obj_pseudo_inter, Annot.rough %in% c("CD4T cells"))
CD4T_scoring=list(CD4T.naive=c("TCF7", "CD4", "CCR7", "IL7R", "FHIT","LEF1", "MAL", "NOSIP", "LDHB", "PIK3IP1"),
                  CD4T.mem=c("IL7R", "TMSB10", "CD4", "ITGB1", "LTB", "TRAC", "AQP3", "LDHB", "IL32", "MAL", "CCL5", "FYB1", "GZMK", "GZMA", "KLRB1"))
CD4T_Ucell=AddModuleScore_UCell(CD4T_subset, features=CD4T_scoring)
score_df=CD4T_Ucell[[]][, colnames(CD4T_Ucell[[]])[grepl("_UCell$", colnames(CD4T_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df=list(CD4T=score_df)
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")
# CD8T
CD8T_subset=subset(test_obj_pseudo_inter, Annot.rough %in% c("CD8T cells"))
CD8T_scoring=list(CD8T.naive=c("CD8B", "S100B", "CCR7", "RGS10", "NOSIP","LINC02446", "LEF1", "CRTAM", "CD8A", "OXNAD1"),
                  CD8T.mem=c("CD8B","ANXA1","CD8A","KRT1","LINC02446","YBX3","IL7R","TRAC","NELL2","LDHB","CCR4","CCR7","CD27","CD3D","CD3G","CD3E","SELL","CD69",
                             "KLRC1","KLRC2","KLRC3","IKZF2","CCL5","GZMH","KLRD1","TRGC2","NKG7","GZMK","CST7","GZMB","FCGR3A","HNRNPH1","XCL1","GATA3",
                             "LYAR","MALAT1","TSHZ2"))
CD8T_Ucell=AddModuleScore_UCell(CD8T_subset, features=CD8T_scoring)
score_df=CD8T_Ucell[[]][, colnames(CD8T_Ucell[[]])[grepl("_UCell$", colnames(CD8T_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df_orig=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")
score_df=c(score_df_orig, list(CD8T=score_df))
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")
# NK
NK_subset=subset(test_obj_pseudo_inter, Annot.rough %in% c("NK cells"))
NK_scoring=list(NK.CD56dim=c("GNLY", "TYROBP", "NKG7", "FCER1G", "GZMB", "TRDC", "PRF1", "FGFBP2", "SPON2", "KLRF1"),
                NK.CD56hi=c("XCL2", "FCER1G", "SPINK2", "TRDC", "KLRC1", "XCL1", "SPTSSB", "PPP1R9A", "NCAM1", "TNFRSF11A"))
NK_Ucell=AddModuleScore_UCell(NK_subset, features=NK_scoring)
score_df=NK_Ucell[[]][, colnames(NK_Ucell[[]])[grepl("_UCell$", colnames(NK_Ucell[[]]))]] %>%
  tibble::rownames_to_column("cellid")
score_df_orig=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")
score_df=c(score_df_orig, list(NK=score_df))
saveRDS(score_df,"~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")

### Arrange the scores for CD4T, CD8T, and NK
score_df=readRDS("~/Project_PBMCage/Tempt_RDS/MapTo.Covid_CD4TCD8TNK_scores.rds")
df_arranged=test_obj_pseudo_inter[[]] %>% 
  subset(Annot.rough %in% c("CD4T cells","CD8T cells","NK cells")) %>%
  dplyr::select(alt.subject.id, condition, age, Annot.rough, age_cut) %>%
  tibble::rownames_to_column("cellid") %>%
  left_join(., score_df[[1]]) %>%
  left_join(., score_df[[2]]) %>%
  left_join(., score_df[[3]])
df_arranged=df_arranged %>% subset(age_cut!="69~97") # 69~97 has only SLE INACT group but no SLE ACT or Control, so remove this

### Check the scores across age_cut and celltypes
# turns out that only CD8T.mem scores (in 48~68) showed a diff across disease status
df_selected=subset(df_arranged, Annot.rough %in% c("CD4T cells", "CD8T cells", "NK cells")) %>% 
  # dplyr::select(-c(CD4T.naive_UCell,CD4T.mem_UCell,NK.CD56dim_UCell,NK.CD56hi_UCell)) %>%
  tidyr::pivot_longer(cols=colnames(.)[grepl("UCell",colnames(.))], values_to="Score", names_to="group") %>%
  mutate(group=gsub("_UCell","",group)) %>%
  mutate(condition=ifelse(condition=="HC","CTRL","COVR"))
df_selected$group=forcats::fct_relevel(df_selected$group, c("CD4T.naive","CD4T.mem","CD8T.naive","CD8T.mem","NK.CD56dim","NK.CD56hi"))
df_selected$condition=forcats::fct_relevel(df_selected$condition, c("CTRL","COVR"))

plot_all=
  ggplot(df_selected, aes(x=Score, color=condition)) +
  facet_grid(group~age_cut) +
  geom_density(linewidth=0.25) +
  scale_color_manual(values=c("grey70","grey20")) +
  scale_x_continuous(breaks=c(0.4,0.8), labels=c(0.4,0.8)) +
  # scale_y_continuous(n.breaks=4) +
  theme_classic() +
  labs(x="UCell score", y="density", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal",
        legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(linewidth=0.5)))

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_CellProportions_CD4T.CD8T.NK.pdf", height=6.5, width=3.5)
plot_all
dev.off()

### Plot CD8T.naive scores across disease status
df_CD8T.naive=subset(df_selected, Annot.rough=="CD8T cells" & group=="CD8T.naive")
plot_=
  ggplot(df_CD8T.naive, aes(x=Score, color=condition)) +
  facet_wrap(~age_cut) +
  geom_density(linewidth=0.25) +
  scale_color_manual(values=c("grey80","grey20")) +
  scale_x_continuous(breaks=c(0.2,0.4,0.6), labels=c(0.2,0.4,0.6)) +
  theme_classic() +
  labs(x="UCell score", y="density", title=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.title=element_blank(),
        title=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(linewidth=0.5)))

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_CellProportions_CD8Tnaive.pdf", height=1.5, width=5)
plot_
dev.off()

#####################################



### Analyze the expression of metabolic- and immune-related hub genes in the CD4T & CD8T from the SLE dataset
#####################################
### 

### Load the metabolic- and immune-related hub genes from the PBMCage WGCNA modules
Modules=readRDS("~/Project_PBMCage/Tempt_RDS/WGCNA_Annot.rough_results.rds")[[4]]
Modules$module=as.character(Modules$module)
Modules=Modules %>%
  subset(module!="grey") %>%
  split(.$module)
Module_genes=Modules %>%
  lapply(., function(df) df %>% 
           arrange(desc(kME)) %>%
           dplyr::select(gene_name) %>% tibble::deframe()
  )
terms_list=c("purple"="cytoplasmic translation",
             "greenyellow"="tRNA modification",
             "yellow"="oxidative phosphorylation",
             "salmon"="chromosome localization/mitotic checkpoint",
             "turquoise"="defense response",
             "blue"="humoral immunity",
             "green"="regulation of T cell differentiation",
             "brown"="leukocyte-mediated cytotoxicity",
             "midnightblue"="regulation of hemopoiesis",
             "black"="antigen processing and presentation",
             "red"="hemostasis")
ordered_terms=c(terms_list[names(Modules)] %>% .[!is.na(.)], 
                names(Module_genes)[!names(Module_genes) %in% names(terms_list)]) %>%
  unname(.)
names(Module_genes)=ordered_terms
Merged_modules=
  list(
    metabolism=c("cytoplasmic translation","tRNA modification","chromosome localization/mitotic checkpoint",
                 "oxidative phosphorylation","regulation of hemopoiesis"),
    immunity=c("defense response","humoral immunity","regulation of T cell differentiation",
               "leukocyte-mediated cytotoxicity","antigen processing and presentation",
               "hemostasis")
  )
Module_genes=lapply(Merged_modules, function(x) Module_genes[x] %>% Reduce(union, .))

### Load the expr data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=Seurat::NormalizeData(test_obj_pseudo_inter)
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter[[]]=test_obj_pseudo_inter[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
test_sle_subset=subset(test_obj_pseudo_inter, Annot.rough %in% c("CD4T cells","CD8T cells"))
DefaultAssay(test_sle_subset)="RNA"

### Arrange the metabolic- and immunity-related df
DF_all=list()
for (i in 1:length(Merged_modules)) {
  genes.selected=Module_genes[[i]]
  expr_data_meta=FetchData(test_sle_subset, vars=c("age", "condition", "age_cut", "alt.subject.id", "Annot.rough"), layer="data")
  expr_data_all=FetchData(test_sle_subset, vars=c(genes.selected), layer="data")
  unique(rownames(expr_data_meta)==rownames(expr_data_all)) # check
  expr_data_all_filter=names(expr_data_all)[colSums(expr_data_all)>0]
  expr_data_all_arranged=expr_data_all %>% dplyr::select(all_of(expr_data_all_filter)) %>% rowMeans(.)
  expr_data_all_arranged=data.frame(cellid=names(expr_data_all_arranged), mean.expr=expr_data_all_arranged) %>%
    right_join(., expr_data_meta %>% tibble::rownames_to_column("cellid"), by="cellid") %>%
    mutate(condition=ifelse(condition=="HC","CTRL","COVR"))
  expr_data_all_arranged$condition=forcats::fct_relevel(expr_data_all_arranged$condition, c("CTRL","COVR"))
  DF_all[[i]]=expr_data_all_arranged
}
names(DF_all)=names(Merged_modules)

### Plot
# plot metabolism
plot_1=
  ggplot(DF_all[[1]], aes(x=condition, y=mean.expr, color=condition)) +
  facet_wrap(~age_cut) +
  geom_boxplot(outlier.size=0.5, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","COVR"))) +
  theme_classic() +
  scale_color_manual(values=c("grey80","grey20")) +
  labs(x=NULL, y="Avg. expression", subtitle="metabolism") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
  scale_y_continuous(expand=c(0, 0.005))
# plot immunity
plot_2=
  ggplot(DF_all[[2]], aes(x=condition, y=mean.expr, color=condition)) +
  facet_wrap(~age_cut) +
  geom_boxplot(outlier.size=0.5, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(group=age_cut, label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","COVR"))) +
  theme_classic() +
  scale_color_manual(values=c("grey80","grey20")) +
  labs(x=NULL, y="Avg. expression", subtitle="immunity") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
  scale_y_continuous(expand=c(0, 0.015))

plot_both=cowplot::plot_grid(plotlist=list(plot_1, plot_2), ncol=1)

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_Metabolism.Immune_hub.genes.expr.pdf", height=5.5, width=5)
plot_both
dev.off()

### If not split age_cut when comparing across disease status
plot_notsplit_1=
  ggplot(DF_all[[1]], aes(x=condition, y=mean.expr, color=condition)) +
  # facet_wrap(~fine_cell_type) +
  geom_boxplot(outlier.size=1, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","COVR"))) +
  theme_classic() +
  scale_color_manual(values=c("grey80","grey20")) +
  labs(x=NULL, y="Avg. expression", subtitle="metabolism") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
  scale_y_continuous(expand=c(0,0.005))

plot_notsplit_2=
  ggplot(DF_all[[2]], aes(x=condition, y=mean.expr, color=condition)) +
  # facet_wrap(~fine_cell_type) +
  geom_boxplot(outlier.size=1, outlier.shape=1, width=0.5, linewidth=0.25) +
  ggpubr::stat_compare_means(aes(label=paste0("p=",after_stat(p.format))), 
                             method="t.test", size=3, tip.length=0.01,
                             comparisons=list(c("CTRL","COVR"))) +
  theme_classic() +
  scale_color_manual(values=c("grey80","grey20")) +
  labs(x=NULL, y="Avg. expression", subtitle="immunity") +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="none",
        legend.title=element_blank(),
        plot.subtitle=element_text(size=11)) +
    scale_y_continuous(expand=c(0,0.015))

plot_both_notsplit=cowplot::plot_grid(plotlist=list(plot_notsplit_1, plot_notsplit_2), ncol=1)

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_Metabolism.Immune_hub.genes.expr_notsplitagecut.pdf", height=5.5, width=2)
plot_both_notsplit
dev.off()

#####################################



### Analyze the expression of the key kinases (JAK1/3, MAPK1/8/MAP3K7, LCK, AKT1) and TFs (NFKB familiy, MYC)
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

### Process the Covid data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter[[]]=test_obj_pseudo_inter[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
test_sle_subset=subset(test_obj_pseudo_inter, Annot.rough %in% c("CD4T cells","CD8T cells"))
test_sle_subset$donor_id=test_sle_subset$alt.subject.id
test_sle_subset_agg=AggregateExpression(test_sle_subset, assays="RNA", return.seurat=T, group.by=c("age","sex","donor_id","condition"))
test_sle_subset_agg$age=as.numeric(unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[1])))
test_sle_subset_agg$sex=unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[2]))
test_sle_subset_agg$donor_id=unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[3]))
test_sle_subset_agg$condition=unlist(lapply(strsplit(colnames(test_sle_subset_agg), split="_"), function(x) x[4])) %>% gsub("\\..*","",.)

### Extract the gene expression
genes_of_interest=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","AKT1","LCK","MYC","NFKB1","NFKB2")
expr_data=FetchData(test_sle_subset_agg, vars=c(genes_of_interest, colnames(test_sle_subset_agg[[]])), layer="data")
expr_data_arranged=expr_data %>%
  dplyr::select(-c(orig.ident)) %>%
  tidyr::pivot_longer(cols=colnames(.)[!grepl("age|sex|donor_id|condition",colnames(.))], names_to="gene", values_to="data")
# sort the gene expr
genes_of_interest_list=list(MYC=c("MYC"), 
                            `Metabolism-related kinases`=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","AKT1"),
                            NFKB2=c("NFKB2"),
                            `Immunity-related kinases`=c("JAK1","JAK3","MAPK1","MAPK8","MAP3K7","LCK"))
mean_expr_data=
  lapply(1:length(genes_of_interest_list), function(idx) expr_data_arranged %>% 
           subset(gene %in% genes_of_interest_list[[idx]]) %>%
           group_by(donor_id, age, condition) %>%
           summarize_at("data", mean) %>%
           mutate(group=names(genes_of_interest_list)[idx])) %>%
  data.table::rbindlist(.) %>%
  mutate(condition=ifelse(condition=="HC","CTRL","COVR")) %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
mean_expr_data$group=forcats::fct_relevel(mean_expr_data$group, c("MYC","NFKB2","Metabolism-related kinases","Immunity-related kinases"))
mean_expr_data$condition=forcats::fct_relevel(mean_expr_data$condition, c("CTRL","COVR"))

### Plot
plot_TFKinase=
  ggplot(mean_expr_data, aes(x=age, y=data, color=condition)) +
  facet_wrap(~group, scales="free_y") +
  geom_point(size=0.5, shape=1) +
  geom_smooth(aes(group=condition), method="lm", se=F, linewidth=0.25, show.legend=F) +
  theme_classic() +
  scale_color_manual(values=c("grey80","grey20")) +
  labs(x="age", y="Avg. expression", subtitle=NULL) +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical",
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.subtitle=element_text(size=11)) +
  guides(color=guide_legend(override.aes=list(size=2, shape=19)))

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_Metabolism.Immune_TFsAndKinases.expr.pdf", height=3.5, width=5)
plot_TFKinase
dev.off()

#####################################



### Analyze the communication among cell types
#####################################
###

library(dplyr)
library(ggplot2)
library(Seurat)

### Load the expr data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=Seurat::NormalizeData(test_obj_pseudo_inter)
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter[[]]=test_obj_pseudo_inter[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97"))))
DefaultAssay(test_obj_pseudo_inter)="RNA"

### Extract the key interactions
B_CD8T=list(`B cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"), 
            `CD8T cells`=c("CD8B"))
allother_NK=list(`B cells,CD4T cells,Monocytes,CD8T cells,NK cells`=c("HLA-E"), 
                 `NK cells`=c("KLRK1","KLRD1","KLRC1"))
CD8T_NK=list(`CD8T cells`=c("CD99"), `NK cells`=c("CD99"))
NK_CD8T=list(`NK cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"), `CD8T cells`=c("CD8A"))
CD4T_CD8T=list(`CD4T cells`=c("MIF"), `CD8T cells`=c("CD44","CD74","CXCR4"))

all_to_CD8T=list(`B cells`=c("HLA-A","HLA-B","HLA-C","HLA-E"), `CD8T cells`=c("CD8B"))
all_to_NK=list(`B cells,CD4T cells,Monocytes,CD8T cells,NK cells`=c("HLA-E"), 
               `NK cells`=c("KLRK1","KLRD1","KLRC1"))

### Plot the signals to CD8T and NK, respectively
molecule_list=list(all_to_CD8T, all_to_NK)
plot_lists=list()
for (j in 1:length(molecule_list)) {
  DF_all=lapply(1:length(molecule_list[[j]]), function(idx) FetchData(test_obj_pseudo_inter %>% 
                                                                        subset(Annot.rough %in% strsplit(names(molecule_list[[j]][idx]), split=",")[[1]]), 
                                                                      vars=c(molecule_list[[j]][[idx]], "age", "age_cut", "condition", "alt.subject.id"), 
                                                                      layer="data"))
  tempt=lapply(1:length(DF_all), function(idx) DF_all[[idx]] %>% dplyr::select(any_of(molecule_list[[j]][[idx]])) %>% rowSums())
  tempt=lapply(1:length(tempt), function(idx) data.frame(cellid=names(tempt[[idx]]), expr=tempt[[idx]]))
  DF_all_mean=lapply(1:length(DF_all), function(idx) {
    tempt_=DF_all[[idx]] %>% 
      dplyr::select(all_of(c("age", "age_cut", "condition", "alt.subject.id"))) %>%
      tibble::rownames_to_column("cellid") %>%
      left_join(., tempt[[idx]]) %>%
      group_by(alt.subject.id, age, age_cut, condition) %>%
      summarize_at("expr", mean)
    colnames(tempt_)[ncol(tempt_)]=names(molecule_list[[j]][idx])
    tempt_
  }) %>%
    Reduce(full_join, .) %>%
    mutate(condition=ifelse(condition=="HC","CTRL","COVR"))
  colnames(DF_all_mean)[c(5,6)]=c("interaction_A","interaction_B")
  DF_all_mean$condition=forcats::fct_relevel(DF_all_mean$condition, c("CTRL","COVR"))
  
  plot_lists[[j]]=
    ggplot(DF_all_mean, aes(x=interaction_A, y=interaction_B, color=condition)) +
    geom_point(size=0.5, shape=1) +
    geom_smooth(aes(group=condition), show.legend=F, se=F, method="lm", linewidth=0.5) +
    ggpmisc::stat_correlation(mapping=ggpmisc::use_label(c("R", "P"), aes(group=condition)), 
                              small.p=T,
                              show.legend=F, size=3, vstep=0.05) +
    theme_classic() +
    scale_color_manual(values=c("grey80","grey20")) +
    labs(x="Avg. expr. of MHC-I", y=c("Avg. expr. of CD8B","Avg. expr. of NK receptor")[j], 
         subtitle=c("Signals to CD8T cells","Signals to NK cells")[j]) +
    theme(axis.text.x=element_text(size=9),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position=c("none","right")[j],
          legend.title=element_blank(),
          legend.text=element_text(size=9),
          plot.subtitle=element_text(size=11)) +
    guides(color=guide_legend(override.aes=list(size=2, shape=19)))
}

plot_both=cowplot::plot_grid(plotlist=plot_lists, align="h", axis="tb", rel_widths=c(0.65,1))

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_CellchatInteractions.expr.pdf", height=5.5, width=5)
plot_both
dev.off()

#####################################



### Analyze the expression patterns
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

### Load the expr data
test_obj_pseudo_inter=readRDS("~/Project_PBMCage/GSE158055/New_covid_flu_Dataset_Pseudobulk_inter.rds")
test_obj_pseudo_inter=Seurat::NormalizeData(test_obj_pseudo_inter)
test_obj_pseudo_inter=subset(test_obj_pseudo_inter, condition!="NA")
test_obj_pseudo_inter[[]]=test_obj_pseudo_inter[[]] %>%
  mutate(age_cut=ifelse(age<=30,"19~30",
                        ifelse(age>30 & age<=47,"31~47",
                               ifelse(age>47 & age<=68, "48~68", "69~97")))) %>%
  mutate(sample_id=ifelse(grepl("HC",rownames(.)), paste0("HC",gsub(".*_HC","",rownames(.))), rownames(.))) %>%
  mutate(sample_id=ifelse(grepl("COVR",rownames(.)), paste0("COVR",gsub(".*_COVR","",rownames(.))), sample_id))
DefaultAssay(test_obj_pseudo_inter)="RNA"
test_sle_aggr=AggregateExpression(test_obj_pseudo_inter, return.seurat=T, group.by=c("sample_id","Annot.rough"))
updated_meta=test_sle_aggr[[]] %>% 
  mutate(sample_id=gsub("_.*","",rownames(.)), Annot.rough=gsub(".*_","",rownames(.))) %>%
  left_join(., test_obj_pseudo_inter[[]] %>% dplyr::select(sample_id, age, age_cut, alt.subject.id, condition, Annot.rough) %>% subset(!duplicated(.))) %>%
  mutate(id=paste0(sample_id,"_",Annot.rough)) %>%
  tibble::column_to_rownames("id")
test_sle_aggr[[]]=updated_meta

### Load the 5 expression patterns from the PBMCage project
pattern_trend_genes_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.txt.gz", sep="\t")
celltype_map=unique(pattern_trend_genes_df$celltype)
pattern_df_modified=pattern_trend_genes_df %>%
  mutate(celltype=ifelse(celltype=="B","B cells",
                         ifelse(celltype=="CD4T","CD4T cells",
                                ifelse(celltype=="CD8T","CD8T cells",
                                       ifelse(celltype=="Mono","Monocytes",
                                              ifelse(celltype=="NK","NK cells",
                                                     ifelse(celltype=="DC","DCs",
                                                            ifelse(celltype=="OtherT","OtherT","Other cells")))))))) %>%
  subset(celltype!="Monocytes") %>% # because the Covid dataset doesn't have monocytes
  split(.$pattern)

### Extract the mean.expr for each pattern
expr_all_patterns=list()
for (i in 1:length(pattern_df_modified)) {
  celltype_list=unique(pattern_df_modified[[i]]$celltype)
  expr_for_celltype_list=list()
  for (j in 1:length(celltype_list)) {
    gene_per_cell_per_pattern=pattern_df_modified[[i]] %>% 
      subset(celltype==celltype_list[j]) %>% dplyr::select(gene) %>% tibble::deframe()
    cellname=strsplit(celltype_list[j], split=",")[[1]]
    expr_data=
      FetchData(test_sle_aggr %>% subset(Annot.rough %in% cellname), vars=c(gene_per_cell_per_pattern, "age","age_cut","condition","alt.subject.id"),
                layer="data")
    # remove the genes with all zeros
    expr_data_filter=colSums(expr_data %>% .[,colnames(.)[!grepl("age|condition|subject",colnames(.))]]) %>% .[.>0] %>% names(.)
    expr_data_cleaned=cbind(expr_data[, c(expr_data_filter)], 
                            expr_data[, c("age","age_cut","condition","alt.subject.id")])
    expr_data_cleaned=expr_data_cleaned %>% group_by(age, age_cut, condition, alt.subject.id) %>%
      summarize_at(colnames(.)[!grepl("age|condition|subject",colnames(.))], mean) %>%
      tidyr::pivot_longer(cols=colnames(.)[!grepl("age|condition|subject",colnames(.))],
                          values_to="mean.expr", names_to="gene") %>%
      mutate(celltype=celltype_list[j])
    expr_for_celltype_list[[j]]=expr_data_cleaned
  }
  all_expr_for_pattern=data.table::rbindlist(expr_for_celltype_list) %>%
    mutate(condition=ifelse(condition=="HC","CTRL","COVR"))
  all_expr_for_pattern$condition=forcats::fct_relevel(all_expr_for_pattern$condition, c("CTRL","COVR"))
  
  expr_all_patterns[[i]]=all_expr_for_pattern
}
saveRDS(expr_all_patterns, "~/Project_PBMCage/Map_to_Covid/ExpressionPatterns_mean.expr.per.genes.rds")

### Plot
expr_all_patterns=readRDS("~/Project_PBMCage/Map_to_Covid/ExpressionPatterns_mean.expr.per.genes.rds")
plot_lists=list()
for (i in 1:length(expr_all_patterns)) {
  quantile_lim=expr_all_patterns[[i]] %>% group_by(age_cut, condition) %>% 
    mutate(median=quantile(mean.expr, 0.5),
           upper_hinge=quantile(mean.expr, 0.75),
           lower_hinge=quantile(mean.expr, 0.25)) %>%
    mutate(IQR=upper_hinge-lower_hinge) %>%
    mutate(upper_whisker=upper_hinge+1.5*IQR,
           lower_whisker=lower_hinge-1.5*IQR) %>%
    dplyr::select(age_cut, condition, upper_hinge, lower_hinge, IQR, upper_whisker, lower_whisker) %>%
    subset(!duplicated(.))
    
  plot_lists[[i]]=
    ggplot(expr_all_patterns[[i]], aes(x=age_cut, y=mean.expr, color=condition)) +
    facet_wrap(~condition) +
    geom_boxplot(outlier.shape=NA, outlier.size=0.5, width=0.5, fill="transparent", coef=0, linewidth=0.25) +
    scale_color_manual(values=c("grey80","grey20")) +
    theme_classic() +
    labs(x=NULL, y="Avg. expression", subtitle=paste0("Pattern",i)) +
    theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position="none",
          legend.title=element_blank(),
          plot.subtitle=element_text(size=11)) +
    coord_cartesian(ylim=c(min(quantile_lim$lower_hinge), max(quantile_lim$upper_hinge)))
}
plot_patterns_all=cowplot::plot_grid(plotlist=plot_lists, axis="tb", align="h", nrow=1)

pdf("~/Project_PBMCage/Plots/Map_to_Covid/MapToCovid_ExprPatterns.pdf", height=1.75, width=18)
plot_patterns_all
dev.off()

#####################################



### Analyze the timing patterns
#####################################
###

library(Seurat)
library(dplyr)
library(ggplot2)

#####################################
