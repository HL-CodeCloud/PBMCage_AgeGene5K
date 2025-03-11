library(dplyr)
library(Seurat)

### Window-sliding to determine the agecut
#####################################
###
###

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
celltypes_r=names(table(THEObj$Annot.rough))
celltypes_i=names(table(THEObj$Annot.inter))
celltypes_d=names(table(THEObj$Annot.detailed))

# Identify the age-associated sig. genes in the whole obj
DEGs_total=list()
names_=c()
Idents(THEObj)="age"
age_=names(table(THEObj$age))
for (i in 1:(length(age_)-1)) {
  tryCatch({
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i])
    DEGs_total=c(DEGs_total, list(markers))
    names_=c(names_, paste0("All", "_", age_[i+1]))
  }, error=function(msg) {print("Skip.")})
}
names(DEGs_total)=names_

# Identify the age-associated sig. genes in celltypes at the rough level
DEGs_rough=list()
names_=c()
for (j in 1:length(celltypes_r)) {
  subset_=subset(THEObj, Annot.rough==celltypes_r[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_rough=c(DEGs_rough, list(markers))
      names_=c(names_, paste0(celltypes_r[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_rough)=names_

# Identify the age-associated sig. genes in celltypes at the inter level
DEGs_inter=list()
names_=c()
for (j in 1:length(celltypes_i)) {
  subset_=subset(THEObj, Annot.inter==celltypes_i[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_inter=c(DEGs_inter, list(markers))
      names_=c(names_, paste0(celltypes_i[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_inter)=names_

# Identify the age-associated sig. genes in celltypes at the detailed level
DEGs_detailed=list()
names_=c()
for (j in 1:length(celltypes_d)) {
  subset_=subset(THEObj, Annot.detailed==celltypes_d[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")

#####################################



### Window-sliding to determine the agecut in females
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
THEObj=subset(THEObj, sex=="female")
celltypes_r=names(table(THEObj$Annot.rough))
celltypes_i=names(table(THEObj$Annot.inter))
celltypes_d=names(table(THEObj$Annot.detailed))

# Identify the age-associated sig. genes in the whole obj
DEGs_total=list()
names_=c()
Idents(THEObj)="age"
age_=names(table(THEObj$age))
for (i in 1:(length(age_)-1)) {
  tryCatch({
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i])
    DEGs_total=c(DEGs_total, list(markers))
    names_=c(names_, paste0("All", "_", age_[i+1]))
  }, error=function(msg) {print("Skip.")})
}
names(DEGs_total)=names_

# Identify the age-associated sig. genes in celltypes at the rough level
DEGs_rough=list()
names_=c()
for (j in 1:length(celltypes_r)) {
  subset_=subset(THEObj, Annot.rough==celltypes_r[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_rough=c(DEGs_rough, list(markers))
      names_=c(names_, paste0(celltypes_r[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_rough)=names_

# Identify the age-associated sig. genes in celltypes at the inter level
DEGs_inter=list()
names_=c()
for (j in 1:length(celltypes_i)) {
  subset_=subset(THEObj, Annot.inter==celltypes_i[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_inter=c(DEGs_inter, list(markers))
      names_=c(names_, paste0(celltypes_i[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_inter)=names_

# Identify the age-associated sig. genes in celltypes at the detailed level
DEGs_detailed=list()
names_=c()
for (j in 1:length(celltypes_d)) {
  subset_=subset(THEObj, Annot.detailed==celltypes_d[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly.rds")

#####################################



### Window-sliding to determine the agecut in males
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
THEObj=subset(THEObj, sex=="male")
celltypes_r=names(table(THEObj$Annot.rough))
celltypes_i=names(table(THEObj$Annot.inter))
celltypes_d=names(table(THEObj$Annot.detailed))

# Identify the age-associated sig. genes in the whole obj
DEGs_total=list()
names_=c()
Idents(THEObj)="age"
age_=names(table(THEObj$age))
for (i in 1:(length(age_)-1)) {
  tryCatch({
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i])
    DEGs_total=c(DEGs_total, list(markers))
    names_=c(names_, paste0("All", "_", age_[i+1]))
  }, error=function(msg) {print("Skip.")})
}
names(DEGs_total)=names_

# Identify the age-associated sig. genes in celltypes at the rough level
DEGs_rough=list()
names_=c()
for (j in 1:length(celltypes_r)) {
  subset_=subset(THEObj, Annot.rough==celltypes_r[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_rough=c(DEGs_rough, list(markers))
      names_=c(names_, paste0(celltypes_r[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_rough)=names_

# Identify the age-associated sig. genes in celltypes at the inter level
DEGs_inter=list()
names_=c()
for (j in 1:length(celltypes_i)) {
  subset_=subset(THEObj, Annot.inter==celltypes_i[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_inter=c(DEGs_inter, list(markers))
      names_=c(names_, paste0(celltypes_i[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_inter)=names_

# Identify the age-associated sig. genes in celltypes at the detailed level
DEGs_detailed=list()
names_=c()
for (j in 1:length(celltypes_d)) {
  subset_=subset(THEObj, Annot.detailed==celltypes_d[j])
  Idents(subset_)="age"
  age_=names(table(subset_$age))
  for (i in 1:(length(age_)-1)) {
    tryCatch({
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i])
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly.rds")

#####################################



### Clean the results
#####################################
###

library(dplyr)

###
The_Whole_List=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")
DEGs_total=The_Whole_List[[1]]
DEGs_rough=The_Whole_List[[2]]
DEGs_inter=The_Whole_List[[3]]
DEGs_detailed=The_Whole_List[[4]]

### Clean DEGs_total
Count_Df=list()
for (i in 1:length(DEGs_total)) {
  result_per_age=DEGs_total[[i]]
  result_per_age_pos_nrow=result_per_age %>%
    subset(avg_log2FC>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_neg_nrow=result_per_age %>%
    subset(avg_log2FC<=(-2) & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_both_nrow=result_per_age %>%
    subset(abs(avg_log2FC)>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  count_df=data.frame(Category=names(DEGs_total)[i],
                      pos_count=result_per_age_pos_nrow,
                      neg_count=result_per_age_neg_nrow,
                      DEG_count=result_per_age_both_nrow)
  count_df$celltype=gsub("_[0-9].*","",count_df$Category)
  count_df$age=as.numeric(gsub("[^0-9]","",count_df$Category))
  Count_Df[[i]]=count_df
}
Count_Df=data.table::rbindlist(Count_Df)

### Clean DEGs_rough
Count_Df_r=list()
for (i in 1:length(DEGs_rough)) {
  result_per_age=DEGs_rough[[i]]
  result_per_age_pos_nrow=result_per_age %>%
    subset(avg_log2FC>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_neg_nrow=result_per_age %>%
    subset(avg_log2FC<=(-2) & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_both_nrow=result_per_age %>%
    subset(abs(avg_log2FC)>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  count_df=data.frame(Category=names(DEGs_rough)[i],
                      pos_count=result_per_age_pos_nrow,
                      neg_count=result_per_age_neg_nrow,
                      DEG_count=result_per_age_both_nrow)
  count_df$celltype=gsub("_[0-9].*","",count_df$Category)
  count_df$age=as.numeric(gsub(".*_","",count_df$Category))
  Count_Df_r[[i]]=count_df
}
Count_Df_r=data.table::rbindlist(Count_Df_r)

### Clean DEGs_inter
Count_Df_i=list()
for (i in 1:length(DEGs_inter)) {
  result_per_age=DEGs_inter[[i]]
  result_per_age_pos_nrow=result_per_age %>%
    subset(avg_log2FC>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_neg_nrow=result_per_age %>%
    subset(avg_log2FC<=(-2) & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_both_nrow=result_per_age %>%
    subset(abs(avg_log2FC)>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  count_df=data.frame(Category=names(DEGs_inter)[i],
                      pos_count=result_per_age_pos_nrow,
                      neg_count=result_per_age_neg_nrow,
                      DEG_count=result_per_age_both_nrow)
  count_df$celltype=gsub("_[0-9].*","",count_df$Category)
  count_df$age=as.numeric(gsub(".*_","",count_df$Category))
  Count_Df_i[[i]]=count_df
}
Count_Df_i=data.table::rbindlist(Count_Df_i)

### Clean DEGs_detailed
Count_Df_d=list()
for (i in 1:length(DEGs_detailed)) {
  result_per_age=DEGs_detailed[[i]]
  result_per_age_pos_nrow=result_per_age %>%
    subset(avg_log2FC>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_neg_nrow=result_per_age %>%
    subset(avg_log2FC<=(-2) & p_val_adj<=0.05) %>%
    nrow(.)
  result_per_age_both_nrow=result_per_age %>%
    subset(abs(avg_log2FC)>=2 & p_val_adj<=0.05) %>%
    nrow(.)
  count_df=data.frame(Category=names(DEGs_detailed)[i],
                      pos_count=result_per_age_pos_nrow,
                      neg_count=result_per_age_neg_nrow,
                      DEG_count=result_per_age_both_nrow)
  count_df$celltype=gsub("_[0-9].*","",count_df$Category)
  count_df$age=as.numeric(gsub(".*_","",count_df$Category))
  Count_Df_d[[i]]=count_df
}
Count_Df_d=data.table::rbindlist(Count_Df_d)

### Save
saveRDS(list(Count_Df, Count_Df_r, Count_Df_i, Count_Df_d), "~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_counts.rds")

#####################################



### Plot the Window sliding results
#####################################
###

library(ggplot2)

###
The_Whole_List=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_counts.rds")
Count_Df=The_Whole_List[[1]]
Count_Df_r=The_Whole_List[[2]]
Count_Df_i=The_Whole_List[[3]]
Count_Df_d=The_Whole_List[[4]]

# Plot PBMC
plot_all=
  ggplot(Count_Df, aes(x=age, y=DEG_count)) +
  geom_point() +
  labs(y="# of DEGs") +
  scale_x_continuous(breaks=seq(0,120,2)) +
  scale_y_continuous(trans='log10') +
  geom_line(aes(group=1), linetype="dashed", show.legend=FALSE, linewidth=0.2, alpha=0.3) +
  geom_smooth(show.legend=FALSE, linewidth=0.75, se=FALSE, color="darkblue") +
  scale_color_d3() +
  theme_classic() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  guides(color=guide_legend(title="Celltype"))

# Plot rough level
plot_rough=
  ggplot(Count_Df_r, aes(x=age, y=DEG_count, color=celltype)) +
    geom_point(size=1, alpha=0.3) +
    labs(y="# of DEGs") +
    scale_x_continuous(breaks=seq(0,120,2)) +
    scale_y_continuous(trans='log10') +
    geom_line(aes(group=celltype), linetype="dashed", show.legend=FALSE, linewidth=0.2, alpha=0.3) +
    geom_smooth(aes(group=celltype, color=celltype), show.legend=FALSE, linewidth=0.75, se=FALSE) +
    scale_color_d3() +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9)) +
    guides(color=guide_legend(title="Celltype", override.aes=list(size=2, alpha=1)))

# Plot inter level
Count_Df_i$celltype_rough=gsub("\\..*","",Count_Df_i$celltype)
plot_inter=
  ggplot(Count_Df_i, aes(x=age, y=DEG_count, color=celltype_rough)) +
    geom_point(size=1, alpha=0.3) +
    labs(y="# of DEGs") +
    scale_x_continuous(breaks=seq(0,120,2)) +
    scale_y_continuous(trans='log10') +
    geom_line(aes(group=celltype), linetype="dashed", show.legend=FALSE, linewidth=0.2, alpha=0.3) +
    geom_smooth(aes(group=celltype_rough, color=celltype_rough), show.legend=FALSE, linewidth=0.75, se=FALSE) +
    scale_color_d3() +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9)) +
    guides(color=guide_legend(title="Celltype", override.aes=list(size=2, alpha=1)))

# Plot detailed level
Count_Df_d$celltype_rough=gsub("\\..*","",Count_Df_d$celltype)
plot_detailed=
  ggplot(Count_Df_d, aes(x=age, y=DEG_count, color=celltype_rough)) +
    geom_point(size=1, alpha=0.3) +
    labs(y="# of DEGs") +
    scale_x_continuous(breaks=seq(0,120,2)) +
    scale_y_continuous(trans='log10') +
    geom_line(aes(group=celltype), linetype="dashed", show.legend=FALSE, linewidth=0.2, alpha=0.3) +
    geom_smooth(aes(group=celltype_rough, color=celltype_rough), method="loess", show.legend=FALSE, linewidth=0.75, se=FALSE) +
    scale_color_d3() +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9)) +
    guides(color=guide_legend(title="Celltype", override.aes=list(size=2, alpha=1)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGsCounts.pdf", height=9, width=8)
plot(plot_all)
plot(plot_rough)
plot(plot_inter)
plot(plot_detailed)
dev.off()

#####################################



### Vote for the agecuts according to the derivatives of loess
#####################################
###

library(ggplot2)

### Vote based on the DEGs
The_Whole_List=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_counts.rds")
Count_Df=The_Whole_List[[1]]
Count_Df_r=The_Whole_List[[2]]
Count_Df_i=The_Whole_List[[3]]
Count_Df_d=The_Whole_List[[4]]

# for PBMC
df=Count_Df
lo_function=loess(DEG_count ~ age, subset(df, celltype=="All"), span=0.5, degree=2)
predictedY=predict(lo_function, c(min(df$age):max(df$age)))
predict_df=data.frame(age=c(min(df$age):max(df$age)),
                      curve_y=predictedY)
predict_df$diff=c(diff(predict_df$curve_y)[1],diff(predict_df$curve_y))
Age_At_Turn=data.frame(celltype=rep("All", nrow(subset(predict_df, abs(diff)<1))),
                       age=subset(predict_df, abs(diff)<1)$age)

# for Rough Level
df=Count_Df_r
age_at_turn_=list()
for (i in 1:length(unique(df$celltype))) {
  tryCatch({
    lo_function=loess(DEG_count ~ age, subset(df, celltype==unique(df$celltype)[i]), span=0.5, degree=2)
    predictedY=predict(lo_function, c(min(df$age):max(df$age)))
    predict_df=data.frame(age=c(min(df$age):max(df$age)),
                          curve_y=predictedY)
    predict_df$diff=c(diff(predict_df$curve_y)[1],diff(predict_df$curve_y))
    predict_df$Sec.diff=c(1, sapply(2:length(predict_df$diff), function(x) predict_df$diff[x]/predict_df$diff[x-1]))
    tempt_=data.frame(celltype=rep(unique(df$celltype)[i], nrow(subset(predict_df, Sec.diff<0))),
                      age=subset(predict_df, Sec.diff<0)$age)
  }, error=function(msg) print("Skip."))
  age_at_turn_=c(age_at_turn_, list(tempt_))
}
Age_At_Turn_rough=data.table::rbindlist(age_at_turn_)

# for inter Level
df=Count_Df_i
age_at_turn_=list()
for (i in 1:length(unique(df$celltype))) {
  tryCatch({
    lo_function=loess(DEG_count ~ age, subset(df, celltype==unique(df$celltype)[i]), span=0.5, degree=2)
    predictedY=predict(lo_function, c(min(df$age):max(df$age)))
    predict_df=data.frame(age=c(min(df$age):max(df$age)),
                          curve_y=predictedY)
    predict_df$diff=c(diff(predict_df$curve_y)[1],diff(predict_df$curve_y))
    predict_df$Sec.diff=c(1, sapply(2:length(predict_df$diff), function(x) predict_df$diff[x]/predict_df$diff[x-1]))
    tempt_=data.frame(celltype=rep(unique(df$celltype)[i], nrow(subset(predict_df, Sec.diff<0))),
                      age=subset(predict_df, Sec.diff<0)$age)
  }, error=function(msg) print("Skip."))
  age_at_turn_=c(age_at_turn_, list(tempt_))
}
Age_At_Turn_inter=data.table::rbindlist(age_at_turn_)

# for detailed Level
df=Count_Df_d
age_at_turn_=list()
for (i in 1:length(unique(df$celltype))) {
  tryCatch({
    lo_function=loess(DEG_count ~ age, subset(df, celltype==unique(df$celltype)[i]), span=0.5, degree=2)
    predictedY=predict(lo_function, c(min(df$age):max(df$age)))
    predict_df=data.frame(age=c(min(df$age):max(df$age)),
                          curve_y=predictedY)
    predict_df$diff=c(diff(predict_df$curve_y)[1],diff(predict_df$curve_y))
    predict_df$Sec.diff=c(1, sapply(2:length(predict_df$diff), function(x) predict_df$diff[x]/predict_df$diff[x-1]))
    tempt_=data.frame(celltype=rep(unique(df$celltype)[i], nrow(subset(predict_df, Sec.diff<0))),
                      age=subset(predict_df, Sec.diff<0)$age)
  }, error=function(msg) print("Skip."))
  age_at_turn_=c(age_at_turn_, list(tempt_))
}
Age_At_Turn_detailed=data.table::rbindlist(age_at_turn_)

# plot the DEG turning points
Age_At_Turn_all=data.table::rbindlist(list(Age_At_Turn, Age_At_Turn_rough, Age_At_Turn_inter, Age_At_Turn_detailed))
x_most=Age_At_Turn_all %>% group_by(age) %>% dplyr::count() %>% ungroup() %>% slice_max(order_by=n, n=20) %>% dplyr::select(age) %>% unlist() %>% unname()
plot_mode_DEGs=
  ggplot(Age_At_Turn_all, aes(x=age, y=celltype)) +
  geom_point() +
  geom_vline(xintercept=x_most) +
  labs(y=NULL, title="DEGs") +
  scale_x_continuous(breaks=seq(0,120,2)) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

###
### Vote based on the celltype freq
pbmc.seu_merged=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly.rds")
ID_=pbmc.seu_merged$donor_id
Age_=pbmc.seu_merged$age
CellType_=pbmc.seu_merged$Annot.detailed
cellNum=table(CellType_, ID_)
cellProp=t(t(cellNum)/rowSums(t(cellNum)))
data=t(cellProp*100)
data=as.data.frame(data)
colnames(data)=c("ID","CellType","Freq")
ID_Age_match=data.frame(ID=ID_, Age=Age_)
ID_Age_match=ID_Age_match[!duplicated(ID_Age_match),]
data=merge(data, ID_Age_match, by="ID")
data$CellType=as.factor(data$CellType)
data=subset(data, !grepl("lin_infidel|pltcontamin", CellType))
data$Age=as.integer(as.character(data$Age))

# for detailed Level
age_at_turn_=list()
for (i in 1:length(unique(data$CellType))) {
  tryCatch({
    lo_function=loess(Freq ~ Age, subset(data, CellType==unique(data$CellType)[i]), span=0.5, degree=2)
    predictedY=predict(lo_function, c(min(data$Age):max(data$Age)))
    predict_df=data.frame(age=c(min(data$Age):max(data$Age)),
                          curve_y=predictedY)
    predict_df$diff=c(diff(predict_df$curve_y)[1],diff(predict_df$curve_y))
    predict_df$Sec.diff=c(1, sapply(2:length(predict_df$diff), function(x) predict_df$diff[x]/predict_df$diff[x-1]))
    tempt_=data.frame(celltype=rep(unique(data$CellType)[i], nrow(subset(predict_df, Sec.diff<0))),
                      age=subset(predict_df, Sec.diff<0)$age)
  }, error=function(msg) print("Skip."))
  age_at_turn_=c(age_at_turn_, list(tempt_))
}
Age_At_Turn_detailed=data.table::rbindlist(age_at_turn_)

# plot the freq turning points
x_most=Age_At_Turn_detailed %>% group_by(age) %>% dplyr::count() %>% ungroup() %>% slice_max(order_by=n, n=20) %>% dplyr::select(age) %>% unlist() %>% unname()
plot_mode_freq=
  ggplot(Age_At_Turn_detailed, aes(x=age, y=celltype)) +
  geom_point() +
  geom_vline(xintercept=x_most) +
  labs(y=NULL, title="Freq") +
  scale_x_continuous(breaks=seq(0,120,2)) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

###
### Vote based on the both the DEGs and celltype freq
Age_At_Turn_both=data.table::rbindlist(list(Age_At_Turn_all, Age_At_Turn_detailed))
# plot the freq turning points
x_most=Age_At_Turn_both %>% group_by(age) %>% dplyr::count() %>% ungroup() %>% slice_max(order_by=n, n=20) %>% dplyr::select(age) %>% unlist() %>% unname()
plot_mode_both=
  ggplot(Age_At_Turn_detailed, aes(x=age, y=celltype)) +
  geom_point() +
  geom_vline(xintercept=x_most) +
  labs(y=NULL, title="DEGs and Freq") +
  scale_x_continuous(breaks=seq(0,120,2)) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9))

###
### Plot the vote
# add colored blocks of PBMCagecuts
AGECUT=list(c(19:24),c(25:30),c(31:42),c(43:51),c(52:58),c(59:64),c(65:73),c(74:81),c(82:88),c(89:91),c(92:97))
colors=paletteer::paletteer_d("ggsci::category20_d3")[1:length(AGECUT)]
codes=""
for (i in 1:length(AGECUT)) {
  codes=paste0(codes, 
               "+ annotate('rect', xmin=",AGECUT[[i]][1],", xmax=",AGECUT[[i]][length(AGECUT[[i]])],", ymin=-Inf, ymax=Inf, fill='",colors[i],"', alpha=0.3)")
}
plot_mode_both_blocks=
  eval(parse(text=paste0("plot_mode_both",codes,
                         "+ scale_x_continuous(breaks=seq(0,120,4))"
                         )))

# start plotting
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGsCounts.and.CelltypeFreq_Mode.pdf", height=10, width=8)
plot(plot_mode_DEGs)
plot(plot_mode_freq)
plot(plot_mode_both)
plot(plot_mode_both_blocks)
dev.off()

#####################################



### Analyze the variation of DEGs across ages
#####################################
###

library(dplyr)

### Load the results
The_Whole_List=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")
The_Whole_List_F=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly.rds")
The_Whole_List_M=readRDS("~/Project_PBMCage/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly.rds")

### Merge the data at the donor_id level
DEGs_total=The_Whole_List[[1]]; DEGs_total_F=The_Whole_List_F[[1]]; DEGs_total_M=The_Whole_List_M[[1]]
DEGs_total_clean=lapply(1:length(DEGs_total), function(i) DEGs_total[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total)[i])), 
                                                                                     gene=rownames(.),
                                                                                     celltype.level="donor_id",
                                                                                     celltypes=gsub("_.*","",names(DEGs_total)[i]),
                                                                                     analysis="all"))
DEGs_total_clean_F=lapply(1:length(DEGs_total_F), function(i) DEGs_total_F[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_F)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="donor_id",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_F)[i]),
                                                                                           analysis="F"))
DEGs_total_clean_M=lapply(1:length(DEGs_total_M), function(i) DEGs_total_M[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_M)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="donor_id",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_M)[i]),
                                                                                           analysis="M"))
DEGs_total_clean=data.table::rbindlist(c(DEGs_total_clean, DEGs_total_clean_F, DEGs_total_clean_M)) %>% 
  select(gene, age, avg_log2FC, p_val_adj, celltype.level, celltypes, analysis) %>% 
  mutate(orientation=ifelse(avg_log2FC>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_log2F=abs(avg_log2FC))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Donorid.csv.gz", sep="\t")

### Merge the data at the rough level
DEGs_total=The_Whole_List[[2]]; DEGs_total_F=The_Whole_List_F[[2]]; DEGs_total_M=The_Whole_List_M[[2]]
unique(names(DEGs_total)) # check
unique(gsub("_.*","",names(DEGs_total))) # check
DEGs_total_clean=lapply(1:length(DEGs_total), function(i) DEGs_total[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total)[i])), 
                                                                                     gene=rownames(.),
                                                                                     celltype.level="rough",
                                                                                     celltypes=gsub("_.*","",names(DEGs_total)[i]),
                                                                                     analysis="all"))
DEGs_total_clean_F=lapply(1:length(DEGs_total_F), function(i) DEGs_total_F[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_F)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="rough",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_F)[i]),
                                                                                           analysis="F"))
DEGs_total_clean_M=lapply(1:length(DEGs_total_M), function(i) DEGs_total_M[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_M)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="rough",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_M)[i]),
                                                                                           analysis="M"))
DEGs_total_clean=data.table::rbindlist(c(DEGs_total_clean, DEGs_total_clean_F, DEGs_total_clean_M)) %>% 
  select(gene, age, avg_log2FC, p_val_adj, celltype.level, celltypes, analysis) %>% 
  mutate(orientation=ifelse(avg_log2FC>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_log2F=abs(avg_log2FC))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")

### Merge the data at the inter level
DEGs_total=The_Whole_List[[3]]; DEGs_total_F=The_Whole_List_F[[3]]; DEGs_total_M=The_Whole_List_M[[3]]
unique(names(DEGs_total)) # check
unique(gsub("_.*","",names(DEGs_total))) # check
DEGs_total_clean=lapply(1:length(DEGs_total), function(i) DEGs_total[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total)[i])), 
                                                                                     gene=rownames(.),
                                                                                     celltype.level="inter",
                                                                                     celltypes=gsub("_.*","",names(DEGs_total)[i]),
                                                                                     analysis="all"))
DEGs_total_clean_F=lapply(1:length(DEGs_total_F), function(i) DEGs_total_F[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_F)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="inter",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_F)[i]),
                                                                                           analysis="F"))
DEGs_total_clean_M=lapply(1:length(DEGs_total_M), function(i) DEGs_total_M[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_M)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="inter",
                                                                                           celltypes=gsub("_.*","",names(DEGs_total_M)[i]),
                                                                                           analysis="M"))
DEGs_total_clean=data.table::rbindlist(c(DEGs_total_clean, DEGs_total_clean_F, DEGs_total_clean_M)) %>% 
  select(gene, age, avg_log2FC, p_val_adj, celltype.level, celltypes, analysis) %>% 
  mutate(orientation=ifelse(avg_log2FC>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_log2F=abs(avg_log2FC))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Inter.csv.gz", sep="\t")

### Merge the data at the detailed level
DEGs_total=The_Whole_List[[4]]; DEGs_total_F=The_Whole_List_F[[4]]; DEGs_total_M=The_Whole_List_M[[4]]
unique(names(DEGs_total)) # check
unique(gsub("_[0-9][0-9].*","",names(DEGs_total))) # check
DEGs_total_clean=lapply(1:length(DEGs_total), function(i) DEGs_total[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total)[i])), 
                                                                                     gene=rownames(.),
                                                                                     celltype.level="detailed",
                                                                                     celltypes=gsub("_[0-9][0-9].*","",names(DEGs_total)[i]),
                                                                                     analysis="all"))
DEGs_total_clean_F=lapply(1:length(DEGs_total_F), function(i) DEGs_total_F[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_F)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="detailed",
                                                                                           celltypes=gsub("_[0-9][0-9].*","",names(DEGs_total_F)[i]),
                                                                                           analysis="F"))
DEGs_total_clean_M=lapply(1:length(DEGs_total_M), function(i) DEGs_total_M[[i]] %>% mutate(age=as.numeric(gsub(".*_","",names(DEGs_total_M)[i])), 
                                                                                           gene=rownames(.),
                                                                                           celltype.level="detailed",
                                                                                           celltypes=gsub("_[0-9][0-9].*","",names(DEGs_total_M)[i]),
                                                                                           analysis="M"))
DEGs_total_clean=data.table::rbindlist(c(DEGs_total_clean, DEGs_total_clean_F, DEGs_total_clean_M)) %>% 
  select(gene, age, avg_log2FC, p_val_adj, celltype.level, celltypes, analysis) %>% 
  mutate(orientation=ifelse(avg_log2FC>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_log2F=abs(avg_log2FC))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Detailed.csv.gz", sep="\t")

### Analyze the variation at the donor_id level
DEGs_Donorid=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Donorid.csv.gz", sep="\t")
DEGs_Donorid_bothsex=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Donorid_bothsex=DEGs_Donorid_bothsex %>% group_by(age, orientation) %>% summarize_at("avg_log2FC", mean)
mixed_sex_plot=
  ggplot(DEGs_Donorid_bothsex, aes(x=age, y=avg_log2FC, color=orientation))+
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c("brown3","dodgerblue3"), labels=c("Increasing", "Decreasing")) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  # geom_segment(aes(x=34, xend=34, y=1.45, yend=1.3), linewidth=0.1, arrow=arrow(length=unit(0.04, "npc")), color="coral4") +
  # geom_segment(aes(x=60, xend=60, y=1.45, yend=1.3), linewidth=0.1, arrow=arrow(length=unit(0.04, "npc")), color="coral4") +
  # geom_segment(aes(x=78, xend=78, y=1.45, yend=1.3), linewidth=0.1, arrow=arrow(length=unit(0.04, "npc")), color="coral4") +
  geom_smooth(aes(group=orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation.pdf", height=4, width=2.5)
plot(mixed_sex_plot)
dev.off()

### Analyze the variation at the donor_id level, comparing F vs. M
DEGs_Donorid_subset=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Donorid_subset=DEGs_Donorid_subset %>% group_by(age, analysis, orientation) %>% summarize_at("avg_log2FC", mean)
DEGs_Donorid_subset=DEGs_Donorid_subset %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))
FvsM_plot=
  ggplot(DEGs_Donorid_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_FvsM.pdf", height=4, width=2.5)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the donor_id level, comparing the 12 superclusters
DEGs_Donorid=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Donorid.csv.gz", sep="\t")
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Donorid_subset=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Donorid_subset=lapply(1:length(Genes_in_SuperClusters), 
                           function(idx) DEGs_Donorid_subset %>% 
                             subset(gene %in% Genes_in_SuperClusters[[idx]]) %>% 
                             group_by(age, analysis, orientation) %>% summarize_at("avg_log2FC", mean) %>% 
                             mutate(supercluster=gsub("\\|.*","",names(Genes_in_SuperClusters))[idx]) %>%
                             mutate(supercluster_orientation=paste0(supercluster,"_",orientation))
                           ) %>%
  data.table::rbindlist(.)
supecluster_=names(table(DEGs_Donorid_subset$supercluster))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
supecluster_reorder=supecluster_[supecluster_reorder_idx]
DEGs_Donorid_subset$supercluster=forcats::fct_relevel(DEGs_Donorid_subset$supercluster, supecluster_reorder)
supercluster_plot=
  ggplot(DEGs_Donorid_subset, aes(x=age, y=avg_log2FC, color=supercluster))+
  facet_wrap(~analysis, labeller=as_labeller(c('F'="Female",'M'="Male"))) +
  geom_point(alpha=0.1, size=0.1) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=supercluster_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL, nrow=12))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_12superclusters.pdf", height=3.5, width=5.5)
plot(supercluster_plot)
dev.off()

### Analyze the variation at the rough level
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_bothsex=DEGs_Rough %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% group_by(age, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
# show only lymphocytes as in the previous text we focused on their shared alterations
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT"))
# DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% subset(!celltypes %in% c("DCs","Other cells")) # have checked, these two do not have enough data points, so excluded
# mixed_sex_plot=
  ggplot(DEGs_Rough_bothsex, aes(x=age, y=avg_log2FC, color=orientation))+
  facet_wrap(~celltypes) +
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c("grey60","grey30"), labels=c("Higher than avg.", "Lower than avg.")) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  scale_x_continuous(breaks=c(30,60,90)) +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Rough.pdf", height=5, width=4.5)
plot(mixed_sex_plot)
dev.off()

### Analyze the variation at the rough level in lymphocytes, comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(analysis!="all") %>% 
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation))
DEGs_Rough_selected=DEGs_Rough %>% 
  subset(analysis!="all" & p_val_adj<0.05) %>% 
  # show only lymphocytes as in the previous text we focused on their shared alterations
  # subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")) %>%
  subset(!(celltypes %in% c("DCs","Other cells"))) %>%
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation))

FvsM_plot_selected=
  ggplot(DEGs_Rough_selected, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, nrow=2) +
  geom_point(alpha=0.1, size=0.1) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.25) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Rough.Lymphocytes_FvsM.pdf", height=4, width=4)
plot(FvsM_plot_selected)
dev.off()
  
FvsM_plot=
  ggplot(DEGs_Rough_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=4) +
  geom_hline(yintercept=c(0.25, -0.25), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change), title="Yazar et al. (2022) dataset") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  # geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  # geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Rough.All_FvsM.pdf", height=3.5, width=6)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the rough level in DCs, Monocytes, and Other cells, comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(analysis!="all")
DEGs_Rough_subset=DEGs_Rough_subset %>% 
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
# select the celltypes of interest
DEGs_Rough_selected=DEGs_Rough_subset %>% subset(celltypes %in% c("CD4T cells","DCs","Monocytes","Other cells"))
DEGs_Rough_selected=DEGs_Rough_selected %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))
FvsM_plot=
  ggplot(DEGs_Rough_selected, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=4) +
  geom_hline(yintercept=c(0.5, -0.5), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Rough.DC.Mono.Othercells_FvsM.pdf", height=2, width=6)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the inter level (celltypes with enough samples for geom_smooth), comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Inter.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(analysis!="all" & p_val_adj<0.05)
DEGs_Rough_subset=DEGs_Rough_subset %>% group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
# take only CD4 and CD8T that are interesting
DEGs_Rough_subset=DEGs_Rough_subset %>% subset(grepl("CD4|CD8",celltypes)) %>% subset(celltypes!="CD4T.Treg")
DEGs_Rough_subset=DEGs_Rough_subset %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))

FvsM_plot=
  ggplot(DEGs_Rough_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=6) +
  geom_hline(yintercept=c(0.5, -0.5), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Inter.CD4CD8T_FvsM.pdf", height=2, width=8)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the inter level (celltypes without enough samples for geom_smooth), comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Inter.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(grepl("CD4T\\.naive|CD4T\\.Treg|CD4T\\.prolif|CD8T\\.CTL",celltypes)) 
DEGs_Rough_pcutoff=DEGs_Rough_subset %>% subset(analysis!="all")
DEGs_Rough_pcutoff=DEGs_Rough_pcutoff %>% group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
DEGs_Rough_pcutoff=DEGs_Rough_pcutoff %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))

FvsM_plot=
  ggplot(DEGs_Rough_pcutoff, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=6) +
  coord_cartesian(ylim=c(-4,4)) +
  geom_hline(yintercept=c(0.5, -0.5), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr.Cor.wAges_CorPlot_Tprolif.Treg.CTL_Inter.pdf", height=2, width=5)
plot(FvsM_plot)
dev.off()


### Analyze the variation at the detailed level, comparing F vs. M
# calculate AUC of avg.log2FC * age
DEGs_Detailed=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Detailed.csv.gz", sep="\t")

AUC_results=DEGs_Detailed %>% subset(p_val_adj<0.05 & analysis=="all") %>% split(.$celltypes) %>%
  lapply(., function(df) sum(diff(df$age[order(df$age)])*zoo::rollmean(abs(df$avg_log2FC)[order(df$age)],2))) %>%
  as.data.frame() %>% t(.) %>% as.data.frame() %>% mutate(celltypes=rownames(.), AUC=V1, analysis="all") %>% select(-V1)

AUC_F=DEGs_Detailed %>% subset(p_val_adj<0.05 & analysis=="F") %>% split(.$celltypes) %>%
  lapply(., function(df) sum(diff(df$age[order(df$age)])*zoo::rollmean(abs(df$avg_log2FC)[order(df$age)],2))) %>%
  as.data.frame() %>% t(.) %>% as.data.frame() %>% mutate(celltypes=rownames(.), AUC=V1, analysis="F") %>% select(-V1)

AUC_M=DEGs_Detailed %>% subset(p_val_adj<0.05 & analysis=="M") %>% split(.$celltypes) %>%
  lapply(., function(df) sum(diff(df$age[order(df$age)])*zoo::rollmean(abs(df$avg_log2FC)[order(df$age)],2))) %>%
  as.data.frame() %>% t(.) %>% as.data.frame() %>% mutate(celltypes=rownames(.), AUC=V1, analysis="M") %>% select(-V1)

AUC_DF=data.table::rbindlist(list(AUC_results, AUC_F, AUC_M))

# clean the data
AUC_DF_plot=AUC_DF %>% 
  tidyr::pivot_wider(names_from="analysis", values_from="AUC") %>% 
  mutate(diff=`F`-`M`, rough=gsub("\\..*","",celltypes), FoverM=ifelse(`F`>`M`,"F>M","F<M")) %>%
  arrange(desc(abs(diff))) %>%
  na.omit(.)
AUC_DF_plot$celltypes=forcats::fct_relevel(AUC_DF_plot$celltypes, unique(AUC_DF_plot$celltypes))
x_length=unlist(lapply(split(AUC_DF_plot$celltypes, AUC_DF_plot$rough), length))
range(AUC_DF_plot$diff)

plotlists=AUC_DF_plot %>% split(~rough) %>% 
  purrr::map(., 
             ~ggplot(., aes(x=celltypes, y=diff, color=FoverM)) +
               geom_point() +
               theme_classic() +
               labs(x=NULL, y="Sex difference\nin AUC") +
               theme(axis.text.y=element_text(size=11),
                     axis.title.y=element_text(size=12),
                     axis.text.x=element_text(size=12, angle=60, vjust=1, hjust=1),
                     legend.position="none",
                     legend.direction="horizontal",
                     plot.margin=margin(t=0.5, l=1.5, r=0.02, unit="cm")) +
               coord_cartesian(ylim=c(-250,250)) +
               scale_y_continuous(labels=abs) +
               geom_hline(yintercept=0, linetype="dashed", linewidth=0.4, color="grey70") +
               scale_color_manual(values=c("dodgerblue3","coral3")) +
               guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
             )
plotlists[c(5,7,4,1)]=plotlists[c(5,7,4,1)] %>% 
  purrr::map(~ .x + theme(axis.text.y=element_blank(),
                          axis.title.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.line.y=element_blank(),
                          plot.margin=margin(t=0.5, l=0.2, r=0.02, unit="cm"))
             ) 
all_plots=
  cowplot::plot_grid(cowplot::plot_grid(plotlist=plotlists[c(2,5)], align="h", axis="tb", rel_widths=c(8,2), ncol=2),
                     cowplot::plot_grid(plotlist=plotlists[c(3,7,4)], align="h", axis="tb", rel_widths=c(9,2,1), ncol=3),
                     cowplot::plot_grid(plotlist=plotlists[c(6,1)], align="h", axis="tb", rel_widths=c(10,7), ncol=2),
                     axis="rl", align="hv", nrow=3, rel_heights=c(1.05,1,1.1))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_Detailed_FvsM.pdf", height=10, width=5)
plot(all_plots)
dev.off()

#####################################



### Analyze the variation of DEGs across ages in the key terms (12 superclusters)
#####################################
###

### Load the genes in the superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")

### Take the supercluster genes in DEGs data
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(analysis=="all" & avg_log2FC<Inf & avg_log2FC>-Inf & p_val_adj<0.05) %>% 
  mutate(celltypes=gsub("\\..*","",celltypes))
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, orientation, celltypes, superclusters) %>% 
                                      summarize_at("avg_log2FC", mean) %>% 
                                      mutate(analysis_orientation=paste0(analysis,"_",orientation))) %>%
  data.table::rbindlist()
DEGs_Rough_subset_final=DEGs_Rough_subset_percluster %>%
  # mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)
supecluster_=names(table(DEGs_Rough_subset_final$superclusters))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
supecluster_reorder=supecluster_[supecluster_reorder_idx]
DEGs_Rough_subset_final$superclusters=forcats::fct_relevel(DEGs_Rough_subset_final$superclusters, supecluster_reorder)

# test
for(i in 1:12) {
  print(paste0("==========",i,"=========="))
  print(
    ggpubr::compare_means(avg_log2FC_abs ~ celltypes,  
                          data=DEGs_Rough_subset_final %>% mutate(avg_log2FC_abs=abs(avg_log2FC)) %>%
                            subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT") & 
                                     superclusters==supecluster_reorder[i]),
                          ref.group="CD4T cells",
                          method="wilcox.test")
  )
  print(
    ggpubr::compare_means(avg_log2FC_abs ~ celltypes,  
                          data=DEGs_Rough_subset_final %>% mutate(avg_log2FC_abs=abs(avg_log2FC)) %>%
                            subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT") & 
                                     superclusters==supecluster_reorder[i]),
                          ref.group="CD8T cells",
                          method="wilcox.test")
  )
}
# ...turns out to be "leuk. devel." and "immune proc." that are sig. diff. between "CD4 & CD8" and others
FvsM_plot=
  ggplot(DEGs_Rough_subset_final %>% subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")), 
         aes(x=agecut, y=abs(avg_log2FC), color=celltypes))+
  facet_wrap(~superclusters) +
  coord_cartesian(ylim=c(0,2)) +
  ggrastr::geom_beeswarm_rast(shape=1, alpha=0.1, size=0.5) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(1,2,3,6,8)]) +
  labs(x=NULL, y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_line(linewidth=0.02, color="grey80"),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="top",
        legend.direction="horizontal") +
  geom_line(data=. %>% group_by(celltypes, superclusters, agecut) %>% summarize_at("avg_log2FC", function(x) mean(abs(x), na.rm=T)),
            aes(group=celltypes), linewidth=0.5, show.legend=F) +
  geom_text(data=. %>% subset(superclusters %in% c("leuk. devel.","immune proc.") & agecut=="19~30" & celltypes=="CD4T cells"), 
            aes(x="31~47", y=1.5, label="p<0.05"), color="black", size=3.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_MapToSupercluster.pdf", height=6, width=8)
FvsM_plot
dev.off()

Immune_plot=
  ggplot(DEGs_Rough_subset_final %>% 
           subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes") &
                    superclusters %in% c("leuk. devel.","immune proc.")), 
         aes(x=agecut, y=abs(avg_log2FC), color=celltypes))+
  facet_wrap(~superclusters) +
  coord_cartesian(ylim=c(0,2)) +
  ggrastr::geom_beeswarm_rast(shape=1, alpha=0.1, size=0.5) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(2,3,5)]) +
  labs(x=NULL, y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_line(linewidth=0.02, color="grey80"),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="top",
        legend.direction="horizontal") +
  geom_line(data=. %>% group_by(celltypes, superclusters, agecut) %>% summarize_at("avg_log2FC", function(x) mean(abs(x), na.rm=T)),
            aes(group=celltypes), linewidth=0.5, show.legend=F) +
  ggpubr::stat_compare_means(data=. %>% 
                               group_by(analysis, orientation, celltypes, superclusters) %>% 
                               summarize_at("avg_log2FC",.funs=function(x) {mean(abs(x))}) %>%
                               mutate(agecut="31~47"),
                             method="anova",
                             aes(group=celltypes), label="p.format", label.y=1.75) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_MapToSupercluster_leukoAndImmune_Only.pdf", height=2.5, width=4.5)
Immune_plot
dev.off()

#####################################



### Analyze the relation between variation of genes within 12 superclusters and Annot.detailed cell counts
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Detailed.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2F, na.rm=T)) %>%
  unlist()
DEGs_df=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                   supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                   sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()

# check the relation between variation and cell counts
# ... turns out that variation in a certain cell type is negatively correlated with its cell count
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
cell_count=meta_info_combined %>% dplyr::select(donor_id, Annot.detailed, N_per_detailed) %>%
  subset(!duplicated(.)) %>%
  group_by(Annot.detailed) %>% summarize_at("N_per_detailed", function(x) sum(x, na.rm=T)) %>%
  rename(celltype=Annot.detailed) %>%
  right_join(., DEGs_df, by=c("celltype"))

cell_count_reorder=names(table(cell_count$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cell_count$supercluster=forcats::fct_relevel(cell_count$supercluster, cell_count_reorder)
plot_=
  ggplot(cell_count, aes(x=N_per_detailed, y=sd, color=supercluster)) +
    geom_point(alpha=0.25) +
    geom_smooth(aes(group=supercluster), se=F, linewidth=0.5) +
    ggpubr::stat_cor(aes(group=1), method="spearman", show.legend=F) +
    scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")) +
    theme_classic() +
    labs(x="Cell count", y=expression(SD~of~Avg.Abs.Log[2]~Fold-change)) +
    theme(axis.text.y=element_text(size=9),
          axis.text.x=element_text(size=9),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          legend.position="right",
          legend.direction="horizontal") +
    guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL, nrow=6))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_sd.vs.cellcounts_Detailed.pdf", height=2.5, width=5)
plot(plot_)
dev.off()

#####################################



### Normalize the variation by cell counts (sd.RMS="Integral of sd across N" / N)
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
# # remove outliers
# outliers=lapply(DEGs_Rough_subset_percluster %>% split(.$superclusters), 
#                 function(df) {boxplot.stats(df$abs.avg_log2F)$out}
# ) %>% unlist() %>% as.data.frame() %>%
#   tibble::rownames_to_column("supercluster") %>%
#   rename(abs.avg_log2FC=".") %>% 
#   mutate(supercluster=gsub("[0-9]$","",supercluster)) %>%
#   mutate(supercluster_abs.avgLOG2FC=paste0(supercluster,"_",abs.avg_log2FC))
# DEGs_Rough_subset_percluster=DEGs_Rough_subset_percluster %>%
#   mutate(supercluster_abs.avgLOG2FC=paste0(superclusters,"_",abs.avg_log2F)) %>%
#   subset(!(supercluster_abs.avgLOG2FC %in% outliers$supercluster_abs.avgLOG2FC))
# calculate sd for abs.avg_log2FC in each supercluster
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2F, na.rm=T)) %>%
  unlist()
DEGs_df=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                   supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                   sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()
# add the cell count info
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
cell_count=meta_info_combined %>% dplyr::select(donor_id, Annot.rough, N_per_rough) %>%
  subset(!duplicated(.)) %>%
  group_by(Annot.rough) %>% summarize_at("N_per_rough", function(x) sum(x, na.rm=T)) %>%
  rename(celltype=Annot.rough) %>%
  right_join(., DEGs_df, by=c("celltype"))
# normalize cell count by sd.RMS
sd.div.N_persupercluster=cell_count %>% 
  subset(!is.na(sd)) %>%
  split(.$supercluster) %>%
  lapply(., 
         function(df) {
           integral_=with(df, 
                          integrate(approxfun(N_per_rough, sd), 
                                    lower=min(N_per_rough), 
                                    upper=max(N_per_rough))) %>%
             .$value
           divident_=max(df$N_per_rough)-min(df$N_per_rough)
           integral_/divident_
         }) %>%
  unlist() %>% as.data.frame() %>%
  dplyr::rename(RMS=".") %>%
  tibble::rownames_to_column("supercluster")

# reorder_=names(table(sd.div.N_persupercluster$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
reorder_=sd.div.N_persupercluster %>% arrange(RMS) %>% select(supercluster) %>% 
  tibble::deframe() %>%
  as.character()
sd.div.N_persupercluster$supercluster=forcats::fct_relevel(sd.div.N_persupercluster$supercluster, reorder_)

plot_=
ggplot(sd.div.N_persupercluster, aes(x=supercluster, y=RMS)) +
  geom_bar(stat="identity", fill="grey80", color="transparent", width=0.75) +
  theme_classic() +
  labs(x="Supercluster", y=expression(Mean~of~SD~of~Avg.Abs.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
  # coord_flip(ylim=c(0.2,0.6)) +
  coord_flip() +
  scale_x_discrete(guide = guide_axis(n.dodge=2))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_RMS_superclusters.pdf", height=2, width=7)
plot(plot_)
dev.off()

#####################################



### Plot the variation across age in 12 superclusters in CD8T
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()

# rearrange the df
DEGs_df=DEGs_Rough_subset_percluster %>%
  subset(celltypes=="CD8T cells") %>%
  # mutate(superclusters_orientation=paste0(superclusters, "_", orientation)) %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                                       "cell activity",
                                       ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                                              "cell reaction")))
supecluster_=names(table(DEGs_df$superclusters))
supecluster_reorder=supecluster_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
DEGs_df$superclusters=forcats::fct_relevel(DEGs_df$superclusters, supecluster_reorder)
DEGs_df$main_group=forcats::fct_relevel(DEGs_df$main_group, 
                                        c("subcellular bioprocesses","cell activity","cell reaction"))

plot_var=
  ggplot(DEGs_df, aes(x=age, y=abs.avg_log2F, color=superclusters, linetype=main_group)) +
  geom_point(size=0.05, alpha=0.05) +
  geom_smooth(aes(group=superclusters), se=F, linewidth=0.5) +
  scale_linetype_manual(values=c("longdash","solid","dotted")) +
  coord_cartesian(ylim=c(0,1.6)) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")) +
  labs(x="age", y=expression(Avg.Abs.Log[2]~Fold-change)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=9, angle=0),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=9, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="right",
        legend.direction="vertical", 
        legend.spacing.y=unit(-0.5, "mm")
        ) +
  # scale_y_continuous(labels={function(x) abs(x)}) +
  guides(color=guide_legend(title=NULL, ncol=1, byrow=F),
         linetype=guide_legend(title=NULL, override.aes=list(color="black")))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_pattern_superclusters.pdf", height=4, width=4)
plot(plot_var)
dev.off()

### Plot the AUC of variation across age in 12 superclusters in CD8T
# calculate the AUC
DEGs_df_AUC=DEGs_df %>% 
  group_by(age, celltypes, superclusters, main_group) %>%
  summarize_at("abs.avg_log2F", mean) %>%
  split(.$superclusters) %>%
  lapply(., function(df)
    {integral_=with(df, 
                   integrate(approxfun(age, abs.avg_log2F), 
                             lower=min(age), 
                             upper=max(age), subdivisions=300)) %>%
      .$value})
# arrange the df
DEGs_df_AUC_arranged=DEGs_df_AUC %>%
  unlist() %>%
  as.data.frame() %>%
  rename(AUC.abs.avg_log2FC=".") %>%
  tibble::rownames_to_column("superclusters") %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                  "cell activity",
                  ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                         "cell reaction")))
# reorder
# reorder_=names(table(DEGs_df_AUC_arranged$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
reorder_=DEGs_df_AUC_arranged %>% arrange(AUC.abs.avg_log2FC) %>% select(superclusters) %>%
  tibble::deframe() %>%
  as.character()
DEGs_df_AUC_arranged$superclusters=forcats::fct_relevel(DEGs_df_AUC_arranged$superclusters, reorder_)

DEGs_df_AUC_arranged$main_group=forcats::fct_relevel(DEGs_df_AUC_arranged$main_group, 
                                                     c("subcellular bioprocesses","cell activity","cell reaction"))
plot_AUC=
ggplot(DEGs_df_AUC_arranged, aes(x=AUC.abs.avg_log2FC, y=superclusters, fill=main_group)) +
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.75) +
  scale_fill_manual(values=c("white","grey80","grey40")) +
  labs(x=NULL, y=expression(AUC~of~Avg.Abs.Log[2]~Fold-change)) +
  theme_classic() +
  theme(plot.background=element_rect(fill="transparent", color="transparent"),
        axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=9, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="bottom",
        legend.direction="horizontal", 
        legend.spacing.y=unit(-0.5, "mm")
  ) +
  guides(fill=guide_legend(title=NULL, override.aes=list(size=1)))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_AUC_superclusters_CD8T.pdf", height=4.2, width=2.5)
plot(plot_AUC)
dev.off()

#####################################



### Plot the variation across age in 12 superclusters in general (regardless of cell types)
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()

# rearrange the df
DEGs_df=DEGs_Rough_subset_percluster %>%
  # subset(celltypes=="CD8T cells") %>%
  # mutate(superclusters_orientation=paste0(superclusters, "_", orientation)) %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                  "cell activity",
                  ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                         "cell reaction")))
supecluster_=names(table(DEGs_df$superclusters))
supecluster_reorder=supecluster_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
DEGs_df$superclusters=forcats::fct_relevel(DEGs_df$superclusters, supecluster_reorder)
DEGs_df$main_group=forcats::fct_relevel(DEGs_df$main_group, 
                                        c("subcellular bioprocesses","cell activity","cell reaction"))

### Plot the AUC of variation across age in 12 superclusters in general
# calculate the AUC
DEGs_df_AUC=DEGs_df %>% 
  group_by(age, celltypes, superclusters, main_group) %>%
  summarize_at("abs.avg_log2F", mean) %>%
  split(.$superclusters) %>%
  lapply(., function(df)
  {integral_=with(df, 
                  integrate(approxfun(age, abs.avg_log2F), 
                            lower=min(age), 
                            upper=max(age), subdivisions=300)) %>%
    .$value})
# arrange the df
DEGs_df_AUC_arranged=DEGs_df_AUC %>%
  unlist() %>%
  as.data.frame() %>%
  rename(AUC.abs.avg_log2FC=".") %>%
  tibble::rownames_to_column("superclusters") %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                  "cell activity",
                  ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                         "cell reaction")))
# reorder
# reorder_=names(table(DEGs_df_AUC_arranged$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
reorder_=DEGs_df_AUC_arranged %>% arrange(AUC.abs.avg_log2FC) %>% select(superclusters) %>%
  tibble::deframe() %>%
  as.character()
DEGs_df_AUC_arranged$superclusters=forcats::fct_relevel(DEGs_df_AUC_arranged$superclusters, reorder_)

DEGs_df_AUC_arranged$main_group=forcats::fct_relevel(DEGs_df_AUC_arranged$main_group, 
                                                     c("subcellular bioprocesses","cell activity","cell reaction"))
plot_AUC=
  ggplot(DEGs_df_AUC_arranged, aes(x=AUC.abs.avg_log2FC, y=superclusters)) +
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.75) +
  # scale_fill_manual(values=c("white","grey80","grey40")) +
  labs(y=NULL, x=expression(AUC~of~Avg.Abs.Log[2]~Fold-change)) +
  theme_void() +
  theme(plot.background=element_rect(fill="transparent", color="transparent"),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10, margin=margin(t=0.5, b=0.5, unit="line")),
        axis.text.y=element_text(size=9, hjust=1, margin=margin(r=-1, unit="line")),
        axis.title.y=element_text(size=10),
        legend.position="bottom",
        legend.direction="horizontal"
  ) +
  guides(fill=guide_legend(title=NULL, override.aes=list(size=1)))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_AUC_superclusters_allAnnot.rough.pdf", height=3.5, width=5)
plot(plot_AUC)
dev.off()

#####################################



### Variation in 12 superclusters in CD8T
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")

### CD8T
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(celltypes=="CD8T cells") %>%
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
# remove outliers
outliers=lapply(DEGs_Rough_subset_percluster %>% split(.$superclusters), 
       function(df) {boxplot.stats(df$abs.avg_log2F)$out}
  ) %>% unlist() %>% as.data.frame() %>%
  tibble::rownames_to_column("supercluster") %>%
  rename(abs.avg_log2FC=".") %>% 
  mutate(supercluster=gsub("[0-9]$","",supercluster)) %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(supercluster,"_",abs.avg_log2FC))
DEGs_Rough_subset_percluster=DEGs_Rough_subset_percluster %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(superclusters,"_",abs.avg_log2F)) %>%
  subset(!(supercluster_abs.avgLOG2FC %in% outliers$supercluster_abs.avgLOG2FC))
# calculate sd for abs.avg_log2FC in each supercluster
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2F, na.rm=T)) %>%
  unlist()
DEGs_df1=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                   supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                   sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()

### CD4T
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(celltypes=="CD4T cells") %>%
  subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
# remove outliers
outliers=lapply(DEGs_Rough_subset_percluster %>% split(.$superclusters), 
                function(df) {boxplot.stats(df$abs.avg_log2F)$out}
) %>% unlist() %>% as.data.frame() %>%
  tibble::rownames_to_column("supercluster") %>%
  rename(abs.avg_log2FC=".") %>% 
  mutate(supercluster=gsub("[0-9]$","",supercluster)) %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(supercluster,"_",abs.avg_log2FC))
DEGs_Rough_subset_percluster=DEGs_Rough_subset_percluster %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(superclusters,"_",abs.avg_log2F)) %>%
  subset(!(supercluster_abs.avgLOG2FC %in% outliers$supercluster_abs.avgLOG2FC))
# calculate sd for abs.avg_log2FC in each supercluster
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2F, na.rm=T)) %>%
  unlist()
DEGs_df2=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                   supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                   sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()

### Merge CD4T and CD8T sd
DEGs_df=rbind(DEGs_df1, DEGs_df2)
# reorder
# reorder_=names(table(DEGs_df$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
reorder_=DEGs_df %>% subset(celltype=="CD8T cells") %>% arrange(sd) %>% select(supercluster) %>% 
  tibble::deframe() %>%
  as.character()
DEGs_df$supercluster=forcats::fct_relevel(DEGs_df$supercluster, reorder_)

plot_=
  ggplot(DEGs_df, aes(x=supercluster, y=sd)) +
  facet_wrap(~celltype) +
  geom_bar(stat="identity", fill="grey80", color="transparent", width=0.75) +
  theme_classic() +
  labs(x="Supercluster", y="SD") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
  # coord_flip(ylim=c(0.2,0.4)) +
  coord_flip() +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  scale_y_continuous(n.breaks=4)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_RNA_DEGs_variation_superclusters_CD8T.pdf", height=2, width=4)
plot(plot_)
dev.off()

#####################################



### Calculate the avg_log2FC with window-sliding on the TSS methylation data
#####################################
###

library(dplyr)
library(Seurat)

### Load the data
beta_meta_merged_ALL=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")
all_celltypes=unique(beta_meta_merged_ALL$celltypes)

### Window-sliding analysis
results_per_age_ALL=list()
for (i in 1:length(all_celltypes)) {
  the_celltype=all_celltypes[i]
  results_per_age_both=list()
  
  # for females
  beta_meta_subset=beta_meta_merged_ALL %>% subset(celltypes==the_celltype & sex=="F")
  beta_meta_subset_wide=beta_meta_subset %>% select(gene, age, beta) %>% group_by(gene, age) %>%
    summarize_at("beta", mean) %>%
    tidyr::pivot_wider(names_from="age", values_from="beta")
  results_per_age=lapply(3:ncol(beta_meta_subset_wide), function(x) {
    df=data.frame(beta_meta_subset_wide[,x-1], beta_meta_subset_wide[,x]); colnames(df)=c("former","latter"); df %>% mutate(avg_log2FC=log2(latter+1e-5)-log2(former+1e-5)) %>% select(avg_log2FC) %>% unlist() %>% unname()
  })
  names(results_per_age)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
  results_per_age_df=as.data.frame(results_per_age)
  rownames(results_per_age_df)=beta_meta_subset_wide$gene
  colnames(results_per_age_df)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
  results_per_age_df=results_per_age_df %>% 
    mutate(gene=rownames(.)) %>%
    tidyr::pivot_longer(cols=colnames(results_per_age_df), names_to="age", values_to="avg_log2FC") %>%
    mutate(celltypes=the_celltype, analysis="F")
  results_per_age_both[[1]]=results_per_age_df
  
  # for males
  beta_meta_subset=beta_meta_merged_ALL %>% subset(celltypes==the_celltype & sex=="M")
  beta_meta_subset_wide=beta_meta_subset %>% select(gene, age, beta) %>% group_by(gene, age) %>%
    summarize_at("beta", mean) %>%
    tidyr::pivot_wider(names_from="age", values_from="beta")
  results_per_age=lapply(3:ncol(beta_meta_subset_wide), function(x) {
    df=data.frame(beta_meta_subset_wide[,x-1], beta_meta_subset_wide[,x]); colnames(df)=c("former","latter"); df %>% mutate(avg_log2FC=log2(latter+1e-5)-log2(former+1e-5)) %>% select(avg_log2FC) %>% unlist() %>% unname()
  })
  names(results_per_age)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
  results_per_age_df=as.data.frame(results_per_age)
  rownames(results_per_age_df)=beta_meta_subset_wide$gene
  colnames(results_per_age_df)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
  results_per_age_df=results_per_age_df %>% 
    mutate(gene=rownames(.)) %>%
    tidyr::pivot_longer(cols=colnames(results_per_age_df), names_to="age", values_to="avg_log2FC") %>%
    mutate(celltypes=the_celltype, analysis="M")
  results_per_age_both[[2]]=results_per_age_df
  
  # merge the sex
  results_per_age_ALL[[i]]=data.table::rbindlist(results_per_age_both)
}

results_ALL=data.table::rbindlist(results_per_age_ALL)
data.table::fwrite(results_ALL, "~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")

#####################################



### Analyze the variation in TSS
#####################################
###

library(dplyr)

### Combine TSS_Cor data and WindowSliding_TSS data
results_ALL=data.table::fread("~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
results_cor_merged=left_join(results_ALL, 
                             TSS_Cor %>% dplyr::select(rho, rho_pval, gene, analysis, celltypes) %>% 
                               subset(analysis!="all") %>% mutate(analysis=ifelse(analysis=="males","M","F")), 
                             by=c("gene","analysis","celltypes"))
correlated_subset=results_cor_merged %>% subset(rho_pval<0.05) %>% subset(!is.na(avg_log2FC)) %>%
  mutate(orientation=ifelse(avg_log2FC<0,"hypo","hyper")) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", function(x) mean(x, na.rm=T))

### Plot for all
correlated_subset_plot=correlated_subset %>% 
  subset(celltypes!="leukocyte") %>% # remove this one as it provides limited information
  mutate(analysis_orientation=paste0(analysis,"_",orientation)) %>%
  mutate(celltypes=ifelse(celltypes=="peripheral blood mononuclear cell","PBMC",celltypes))
# FvsM_plot=
  ggplot(correlated_subset_plot, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes) +
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/TSS.DMPs_Correlation_wAges_variation_FvsM.pdf", height=5, width=4.5)
plot(FvsM_plot)
dev.off()

### Plot for CD4T and CD8T only, and Mono as control
correlated_subset_plot=correlated_subset %>% 
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells",
                                 ifelse(celltypes=="CD14+ monocyte","Monocytes",celltypes)))) %>%
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  mutate(analysis_orientation=paste0(analysis,"_",orientation)) %>%
  mutate(panel="TSS")
# add the RNA variation below
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/PBMC_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(analysis!="all") %>% 
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation))
DEGs_Rough_selected=DEGs_Rough %>% 
  subset(analysis!="all" & p_val_adj<0.05) %>% 
  subset(celltypes %in% c("CD4T cells","CD8T cells","Monocytes")) %>%
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation)) %>%
  mutate(panel="RNA")

correlated_subset_plot=data.table::rbindlist(list(correlated_subset_plot, DEGs_Rough_selected))
plot_=
  ggplot(correlated_subset_plot %>% subset(celltypes!="Monocytes"), aes(x=age, y=avg_log2FC, color=analysis))+
  # facet_grid(panel~celltypes) +
  facet_wrap(~celltypes, nrow=1) +
  # geom_point(alpha=0.1, size=0.1) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  # coord_cartesian(ylim=c(-3.2,3.2)) +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=analysis_orientation, linetype=panel), method="loess", formula=y~log(x), se=F, linewidth=0.25) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL),
         linetype=guide_legend(title=NULL, override.aes=list(color="black")))

pdf("~/Project_PBMCage/Other_omics/Methy_plots/TSS.DMPs_Correlation_wAges_variation_and_RNAvariation_Rough_CD4T.CD8T_FvsM.pdf", height=2.5, width=4)
plot(plot_)
dev.off()

#####################################



### Analyze the variation of TSSs across ages in the key terms (12 superclusters)
#####################################
###

### Load the genes in the superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")

### Take the supercluster genes in DEGs data
DEGs_Rough=data.table::fread("~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(celltypes %in% c("CD4+ T cell","CD8+ T cell","CD14+ monocyte")) %>% 
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells","Monocytes")))
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      mutate(orientation=ifelse(avg_log2FC>0,"hyper.","hypo.")) %>%
                                      group_by(age, analysis, orientation, celltypes, superclusters) %>% 
                                      summarize_at("avg_log2FC", function(x) mean(x, na.rm=T)) %>% 
                                      .[!is.na(.$avg_log2FC),] %>%
                                      mutate(analysis_orientation=paste0(analysis,"_",orientation))) %>%
  data.table::rbindlist() %>%
  mutate(age=round(age))
DEGs_Rough_subset_final=DEGs_Rough_subset_percluster %>%
  # mutate(agecut.half=ifelse(age %in% c(19:25), "19~25", age)) %>%
  # mutate(agecut.half=ifelse(age %in% c(26:30), "26~30", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(31:35), "31~35", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(36:40), "36~40", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(41:45), "41~45", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(45:47), "45~47", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(48:50), "48~50", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(51:55), "51~55", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(56:60), "56~60", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(61:65), "61~65", agecut.half)) %>%
  # mutate(agecut.half=ifelse(age %in% c(66:68), "66~68", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(69:70), "69~70", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(71:75), "71~75", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(76:80), "76~80", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(81:85), "81~85", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(86:90), "86~90", agecut.half)) %>%
# mutate(agecut.half=ifelse(age %in% c(91:99), "91~99", agecut.half)) %>%
mutate(agecut.half=ifelse(age %in% c(19:30), "19~30", age)) %>%
  mutate(agecut.half=ifelse(age %in% c(31:47), "31~47", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(48:68), "48~68", agecut.half)) %>%
  mutate(agecut.half=ifelse(age %in% c(69:97), "69~97", agecut.half)) %>%
  dplyr::rename(agecut=agecut.half)
supecluster_=names(table(DEGs_Rough_subset_final$superclusters))
supecluster_reorder_idx=c(9,11,12,5, 6,3,2,10, 4,1,8,7)
supecluster_reorder=supecluster_[supecluster_reorder_idx]
DEGs_Rough_subset_final$superclusters=forcats::fct_relevel(DEGs_Rough_subset_final$superclusters, supecluster_reorder)
# test
for(i in 1:12) {
  print(paste0("==========",i,"=========="))
  print(
    ggpubr::compare_means(avg_log2FC_abs ~ celltypes,  
                          data=DEGs_Rough_subset_final %>% mutate(avg_log2FC_abs=abs(avg_log2FC)) %>%
                            subset(superclusters==supecluster_reorder[i]),
                          ref.group="CD4T cells",
                          method="wilcox.test")
  )
  print(
    ggpubr::compare_means(avg_log2FC_abs ~ celltypes,  
                          data=DEGs_Rough_subset_final %>% mutate(avg_log2FC_abs=abs(avg_log2FC)) %>%
                            subset(superclusters==supecluster_reorder[i]),
                          ref.group="CD8T cells",
                          method="wilcox.test")
  )
}
# focus on "leuk. devel." and "immune proc." which have been identified in RNA expr
FvsM_plot=
  ggplot(DEGs_Rough_subset_final, 
         aes(x=agecut, y=abs(avg_log2FC), color=celltypes))+
  facet_wrap(~superclusters) +
  coord_cartesian(ylim=c(0,2)) +
  ggrastr::geom_beeswarm_rast(shape=1, alpha=0.1, size=0.5) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(2,3,5)]) +
  labs(x=NULL, y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_line(linewidth=0.02, color="grey80"),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="top",
        legend.direction="horizontal") +
  geom_line(data=. %>% group_by(celltypes, superclusters, agecut) %>% summarize_at("avg_log2FC", function(x) mean(abs(x), na.rm=T)),
            aes(group=celltypes), linewidth=0.5, show.legend=F) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_Meyth.TSSs_variation_MapToSupercluster.pdf", height=6, width=8)
FvsM_plot
dev.off()

Immune_plot=
  ggplot(DEGs_Rough_subset_final %>% subset(superclusters %in% c("leuk. devel.","immune proc.")), 
         aes(x=agecut, y=abs(avg_log2FC), color=celltypes))+
  facet_wrap(~superclusters) +
  coord_cartesian(ylim=c(0,2)) +
  ggrastr::geom_beeswarm_rast(shape=1, alpha=0.1, size=0.5) +
  theme_classic() +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3")[c(2,3,5)]) +
  labs(x=NULL, y=expression(Avg.Log[2]~Fold-change)) +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        panel.grid.major=element_line(linewidth=0.02, color="grey80"),
        title=element_text(size=11),
        strip.text=element_text(size=10),
        legend.text=element_text(size=9),
        legend.position="top",
        legend.direction="horizontal") +
  geom_line(data=. %>% group_by(celltypes, superclusters, agecut) %>% summarize_at("avg_log2FC", function(x) mean(abs(x), na.rm=T)),
            aes(group=celltypes), linewidth=0.5, show.legend=F) +
    ggpubr::stat_compare_means(data=. %>% 
                                 group_by(analysis, orientation, celltypes, superclusters) %>% 
                                 summarize_at("avg_log2FC",.funs=function(x) {mean(abs(x))}) %>%
                                 mutate(agecut="31~47"),
                               method="anova",
                               aes(group=celltypes), label="p.format", label.y=1) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title=NULL))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_Meyth.TSSs_variation_MapToSupercluster_leukoAndImmune_Only.pdf", height=2.5, width=4.5)
Immune_plot
dev.off()

#####################################



### Variation in TSSs of 12 superclusters in CD8T
#####################################
###

### Load the genes in the superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")

### Take the supercluster genes in DEGs data
DEGs_Rough=data.table::fread("~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
DEGs_Rough=DEGs_Rough %>% 
  subset(celltypes %in% c("CD4+ T cell","CD8+ T cell","CD14+ monocyte")) %>% 
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells","Monocytes")))

### CD8T
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(celltypes=="CD8T cells") %>%
  subset(analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      mutate(abs.avg_log2FC=abs(avg_log2FC)) %>%
                                      summarize_at("abs.avg_log2FC", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
# remove outliers
outliers=lapply(DEGs_Rough_subset_percluster %>% split(.$superclusters),
                function(df) {boxplot.stats(df$abs.avg_log2FC)$out}
) %>% unlist() %>% as.data.frame() %>%
  tibble::rownames_to_column("supercluster") %>%
  rename(abs.avg_log2FC=".") %>%
  mutate(supercluster=gsub("[0-9]$","",supercluster)) %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(supercluster,"_",abs.avg_log2FC))
DEGs_Rough_subset_percluster=DEGs_Rough_subset_percluster %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(superclusters,"_",abs.avg_log2FC)) %>%
  subset(!(supercluster_abs.avgLOG2FC %in% outliers$supercluster_abs.avgLOG2FC))
# calculate sd for abs.avg_log2FC in each supercluster
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2FC", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2FC, na.rm=T)) %>%
  unlist()
DEGs_df1=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                    supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                    sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()

### CD4T
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(celltypes=="CD4T cells") %>%
  subset(analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      mutate(abs.avg_log2FC=abs(avg_log2FC)) %>%
                                      summarize_at("abs.avg_log2FC", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()
# remove outliers
outliers=lapply(DEGs_Rough_subset_percluster %>% split(.$superclusters),
                function(df) {boxplot.stats(df$abs.avg_log2FC)$out}
) %>% unlist() %>% as.data.frame() %>%
  tibble::rownames_to_column("supercluster") %>%
  rename(abs.avg_log2FC=".") %>%
  mutate(supercluster=gsub("[0-9]$","",supercluster)) %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(supercluster,"_",abs.avg_log2FC))
DEGs_Rough_subset_percluster=DEGs_Rough_subset_percluster %>%
  mutate(supercluster_abs.avgLOG2FC=paste0(superclusters,"_",abs.avg_log2FC)) %>%
  subset(!(supercluster_abs.avgLOG2FC %in% outliers$supercluster_abs.avgLOG2FC))
# calculate sd for abs.avg_log2FC in each supercluster
DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
  mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
  group_by(celltype_cl, age) %>%
  summarize_at("abs.avg_log2FC", function(x) mean(x, na.rm=T)) %>%
  ungroup() %>%
  split(.$celltype_cl) %>%
  lapply(., function(df) sd(df$abs.avg_log2FC, na.rm=T)) %>%
  unlist()
DEGs_df2=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
                    supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
                    sd=DEGs_summarize.at.age) %>%
  tibble::remove_rownames()

### Merge CD4T and CD8T sd
DEGs_df=rbind(DEGs_df1, DEGs_df2)
# reorder
# reorder_=names(table(DEGs_df$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
reorder_=DEGs_df %>% subset(celltype=="CD8T cells") %>% arrange(sd) %>% select(supercluster) %>% 
  tibble::deframe() %>%
  as.character()
DEGs_df$supercluster=forcats::fct_relevel(DEGs_df$supercluster, reorder_)

plot_=
  ggplot(DEGs_df, aes(x=supercluster, y=sd)) +
  facet_wrap(~celltype) +
  geom_bar(stat="identity", fill="grey80", color="transparent", width=0.75) +
  theme_classic() +
  labs(x="Supercluster", y="SD") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10)) +
  # coord_flip(ylim=c(0.2,0.4)) +
  coord_flip() +
  scale_x_discrete(guide=guide_axis(n.dodge=2)) +
  scale_y_continuous(n.breaks=4)

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_Meyth.TSSs_variation_superclusters_CD8T.CD4T.pdf", height=2, width=4)
plot(plot_)
dev.off()

#####################################



### Plot the variation in TSS across age in 12 superclusters in CD8T
#####################################
###

Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
DEGs_Rough=data.table::fread("~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
DEGs_Rough=DEGs_Rough %>% 
  subset(celltypes %in% c("CD4+ T cell","CD8+ T cell","CD14+ monocyte")) %>% 
  mutate(celltypes=ifelse(celltypes=="CD4+ T cell","CD4T cells",
                          ifelse(celltypes=="CD8+ T cell","CD8T cells","Monocytes")))

DEGs_Rough_subset=DEGs_Rough %>% subset(celltypes=="CD8T cells")
DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
                                    function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
                                      mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
                                      mutate(abs.avg_log2FC=abs(avg_log2FC)) %>%
                                      group_by(age, analysis, celltypes, superclusters) %>% 
                                      summarize_at("abs.avg_log2FC", function(x) mean(x, na.rm=T))) %>%
  data.table::rbindlist()

# rearrange the df
DEGs_df=DEGs_Rough_subset_percluster %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                  "cell activity",
                  ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                         "cell reaction")))
supecluster_=names(table(DEGs_df$superclusters))
supecluster_reorder=supecluster_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
DEGs_df$superclusters=forcats::fct_relevel(DEGs_df$superclusters, supecluster_reorder)
DEGs_df$main_group=forcats::fct_relevel(DEGs_df$main_group, 
                                        c("subcellular bioprocesses","cell activity","cell reaction"))

# plot_var=
  ggplot(DEGs_df, aes(x=age, y=abs.avg_log2FC, color=superclusters, linetype=main_group)) +
  geom_point(size=0.05, alpha=0.05) +
  geom_smooth(aes(group=superclusters), se=F, linewidth=0.5) +
  scale_linetype_manual(values=c("longdash","solid","dotted")) +
  coord_cartesian(ylim=c(0,0.6)) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=9, angle=0),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=9, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="right",
        legend.direction="vertical", 
        legend.spacing.y=unit(-0.5, "mm")
  ) +
  scale_x_continuous(limits=c(19,97)) +
  guides(color=guide_legend(title=NULL, ncol=1, byrow=F),
         linetype=guide_legend(title=NULL, override.aes=list(color="black"))) +
    geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1")

### Plot the AUC of variation across age in 12 superclusters in CD8T
# calculate the AUC
DEGs_df_AUC=DEGs_df %>% 
  group_by(age, celltypes, superclusters, main_group) %>%
  summarize_at("abs.avg_log2FC", mean) %>%
  split(.$superclusters) %>%
  lapply(., function(df)
  {integral_=with(df, 
                  integrate(approxfun(age, abs.avg_log2FC), 
                            lower=min(age), 
                            upper=max(age), subdivisions=300)) %>%
    .$value})
# arrange the df
DEGs_df_AUC_arranged=DEGs_df_AUC %>%
  unlist() %>%
  as.data.frame() %>%
  rename(AUC.abs.avg_log2FC=".") %>%
  tibble::rownames_to_column("superclusters") %>%
  mutate(main_group=
           ifelse(superclusters %in% c("energy metab.","cell division","autophagy","prog. death"),
                  "cell activity",
                  ifelse(superclusters %in% c("nucl. metab.","prot. metab.","ribosome syn.","cellular org."), "subcellular bioprocesses",
                         "cell reaction")))
# reorder
# reorder_=names(table(DEGs_df_AUC_arranged$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
# reorder_=DEGs_df_AUC_arranged %>% arrange(AUC.abs.avg_log2FC) %>% select(superclusters) %>%
#   tibble::deframe() %>%
#   as.character()
reorder_=c("ribosome syn.", "energy metab.", "prot. metab.", "autophagy",
           "anti-virus", "nucl. metab.", "prog. death", "cellular org.", 
           "cell division", "leuk. devel.", "cell response", "immune proc.")
DEGs_df_AUC_arranged$superclusters=forcats::fct_relevel(DEGs_df_AUC_arranged$superclusters, reorder_)

DEGs_df_AUC_arranged$main_group=forcats::fct_relevel(DEGs_df_AUC_arranged$main_group, 
                                                     c("subcellular bioprocesses","cell activity","cell reaction"))
plot_AUC=
  ggplot(DEGs_df_AUC_arranged, aes(x=superclusters, y=AUC.abs.avg_log2FC, fill=main_group)) +
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.75) +
  scale_fill_manual(values=c("white","grey80","grey40"),
                    labels=function(x) {gsub("subcellular ","subcellular\n",x)}) +
  labs(x=NULL, y=expression(atop(AUC~of,Avg.Abs.Log[2]~Fold-change))) +
  theme_classic() +
  theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
        axis.title.x=element_text(size=10),
        axis.text.y=element_text(size=9),
        axis.title.y=element_text(size=10),
        strip.text=element_text(size=9, color="black"),
        strip.background=element_rect(fill="white", color="black"),
        legend.position="top",
        legend.direction="horizontal", 
        legend.spacing.y=unit(-0.5, "mm")
  ) +
  guides(fill=guide_legend(title=NULL, override.aes=list(size=1)))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/WindowSliding_Meyth.TSSs_variation_AUC_superclusters.pdf", height=2.5, width=4)
plot(plot_AUC)
dev.off()

#####################################

