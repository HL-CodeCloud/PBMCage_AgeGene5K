library(dplyr)
library(Seurat)

### Window-sliding to determine the agecut
#####################################
###
###

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
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
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i], slot="data")
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
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")

#####################################



### Window-sliding to determine the agecut in females
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
THEObj=subset(THEObj, Sex=="Female")
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
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly.rds")

#####################################



### Window-sliding to determine the agecut in males
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
THEObj=subset(THEObj, Sex=="Male")
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
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly.rds")

#####################################



### Clean the results
#####################################
###

library(dplyr)

###
The_Whole_List=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")
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
saveRDS(list(Count_Df, Count_Df_r, Count_Df_i, Count_Df_d), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_counts.rds")

#####################################



### Plot the Window sliding results
#####################################
###

library(ggplot2)

###
The_Whole_List=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_counts.rds")
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
  ggsci::scale_color_d3() +
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
    ggsci::scale_color_d3() +
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
    ggsci::scale_color_d3() +
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
    ggsci::scale_color_d3() +
    theme_classic() +
    theme(axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          axis.title.x=element_text(size=9),
          axis.title.y=element_text(size=9),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9)) +
    guides(color=guide_legend(title="Celltype", override.aes=list(size=2, alpha=1)))

# pdf("   ", height=9, width=8)
# plot(plot_all)
# plot(plot_rough)
# plot(plot_inter)
# plot(plot_detailed)
# dev.off()

#####################################



### Analyze the variation of DEGs across ages
#####################################
###

library(dplyr)

### Load the results
The_Whole_List=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs.rds")
The_Whole_List_F=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly.rds")
The_Whole_List_M=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly.rds")

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
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Donorid.csv.gz", sep="\t")

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
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")

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
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Inter.csv.gz", sep="\t")

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
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Detailed.csv.gz", sep="\t")

### Analyze the variation at the donor_id level
DEGs_Donorid=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Donorid.csv.gz", sep="\t")
DEGs_Donorid_bothsex=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all" & !(avg_log2FC %in% c(Inf, -Inf)))
DEGs_Donorid_bothsex=DEGs_Donorid_bothsex %>% group_by(age, orientation) %>% summarize_at("avg_log2FC", mean)
# mixed_sex_plot=
  ggplot(DEGs_Donorid_bothsex, aes(x=age, y=avg_log2FC, color=orientation))+
    geom_hline(yintercept=0) +
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
  geom_smooth(aes(group=orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  geom_line(linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation.pdf", height=4, width=2.5)
plot(mixed_sex_plot)
dev.off()

### Analyze the variation at the donor_id level, comparing F vs. M
DEGs_Donorid_subset=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
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
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_FvsM.pdf", height=4, width=2.5)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the rough level
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_bothsex=DEGs_Rough %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% group_by(age, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% subset(!celltypes %in% c("DCs","Other cells")) # have checked, these two do not have enough data points, so excluded
mixed_sex_plot=
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

### Analyze the variation at the rough level, comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% 
  subset(analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf) %>% 
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation))
DEGs_Rough_selected=DEGs_Rough %>% 
  subset(analysis!="all" & p_val_adj<0.05 & avg_log2FC<Inf & avg_log2FC>-Inf) %>% 
  # show only lymphocytes as in the previous text we focused on their shared alterations
  subset(celltypes %in% c("B cells","CD4T cells","CD8T cells","NK cells","OtherT")) %>%
  mutate(celltypes=gsub("\\..*","",celltypes)) %>%
  group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean) %>% 
  mutate(analysis_orientation=paste0(analysis,"_",orientation))

# FvsM_plot_selected=
  ggplot(DEGs_Rough_selected, aes(x=age, y=avg_log2FC))+
  facet_wrap(~celltypes, ncol=1) +
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
  # geom_smooth(aes(group=orientation), method="loess", formula=y~log(x), se=F, linewidth=0.25) +
  geom_line(aes(group=orientation), linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

FvsM_plot=
  ggplot(DEGs_Rough_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=4) +
  geom_hline(yintercept=c(0.45, -0.45), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y=expression(Avg.Log[2]~Fold-change), title="Terekhova et al. (2023) dataset") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        legend.position="right",
        legend.direction="vertical") +
  scale_x_continuous(breaks=c(30,60,90)) +
  # geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  # geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_Rough.All_FvsM.pdf", height=3.5, width=6)
plot(FvsM_plot)
dev.off()

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_Rough_FvsM.pdf", height=5, width=8)
plot(FvsM_plot_selected)
dev.off()

### Analyze the variation at the inter level, comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Inter.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset=DEGs_Rough_subset %>% group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
# have checked, B.pbpc, B.trans, CD4T.Treg, DC.cDC2, DC.pDC, NK.prolif, Other.progenitor, and Other.MAIT, otherT.NKT have too few data samples,
#... not really meaningful, therefore removed
DEGs_Rough_subset=DEGs_Rough_subset %>% subset(!celltypes %in% c("B.pbpc","B.trans","CD4T.Treg","CD8T.Trm","DC.cDC2","DC.pDC","NK.prolif","Other.progenitor","OtherT.MAIT","OtherT.NKT","OtherT.dnT"))
# take only CD4 and CD8T that are interesting
DEGs_Rough_subset=DEGs_Rough_subset %>% subset(grepl("CD4|CD8",celltypes))
DEGs_Rough_subset=DEGs_Rough_subset %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))

FvsM_plot=
  ggplot(DEGs_Rough_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=3) +
  geom_hline(yintercept=c(5, -5), linetype="dashed", linewidth=0.25) +
  geom_point(alpha=0.25, size=0.75) +
  coord_cartesian(ylim=c(-100,100)) +
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

pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_Inter.CD4CD8T_FvsM.pdf", height=3, width=4.5)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the detailed level, comparing F vs. M
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Detailed.csv.gz", sep="\t")
DEGs_Rough_subset=DEGs_Rough %>% subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
DEGs_Rough_subset=DEGs_Rough_subset %>% group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", mean)
# take only CD4 and CD8T that are interesting
DEGs_Rough_subset=DEGs_Rough_subset %>% subset(grepl("CD4|CD8",celltypes))
DEGs_Rough_subset=DEGs_Rough_subset %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))

FvsM_plot=
  ggplot(DEGs_Rough_subset, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=3) +
  coord_cartesian(ylim=c(-100,100)) +
  geom_point(alpha=0.25) +
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
  geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

#... turns out to be not good because some celltypes reduced numbers significantly thereby influencing gene expression

#####################################
  


### Analyze the variation of DEGs across ages in the key terms (12 superclusters)
#####################################
###

### Load the genes in the superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")

### Take the supercluster genes in DEGs data
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
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

# FvsM_plot_selected=
ggplot(DEGs_Rough_selected, aes(x=age, y=avg_log2FC, color=analysis))+
  facet_wrap(~celltypes, ncol=5) +
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
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

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
  coord_cartesian(ylim=c(0,50)) +
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
            aes(x="31~47", y=45, label="p<0.05"), color="black", size=3.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1, shape=19), title=NULL))
  
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_MapToSupercluster.pdf", height=6, width=8)
FvsM_plot
dev.off()
  
#####################################



# ### Normalize the variation by cell counts (sd.RMS="Integral of sd across N" / N)
# # ... not consistent with PBMCage dataset, for undetermined reason
# #####################################
# ###
# 
# Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
# DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough.csv.gz", sep="\t")
# DEGs_Rough_subset=DEGs_Rough %>% 
#   subset(p_val_adj<0.05 & analysis!="all" & avg_log2FC<Inf & avg_log2FC>-Inf)
# DEGs_Rough_subset_percluster=lapply(1:length(Genes_in_SuperClusters), 
#                                     function(idx) DEGs_Rough_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]) %>%
#                                       mutate(superclusters=gsub("\\|.*","",names(Genes_in_SuperClusters)[idx])) %>%
#                                       group_by(age, analysis, celltypes, superclusters) %>% 
#                                       summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T))) %>%
#   data.table::rbindlist()
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
# # calculate sd for abs.avg_log2FC in each supercluster
# DEGs_summarize.at.age=DEGs_Rough_subset_percluster %>%
#   mutate(celltype_cl=paste0(celltypes,"__",superclusters)) %>%
#   group_by(celltype_cl, age) %>%
#   summarize_at("abs.avg_log2F", function(x) mean(x, na.rm=T)) %>%
#   ungroup() %>%
#   split(.$celltype_cl) %>%
#   lapply(., function(df) sd(df$abs.avg_log2F, na.rm=T)) %>%
#   unlist()
# DEGs_df=data.frame(celltype=gsub("__.*","",names(DEGs_summarize.at.age)),
#                    supercluster=gsub(".*__","",names(DEGs_summarize.at.age)),
#                    sd=DEGs_summarize.at.age) %>%
#   tibble::remove_rownames()
# # add the cell count info
# meta_info_combined=read.delim("~/Project_PBMCage/Results/Immunity_results/Cell_frequency_values.txt", sep="\t")
# cell_count=meta_info_combined %>% dplyr::select(donor_id, Annot.rough, N_per_rough) %>%
#   subset(!duplicated(.)) %>%
#   group_by(Annot.rough) %>% summarize_at("N_per_rough", function(x) sum(x, na.rm=T)) %>%
#   rename(celltype=Annot.rough) %>%
#   right_join(., DEGs_df, by=c("celltype"))
# # normalize cell count by sd.RMS
# sd.div.N_persupercluster=cell_count %>% 
#   subset(!is.na(sd)) %>%
#   split(.$supercluster) %>%
#   lapply(., 
#          function(df) {
#            integral_=with(df, 
#                           integrate(approxfun(N_per_rough, sd), 
#                                     lower=min(N_per_rough), 
#                                     upper=max(N_per_rough))) %>%
#              .$value
#            divident_=max(df$N_per_rough)-min(df$N_per_rough)
#            integral_/divident_
#          }) %>%
#   unlist() %>% as.data.frame() %>%
#   dplyr::rename(RMS=".") %>%
#   tibble::rownames_to_column("supercluster")
# 
# # reorder_=names(table(sd.div.N_persupercluster$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
# reorder_=sd.div.N_persupercluster %>% arrange(RMS) %>% select(supercluster) %>% 
#   tibble::deframe() %>%
#   as.character()
# sd.div.N_persupercluster$supercluster=forcats::fct_relevel(sd.div.N_persupercluster$supercluster, reorder_)
# 
# # plot_=
#   ggplot(sd.div.N_persupercluster, aes(x=supercluster, y=RMS)) +
#   geom_bar(stat="identity", fill="grey80", color="transparent", width=0.75) +
#   theme_classic() +
#   labs(x="Supercluster", y="Mean of SD") +
#   theme(axis.text.y=element_text(size=9),
#         axis.text.x=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10)) +
#   # coord_flip(ylim=c(0.2,0.6)) +
#   coord_flip() +
#   scale_x_discrete(guide = guide_axis(n.dodge=2))
# 
# #####################################



# ### Calculate the avg_log2FC with window-sliding on the TSS methylation data
# #####################################
# ###
# 
# library(dplyr)
# library(Seurat)
# 
# ### Load the data
# beta_meta_merged_ALL=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets.csv.gz", sep="\t")
# all_celltypes=unique(beta_meta_merged_ALL$celltypes)
# 
# ### Window-sliding analysis
# results_per_age_ALL=list()
# for (i in 1:length(all_celltypes)) {
#   the_celltype=all_celltypes[i]
#   results_per_age_both=list()
#   
#   # for females
#   beta_meta_subset=beta_meta_merged_ALL %>% subset(celltypes==the_celltype & sex=="F")
#   beta_meta_subset_wide=beta_meta_subset %>% select(gene, age, beta) %>% group_by(gene, age) %>%
#     summarize_at("beta", mean) %>%
#     tidyr::pivot_wider(names_from="age", values_from="beta")
#   results_per_age=lapply(3:ncol(beta_meta_subset_wide), function(x) {
#     df=data.frame(beta_meta_subset_wide[,x-1], beta_meta_subset_wide[,x]); colnames(df)=c("former","latter"); df %>% mutate(avg_log2FC=log2(latter+1e-5)-log2(former+1e-5)) %>% select(avg_log2FC) %>% unlist() %>% unname()
#   })
#   names(results_per_age)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
#   results_per_age_df=as.data.frame(results_per_age)
#   rownames(results_per_age_df)=beta_meta_subset_wide$gene
#   colnames(results_per_age_df)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
#   results_per_age_df=results_per_age_df %>% 
#     mutate(gene=rownames(.)) %>%
#     tidyr::pivot_longer(cols=colnames(results_per_age_df), names_to="age", values_to="avg_log2FC") %>%
#     mutate(celltypes=the_celltype, analysis="F")
#   results_per_age_both[[1]]=results_per_age_df
#   
#   # for males
#   beta_meta_subset=beta_meta_merged_ALL %>% subset(celltypes==the_celltype & sex=="M")
#   beta_meta_subset_wide=beta_meta_subset %>% select(gene, age, beta) %>% group_by(gene, age) %>%
#     summarize_at("beta", mean) %>%
#     tidyr::pivot_wider(names_from="age", values_from="beta")
#   results_per_age=lapply(3:ncol(beta_meta_subset_wide), function(x) {
#     df=data.frame(beta_meta_subset_wide[,x-1], beta_meta_subset_wide[,x]); colnames(df)=c("former","latter"); df %>% mutate(avg_log2FC=log2(latter+1e-5)-log2(former+1e-5)) %>% select(avg_log2FC) %>% unlist() %>% unname()
#   })
#   names(results_per_age)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
#   results_per_age_df=as.data.frame(results_per_age)
#   rownames(results_per_age_df)=beta_meta_subset_wide$gene
#   colnames(results_per_age_df)=colnames(beta_meta_subset_wide)[3:ncol(beta_meta_subset_wide)]
#   results_per_age_df=results_per_age_df %>% 
#     mutate(gene=rownames(.)) %>%
#     tidyr::pivot_longer(cols=colnames(results_per_age_df), names_to="age", values_to="avg_log2FC") %>%
#     mutate(celltypes=the_celltype, analysis="M")
#   results_per_age_both[[2]]=results_per_age_df
#   
#   # merge the sex
#   results_per_age_ALL[[i]]=data.table::rbindlist(results_per_age_both)
# }
# 
# results_ALL=data.table::rbindlist(results_per_age_ALL)
# data.table::fwrite(results_ALL, "~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
# 
# #####################################
# 
# 
# 
# ### Analyze the variation in TSS
# #####################################
# ###
# 
# ### Combine TSS_Cor data and WindowSliding_TSS data
# results_ALL=data.table::fread("~/Project_PBMCage/Other_omics/Methy_results/WindowSliding_TSS.csv.gz", sep="\t")
# TSS_Cor=data.table::fread("~/Project_PBMCage/Other_omics/MethyTSS_subsets/AgeMeth_data_TSS.subsets_ageCorrelationResult.csv.gz", sep="\t")
# results_cor_merged=left_join(results_ALL, 
#                              TSS_Cor %>% select(rho, rho_pval, gene, analysis, celltypes) %>% subset(analysis!="all") %>% mutate(analysis=ifelse(analysis=="males","M","F")), 
#                              by=c("gene","analysis","celltypes"))
# correlated_subset=results_cor_merged %>% subset(rho_pval<0.05) %>% subset(!is.na(avg_log2FC)) %>%
#   mutate(orientation=ifelse(avg_log2FC<0,"hypo","hyper")) %>%
#   group_by(age, analysis, orientation, celltypes) %>% summarize_at("avg_log2FC", function(x) mean(x, na.rm=T))
# 
# ### Plot
# correlated_subset=correlated_subset %>% 
#   subset(celltypes!="leukocyte") %>% # remove this one as it provides limited information
#   mutate(analysis_orientation=paste0(analysis,"_",orientation)) %>%
#   mutate(celltypes=ifelse(celltypes=="peripheral blood mononuclear cell","PBMC",celltypes))
# FvsM_plot=
#   ggplot(correlated_subset, aes(x=age, y=avg_log2FC, color=analysis))+
#   facet_wrap(~celltypes) +
#   geom_point(alpha=0.25) +
#   theme_classic() +
#   scale_color_manual(values=c("coral3","deepskyblue3")) +
#   labs(x="age", y=expression(Avg.Log[2]~Fold-change)) +
#   theme(axis.text.y=element_text(size=9),
#         axis.text.x=element_text(size=9),
#         axis.title.x=element_text(size=10),
#         axis.title.y=element_text(size=10),
#         legend.position="top",
#         legend.direction="horizontal") +
#   # geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
#   geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
#   guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
# 
# pdf("~/Project_PBMCage/Other_omics/Methy_plots/TSS.DMPs_Correlation_wAges_variation_FvsM.pdf", height=5, width=4.5)
# plot(FvsM_plot)
# dev.off()
# 
# #####################################
# 
# 
# 











###############################################################################################################
##################################### REANALYZE THE WINDOW-SLIDING WITH SCALE.DATA #####################################
###############################################################################################################

library(dplyr)
library(Seurat)

### Window-sliding to determine the agecut
#####################################
###
###

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
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
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_wScaleData.rds")

#####################################



### Window-sliding to determine the agecut in females
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
THEObj=subset(THEObj, Sex=="Female")
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
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly_wScaleData.rds")

#####################################



### Window-sliding to determine the agecut in males
#####################################
###

library(dplyr)
library(Seurat)

# Load
THEObj=readRDS("~/Project_PBMCage/Immunity/Immunity_RawCountObj_OneSamplePerDonor.rds")
THEObj[[]]=THEObj[[]] %>% dplyr::rename(age=Age)
THEObj=subset(THEObj, Sex=="Male")
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
    markers=FindMarkers(THEObj, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
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
      markers=FindMarkers(subset_, ident.1=age_[i+1], ident.2=age_[i], slot="scale.data")
      DEGs_detailed=c(DEGs_detailed, list(markers))
      names_=c(names_, paste0(celltypes_d[j], "_", age_[i+1]))
    }, error=function(msg) {print("Skip.")})
  }
}
names(DEGs_detailed)=names_

# save
saveRDS(list(DEGs_total, DEGs_rough, DEGs_inter, DEGs_detailed), "~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly_wScaleData.rds")

#####################################



### Analyze the variation of DEGs across ages
#####################################
###

library(dplyr)

### Load the results
The_Whole_List=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_wScaleData.rds")
The_Whole_List_F=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_FemalesOnly_wScaleData.rds")
The_Whole_List_M=readRDS("~/Project_PBMCage/Immunity/Tempt_RDS/WindowSliding_RNAseq_DEGs_MalesOnly_wScaleData.rds")

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
  select(gene, age, avg_diff, p_val_adj, celltype.level, celltypes, analysis) %>%
  mutate(orientation=ifelse(avg_diff>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_diff=abs(avg_diff))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Donorid_wScaleData.csv.gz", sep="\t")

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
  select(gene, age, avg_diff, p_val_adj, celltype.level, celltypes, analysis) %>%
  mutate(orientation=ifelse(avg_diff>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_diff=abs(avg_diff))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough_wScaleData.csv.gz", sep="\t")

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
  select(gene, age, avg_diff, p_val_adj, celltype.level, celltypes, analysis) %>%
  mutate(orientation=ifelse(avg_diff>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_diff=abs(avg_diff))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Inter_wScaleData.csv.gz", sep="\t")

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
  select(gene, age, avg_diff, p_val_adj, celltype.level, celltypes, analysis) %>%
  mutate(orientation=ifelse(avg_diff>0, "higher.than.avg", "lower.than.avg"),
         abs.avg_diff=abs(avg_diff))
data.table::fwrite(DEGs_total_clean, "~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Detailed_wScaleData.csv.gz", sep="\t")

### Analyze the variation at the donor_id level
DEGs_Donorid=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Donorid_wScaleData.csv.gz", sep="\t")
DEGs_Donorid_bothsex=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all" )
DEGs_Donorid_bothsex=DEGs_Donorid_bothsex %>% group_by(age, orientation) %>% summarize_at("avg_diff", mean)
mixed_sex_plot=
ggplot(DEGs_Donorid_bothsex, aes(x=age, y=avg_diff, color=orientation))+
  # geom_hline(yintercept=0) +
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c("brown3","dodgerblue3"), labels=c("Increasing", "Decreasing")) +
  labs(x="age", y="Avg.Diff") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  geom_line(linewidth=0.5) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_wScaleData.pdf", height=4, width=2.5)
plot(mixed_sex_plot)
dev.off()

### Analyze the variation at the donor_id level, comparing F vs. M
DEGs_Donorid_subset=DEGs_Donorid %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Donorid_subset=DEGs_Donorid_subset %>% group_by(age, analysis, orientation) %>% summarize_at("avg_diff", mean)
DEGs_Donorid_subset=DEGs_Donorid_subset %>% mutate(analysis_orientation=paste0(analysis,"_",orientation))
FvsM_plot=
  ggplot(DEGs_Donorid_subset, aes(x=age, y=avg_diff, color=analysis))+
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c(scales::hue_pal()(2)[1],scales::hue_pal()(2)[2])) +
  labs(x="age", y="Avg.Diff") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top") +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  # geom_smooth(aes(group=analysis_orientation), method="loess", formula=y~log2(x), se=F, linewidth=0.5) +
  geom_line(aes(linetype=orientation), show.legend=F) +
  scale_linetype_manual(values=c("solid","solid")) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))
pdf("~/Project_PBMCage/Plots/Immunity/Immunity_WindowSliding_RNA_DEGs_variation_FvsM_wScaleData.pdf", height=4, width=2.5)
plot(FvsM_plot)
dev.off()

### Analyze the variation at the rough level
DEGs_Rough=data.table::fread("~/Project_PBMCage/Results/Immunity_results/WindowSliding_DEGs_Rough_wScaleData.csv.gz", sep="\t")
DEGs_Rough_bothsex=DEGs_Rough %>% subset(p_val_adj<0.05 & analysis!="all")
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% group_by(age, orientation, celltypes) %>% summarize_at("avg_diff", mean)
DEGs_Rough_bothsex=DEGs_Rough_bothsex %>% subset(!celltypes %in% c("DCs","Other cells")) # have checked, these two do not have enough data points, so excluded
# mixed_sex_plot=
  ggplot(DEGs_Rough_bothsex, aes(x=age, y=avg_diff, color=orientation))+
  facet_wrap(~celltypes) +
  geom_point(alpha=0.25) +
  theme_classic() +
  scale_color_manual(values=c("grey60","grey30"), labels=c("Higher than avg.", "Lower than avg.")) +
  labs(x="age", y="Avg.Diff") +
  theme(axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        legend.position="top",
        legend.direction="horizontal") +
  scale_x_continuous(breaks=c(30,60,90)) +
  geom_vline(xintercept=c(30,47,68), linetype="dashed", linewidth=0.2, color="coral1") +
  geom_smooth(aes(group=orientation), method="loess", formula=y~log(x), se=F, linewidth=0.5) +
  geom_line(aes(linetype=orientation)) +
  scale_linetype_manual(values=c("solid","solid")) +
  guides(color=guide_legend(override.aes=list(size=2, alpha=1), title=NULL))

