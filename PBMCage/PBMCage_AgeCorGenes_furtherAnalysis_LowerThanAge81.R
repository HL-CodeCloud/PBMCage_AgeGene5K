


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



### Enrich age-correlated genes across all celltypes (common ones)
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough_LowerThanAge81.csv.gz", sep="\t")
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0 & cor_pval<0.05 & cor<0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough_LowerThanAge81.rds")

upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0 & cor_pval<0.05 & cor>0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Rough_LowerThanAge81.rds")

#####################################



### Enrich age-correlated genes after exclusion of ribosome-related ones
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough_LowerThanAge81.csv.gz", sep="\t")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough_LowerThanAge81.rds")$ego
ego_results=ego_cluster@compareClusterResult %>% 
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]

### Exclude them and do EnrichGO
# downregulated ones
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0 & cor_pval<0.05 & cor<0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check
# ... turns out to be a part of metabolic processes shared by some celltypes while another part shared by other cell types...

# downregulated ones
upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0 & cor_pval<0.05 & cor>0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check
# ... turns out that cytoplasmic translation/RNA splicing/etc. are common, while leukocyte-mediated immunity is shared by some cell types...

# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough_LowerThanAge81.rds")

#####################################



### Get all the genes annotated to the key terms enriched based on the GO gaf file
#####################################
###

library(dplyr)

### Map age-corr. genes in each celltypes to the SuperCluster genes
# get the downreg./upreg. age-correlated genes in each celltype
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Rough_LowerThanAge81.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
Cor_analysis_DF_subset=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  mutate(orientation=ifelse(rho<0,"downreg.","upreg.")) %>%
  arrange(desc(abs(rho))) %>%
  split(.$orientation) %>%
  .[c("downreg.","upreg.")] %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% .[celltypes_] %>%
           lapply(., function(dff) dff %>% dplyr::select(gene, rho, rho_pval, celltypes, orientation)))
# map to the genes in superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=list()
for (o in 1:length(Cor_analysis_DF_subset)) {
  orientation_subset=Cor_analysis_DF_subset[[o]]
  exist_genes_per_celltype=list()
  for (i in 1:length(orientation_subset)) {
    celltype_subset=orientation_subset[[i]]
    exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters), function(idx) celltype_subset %>% subset(gene %in% Genes_in_SuperClusters[[idx]]))
    names(exist_genes_per_supercluster)=names(Genes_in_SuperClusters)
    exist_genes_per_celltype[[i]]=exist_genes_per_supercluster
  }
  names(exist_genes_per_celltype)=names(orientation_subset)
  exist_genes_per_orientation[[o]]=exist_genes_per_celltype
}
names(exist_genes_per_orientation)=c("downreg.","upreg.")

saveRDS(exist_genes_per_orientation, "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_Genes_in_GOSuperClusters_LowerThanAge81.rds")

#####################################



### Enrichment term scoring of the 12 superclusters in each celltype at the rough level
#####################################
###

### Load the GO terms for each supercluster
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
GOs_both_df=
  as.data.frame(unlist(GOs_both)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_both)') %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_pos=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[1]]
names(GOs_pos)=gsub("\\|.*","",names(GOs_pos))
GOs_pos_df=
  as.data.frame(unlist(GOs_pos)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_pos)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_neg=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[2]]
names(GOs_neg)=gsub("\\|.*","",names(GOs_neg))
GOs_neg_df=
  as.data.frame(unlist(GOs_neg)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_neg)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_all=data.table::rbindlist(list(GOs_pos_df, GOs_neg_df, GOs_both_df))

### Load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Rough_LowerThanAge81.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Rough_LowerThanAge81.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough_LowerThanAge81.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Rough_LowerThanAge81.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)

### Map the GO IDs to the 12 supercluster
merged_df=full_join(ego_all, GOs_all, by=c("GO_ID")) %>% 
  subset(!is.na(Cluster) & !is.na(supercluster)) %>%
  mutate(p.adjust_wSign=ifelse(orientation=="upreg",p.adjust,-p.adjust)) %>%
  group_by(Cluster, supercluster, orientation) %>%
  slice_max(p.adjust, n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(Cluster, supercluster) %>%
  summarize_at("p.adjust_wSign", sum) %>%
  left_join(data.frame(Cluster='B cells', supercluster=NA, p.adjust_wSign=NA)) %>%
  tidyr::pivot_wider(names_from="supercluster", id_cols="Cluster", values_from="p.adjust_wSign") %>%
  tibble::column_to_rownames("Cluster") %>%
  replace(is.na(.), 0)

### Plot
library(ComplexHeatmap)
mtx_=scale(merged_df, scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]

range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-3, 0, 3), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(mtx_,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Scaled enrichment index",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Yazar et al. (2022) dataset (age<=80)", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(6, "mm"),
          height=nrow(mtx_)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_LowerThanAge81.pdf", height=3, width=4.5)
draw(ht_final)
dev.off()

#####################################






###### ====== Use JapanSC dataset to prove the ribosomal changes in lymphocytes ======


### Get the DEGs of SC vs. CT
#####################################
###

library(Seurat)

### Load the object
THEObj_sc=readRDS("~/Project_PBMCage/Japan supercentenarians/JapanSC_annotated.rds")

### DEGs of SC vs. CT
DEG_lists=list()
celltypes=names(table(THEObj_sc$Annot.rough))
for (i in 1:length(celltypes)) {
  subset_=subset(THEObj_sc, Annot.rough==celltypes[i])
  Idents(subset_)="agecut"
  DEG_lists[[i]]=FindMarkers(subset_, ident.1="SC", ident.2="CT")
}
names(DEG_lists)=celltypes
saveRDS(DEG_lists, "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")

### DEGs of over vs. below 80
DEG_lists=list()
celltypes=names(table(THEObj_sc$Annot.rough))
for (i in 1:length(celltypes)) {
  subset_=subset(THEObj_sc, Annot.rough==celltypes[i])
  subset_[[]]=subset_[[]] %>% mutate(agecut_80=ifelse(age %in% c("50","60","70"),"below80","over80"))
  Idents(subset_)="agecut_80"
  DEG_lists[[i]]=FindMarkers(subset_, ident.1="over80", ident.2="below80")
}
names(DEG_lists)=celltypes
saveRDS(DEG_lists, "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGsAt80_perAnnot.rough.rds")

#####################################



### Enrich DEGs across all celltypes (common ones)
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the DEGs of SC vs. CT and run enrichment across celltypes
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")
downreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC<0 & p_val_adj<0.05) %>% rownames(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedDownregGenes.EnrichResults_SCvsCT.rds")
upreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC>0 & p_val_adj<0.05) %>% rownames(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_SCvsCT.rds")

### Load the DEGs of over80 vs. below80 and run enrichment across celltypes
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGsAt80_perAnnot.rough.rds")
downreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC<0 & p_val_adj<0.05) %>% rownames(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedDownregGenes.EnrichResults_At80.rds")
upreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC>0 & p_val_adj<0.05) %>% rownames(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_At80.rds")

#####################################



### Enrich DEGs after exclusion of ribosome-related ones
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the DEGs of SC vs. CT and run enrichment across celltypes
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_SCvsCT.rds")$ego
ego_results=ego_cluster@compareClusterResult %>% 
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]
# exclude them and do EnrichGO
downreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC<0 & p_val_adj<0.05 & !(rownames(.) %in% rb_related_genes)) %>% rownames(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check
upreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC>0 & p_val_adj<0.05 & !(rownames(.) %in% rb_related_genes)) %>% rownames(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check
# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_SCvsCT.rds")

### Load the DEGs of over80 vs. below80 and run enrichment across celltypes
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGsAt80_perAnnot.rough.rds")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_At80.rds")$ego
ego_results=ego_cluster@compareClusterResult %>% 
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]
# exclude them and do EnrichGO
downreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC<0 & p_val_adj<0.05 & !(rownames(.) %in% rb_related_genes)) %>% rownames(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check
upreg_genes=lapply(DEG_lists, function(df) df %>% subset(avg_log2FC>0 & p_val_adj<0.05 & !(rownames(.) %in% rb_related_genes)) %>% rownames(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check
# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_At80.rds")

#####################################



### Get all the genes annotated to the key terms enriched based on the GO gaf file
#####################################
###

library(dplyr)

### Map DEGs of SC vs. CT in each celltypes to the SuperCluster genes
# get the downreg./upreg. age-correlated genes in each celltype
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")
celltypes_=names(DEG_lists)
orientation_subset_down=orientation_subset_up=list()
for (i in 1:length(celltypes_)) {
  tempt_=
    DEG_lists[[i]] %>% subset(p_val_adj<0.05) %>%
    mutate(orientation=ifelse(avg_log2FC<0,"downreg.","upreg.")) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    split(.$orientation) %>%
    .[c("downreg.","upreg.")]
  orientation_subset_down=c(orientation_subset_down, list(tempt_[["downreg."]]))
  orientation_subset_up=c(orientation_subset_up, list(tempt_[["downreg."]]))
}
names(orientation_subset_down)=celltypes_
names(orientation_subset_up)=celltypes_
Cor_analysis_DF_subset=list(orientation_subset_down, orientation_subset_up)
# map to the genes in superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=list()
for (o in 1:length(Cor_analysis_DF_subset)) {
  orientation_subset=Cor_analysis_DF_subset[[o]]
  exist_genes_per_celltype=list()
  for (i in 1:length(orientation_subset)) {
    celltype_subset=orientation_subset[[i]]
    exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters), function(idx) celltype_subset %>% subset(rownames(.) %in% Genes_in_SuperClusters[[idx]]))
    names(exist_genes_per_supercluster)=names(Genes_in_SuperClusters)
    exist_genes_per_celltype[[i]]=exist_genes_per_supercluster
  }
  names(exist_genes_per_celltype)=names(orientation_subset)
  exist_genes_per_orientation[[o]]=exist_genes_per_celltype
}
names(exist_genes_per_orientation)=c("downreg.","upreg.")
saveRDS(exist_genes_per_orientation, "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_Genes_in_GOSuperClusters_SCvsCT.rds")

### Map DEGs of over80 vs. below80 in each celltypes to the SuperCluster genes
# get the downreg./upreg. age-correlated genes in each celltype
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGsAt80_perAnnot.rough.rds")
celltypes_=names(DEG_lists)
orientation_subset_down=orientation_subset_up=list()
for (i in 1:length(celltypes_)) {
  tempt_=
    DEG_lists[[i]] %>% subset(p_val_adj<0.05) %>%
    mutate(orientation=ifelse(avg_log2FC<0,"downreg.","upreg.")) %>%
    arrange(desc(abs(avg_log2FC))) %>%
    split(.$orientation) %>%
    .[c("downreg.","upreg.")]
  orientation_subset_down=c(orientation_subset_down, list(tempt_[["downreg."]]))
  orientation_subset_up=c(orientation_subset_up, list(tempt_[["downreg."]]))
}
names(orientation_subset_down)=celltypes_
names(orientation_subset_up)=celltypes_
Cor_analysis_DF_subset=list(orientation_subset_down, orientation_subset_up)
# map to the genes in superclusters
Genes_in_SuperClusters=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
exist_genes_per_orientation=list()
for (o in 1:length(Cor_analysis_DF_subset)) {
  orientation_subset=Cor_analysis_DF_subset[[o]]
  exist_genes_per_celltype=list()
  for (i in 1:length(orientation_subset)) {
    celltype_subset=orientation_subset[[i]]
    exist_genes_per_supercluster=lapply(1:length(Genes_in_SuperClusters), function(idx) celltype_subset %>% subset(rownames(.) %in% Genes_in_SuperClusters[[idx]]))
    names(exist_genes_per_supercluster)=names(Genes_in_SuperClusters)
    exist_genes_per_celltype[[i]]=exist_genes_per_supercluster
  }
  names(exist_genes_per_celltype)=names(orientation_subset)
  exist_genes_per_orientation[[o]]=exist_genes_per_celltype
}
names(exist_genes_per_orientation)=c("downreg.","upreg.")
saveRDS(exist_genes_per_orientation, "~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_Genes_in_GOSuperClusters_At80.rds")

#####################################



### Enrichment term (derived from DEGs of SC vs. CT) scoring of the 12 superclusters 
### ...in each celltype at the rough level
### *** NOT USED FOR PUBLICAION
#####################################
###

### Load the GO terms for each supercluster
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
GOs_both_df=
  as.data.frame(unlist(GOs_both)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_both)') %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_pos=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[1]]
names(GOs_pos)=gsub("\\|.*","",names(GOs_pos))
GOs_pos_df=
  as.data.frame(unlist(GOs_pos)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_pos)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_neg=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[2]]
names(GOs_neg)=gsub("\\|.*","",names(GOs_neg))
GOs_neg_df=
  as.data.frame(unlist(GOs_neg)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_neg)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_all=data.table::rbindlist(list(GOs_pos_df, GOs_neg_df, GOs_both_df))

### Load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedDownregGenes.EnrichResults_SCvsCT.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_SCvsCT.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_SCvsCT.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_SCvsCT.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)

### Map the GO IDs to the 12 supercluster
merged_df=full_join(ego_all, GOs_all, by=c("GO_ID")) %>% 
  subset(!is.na(Cluster) & !is.na(supercluster)) %>%
  mutate(p.adjust_wSign=ifelse(orientation=="upreg",p.adjust,-p.adjust)) %>%
  group_by(Cluster, supercluster, orientation) %>%
  slice_max(p.adjust, n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(Cluster, supercluster) %>%
  summarize_at("p.adjust_wSign", sum) %>%
  left_join(data.frame(Cluster='B cells', supercluster=NA, p.adjust_wSign=NA)) %>%
  tidyr::pivot_wider(names_from="supercluster", id_cols="Cluster", values_from="p.adjust_wSign") %>%
  tibble::column_to_rownames("Cluster") %>%
  replace(is.na(.), 0)

### Plot
library(ComplexHeatmap)
mtx_=scale(merged_df, scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(8,10,11,4, 5,2,9, 3,1,7,6)]
# "cell division" is not there
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]

range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-3, 0, 3), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(mtx_,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Scaled enrichment index",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(6, "mm"),
          height=nrow(mtx_)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

#####################################



### Enrichment term (derived from DEGs of over80 vs. below80) scoring of the 12 superclusters 
### ...in each celltype at the rough level
### *** NOT USED FOR PUBLICAION
#####################################
###

library(dplyr)

### Load the GO terms for each supercluster
GOs_both=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster.rds")
names(GOs_both)=gsub("\\|.*","",names(GOs_both))
GOs_both_df=
  as.data.frame(unlist(GOs_both)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_both)') %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_pos=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[1]]
names(GOs_pos)=gsub("\\|.*","",names(GOs_pos))
GOs_pos_df=
  as.data.frame(unlist(GOs_pos)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_pos)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_neg=readRDS("~/Project_PBMCage/Tempt_RDS/GOterms_inEachSuperCluster_posAndnegRegulate.rds")[[2]]
names(GOs_neg)=gsub("\\|.*","",names(GOs_neg))
GOs_neg_df=
  as.data.frame(unlist(GOs_neg)) %>% 
  tibble::rownames_to_column("supercluster") %>% dplyr::rename(GO_ID='unlist(GOs_neg)') %>%
  # mutate(orientation="downreg") %>%
  mutate(supercluster=gsub("[0-9]","",supercluster))

GOs_all=data.table::rbindlist(list(GOs_pos_df, GOs_neg_df, GOs_both_df))

### Load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedDownregGenes.EnrichResults_At80.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_SharedUpregGenes.EnrichResults_At80.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_At80.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough_RbExcluded.EnrichResults_At80.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)

### Map the GO IDs to the 12 supercluster
merged_df=full_join(ego_all, GOs_all, by=c("GO_ID")) %>% 
  subset(!is.na(Cluster) & !is.na(supercluster)) %>%
  mutate(p.adjust_wSign=ifelse(orientation=="upreg",p.adjust,-p.adjust)) %>%
  group_by(Cluster, supercluster, orientation) %>%
  slice_max(p.adjust, n=1, with_ties=FALSE) %>%
  ungroup() %>% group_by(Cluster, supercluster) %>%
  summarize_at("p.adjust_wSign", sum) %>%
  left_join(data.frame(Cluster='B cells', supercluster=NA, p.adjust_wSign=NA)) %>%
  tidyr::pivot_wider(names_from="supercluster", id_cols="Cluster", values_from="p.adjust_wSign") %>%
  tibble::column_to_rownames("Cluster") %>%
  replace(is.na(.), 0)

### Plot
library(ComplexHeatmap)
mtx_=scale(merged_df, scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(8,10,11,4, 5,2,9, 3,1,7,6)]
# "cell division" is not there
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]

range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-3, 0, 3), c("dodgerblue3", "white", "brown3"))
ht_=
  Heatmap(mtx_,
          name="mat", show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Scaled enrichment index",
            legend_height=unit(4, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9),
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          border_gp=gpar(col="black", lty=2, lwd=0.5),
          rect_gp=gpar(col="white", lwd=0.5),
          col=col_fun_sub,
          cluster_rows=F, cluster_columns=F, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(6, "mm"),
          height=nrow(mtx_)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

#####################################



### Prove that NK cells show certain downregulation in SC vs. CT
#####################################
###

library(clusterProfiler)
library(ggplot2)

### Enrichment analysis on downreg./upreg. DEGs of SC vs. CT
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")
NK_downreg_genes=DEG_lists[["NK cells"]] %>% arrange(desc(avg_log2FC)) %>% subset(avg_log2FC<0)
NK_upreg_genes=DEG_lists[["NK cells"]] %>% arrange(desc(avg_log2FC)) %>% subset(avg_log2FC>0)
enrich_res_down=enrichGO(rownames(NK_downreg_genes), ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
dotplot(enrich_res_down)
# enrich_res_up=enrichGO(rownames(NK_upreg_genes), ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
# dotplot(enrich_res_up)

### Plot
module_=enrich_res_down@result %>% 
  .[1:20,] %>%
  mutate(GeneRatio_up=gsub("\\/.*","",GeneRatio), GeneRatio_down=gsub(".*\\/","",GeneRatio)) %>%
  mutate(GeneRatio=as.numeric(GeneRatio_up)/as.numeric(GeneRatio_down)) %>%
  dplyr::select(-c(GeneRatio_up, GeneRatio_down)) %>%
  arrange(desc(GeneRatio))
module_$Description=forcats::fct_relevel(module_$Description, module_$Description)
plot_NK=
  ggplot(module_, aes(y=GeneRatio, x=Description, color=p.adjust, size=Count)) +
  geom_point() +
  scale_x_discrete(labels=function(x) stringr::str_wrap(x, 40)) +
  scale_colour_gradient(limits=quantile(module_$p.adjust, na.rm=T)[c(1,5)],
                        low="brown4", high="pink",
                        labels=sprintf(fmt="%0.01e", quantile(module_$p.adjust, na.rm=T)[c(1,5)]),
                        breaks=quantile(module_$p.adjust, na.rm=T)[c(1,5)]) +
  scale_size_continuous(range=c(0.5,3), labels=c(60,120,180), breaks=c(60,120,180)) +
  labs(title="Downreg. genes in NKs from SC", y="GeneRatio", x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5, lineheight=0.8),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.4,"cm"),
        # legend.margin=margin(t=3, b=0),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical"
        
  ) +
  guides(color=guide_colorbar(label.theme=element_text(size=9),
                              title.theme=element_text(size=10), order=1, title="\n\n\n\n\np.adjust"),
         size=guide_legend(label.theme=element_text(size=9),
                           title.theme=element_text(size=10), order=2))

pdf("~/Project_PBMCage/Plots/JapanSC_Plots/DEGs_SC.vs.CT_EnrichGO_downreg.enrichment_NK.pdf", height=4, width=10)
plot_NK
dev.off()

#####################################



### Prove that CD4T cells show certain downregulation in SC vs. CT
#####################################
###

library(clusterProfiler)
library(ggplot2)

### Enrichment analysis on downreg./upreg. DEGs of SC vs. CT
DEG_lists=readRDS("~/Project_PBMCage/Japan supercentenarians/Tempt_RDS/scRNAseqAnalysis_DEGs_perAnnot.rough.rds")
CD4T_downreg_genes=DEG_lists[["CD4T cells"]] %>% arrange(desc(avg_log2FC)) %>% subset(avg_log2FC<0)
CD4T_upreg_genes=DEG_lists[["CD4T cells"]] %>% arrange(desc(avg_log2FC)) %>% subset(avg_log2FC>0)
enrich_res_down_CD4T=enrichGO(rownames(CD4T_downreg_genes), ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
dotplot(enrich_res_down_CD4T)
# enrich_res_up=enrichGO(rownames(NK_upreg_genes), ont="BP", OrgDb="org.Hs.eg.db", keyType="SYMBOL")
# dotplot(enrich_res_up)

### Plot
module_CD4T=enrich_res_down_CD4T@result %>% 
  .[1:5,] %>%
  mutate(GeneRatio_up=gsub("\\/.*","",GeneRatio), GeneRatio_down=gsub(".*\\/","",GeneRatio)) %>%
  mutate(GeneRatio=as.numeric(GeneRatio_up)/as.numeric(GeneRatio_down)) %>%
  dplyr::select(-c(GeneRatio_up, GeneRatio_down)) %>%
  arrange(desc(GeneRatio))
module_CD4T$Description=forcats::fct_relevel(module_CD4T$Description, module_CD4T$Description)
plot_CD4T=
  ggplot(module_CD4T, aes(y=GeneRatio, x=Description, color=p.adjust, size=Count)) +
  geom_point() +
  scale_x_discrete(labels=function(x) stringr::str_wrap(x, 40)) +
  scale_colour_gradient(limits=quantile(module_CD4T$p.adjust, na.rm=T)[c(1,5)],
                        low="brown4", high="pink",
                        labels=sprintf(fmt="%0.01e", quantile(module_CD4T$p.adjust, na.rm=T)[c(1,5)]),
                        breaks=quantile(module_CD4T$p.adjust, na.rm=T)[c(1,5)]) +
  scale_size_continuous(range=c(0.5,3), labels=c(70,90,110), breaks=c(70,90,110)) +
  labs(title="Downreg. genes in CD4T from SC", y="GeneRatio", x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=10, angle=90, hjust=1, vjust=0.5, lineheight=0.8),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="panel",
        plot.title=element_text(hjust=0.5),
        legend.key.width=unit(0.3,"cm"),
        legend.key.height=unit(0.4,"cm"),
        legend.box.spacing=unit(2,"cm"),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical"
        
  ) +
  guides(color=guide_colorbar(label.theme=element_text(size=9),
                              title.theme=element_text(size=10), order=1, title="\n\n\n\n\np.adjust"),
         size=guide_legend(label.theme=element_text(size=9),
                           title.theme=element_text(size=10), order=2))

pdf("~/Project_PBMCage/Plots/JapanSC_Plots/DEGs_SC.vs.CT_EnrichGO_downreg.enrichment_CD4T.pdf", height=3.8, width=4)
plot_CD4T
dev.off()

#####################################


