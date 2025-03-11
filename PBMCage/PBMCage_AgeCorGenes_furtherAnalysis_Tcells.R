
### Enrich age-correlated genes across CD4T and CD8T subpopulations (common ones)
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & grepl("CD4T\\.|CD8T\\.",celltypes) & rho_pval<0.05 & rho<0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_CD4CD8Detailed.rds")

upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & grepl("CD4T\\.|CD8T\\.",celltypes) & rho_pval<0.05 & rho>0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_CD4CD8Detailed.rds")

#####################################



### Enrich age-correlated genes after exclusion of ribosome-related ones
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_CD4CD8Detailed.rds")$ego
ego_results=ego_cluster@compareClusterResult %>% 
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]

### Exclude them and do EnrichGO
# downregulated ones
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0 & grepl("CD4T\\.|CD8T\\.",celltypes) & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check

# upregulated ones
upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0 & grepl("CD4T\\.|CD8T\\.",celltypes) & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check

# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_CD4CD8Detailed.rds")

#####################################



### Enrichment term scoring of the 12 superclusters in each T subpopulation at the detailed level
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
shared_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_CD4CD8Detailed.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_CD4CD8Detailed.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_CD4CD8Detailed.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_CD4CD8Detailed.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
# ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
#   mutate(p.adjust=-log10(p.adjust)) %>%
#   dplyr::rename(GO_ID=ID)
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df)) %>%
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
mtx_=scale(merged_df, scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
# mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
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
          column_title="Yazar et al. (2022) dataset", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(6, "mm"),
          height=nrow(mtx_)*unit(5.5, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_final=draw(ht_, heatmap_legend_side="left")

# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed.pdf", height=3, width=4.5)
# draw(ht_final)
# dev.off()

#####################################



### Correlation between Scaled Enrichment Index of the 12 superclusters in each T subpopulation at the detailed level,
### ...add age-correlated rho of cell proportions
#####################################
###

### Load the cell freq
# age-cor rho
COR_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_correlation.txt", sep="\t")
COR_detailed_CD4CD8T=COR_combined %>% subset(analysis=="all" & celltype.level=="Annot.detailed" & grepl("CD4T\\.|CD8T\\.",celltypes))
# freq data
# meta_info_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")

### Get the Scaled Enrichment Index
# load the GO terms for each supercluster
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
# load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_CD4CD8Detailed.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_CD4CD8Detailed.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_CD4CD8Detailed.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_CD4CD8Detailed.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
# ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
#   mutate(p.adjust=-log10(p.adjust)) %>%
#   dplyr::rename(GO_ID=ID)
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)
# map the GO IDs to the 12 supercluster
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

### Merge the enrichment index and cell freq
cor_between_EI_and_freq=dplyr::left_join(COR_detailed_CD4CD8T, 
                                         merged_df %>% tibble::rownames_to_column("celltypes"), 
                                         by="celltypes") %>%
  na.omit(.) %>% # remove those celltypes with NA in enrichment indexes
  tidyr::pivot_longer(cols=colnames(merged_df), names_to="Supercluster", values_to="Enrichment index")

### Plot general
ggplot(cor_between_EI_and_freq, aes(x=abs(cor), y=abs(`Enrichment index`), color=Supercluster)) +
  geom_point() +
  geom_smooth(aes(group=Supercluster))

### Plot by down/up-reg.
ggplot(cor_between_EI_and_freq %>% mutate(cor_orientation=ifelse(cor>0,"upreg.","downreg.")), 
       aes(x=abs(cor), y=`Enrichment index`, color=Supercluster)) +
  facet_wrap(~cor_orientation) +
  geom_point() +
  geom_smooth(aes(group=Supercluster))

### Plot by lympho/innate
ggplot(cor_between_EI_and_freq %>% mutate(cellgroup=ifelse(cor>0,"upreg.","downreg.")), 
       aes(x=abs(cor), y=`Enrichment index`, color=Supercluster)) +
  facet_wrap(~cor_orientation) +
  geom_point() +
  geom_smooth(aes(group=Supercluster))











### ===== Analysis of all celltypes at the detailed level =====



### Enrich age-correlated genes across Annot.detailed (common ones)
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed_updated.csv.gz", sep="\t")
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=downreg_genes, ego=ego_cluster),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Detailed_updated.rds")

upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe(.))
ego_cluster_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
saveRDS(list(genes=upreg_genes, ego=ego_cluster_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Detailed_updated.rds")

#####################################



### Enrich age-correlated genes after exclusion of ribosome-related ones
#####################################
###

library(clusterProfiler)
library(dplyr)
library(ggplot2)

### Load the correlation df
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed_updated.csv.gz", sep="\t")
# extract the enriched ribosome-related age-cor. genes
ego_cluster=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Detailed_updated.rds")$ego
ego_results=ego_cluster@compareClusterResult %>% 
  # subset(p.adjust<0.05) %>%
  subset(grepl("ribo|rRNA",Description))
ego_genes=sapply(ego_results$geneID, function(x) {strsplit(x, split="/")[[1]]})
names(ego_genes)=ego_results$Description
rb_related_genes=ego_genes %>% unlist() %>% unname() %>% .[!duplicated(.)]

### Exclude them and do EnrichGO
# downregulated ones
downreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho<0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_down=compareCluster(downreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_down) # check

# upregulated ones
upreg_genes=Cor_analysis_DF %>% subset(analysis=="all" & rho_pval<0.05 & rho>0 & !(gene %in% rb_related_genes)) %>%
  split(.$celltypes) %>%
  lapply(., function(df) df %>% select(gene) %>% tibble::deframe(.))
ego_up=compareCluster(upreg_genes, fun="enrichGO", OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
dotplot(ego_up) # check

# save
saveRDS(list(ego_down=ego_down, ego_up=ego_up),
        "~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Detailed_updated.rds")

#####################################



### Enrichment term scoring of the 12 superclusters in each celltype at the detailed level
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
shared_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Detailed_updated.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Detailed_updated.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Detailed_updated.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Detailed_updated.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
# ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
#   mutate(p.adjust=-log10(p.adjust)) %>%
#   dplyr::rename(GO_ID=ID)
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df)) %>%
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


width_unit=4
height_unit=4.5
### Plot for CD4T and CD8T cells
mtx_=scale(merged_df[grepl("CD4T\\.|CD8T\\.",rownames(merged_df)),], scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
mtx_row_reorder=c(mtx_row_reorder[grepl("naive",mtx_row_reorder)],
                  mtx_row_reorder[grepl("Tcm",mtx_row_reorder)],
                  mtx_row_reorder[grepl("Tem",mtx_row_reorder)],
                  mtx_row_reorder[grepl("Treg",mtx_row_reorder)],
                  mtx_row_reorder[grepl("prolif",mtx_row_reorder)])
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]
range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "brown3"))
mtx_=t(mtx_)
ht_Tcells=
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
          cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="T cells", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(width_unit, "mm"),
          height=nrow(mtx_)*unit(height_unit, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_Tcells_final=draw(ht_Tcells, heatmap_legend_side="left")



### Plot for NKs
mtx_=scale(merged_df[grepl("NK\\.",rownames(merged_df)),], scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
# mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]
range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "brown3"))
mtx_=t(mtx_)
ht_NKs=
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
          cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="NK cells", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(width_unit, "mm"),
          height=nrow(mtx_)*unit(height_unit, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_NKs_final=draw(ht_NKs, heatmap_legend_side="left")



### Plot for DCs and Monocytes
mtx_=scale(merged_df[grepl("DC\\.|Mono\\.",rownames(merged_df)),], scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
# mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]
range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "brown3"))
mtx_=t(mtx_)
ht_DCandMono=
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
          cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="DCs & Monocytes", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(width_unit, "mm"),
          height=nrow(mtx_)*unit(height_unit, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_DCandMono_final=draw(ht_DCandMono, heatmap_legend_side="left")



### Plot for B cells
mtx_=scale(merged_df[grepl("B\\.",rownames(merged_df)),], scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
# mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]
range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "brown3"))
mtx_=t(mtx_)
ht_Bcells=
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
          cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="B cells", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(width_unit, "mm"),
          height=nrow(mtx_)*unit(height_unit, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_Bcells_final=draw(ht_Bcells, heatmap_legend_side="left")



### Plot for OtherT cells
mtx_=scale(merged_df[grepl("OtherT\\.",rownames(merged_df)),], scale=T, center=F)
# reorder
mtx_row_reorder=names(table(rownames(mtx_)))
# mtx_row_reorder=mtx_row_reorder[c(1,2,3,6,8,4,5,7)]
mtx_col_reorder=names(table(colnames(mtx_)))
mtx_col_reorder=mtx_col_reorder[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
mtx_=mtx_[mtx_row_reorder, mtx_col_reorder]
range(mtx_) # check
col_fun_sub=circlize::colorRamp2(c(-4, 0, 4), c("dodgerblue3", "white", "brown3"))
mtx_=t(mtx_)
ht_OtherTs=
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
          cluster_rows=F, cluster_columns=T, show_column_dend=F, show_row_dend=F,
          column_names_rot=90,
          show_column_names=T,
          column_names_gp=gpar(fontsize=10), row_names_gp=gpar(fontsize=10),
          column_title="OtherT", column_title_gp=gpar(fontsize=11),
          width=ncol(mtx_)*unit(width_unit, "mm"),
          height=nrow(mtx_)*unit(height_unit, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3, "mm")
ht_OtherTs_final=draw(ht_OtherTs, heatmap_legend_side="left")



ht_final=draw(ht_Bcells+ ht_Tcells + ht_NKs + ht_DCandMono + ht_OtherTs, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_Annot.detailed.pdf", height=5, width=12)
draw(ht_final)
dev.off()

#####################################



### Correlation between Scaled Enrichment Index of the 12 superclusters in each celltype at the detailed level,
### ...add age-correlated rho of cell proportions
#####################################
###

library(dplyr)
library(ggplot2)

### Load the cell freq
# age-cor rho
COR_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_correlation.txt", sep="\t")
COR_detailed=COR_combined %>% subset(analysis=="all" & celltype.level=="Annot.detailed")
# # freq data (not good to distinguish between effector and supportive cells)
# meta_info_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
# meta_info_combined_arranged=meta_info_combined %>%
#   split(.$Annot.detailed) %>%
#   lapply(., function(df) {
#     tempt_=lm(percent.detailed ~ age, data=df)
#     tempt_$coefficients[2]
#   }) %>% 
#   as.data.frame() %>%
#   tidyr::pivot_longer(cols=colnames(.), names_to="Annot.detailed", values_to="Beta")

### Get the Scaled Enrichment Index
# load the GO terms for each supercluster
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
# load the enrichment results
shared_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedDownregGenes.EnrichResults_Detailed.rds")$ego
shared_down_df=shared_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
shared_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_SharedUpregGenes.EnrichResults_Detailed.rds")$ego
shared_up_df=shared_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
ego_down=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Detailed.rds")[[1]]
ego_down_df=ego_down@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="downreg")
ego_up=readRDS("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr.Cor.wAges_RbExcluded.EnrichResults_Detailed.rds")[[2]]
ego_up_df=ego_up@compareClusterResult %>% 
  dplyr::select(any_of(c("Cluster","ID","GeneRatio","p.adjust"))) %>%
  mutate(orientation="upreg")
# ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df, ego_down_df, ego_up_df)) %>%
#   mutate(p.adjust=-log10(p.adjust)) %>%
#   dplyr::rename(GO_ID=ID)
ego_all=data.table::rbindlist(list(shared_down_df, shared_up_df)) %>%
  mutate(p.adjust=-log10(p.adjust)) %>%
  dplyr::rename(GO_ID=ID)
# map the GO IDs to the 12 supercluster
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

### Merge the enrichment index and cell freq
cor_between_EI_and_freq=dplyr::left_join(COR_detailed, 
                                         merged_df %>% tibble::rownames_to_column("celltypes"), 
                                         by="celltypes") %>%
  na.omit(.) %>% # remove those celltypes with NA in enrichment indexes
  tidyr::pivot_longer(cols=colnames(merged_df), names_to="Supercluster", values_to="Enrichment index") %>%
  subset(!grepl("Other\\.Mast|Other\\.Eryth",celltypes)) %>%
  mutate(Annot.rough=gsub("\\..*","",celltypes),
         cor_orientation=ifelse(rho>0,"increased","declined")) %>%
  mutate(cell_group=ifelse(Annot.rough %in% c("CD4T","CD8T","NK","OtherT"), 
                           "Effector lymphocytes", 
                           "Supportive cells")) %>%
  mutate(Annot.rough=paste0(Annot.rough, " cells")) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="DC cells","DCs",
                            ifelse(Annot.rough=="Mono cells","Monocytes",
                                   ifelse(Annot.rough=="OtherT cells","OtherT",Annot.rough))))

### Plot general
plot_general_annot=
  ggplot(cor_between_EI_and_freq, aes(x=abs(rho), y=abs(`Enrichment index`), color=Annot.rough)) +
  geom_point(alpha=0.5, size=0.5, shape=20) +
  geom_smooth(aes(group=1), method="lm", linewidth=1.5, fill="gray93", color="black", show.legend=F) +
  ggpubr::stat_cor(aes(group=1)) +
  ggrepel::geom_text_repel(data=. %>% 
                             subset(cor_orientation=="downreg."|Annot.rough=="NK cells") %>%
                             group_by(celltypes) %>% 
                             slice_max(abs(`Enrichment index`)) %>%
                             ungroup() %>% slice_max(abs(`Enrichment index`)*abs(rho), n=6),
                           aes(label=celltypes), size=3, 
                           # nudge_y=2, 
                           min.segment.length=0, 
                           show.legend=F) +
  geom_smooth(aes(group=Annot.rough), method="lm", linewidth=0.5, se=F, show.legend=F) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1))) +
  labs(x=expression(Abs.Spearman~rho), y="Abs.Enrichment index") +
  theme_classic() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_EnrichmentIndex_vs_CellProportion_cor.pdf", height=4, width=4)
plot(plot_general_annot)
dev.off()

### Plot by down/up-reg.
plot_byOrientation=
  ggplot(cor_between_EI_and_freq, 
         aes(x=abs(rho), y=abs(`Enrichment index`), color=Annot.rough)) +
  facet_wrap(~cor_orientation, nrow=2) +
  geom_point(alpha=0.5, size=0.5, shape=20) +
  ggrepel::geom_text_repel(data=. %>% 
                             subset(cor_orientation=="declined"|Annot.rough=="NK cells") %>%
                             subset(grepl("naive|CD56dim",celltypes)) %>%
                             group_by(celltypes) %>% 
                             slice_max(abs(`Enrichment index`)) %>%
                             ungroup() %>% slice_max(abs(`Enrichment index`)*abs(rho), n=4),
                           aes(label=celltypes), size=3,
                           min.segment.length=0, 
                           show.legend=F) +
  geom_smooth(aes(group=Annot.rough), method="lm", linewidth=0.5, se=F, show.legend=F) +
  guides(color=guide_legend(title=NULL, override.aes=list(size=2, alpha=1))) +
  labs(x=expression(Abs.Spearman~rho), y="Abs.Enrichment index") +
  theme_classic() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.text=element_text(size=9, margin=margin(l=-6)),
        legend.box.spacing=unit(0,"line")) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_EnrichmentIndex_vs_CellProportion_cor_byOrientation.pdf", height=4, width=3.7)
plot(plot_byOrientation)
dev.off()

### Plot by grouping the celltypes
plot_byCellGroups=
  ggplot(cor_between_EI_and_freq, 
         aes(x=abs(rho), y=abs(`Enrichment index`), color=Annot.rough)) +
  facet_wrap(~cell_group, nrow=2) +
  geom_point(alpha=0.5, size=0.5, shape=20) +
  geom_smooth(aes(group=Annot.rough), method="lm", linewidth=0.5, se=F, show.legend=F) +
  guides(color="none") +
  labs(x=expression(Abs.Spearman~rho), y="Abs.Enrichment index") +
  theme_classic() +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.text=element_text(size=9)) +
  scale_color_manual(values=paletteer::paletteer_d("ggsci::category20_d3"))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Scoring_Phenotye_12Keyterms_basedOn_Pseudobulk_detailed_EnrichmentIndex_vs_CellProportion_cor_byCellGroup.pdf", height=4, width=2)
plot(plot_byCellGroups)
dev.off()

#####################################



### Correlation between upreg./downreg. gene counts in each celltype at the detailed level,
### ...and age-correlated rho of cell proportions
#####################################
###

### Load the cell freq
# age-cor rho
COR_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_correlation.txt", sep="\t")
COR_detailed=COR_combined %>% subset(analysis=="all" & celltype.level=="Annot.detailed")
# freq data
# meta_info_combined=read.table("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")

### Load the age-corr. genes
Cor_analysis_DF=data.table::fread("~/Project_PBMCage/Tempt_RDS/Bulkanalysis_RNAexpr_Correlation_wAges_Detailed.csv.gz", sep="\t")
celltypes_=names(table(Cor_analysis_DF$celltypes))
geneNO=Cor_analysis_DF %>%
  subset(rho_pval<0.05 & analysis=="all") %>%
  mutate(orientation=ifelse(rho<0,"downreg.","upreg.")) %>%
  split(.$orientation) %>%
  .[c("downreg.","upreg.")] %>%
  lapply(., function(df) df %>% split(.$celltypes) %>% .[celltypes_] %>%
           lapply(., function(dff) dff$gene %>% unique(.) %>% length(.)))
geneNO_arranged=lapply(geneNO, function(list_) as.data.frame(list_) %>% 
                         tidyr::pivot_longer(cols=colnames(.),names_to="Annot.detailed",values_to="geneN"))
geneNO_arranged[[1]]=geneNO_arranged[[1]] %>% mutate(orientation="downreg.")
geneNO_arranged[[2]]=geneNO_arranged[[2]] %>% mutate(orientation="upreg.")
geneNO_arranged_df=data.table::rbindlist(geneNO_arranged)

### Make the correlation df
COR_df=COR_detailed %>% dplyr::rename(Annot.detailed=celltypes) %>%
  full_join(., geneNO_arranged_df, by="Annot.detailed") %>%
  mutate(Annot.rough=gsub("\\..*","",Annot.detailed)) %>%
  mutate(cell_group=ifelse(Annot.rough %in% c("CD4T","CD8T","NK","OtherT"), 
                           "Effector lymphocytes", 
                           "Supportive cells")) %>%
  mutate(Annot.rough=paste0(Annot.rough, " cells")) %>%
  mutate(Annot.rough=ifelse(Annot.rough=="DC cells","DCs",
                            ifelse(Annot.rough=="Mono cells","Monocytes",
                                   ifelse(Annot.rough=="OtherT cells","OtherT",Annot.rough))))

### Plot
# general
COR_df_upAnddown=COR_df %>% 
  group_by_if(!grepl("geneN|orientation",colnames(.))) %>%
  summarize_at("geneN", sum)
ggplot(COR_df_upAnddown, aes(x=abs(cor), y=geneN, color=Annot.rough)) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=Annot.detailed))

# by upreg. and downreg.
ggplot(COR_df, aes(x=abs(cor), y=geneN, color=Annot.rough)) +
  facet_wrap(~orientation) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=Annot.detailed))

# by cell characteristics
ggplot(COR_df, aes(x=abs(cor), y=geneN, color=Annot.rough)) +
  facet_wrap(~cell_group) +
  geom_point() +
  ggrepel::geom_text_repel(aes(label=Annot.detailed)) +
  geom_smooth(aes(group=Annot.rough), method="lm")
