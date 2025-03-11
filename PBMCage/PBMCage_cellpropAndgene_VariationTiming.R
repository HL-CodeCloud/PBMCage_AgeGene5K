
### Cluster the variation of genes across ages in pseudobulk.donor_id, regardless of sex
#####################################
###

library(dplyr)
library(Seurat)
library(ClusterGVis)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_donor.rds")
pseudobulkobj_all$age=as.numeric(colnames(pseudobulkobj_all))

### Run ClusterGVis
pseudo_data=LayerData(pseudobulkobj_all, layer="data")
# run getClusters to determine the best cluster.num
pdf("~/Project_PBMCage/Plots/Tempt.pdf", height=15, width=15)
getClusters(exp=pseudo_data) # so take cluster.num=4 or 5
dev.off()

ck=clusterData(exp=pseudo_data,
               cluster.method="kmeans",
               cluster.num=4)
saveRDS(ck, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_CorGenes.CalculatedonDetailedcelltype_BasedOnPseudobulkDonorid_bothSex.rds")

#####################################



### Cluster the variation of genes across ages in pseudobulk.rough, regardless of sex
#####################################
###

library(dplyr)
library(Seurat)
library(ClusterGVis)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.rough))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.rough==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Run ClusterGVis
CKs=list()
for (i in 1:length(all_celltypes)) {
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # remove the genes with two small expr or variation
  gene_expr_max=apply(pseudo_data, 1, max)
  pseudo_data=pseudo_data[gene_expr_max>0.1, ]
  # run getClusters to determine the best cluster.num
  nbclust_res=factoextra::fviz_nbclust(pseudo_data, FUNcluster=kmeans, method="wss", k.max=10) # check the optimal number of clusters
  slope=sapply(2:length(-diff(nbclust_res$data$y)), function(idx) -diff(nbclust_res$data$y)[idx-1]/-diff(nbclust_res$data$y)[idx])
  k_best=which(slope<2)[1]
  ck_save=function(data) {
    library(ClusterGVis)
    tryCatch({
      ck_=clusterData(exp=data,
                      cluster.method="kmeans",
                      cluster.num=k_best)
      return(ck_)
    }, error=function(msg) {
      print("Oversize.")
      return(NA)
    })
  }
  CKs[[i]]=ck_save(pseudo_data)
}
names(CKs)=all_celltypes
saveRDS(CKs, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex.rds")

#####################################



### Pattern clustering of the pseudobulk.rough_both sex's aging trends
#####################################
###

### Load the aging trends
CKs=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex.rds")

### Compare the trends across celltypes
# get the median as the trends
celltype_shortname=c("B_","CD4T_","CD8T_","DC_","Mono_","NK_","Other_","OtherT_")
CK_medians=c()
for (i in 1:length(CKs)) {
  ck_median=CKs[[i]]$long.res %>% group_by(cluster, cell_type) %>% summarize(median=median(norm_value, rm.na=TRUE)) %>%
    ungroup() %>%
    split(.$cluster) %>% lapply(., function(df) df %>% arrange(cell_type) %>% dplyr::select(cell_type, median) %>% tibble::deframe())
  names(ck_median)=paste0(celltype_shortname[i],names(ck_median))
  CK_medians=c(CK_medians, ck_median)
}
# arrange the df
CK_medians_df=lapply(1:length(CK_medians), function(idx) {df=data.frame(age=names(CK_medians[[idx]]), median=CK_medians[[idx]]); colnames(df)=c("age",names(CK_medians)[idx]); df})
CK_medians_df=Reduce(full_join, CK_medians_df)
# compare
median_names=colnames(CK_medians_df)[2:ncol(CK_medians_df)]
tau=tau_names=pval=c()
for (i in 1:length(median_names)) {
  for (j in 1:length(median_names)) {
    res=cor.test(CK_medians_df[[median_names[i]]], CK_medians_df[[median_names[j]]], method="kendall")
    tau=c(tau, res$estimate %>% unname())
    pval=c(pval, res$p.value)
    tau_names=c(tau_names, paste0(median_names[i],"_vs_",median_names[j]))
  }
}
# rearrange the correlation results
cor_mtx=data.frame(pair=tau_names, tau_val=tau) %>%
  mutate(group1=gsub("_vs_.*","",pair), group2=gsub(".*_vs_","",pair)) %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="group2", values_from="tau_val") %>%
  tibble::column_to_rownames("group1")
pval_mtx=data.frame(pair=tau_names, p_val=pval) %>%
  mutate(group1=gsub("_vs_.*","",pair), group2=gsub(".*_vs_","",pair)) %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="group2", values_from="p_val") %>%
  tibble::column_to_rownames("group1")

### Plot
quantile(as.vector(cor_mtx) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-1, 0, 1), c("dodgerblue3", "white", "brown3"))
# cluster
pa=cluster::pam(cor_mtx, k=5)
ht_cormtx=
  Heatmap(cor_mtx,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title=expression(Kendall~tau),
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun_cor,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=T, cluster_rows=T,
          row_names_gp=gpar(fontsize=10), show_row_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          height=unit(4.8, "mm")*nrow(cor_mtx),
          width=unit(4.8, "mm")*ncol(cor_mtx),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm"),
          row_split=paste0("pattern", pa$clustering),
          column_split=paste0("pattern", pa$clustering)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_cormtx=draw(ht_cormtx, heatmap_legend_side="top")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr_AgingTrendPatterns_Rough_CorrelationHeatmap.pdf", height=8.5, width=8)
ht_cormtx
dev.off()
# save the column/row order
order_=ComplexHeatmap::column_order(ht_cormtx)
order_names=lapply(order_, function(x) colnames(cor_mtx)[x])
saveRDS(order_names, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_Patterns_CorrelationHeatmapOrder.rds")

### Extract the aging trends within the 5 patterns
cluster_df=as.data.frame(pa$clustering) %>% tibble::rownames_to_column("agingtrend") %>% rename(cluster=`pa$clustering`)
pattern_clusters=split(cluster_df$agingtrend, cluster_df$cluster)
names(pattern_clusters)=paste0("pattern",names(pattern_clusters))

### Plot the aging trends within each pattern
ages=CK_medians_df[["age"]]
CK_medians_df_perpattern=lapply(pattern_clusters, 
                                function(p) {
                                  medians=CK_medians_df[,p]
                                  medians$age=ages
                                  medians=medians %>% tidyr::pivot_longer(cols=colnames(.)[!grepl("age",colnames(.))], names_to="agingtrend", values_to="median")
                                  medians
                                })
CK_medians_df_perpattern=lapply(1:length(CK_medians_df_perpattern), function(idx) CK_medians_df_perpattern[[idx]] %>% mutate(pattern=names(CK_medians_df_perpattern)[idx]))
median_pattern_all=data.table::rbindlist(CK_medians_df_perpattern)
median_pattern_all$age=as.integer(median_pattern_all$age)
pattern_plot=
  ggplot(median_pattern_all, aes(x=age, y=median, color=agingtrend)) +
  facet_wrap(~pattern, ncol=1, strip.position="right") +
  geom_line(aes(group=agingtrend), linewidth=0.25, alpha=0.5) +
  geom_smooth(aes(group=1), color="grey50", linewidth=0.5) +
  labs(y="Median of norm. expr.") +
  theme_bw() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=10),
        legend.position="top",
        legend.spacing.x=unit(0.1,"cm"),
        legend.direction="horizontal") +
  guides(color="none") +
  scale_x_continuous(labels=c(30,60,90), breaks=c(30,60,90))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Bulkanalysis_RNAexpr_AgingTrendPatterns_Rough_ExprAcrossAging.LinePlot.pdf", height=5, width=1.5)
plot(pattern_plot)
dev.off()

### Save the patterns identified based on the aging trends
saveRDS(median_pattern_all, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_Patterns.rds")

#####################################



### Extract the genes within each agingtrend belong to each pattern
#####################################
###

### Load the aging trends
CKs=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex.rds")
# extract the genes from each aging trend
celltype_shortname=c("B_","CD4T_","CD8T_","DC_","Mono_","NK_","Other_","OtherT_")
ck_genes_all=c()
for (i in 1:length(CKs)) {
  ck_=CKs[[i]]
  ck_genes=ck_$long.res %>% 
    dplyr::select(cluster, gene) %>%
    subset(!duplicated(.)) %>%
    split(.$cluster) %>% lapply(., function(df) df %>% dplyr::select(gene) %>% tibble::deframe())
  names(ck_genes)=paste0(celltype_shortname[i], names(ck_genes))
  ck_genes_all=c(ck_genes_all, ck_genes)
}
# group them into the 5 patterns
pattern_clusters=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_Patterns.rds")
pattern_clusters=pattern_clusters %>% dplyr::select(agingtrend, pattern) %>% subset(!duplicated(.))
pattern_clusters=split(pattern_clusters$agingtrend, pattern_clusters$pattern)
pattern_genes=lapply(pattern_clusters, function(cl) ck_genes_all[cl])
saveRDS(pattern_genes, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")

### Rearrange to a dataframe
pattern_trend_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
# list all the genes from each cell type within each pattern
pattern_trend_genes_df=data.frame(group=names(unlist(pattern_trend_genes)), gene=unlist(pattern_trend_genes))
pattern_trend_genes_df=pattern_trend_genes_df %>%
  mutate(pattern=gsub("\\..*","",group)) %>%
  mutate(celltype=gsub("^pattern[0-9]\\.|_.*","",group)) %>%
  mutate(agingtrend=gsub("^\\.","",stringr::str_match(group, "\\..+_[1-5]")[,1]))
table(pattern_trend_genes_df$agingtrend) # check
table(pattern_trend_genes_df$pattern)
table(pattern_trend_genes_df$celltype)
data.table::fwrite(pattern_trend_genes_df, "~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.txt.gz", sep="\t")

#####################################



### Calculate the proportion of the genes (from each celltype within each pattern) in 12 superclusters
#####################################
###

### Map the genes to the superclusters
pattern_trend_genes_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.txt.gz", sep="\t")
gene_supercluster=readRDS("~/Project_PBMCage/Tempt_RDS/GOSuperClusters_genes_inEachSuperCluster.rds")
names(gene_supercluster)=gsub("\\|.*","",names(gene_supercluster))
gene_supercluster_df=lapply(1:length(gene_supercluster), function(idx) data.frame(supercluster=names(gene_supercluster)[idx], gene=gene_supercluster[[idx]])) %>%
  data.table::rbindlist()
pattern_gene_mapTOsupercluster=
  full_join(pattern_trend_genes_df, gene_supercluster_df, by="gene") %>%
  subset(!is.na(pattern) & !is.na(supercluster))

### Arrange and calculate the numbers of genes in each supercluster
pattern_gene_arranged=pattern_gene_mapTOsupercluster %>%
  mutate(pattern.celltype=paste0(pattern,"_|_",celltype)) %>%
  split(.$pattern.celltype)
pattern_gene_arranged=lapply(1:length(pattern_gene_arranged), function(idx) pattern_gene_arranged[[idx]] %>% 
                               mutate(n.total=nrow(.)) %>% group_by(n.total, supercluster) %>% summarize(n.persupercluster=n()) %>%
                               mutate(pattern.celltype=names(pattern_gene_arranged)[idx])) %>%
  data.table::rbindlist() %>%
  mutate(freq=n.persupercluster/n.total*100) %>%
  mutate(pattern=gsub("_\\|_.*","",pattern.celltype), celltype=gsub(".*_\\|_","",pattern.celltype))

reorder_=names(table(pattern_gene_arranged$supercluster))[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
pattern_gene_arranged$supercluster=forcats::fct_relevel(pattern_gene_arranged$supercluster, reorder_)
piechart_=
  ggplot(pattern_gene_arranged, aes(x="", y=n.persupercluster, fill=supercluster)) +
  facet_grid(pattern~celltype) +
  geom_col(color="white", linewidth=0.1) +
  scale_fill_manual(values=paletteer::paletteer_d("ggthemes::Classic_Green_Orange_12")) +
  coord_polar(theta="y") +
  theme_minimal() +
  theme(legend.position="left",
        legend.direction="vertical",
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        strip.text=element_text(size=9),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        panel.spacing.x=unit(1,"mm"),
        panel.spacing.y=unit(10,"mm")) +
  labs(x=NULL, y=NULL) +
  guides(fill=guide_legend(title=NULL))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_Piechart.pdf", width=6, height=5)
piechart_
dev.off()

#####################################



### Enrich the genes from each celltype within each pattern
#####################################
###

### Load the df and extract the genes
pattern_trend_genes_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.txt.gz", sep="\t")
pattern_trend_genes=pattern_trend_genes_df %>%
  mutate(pattern.celltype=paste0(pattern,"_|_",celltype)) %>%
  split(.$pattern.celltype) %>%
  lapply(., function(df) df$gene %>% .[!duplicated(.)])

### Enrich
goenrich_res=lapply(pattern_trend_genes, function(x) clusterProfiler::enrichGO(x, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP"))
names(goenrich_res)=names(pattern_trend_genes)
goenrich_res_df=lapply(1:length(goenrich_res), function(idx) goenrich_res[[idx]]@result %>% 
                         subset(p.adjust<0.05) %>% dplyr::select(ID, Description, GeneRatio, p.adjust, Count) %>%
                         mutate(group=names(goenrich_res)[idx])) %>%
  data.table::rbindlist() %>%
  mutate(pattern=gsub("_\\|_.*","",group), celltype=gsub(".*_\\|_","",group))
goenrich_res_df$GeneRatio=lapply(strsplit(goenrich_res_df$GeneRatio, split="/"), function(x) as.numeric(x[1])/as.numeric(x[2])) %>%
  unlist()
data.table::fwrite(goenrich_res_df, "~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_enrichGOres.txt.gz", sep="\t")

### Plot
goenrich_res_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_enrichGOres.txt.gz", sep="\t")
# remove the replicated terms
goenrich_res_df_rep=goenrich_res_df %>% split(.$group) %>% lapply(., function(df) df$ID[1:30])
goenrich_res_specificGO=lapply(1:length(goenrich_res_df_rep), function(idx) data.frame(group=names(goenrich_res_df_rep)[idx], ID=goenrich_res_df_rep[[idx]])) %>%
  data.table::rbindlist() %>%
  group_by(ID) %>%
  summarize(replicates=n()) %>%
  subset(replicates<=3) %>%
  dplyr::select(ID) %>% unlist() %>% unname()
goenrich_res_df_rep.rm=goenrich_res_df %>%
  subset(ID %in% goenrich_res_specificGO) %>%
  slice_max(GeneRatio, n=3, with_ties=F, by="group") %>%
  group_by(group) %>%
  mutate(id=row_number())
plot_=
  ggplot(goenrich_res_df_rep.rm, 
         aes(x=id, y=0, color=p.adjust)) +
    facet_grid(pattern~celltype) +
    geom_point(size=2) +
    coord_flip(ylim=c(-0.01,5), xlim=c(0.5,3.5)) +
    geom_text(aes(label=stringr::str_wrap(Description, width=25, indent=3, exdent=3), y=0), 
              hjust=0, size=3, lineheight=0.75, color="black") +
    theme_void() +
    scale_color_continuous(limits=
                            c(min(goenrich_res_df_rep.rm$p.adjust), max(goenrich_res_df_rep.rm$p.adjust)),
                          low="lightpink", high="brown3") +
    theme(legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          strip.text.y.right=element_text(angle=-90, hjust=0.5, vjust=0.5, size=9),
          strip.text.x=element_text(size=9),
          strip.background=element_rect(color="grey50"),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.key.height=unit(1,"mm"),
          legend.key.width=unit(3,"cm"),
          axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA, linewidth=0.5))
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_EnrichPlot.pdf", width=12.5, height=6)
plot_
dev.off()

#####################################



### Similarity analysis on the genes across aging trends in each clustered pattern based on the pseudobulk.rough_both sex
#####################################
###

library(ComplexHeatmap)

### Match the column/row order between similarity heatmap and pattern heatmap
order_names=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_Patterns_CorrelationHeatmapOrder.rds")
order_names=order_names[sort(names(order_names))]

### Jaccard similarity analysis
pattern_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
# the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
# analyze
Jaccard_SIM_all=list()
for (p in 1:length(pattern_genes)) {
  pattern=pattern_genes[[p]]
  # similarity between each pair of aging trend within this pattern
  SIM=SIM_Name=c()
  for (i in 1:length(pattern)) {
    for (j in 1:length(pattern)) {
      similarity_result=jaccard(pattern[[i]], pattern[[j]])
      SIM=c(SIM, similarity_result)
      SIM_Name=c(SIM_Name, paste0(names(pattern)[i],"_vs_",names(pattern)[j]))
    }
  }
  names(SIM)=SIM_Name
  # arrange the similarity matrix
  Jaccard_sim_=data.frame(pair=names(SIM), jaccard=SIM) %>%
    mutate(trend1=gsub("_vs_.*","",pair), trend2=gsub(".*_vs_","",pair)) %>%
    dplyr::select(-pair) %>%
    tidyr::pivot_wider(names_from="trend2", values_from="jaccard") %>%
    tibble::column_to_rownames("trend1")
    
  Jaccard_SIM_all[[p]]=Jaccard_sim_
}
names(Jaccard_SIM_all)=names(pattern_genes)

### Plot the similarity heatmap
lapply(Jaccard_SIM_all, range) # check
col_fun=circlize::colorRamp2(c(0, 1), c("white","brown3"))
ht_genesims=list()
for (i in 1:length(Jaccard_SIM_all)) {
  ht_genesims[[i]]=
    Heatmap(Jaccard_SIM_all[[i]],
            name=" ",
            show_heatmap_legend=c(rep(FALSE,3), TRUE, FALSE)[i],
            heatmap_legend_param=list(
              title="Gene-level similarity",
              legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
              direction="horizontal", 
              title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
            ),
            border_gp=gpar(col="black", lty=2),
            col=col_fun,
            rect_gp=gpar(col="white", lwd=0.5),
            cluster_columns=F, cluster_rows=F,
            column_order=match(order_names[[i]], colnames(Jaccard_SIM_all[[i]])), 
            row_order=match(order_names[[i]], rownames(Jaccard_SIM_all[[i]])),
            row_names_gp=gpar(fontsize=10), show_row_names=T,
            row_title=names(Jaccard_SIM_all)[i],
            row_title_gp=gpar(fontsize=11),
            show_column_names=F,
            column_names_side="bottom",
            column_names_rot=90,
            column_names_gp=gpar(fontsize=10),
            height=unit(5, "mm")*nrow(Jaccard_SIM_all[[i]]),
            width=unit(5, "mm")*ncol(Jaccard_SIM_all[[i]]))
}
ht_genesims=lapply(ht_genesims, function(x) draw(x, heatmap_legend_side="top"))

pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_IndividualGeneSim_withinPattern.pdf", height=3, width=2.5)
ht_genesims[[1]]
ht_genesims[[2]]
ht_genesims[[3]]
ht_genesims[[4]]
ht_genesims[[5]]
dev.off()

### Perform functional similarity analysis (mclusterSim) on aging trends within each pattern
library(dplyr)
pattern_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
pattern_genes_entrez=
  lapply(pattern_genes, 
         function(pattern) lapply(pattern, function(genelist) clusterProfiler::bitr(genelist, fromType="SYMBOL", OrgDb="org.Hs.eg.db", toType="ENTREZID") %>% .[["ENTREZID"]]))
names(pattern_genes_entrez)=names(pattern_genes)
# analyze functional similarity with mclusterSim
d=GOSemSim::godata('org.Hs.eg.db', ont="BP")
res=lapply(pattern_genes_entrez, function(agingtrend) GOSemSim::mclusterSim(agingtrend, semData=d, measure="Wang", combine="BMA"))
saveRDS(res, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_mclusterSim_withinPattern.rds")

### Plot
SIM=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_mclusterSim_withinPattern.rds")
range(SIM)
col_fun=circlize::colorRamp2(c(0.9, 1), c("white","brown3"))
ht_genesims=list()
for (i in 1:length(SIM)) {
  ht_genesims[[i]]=
    Heatmap(SIM[[i]],
            name=" ",
            show_heatmap_legend=c(rep(FALSE,3), TRUE, FALSE)[i],
            heatmap_legend_param=list(
              title="Cluster-level functional similarity",
              legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
              direction="horizontal", 
              title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
            ),
            border_gp=gpar(col="black", lty=2),
            col=col_fun,
            rect_gp=gpar(col="white", lwd=0.5),
            cluster_columns=F, cluster_rows=F,
            column_order=match(order_names[[i]], colnames(SIM[[i]])), 
            row_order=match(order_names[[i]], rownames(SIM[[i]])),
            row_names_gp=gpar(fontsize=10), show_row_names=T,
            row_title=names(SIM)[i],
            row_title_gp=gpar(fontsize=11),
            show_column_names=F,
            column_names_side="bottom",
            column_names_rot=90,
            column_names_gp=gpar(fontsize=10),
            height=unit(5, "mm")*nrow(SIM[[i]]),
            width=unit(5, "mm")*ncol(SIM[[i]]))
}
ht_genesims=lapply(ht_genesims, function(x) draw(x, heatmap_legend_side="top"))

pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_mclusterSim_withinPattern.pdf", height=3, width=2.5)
ht_genesims[[1]]
ht_genesims[[2]]
ht_genesims[[3]]
ht_genesims[[4]]
ht_genesims[[5]]
dev.off()

#####################################



### Similarity analysis on the genes across 5 patterns based on the pseudobulk.rough_both sex
#####################################
###

library(ComplexHeatmap)

### Get the genes (across celltypes in general) within each pattern
pattern_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
genes_per_pattern=lapply(pattern_genes, function(item) {
  names_=names(item)
  names_celltypes=gsub("_[0-9]","",names_)
  match_=match(names_celltypes, unique(names_celltypes))
  merged=sapply(1:max(match_), function(idx) item[match_==idx] %>% Reduce(union, .))
  names(merged)=unique(names_celltypes)
  merged
})
genes_per_pattern_intersect=lapply(genes_per_pattern, function(item) Reduce(intersect, item))

### Jaccard similarity analysis
# the jaccard similarity function
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
# similarity between each pair of patterns
SIM=SIM_Name=c()
for (i in 1:length(genes_per_pattern_intersect)) {
  for (j in 1:length(genes_per_pattern_intersect)) {
    similarity_result=jaccard(genes_per_pattern_intersect[[i]], genes_per_pattern_intersect[[j]])
    SIM=c(SIM, similarity_result)
    SIM_Name=c(SIM_Name, paste0(names(genes_per_pattern_intersect)[i],"_vs_",names(genes_per_pattern_intersect)[j]))
  }
}
names(SIM)=SIM_Name
# arrange the similarity matrix
Jaccard_sim_=data.frame(pair=names(SIM), jaccard=SIM) %>%
  mutate(trend1=gsub("_vs_.*","",pair), trend2=gsub(".*_vs_","",pair)) %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="trend2", values_from="jaccard") %>%
  tibble::column_to_rownames("trend1")

### Plot the similarity heatmap
lapply(Jaccard_sim_, range) # check
col_fun=circlize::colorRamp2(c(0, 1), c("white","brown3"))
ht_genesim_=
  Heatmap(Jaccard_sim_,
          name=" ",
          show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="Individual gene similarity",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T,
          row_title_gp=gpar(fontsize=11),
          column_names_side="bottom", column_names_rot=90,
          column_names_gp=gpar(fontsize=10), show_column_names=T,
          height=unit(5, "mm")*nrow(Jaccard_sim_),
          width=unit(5, "mm")*ncol(Jaccard_sim_))
ht_genesim_=draw(ht_genesim_, heatmap_legend_side="top")

pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_IndividualGeneJaccardSim_AcrossPattern.pdf", height=2.5, width=2)
ht_genesim_
dev.off()

### Perform functional similarity analysis on across-celltype genes within each pattern
enrichs=
  lapply(genes_per_pattern_intersect,
         function(genelist) clusterProfiler::enrichGO(genelist, keyType="SYMBOL", OrgDb="org.Hs.eg.db", ont="BP"))
# # analyze functional similarity with mgosim
# go_terms=lapply(enrichs, function(x) (x@result %>% subset(p.adjust<0.05))$ID)
# d=GOSemSim::godata('org.Hs.eg.db', ont="BP")
# GoSIM=GoSIM_Name=c()
# for (i in 1:length(go_terms)) {
#   for (j in 1:length(go_terms)) {
#     similarity_result=GOSemSim::mgoSim(go_terms[[i]], go_terms[[j]], semData=d, measure="Wang", combine="BMA")
#     GoSIM=c(GoSIM, similarity_result)
#     GoSIM_Name=c(GoSIM_Name, paste0(names(go_terms)[i],"_vs_",names(go_terms)[j]))
#   }
# }
# names(GoSIM)=GoSIM_Name
# jaccard similarity analysis on go terms
jaccard=function(a, b) {
  intersection=length(intersect(a, b))
  union=length(a)+length(b)-intersection
  return (intersection/union)
}
GoSIM=GoSIM_Name=c()
go_terms=lapply(enrichs, function(x) (x@result %>% subset(p.adjust<0.05))$ID)
for (i in 1:length(go_terms)) {
  for (j in 1:length(go_terms)) {
    similarity_result=jaccard(go_terms[[i]], go_terms[[j]])
    GoSIM=c(GoSIM, similarity_result)
    GoSIM_Name=c(GoSIM_Name, paste0(names(go_terms)[i],"_vs_",names(go_terms)[j]))
  }
}
names(GoSIM)=GoSIM_Name
Go_sim_=data.frame(pair=names(GoSIM), similarity=GoSIM) %>%
  mutate(trend1=gsub("_vs_.*","",pair), trend2=gsub(".*_vs_","",pair)) %>%
  dplyr::select(-pair) %>%
  tidyr::pivot_wider(names_from="trend2", values_from="similarity") %>%
  tibble::column_to_rownames("trend1")

### Plot
range(Go_sim_)
col_fun=circlize::colorRamp2(c(0, 1), c("white","brown3"))
ht_gosim_=
  Heatmap(Go_sim_,
          name=" ",
          show_heatmap_legend=T,
          heatmap_legend_param=list(
            title="GO term similarity",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T,
          row_title_gp=gpar(fontsize=11),
          column_names_side="bottom", column_names_rot=90,
          column_names_gp=gpar(fontsize=10), show_column_names=T,
          height=unit(5, "mm")*nrow(Go_sim_),
          width=unit(5, "mm")*ncol(Go_sim_))
ht_gosim_=draw(ht_gosim_, heatmap_legend_side="top")

pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_GOJaccardSim_AcrossPattern.pdf", height=2.5, width=2)
ht_gosim_
dev.off()

#####################################



### GO enrichments across 5 patterns based on the pseudobulk.rough_both sex
#####################################
###

library(dplyr)

### Get the genes (across celltypes in general) within each pattern
pattern_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
genes_per_pattern=lapply(pattern_genes, function(item) {
  names_=names(item)
  names_celltypes=gsub("_[0-9]","",names_)
  match_=match(names_celltypes, unique(names_celltypes))
  merged=sapply(1:max(match_), function(idx) item[match_==idx] %>% Reduce(union, .))
  names(merged)=unique(names_celltypes)
  merged
})
genes_per_pattern_intersect=lapply(genes_per_pattern, function(item) Reduce(union, item))

### Enrich
enrichs=
  lapply(genes_per_pattern_intersect,
         function(genelist) clusterProfiler::enrichGO(genelist, keyType="SYMBOL", OrgDb="org.Hs.eg.db", ont="BP"))
clusterProfiler::dotplot(enrichs[[1]])
enrichs_compare=clusterProfiler::compareCluster(genes_per_pattern_intersect, fun="enrichGO", keyType="SYMBOL", OrgDb="org.Hs.eg.db", ont="BP")
clusterProfiler::dotplot(enrichs_compare, showCategory=5)

### One hot encoder
goenrich_res_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_enrichGOres.txt.gz", sep="\t")
goenrich_res_df_arranged=goenrich_res_df %>%
  group_by(group) %>%
  slice_max(GeneRatio, n=50, with_ties=F) %>%
  mutate(ID_rank=row_number()) %>%
  mutate(ID_rank=paste0("rank",ID_rank)) %>%
  dplyr::select(ID, group, ID_rank) %>%
  tidyr::pivot_wider(id_cols="group", names_from="ID_rank", values_from="ID") %>%
  tibble::column_to_rownames("group")
library(caret)
dmy=caret::dummyVars(" ~ .", data=goenrich_res_df_arranged)
trsf=data.frame(predict(dmy, newdata=goenrich_res_df_arranged))
trsf_arranged=apply(trsf, 1, function(x) paste0(x, collapse=""))
trsf_arranged_convert=Gmisc::baseConvert(trsf_arranged, target=10, base=2)
trsf_arranged_convert=data.frame(group=names(trsf_arranged), trsf_arranged_convert) %>%
  mutate(celltype=gsub(".*_\\|_","",group), pattern=gsub("_\\|_.*","",group)) %>%
  dplyr::select(-group) %>%
  tidyr::pivot_wider(id_cols="pattern", names_from="celltype", values_from="trsf_arranged_convert") %>%
  tibble::column_to_rownames("pattern")

### SVD
trsf_mtx=as.matrix(trsf_arranged_convert)
trsf_mtx[is.na(trsf_mtx)]=0
trsf_mtx.svd=svd(trsf_mtx)

### Convert back to base2
tempt_convert=Gmisc::baseConvert(trsf_mtx.svd$d, target=2, base=10)
lacking_zeros=ncol(trsf)-nchar(tempt_convert)
svd_d_convert=
  sapply(1:length(lacking_zeros), function(idx) ifelse(lacking_zeros[idx]<0, gsub("^[01]","",tempt_convert[idx]),
                                                       ifelse(lacking_zeros[idx]==0, tempt_convert[idx], 
                                                              paste0(rep("0",lacking_zeros[[idx]]), tempt_convert[idx]))))
# map to the GO terms
base2_convertTo_GO=lapply(1:length(svd_d_convert), function(idx) {
  go_terms=names(trsf)[which(as.integer(strsplit(svd_d_convert[[idx]], split="")[[1]])==1)]
  data.frame(pattern_svd=idx, rank_go=go_terms)
}) %>%
  data.table::rbindlist() %>%
  mutate(rank=gsub("rank|GO.*","",rank_go), ID=gsub("rank[0-9]+","",rank_go)) %>%
  dplyr::select(-rank_go) %>%
  mutate(ID=gsub("\\.",":",ID), pattern_svd=paste0("pattern",pattern_svd))
# add annotation
svd_d_GOterms_annot=right_join(goenrich_res_df %>% dplyr::select(ID, Description) %>% subset(!duplicated(.)), 
                               base2_convertTo_GO) %>%
  subset(!duplicated(paste0(ID,"_",pattern_svd)))

### Plot
plot_all=
  ggplot(svd_d_GOterms_annot, aes(x=rank, y=Description, color=pattern_svd, size=pattern_svd)) +
  geom_point(shape=1) +
  scale_size_manual(values=c(1,3,5,7,9)) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_10")) +
  labs(x="rank", y=NULL, subtitle="\n") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        plot.margin=margin(t=10, l=5, unit="pt"),
        # legend.position="top",
        legend.direction="horizontal",
        legend.position=c(-1,1.075),
        legend.justification='left',
        # legend.box="horizontal"
  )
pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_PatternSpecificGOenrichments_all.pdf", height=4, width=7)
plot_all
dev.off()

plot_rm.pattern5=
  ggplot(svd_d_GOterms_annot %>% subset(pattern_svd!="pattern5"), 
         aes(x=rank, y=Description, color=pattern_svd, size=pattern_svd)) +
  geom_point(shape=1) +
  scale_size_manual(values=c(1,3,5,7)) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_10")) +
  labs(x="rank", y=NULL, subtitle="\n") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        # legend.position="top",
        legend.direction="horizontal",
        legend.position=c(-1,1.075), 
        legend.justification='left',
        # legend.box="horizontal",
        legend.margin=margin(0)
  )
pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_PatternSpecificGOenrichments_excludepattern5.pdf", height=4, width=7)
plot_rm.pattern5
dev.off()

plot_pattern5.only=
  ggplot(svd_d_GOterms_annot %>% subset(pattern_svd=="pattern5"), 
         aes(x=rank, y=Description, color=pattern_svd)) +
  geom_point(shape=1, size=9) +
  scale_color_manual(values=paletteer::paletteer_d("ggthemes::Classic_10")[5]) +
  labs(x="rank", y=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        legend.position="none",
        legend.title=element_blank(),
        legend.text=element_text(size=9),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        # legend.position="top",
        # legend.direction="horizontal",
        # legend.position=c(-1,1.075), 
        # legend.justification='left',
        # legend.box="horizontal",
        legend.margin=margin(0)
  )
pdf( "~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_PatternSpecificGOenrichments_pattern5only.pdf", height=4, width=7)
plot_pattern5.only
dev.off()

#####################################



### Analyze the timing (lagged correlation) in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)
library(Seurat)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.rough))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.rough==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Run CCF for each celltype
celltype_shortname=c("B","CD4T","CD8T","DC","Mono","NK","Other","OtherT")
for (i in 1:length(all_celltypes)) {
  
  ### Load the data for the celltype
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # remove the genes with two small expr or variation
  gene_expr_max=apply(pseudo_data, 1, max)
  pseudo_data=pseudo_data[gene_expr_max>0.2, ]
  pseudo_data.t=t(pseudo_data) %>% as.data.frame()
  
  ### Run lag correlation
  Lag_Cor=function(x, y) {
    corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
    corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
    max.cor_=max(corLag_clean)
    lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
    return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
  } 
  
  Lag_Cor_Res=list()
  for (n in 1:ncol(pseudo_data.t)) {
    Lag_Cor_Res[[n]]=lapply(1:ncol(pseudo_data.t), function(j) Lag_Cor(pseudo_data.t[[n]], pseudo_data.t[[j]]))
    names(Lag_Cor_Res[[n]])=colnames(pseudo_data.t)
    print(paste0(n," out of ",ncol(pseudo_data.t)," is done."))
  }
  names(Lag_Cor_Res)=colnames(pseudo_data.t)
  
  saveRDS(Lag_Cor_Res, paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],".rds"))
}

### Save the results as dataframes
for (i in 1:length(celltype_shortname)) {
  xxx=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],".rds"))
  gene1=rep(names(xxx), each=length(names(xxx)))
  gene2=lapply(xxx, function(x) names(x)) %>% Reduce(c, .)
  lag.time=lapply(xxx, function(x) lapply(x, function(y) y$s2.later.than.s1) %>% Reduce(c, .)) %>% Reduce(c, .)
  max.cor=lapply(xxx, function(x) lapply(x, function(y) y$max.correlation) %>% Reduce(c, .)) %>% Reduce(c, .)
  df_arranged=data.table::data.table(gene1=gene1, gene2=gene2, lag.time=lag.time, max.cor=max.cor)
  data.table::fwrite(df_arranged, paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],".txt.gz", sep="\t"))
}

# ### Remove the rds
# for (i in 1:length(celltype_shortname)) {
#   unlink(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],".rds"))
# }

#####################################



### Analyze the timing (lagged correlation) in B cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_B.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:20,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:20, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#     geom_line() + geom_point() +
#     scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
  Heatmap(df_arranged_lag.time_selected,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title=expression(Kendall~tau),
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun_cor,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=T, cluster_rows=T,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
          column_title_gp=gpar(fontsize=11),
          row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=10),
          row_split=paste0("timing\np.", pa$clustering),
          column_split=paste0("timing\np.", pa$clustering),
          height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
          width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
# p2>p1 ('>' means later)
  
### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(2, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="B cells", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
saveRDS(list(B=mean_of_the_cell), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2>p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=3))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_Bcells.pdf", width=2.05, height=1.2)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in CD4T cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_CD4T.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=3) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title=expression(Kendall~tau),
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        row_split=paste0("timing\np.", pa$clustering),
        column_split=paste0("timing\np.", pa$clustering),
        height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p2>p3>(=?)p1 ('>' means later)

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="CD4T cells", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(CD4T=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2>p3>(=?)p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,3,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=3))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_CD4Tcells.pdf", width=4, height=1.5)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in CD8T cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_CD8T.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=3) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title=expression(Kendall~tau),
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        row_split=paste0("timing\np.", pa$clustering),
        column_split=paste0("timing\np.", pa$clustering),
        height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p3>p2>p1 ('>' means later)

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="CD8T cells", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(CD8T=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p3>p2>p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2,3)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=2))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_CD8Tcells.pdf", width=3.7, height=1.3)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in DC cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_DC.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title="lags",
          legend_width=unit(16, "cm"), grid_height=unit(0.4, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=18), title_gp=gpar(fontsize=20, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=20), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=22),
        row_title_gp=gpar(fontsize=22),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=20),
        row_split=paste0("timing\npattern ", pa$clustering),
        column_split=paste0("timing\npattern ", pa$clustering),
        height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p1=p2 ('>' means later)
ht_cormtx=draw(ht_cormtx, heatmap_legend_side="top")

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="DCs", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(DC=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2=p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="ALL", OrgDb="org.Hs.eg.db")
saveRDS(enrich.compare, "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_DCs_enrichGO.rds")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=5))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10, lineheight=0.9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right", labels=function(x) stringr::str_wrap(x, width=30))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_DCs.pdf", width=3, height=3.2)
plot_
dev.off()

### Plot the gene-by-gene lags with fewer genes to reduce the pdf size
dim(df_arranged_max.cor)
cor_thr=0.85
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time (n.of.genes should be arround 133 to be consistent with the NK figure)
lag.time_thr=0.15
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist(), c(0.2,0.4,0.6,0.8,1))
col_fun_cor=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3"))
ht_cormtx=
  Heatmap(df_arranged_lag.time_selected,
          name=" ",
          show_heatmap_legend=TRUE,
          heatmap_legend_param=list(
            title="[DCs] lags",
            legend_width=unit(4.5, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=14), title_gp=gpar(fontsize=15, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=col_fun_cor,
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=T, cluster_rows=T,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=14), show_row_names=F, show_column_names=F,
          column_title_gp=gpar(fontsize=15),
          row_title_gp=gpar(fontsize=15),
          column_names_side="bottom",
          column_names_rot=90,
          column_names_gp=gpar(fontsize=14),
          row_split=paste0("timing\npattern ", pa$clustering),
          column_split=paste0("timing\npattern ", pa$clustering),
          height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
          width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
# p1=p2 ('>' means later)
ht_cormtx=draw(ht_cormtx, heatmap_legend_side="top")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPattern.AllLags_DCs.pdf", width=4, height=4.5)
ht_cormtx
dev.off()

#####################################



### Analyze the timing (lagged correlation) in Monocytes
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_Mono.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title=expression(Kendall~tau),
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        row_split=paste0("timing\np.", pa$clustering),
        column_split=paste0("timing\np.", pa$clustering),
        height=unit(1, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(1, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p2>p1 ('>' means later)

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="Monocytes", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(Monocyte=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2>p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=3))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10, lineheight=0.75),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right", labels=function(x) stringr::str_wrap(x, width=35))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_Monocytes.pdf", width=2.9, height=1.8)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in NK cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_NK.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3"))
ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title="[NK cells] lags",
          legend_width=unit(4.5, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=14), title_gp=gpar(fontsize=15, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=14), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=15),
        row_title_gp=gpar(fontsize=15),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=14),
        row_split=paste0("timing\npattern", pa$clustering),
        column_split=paste0("timing\npattern", pa$clustering),
        height=unit(0.6, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.6, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p2=p1 ('>' means later)
ht_cormtx=draw(ht_cormtx, heatmap_legend_side="top")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPattern.AllLags_NKs.pdf", width=4, height=4.5)
ht_cormtx
dev.off()

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="NK cells", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(NK=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2=p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db")
saveRDS(enrich.compare, "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_NKs_enrichGO.rds")
data_=(clusterProfiler::dotplot(enrich.compare, showCategory=5))$data

plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10, lineheight=0.9),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right", labels=function(x) stringr::str_wrap(x, width=30))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_NK.pdf", width=3, height=3.2)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in Other cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_Other.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=2) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
png("~/Project_PBMCage/tempt.png", width = 4000, height=4000, res = 100)
# ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title=expression(Kendall~tau),
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        row_split=paste0("timing\np.", pa$clustering),
        column_split=paste0("timing\np.", pa$clustering),
        height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
dev.off()
# p2>p1 ('>' means later)

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="NK cells", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(Other=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p2>p1 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(1,2)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db") # pvalueCutoff=0.05

data_=(clusterProfiler::dotplot(enrich.compare, showCategory=3))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10, lineheight=0.75),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right", labels=function(x) stringr::str_wrap(x, width=35))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_Othercells.pdf", width=2.8, height=1.2)
plot_
dev.off()

#####################################



### Analyze the timing (lagged correlation) in OtherT cells
#####################################
###

df_arranged=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_OtherT.txt.gz", sep=",")

df_arranged_lag.time=df_arranged %>% 
  dplyr::select(-max.cor) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="lag.time") %>%
  tibble::column_to_rownames("gene1")
df_arranged_max.cor=df_arranged %>% 
  dplyr::select(-lag.time) %>%
  tidyr::pivot_wider(names_from="gene2", values_from="max.cor") %>%
  tibble::column_to_rownames("gene1")

# remove too low correlation
dim(df_arranged_max.cor)
cor_thr=0.8
df_deselected=(apply(df_arranged_max.cor, 1, function(x) max(x[x!=1]))<cor_thr) | (apply(df_arranged_max.cor, 2, function(x) max(x[x!=1]))<cor_thr) # to reduce num of genes, increase the threshold
table(df_deselected)
df_arranged_selected=df_arranged_max.cor[!df_deselected, !df_deselected]
dim(df_arranged_selected)
# remove too many zero-lag-time
lag.time_thr=0.2
df_deselected2=(apply(df_arranged_lag.time, 1, function(x) length(which(x==0))>length(x)*lag.time_thr)) | (apply(df_arranged_lag.time, 2, function(x) length(which(x==0))>length(x)*lag.time_thr)) # to reduce num of genes, decrease the threshold
table(df_deselected|df_deselected2)
df_arranged_lag.time_selected=df_arranged_lag.time[!(df_deselected|df_deselected2), !(df_deselected|df_deselected2)]
dim(df_arranged_lag.time_selected) # better down to ~500 due to heatmap visualization

### Cluster
# tot_withinss=purrr::map_dbl(1:30,  function(k){
#   model=kmeans(x=df_arranged_lag.time_selected, centers=k)
#   model$tot.withinss
# })
# elbow_df=data.frame(k=1:30, tot_withinss=tot_withinss)
# nbclust_res=
#   ggplot(elbow_df, aes(x=k, y=tot_withinss)) +
#   geom_line() + geom_point() +
#   scale_x_continuous(breaks=1:30)
# slope=sapply(2:length(-diff(nbclust_res$data$tot_withinss)), function(idx) -diff(nbclust_res$data$tot_withinss)[idx-1]/-diff(nbclust_res$data$tot_withinss)[idx])
# k_best=which(slope<1)[1] # modify the threshold according to the plot
pa=cluster::pam(df_arranged_lag.time_selected, k=3) # use k_best or determine based on the final heatmap (staring with 3)
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
split(pa_clustering$gene, pa_clustering$cluster)

### Plot
quantile(as.vector(df_arranged_lag.time_selected) %>% unlist())
col_fun_cor=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
Heatmap(df_arranged_lag.time_selected,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title=expression(Kendall~tau),
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=F, show_column_dend=F,
        row_names_gp=gpar(fontsize=10), show_row_names=F, show_column_names=F,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        row_split=paste0("timing\np.", pa$clustering),
        column_split=paste0("timing\np.", pa$clustering),
        height=unit(0.35, "mm")*nrow(df_arranged_lag.time_selected), # 0.35mm
        width=unit(0.35, "mm")*ncol(df_arranged_lag.time_selected),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)
# p1>p2>p3 ('>' means later)

### Calculate the median lags between each pair of timing patterns
pa_clustering=data.frame(gene=names(pa$clustering), cluster=pa$clustering)
pa_clustering_split=split(pa_clustering$gene, pa_clustering$cluster)
mean_of_the_cell=matrix(NA, nrow=length(pa_clustering_split), ncol=length(pa_clustering_split))
rownames(mean_of_the_cell)=colnames(mean_of_the_cell)=paste0("T",names(pa_clustering_split))
for (i in 1:length(pa_clustering_split)) {
  for (j in 1:length(pa_clustering_split)) {
    gene_c1=pa_clustering_split[[i]]
    gene_c2=pa_clustering_split[[j]]
    mean_of_the_cell[i,j]=apply(df_arranged_lag.time_selected[gene_c1, gene_c2], c(1,2), as.numeric) %>% median()
  }
}
# plot the median lags
ht_median=
  Heatmap(mean_of_the_cell,
          name=" ",
          show_heatmap_legend=FALSE,
          heatmap_legend_param=list(
            title="Median lags",
            legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
            direction="horizontal", 
            title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
          ),
          border_gp=gpar(col="black", lty=2),
          col=circlize::colorRamp2(c(-5, 0, 5), c("dodgerblue3", "white", "brown3")),
          rect_gp=gpar(col="white", lwd=0.5),
          cluster_columns=F, cluster_rows=F,
          show_row_dend=F, show_column_dend=F,
          row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
          column_title_gp=gpar(fontsize=11),
          row_title="OtherT", row_title_gp=gpar(fontsize=11),
          column_names_side="bottom",
          column_names_rot=0,
          column_names_gp=gpar(fontsize=10), column_names_centered=T,
          height=unit(5, "mm")*nrow(mean_of_the_cell),
          width=unit(5, "mm")*ncol(mean_of_the_cell),
          row_dend_width=unit(2, "cm"),
          column_dend_height=unit(2, "cm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_median=draw(ht_median, heatmap_legend_side="top")

### Save the median lag
tempt=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
saveRDS(c(tempt, list(OtherT=mean_of_the_cell)), "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")

### Compare the timing patterns
# p1>p2>p3 ('>' means later)
t.patterns=split(pa_clustering$gene, pa_clustering$cluster)
t.patterns=t.patterns[c(3,2,1)]
enrich.compare=clusterProfiler::compareCluster(t.patterns, fun="enrichGO", keyType="SYMBOL", ont="BP", OrgDb="org.Hs.eg.db") # pvalueCutoff=0.05

data_=(clusterProfiler::dotplot(enrich.compare, showCategory=2))$data
plot_=
  ggplot(data_, aes(x=Cluster, y=Description)) +
  geom_point() +
  labs(title=NULL, y=NULL, x=NULL) +
  theme_minimal() +
  theme(axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=10, lineheight=0.75),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        title=element_text(size=11),
        plot.title.position="plot",
        plot.title=element_text(hjust=0.5),
        legend.position="right",
        legend.box="vertical",
        legend.direction="vertical",
        legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
        legend.spacing.y=unit(1,"mm"),
        legend.margin=margin(t=10, b=10)
  ) +
  scale_x_discrete(labels=paste0("T",names(t.patterns))) +
  scale_y_discrete(position="right", labels=function(x) stringr::str_wrap(x, width=35))

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatternEnrich_OtherT.pdf", width=3, height=1.2)
plot_
dev.off()

#####################################



### Plot median lags between each pair of timing patterns for rough celltypes
#####################################
### 

medianlags=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_meanLag_timingpatterns_rough.rds")
celltype_names=c("B cells","CD4T cells","CD8T cells","DCs","Monocytes","NK cells","Other cells","OtherT")

# plot the median lags
ht_medians=
  lapply(1:length(medianlags), function(idx) {
    Heatmap(medianlags[[idx]],
            name=" ",
            show_heatmap_legend=c(TRUE, rep(FALSE,7))[idx],
            heatmap_legend_param=list(
              title="Median lags",
              legend_width=unit(2.5, "cm"), grid_height=unit(0.2, "cm"), 
              direction="horizontal", 
              title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
            ),
            border_gp=gpar(col="black", lty=2),
            col=circlize::colorRamp2(c(-10, 0, 10), c("dodgerblue3", "white", "brown3")),
            rect_gp=gpar(col="white", lwd=0.5),
            cluster_columns=F, cluster_rows=F,
            show_row_dend=F, show_column_dend=F,
            row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
            column_title_gp=gpar(fontsize=11),
            row_title=celltype_names[idx], row_title_gp=gpar(fontsize=11),
            column_names_side="bottom",
            column_names_rot=0,
            column_names_gp=gpar(fontsize=10), column_names_centered=T,
            height=unit(6.5, "mm")*nrow(mean_of_the_cell),
            width=unit(6.5, "mm")*ncol(mean_of_the_cell),
            row_dend_width=unit(2, "cm"),
            column_dend_height=unit(2, "cm"))
  })
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")

### Save the plots
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_TimingPatterns_medianLags.pdf", width=1.5, height=1.5)
for (i in 1:length(ht_medians)) {
  draw(ht_medians[[i]], heatmap_legend_side="top")
}
dev.off()

#####################################




### Analyze the genes in each aging trend/cell type within each pattern
### 1) Same gene in different celltypes is grouped into different patterns
# ... i.e., same genes but in different patterns (because one gene has to be from diff. celltypes to be grouped into diff. pattern)
#####################################
###

### Load to analyze the same genes across aging trends
pattern_trend_genes=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes.rds")
# rearrange to have all the genes in each pattern
pattern_allgenes=lapply(pattern_trend_genes, function(x) Reduce(union, x))

### Plot upset plot of intersects
library(ComplexHeatmap)
library(ggplot2)
m=make_comb_mat(pattern_allgenes, mode="intersect")
m=normalize_comb_mat(m, full_comb_sets=T)
set_name(m) # check
# plot
upset_plot_=
  UpSet(m,
        top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black")),
        right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F),
        row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9))

# ### Enrich the intersects
# # ... DON'T RUN (because all the enrichments are the same)
# comb_name(m) # check
# pattern_intersect_genes=lapply(comb_name(m), function(x) extract_comb(m, x))
# goenrich_res=lapply(pattern_intersect_genes, function(x) clusterProfiler::enrichGO(x, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP"))
# comb_names=lapply(comb_name(m), function(x) paste0(set_name(m)[as.logical(as.integer(strsplit(x, split="")[[1]]))], collapse="_"))
# comb_names=unlist(comb_names)
# names(goenrich_res)=comb_names
# saveRDS(goenrich_res, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_IntersectsBetweenPattern.EnrichGO.rds")

### Plot upset plot of distincts
library(ComplexHeatmap)
library(ggplot2)
m=make_comb_mat(pattern_allgenes, mode="distinct")
m=normalize_comb_mat(m, full_comb_sets=T)
set_name(m) # check
comb_name(m) # check
# plot the upset
upset_plot_=
  UpSet(m,
        top_annotation=upset_top_annotation(m, add_numbers=TRUE, numbers_rot=90, annotation_name_rot=90, annotation_name_gp=gpar(fontsize=10), gp=gpar(col="black"), 
                                            bar_width=unit(2, "mm")),
        right_annotation=upset_right_annotation(m, add_numbers=TRUE, gp=gpar(fill="gray50"), show_annotation_name=F, axis_param=list(at=0, labels=c())),
        row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
        pt_size=unit(2, "mm"), lwd=0.5
        )

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_DinstinctsBetweenPattern.Upset.pdf", width=5.5, height=2.5)
upset_plot_
dev.off()

### Enrich the distincts (e.g., genes have both pattern1 and 2, but not show pattern 3,4,5)
pattern_distinct_genes=lapply(comb_name(m), function(x) extract_comb(m, x))
goenrich_res=lapply(pattern_distinct_genes, function(x) clusterProfiler::enrichGO(x, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP"))
comb_names=lapply(comb_name(m), function(x) paste0(set_name(m)[as.logical(as.integer(strsplit(x, split="")[[1]]))], collapse="_"))
comb_names=unlist(comb_names)
names(goenrich_res)=comb_names
# remove empty enrichments
goenrich_res_rm.null=goenrich_res[lapply(goenrich_res, function(x) !is.null(x)) %>% unlist() %>% unname()]
goenrich_res_rm.null_more=lapply(goenrich_res_rm.null, function(x) (subset(x@result, p.adjust<0.05) %>% nrow(.))!=0) %>%
  unlist() %>% unname()
goenrich_res_rm.null=goenrich_res_rm.null[goenrich_res_rm.null_more]
names(goenrich_res_rm.null) # check
saveRDS(goenrich_res_rm.null, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_DinstinctsBetweenPattern.EnrichGO.rds")
# dotplot
goenrich_res_rm.null=readRDS("~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_DinstinctsBetweenPattern.EnrichGO.rds")
tempts=lapply(goenrich_res_rm.null, function(x) clusterProfiler::dotplot(x, showCategory=5)) # take 5 to save space
names(goenrich_res_rm.null)
cowplot::plot_grid(plotlist=tempts[1:5]) # genes fallen into 5 or 4 different patterns
cowplot::plot_grid(plotlist=tempts[6:9])
cowplot::plot_grid(plotlist=tempts[10:15])
cowplot::plot_grid(plotlist=tempts[16:18])
# arrange the enrichment results with >10 enriched genes
goenrich_res_arranged=lapply(goenrich_res_rm.null, function(x) subset(x@result, p.adjust<0.05 & Count>10))
goenrich_res_arranged=lapply(1:length(goenrich_res_arranged), function(idx) {
  goenrich_res_arranged[[idx]] %>% mutate(distinct_among_patterns=names(goenrich_res_arranged)[idx])
}) %>% data.table::rbindlist()

### Plot the dotplots with genes exhibiting 4-5 patterns and pattern5-specific
# get the sorted data by clusterprofiler plotting
plot_4or5=tempts[c(1:5,18)]
names(plot_4or5)=paste0("pattern ",gsub("pattern","",gsub("_","\\/",names(plot_4or5))))
plot_data=lapply(plot_4or5, function(x) x$data)
# plot with ggplot2
MODULE_ENRICH_Plots=list()
for (i in 1:length(plot_data)) {
  module_=plot_data[[i]] %>% dplyr::select(Description, GeneRatio, p.adjust, Count)
  module_descrip_levels=levels(module_$Description)
  module_$Description=as.character(module_$Description)
  module_$Description=sapply(module_$Description, function(x) formatters::wrap_string(x, width=40, collapse="\n"))
  module_$Description=forcats::fct_relevel(module_$Description, formatters::wrap_string(module_descrip_levels, width=40, collapse="\n"))
  
  x_axis_expansion=mean(range(module_$GeneRatio, na.rm=T))/10
  
  plot_=
    ggplot(module_, aes(x=GeneRatio, y=Description, color=p.adjust, size=Count)) +
    geom_point() +
    scale_colour_gradient(limits=quantile(module_$p.adjust, na.rm=T)[c(1,5)],
                          low="brown4", high="pink",
                          labels=sprintf(fmt="%0.01e", quantile(module_$p.adjust, na.rm=T)[c(1,5)]),
                          breaks=quantile(module_$p.adjust, na.rm=T)[c(1,5)]) +
    scale_size_continuous(breaks=c(min(module_$Count, na.rm=T),
                                   round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                   max(module_$Count, na.rm=T)),
                          labels=c(min(module_$Count, na.rm=T),
                                   round((min(module_$Count, na.rm=T)+max(module_$Count, na.rm=T))/2),
                                   max(module_$Count, na.rm=T)),
                          range=c(2,6)) +
    guides(size=guide_legend(order=1)) +
    labs(title=names(plot_data)[i], y=NULL, x=NULL) +
    theme_minimal() +
    theme(axis.text.x=element_text(size=9, angle=90, hjust=1, vjust=0.5),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=10),
          axis.title.y=element_text(size=10),
          title=element_text(size=11),
          plot.title.position="plot",
          plot.title=element_text(hjust=0.5),
          legend.position="right",
          legend.box="vertical",
          legend.direction="vertical",
          legend.key.width=unit(0.5,"cm"), legend.key.height=unit(0.4,"cm"),
          legend.spacing.y=unit(1,"mm"),
          legend.margin=margin(t=10, b=10)
    ) +
    guides(color=guide_colorbar(title="p.adj",
                                label.theme=element_text(angle=90, size=9, vjust=0.5, hjust=0.5),
                                title.theme=element_text(angle=90, size=10, vjust=1, hjust=0), order=1),
           size=guide_legend(label.theme=element_text(size=9),
                             title.theme=element_text(size=10), order=2)) +
    scale_x_continuous(expand=c(x_axis_expansion, x_axis_expansion))
  
  MODULE_ENRICH_Plots=c(MODULE_ENRICH_Plots, list(plot_))
}
plot_all=
  cowplot::plot_grid(plotlist=MODULE_ENRICH_Plots, align="hv", nrow=3)
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Rough_bothSex_PatternGenes_DinstinctsBetweenPattern.EnrichGO_distinctAcross4or5patterns.pdf", width=9, height=7.5)
plot_all
dev.off()

#####################################



### Analyze the timing (lagged correlation) across celltypes based on pseudobulk_rough
#####################################
###

library(dplyr)
library(Seurat)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.rough))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.rough==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Annotate each gene with its celltype origin and merge the expr data
pseudo_data.arranged=list()
for (i in 1:length(all_celltypes)) {
  
  ### Load the data for the celltype
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # remove the genes with two small expr or variation
  gene_expr_max=apply(pseudo_data, 1, max)
  pseudo_data=pseudo_data[gene_expr_max>quantile(gene_expr_max, 0.975), ]
  
  pseudo_data.arranged[[i]]=pseudo_data %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("celltype.gene") %>%
    mutate(celltype.gene=paste0(all_celltypes[i],"|||",celltype.gene))
}
pseudo_data.merged=data.table::rbindlist(pseudo_data.arranged, fill=TRUE) %>%
  na.omit()

### Run CCF
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 

Lag_Cor_Res=list()
pseudo_data.t=pseudo_data.merged %>% tibble::column_to_rownames("celltype.gene") %>%
  t() %>% as.data.frame()
for (n in 1:ncol(pseudo_data.t)) {
  Lag_Cor_Res[[n]]=lapply(1:ncol(pseudo_data.t), function(j) Lag_Cor(pseudo_data.t[[n]], pseudo_data.t[[j]]))
  names(Lag_Cor_Res[[n]])=colnames(pseudo_data.t)
  print(paste0(n," out of ",ncol(pseudo_data.t)," is done."))
}
names(Lag_Cor_Res)=colnames(pseudo_data.t)
saveRDS(Lag_Cor_Res, "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_all.rds")

### Save the results as dataframes
xxx=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_all.rds")
gene1=rep(names(xxx), each=length(names(xxx)))
gene2=lapply(xxx, function(x) names(x)) %>% Reduce(c, .)
lag.time=lapply(xxx, function(x) lapply(x, function(y) y$s2.later.than.s1) %>% Reduce(c, .)) %>% Reduce(c, .)
max.cor=lapply(xxx, function(x) lapply(x, function(y) y$max.correlation) %>% Reduce(c, .)) %>% Reduce(c, .)
df_arranged=data.table::data.table(gene1=gene1, gene2=gene2, lag.time=lag.time, max.cor=max.cor)
data.table::fwrite(df_arranged, "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_all.txt.gz", sep="\t")

# ### Remove the rds
# unlink("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_all.rds")

#####################################






### Cluster the variation of genes across ages in pseudobulk.inter, regardless of sex
#####################################
###

library(dplyr)
library(Seurat)
library(ClusterGVis)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.inter))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.inter==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Run ClusterGVis
CKs=list()
for (i in 1:length(all_celltypes)) {
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # remove the genes with two small expr or variation
  gene_expr_max=apply(pseudo_data, 1, max)
  pseudo_data=pseudo_data[gene_expr_max>0.1, ]
  # run getClusters to determine the best cluster.num
  nbclust_res=factoextra::fviz_nbclust(pseudo_data, FUNcluster=kmeans, method="wss", k.max=10) # check the optimal number of clusters
  slope=sapply(2:length(-diff(nbclust_res$data$y)), function(idx) -diff(nbclust_res$data$y)[idx-1]/-diff(nbclust_res$data$y)[idx])
  k_best=which(slope<2)[1]
  ck_save=function(data) {
    library(ClusterGVis)
    tryCatch({
      ck_=clusterData(exp=data,
                      cluster.method="kmeans",
                      cluster.num=k_best)
      return(ck_)
    }, error=function(msg) {
      print("Oversize.")
      return(NA)
    })
  }
  CKs[[i]]=ck_save(pseudo_data)
}
names(CKs)=all_celltypes
saveRDS(CKs, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Inter_bothSex.rds")

#####################################



### Cluster the variation of genes across ages in pseudobulk.detailed, regardless of sex
#####################################
###

library(dplyr)
library(Seurat)
library(ClusterGVis)

### Load the data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.detailed))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.detailed==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Run ClusterGVis
CKs=list()
for (i in 1:length(all_celltypes)) {
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # remove the genes with two small expr or variation
  gene_expr_max=apply(pseudo_data, 1, max)
  pseudo_data=pseudo_data[gene_expr_max>0.1, ]
  # run getClusters to determine the best cluster.num
  nbclust_res=factoextra::fviz_nbclust(pseudo_data, FUNcluster=kmeans, method="wss", k.max=10) # check the optimal number of clusters
  slope=sapply(2:length(-diff(nbclust_res$data$y)), function(idx) -diff(nbclust_res$data$y)[idx-1]/-diff(nbclust_res$data$y)[idx])
  k_best=which(slope<2)[1]
  ck_save=function(data) {
    library(ClusterGVis)
    tryCatch({
      ck_=clusterData(exp=data,
                      cluster.method="kmeans",
                      cluster.num=k_best)
      return(ck_)
    }, error=function(msg) {
      print("Oversize.")
      return(NA)
    })
  }
  CKs[[i]]=ck_save(pseudo_data)
}
names(CKs)=all_celltypes
saveRDS(CKs, "~/Project_PBMCage/Tempt_RDS/ClusterGVis_Allgenes_Clustering.BasedOn.Pseudobulk.Detailed_bothSex.rds")

#####################################



### Analyze the timing (lagged correlation) between cell proportion at the rough, inter, and detailed levels respectively
#####################################
###

library(dplyr)

### Load the cell proportion (use the monocyte-modified version to be consistent with the annotation to be published)
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values_updated.txt", sep="\t")
  
freq_rough=meta_info_combined %>% 
  group_by(Annot.rough, age) %>% summarize_at("percent.rough", mean) %>%
  tidyr::pivot_wider(id_cols="age", names_from="Annot.rough", values_from="percent.rough") %>%
  tibble::column_to_rownames("age")
freq_rough[is.na(freq_rough)]=0

freq_inter=meta_info_combined %>% 
  group_by(Annot.inter, age) %>% summarize_at("percent.inter", mean) %>%
  tidyr::pivot_wider(id_cols="age", names_from="Annot.inter", values_from="percent.inter") %>%
  tibble::column_to_rownames("age")
freq_inter[is.na(freq_inter)]=0

freq_detailed=meta_info_combined %>% 
  group_by(Annot.detailed, age) %>% summarize_at("percent.detailed", mean) %>%
  tidyr::pivot_wider(id_cols="age", names_from="Annot.detailed", values_from="percent.detailed") %>%
  tibble::column_to_rownames("age")
freq_detailed[is.na(freq_detailed)]=0

### Run CCF for each level
freqs_all=list(freq_rough, freq_inter, freq_detailed)
DFs_list=list()
for (i in 1:length(freqs_all)) {
  
  ### Run lag correlation
  Lag_Cor=function(x, y) {
    corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
    corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
    max.cor_=max(corLag_clean)
    lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
    return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
  } 
  Lag_Cor_Res=lapply(1:ncol(freqs_all[[i]]), function(j) {
    tempt_=apply(freqs_all[[i]], 2, function(celltype) Lag_Cor(celltype, freqs_all[[i]][[j]]))
    names(tempt_)=colnames(freqs_all[[i]])
    tempt_
  })
  names(Lag_Cor_Res)=colnames(freqs_all[[i]])
  
  ### Arrange the df
  celltype1=rep(colnames(freqs_all[[i]]), each=ncol(freqs_all[[i]])*2)
  celltype2=rep(rep(colnames(freqs_all[[i]]), each=2), times=ncol(freqs_all[[i]]))
  parameter=rep(c("s2.later.than.s1","max.correlation"), times=ncol(freqs_all[[i]])*ncol(freqs_all[[i]]))
  df_arranged=data.frame(celltype1=celltype1, celltype2=celltype2, parameter=parameter, value=unlist(Lag_Cor_Res)) %>%
    tibble::remove_rownames() %>%
    tidyr::pivot_wider(id_cols=c("celltype1","celltype2"), names_from="parameter", values_from="value")
  
  ### Save all the 3 levels
  DFs_list[[i]]=df_arranged
}
names(DFs_list)=c("Annot.rough","Annot.inter","Annot.detailed")
saveRDS(DFs_list, "~/Project_PBMCage/Tempt_RDS/Timing.Cellprop_All3Levels.rds")

### Plot the rough
DFs_list=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.Cellprop_All3Levels.rds")
cor_mtx=DFs_list[[1]] %>% dplyr::select(celltype1, celltype2, max.correlation) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="max.correlation") %>%
  tibble::column_to_rownames("celltype1")
lag_mtx=DFs_list[[1]] %>% dplyr::select(celltype1, celltype2, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("celltype1")
range(lag_mtx)
# rownames as celltype1, colnames as celltype2, so las refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags",
            legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            grid.rect(x=x, y=y, width=width, height=height,
                      gp=gpar(col="grey50", fill=NA, lwd=0.1))
            if (i!=j) {
              grid.circle(x=x, y=y, r=ifelse(cor_mtx[i, j]>0, cor_mtx[i, j], 0)/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(6,"mm")*nrow(lag_mtx),
          width=unit(6,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="top")


### Plot the inter level
DFs_list=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.Cellprop_All3Levels.rds")
cor_mtx=DFs_list[[2]] %>% dplyr::select(celltype1, celltype2, max.correlation) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="max.correlation") %>%
  tibble::column_to_rownames("celltype1")
lag_mtx=DFs_list[[2]] %>% dplyr::select(celltype1, celltype2, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("celltype1")
# remove Other.Eryth and Other.Mast
cor_mtx=cor_mtx[!grepl("Eryth|Mast",rownames(cor_mtx)),!grepl("Eryth|Mast",colnames(cor_mtx))]
lag_mtx=lag_mtx[!grepl("Eryth|Mast",rownames(lag_mtx)),!grepl("Eryth|Mast",colnames(lag_mtx))]
range(lag_mtx)
# rownames as celltype1, colnames as celltype2, so las refer to colnames later than rownames

col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_inter=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              # highlight the CD4T.naive & CD8T.naive on the figure
              if (i==which(rownames(lag_mtx)=="CD4T.naive") & j==which(rownames(lag_mtx)=="CD8T.naive")) {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="black", fill=NA, lwd=1, lty=2))
              } else {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="grey80", fill=NA, lwd=0.1))
              }
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor and not at the diagnal
            if (i!=j & cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_inter=draw(ht_inter, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing.Cellprop_inter.pdf", width=6.6, height=6.1)
ht_inter
dev.off()

### Plot the detailed level
DFs_list=readRDS("~/Project_PBMCage/Tempt_RDS/Timing.Cellprop_All3Levels.rds")
cor_mtx=DFs_list[[3]] %>% dplyr::select(celltype1, celltype2, max.correlation) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="max.correlation") %>%
  tibble::column_to_rownames("celltype1")
lag_mtx=DFs_list[[3]] %>% dplyr::select(celltype1, celltype2, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="celltype1", names_from="celltype2", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("celltype1")
range(lag_mtx)
# rownames as celltype1, colnames as celltype2, so las refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags",
            legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            grid.rect(x=x, y=y, width=width, height=height,
                      gp=gpar(col="grey50", fill=NA, lwd=0.1))
            if (i!=j) {
              grid.circle(x=x, y=y, r=ifelse(cor_mtx[i, j]>0, cor_mtx[i, j], 0)/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(6,"mm")*nrow(lag_mtx),
          width=unit(6,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="top")

#####################################









##########################################################################
###
### Lag correlation between cell proportions and the mean expr of all genes
###
##########################################################################



### Analyze the timing (lagged correlation) between cell proportion and gene expression in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)
library(Seurat)

### Load the gene expression data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_rough.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.rough))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.rough==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Load the cell proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% 
  group_by(Annot.rough, age) %>% summarize_at("percent.rough", mean) %>%
  split(.$Annot.rough) %>% 
  .[all_celltypes]
  
### Run CCF for each celltype
celltype_shortname=c("B","CD4T","CD8T","DC","Mono","NK","Other","OtherT")
for (i in 1:length(all_celltypes)) {
  
  ### Load the data for the celltype
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # # remove the genes with two small expr or variation
  # gene_expr_max=apply(pseudo_data, 1, max)
  # pseudo_data=pseudo_data[gene_expr_max>0.2, ]
  pseudo_data.t=t(pseudo_data) %>% as.data.frame() %>% tibble::rownames_to_column("age") %>%
    mutate(age=as.integer(age)) %>%
    left_join(., meta_info_combined[[i]], by="age")
  expr.data=pseudo_data.t %>% .[,!grepl("age|percent\\.rough|Annot\\.rough",colnames(.))]
  percent.data=pseudo_data.t[["percent.rough"]]
  
  ### Run lag correlation
  Lag_Cor=function(x, y) {
    corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
    corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
    max.cor_=max(corLag_clean)
    lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
    return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
  } 
  
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, percent.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  print(paste0(i, " out of ", length(all_celltypes), " is done."))
  
  saveRDS(Lag_Cor_Res, paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

### Save the results as dataframes
df_arranged=list()
for (i in 1:length(celltype_shortname)) {
  xxx=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
  genes=names(xxx)
  lag.time=lapply(xxx, function(y) y$s2.later.than.s1) %>% Reduce(c, .)
  max.cor=lapply(xxx, function(y) y$max.correlation) %>% Reduce(c, .)
  df_arranged[[i]]=data.table::data.table(gene.vs.cellprop=genes, lag.time=lag.time, max.cor=max.cor, Annot.rough=all_celltypes[i])
}
data.table::fwrite(data.table::rbindlist(df_arranged), 
                   "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_genesVScellprop_all.txt.gz", 
                   sep="\t")
tempt_=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_genesVScellprop_all.txt.gz", sep="\t")

### Remove the rds
for (i in 1:length(celltype_shortname)) {
  unlink(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

#####################################



### Analyze the timing (lagged correlation) between cell proportion and gene expression in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)
library(Seurat)

### Load the gene expression data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_inter.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.inter))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.inter==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Load the cell proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% 
  group_by(Annot.inter, age) %>% summarize_at("percent.inter", mean) %>%
  split(.$Annot.inter) %>% 
  .[all_celltypes]

### Run CCF for each celltype
celltype_shortname=gsub("\\.| ","",all_celltypes)
for (i in 1:length(all_celltypes)) {
  
  ### Load the data for the celltype
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # # remove the genes with two small expr or variation
  # gene_expr_max=apply(pseudo_data, 1, max)
  # pseudo_data=pseudo_data[gene_expr_max>0.2, ]
  pseudo_data.t=t(pseudo_data) %>% as.data.frame() %>% tibble::rownames_to_column("age") %>%
    mutate(age=as.integer(age)) %>%
    left_join(., meta_info_combined[[i]], by="age")
  expr.data=pseudo_data.t %>% .[,!grepl("age|percent\\.inter|Annot\\.inter",colnames(.))]
  percent.data=pseudo_data.t[["percent.inter"]]
  
  ### Run lag correlation
  Lag_Cor=function(x, y) {
    corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
    corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
    max.cor_=max(corLag_clean)
    lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
    return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
  } 
  
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, percent.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  print(paste0(i, " out of ", length(all_celltypes), " is done."))
  
  saveRDS(Lag_Cor_Res, paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

### Save the results as dataframes
df_arranged=list()
for (i in 1:length(celltype_shortname)) {
  xxx=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
  genes=names(xxx)
  lag.time=lapply(xxx, function(y) y$s2.later.than.s1) %>% Reduce(c, .)
  max.cor=lapply(xxx, function(y) y$max.correlation) %>% Reduce(c, .)
  df_arranged[[i]]=data.table::data.table(gene.vs.cellprop=genes, lag.time=lag.time, max.cor=max.cor, Annot.inter=all_celltypes[i])
}
data.table::fwrite(data.table::rbindlist(df_arranged), 
                   "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_genesVScellprop_all.txt.gz", 
                   sep="\t")
tempt_=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_genesVScellprop_all.txt.gz", sep="\t")

### Remove the rds
for (i in 1:length(celltype_shortname)) {
  unlink(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

#####################################



### Analyze the timing (lagged correlation) between cell proportion and gene expression in each celltype based on pseudobulk_detailed
#####################################
###

library(dplyr)
library(Seurat)

### Load the gene expression data
pseudobulkobj_all=readRDS("~/Project_PBMCage/raw_data/PBMCage_annotated_AssignedOnly_Pseudobulk_detailed.rds")
all_celltypes=names(table(pseudobulkobj_all$Annot.detailed))
pseudobulkobj_all_subsets=lapply(all_celltypes, function(celltype_) subset(pseudobulkobj_all, Annot.detailed==celltype_))
pseudobulkobj_subsets_ageMerged=lapply(pseudobulkobj_all_subsets, function(x) Seurat::AggregateExpression(x, group.by="age", return.seurat=TRUE))
names(pseudobulkobj_subsets_ageMerged)=all_celltypes

### Load the cell proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
meta_info_combined=meta_info_combined %>% 
  group_by(Annot.detailed, age) %>% summarize_at("percent.detailed", mean) %>%
  split(.$Annot.detailed) %>% 
  .[all_celltypes]

### Run CCF for each celltype
celltype_shortname=gsub("\\.| ","",all_celltypes)
for (i in 1:length(all_celltypes)) {
  
  ### Load the data for the celltype
  pseudo_data=LayerData(pseudobulkobj_subsets_ageMerged[[i]], layer="data")
  # # remove the genes with two small expr or variation
  # gene_expr_max=apply(pseudo_data, 1, max)
  # pseudo_data=pseudo_data[gene_expr_max>0.2, ]
  pseudo_data.t=t(pseudo_data) %>% as.data.frame() %>% tibble::rownames_to_column("age") %>%
    mutate(age=as.integer(age)) %>%
    left_join(., meta_info_combined[[i]], by="age")
  expr.data=pseudo_data.t %>% .[,!grepl("age|percent\\.detailed|Annot\\.detailed",colnames(.))]
  percent.data=pseudo_data.t[["percent.detailed"]]
  
  ### Run lag correlation
  Lag_Cor=function(x, y) {
    corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
    corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
    max.cor_=max(corLag_clean)
    lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
    return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
  } 
  
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, percent.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  print(paste0(i, " out of ", length(all_celltypes), " is done."))
  
  saveRDS(Lag_Cor_Res, paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

### Save the results as dataframes
df_arranged=list()
for (i in 1:length(celltype_shortname)) {
  xxx=readRDS(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
  genes=names(xxx)
  lag.time=lapply(xxx, function(y) y$s2.later.than.s1) %>% Reduce(c, .)
  max.cor=lapply(xxx, function(y) y$max.correlation) %>% Reduce(c, .)
  df_arranged[[i]]=data.table::data.table(gene.vs.cellprop=genes, lag.time=lag.time, max.cor=max.cor, Annot.detailed=all_celltypes[i])
}
data.table::fwrite(data.table::rbindlist(df_arranged), 
                   "~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_genesVScellprop_all.txt.gz", 
                   sep="\t")
tempt_=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_genesVScellprop_all.txt.gz", sep="\t")

### Remove the rds
for (i in 1:length(celltype_shortname)) {
  unlink(paste0("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_",celltype_shortname[i],"_genesVScellprop.rds"))
}

#####################################



### Plot the timing (lagged correlation) between cell proportion and gene expression in each celltype based on pseudobulk_rough, inter, detailed
#####################################
###

### Load the timing data of rough
timing_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_genesVScellprop_all.txt.gz", sep="\t")
timing_data=timing_data %>% split(.$Annot.rough)

### Analyze per cell type at rough
n_all=c()
for (i in 1:length(timing_data)) {
  # sort only correlation>0
  timing_data_celltypes=timing_data[[i]] %>%
    dplyr::select(gene.vs.cellprop, lag.time, max.cor) %>%
    tibble::column_to_rownames("gene.vs.cellprop") %>%
    subset(max.cor>0)
  # count the lag.time that is pos, neg, and zero
  n.pos=timing_data_celltypes %>% subset(lag.time>0) %>% nrow()
  n.neg=timing_data_celltypes %>% subset(lag.time<0) %>% nrow()
  n.zero=timing_data_celltypes %>% subset(lag.time==0) %>% nrow()
  # summarize
  n_all[[i]]=
    data.frame(pos.lags=n.pos/nrow(timing_data_celltypes),
               neg.lags=n.neg/nrow(timing_data_celltypes),
               zero.lags=n.zero/nrow(timing_data_celltypes),
               celltype=names(timing_data)[i])
}
lags.counts_rough=data.table::rbindlist(n_all) %>% mutate(celltype.level="Annot.rough")
matrix_df=
  lags.counts_rough %>% 
  dplyr::select(-celltype.level) %>% tibble::column_to_rownames("celltype") %>% na.omit()

col_fun_cor=circlize::colorRamp2(c(0,1), c("white", "brown3"))
# ht_cormtx=
Heatmap(matrix_df,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title="Percent",
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=T, show_column_dend=T,
        row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        height=unit(6, "mm")*nrow(matrix_df), # 0.35mm
        width=unit(6, "mm")*ncol(matrix_df)
)

### Analyze all the Annot.rough proportions' lags from genes
timing_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Rough_bothSex_genesVScellprop_all.txt.gz", sep="\t")
cor.data=timing_data %>%
  tidyr::pivot_wider(id_cols="gene.vs.cellprop", names_from="Annot.rough", values_from="max.cor") %>%
  tibble::column_to_rownames("gene.vs.cellprop") %>%
  na.omit()
lag.data=timing_data %>%
  tidyr::pivot_wider(id_cols="gene.vs.cellprop", names_from="Annot.rough", values_from="lag.time") %>%
  tibble::column_to_rownames("gene.vs.cellprop") %>% 
  .[rownames(cor.data), ]
# enrich the genes with (at least 1 lag>0) and (all cor>=0.1)
lag_keep=lag.data %>% apply(., 1, max) %>% .[.>0] # select the genes earlier than prop
cor_keep=cor.data %>% apply(., 1, min) %>% .[.>=0.1] # select the genes correlated than prop
both_cor_lag_keep=intersect(names(lag_keep), names(cor_keep))
enrich_res=clusterProfiler::enrichGO(both_cor_lag_keep, OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
# draw heatmap
lag_keep_less=lag.data %>% apply(., 1, max) %>% .[.>0] # to reduce no of genes for plotting
cor_keep_less=cor.data %>% apply(., 1, min) %>% .[.>=0.2] # to reduce no of genes for plotting
both_cor_lag_keep_less=intersect(names(lag_keep_less), names(cor_keep_less))
lag.data_selected=lag.data[both_cor_lag_keep_less,]
cor.data_selected=cor.data[both_cor_lag_keep_less,]
lag.data_selected.t=t(lag.data_selected)
marking_genes=enrich_res@result %>% .[c(1:8),] %>% dplyr::select(geneID) %>% 
  tibble::deframe() %>% 
  lapply(., function(x) strsplit(x, split="/")[[1]]) %>% Reduce(union,.) %>% 
  na.omit() %>%
  intersect(., both_cor_lag_keep_less)
marking_indices=sapply(marking_genes, function(x) which(both_cor_lag_keep_less==x)) %>% unname()
topanno=HeatmapAnnotation(foo=anno_mark(at=marking_indices,
                                        labels=marking_genes,
                                        which="column", side="left"),
                                   gp=gpar(fontsize=3))
range(lag.data_selected.t)
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag.data_selected.t, 
          name=" ",
          heatmap_legend_param=list(
            title="lags",
            legend_width=unit(5, "cm"), grid_height=unit(0.2, "cm"),
            direction="horizontal", title_position="topcenter", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun,
          rect_gp=gpar(col="white", lwd=0.1),
          bottom_annotation=topanno,
          cluster_rows=T, cluster_columns=T,
          show_row_names=T, show_column_names=F,
          show_row_dend=T, show_column_dend=F,
          row_names_side="right", row_names_gp=gpar(fontsize=10),
          height=unit(5, "mm")*nrow(lag.data_selected.t),
          width=unit(0.5, "mm")*ncol(lag.data_selected.t)
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="top")
# pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing.Genes_VS_Cellprop_rough_heatmap.pdf", width=20, height=3.5)
# ht_rough
# dev.off()


### Load the timing data of inter
timing_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_genesVScellprop_all.txt.gz", sep="\t")
timing_data=timing_data %>% split(.$Annot.inter)

### Analyze per cell type at inter (numbers of lag<0, lag==0, lag>0)
n_all=c()
for (i in 1:length(timing_data)) {
  # sort only correlation>0
  timing_data_celltypes=timing_data[[i]] %>%
    dplyr::select(gene.vs.cellprop, lag.time, max.cor) %>%
    tibble::column_to_rownames("gene.vs.cellprop") %>%
    subset(max.cor>quantile(max.cor, 0.75, na.rm=T)[1])
  # count the lag.time that is pos, neg, and zero
  n.pos=timing_data_celltypes %>% subset(lag.time>0) %>% nrow()
  n.neg=timing_data_celltypes %>% subset(lag.time<0) %>% nrow()
  n.zero=timing_data_celltypes %>% subset(lag.time==0) %>% nrow()
  # summarize
  n_all[[i]]=
    data.frame(pos.lags=n.pos/nrow(timing_data_celltypes),
               neg.lags=n.neg/nrow(timing_data_celltypes),
               zero.lags=n.zero/nrow(timing_data_celltypes),
               celltype=names(timing_data)[i])
}
lags.counts_inter=data.table::rbindlist(n_all) %>% mutate(celltype.level="Annot.inter")
matrix_df=
  lags.counts_inter %>% 
  dplyr::select(-celltype.level) %>% tibble::column_to_rownames("celltype") %>% na.omit()

col_fun_cor=circlize::colorRamp2(c(0,1), c("white", "brown3"))
# ht_cormtx=
Heatmap(matrix_df,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title="Percent",
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=T, show_column_dend=T,
        row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        height=unit(6, "mm")*nrow(matrix_df), # 0.35mm
        width=unit(6, "mm")*ncol(matrix_df)
)

### Analyze all the Annot.inter proportions' lags from genes
timing_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Inter_bothSex_genesVScellprop_all.txt.gz", sep="\t")
Annot.inter_selected=unique(timing_data$Annot.inter) %>% .[grepl("CD4T|CD8T|NK\\.CD56dim|MAIT|B\\.naive|Mono\\.nonclassical|OtherT\\.NKT",.)]
cor.data=timing_data %>%
  subset(Annot.inter %in% Annot.inter_selected) %>%
  tidyr::pivot_wider(id_cols="gene.vs.cellprop", names_from="Annot.inter", values_from="max.cor") %>%
  tibble::column_to_rownames("gene.vs.cellprop") %>%
  na.omit()
lag.data=timing_data %>%
  subset(Annot.inter %in% Annot.inter_selected) %>%
  tidyr::pivot_wider(id_cols="gene.vs.cellprop", names_from="Annot.inter", values_from="lag.time") %>%
  tibble::column_to_rownames("gene.vs.cellprop") %>% 
  .[rownames(cor.data), ]
# enrich the genes with (at least 1 lag>0) and (all cor>=0.1)
lag_keep=lag.data %>% apply(., 1, max) %>% .[.>0] # select the genes earlier than prop
# cor_keep=cor.data %>% apply(., 1, min) %>% .[.>=0.1] # select the genes correlated than prop
both_cor_lag_keep=intersect(names(lag_keep), names(cor_keep))
row_dend=as.dendrogram(hclust(dist(t(lag.data[names(lag_keep),]))))
enrich_res=clusterProfiler::enrichGO(names(lag_keep), OrgDb="org.Hs.eg.db", keyType="SYMBOL", ont="BP")
# draw heatmap
lag_keep_less=lag.data %>% apply(., 1, max) %>% .[.>0] # to reduce no of genes for plotting
cor_keep_less=cor.data %>% apply(., 1, min) %>% .[.>=0.15] # to reduce no of genes for plotting
both_cor_lag_keep_less=intersect(names(lag_keep_less), names(cor_keep_less))
lag.data_selected=lag.data[both_cor_lag_keep_less,]
cor.data_selected=cor.data[both_cor_lag_keep_less,]
lag.data_selected.t=t(lag.data_selected)
per_description=enrich_res@result %>% subset(p.adjust<0.05) %>% .$geneID; names(per_description)=enrich_res@result %>% subset(p.adjust<0.05) %>% .$Description
no_of_marking_genes_inplot=lapply(per_description, function(x) strsplit(x, split="/")[[1]] %>% intersect(., both_cor_lag_keep_less) %>% length(.)) %>% 
  unlist()
View(no_of_marking_genes_inplot)
marking_genes=enrich_res@result %>% subset(p.adjust<0.05) %>% 
  subset(Description %in% names(no_of_marking_genes_inplot)[1:30]) %>%
  dplyr::select(geneID) %>% 
  tibble::deframe() %>% 
  lapply(., function(x) strsplit(x, split="/")[[1]]) %>% Reduce(union,.) %>% 
  na.omit() %>%
  intersect(., both_cor_lag_keep_less)
marking_indices=sapply(marking_genes, function(x) which(both_cor_lag_keep_less==x)) %>% unname()
topanno=HeatmapAnnotation(foo=anno_mark(at=marking_indices,
                                        labels=marking_genes,
                                        which="column", side="left",
                                        labels_gp=gpar(fontsize=9))
                          )
range(lag.data_selected.t)
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_inter=
  Heatmap(lag.data_selected.t, 
          name=" ",
          heatmap_legend_param=list(
            title="lags",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun,
          # rect_gp=gpar(col="white", lwd=3.4),
          cell_fun=function(j, i, x, y, w, h, fill) {
            grid.lines(x=c(x-w/2, x+w/2), y=c(y+h/2, y+h/2), gp=gpar(col='white', lwd=6))
          },
          bottom_annotation=topanno,
          cluster_columns=T, cluster_rows=row_dend,
          show_row_names=T, show_column_names=F, 
          show_row_dend=T, show_column_dend=F,
          row_names_side="right", row_names_gp=gpar(fontsize=10),
          height=unit(4, "mm")*nrow(lag.data_selected.t),
          width=unit(0.55, "mm")*ncol(lag.data_selected.t),
          row_dend_width=unit(3, "mm")
  )
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_inter=draw(ht_inter, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing.Genes_VS_Cellprop_inter_heatmap.pdf", width=6.6, height=3.4)
ht_inter
dev.off()


### Load the timing data of detailed
timing_data=data.table::fread("~/Project_PBMCage/Tempt_RDS/Timing.BasedOn.Pseudobulk.Detailed_bothSex_genesVScellprop_all.txt.gz", sep="\t")
timing_data=timing_data %>% split(.$Annot.detailed)

### Analyze per cell type at detailed
n_all=c()
for (i in 1:length(timing_data)) {
  # sort only correlation>0
  timing_data_celltypes=timing_data[[i]] %>%
    dplyr::select(gene.vs.cellprop, lag.time, max.cor) %>%
    tibble::column_to_rownames("gene.vs.cellprop") %>%
    subset(max.cor>quantile(max.cor, 0.75, na.rm=T)[1])
  # count the lag.time that is pos, neg, and zero
  n.pos=timing_data_celltypes %>% subset(lag.time>0) %>% nrow()
  n.neg=timing_data_celltypes %>% subset(lag.time<0) %>% nrow()
  n.zero=timing_data_celltypes %>% subset(lag.time==0) %>% nrow()
  # summarize
  n_all[[i]]=
    data.frame(pos.lags=n.pos/nrow(timing_data_celltypes),
               neg.lags=n.neg/nrow(timing_data_celltypes),
               zero.lags=n.zero/nrow(timing_data_celltypes),
               celltype=names(timing_data)[i])
}
lags.counts_detailed=data.table::rbindlist(n_all) %>% mutate(celltype.level="Annot.detailed")
matrix_df=
  lags.counts_detailed %>% 
  dplyr::select(-celltype.level) %>% tibble::column_to_rownames("celltype") %>% na.omit()

col_fun_cor=circlize::colorRamp2(c(0,0.5,1), c("dodgerblue3", "white", "brown3"))
# ht_cormtx=
Heatmap(matrix_df,
        name=" ",
        show_heatmap_legend=TRUE,
        heatmap_legend_param=list(
          title="Percent",
          legend_width=unit(4, "cm"), grid_height=unit(0.2, "cm"), 
          direction="horizontal", 
          title_position="topcenter", labels_gp=gpar(fontsize=9), title_gp=gpar(fontsize=10, fontface="plain")
        ),
        border_gp=gpar(col="black", lty=2),
        col=col_fun_cor,
        rect_gp=gpar(col="white", lwd=0.5),
        cluster_columns=T, cluster_rows=T,
        show_row_dend=T, show_column_dend=T,
        row_names_gp=gpar(fontsize=10), show_row_names=T, show_column_names=T,
        column_title_gp=gpar(fontsize=11),
        row_title_gp=gpar(fontsize=11),
        column_names_side="bottom",
        column_names_rot=90,
        column_names_gp=gpar(fontsize=10),
        height=unit(6, "mm")*nrow(matrix_df), # 0.35mm
        width=unit(6, "mm")*ncol(matrix_df),
        row_dend_width=unit(2, "cm"),
        column_dend_height=unit(2, "cm")
)

#####################################









##########################################################################
###
### Lag correlation between cell proportions and Superclusters' mean expr
###
##########################################################################



### Analyze the timing (lagged correlation) between cell proportions, 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)

### Load Annot.rough proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
prop_arranged=meta_info_combined %>%
  dplyr::select(age, Annot.rough, N_per_rough, N) %>%
  group_by(age, Annot.rough) %>%
  summarize(N_per_rough=sum(N_per_rough, na.rm=T),
            N=sum(N, na.rm=T)) %>%
  mutate(percent.rough=N_per_rough/N) %>%
  ungroup() %>%
  dplyr::select(-c(N_per_rough, N))

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.rough.txt.gz")
All_df_arranged=All_df %>%
  group_by(Annot.rough, age, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.rough","age"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(prop_arranged, All_df_arranged, by=c("age","Annot.rough")) %>%
  split(.$Annot.rough)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  percent.data=DF_merged_allAnnot.rough[[i]]$percent.rough
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, percent.data)) # freq is later than expr
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.rough=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CellFreq.vs.SuperclusterMeanExpr_Annot.rough.txt.gz", sep="\t")

### Plot lag correlation between cell freq and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CellFreq.vs.SuperclusterMeanExpr_Annot.rough.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
range(lag_mtx)
# rownames as Supercluster, colnames as Cell proportion, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Cell proportion vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_Cell.freq.vs.Supercluster.mean.expr_rough.pdf", width=3, height=2.85)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell proportions, 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)

### Load Annot.inter proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
prop_arranged=meta_info_combined %>%
  dplyr::select(age, Annot.inter, N_per_inter, N) %>%
  group_by(age, Annot.inter) %>%
  summarize(N_per_inter=sum(N_per_inter, na.rm=T),
            N=sum(N, na.rm=T)) %>%
  mutate(percent.inter=N_per_inter/N) %>%
  ungroup() %>%
  dplyr::select(-c(N_per_inter, N))

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.inter.txt.gz")
All_df_arranged=All_df %>%
  group_by(Annot.inter, age, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.inter","age"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(prop_arranged, All_df_arranged, by=c("age","Annot.inter")) %>%
  split(.$Annot.inter)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  percent.data=DF_merged_allAnnot.rough[[i]]$percent.inter
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, percent.data)) # freq is later than expr
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.inter=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CellFreq.vs.SuperclusterMeanExpr_Annot.inter.txt.gz", sep="\t")

### Plot lag correlation between cell freq and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CellFreq.vs.SuperclusterMeanExpr_Annot.inter.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
# remove Other.Eryth and Other.Mast
cor_mtx=cor_mtx[!grepl("Eryth|Mast",rownames(cor_mtx)),!grepl("Eryth|Mast",colnames(cor_mtx))]
lag_mtx=lag_mtx[!grepl("Eryth|Mast",rownames(lag_mtx)),!grepl("Eryth|Mast",colnames(lag_mtx))]
range(lag_mtx)
# rownames as Supercluster, colnames as Cell proportion, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-15, 0, 15), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Cell proportion vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_Cell.freq.vs.Supercluster.mean.expr_inter.pdf", width=6.5, height=3.25)
ht_rough
dev.off()

#####################################









##########################################################################
###
### Lag correlation among Superclusters' mean expr
###
##########################################################################



### Analyze the timing (lagged correlation) among the 12 superclusters' genes mean expr
### ... in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.rough.txt.gz")
All_df_arranged=All_df %>%
  group_by(Annot.rough, age, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.rough","age"), names_from="Supercluster", values_from="merged") %>%
  split(.$Annot.rough)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL.rough=list()
for (i in 1:length(All_df_arranged)) {
  expr.data=All_df_arranged[[i]] %>% .[,3:ncol(All_df_arranged[[i]])]
  Lag_Cor_Res_alljs=list()
  for (j in 1:ncol(expr.data)) {
    gene.series_j=expr.data[,j]
    Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series_j, gene.series)) # all the other superclusters later than the jth
    Lag_Cor_Res_alljs[[j]]=
      lapply(1:length(Lag_Cor_Res), function(idx) {
        data.frame(s2.later.than.s1=Lag_Cor_Res[[idx]][[1]], 
                   max.correlation=Lag_Cor_Res[[idx]][[2]],
                   Later=names(Lag_Cor_Res)[idx])
      }) %>%
      data.table::rbindlist(.) %>%
      mutate(Earlier=colnames(expr.data)[j])
  }
  Lag_Cor_Res_alljs=data.table::rbindlist(Lag_Cor_Res_alljs)
  Lag_Cor_Res_ALL.rough[[i]]=Lag_Cor_Res_alljs
  print(paste0(i, " out of ", length(All_df_arranged), " is done."))
}

### Plot lag correlation of mean expr among various Superclusters, across Annot.rough
library(ComplexHeatmap)
ht_=list()
for (i in 1:length(Lag_Cor_Res_ALL.rough)) {
  celltype_=names(Lag_Cor_Res_ALL.rough)[i]
  df_=Lag_Cor_Res_ALL.rough[[i]]
  
  cor_mtx=df_ %>% dplyr::select(Earlier, Later, max.correlation) %>%
    tidyr::pivot_wider(id_cols="Earlier", names_from="Later", values_from="max.correlation") %>%
    tibble::column_to_rownames("Earlier")
  superclusters_=names(table(rownames(cor_mtx)))
  superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
  cor_mtx=cor_mtx[superclusters_reorder,superclusters_reorder]
  lag_mtx=df_ %>% dplyr::select(Earlier, Later, s2.later.than.s1) %>%
    tidyr::pivot_wider(id_cols="Earlier", names_from="Later", values_from="s2.later.than.s1") %>%
    tibble::column_to_rownames("Earlier")
  superclusters_=names(table(rownames(lag_mtx)))
  superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
  lag_mtx=lag_mtx[superclusters_reorder,superclusters_reorder]
  range(lag_mtx)
  
  # rownames as Earlier, colnames as Later, so lag refer to colnames later than rownames
  col_fun=circlize::colorRamp2(c(-50, 0, 50), c("dodgerblue3", "white", "brown3"))
  ht_[[i]]=
    Heatmap(lag_mtx, 
            name=" ",
            heatmap_legend_param=list(
              title="lags",
              legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
              title_gp=gpar(fontface="plain", fontsize=10)
            ),
            col=col_fun, rect_gp=gpar(type="none"), 
            cell_fun=function(j, i, x, y, width, height, fill) {
              # darken the cells that do not correlate (cor<=0)
              if (cor_mtx[i, j]>0) {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="grey80", fill=NA, lwd=0.1))
              } else {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="grey80", fill="grey93", lwd=0.1))
              }
              # draw lags only in the cells with pos.cor
              if (cor_mtx[i, j]>0) {
                grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                            gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
              }
            }, 
            cluster_rows=FALSE, cluster_columns=FALSE,
            show_row_names=TRUE, show_column_names=FALSE,
            row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
            column_title=names(All_df_arranged)[i],
            column_title_gp=gpar(fontsize=10),
            height=unit(3,"mm")*nrow(lag_mtx),
            width=unit(2,"mm")*ncol(lag_mtx))
}
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(Reduce(`+`,ht_), heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SuperclusterMeanExpr.vs.SuperclusterMeanExpr_rough.pdf", width=9.6, height=1.75)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) among the 12 superclusters' genes mean expr
### ... in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.inter.txt.gz")
All_df_arranged=All_df %>%
  group_by(Annot.inter, age, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.inter","age"), names_from="Supercluster", values_from="merged") %>%
  split(.$Annot.inter)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL.rough=list()
for (i in 1:length(All_df_arranged)) {
  expr.data=All_df_arranged[[i]] %>% .[,3:ncol(All_df_arranged[[i]])]
  Lag_Cor_Res_alljs=list()
  for (j in 1:ncol(expr.data)) {
    gene.series_j=expr.data[,j]
    Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series_j, gene.series)) # all the other superclusters later than the jth
    Lag_Cor_Res_alljs[[j]]=
      lapply(1:length(Lag_Cor_Res), function(idx) {
        data.frame(s2.later.than.s1=Lag_Cor_Res[[idx]][[1]], 
                   max.correlation=Lag_Cor_Res[[idx]][[2]],
                   Later=names(Lag_Cor_Res)[idx])
      }) %>%
      data.table::rbindlist(.) %>%
      mutate(Earlier=colnames(expr.data)[j])
  }
  Lag_Cor_Res_alljs=data.table::rbindlist(Lag_Cor_Res_alljs)
  Lag_Cor_Res_ALL.rough[[i]]=Lag_Cor_Res_alljs
  print(paste0(i, " out of ", length(All_df_arranged), " is done."))
}

### Plot lag correlation of mean expr among various Superclusters, across Annot.inter
library(ComplexHeatmap)
ht_=list()
for (i in 1:length(Lag_Cor_Res_ALL.rough)) {
  celltype_=names(Lag_Cor_Res_ALL.rough)[i]
  df_=Lag_Cor_Res_ALL.rough[[i]]
  
  cor_mtx=df_ %>% dplyr::select(Earlier, Later, max.correlation) %>%
    tidyr::pivot_wider(id_cols="Earlier", names_from="Later", values_from="max.correlation") %>%
    tibble::column_to_rownames("Earlier")
  superclusters_=names(table(rownames(cor_mtx)))
  superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
  cor_mtx=cor_mtx[superclusters_reorder,superclusters_reorder]
  lag_mtx=df_ %>% dplyr::select(Earlier, Later, s2.later.than.s1) %>%
    tidyr::pivot_wider(id_cols="Earlier", names_from="Later", values_from="s2.later.than.s1") %>%
    tibble::column_to_rownames("Earlier")
  superclusters_=names(table(rownames(lag_mtx)))
  superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
  lag_mtx=lag_mtx[superclusters_reorder,superclusters_reorder]
  range(lag_mtx)
  
  # rownames as Earlier, colnames as Later, so lag refer to colnames later than rownames
  col_fun=circlize::colorRamp2(c(-50, 0, 50), c("dodgerblue3", "white", "brown3"))
  ht_[[i]]=
    Heatmap(lag_mtx, 
            name=" ",
            heatmap_legend_param=list(
              title="lags",
              legend_height=unit(3, "cm"), grid_width=unit(0.2, "cm"),
              direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
              title_gp=gpar(fontface="plain", fontsize=10)
            ),
            col=col_fun, rect_gp=gpar(type="none"), 
            cell_fun=function(j, i, x, y, width, height, fill) {
              # darken the cells that do not correlate (cor<=0)
              if (cor_mtx[i, j]>0) {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="grey80", fill=NA, lwd=0.1))
              } else {
                grid.rect(x=x, y=y, width=width, height=height,
                          gp=gpar(col="grey80", fill="grey93", lwd=0.1))
              }
              # draw lags only in the cells with pos.cor
              if (cor_mtx[i, j]>0) {
                grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                            gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
              }
            }, 
            cluster_rows=FALSE, cluster_columns=FALSE,
            show_row_names=TRUE, show_column_names=FALSE,
            row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
            column_title=names(All_df_arranged)[i],
            column_title_gp=gpar(fontsize=10),
            height=unit(3,"mm")*nrow(lag_mtx),
            width=unit(2,"mm")*ncol(lag_mtx))
}
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(Reduce(`+`,ht_), heatmap_legend_side="left")

pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SuperclusterMeanExpr.vs.SuperclusterMeanExpr_inter.pdf", width=20, height=1.75)
ht_rough
dev.off()

#####################################







##########################################################################
###
### Lag correlation between cell communication, proportions, metabolism, and functional changes
###
##########################################################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (source), 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange
CELLCHAT.WEIGHT_arranged=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Source) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.rough"="Source")
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.rough.txt.gz")
All_df_arranged=All_df
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(-age) %>%
  group_by(Annot.rough, Agecut, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.rough","Agecut"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(CELLCHAT.WEIGHT_arranged, All_df_arranged, by=c("Agecut","Annot.rough")) %>%
  split(.$Annot.rough)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  weight.data=DF_merged_allAnnot.rough[[i]]$Weight
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, weight.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.rough=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Source.vs.SuperclusterMeanExpr_Annot.rough.txt.gz", sep="\t")

### Plot lag correlation between communication strength and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Source.vs.SuperclusterMeanExpr_Annot.rough.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
range(lag_mtx)
# rownames as Supercluster, colnames as Communication.weight, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Sending signal strength vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SendingSignalStrength.vs.Supercluster.mean.expr_rough.pdf", width=3, height=2.85)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (target), 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange
CELLCHAT.WEIGHT_arranged=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Target) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.rough"="Target")
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load gene expr in 12 superclusters
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.rough.txt.gz")
All_df_arranged=All_df
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(-age) %>%
  group_by(Annot.rough, Agecut, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.rough","Agecut"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(CELLCHAT.WEIGHT_arranged, All_df_arranged, by=c("Agecut","Annot.rough")) %>%
  split(.$Annot.rough)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  weight.data=DF_merged_allAnnot.rough[[i]]$Weight
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, weight.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.rough=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Target.vs.SuperclusterMeanExpr_Annot.rough.txt.gz", sep="\t")

### Plot lag correlation between communication strength and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Target.vs.SuperclusterMeanExpr_Annot.rough.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.rough, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.rough", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
range(lag_mtx)
# rownames as Supercluster, colnames as Communication.weight, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Receiving signal strength vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_ReceivingSignalStrength.vs.Supercluster.mean.expr_rough.pdf", width=3, height=2.85)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (source), 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_inter_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange
CELLCHAT.WEIGHT_arranged=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Source) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.inter"="Source")
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load cell proportion
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.inter.txt.gz")
All_df_arranged=All_df
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(-age) %>%
  group_by(Annot.inter, Agecut, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.inter","Agecut"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(CELLCHAT.WEIGHT_arranged, All_df_arranged, by=c("Agecut","Annot.inter")) %>%
  split(.$Annot.inter)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  weight.data=DF_merged_allAnnot.rough[[i]]$Weight
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, weight.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.inter=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Source.vs.SuperclusterMeanExpr_Annot.inter.txt.gz", sep="\t")

### Plot lag correlation between communication strength and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Source.vs.SuperclusterMeanExpr_Annot.inter.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
# remove Other.Eryth and Other.Mast
cor_mtx=cor_mtx[!grepl("Eryth|Mast",rownames(cor_mtx)),!grepl("Eryth|Mast",colnames(cor_mtx))]
lag_mtx=lag_mtx[!grepl("Eryth|Mast",rownames(lag_mtx)),!grepl("Eryth|Mast",colnames(lag_mtx))]
range(lag_mtx)
# rownames as Supercluster, colnames as Communication.weight, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Sending signal strength vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SendingSignalStrength.vs.Supercluster.mean.expr_inter.pdf", width=6.5, height=3.25)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (target), 
### ... and gene expression in the 12 superclusters,
### ... in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_inter_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange
CELLCHAT.WEIGHT_arranged=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Target) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.inter"="Target")
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load cell proportion
All_df=data.table::fread("~/Project_PBMCage/Results/PBMC_results/Supercluster.genes_meanExpr.across.age_Annot.inter.txt.gz")
All_df_arranged=All_df
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(-age) %>%
  group_by(Annot.inter, Agecut, Supercluster) %>%
  summarize_at("merged", mean) %>%
  ungroup() %>%
  tidyr::pivot_wider(id_cols=c("Annot.inter","Agecut"), names_from="Supercluster", values_from="merged")

### Merge communication interaction and supercluster mean expr
DF_merged_allAnnot.rough=
  dplyr::full_join(CELLCHAT.WEIGHT_arranged, All_df_arranged, by=c("Agecut","Annot.inter")) %>%
  split(.$Annot.inter)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  expr.data=DF_merged_allAnnot.rough[[i]] %>% .[,4:ncol(DF_merged_allAnnot.rough[[i]])]
  weight.data=DF_merged_allAnnot.rough[[i]]$Weight
  Lag_Cor_Res=apply(expr.data, 2, function(gene.series) Lag_Cor(gene.series, weight.data))
  names(Lag_Cor_Res)=colnames(expr.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Supercluster=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.inter=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Target.vs.SuperclusterMeanExpr_Annot.inter.txt.gz", sep="\t")

### Plot lag correlation between communication strength and Superclusters' mean expr
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength_Target.vs.SuperclusterMeanExpr_Annot.inter.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="max.correlation") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(cor_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
cor_mtx=cor_mtx[superclusters_reorder,]
lag_mtx=DF_RES %>% dplyr::select(Supercluster, Annot.inter, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Supercluster", names_from="Annot.inter", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Supercluster")
superclusters_=names(table(rownames(lag_mtx)))
superclusters_reorder=superclusters_[c(9,11,12,5, 6,3,2,10, 4,1,8,7)]
lag_mtx=lag_mtx[superclusters_reorder,]
# remove Other.Eryth and Other.Mast
cor_mtx=cor_mtx[!grepl("Eryth|Mast",rownames(cor_mtx)),!grepl("Eryth|Mast",colnames(cor_mtx))]
lag_mtx=lag_mtx[!grepl("Eryth|Mast",rownames(lag_mtx)),!grepl("Eryth|Mast",colnames(lag_mtx))]
range(lag_mtx)
# rownames as Supercluster, colnames as Communication.weight, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Receiving signal strength vs. Expr.)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_ReceivingSignalStrength.vs.Supercluster.mean.expr_inter.pdf", width=6.5, height=3.25)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (both sending and receiving), 
### ... and Annot.rough proportions,
### ... in each celltype based on pseudobulk_rough
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_rough_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange of sending signals
CELLCHAT.WEIGHT_sending=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Source) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.rough"="Source")
# arrange of receiving signals
CELLCHAT.WEIGHT_receiving=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Target) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.rough"="Target")
# combine
CELLCHAT.WEIGHT_both=CELLCHAT.WEIGHT_sending %>% dplyr::rename("Sending"="Weight") %>%
  full_join(CELLCHAT.WEIGHT_receiving %>% dplyr::rename("Receiving"="Weight"), by=c("Agecut","Annot.rough"))
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load Annot.rough proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
All_df_arranged=meta_info_combined
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(Agecut, Annot.rough, N_per_rough, N) %>%
  group_by(Agecut, Annot.rough) %>%
  summarize(N_per_rough=sum(N_per_rough, na.rm=T),
            N=sum(N, na.rm=T)) %>%
  mutate(percent.rough=N_per_rough/N) %>%
  ungroup() %>%
  dplyr::select(-c(N_per_rough, N))

### Merge communication interaction and Annot.rough proportions
DF_merged_allAnnot.rough=
  dplyr::full_join(CELLCHAT.WEIGHT_both, All_df_arranged, by=c("Agecut","Annot.rough")) %>%
  split(.$Annot.rough)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.rough)) {
  freq.data=DF_merged_allAnnot.rough[[i]]$percent.rough
  weight.data=DF_merged_allAnnot.rough[[i]] %>% .[,3:4]
  Lag_Cor_Res=apply(weight.data, 2, function(orientation) Lag_Cor(freq.data, orientation)) # communication later than proportion
  names(Lag_Cor_Res)=colnames(weight.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.rough), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.rough)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Communication_orientation=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.rough=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength.vs.Cellfreq_Annot.rough.txt.gz", sep="\t")

### Plot lag correlation between communication strength and cell proportion
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength.vs.Cellfreq_Annot.rough.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Annot.rough, Communication_orientation, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Annot.rough", names_from="Communication_orientation", values_from="max.correlation") %>%
  tibble::column_to_rownames("Annot.rough")
lag_mtx=DF_RES %>% dplyr::select(Annot.rough, Communication_orientation, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Annot.rough", names_from="Communication_orientation", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Annot.rough")
range(lag_mtx)
# rownames as Annot.rough, colnames as Communication.orientation, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Signal strength vs. Cell proportion)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SignalStrength.vs.CellFreq_rough.pdf", width=2, height=2.25)
ht_rough
dev.off()

#####################################



### Analyze the timing (lagged correlation) between cell communication interaction strength (both sending and receiving), 
### ... and Annot.rough proportions,
### ... in each celltype based on pseudobulk_inter
#####################################
###

library(dplyr)

### Load cell communication interaction strength
CELLCHAT.WEIGHT=readRDS("~/Project_PBMCage/Tempt_RDS/Cellchat_inter_allresults.rds")[["CELLCHAT.WEIGHT"]]
# arrange of sending signals
CELLCHAT.WEIGHT_sending=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Source) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.inter"="Source")
# arrange of receiving signals
CELLCHAT.WEIGHT_receiving=CELLCHAT.WEIGHT %>%
  filter(Source!=Target) %>%
  group_by(Agecut, Target) %>%
  summarize_at("Weight", sum) %>%
  dplyr::rename("Annot.inter"="Target")
# combine
CELLCHAT.WEIGHT_both=CELLCHAT.WEIGHT_sending %>% dplyr::rename("Sending"="Weight") %>%
  full_join(CELLCHAT.WEIGHT_receiving %>% dplyr::rename("Receiving"="Weight"), by=c("Agecut","Annot.inter"))
agecuts_all=names(table(CELLCHAT.WEIGHT$Agecut))
agecuts_num=lapply(agecuts_all, function(x) c(as.integer(gsub("~.*","",x)):as.integer(gsub(".*~","",x))))

### Load Annot.inter proportion
meta_info_combined=read.delim("~/Project_PBMCage/Results/PBMC_results/Cell_frequency_values.txt", sep="\t")
All_df_arranged=meta_info_combined
All_df_arranged$Agecut=All_df_arranged$age
# arrange
for (i in 1:length(agecuts_all)) {
  All_df_arranged=All_df_arranged %>%
    mutate(Agecut=ifelse(Agecut %in% agecuts_num[[i]], agecuts_all[i], Agecut))
}
All_df_arranged=All_df_arranged %>%
  dplyr::select(Agecut, Annot.inter, N_per_inter, N) %>%
  group_by(Agecut, Annot.inter) %>%
  summarize(N_per_inter=sum(N_per_inter, na.rm=T),
            N=sum(N, na.rm=T)) %>%
  mutate(percent.inter=N_per_inter/N) %>%
  ungroup() %>%
  dplyr::select(-c(N_per_inter, N))

### Merge communication interaction and Annot.inter proportions
DF_merged_allAnnot.inter=
  dplyr::full_join(CELLCHAT.WEIGHT_both, All_df_arranged, by=c("Agecut","Annot.inter")) %>%
  split(.$Annot.inter)

### Run lag correlation
Lag_Cor=function(x, y) {
  corLag=ccf(x, y, lag.max=nrow(x), pl=FALSE)
  corLag_clean=corLag$acf[,,1]; names(corLag_clean)=corLag$lag[,,1]
  max.cor_=max(corLag_clean)
  lag_=as.integer(names(which.max(corLag_clean))) # time.series2 is [lag_] later than time.series1
  return(list(s2.later.than.s1=lag_, max.correlation=max.cor_))
} 
Lag_Cor_Res_ALL=list()
for (i in 1:length(DF_merged_allAnnot.inter)) {
  freq.data=DF_merged_allAnnot.inter[[i]]$percent.inter
  weight.data=DF_merged_allAnnot.inter[[i]] %>% .[,3:4]
  Lag_Cor_Res=apply(weight.data, 2, function(orientation) Lag_Cor(freq.data, orientation)) # communication later than proportion
  names(Lag_Cor_Res)=colnames(weight.data)
  Lag_Cor_Res_ALL[[i]]=Lag_Cor_Res
  print(paste0(i, " out of ", length(DF_merged_allAnnot.inter), " is done."))
}

### Arrange the Lag_Cor results
DF_RES=list()
names(Lag_Cor_Res_ALL)=names(DF_merged_allAnnot.inter)
for (i in 1:length(Lag_Cor_Res_ALL)) {
  tempt_df=Lag_Cor_Res_ALL[[i]] %>%
    lapply(., function(x) {
      data.frame(s2.later.than.s1=x[[1]], max.correlation=x[[2]])
    })
  tempt_df=
    lapply(1:length(tempt_df), function(idx) tempt_df[[idx]]=tempt_df[[idx]] %>% mutate(Communication_orientation=names(tempt_df)[idx])) %>%
    data.table::rbindlist() %>%
    mutate(Annot.inter=names(Lag_Cor_Res_ALL)[i])
  DF_RES[[i]]=tempt_df
}
DF_RES=data.table::rbindlist(DF_RES)
data.table::fwrite(DF_RES, "~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength.vs.Cellfreq_Annot.inter.txt.gz", sep="\t")

### Plot lag correlation between communication strength and cell proportion
library(ComplexHeatmap)
DF_RES=data.table::fread("~/Project_PBMCage/Results/PBMC_results/LagCor_CommunicationStrength.vs.Cellfreq_Annot.inter.txt.gz")
cor_mtx=DF_RES %>% dplyr::select(Annot.inter, Communication_orientation, max.correlation) %>%
  tidyr::pivot_wider(id_cols="Annot.inter", names_from="Communication_orientation", values_from="max.correlation") %>%
  tibble::column_to_rownames("Annot.inter")
lag_mtx=DF_RES %>% dplyr::select(Annot.inter, Communication_orientation, s2.later.than.s1) %>%
  tidyr::pivot_wider(id_cols="Annot.inter", names_from="Communication_orientation", values_from="s2.later.than.s1") %>%
  tibble::column_to_rownames("Annot.inter")
# remove Other.Eryth and Other.Mast
cor_mtx=cor_mtx[!grepl("Eryth|Mast",rownames(cor_mtx)),!grepl("Eryth|Mast",colnames(cor_mtx))]
lag_mtx=lag_mtx[!grepl("Eryth|Mast",rownames(lag_mtx)),!grepl("Eryth|Mast",colnames(lag_mtx))]
range(lag_mtx)
# rownames as Annot.inter, colnames as Communication.orientation, so lag refer to colnames later than rownames
col_fun=circlize::colorRamp2(c(-9, 0, 9), c("dodgerblue3", "white", "brown3"))
ht_rough=
  Heatmap(lag_mtx, 
          name=" ",
          heatmap_legend_param=list(
            title="lags\n(Signal strength vs. Cell proportion)",
            legend_height=unit(5, "cm"), grid_width=unit(0.2, "cm"),
            direction="vertical", title_position="leftcenter-rot", labels_gp=gpar(fontsize=9), 
            title_gp=gpar(fontface="plain", fontsize=10)
          ),
          col=col_fun, rect_gp=gpar(type="none"), 
          cell_fun=function(j, i, x, y, width, height, fill) {
            # darken the cells that do not correlate (cor<=0)
            if (cor_mtx[i, j]>0) {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill=NA, lwd=0.1))
            } else {
              grid.rect(x=x, y=y, width=width, height=height,
                        gp=gpar(col="grey80", fill="grey93", lwd=0.1))
            }
            # draw lags only in the cells with pos.cor
            if (cor_mtx[i, j]>0) {
              grid.circle(x=x, y=y, r=cor_mtx[i, j]/2 * min(unit.c(width, height)), 
                          gp=gpar(fill=col_fun(lag_mtx[i, j]), col=NA))
            }
          }, 
          cluster_rows=FALSE, cluster_columns=FALSE,
          show_row_names=TRUE, show_column_names=TRUE,
          row_names_gp=gpar(fontsize=9), column_names_gp=gpar(fontsize=9),
          height=unit(4,"mm")*nrow(lag_mtx),
          width=unit(4,"mm")*ncol(lag_mtx))
ht_opt$HEATMAP_LEGEND_PADDING=unit(3,"mm")
ht_rough=draw(ht_rough, heatmap_legend_side="left")
pdf("~/Project_PBMCage/Plots/PBMCage_Plots/Timing_SignalStrength.vs.CellFreq_inter.pdf", width=2.25, height=5.5)
ht_rough
dev.off()

#####################################